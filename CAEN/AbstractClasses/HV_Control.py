#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
from pykeyboard import PyKeyboard
from ConfigParser import ConfigParser
import subprocess as subp
import struct
import ROOT as ro
import shutil
from copy import deepcopy


# from DataAcquisition import DataAcquisition

class HV_Control:
	def __init__(self, settings):
		self.settings = settings
		self.hv_supply, self.ch, self.bias, self.current_limit, self.hot_start = settings.hv_supply, settings.hv_ch, settings.bias, settings.current_limit, settings.hot_start
		self.filename, self.dut = settings.filename, settings.dut
		self.Pics_folder_path = settings.pics_folder_path
		self.doControlHV = False if self.hv_supply == '' else True
		if self.Pics_folder_path == '':
			print 'Cannot control voltage because Pics folder (Micha) was not found XD'
			self.doControlHV = False
		self.process = None
		if self.doControlHV:
			self.LinkPicsFolder()
			self.CreateConfigFiles()
			if self.hot_start:
				self.process = subp.Popen(['HVClient.py', '-H'], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
				settings.Delay(3)
				self.process.stdin.write('yes\n')
				self.process.stdin.flush()
			else:
				print 'Only hot start has been implemented :P'
				# self.process = subp.Popen(['HVClient.py'], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)

	def LinkPicsFolder(self):
		if not os.path.isdir('Pics'):
			os.symlink(self.Pics_folder_path, 'Pics')

	def CreateConfigFiles(self):
		if not os.path.isdir('config'):
			os.mkdir('config')
		self.CreateHVClientConfig()

	def CreateHVClientConfig(self):
		num_supplies = 7
		supplies_num = range(1, num_supplies + 1)
		supplies_ids = {1: 'Keithley1', 2: 'Keithley2', 3: 'Keithley237', 4: 'Keithley6517', 5: '', 6: 'Keithley2657A', 7: 'ISEG-NHS-6220x'}
		conf_file = open('config/hv_{f}.cfg'.format(f=self.filename), 'w')

		conf_file.write('[Main]\n')
		conf_file.write('devices = {n}.\n'.format(n=num_supplies))
		conf_file.write('testbeam_name: = {f}\n'.format(f=self.filename))

		conf_file.write('\n[Names]\n')
		for s in supplies_num:
			if supplies_ids[s] == self.hv_supply:
				conf_file.write('CH{s}: {d}\n'.format(s=s, d=self.dut))
			else:
				conf_file.write('CH{s}: None\n'.format(s=s))

		for s in supplies_num:
			if supplies_ids[s] == self.hv_supply:
				if self.hv_supply == 'ISEG-NHS-6220x':
					conf_file.write('\n[HV{s}]\n'.format(s=s))
					conf_file.write('name: {n}\n'.format(n=supplies_ids[s]))
					conf_file.write('model: NHS-6220x\n')
					conf_file.write('module_name: ISEG\n')
					conf_file.write('nChannels: 6\n')
					conf_file.write('active_channels: [{ch}]\n'.format(ch=self.ch))
					conf_file.write('address: /dev/iseg\n')
					conf_file.write('# in V/s\n')
					conf_file.write('ramp: 10\n')
					conf_file.write('config_file: iseg.cfg\n')
			# TODO: write the other cases for the other supplies

		conf_file.close()
		self.UnlinkConfigFile('keithley.cfg')
		os.symlink('config/hv_{f}.cfg'.format(f=self.filename), 'keithley.cfg')

		if self.hv_supply == 'ISEG-NHS-6220x':
			self.CreateIsegConfigFile()

	def CreateIsegConfigFile(self):
		iseg_chs = 6
		conf_file = open('config/iseg_{f}.cfg'.format(f=self.filename), 'w')
		conf_file.write('[Names]\n')
		for ch in xrange(iseg_chs):
			if ch == self.ch:
				conf_file.write('CH{ch}: {d}\n'.format(ch=ch, d=self.dut))
			else:
				conf_file.write('CH{ch}: None\n'.format(ch=ch))
		for ch in xrange(iseg_chs):
			conf_file.write('\n[CH{ch}]\n'.format(ch=ch))
			conf_file.write('name: CH{ch}\n'.format(ch=ch))
			compliance = self.current_limit if ch == self.ch else 10e-6
			conf_file.write('compliance: {c}\n'.format(c=compliance))
			conf_file.write('measure_range: 10e-6\n')
			conf_file.write('bias: 0\n')
			min_bias = -1 if self.bias >= 0 else self.bias - 10
			max_bias = self.bias + 10 if self.bias >= 0 else 1
			if ch == self.ch:
				conf_file.write('min_bias: {m}\n'.format(m=min_bias))
				conf_file.write('max_bias: {m}\n'.format(m=max_bias))
			else:
				conf_file.write('min_bias: -1\n')
				conf_file.write('max_bias: 1\n')
		conf_file.close()
		self.UnlinkConfigFile('iseg.cfg')
		os.symlink('config/iseg_{f}.cfg'.format(f=self.filename), 'config/iseg.cfg')

	def UnlinkConfigFile(self, name):
		if os.path.isfile('config/{n}'.format(n=name)):
			if os.path.islink('config/{n}'.format(n=name)):
				os.unlink('config/{n}'.format(n=name))
			else:
				os.remove('config/{n}'.format(n=name))

	def SetupRootFile(self, arguments):
		if self.simultaneous_conversion:
			print '\n', str(arguments)
			print 'Start creating root file simultaneously with data taking'
		else:
			print str(arguments)
			print 'Start creating root file'
		self.raw_file = ro.TFile('{wd}/{d}/Runs/{r}.root'.format(wd=self.working_dir_location, d=self.run_dir_location, r=self.filename), 'RECREATE')
		self.raw_tree = ro.TTree(self.filename, self.filename)
		if self.doVeto:
			self.vetoBra = np.zeros(self.points, 'f8')
			self.vetoedBra = np.zeros(1, '?')
		self.eventBra = np.zeros(1, 'I')
		self.voltBra = np.zeros(self.points, 'f8')
		self.trigBra = np.zeros(self.points, 'f8')
		self.timeBra = np.zeros(self.points, 'f8')
		self.badShapeBra = np.zeros(1, dtype=np.dtype('int8'))  # signed char
		self.badPedBra = np.zeros(1, '?')
		self.raw_tree.Branch('event', self.eventBra, 'event/i')
		self.raw_tree.Branch('time', self.timeBra, 'time[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageSignal', self.voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.raw_tree.Branch('voltageTrigger', self.trigBra, 'voltageTrigger[{s}]/D'.format(s=self.points))
		if self.doVeto:
			self.raw_tree.Branch('voltageVeto', self.vetoBra, 'voltageVeto[{s}]/D'.format(s=self.points))
			self.raw_tree.Branch('vetoedEvent', self.vetoedBra, 'vetoedEvent/O')
		self.raw_tree.Branch('badShape', self.badShapeBra, 'badShape/B')  # signed char
		self.raw_tree.Branch('badPedestal', self.badPedBra, 'badPedestal/O')

	def OpenRawBinaries(self):
		self.fs = open('{wd}/raw_wave{s}.dat'.format(wd=self.working_dir_location, s=self.signal_ch), 'rb')
		self.ft = open('{wd}/raw_wave{t}.dat'.format(wd=self.working_dir_location, t=self.trigger_ch), 'rb')
		if self.doVeto:
			self.fa = open('{wd}/raw_wave{a}.dat'.format(wd=self.working_dir_location, a=self.anti_co_ch), 'rb')

	def GetBinariesWrittenEvents(self):
		self.signal_written_events = int(round(os.path.getsize('{d}/raw_wave{s}.dat'.format(d=self.working_dir_location, s=self.signal_ch)) / self.struct_len))
		self.trigger_written_events = int(round(os.path.getsize('{d}/raw_wave{t}.dat'.format(d=self.working_dir_location, t=self.trigger_ch)) / self.struct_len))
		if self.doVeto:
			self.anti_co_written_events = int(round(os.path.getsize('{d}/raw_wave{a}.dat'.format(d=self.working_dir_location, a=self.anti_co_ch)) / self.struct_len))

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)

	def CheckFilesSizes(self, ev):
		self.wait_for_data = True if (self.signal_written_events <= ev or self.trigger_written_events <= ev) else False
		if self.doVeto:
			self.wait_for_data = bool(self.wait_for_data or (self.anti_co_written_events <= ev))

	def WaitForData(self, ev, t1):
		while self.wait_for_data:
			if self.simultaneous_conversion:
				time_break = int(np.ceil(self.time_recal + 20))
				if time.time() - t1 > time_break:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=time_break)
					exit()
				if not self.fs.closed:
					self.fs.close()
				if not self.ft.closed:
					self.ft.close()
				if self.doVeto:
					if not self.fa.closed:
						self.fa.close()
				self.GetBinariesWrittenEvents()
				self.CheckFilesSizes(ev)
				if not self.wait_for_data:
					self.OpenRawBinaries()
			else:
				print 'The data is corrupted... exiting'
				exit()

	def ReadData(self, ev):
		self.fs.seek(ev * self.struct_len, 0)
		self.datas = self.fs.read(self.struct_len)
		self.ft.seek(ev * self.struct_len, 0)
		self.datat = self.ft.read(self.struct_len)
		if self.doVeto:
			self.fa.seek(ev * self.struct_len, 0)
			self.dataa = self.fa.read(self.struct_len)

	def CheckData(self):
		if not self.datas or not self.datat:
			print 'No event in signal or trigger files... exiting'
			exit()
		elif self.doVeto:
			if not self.dataa:
				print 'No event in veto trigger file... exiting'
				exit()

	def FillBranches(self, ev):
		self.eventBra.fill(ev)
		np.putmask(self.timeBra, np.bitwise_not(np.zeros(self.points, '?')), self.timeVect)
		np.putmask(self.voltBra, np.bitwise_not(np.zeros(self.points, '?')), self.sigVolts)
		np.putmask(self.trigBra, np.bitwise_not(np.zeros(self.points, '?')), self.trigVolts)
		if self.doVeto:
			np.putmask(self.vetoBra, np.bitwise_not(np.zeros(self.points, '?')), self.vetoVolts)
			self.vetoedBra.fill(self.vetoed_event)
		self.badShapeBra.fill(self.bad_shape_event)
		self.badPedBra.fill(self.bad_pedstal_event)

	def ConvertEvents(self):
		self.bar.start()

		for ev in xrange(self.num_events):
			t1 = time.time()
			self.CheckFilesSizes(ev)
			self.WaitForData(ev, t1)
			self.ReadData(ev)
			self.CheckData()
			s = struct.Struct(self.struct_fmt).unpack_from(self.datas)
			self.sigADC = np.array(s, 'H')
			self.sigVolts = self.ADC_to_Volts('signal')
			t = struct.Struct(struct_fmt).unpack_from(self.datat)
			self.trigADC = np.array(t, 'H')
			self.trigVolts = self.ADC_to_Volts('trigger')
			self.LookForTime0()
			self.timeVect = np.linspace(-self.trigPos * self.time_res, self.time_res * (self.points - 1 - self.trigPos), self.points, dtype='f8')
			if self.doVeto:
				ac = struct.Struct(struct_fmt).unpack_from(self.dataa)
				self.vetoADC = np.array(ac, 'H')
				self.vetoVolts = self.ADC_to_Volts('veto')
				self.vetoed_event = self.IsEventVetoed()
			self.DefineSignalBaseLineAndPeakPosition()
			self.bad_shape_event = self.IsEventBadShape()
			self.bad_pedstal_event = self.IsPedestalBad()
			self.FillBranches(ev)
			numFil = self.raw_tree.Fill()
			self.bar.update(ev + 1)

	def CloseAll(self):
		self.bar.finish()
		self.raw_file.Write()
		self.raw_file.Close()
		self.fs.close()
		self.ft.close()
		if self.doVeto:
			self.fa.close()
		self.t0 = time.time() - self.t0
		print 'Time creating root tree:', self.t0, 'seconds'
		exit()

	def LookForTime0(self):
		# ipdb.set_trace(context=7)
		guess_pos = int(round(self.points * (100.0 - self.post_trig_percent)/100.0))
		condition_trigg = np.array(np.abs(self.array_points - guess_pos) <= int(round(0.1e-6/self.time_res)), dtype='?')
		condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
		mean = np.extract(condition_no_trigg, self.trigVolts).mean()
		sigma = np.extract(condition_no_trigg, self.trigVolts).std()
		temp_trig_volts = np.copy(self.trigVolts)
		np.putmask(temp_trig_volts, condition_no_trigg, 100)
		volt_min_pos = temp_trig_volts.argmin()
		condition_trigg = np.bitwise_and(condition_trigg, np.array(self.array_points <= volt_min_pos))
		np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), 100)
		self.trigPos = np.abs(temp_trig_volts - self.trig_value).argmin()

	def IsEventVetoed(self):
		window_around_trigg = 50e-9
		condition_veto_base_line = np.array(np.abs(self.array_points - self.trigPos) > int(round(window_around_trigg/float(self.time_res))), dtype='?')
		condition_search = np.array(1 - condition_veto_base_line, dtype='?')
		mean = np.extract(condition_veto_base_line, self.vetoADC).mean()
		sigma = np.extract(condition_veto_base_line, self.vetoADC).std()
		vetoValNew = 4 * sigma if self.veto_value < 0.9 * 4 * sigma else self.veto_value
		veto_event = bool((np.extract(condition_search, self.vetoADC) - mean + vetoValNew).min() <= 0)
		return veto_event

	def IsPedestalBad(self):
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		sigma_volts = sigma * self.adc_res
		diff_volts = abs(int(self.sigADC[0]) - int(self.sigADC[self.trigPos])) * self.adc_res
		if sigma_volts >= 2e-3 or diff_volts >= 15e-3:
			return True
		else:
			return False

	def DefineSignalBaseLineAndPeakPosition(self):
		self.condition_base_line = np.array(self.array_points <= self.trigPos, dtype='?')
		# values 2.2 and 0.5 us come from shape of signal
		self.condition_peak_pos = np.array(np.abs(self.array_points - (2.2e-6/float(self.time_res) + self.trigPos)) <= 0.6e-6/float(self.time_res), dtype='?')

	def IsEventBadShape(self):
		# mean = np.extract(self.condition_base_line, self.sigADC).mean()
		sigma = np.extract(self.condition_base_line, self.sigADC).std()
		lim_inf = self.condition_peak_pos.argmax()
		lim_sup = self.points - self.condition_peak_pos[::-1].argmax() - 1
		peak_pos = self.sigADC.argmin()
		if lim_inf < peak_pos < lim_sup:
			# The event has a good shape
			return 0
		else:
			modified_adc = self.sigADC - sigma
			modified_adc[lim_inf] += 2*sigma
			modified_adc[lim_sup] += 2*sigma
			peak_pos = modified_adc.argmin()
			if lim_inf < peak_pos < lim_sup:
				# Can't tell if the event has a bad shape
				return -1
			else:
				# Event has bad shape
				return 1

	def IsPedestalNotFlat(self, signalADC, points, trigPos, time_res):
		array_points = np.arange(points, dtype=np.dtype('int32'))
		condition_base_line = np.array(array_points - trigPos <= 0, dtype='?')

	def ADC_to_Volts(self, sig_type):
		if sig_type == 'signal':
			adcs = self.sigADC
			offset = self.sig_offset
		elif sig_type == 'trigger':
			adcs = self.trigADC
			offset = self.trig_offset
		elif sig_type == 'veto':
			adcs = self.vetoADC
			offset = self.anti_co_offset
		else:
			print 'Wrong type. Exiting'
			exit()
		result = np.multiply(self.adc_res, np.add(adcs, np.multiply(2 ** self.dig_bits - 1.0, offset / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')
		return result

if __name__ == '__main__':
	run_dir_location = str(sys.argv[1])  # output directory location
	working_dir_location = str(sys.argv[2])  # current directory location, where the raw_wave files are
	filename = str(sys.argv[3])  # file name for the files
	signal_ch = int(sys.argv[4])  # caen channel for the ccd signal
	trigger_ch = int(sys.argv[5])  # caen channel for the trigger signal
	anti_co_ch = int(sys.argv[6])  # caen channel for the veto scintillator signal
	points = int(sys.argv[7])  # caen number of points in each event
	num_events = int(sys.argv[8])  # number of events to convert
	struct_len = int(sys.argv[9])  # structure length (obtained from struct_fmt and struct.calcsize(self.struct_fmt)
	struct_fmt = str(sys.argv[10])  # structure format of the data per event i.e. '@{p}H'.format(p=self.points)
	adc_res = np.double(sys.argv[11])  # adc resolution of the digitizer
	sig_offset = float(sys.argv[12])  # percentage offset for ccd signal (value within [-50, 50]) look at wavedump config file
	trig_offset = float(sys.argv[13])  # percentage offset for trigger signal (value within [-50, 50]) look at wavedump config file
	anti_co_offset = float(sys.argv[14])  # percentage offset for veto signal (value within [-50, 50]) look at wavedump config file
	time_res = np.double(sys.argv[15])  # time resolution of digitizer. i.e. 2e-9
	post_trig_percent = float(sys.argv[16])  # percentage of the acquired data after trigger (look at wave dump config file)
	trig_value = np.double(sys.argv[17])  # volts below base line for trigger
	veto_value = int(sys.argv[18])  # counts below base line on the veto signal for vetoing
	dig_bits = int(sys.argv[19])  # number of bits of the ADC i.e. 14
	simultaneous_conversion = bool(sys.argv[20] != '0')  # whether or not to do the simultaneous conversion while taking data
	time_recal = float(sys.argv[21])  # time between digitiser recalibrations

	converter = Converter_Caen(run_dir_location, working_dir_location, filename, signal_ch, trigger_ch, anti_co_ch, points,
	                           num_events, struct_len, struct_fmt, adc_res, sig_offset, trig_offset, anti_co_offset,
	                           time_res, post_trig_percent, trig_value, veto_value, dig_bits, simultaneous_conversion, time_recal)

	converter.SetupRootFile(sys.argv)
	converter.GetBinariesWrittenEvents()
	converter.OpenRawBinaries()
	converter.CreateProgressBar(converter.num_events)
	converter.ConvertEvents()
	converter.CloseAll()




