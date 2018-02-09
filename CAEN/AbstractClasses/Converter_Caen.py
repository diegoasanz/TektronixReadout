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

class Converter_Caen:
	def __init__(self, rundir, workingdir, filename, sigCh, trigCh, vetoCh, points, nevents, struct_len, struct_fmt,
	             adc_res, sig_dcop, trig_dcop, veto_dcop, time_res, post_trig, trig_val, veto_val, dig_bits, sim_conv, time_recal):
		self.run_dir_location = rundir  # output directory location
		self.working_dir_location = workingdir  # current directory location, where the raw_wave files are
		self.filename = filename  # file name for the files
		self.signal_ch = sigCh  # caen channel for the ccd signal
		self.trigger_ch = trigCh  # caen channel for the trigger signal
		self.anti_co_ch = vetoCh  # caen channel for the veto scintillator signal
		self.points = points  # caen number of points in each event
		self.num_events = nevents  # number of events to convert
		self.struct_len = struct_len  # structure length (obtained from struct_fmt and struct.calcsize(self.struct_fmt)
		self.struct_fmt = struct_fmt  # structure format of the data per event i.e. '@{p}H'.format(p=self.points)
		self.adc_res = adc_res  # adc resolution of the digitizer
		self.sig_offset = sig_dcop  # percentage offset for ccd signal (value within [-50, 50]) look at wavedump config file
		self.trig_offset = trig_dcop  # percentage offset for trigger signal (value within [-50, 50]) look at wavedump config file
		self.anti_co_offset = veto_dcop  # percentage offset for veto signal (value within [-50, 50]) look at wavedump config file
		self.time_res = time_res  # time resolution of digitizer. i.e. 2e-9
		self.post_trig_percent = post_trig  # percentage of the acquired data after trigger (look at wave dump config file)
		self.trig_value = trig_val  # volts below base line for trigger
		self.veto_value = veto_val  # counts below base line on the veto signal for vetoing
		self.dig_bits = dig_bits  # number of bits of the ADC i.e. 14
		self.simultaneous_conversion = sim_conv  # whether or not to do the simultaneous conversion while taking data
		self.time_recal = time_recal  # time between digitiser recalibrations

		self.doVeto = True if self.anti_co_ch != -1 else False
		self.array_points = np.arange(self.points, dtype=np.dtype('int32'))

		self.raw_file = None
		self.raw_tree = None
		self.eventBra = self.voltBra = self.trigBra = self.vetoBra = self.timeBra = self.vetoedBra = self.badShapeBra = self.badPedBra = None
		self.t0 = time.time()

		self.signal_written_events = self.trigger_written_events = self.anti_co_written_events = None
		self.fs = self.ft = self.fa = None
		self.wait_for_data = None

		self.datas = self.datat = self.dataa = None
		self.sigADC = self.trigADC = self.vetoADC = None
		self.sigVolts = self.trigVolts = self.vetoVolts = None
		self.trigPos = None
		self.timeVect = None
		self.vetoed_event = None
		self.bad_shape_event = None
		self.bad_pedstal_event = None
		self.condition_base_line = None
		self.condition_peak_pos = None

		self.bar = None

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




