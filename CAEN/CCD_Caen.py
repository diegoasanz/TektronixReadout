#!/usr/bin/env python
import os
import shutil
import struct
import subprocess as subp
import sys
import time
from ConfigParser import ConfigParser
from optparse import OptionParser

import ROOT as ro
import numpy as np
import progressbar
import ipdb

from AbstractClasses.Channel_Caen import Channel_Caen
from AbstractClasses.Settings_Caen import Settings_Caen


# from DataAcquisition import DataAcquisition


class CCD_Caen:
	def __init__(self, infile='None', verbose=False):
		print 'Starting CCD program ...'
		self.infile = infile
		self.verb = verbose
		self.settings = Settings_Caen(self.infile, self.verb)
		self.settings.ReadInputFile()
		self.settings.Get_Calibration_Constants()
		self.settings.SetOutputFiles()

		self.signal = Channel_Caen(self.settings.sigCh, 'signal', self.verb)
		self.signal.Set_Channel(self.settings)
		self.trigger = Channel_Caen(self.settings.trigCh, 'trigger', self.verb)
		self.trigger.Set_Channel(self.settings)
		self.anti_co = None
		if self.settings.ac_enable:
			self.anti_co = Channel_Caen(self.settings.acCh, 'veto', self.verb)
			self.anti_co.Set_Channel(self.settings)

		self.fs0, self.ft0, self.fa0 = None, None, None

	def GetBaseLines(self):
		self.settings.SetupDigitiser(doBaseLines=True, signal=self.signal, trigger=self.trigger, ac=self.anti_co)
		p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE)
		t0 = time.time()
		self.CreateEmptyFiles()
		events_written = self.GetWaveforms(p, events=1)
		if events_written >= 1:
			self.CloseFiles()
			self.ReadBaseLines()
			t0 = time.time() - t0
			self.settings.RemoveBinaries()
		print 'Time getting base lines: {t} seconds'.format(t=t0)

	def CreateEmptyFiles(self):
		self.ft0 = file('raw_wave{t}.dat'.format(t=self.trigger.ch), 'wb')
		# self.ft0.close()
		self.fs0 = file('raw_wave{s}.dat'.format(s=self.signal.ch), 'wb')
		# self.fs0.close()
		if self.settings.ac_enable:
			self.fa0 = file('raw_wave{a}.dat'.format(a=self.anti_co.ch), 'wb')
			# self.fa0.close()

	def CloseFiles(self):
		self.ft0.close()
		self.fs0.close()
		if self.settings.ac_enable:
			self.fa0.close()

	def GetWaveforms(self, p, events=1, sig_written=0, trg_written=0, aco_written=0):
		t1 = time.time()
		if events == 1:
			# while p.poll() is None:
			self.settings.Delay(1)
			p.stdin.write('c')
			p.stdin.flush()
			self.settings.Delay(1)
			p.stdin.write('s')
			p.stdin.flush()
			# self.settings.Delay(1)
			p.stdin.write('P')
			p.stdin.flush()
			# self.settings.Delay(1)
			p.stdin.write('w')
			p.stdin.flush()
			# self.settings.Delay(1)
			p.stdin.write('t')
			p.stdin.flush()
			self.settings.Delay(1)
			p.stdin.write('s')
			p.stdin.flush()
			self.settings.Delay(1)
			p.stdin.write('q')
			p.stdin.flush()
			while p.poll() is None:
				continue
			written_events_sig, written_events_trig, written_events_aco = self.ConcatenateBinaries2(0, 0, 0, sig_written, trg_written, aco_written)
		else:
			self.settings.Delay(1)
			p.stdin.write('c')
			p.stdin.flush()
			self.settings.Delay(1)
			p.stdin.write('W')
			p.stdin.flush()
			# self.settings.Delay(1)
			p.stdin.write('P')
			p.stdin.flush()
			# self.settings.Delay(1)
			p.stdin.write('s')
			p.stdin.flush()
			written_events_sig = written_events_trig = written_events_aco = 0
			# self.settings.Delay(1)
			while p.poll() is None:
				if time.time() - t1 >= self.settings.time_calib:
					p.stdin.write('s')
					p.stdin.flush()
					p.stdin.write('c')
					p.stdin.flush()
					p.stdin.write('q')
					p.stdin.flush()
					self.settings.Delay(1)
				if written_events_sig + sig_written >= events:
					p.stdin.write('s')
					p.stdin.flush()
					p.stdin.write('q')
					p.stdin.flush()
					self.settings.Delay(1)
				written_events_sig, written_events_trig, written_events_aco = self.ConcatenateBinaries2(written_events_sig, written_events_trig, written_events_aco, sig_written, trg_written, aco_written)
				if not self.settings.simultaneous_conversion:
					self.settings.bar.update(int(min(written_events_sig + sig_written, self.settings.num_events)))
			self.settings.Delay(1)
			written_events_sig, written_events_trig, written_events_aco = self.ConcatenateBinaries2(written_events_sig, written_events_trig, written_events_aco, sig_written, trg_written, aco_written)
		total_events_sig = self.CalculateEventsWritten(self.signal.ch)
		total_events_trig = self.CalculateEventsWritten(self.trigger.ch)
		if self.settings.ac_enable:
			total_events_aco = self.CalculateEventsWritten(self.anti_co.ch)
		if total_events_sig == total_events_trig:
			total_events = total_events_sig
			if self.settings.ac_enable:
				if total_events_sig != total_events:
					print 'Events written are of different sizes. Missmatch'
					exit()
			return total_events
		else:
			print 'Written events are of different sizes. Missmatch. quitting...'
			exit()

	def ConcatenateBinaries(self, fbase, fadd):
		fbin3 = file('tempMerge.dat', 'wb')
		f1 = file(fbase, 'rb')
		f2 = file(fadd, 'rb')
		fbin3.write(f1.read() + f2.read())
		fbin3.close()
		f1.close()
		f2.close()
		sys.stdout.flush()
		shutil.move('tempMerge.dat', fbase)

	def ConcatenateBinaries2(self, written_events_sig, written_events_trig, written_events_aco, sig_written, trg_written, aco_written):
		self.fs0.flush()
		self.ft0.flush()
		if self.settings.ac_enable:
			self.fa0.flush()
		measured_signal_data = measured_trigger_data = measured_antico_data = 0
		if os.path.isfile('wave{s}.dat'.format(s=self.signal.ch)) and os.path.isfile('wave{t}.dat'.format(t=self.trigger.ch)):
			measured_signal_data = int(os.path.getsize('wave{s}.dat'.format(s=self.signal.ch)))
			measured_trigger_data = int(os.path.getsize('wave{t}.dat'.format(t=self.trigger.ch)))
			if self.settings.ac_enable and os.path.isfile('wave{a}.dat'.format(a=self.anti_co.ch)):
				measured_antico_data = int(os.path.getsize('wave{a}.dat'.format(a=self.anti_co.ch)))

		merged_signal_data = int(self.fs0.tell())
		merged_trigger_data = int(self.ft0.tell())
		doMerge = (measured_signal_data + sig_written * self.settings.struct_len > merged_signal_data) and (measured_trigger_data + trg_written * self.settings.struct_len > merged_trigger_data)
		if self.settings.ac_enable:
			merged_antico_data = int(self.fa0.tell())
			doMerge = doMerge and (measured_antico_data + aco_written * self.settings.struct_len > merged_antico_data)
		if doMerge:
			min_measured_data = min(measured_signal_data, measured_trigger_data)
			if self.settings.ac_enable:
				min_measured_data = min(min_measured_data, measured_antico_data)
			data_to_write_sig = min_measured_data - merged_signal_data + sig_written * self.settings.struct_len
			data_to_write_trg = min_measured_data - merged_trigger_data + trg_written * self.settings.struct_len
			min_data_to_write = min(data_to_write_sig, data_to_write_trg)
			if self.settings.ac_enable:
				data_to_write_aco = min_measured_data - merged_antico_data + aco_written * self.settings.struct_len
				min_data_to_write = min(min_data_to_write, data_to_write_aco)
			events_to_write = int(np.floor(min_data_to_write / float(self.settings.struct_len)))
			read_size = events_to_write * self.settings.struct_len

			fins = file('wave{s}.dat'.format(s=self.signal.ch), 'rb')
			fins.seek(written_events_sig * self.settings.struct_len, 0)
			datas = fins.read(read_size)
			self.fs0.write(datas)
			self.fs0.flush()
			fins.close()

			fint = file('wave{t}.dat'.format(t=self.trigger.ch), 'rb')
			fint.seek(written_events_trig * self.settings.struct_len, 0)
			datat = fint.read(read_size)
			self.ft0.write(datat)
			self.ft0.flush()
			fint.close()

			if self.settings.ac_enable:
				fina = file('wave{a}.dat'.format(a=self.anti_co.ch), 'rb')
				fina.seek(written_events_aco * self.settings.struct_len, 0)
				dataa = fina.read(read_size)
				self.fa0.write(dataa)
				self.fa0.flush()
				fina.close()
			written_events_sig += int(round(read_size/float(self.settings.struct_len)))
			written_events_trig += int(round(read_size/float(self.settings.struct_len)))
			written_events_aco += int(round(read_size/float(self.settings.struct_len)))
		self.fs0.flush()
		self.ft0.flush()
		if self.settings.ac_enable:
			self.fa0.flush()
		return written_events_sig, written_events_trig, written_events_aco

	def CalculateEventsWritten(self, ch):
		events_written = int(round(float(os.path.getsize('raw_wave{c}.dat'.format(c=ch))) / float(self.settings.struct_len)))
		return events_written

	def ReadBaseLines(self):
		ft = open('raw_wave{t}.dat'.format(t=self.trigger.ch), 'rb')
		ft.seek(0)
		data_t = ft.read(self.settings.struct_len)
		t = struct.Struct(self.settings.struct_fmt).unpack_from(data_t)
		triggADCs = np.array(t, 'H')
		mean_t = triggADCs.mean()
		std_t = triggADCs.std()
		if self.settings.ac_enable:
			fac = open('raw_wave{ac}.dat'.format(ac=self.anti_co.ch), 'rb')
			fac.seek(0)
			data_ac = fac.read(self.settings.struct_len)
			ac = struct.Struct(self.settings.struct_fmt).unpack_from(data_ac)
			acADCs = np.array(ac, 'H')
			mean_ac = acADCs.mean()
			std_ac = acADCs.std()
			# ipdb.set_trace(context=5)
		for i in xrange(10):
			condition_t = (np.abs(triggADCs - mean_t) < 3 * std_t)
			mean_t = np.extract(condition_t, triggADCs).mean()
			std_t = np.extract(condition_t, triggADCs).std()
			if self.settings.ac_enable:
				condition_ac = (np.abs(acADCs - mean_ac) < 3 * std_ac)
				mean_ac = np.extract(condition_ac, acADCs).mean()
				std_ac = np.extract(condition_ac, acADCs).std()
				# ipdb.set_trace(context=5)
		# self.trigger.Correct_Base_Line(mean_volts=self.settings.ADC_to_Volts(mean_t, self.trigger), sigma_counts=std_t, settings=self.settings)
		self.trigger.Correct_Base_Line2(mean_adc=mean_t, sigma_adc=std_t, settings=self.settings)
		self.trigger.Correct_Threshold(sigma=std_t)
		if self.settings.ac_enable:
			# self.anti_co.Correct_Base_Line(mean_volts=self.settings.ADC_to_Volts(mean_ac, self.anti_co), sigma_counts=std_ac, settings=self.settings)
			self.anti_co.Correct_Base_Line2(mean_adc=mean_ac, sigma_adc=std_ac, settings=self.settings)
			self.anti_co.Correct_Threshold(sigma=std_ac)

	def GetData(self):
		t0 = time.time()
		self.CreateEmptyFiles()
		total_events = 0
		print 'Getting {n} events...'.format(n=self.settings.num_events)
		if self.settings.simultaneous_conversion:
			pconv = self.CreateRootFile()
		else:
			self.settings.CreateProgressBar(self.settings.num_events)
			self.settings.bar.start()
		while total_events < self.settings.num_events:
			print "\nCalibrating ADC's..."
			if self.settings.ac_enable:
				self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal, trigger=self.trigger, ac=self.anti_co, events_written=total_events)
			else:
				self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal, trigger=self.trigger, events_written=total_events)
			sig_written = self.CalculateEventsWritten(self.signal.ch)
			trg_written = self.CalculateEventsWritten(self.trigger.ch)
			aco_written = 0
			if self.settings.ac_enable:
				aco_written = self.CalculateEventsWritten(self.anti_co.ch)
			# p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE)
			p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE)
			total_events = self.GetWaveforms(p, self.settings.num_events, sig_written, trg_written, aco_written)
		self.CloseFiles()
		t0 = time.time() - t0
		if not self.settings.simultaneous_conversion:
			print 'Time getting {n} events: {t} seconds'.format(n=total_events, t=t0)
			self.settings.bar.finish()
		else:
			while pconv.poll() is None:
				continue
		return total_events

	def CreateRootFile(self):
		ac_ch = self.anti_co.ch if self.settings.ac_enable else -1
		ac_offset = self.anti_co.dc_offset_percent if self.settings.ac_enable else -1
		trig_th_in_volts = self.settings.ADC_to_Volts(self.settings.GetTriggerValueADCs(self.trigger), self.trigger)
		veto_value = 0 if self.settings.ac_enable else self.anti_co.thr_counts
		p = subp.Popen(['python', 'AbstractClasses/Converter_Caen.py', self.settings.outdir, os.getcwd(), self.settings.filename,
		                str(self.signal.ch), str(self.trigger.ch), str(ac_ch), str(self.settings.points),
		                str(self.settings.num_events),str(self.settings.struct_len), self.settings.struct_fmt,
		                str(self.settings.sigRes), str(self.signal.dc_offset_percent), str(self.trigger.dc_offset_percent),
		                str(ac_offset), str(self.settings.time_res), str(self.settings.post_trig_percent), str(trig_th_in_volts),
		                str(veto_value), str(self.settings.dig_bits), str(int(self.settings.simultaneous_conversion)),
		                str(self.settings.time_calib)])
		return p

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


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i', '--infile', dest='infile', default='', type='string',
	                  help='Input configuration file. e.g. CAENCalibration.cfg')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')

	(options, args) = parser.parse_args()
	infile = str(options.infile)
	auto = bool(options.auto)
	verb = bool(options.verb)
	ccd = CCD_Caen(infile, verb)
	ccd.GetBaseLines()
	written_events = ccd.GetData()
	ccd.settings.num_events = written_events
	if auto:
		if not ccd.settings.simultaneous_conversion:
			pconv = ccd.CreateRootFile()
			while pconv.poll() is None:
				continue
		ccd.settings.MoveBinaryFiles()

	# ccd.SetOutputFilesNames()
	# ccd.TakeTwoWaves()
	print 'Finished :)'
	sys.stdout.write('\a\a\a')
	sys.stdout.flush()
