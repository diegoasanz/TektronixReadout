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
import pickle

from AbstractClasses.Channel_Caen import Channel_Caen
from AbstractClasses.Settings_Caen import Settings_Caen
from AbstractClasses.HV_Control import HV_Control
from AbstractClasses.Utils import *
from memory_profiler import profile


trig_rand_time = 0.2
wait_time_hv = 7

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
		self.hv_control = None
		self.utils = Utils()
		self.RemoveFiles()
		self.t0, self.t1, self.t2 = None, None, None
		self.p, self.pconv = None, None
		self.total_events = None
		self.written_events_sig, self.written_events_trig, self.written_events_aco = 0, 0, 0
		self.total_events_sig, self.total_events_trig, self.total_events_aco = 0, 0, 0
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_aco = 0, 0, 0
		self.total_merged_data_sig, self.total_merged_data_trig, self.total_merged_data_aco = 0, 0, 0
		self.doMerge = False
		self.min_measured_data = 0
		self.min_data_to_write = 0
		self.events_to_write = 0
		self.read_size = 0
		self.sig_written, self.trg_written, self.aco_written = 0, 0, 0
		self.session_written_events_sig, self.session_written_events_trg, self.session_written_events_aco = 0, 0, 0
		self.fins, self.fint, self.fina = None, None, None
		self.datas, self.datat, self.dataa = None, None, None

	def StartHVControl(self):
		wait_time_hv = 7
		if self.settings.do_hv_control:
			self.hv_control = HV_Control(self.settings)
			print 'Waiting {t} seconds for the HVClient to start...'.format(t=wait_time_hv)
			time.sleep(wait_time_hv)
			self.hv_control.CheckVoltage()

	def GetBaseLines(self):
		self.settings.SetupDigitiser(doBaseLines=True, signal=self.signal, trigger=self.trigger, ac=self.anti_co)
		self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		t0 = time.time()
		self.CreateEmptyFiles()
		self.GetWaveforms(events=1, stdin=True, stdout=False)
		if self.total_events >= 1:
			self.ReadBaseLines()
			t0 = time.time() - t0
			self.settings.RemoveBinaries()
			self.RemoveFiles()
		print 'Time getting base lines: {t} seconds'.format(t=t0)
		del t0

	def CreateEmptyFiles(self):
		self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger.ch), 'wb')
		# self.ft0.close()
		self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal.ch), 'wb')
		# self.fs0.close()
		if self.settings.ac_enable:
			self.fa0 = open('raw_wave{a}.dat'.format(a=self.anti_co.ch), 'wb')
			# self.fa0.close()

	def RemoveFiles(self):
		# used, for example, to remove old files that may have stayed due to crashes
		if self.fs0:
			if not self.fs0.closed:
				self.fs0.close()
				del self.fs0
		if self.ft0:
			if not self.ft0.closed:
				self.ft0.close()
				del self.ft0
		if self.settings.ac_enable:
			if self.fa0:
				if not self.fa0.closed:
					self.fa0.close()
					del self.fa0
		channels = [self.signal.ch, self.trigger.ch, self.anti_co.ch] if self.settings.ac_enable else [self.signal.ch, self.trigger.ch]
		for ch in channels:
			if os.path.isfile('raw_waves{c}.dat'.format(c=ch)):
				os.remove('raw_waves{c}.dat'.format(c=ch))
			if os.path.isfile('waves{c}.dat'.format(c=ch)):
				os.remove('waves{c}.dat'.format(c=ch))
		del channels

	def OpenFiles(self, mode='rb'):
		if not self.fs0:
			self.fs0 = open('raw_wave{s}.dat'.format(s=self.signal.ch), mode)
		if not self.ft0:
			self.ft0 = open('raw_wave{t}.dat'.format(t=self.trigger.ch), mode)
		if self.settings.ac_enable:
			if not self.fa0:
				self.fa0 = open('raw_wave{a}.dat'.format(a=self.anti_co.ch), mode)

	def CloseFiles(self):
		if self.ft0:
			self.ft0.close()
			if self.ft0.closed:
				del self.ft0
				self.ft0 = None
		if self.fs0:
			self.fs0.close()
			if self.fs0.closed:
				del self.fs0
				self.fs0 = None
		if self.settings.ac_enable:
			if self.fa0:
				self.fa0.close()
				if self.fa0.closed:
					del self.fa0
					self.fa0 = None

	def GetWaveforms(self, events=1, stdin=False, stdout=False):
		self.t1 = time.time()
		if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
		if events == 1:
			# while self.p.poll() is None:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('w')
			self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('t')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('q')
			self.p.stdin.flush()
			while self.p.poll() is None:
				continue
			if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
			self.ConcatenateBinaries2()
			self.CloseSubprocess('p', stdin=stdin, stdout=stdout)
			self.settings.RemoveBinaries()
		else:
			time.sleep(1)
			self.p.stdin.write('c')
			self.p.stdin.flush()
			time.sleep(1)
			self.p.stdin.write('W')
			self.p.stdin.flush()
			if self.settings.plot_waveforms:
				# time.sleep(1)
				self.p.stdin.write('P')
				self.p.stdin.flush()
			# time.sleep(1)
			self.p.stdin.write('s')
			self.p.stdin.flush()
			self.written_events_sig, self.written_events_trig, self.written_events_aco = 0, 0, 0
			# time.sleep(1)
			self.t2 = time.time()
			while self.p.poll() is None:
				if time.time() - self.t1 >= self.settings.time_calib:
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('c')
					self.p.stdin.flush()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
					if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
				elif self.written_events_sig + self.sig_written >= events:
					self.p.stdin.write('s')
					self.p.stdin.flush()
					self.settings.RemoveBinaries()
					self.p.stdin.write('q')
					self.p.stdin.flush()
					time.sleep(1)
					if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
				else:
					if self.settings.random_test and (time.time() - self.t2 > trig_rand_time):
						self.p.stdin.write('t')
						self.p.stdin.flush()
						self.t2 = time.time()
					self.ConcatenateBinaries2()
					if self.settings.do_hv_control:
						self.hv_control.UpdateHVFile()
					if not self.settings.simultaneous_conversion:
						self.settings.bar.update(int(min(self.written_events_sig + self.sig_written, self.settings.num_events)))
			del self.t1
			self.t1 = None
			self.CloseSubprocess('p', stdin=stdin, stdout=stdout)
			time.sleep(1)
			del self.t2
			self.t2 = None
			if self.settings.do_hv_control: self.hv_control.UpdateHVFile()
			# written_events_sig, written_events_trig, written_events_aco = self.ConcatenateBinaries2(written_events_sig, written_events_trig, written_events_aco, sig_written, trg_written, aco_written)
		self.total_events_sig = self.CalculateEventsWritten(self.signal.ch)
		self.total_events_trig = self.CalculateEventsWritten(self.trigger.ch)
		if self.settings.ac_enable:
			self.total_events_aco = self.CalculateEventsWritten(self.anti_co.ch)
		if self.total_events_sig == self.total_events_trig:
			self.total_events = self.total_events_sig
			if self.settings.ac_enable:
				if self.total_events_aco != self.total_events:
					print 'Events written are of different sizes. Missmatch'
					exit()
			del self.total_events_sig, self.total_events_trig, self.total_events_aco
			self.total_events_sig, self.total_events_trig, self.total_events_aco = None, None, None
		else:
			print 'Written events are of different sizes. Missmatch. quitting...'
			exit()

	def CloseSubprocess(self, pname, stdin=False, stdout=False):
		if pname == 'p':
			p = self.p
		elif pname == 'pconv':
			p = self.pconv

		pid = p.pid
		if stdin:
			p.stdin.close()
		if stdout:
			p.stdout.close()
		if p.wait() is None:
			print 'Could not terminate subprocess... forcing termination'
			p.kill()
			if p.wait() is None:
				print 'Could not kill subprocess... quitting'
				exit()
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print 'The subprocess is still running. Killing it with os.kill'
			os.kill(pid, 15)
			try:
				os.kill(pid, 0)
			except OSError:
				pass
			else:
				print 'The process does not die... quitting program'
				exit()
		del p, pid

		if pname == 'p':
			del self.p
			self.p = None
		elif pname == 'pconv':
			del self.pconv
			self.pconv = None

	def ConcatenateBinaries(self, fbase, fadd):
		fbin3 = open('tempMerge.dat', 'wb')
		f1 = open(fbase, 'rb')
		f2 = open(fadd, 'rb')
		fbin3.write(f1.read() + f2.read())
		fbin3.close()
		f1.close()
		f2.close()
		sys.stdout.flush()
		shutil.move('tempMerge.dat', fbase)

	def ConcatenateBinaries2(self):
		self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_aco = 0, 0, 0
		if os.path.isfile('wave{s}.dat'.format(s=self.signal.ch)) and os.path.isfile('wave{t}.dat'.format(t=self.trigger.ch)):
			self.session_measured_data_sig = int(os.path.getsize('wave{s}.dat'.format(s=self.signal.ch)))
			self.session_measured_data_trig = int(os.path.getsize('wave{t}.dat'.format(t=self.trigger.ch)))
			if self.settings.ac_enable and os.path.isfile('wave{a}.dat'.format(a=self.anti_co.ch)):
				self.session_measured_data_aco = int(os.path.getsize('wave{a}.dat'.format(a=self.anti_co.ch)))

		self.total_merged_data_sig = int(os.path.getsize('raw_wave{s}.dat'.format(s=self.signal.ch)))
		self.total_merged_data_trig = int(os.path.getsize('raw_wave{t}.dat'.format(t=self.trigger.ch)))
		if self.settings.ac_enable:
			self.total_merged_data_aco = int(os.path.getsize('raw_wave{a}.dat'.format(a=self.anti_co.ch)))
			self.doMerge = (self.session_measured_data_sig + self.sig_written * self.settings.struct_len > self.total_merged_data_sig) and (self.session_measured_data_trig + self.trg_written * self.settings.struct_len > self.total_merged_data_trig) and (self.session_measured_data_aco + self.aco_written * self.settings.struct_len > self.total_merged_data_aco)
		else:
			self.doMerge = (self.session_measured_data_sig + self.sig_written * self.settings.struct_len > self.total_merged_data_sig) and (self.session_measured_data_trig + self.trg_written * self.settings.struct_len > self.total_merged_data_trig)
		if self.doMerge:
			# self.OpenFiles(mode='ab')
			if self.settings.ac_enable:
				self.min_measured_data = min(self.session_measured_data_sig, self.session_measured_data_trig, self.session_measured_data_aco)
			else:
				self.min_measured_data = min(self.session_measured_data_sig, self.session_measured_data_trig)
			data_to_write_sig = self.min_measured_data - self.total_merged_data_sig + self.sig_written * self.settings.struct_len
			data_to_write_trg = self.min_measured_data - self.total_merged_data_trig + self.trg_written * self.settings.struct_len
			if self.settings.ac_enable:
				data_to_write_aco = self.min_measured_data - self.total_merged_data_aco + self.aco_written * self.settings.struct_len
				self.min_data_to_write = min(data_to_write_sig, data_to_write_trg, data_to_write_aco)
				del data_to_write_sig, data_to_write_trg, data_to_write_aco
			else:
				self.min_data_to_write = min(data_to_write_sig, data_to_write_trg)
				del data_to_write_sig, data_to_write_trg
			self.events_to_write = int(np.floor(self.min_data_to_write / float(self.settings.struct_len)))
			self.read_size = self.events_to_write * self.settings.struct_len

			with open('wave{s}.dat'.format(s=self.signal.ch), 'rb') as self.fins:
				self.fins.seek(self.written_events_sig * self.settings.struct_len, 0)
				self.datas = self.fins.read(self.read_size)
				# while not self.fins.closed:
				# 	self.fins.close()
			del self.fins
			self.fins = None
			with open('raw_wave{s}.dat'.format(s=self.signal.ch), 'ab') as self.fs0:
				self.fs0.write(self.datas)
				self.fs0.flush()
				# while not self.fs0.closed:
				# 	self.fs0.close()
			del self.fs0, self.datas
			self.fs0, self.datas = None, None

			with open('wave{t}.dat'.format(t=self.trigger.ch), 'rb') as self.fint:
				self.fint.seek(self.written_events_trig * self.settings.struct_len, 0)
				self.datat = self.fint.read(self.read_size)
				# while not self.fint.closed:
				# 	self.fint.close()
			del self.fint
			self.fint = None
			with open('raw_wave{t}.dat'.format(t=self.trigger.ch), 'ab') as self.ft0:
				self.ft0.write(self.datat)
				self.ft0.flush()
				# while not self.ft0.closed:
				# 	self.ft0.close()
			del self.ft0, self.datat
			self.ft0, self.datat = None, None

			if self.settings.ac_enable:
				with open('wave{a}.dat'.format(a=self.anti_co.ch), 'rb') as self.fina:
					self.fina.seek(self.written_events_aco * self.settings.struct_len, 0)
					self.dataa = self.fina.read(self.read_size)
					# while not self.fina.closed:
					# 	self.fina.close()
				del self.fina
				self.fina = None
				with open('raw_wave{a}.dat'.format(a=self.anti_co.ch), 'ab') as self.fa0:
					self.fa0.write(self.dataa)
					self.fa0.flush()
					# while not self.fa0.closed:
					# 	self.fa0.close()
				del self.fa0, self.dataa
				self.fa0, self.dataa = None, None

			self.written_events_sig += int(self.events_to_write)
			self.written_events_trig += int(self.events_to_write)
			if self.settings.ac_enable:
				self.written_events_aco += int(self.events_to_write)

			del self.events_to_write, self.read_size
			self.events_to_write, self.read_size = None, None

			# del self.datas, self.datat
			# self.datas, self.datat = None, None
			# if self.settings.ac_enable:
			# 	del self.dataa
			# 	self.dataa = None
			# self.CloseFiles()
		del self.doMerge
		self.doMerge = None

	def CalculateEventsWritten(self, ch):
		return int(round(float(os.path.getsize('raw_wave{c}.dat'.format(c=ch))) / float(self.settings.struct_len)))

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
		for i in xrange(10):
			condition_t = (np.abs(triggADCs - mean_t) < 3 * std_t)
			mean_t = np.extract(condition_t, triggADCs).mean()
			std_t = np.extract(condition_t, triggADCs).std()
			if self.settings.ac_enable:
				condition_ac = (np.abs(acADCs - mean_ac) < 3 * std_ac)
				mean_ac = np.extract(condition_ac, acADCs).mean()
				std_ac = np.extract(condition_ac, acADCs).std()
		# self.trigger.Correct_Base_Line(mean_volts=self.settings.ADC_to_Volts(mean_t, self.trigger), sigma_counts=std_t, settings=self.settings)
		self.trigger.Correct_Base_Line2(mean_adc=mean_t, sigma_adc=std_t, settings=self.settings)
		self.trigger.Correct_Threshold(sigma=std_t)
		self.settings.trig_base_line = np.multiply(self.trigger.base_line_u_adcs, self.settings.sigRes, dtype='f8')
		self.settings.trig_thr_counts = self.trigger.thr_counts
		if self.settings.ac_enable:
			# self.anti_co.Correct_Base_Line(mean_volts=self.settings.ADC_to_Volts(mean_ac, self.anti_co), sigma_counts=std_ac, settings=self.settings)
			self.anti_co.Correct_Base_Line2(mean_adc=mean_ac, sigma_adc=std_ac, settings=self.settings)
			self.anti_co.Correct_Threshold(sigma=std_ac)
			self.settings.ac_base_line = np.multiply(self.anti_co.base_line_u_adcs, self.settings.sigRes, dtype='f8')
			self.settings.ac_thr_counts = self.anti_co.thr_counts

		del ft, data_t, t, triggADCs, mean_t, std_t
		if self.settings.ac_enable:
			del fac, data_ac, ac, acADCs, mean_ac, std_ac

	# @profile(precision=12)
	def GetData(self):
		self.t0 = time.time()
		self.CreateEmptyFiles()
		self.CloseFiles()
		self.total_events = 0
		print 'Getting {n} events...'.format(n=self.settings.num_events)
		if self.settings.simultaneous_conversion:
			self.CreateRootFile()
		else:
			self.settings.CreateProgressBar(self.settings.num_events)
			self.settings.bar.start()
		if self.settings.ac_enable:
			self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal, trigger=self.trigger, ac=self.anti_co, events_written=self.total_events)
		else:
			self.settings.SetupDigitiser(doBaseLines=False, signal=self.signal, trigger=self.trigger, events_written=self.total_events)
		while self.total_events < self.settings.num_events:
			print "Calibrating ADC's...\r", ; sys.stdout.flush()
			self.sig_written = self.CalculateEventsWritten(self.signal.ch)
			self.trg_written = self.CalculateEventsWritten(self.trigger.ch)
			self.aco_written = 0
			if self.settings.ac_enable:
				self.aco_written = self.CalculateEventsWritten(self.anti_co.ch)
			# p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, close_fds=True)
			self.p = subp.Popen(['{p}/wavedump'.format(p=self.settings.wavedump_path), '{d}/WaveDumpConfig_CCD.txt'.format(d=self.settings.outdir)], bufsize=-1, stdin=subp.PIPE, stdout=subp.PIPE, close_fds=True)
			self.GetWaveforms(self.settings.num_events, stdin=True, stdout=True)
			# if time.time() - self.t0 >= 165:
			# 	self.total_events = self.settings.num_events
		self.CloseFiles()
		self.t0 = time.time() - self.t0
		if not self.settings.simultaneous_conversion:
			print 'Time getting {n} events: {t} seconds'.format(n=self.total_events, t=self.t0)
			self.settings.bar.finish()
		else:
			while self.pconv.poll() is None:
				continue
			self.CloseSubprocess('pconv', stdin=False, stdout=False)
		return self.total_events

	def CreateRootFile2(self):
		ac_ch = self.anti_co.ch if self.settings.ac_enable else -1
		ac_offset = self.anti_co.dc_offset_percent if self.settings.ac_enable else -1
		trig_th_in_volts = self.settings.ADC_to_Volts(self.settings.GetTriggerValueADCs(self.trigger), self.trigger)
		veto_value = self.anti_co.thr_counts if self.settings.ac_enable else 0
		control_hv = 1 if self.settings.do_hv_control and self.settings.simultaneous_conversion else 0
		settings_bin_path = self.settings.outdir + '/Runs/{f}.settings'.format(f=self.settings.filename)
		signal_bin_path = self.settings.outdir + '/Runs/{f}.signal'.format(f=self.settings.filename)
		trigger_bin_path = self.settings.outdir + '/Runs/{f}.trigger'.format(f=self.settings.filename)
		veto_bin_path = self.settings.outdir + '/Runs/{f}.veto'.format(f=self.settings.filename) if self.settings.ac_enable else ''
		bin_pat_sig = 'raw_wave{s}.dat'.format(s=self.signal.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_signal.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		bin_pat_trig = 'raw_wave{t}.dat'.format(t=self.trigger.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_trigger.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		bin_pat_ac = ''
		polarity = 1 if self.settings.bias >= 0 else -1
		if self.settings.ac_enable:
			bin_pat_ac = 'raw_wave{a}.dat'.format(a=self.anti_co.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_veto.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		self.pconv = subp.Popen(['python', 'AbstractClasses/Converter_Caen2.py', os.getcwd(), settings_bin_path, signal_bin_path, trigger_bin_path, veto_bin_path, bin_pat_sig, bin_pat_trig, bin_pat_ac], close_fds=True)
		del ac_ch, ac_offset, trig_th_in_volts, veto_value, control_hv, bin_pat_sig, bin_pat_trig, bin_pat_ac

	def CreateRootFile(self):
		ac_ch = self.anti_co.ch if self.settings.ac_enable else -1
		ac_offset = self.anti_co.dc_offset_percent if self.settings.ac_enable else -1
		trig_th_in_volts = self.settings.ADC_to_Volts(self.settings.GetTriggerValueADCs(self.trigger), self.trigger)
		veto_value = self.anti_co.thr_counts if self.settings.ac_enable else 0
		control_hv = 1 if self.settings.do_hv_control and self.settings.simultaneous_conversion else 0
		bin_pat_sig = 'raw_wave{s}.dat'.format(s=self.signal.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_signal.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		bin_pat_trig = 'raw_wave{t}.dat'.format(t=self.trigger.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_trigger.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		bin_pat_ac = ''
		polarity = 1 if self.settings.bias >= 0 else -1
		if self.settings.ac_enable:
			bin_pat_ac = 'raw_wave{a}.dat'.format(a=self.anti_co.ch) if self.settings.simultaneous_conversion else '{d}/Runs/{f}_veto.dat'.format(wd=os.getcwd(), d=self.settings.outdir, f=self.settings.filename)
		self.pconv = subp.Popen(['python', 'AbstractClasses/Converter_Caen.py', self.settings.outdir, os.getcwd(), self.settings.filename,
		                str(bin_pat_sig), str(bin_pat_trig), str(bin_pat_ac), str(self.settings.points),
		                str(self.settings.num_events),str(self.settings.struct_len), self.settings.struct_fmt,
		                str(self.settings.sigRes), str(self.signal.dc_offset_percent), str(self.trigger.dc_offset_percent),
		                str(ac_offset), str(self.settings.time_res), str(self.settings.post_trig_percent), str(trig_th_in_volts),
		                str(veto_value), str(self.settings.dig_bits), str(int(self.settings.simultaneous_conversion)),
		                str(self.settings.time_calib), str(control_hv), str(polarity)], close_fds=True)
		del ac_ch, ac_offset, trig_th_in_volts, veto_value, control_hv, bin_pat_sig, bin_pat_trig, bin_pat_ac

	def CloseHVClient(self):
		if self.settings.do_hv_control:
			self.hv_control.CloseClient()

	def SavePickles(self):
		with open('{d}/Runs/{f}.settings'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fs:
			pickle.dump(self.settings, fs, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}.signal'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fsig:
			pickle.dump(self.signal, fsig, pickle.HIGHEST_PROTOCOL)
		with open('{d}/Runs/{f}.trigger'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as ft:
			pickle.dump(self.trigger, ft, pickle.HIGHEST_PROTOCOL)
		if self.settings.ac_enable:
			with open('{d}/Runs/{f}.veto'.format(d=self.settings.outdir, f=self.settings.filename), 'wb') as fv:
				pickle.dump(self.anti_co, fv, pickle.HIGHEST_PROTOCOL)

	def PrintPlotLimits(self, ti=-5.12e-7, tf=4.606e-6, vmin=-0.7, vmax=0.05):
		print np.double([(tf-ti)/float(self.settings.time_res) +1, ti-self.settings.time_res/2.0,
		                      tf+self.settings.time_res/2.0, (vmax-vmin)/self.settings.sigRes, vmin, vmax])

def main():
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
	if auto:
		ccd.StartHVControl()
		ccd.GetBaseLines()
		ccd.SavePickles()
		written_events = ccd.GetData()
		ccd.settings.num_events = written_events
		ccd.SavePickles()
		ccd.settings.MoveBinaryFiles()
		ccd.settings.RenameDigitiserSettings()
		ccd.CloseHVClient()
		if not ccd.settings.simultaneous_conversion:
			ccd.CreateRootFile()
			while ccd.pconv.poll() is None:
				continue
			ccd.CloseSubprocess('pconv', stdin=False, stdout=False)

	# ccd.SetOutputFilesNames()
	# ccd.TakeTwoWaves()
	print 'Finished :)'
	sys.stdout.write('\a\a\a')
	sys.stdout.flush()

if __name__ == '__main__':
	main()
	# parser = OptionParser()
	# parser.add_option('-i', '--infile', dest='infile', default='', type='string',
	#                   help='Input configuration file. e.g. CAENCalibration.cfg')
	# parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	# parser.add_option('-a', '--automatic', dest='auto', default=False, help='Toggles automatic conversion and analysis afterwards', action='store_true')
	#
	# (options, args) = parser.parse_args()
	# infile = str(options.infile)
	# auto = bool(options.auto)
	# verb = bool(options.verb)
	# ccd = CCD_Caen(infile, verb)
	# ccd.StartHVControl()
	# ccd.GetBaseLines()
	# written_events = ccd.GetData()
	# ccd.settings.num_events = written_events
	# ccd.settings.MoveBinaryFiles()
	# ccd.CloseHVClient()
	# if auto:
	# 	if not ccd.settings.simultaneous_conversion:
	# 		ccd.pconv = ccd.CreateRootFile()
	# 		while ccd.pconv.poll() is None:
	# 			continue
	# 		ccd.CloseSubprocess(ccd.pconv, stdin=False, stdout=False)
	#
	# # ccd.SetOutputFilesNames()
	# # ccd.TakeTwoWaves()
	# print 'Finished :)'
	# sys.stdout.write('\a\a\a')
	# sys.stdout.flush()
