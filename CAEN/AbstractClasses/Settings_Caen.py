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
# from DataAcquisition import DataAcquisition


class Settings_Caen:
	def __init__(self, infile='None', verbose=False):
		self.infile = infile
		self.verb = verbose
		self.optlink = 1
		self.node = 0
		self.vme_b_addr = 32100000
		self.dig_bits = 14
		self.points = 2560
		self.post_trig_percent = 90
		self.num_events = 10
		self.time_calib = 300
		self.bias = 0
		self.input_range = 2.15
		self.calib_path = ''
		self.simultaneous_conversion = False
		self.time_res = 2e-9
		self.sigCh = 0
		self.trigCh = 1
		self.trig_base_line = -0.08
		self.trig_thr_counts = 35
		self.ac_enable = False
		self.acCh = 2
		self.ac_base_line = -0.08
		self.ac_thr_counts = 15
		self.outdir = '.'
		self.prefix = 'waves'
		self.suffix = 'default'
		self.sigRes = 0
		self.UpdateSignalResolution()
		self.filename = ''

		self.fit_signal_vcal_params = np.array([-1.21888705e-04, -8.96215025e-01], dtype='f8')
		self.fit_signal_vcal_params_errors = np.array([0.0001923, 0.00559264], dtype='f8')
		self.fit_charge_signal_params = np.array([-1.65417060e+01, -1.35735246e+05], dtype='f8')
		self.fit_charge_signal_params_errors = np.array([26.10228586,  847.0342207], dtype='f8')
		self.fit_vcal_signal_params, self.fit_vcal_signal_params_errors = None, None
		self.UpdateVcalVsSignal()

		self.struct_fmt = '@{p}H'.format(p=self.points)
		self.struct_len = struct.calcsize(self.struct_fmt)

		self.bar = None

	def ReadInputFile(self):
		parser = ConfigParser()
		if self.infile != 'None':
			if os.path.isfile(self.infile):
				print 'Reading input file: {f} ...'.format(f=self.infile)
				parser.read(self.infile)

				if parser.has_section('OPTILINK'):
					if parser.has_option('OPTILINK', 'link'):
						self.optlink = parser.getint('OPTILINK', 'link')
					if parser.has_option('OPTILINK', 'node'):
						self.node = parser.getint('OPTILINK', 'node')
					if parser.has_option('OPTILINK', 'vme_base_address'):
						self.vme_b_addr = parser.getint('OPTILINK', 'vme_base_address')

				if parser.has_section('RUN'):
					if parser.has_option('RUN', 'time'):
						self.points = int(np.ceil(parser.getfloat('RUN', 'time') * 1e-6 / self.time_res))
					if parser.has_option('RUN', 'post_trigger_percent'):
						self.post_trig_percent = parser.getint('RUN', 'post_trigger_percent')
					if parser.has_option('RUN', 'num_events'):
						self.num_events = parser.getint('RUN', 'num_events')
					if parser.has_option('RUN', 'time_calib'):
						self.time_calib = parser.getfloat('RUN', 'time_calib')
					if parser.has_option('RUN', 'sample_voltage'):
						self.bias = parser.getfloat('RUN', 'sample_voltage')
					if parser.has_option('RUN', 'input_range'):
						self.input_range = parser.getfloat('RUN', 'input_range')
					if parser.has_option('RUN', 'calib_path'):
						self.calib_path = parser.get('RUN', 'calib_path')
					if parser.has_option('RUN', 'simultaneous_conversion'):
						self.simultaneous_conversion = bool(parser.getint('RUN', 'simultaneous_conversion'))

				if parser.has_section('SIGNAL'):
					if parser.has_option('SIGNAL', 'channel'):
						self.sigCh = parser.getint('SIGNAL', 'channel')

				if parser.has_section('TRIGGER'):
					if parser.has_option('TRIGGER', 'channel'):
						self.trigCh = parser.getint('TRIGGER', 'channel')
					if parser.has_option('TRIGGER', 'base_line'):
						self.trig_base_line = parser.getfloat('TRIGGER', 'base_line')
					if parser.has_option('TRIGGER', 'thr_counts'):
						self.trig_thr_counts = parser.getint('TRIGGER', 'thr_counts')

				if parser.has_section('ANTICOINCIDENCE'):
					if parser.has_option('ANTICOINCIDENCE', 'channel'):
						self.ac_enable = True
						self.acCh = parser.getint('ANTICOINCIDENCE', 'channel')
					if parser.has_option('ANTICOINCIDENCE', 'base_line'):
						self.ac_base_line = parser.getfloat('ANTICOINCIDENCE', 'base_line')
					if parser.has_option('ANTICOINCIDENCE', 'thr_counts'):
						self.ac_thr_counts = parser.getint('ANTICOINCIDENCE', 'thr_counts')

				if parser.has_section('OUTPUT'):
					if parser.has_option('OUTPUT', 'dir'):
						self.outdir = parser.get('OUTPUT', 'dir')
					if parser.has_option('OUTPUT', 'prefix'):
						self.prefix = parser.get('OUTPUT', 'prefix')
					if parser.has_option('OUTPUT', 'suffix'):
						self.suffix = parser.get('OUTPUT', 'suffix')
				self.UpdateSignalResolution()

	def UpdateSignalResolution(self):
		self.sigRes = np.double(np.divide(np.double(self.input_range), (np.power(2.0, 14.0, dtype='f8') - 1)))

	def Get_Calibration_Constants(self):
		if self.calib_path == '':
			pass
		else:
			if os.path.isfile(self.calib_path):
				tempf = ro.TFile(self.calib_path, 'READ')
				fit_signal_vcal = tempf.Get('TFitResult-Signal_vs_Vcal-sig_vcal_fit')
				fit_charge_signal = tempf.Get('TFitResult-Charge_vs_Signal-q_sig_fit')
				self.fit_signal_vcal_params = np.array(fit_signal_vcal.Parameters(), dtype='f8')
				self.fit_charge_signal_params = np.array(fit_charge_signal.Parameters(), dtype='f8')
				self.UpdateVcalVsSignal()

	def UpdateVcalVsSignal(self):
		self.fit_vcal_signal_params = np.array([np.divide(-self.fit_signal_vcal_params[0], self.fit_signal_vcal_params[1], dtype='f8'),
		                                         np.divide(1.0, self.fit_signal_vcal_params[1], dtype='f8')], dtype='f8')
		self.fit_vcal_signal_params_errors = np.array([np.sqrt((self.fit_signal_vcal_params_errors[0] / self.fit_signal_vcal_params[1])**2 + (self.fit_signal_vcal_params[0] * self.fit_signal_vcal_params_errors[1] / (self.fit_signal_vcal_params[1] ** 2))**2)], dtype='f8')

	def SetOutputFiles(self):
		if not os.path.isdir(self.outdir):
			os.makedirs(self.outdir)
		if not os.path.isdir('{d}/Runs'.format(d=self.outdir)):
			os.makedirs('{d}/Runs'.format(d=self.outdir))

		def AddSuffix(string1):
			string1 += '_Pos' if self.bias >= 0 else '_Neg'
			string1 += '_{b}V'.format(b=abs(self.bias))
			if self.suffix != '':
				string1 += '_{s}'.format(s=self.suffix)
			return string1

		self.filename = '{p}_ccd'.format(p=self.prefix)
		self.filename = AddSuffix(self.filename)

	def Delay(self, ti=1.0):
		t0 = time.time()
		while time.time() - t0 < ti:
			continue
		return

	def SetupDigitiser(self, doBaseLines=False, signal=None, trigger=None, ac=None, events_written=0):
		print 'Creating digitiser CAEN V1730D configuration file... ', ; sys.stdout.flush()
		if doBaseLines:
			rfile = open('{d}/WaveDumpConfig_CCD_BL.txt'.format(d=self.outdir), 'w')
		else:
			rfile = open('{d}/WaveDumpConfig_CCD.txt'.format(d=self.outdir), 'w')
		rfile.write('[COMMON]')
		rfile.write('\n\n# open the digitezer')
		rfile.write('\nOPEN PCI {ol} {n} {ba}'.format(ol=int(self.optlink), n=int(self.node), ba=int(self.vme_b_addr)))
		rfile.write('\n\n# GNUPLOT path, normally /usr/bin/')
		rfile.write('\nGNUPLOT_PATH\t"/usr/bin/"')
		rfile.write('\n\n# output format can be BINARY or ASCII')
		rfile.write('\nOUTPUT_FILE_FORMAT\tBINARY')
		rfile.write('\n\n# if OUTPUT_FILE_HEADER is YES, the structure of the event has to be changed in self.struct_fmt to include the header')
		rfile.write('\nOUTPUT_FILE_HEADER NO')
		rfile.write('\n\n# specify the amount of samples to save. This defines the event window')
		rfile.write('\nRECORD_LENGTH\t{p}'.format(p=int(self.points)))
		rfile.write('\n\n# number of events to save in the file')
		if doBaseLines:
			rfile.write('\nMAX_NUM_EVENTS\t{n}'.format(n=1))
		else:
			rfile.write('\nMAX_NUM_EVENTS\t{n}'.format(n=self.num_events - events_written))
		rfile.write('\n\nTEST_PATTERN\tNO')
		rfile.write('\n\nENABLE_DES_MODE\tNO')
		rfile.write('\n\n# use external trigger. Options are: DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT')
		rfile.write('\nEXTERNAL_TRIGGER\tDISABLED')
		rfile.write('\n\n# specify maximum number of events to read out in one Block Transfer. Must be between 1 and 1023')
		rfile.write('\nMAX_NUM_EVENTS_BLT\t100')
		rfile.write('\n\n# the percentage of the amount of data stored after the trigger in the event window. Has an offset of ~1.6... only accepts integers')
		rfile.write('\nPOST_TRIGGER\t{pt}'.format(pt=int(round(self.post_trig_percent*0.9996 - 1.6384))))
		rfile.write('\n\n# number of events that have to be ready before readout when the IRQ is asserted. 0 means run continuously. 1023 is the maximum')
		rfile.write('\nUSE_INTERRUPT\t0')
		rfile.write('\n\n# type of the front panel LEMO connectors: NIM, TTL')
		rfile.write('\n\nFPIO_LEVEL\tNIM')
		rfile.write('\n\nSKIP_STARTUP_CALIBRATION\tNO')
		rfile.write('\n\nCHANNEL_TRIGGER\tDISABLED')

		sig_polarity = 'POSITIVE' if self.bias < 0 else 'NEGATIVE'

		rfile.write('\n\n# configuration for each channel [0] to [16], although it only has 8 channels ;)')
		for ch in xrange(16):
			rfile.write('\n\n[{ch}]'.format(ch=ch))
			if ch == signal.ch or ch == trigger.ch:
				rfile.write('\nENABLE_INPUT\tYES')
			elif not not ac:
				if ch == ac.ch:
					rfile.write('\nENABLE_INPUT\tYES')
				else:
					rfile.write('\nENABLE_INPUT\tNO')
			else:
				rfile.write('\nENABLE_INPUT\tNO')
			if ch == signal.ch:
				rfile.write('\nPULSE_POLARITY\t{sp}'.format(sp=sig_polarity))
				rfile.write('\nDC_OFFSET\t{o}'.format(o=signal.dc_offset_percent))
				rfile.write('\nCHANNEL_TRIGGER\tDISABLED')
			if ch == self.trigCh:
				rfile.write('\nPULSE_POLARITY\tNEGATIVE')
				rfile.write('\nDC_OFFSET\t{o}'.format(o=trigger.dc_offset_percent))
				if doBaseLines:
					rfile.write('\nCHANNEL_TRIGGER\tDISABLED')
				else:
					rfile.write('\nCHANNEL_TRIGGER\tACQUISITION_ONLY')
					rfile.write('\nTRIGGER_THRESHOLD\t{th}'.format(th=self.GetTriggerValueADCs(trigger.dc_offset_percent)))
			if not not ac:
				if ch == ac.ch:
					rfile.write('\nPULSE_POLARITY\tNEGATIVE')
					rfile.write('\nDC_OFFSET\t{o}'.format(o=ac.dc_offset_percent))
					rfile.write('\nCHANNEL_TRIGGER\tDISABLED')
		rfile.write('\n')
		rfile.close()
		print 'Done'

	def ADC_to_Volts(self, adcs, channel):
		return channel.ADC_to_Volts(adcs, self.sigRes)

	def GetTriggerValueADCs(self, offset):
		return int(round(np.divide(self.trig_thr_counts + 1 - offset / 50.0, self.sigRes)))

	def MoveBinaryFiles(self):
		print 'Moving binary files... ', ; sys.stdout.flush()
		shutil.move('raw_wave{chs}.dat'.format(chs=self.sigCh), '{d}/Runs/{f}_signal.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('raw_wave{cht}.dat'.format(cht=self.trigCh), '{d}/Runs/{f}_trigger.dat'.format(d=self.outdir, f=self.filename))
		if self.ac_enable:
			shutil.move('raw_wave{cha}.dat'.format(cht=self.acCh), '{d}/Runs/{f}_veto.dat'.format(d=self.outdir, f=self.filename))
		print 'Done'

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
	print 'bla'