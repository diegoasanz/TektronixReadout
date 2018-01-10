#!/usr/bin/env python
import visa
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ipdb
# from DataAcquisition import DataAcquisition


class CalibrationRuns:
	def __init__(self, ip='192.168.1.13', outdir='./', filename='', sig=2, sigres=0.21, trig=1, meas=100, calv=False, atten=0.0, bias=-400.0, verb=False):
		self.ip = ip
		self.outdir = outdir
		self.filename = filename
		self.sigCh = sig
		self.trigCh = trig
		self.calv = calv
		self.atten = atten
		self.bias = bias * 10.0 ** (-self.atten/20.0)
		self.sigRes = 0.21
		self.trigVal = 1 #self.bias/2.0 if self.atten > 0 else self.bias/10.0
		self.meas = meas
		self.verb = verb
		self.points = 10000
		self.windowLow = 500E-9
		self.windowHigh = 5E-6
		try:
			self.rm = visa.ResourceManager('@py')
			self.inst = self.rm.open_resource('TCPIP0::{IP}::inst0::INSTR'.format(IP=self.ip))
			print 'Connected to:', self.inst.ask('*idn?')
		except:
			print 'Check the IP address!'
		self.RoRate = 30
		self.inst.write('HEADer 1')
		self.inst.write('VERBose 1')
		self.wfmo, self.nrpt, self.xincr, self.xunit, self.xzero, self.ymult, self.yoffs,self.yunit, self.yzero = None, None, None, None, None, None, None, None, None
		# self.SetOutputFormatTwo()
		self.bindata1 = None
		self.bindata2 = None
		self.volts1, self.time1 = np.zeros(self.points, 'f'), np.zeros(self.points, 'f')
		self.volts2, self.time2 = np.zeros(self.points, 'f'), np.zeros(self.points, 'f')
		self.outString1 = ''
		self.outString2 = ''
		self.fileWaves = None
		self.fileWavesWriter = None
		self.iteration = 0
		# self.daq = DataAcquisition(self, self.verb)
		self.peak_values_waves = np.empty(0,'f')
		self.peak_values_measu = np.empty(0,'f')

	def WaitForSRQ(self):
		self.inst.write('ACQuire:STATE OFF')
		self.inst.write('ACQuire:STOPAfter SEQuence')
		self.inst.write('DESE 1') # enable device event status enable register
		self.inst.write('*ESE 1') # enable event status enable register
		self.inst.write('*SRE 32') # Generate a service request from the oscilloscope
		self.inst.write('ACQuire:STATE ON')
		self.inst.write('*OPC') # wait for interruption
		while True:
			ev = int(self.inst.ask('*ESR?')) # event in queue
			if int(self.inst.ask('EVENT?')) == 402: # message code for completion
				break
		return

	def SetOutputFilesNames(self):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))

		def AddSuffix(string1):
			string1 += '_Input' if self.calv else '_Output'
			string1 += '_Pos' if self.bias >= 0 else '_Neg'
			string1 += '_{b}mV'.format(b=abs(1000 * self.bias))
			string1 += '.csv'
			return string1

		self.outString3 = '{dir}/Runs/{f}'.format(dir=self.outdir, f=self.filename)
		self.outString3 = AddSuffix(self.outString3)

	def Quit(self):
		self.inst.before_close()
		self.inst.close()
		# sys.exit(0)
		exit()

	def SetOutputFormatTwo(self):
		self.points = 10000
		self.inst.write('DLC')
		self.inst.write('*CLS')
		# Trigger
		self.inst.write('SELect:CH{t} ON'.format(t=self.trigCh))
		self.inst.write('SELect:CH{s} ON'.format(s=self.sigCh))
		self.inst.write('TRIGger:A:TYPe EDG')
		self.inst.write('TRIGger:A:MODe NORMal')
		self.inst.write('TRIGger:A:EDGE:COUP DC')
		# if self.bias < 0:
		# 	self.inst.write('TRIGger:A:EDGE:SLOpe FALL')
		# else:
		self.inst.write('TRIGger:A:EDGE:SLOpe RISe')
		self.inst.write('TRIGger:A:EDGE:SOUrce CH{t}'.format(t=self.trigCh))
		self.inst.write('TRIGger:A:LEVel:CH{t} {trVal}'.format(t=self.trigCh, trVal=self.trigVal))
		# Horizontal
		self.inst.write('HORizontal:DELay:MODe OFF')
		self.inst.write('HORizontal:POSition 10')  # Sets the horizontal position trigger to 10%
		self.inst.write('HORizontal:RECOrdlength {p}'.format(p=self.points))
		self.inst.write('HORizontal:SCAle 1E-6')
		# Vertical
		# self.inst.write('CH{t}:TERmination 1E+6'.format(t=self.trigCh))
		self.inst.write('CH{t}:TERmination 50'.format(t=self.trigCh))
		self.inst.write('CH{s}:TERmination 50'.format(s=self.sigCh))
		if self.bias < 0:
			if not self.calv:
				print 'Negative bias: Positive signals'
				# if self.bias > -0.05:
				# 	self.inst.write('CH{s}:POSition -2'.format(s=self.sigCh))
				# elif self.bias > -0.5:
				# 	self.inst.write('CH{s}:POSition -3'.format(s=self.sigCh))
				# else:
				# 	self.inst.write('CH{s}:POSition -4'.format(s=self.sigCh))
				self.inst.write('CH{s}:POSition -4.5'.format(s=self.sigCh))
			else:
				self.inst.write('CH{s}:POSition 3.9'.format(s=self.sigCh))
# self.inst.write('CH{t}:POSition 4'.format(t=self.trigCh))
		else:
			if not self.calv:
				print 'Positive bias: Negative signals'
				# if self.bias < 0.05:
				# 	self.inst.write('CH{s}:POSition 2'.format(s=self.sigCh))
				# elif self.bias < 0.5:
				# 	self.inst.write('CH{s}:POSition 3'.format(s=self.sigCh))
				# else:
				# 	self.inst.write('CH{s}:POSition 4'.format(s=self.sigCh))
				self.inst.write('CH{s}:POSition 4.5'.format(s=self.sigCh))
			else:
				self.inst.write('CH{s}:POSition -3.9'.format(s=self.sigCh))
		self.inst.write('CH{t}:POSition -4'.format(t=self.trigCh))
		# if self.bias >= 0:
		self.inst.write('CH{t}:COUPling DC'.format(t=self.trigCh))
		self.inst.write('CH{s}:COUPling DC'.format(s=self.sigCh))
		# self.inst.write('CH{t}:SCAle {v}'.format(t=self.trigCh, v=abs(self.bias/4)))
		self.inst.write('CH{t}:SCAle 0.25'.format(t=self.trigCh))
		self.inst.write('CH{s}:SCAle {v}'.format(s=self.sigCh, v=abs(self.sigRes)))
		# if abs(self.bias) >= 1:
		# 	value = int(self.bias * 1000.0 / 4.0)
		# 	self.inst.write('CH{t}:SCAle {v}E-3'.format(t=self.trigCh, v=value)) # 50E-3 # 10E-3 for attenuated calibrations
		# else:
		# 	self.inst.write('CH{t}:SCAle {v:1.2f}'.format(t=self.trigCh, v=self.bias)) # 50E-3 # 10E-3 for attenuated calibrations
		# v_sig = int(self.sigRes * 1000)
		# self.inst.write('CH{s}:SCAle {r}E-3'.format(s=self.sigCh, r=v_sig)) # 50E-3 # 10E-3 for attenuated calibrations
		self.inst.write('CH{t}:YUNits "V"'.format(t=self.trigCh))
		self.inst.write('CH{s}:YUNits "V"'.format(s=self.sigCh))
		self.inst.write('CH{t}:LABel "CalTrigger"'.format(t=self.trigCh))
		self.inst.write('CH{s}:LABel "Signal"'.format(s=self.sigCh))
		# Aqcuisition
		self.inst.write('ACQuire:MODe HIRes') # SAMple
		self.inst.write('ACQuire:NUMEnv 1')
		self.inst.write('ACQuire:STOPAfter RUNSTOP')
		self.inst.write('ACQuire:STATE ON')
		self.nrpt = np.zeros(4, 'I')
		self.xzero = np.zeros(4, 'f8')
		self.xincr = np.zeros(4, 'd')
		self.xunit = np.empty(4, 'U')
		self.yzero = np.zeros(4, 'f8')
		self.ymult = np.zeros(4, 'f8')
		self.yoffs = np.zeros(4, 'f8')
		self.yunit = np.empty(4, 'U')

	def SetupOutputChannel(self, ch=1):
		self.inst.write('DATa INIT')
		self.inst.write('DATa:SOUrce CH{s}'.format(s=ch))
		self.inst.write('DATa:STARt 1')
		self.inst.write('DATa:STOP {p}'.format(p=self.points))
		self.inst.write('DATa:ENCdg SRPbinary') # RPBinary, SRIbinary,SRPbinary, FAStest, FPbinary
		self.inst.write('DATa:WIDTH 2')
		self.wfmo = str(self.inst.ask('WFMOutpre?'))
		self.inst.write('HEADer 0')
		self.nrpt = int(self.inst.ask('WFMOutpre:NR_pt?'))
		self.xzero[ch-1] = float(self.inst.ask('WFMOutpre:XZEro?'))
		self.xincr[ch-1] = float(self.inst.ask('WFMOutpre:XINcr?'))
		self.xunit[ch-1] = str(self.inst.ask('WFMOutpre:XUNit?'))
		self.yzero[ch-1] = float(self.inst.ask('WFMOutpre:YZEro?'))
		self.ymult[ch-1] = float(self.inst.ask('WFMOutpre:YMUlt?'))
		self.yoffs[ch-1] = float(self.inst.ask('WFMOutpre:YOFf?'))
		self.yunit[ch-1] = str(self.inst.ask('WFMOutpre:YUNit?'))

	def TakeTwoWaves(self):
		self.SetOutputFormatTwo()
		self.SetupOutputChannel(1)
		widgets = [
			progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			# ' ', progressbar.AdaptativeTransferSpeed()
		]
		bar = progressbar.ProgressBar(widgets=widgets, maxval=self.meas)
		bar.start()
		t0 = time.time()
		for iter in xrange(self.meas):
			self.DoMeasurementTwo()
			self.WriteCSVRowTwo(iter)
			# self.peak_values = np.append(self.peak_values, self.peak_val)
			# self.peak_times = np.append(self.peak_times, self.peak_pos)
			bar.update(iter+1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1-t0), 's'

	def DoMeasurementTwo(self):
		self.ReadTwoData()
		self.FormatTwoData()

	def ReadTwoData(self):
		self.inst.write('HEADer 0')
		self.WaitForSRQ()
		self.SetupOutputChannel(int(self.trigCh))
		self.bindata1 = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.SetupOutputChannel(int(self.sigCh))
		self.bindata2 = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.inst.write('HEADer 1')
	
	def FormatTwoData(self):
		tp = self.trigCh
		sp = self.sigCh
		self.volts1 = self.yzero[tp-1] + self.ymult[tp-1] * (self.bindata1 - self.yoffs[tp-1])
		self.volts2 = self.yzero[sp-1] + self.ymult[sp-1] * (self.bindata2 - self.yoffs[sp-1])
		self.time1 = np.linspace(self.xzero[tp-1], self.xzero[tp-1] + self.xincr[tp-1] * (self.nrpt - 1), self.nrpt)
		self.time2 = np.linspace(self.xzero[sp-1], self.xzero[sp-1] + self.xincr[sp-1] * (self.nrpt - 1), self.nrpt)

	def WriteCSVRowTwo(self, iteration):
		iterArray = np.full(len(self.volts1), iteration)
		# ipdb.set_trace()
		event = np.array((iterArray, self.time1, self.volts1, self.volts2), 'f').transpose()
		self.fileWavesWriter.writerows(event)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i', '--ip', dest='ip', default='169.254.2.226', type='string', help='IP address of the instrument. e.g. 192.168.1.13')
	parser.add_option('-f', '--filename', dest='filename', default='waves_cal', help='Stem of the name for the files. e.g. waves_cal')
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string', help='Relative path to output directory where the Runs folder is e.g. ./')
	parser.add_option('-s', '--signal', dest='sig', default=1, type='int', help='The channel connected to the signal. e.g. 1')
	parser.add_option('-r', '--sigres', dest='sigres', default=0.21, type='float', help='The value of the voltage per division. No longer used e.g. 0.210')
	parser.add_option('-t', '--trigger', dest='trig', default=2, type='int', help='The channel connected to the TTL trigger (has to be attenuated with 6dB) e.g. 2')
	parser.add_option('-m', '--measurements', dest='meas', default=500, type='int', help='The number of measurements to make e.g. 100')
	parser.add_option('-c', '--calvoltage', dest='calv', default=False, action='store_true', help='Saving the calibration pulse from the pulser')
	parser.add_option('-b', '--bias', dest='bias', default=-400, type='float', help='The calibration voltage set in the pulser in volts used e.g. 0.05024')
	parser.add_option('-a', '--attenuator', dest='atten', default=0, type='int', help='The total attenuation used in dB e.g. 20')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	ip = str(options.ip)
	filename = str(options.filename)
	outdir = str(options.outdir)
	sig = int(options.sig)
	sigres = float(options.sigres)
	trig = int(options.trig)
	meas = int(options.meas)
	bias = float(options.bias)
	atten = float(options.atten)
	calv = bool(options.calv)
	verb = bool(options.verb)
	z = CalibrationRuns(ip, outdir, filename, sig, sigres, trig, meas, calv, atten, bias, verb)
	z.SetOutputFilesNames()
	z.fileWaves = open(z.outString3, 'wb')
	z.fileWavesWriter = csv.writer(z.fileWaves, delimiter=',')
	z.TakeTwoWaves()
	print 'Finished :)'
	z.fileWaves.close()
	z.inst.before_close()
	z.inst.close()
	sys.stdout.write('\a\a\a')
	sys.stdout.flush()
