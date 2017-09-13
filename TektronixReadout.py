#!/usr/bin/env python
import visa
import numpy as np
from struct import unpack
import time, os
from optparse import OptionParser
import progressbar
from DataAcquisition import DataAcquisition
from DeepAnalysis import DeepAnalysis
from SaveMeasurements import SaveMeasurements


class TektronixReadout:
	def __init__(self, ip='192.168.1.13', outdir='./', filename='', sig=2, trig=1, trigVal=-0.040, meas=100, waves=False, cont=False, pulser=False, ped=None, deep=False, verb=False):
		self.ip = ip
		self.outdir = outdir
		self.filename = filename
		self.sigCh = sig
		self.trigCh = trig
		self.trigVal = trigVal
		self.meas = meas
		self.doWaves = waves
		self.doCont = cont
		self.isPulser = pulser
		self.ped = ped
		self.verb = verb
		self.deep = deep
		# self.calcPed = None
		try:
			self.rm = visa.ResourceManager('@py')
			self.inst = self.rm.open_resource('TCPIP0::{IP}::inst0::INSTR'.format(IP=self.ip))
			print 'Connected to:', self.inst.ask('*idn?')
		except:
			print 'Check the IP address!'
		self.RoRate = 30
		self.inst.write('HEADer 1')
		self.inst.write('VERBose 1')
		self.wfmo, self.immed, self.nrpt, self.xincr, self.xunit, self.xzero, self.ymult, self.yoffs,self.yunit, self.yzero = None, None, None, None, None, None, None, None, None, None
		self.SetChAndTrigger()
		self.SetOutputFormat()
		self.daq = DataAcquisition(self) if not deep else DeepAnalysis(self)
		self.save = SaveMeasurements(self)
		if self.ped is not None:
			self.meas_points = 1000
			print 'Ready to take measurements...'
		else:
			self.meas_points = 90
			print 'Ready to measure the pedestal...'

	def SetChAndTrigger(self):
		# Reset registers
		self.inst.write('DLC')
		self.inst.write('*CLS')
		# Trigger
		self.inst.write('TRIGger:A:TYPe EDG')
		self.inst.write('TRIGger:A:MODe NORMal')
		self.inst.write('TRIGger:A:EDGE:COUP DC')
		self.inst.write('TRIGger:A:EDGE:SLOpe FALL')
		self.inst.write('TRIGger:A:EDGE:SOUrce CH{t}'.format(t=self.trigCh))
		self.inst.write('TRIGger:A:LEVel:CH{t} {trVal}'.format(t=self.trigCh, trVal=self.trigVal))
		# Horizontal
		self.inst.write('HORizontal:DELay:MODe OFF')
		self.inst.write('HORizontal:POSition 10')  # Sets the horizontal position trigger to 10%
		self.inst.write('HORizontal:RECOrdlength 1000')
		self.inst.write('HORizontal:SCAle 1E-6')
		# Vertical
		if self.isPulser:
			self.inst.write('CH{t}:TERmination 1E+6'.format(t=self.trigCh))
		else:
			self.inst.write('CH{t}:TERmination 50'.format(t=self.trigCh))
		self.inst.write('CH{t}:POSition 4'.format(t=self.trigCh))
		self.inst.write('CH{t}:COUPling DC'.format(t=self.trigCh))
		self.inst.write('CH{t}:SCAle 50E-3'.format(t=self.trigCh)) # 50E-3 # 10E-3 for attenuated calibrations
		self.inst.write('CH{t}:YUNits "V"'.format(t=self.trigCh))
		if not self.isPulser:
			self.inst.write('CH{t}:LABel "Trigger"'.format(t=self.trigCh))
		else:
			self.inst.write('CH{t}:LABel "Pulser"'.format(t=self.trigCh))
		self.inst.write('CH{s}:TERmination 50'.format(s=self.sigCh))
		self.inst.write('CH{s}:COUPling DC'.format(s=self.sigCh))
		if self.ped is None:
			self.inst.write('CH{s}:POSition -3'.format(s=self.sigCh))
			self.inst.write('CH{s}:SCAle 10E-3'.format(s=self.sigCh)) # 27.5E-3 #50E-3 #2E-3 for attenuated calibrations
		else:
			self.inst.write('CH{s}:POSition -4'.format(s=self.sigCh))
			self.inst.write('CH{s}:SCAle 27.5E-3'.format(s=self.sigCh)) # 27.5E-3 #50E-3 #2E-3 for attenuated calibrations
		self.inst.write('CH{s}:YUNits "V"'.format(s=self.sigCh))
		self.inst.write('CH{s}:LABel "Signal"'.format(s=self.sigCh))
		# Aqcuisition
		self.inst.write('ACQuire:MODe HIRes') # SAMple
		self.inst.write('ACQuire:NUMEnv 1')
		if not self.doCont:
			self.inst.write('ACQuire:STOPAfter SEQuence')
		else:
			self.inst.write('ACQuire:STOPAfter RUNSTOP')
		self.inst.write('ACQuire:STATE ON')
		# Measurement
		if not self.doWaves:
			self.inst.write('MEASUrement:CLEARSNapshot')
			self.inst.write('MEASUrement:GATing CURSor')
			self.inst.write('MEASUrement:IMMed:SOUrce CH{s}'.format(s=self.sigCh))
			self.inst.write('MEASUrement:IMMed:TYPe MAXimum')
			# Cursors
			self.inst.write('CURSor:FUNCtion WAVEform')
			self.inst.write('CURSor:SOUrce CH{s}'.format(s=self.sigCh))
			self.inst.write('CURSor:MODe TRACk')
			self.inst.write('CURSor:VBArs:UNIts SEConds')
			self.inst.write('CURSor:VBArs:POSITION1 500E-9')
			self.inst.write('CURSor:VBArs:POSITION2 5.00E-6')
			self.immed = str(self.inst.ask('MEASUrement:IMMed?'))
			if self.verb: print 'Measurement parameters from device:\n{mf}'.format(self.immed)

	def SetOutputFormat(self):
		self.inst.write('DATa INIT')
		self.inst.write('DATa:SOUrce CH{s}'.format(s=self.sigCh))
		self.inst.write('DATa:STARt 1')
		self.inst.write('DATa:STOP {ep}'.format(ep=self.meas_points))
		self.inst.write('DATa:ENCdg SRPbinary') # RPBinary, SRIbinary,SRPbinary, FAStest, FPbinary
		self.inst.write('DATa:WIDTH 2')
		self.wfmo = str(self.inst.ask('WFMOutpre?'))
		if self.verb and self.doWaves: print 'Measurement wave format:\n{wf}'.format(wf=self.wfmo)
		self.inst.write('HEADer 0')
		self.nrpt = int(self.inst.ask('WFMOutpre:NR_pt?'))
		self.xzero = float(self.inst.ask('WFMOutpre:XZEro?'))
		self.xincr = float(self.inst.ask('WFMOutpre:XINcr?'))
		self.xunit = str(self.inst.ask('WFMOutpre:XUNit?'))
		self.yzero = float(self.inst.ask('WFMOutpre:YZEro?'))
		self.ymult = float(self.inst.ask('WFMOutpre:YMUlt?'))
		self.yoffs = float(self.inst.ask('WFMOutpre:YOFf?'))
		self.yunit = str(self.inst.ask('WFMOutpre:YUNit?'))
		# for Curve
		# self.inst.values_format.is_binary = True
		# self.inst.values_format.datatype = 'H' # h -> SRI(2), H -> SRP(2), b -> SRI(1), B -> SRP(1)
		# self.inst.values_format.is_big_endian = False
		# self.inst.values_format.container = np.array
		# self.inst.write('HEADer 1')

	def TakeMeasurements(self):
		self.daq.TakeMeasurements()

	def SaveMeasurements(self):
		self.save.SaveMeasurements(self.daq.peak_values)

	def Quit(self):
		self.inst.before_close()
		self.inst.close()
		# sys.exit(0)
		exit()

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i', '--ip', dest='ip', default='192.168.1.13', type='string', help='IP address of the instrument. e.g. 192.168.1.13')
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string', help='Relative path to output directory. e.g. ./')
	parser.add_option('-f', '--filename', dest='filename', default='', type='string', help='Filename to save the data. e.g. output_2017_06')
	parser.add_option('-s', '--signal', dest='sig', default=2, type='int', help='The channel connected to the signal. e.g. 2')
	parser.add_option('-t', '--trigger', dest='trig', default=1, type='int', help='The channel connected to the trigger (e.g. scintillator) e.g. 1')
	parser.add_option('-l', '--trigval', dest='trigval', default=-0.07, type='float', help='The value of the trigger e.g. -0.04')
	parser.add_option('-m', '--measurements', dest='meas', default=100, type='int', help='The number of measurements to make e.g. 100')
	parser.add_option('-w', '--waveforms', dest='waves', default=False, action='store_true', help='Toggles full waveforms calculation instead of measurement acquisition. Overridden by deep analysis')
	parser.add_option('-c', '--continuous', dest='cont', default=False, action='store_true', help='Toggles continuous data sampling instead of successive single measurements. Overridden by deep analysis')
	parser.add_option('-p', '--pulser', dest='pulser', default=False, action='store_true', help='Toggles trigger from pulser')
	parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true', help='Toggles automatic analysis')
	parser.add_option('-d', '--deep', dest='deep', default=False, help='Toggles deep automatic analysis', action='store_true')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	ip = str(options.ip)
	outdir = str(options.outdir)
	filename = str(options.filename)
	sig = int(options.sig)
	trig = int(options.trig)
	trigval = float(options.trigval)
	meas = int(options.meas)
	deep = bool(options.complex)
	waves = bool(options.waves) if not deep else True
	cont = bool(options.cont) if not deep else False
	pulser = bool(options.pulser)
	automatic = bool(options.automatic)
	verb = bool(options.verb)
	pedMeas = 500
	if automatic or deep:
		zp = TektronixReadout(ip, outdir, filename, sig, trig, trigval, pedMeas, waves=True, cont=True, pulser=pulser, ped=None, deep=deep, verb=verb)
		zp.TakeMeasurements()
		ped = {'mean': zp.daq.ped_values.mean(axis=0, dtype='f'), 'sigma': zp.daq.ped_values.std(axis=0, dtype='f', ddof=1)}
		print 'Measured pedestal: Mean = {m}; Sigma = {s} for {pm} pedestal measurements'.format(m=ped['mean'], s=ped['sigma'], pm=pedMeas)
		z = TektronixReadout(ip, outdir, filename, sig, trig, trigval, meas, waves, cont, pulser, ped, deep, verb)
		z.TakeMeasurements()
		print 'Results:', str(z.daq.peak_values.mean(axis=0, dtype='f')), '+/-', str(z.daq.peak_values.std(axis=0, dtype='f', ddof=1)), 'V. For ', str(z.meas), 'measurements :)'
		z.SaveMeasurements()
		z.inst.before_close()
		z.inst.close()
