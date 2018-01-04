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


class TektronixReadout:
	def __init__(self, ip='192.168.1.13', outdir='./', filename='', sig=2, trig=1, trigVal=-0.040, meas=100, waves=False, cont=False, pulser=False, bias=-400, verb=False):
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
		self.bias = bias
		self.verb = verb
		self.points = 1000
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
		self.SetChAndTrigger()
		self.SetOutputFormat()
		self.bindata = None
		self.bindata1 = None
		self.bindata2 = None
		self.bindata3 = None
		self.volts, self.time, self.peak_val, self.peak_pos = np.zeros(self.points, 'f'), np.zeros(self.points, 'f'), None, None
		self.volts1, self.time1 = np.zeros(self.points, 'f'), np.zeros(self.points, 'f')
		self.volts2, self.time2 = np.zeros(self.points, 'f'), np.zeros(self.points, 'f')
		self.volts3, self.time3 = np.zeros(self.points, 'f'), np.zeros(self.points, 'f')
		self.peak_values = np.empty(0, 'f')
		self.peak_times = np.empty(0, 'f')
		self.outString1 = ''
		self.outString2 = ''
		self.outString3 = ''
		self.fileWaves = None
		self.fileWavesWriter = None
		self.iteration = 0
		# self.daq = DataAcquisition(self, self.verb)
		self.peak_val2 = None
		self.peak_values_waves = np.empty(0,'f')
		self.peak_values_measu = np.empty(0,'f')

	def SetChAndTrigger(self):
		# Reset registers
		self.inst.write('DLC')
		self.inst.write('*CLS')
		# Trigger
		self.inst.write('SELect:CH{s} ON'.format(s=self.sigCh))
		self.inst.write('SELect:CH{t} ON'.format(t=self.trigCh))
		self.inst.write('TRIGger:A:TYPe EDG')
		self.inst.write('TRIGger:A:MODe NORMal')
		self.inst.write('TRIGger:A:EDGE:COUP DC')
		self.inst.write('TRIGger:A:EDGE:SLOpe FALL')
		self.inst.write('TRIGger:A:EDGE:SOUrce CH{t}'.format(t=self.trigCh))
		self.inst.write('TRIGger:A:LEVel:CH{t} {trVal}'.format(t=self.trigCh, trVal=self.trigVal))
		# Horizontal
		self.inst.write('HORizontal:DELay:MODe OFF')
		self.inst.write('HORizontal:POSition 10')  # Sets the horizontal position trigger to 10%
		self.inst.write('HORizontal:RECOrdlength {p}'.format(p=self.points))
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
		self.inst.write('CH{t}:LABel "Trigger"'.format(t=self.trigCh))
		self.inst.write('CH{s}:TERmination 50'.format(s=self.sigCh))
		if self.bias < 0:
			print 'Positive waves'
			self.inst.write('CH{s}:POSition -4'.format(s=self.sigCh))
		else:
			print 'Negative waves'
			self.inst.write('CH{s}:POSition 4'.format(s=self.sigCh))
		self.inst.write('CH{s}:COUPling DC'.format(s=self.sigCh))
		self.inst.write('CH{s}:SCAle 60E-3'.format(s=self.sigCh))  # 60E-3# 27.5E-3 #50E-3 #2E-3 for attenuated calibrations
		self.inst.write('CH{s}:YUNits "V"'.format(s=self.sigCh))
		self.inst.write('CH{s}:LABel "Signal"'.format(s=self.sigCh))
		# if self.bias >= 0:
		# 	self.inst.write('CH{s}:INVert ON'.format(s=self.sigCh))
		# 	print 'Inverting signals to get positive waves'
		# else:
		# 	self.inst.write('CH{s}:INVert OFF'.format(s=self.sigCh))
		# 	print 'Not inverting signals. They are already positive'
		# Aqcuisition
		self.inst.write('ACQuire:MODe HIRes') # SAMple
		self.inst.write('ACQuire:NUMEnv 1')
		if not self.doCont:
			self.inst.write('ACQuire:STOPAfter SEQuence')
		else:
			self.inst.write('ACQuire:STOPAfter RUNSTOP')
		self.inst.write('ACQuire:STATE ON')
		# Measurement
		self.inst.write('MEASUrement:CLEARSNapshot')
		self.inst.write('MEASUrement:GATing CURSor')
		self.inst.write('MEASUrement:IMMed:SOUrce CH{s}'.format(s=self.sigCh))
		self.inst.write('MEASUrement:IMMed:TYPe MAXimum')
		# Cursors
		self.inst.write('CURSor:FUNCtion WAVEform')
		self.inst.write('CURSor:SOUrce CH{s}'.format(s=self.sigCh))
		self.inst.write('CURSor:MODe TRACk')
		self.inst.write('CURSor:VBArs:UNIts SEConds')
		self.inst.write('CURSor:VBArs:POSITION1 {l}'.format(l=self.windowLow))
		self.inst.write('CURSor:VBArs:POSITION2 {h}'.format(h=self.windowHigh))


	def SetOutputFormat(self):
		self.inst.write('DATa INIT')
		self.inst.write('DATa:SOUrce CH{s}'.format(s=self.sigCh))
		self.inst.write('DATa:STARt 1')
		self.inst.write('DATa:STOP {p}'.format(p=self.points))
		self.inst.write('DATa:ENCdg SRPbinary') # RPBinary, SRIbinary,SRPbinary, FAStest, FPbinary
		self.inst.write('DATa:WIDTH 2')
		self.wfmo = str(self.inst.ask('WFMOutpre?'))
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

	def ReadData(self):
		self.inst.write('HEADer 0')
		if not self.doCont:
			self.WaitForSRQ()
		self.bindata = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.inst.write('HEADer 1')

	def ReadMeasurementFromDevice(self):
		self.inst.write('HEADer 0')
		if not self.doCont:
			self.WaitForSRQ()
		self.peak_val = float(self.inst.ask('MEASUrement:IMMed:VALue?'))
		self.inst.write('HEADer 1')

	def FormatData(self):
		self.volts = self.yzero + self.ymult * (self.bindata - self.yoffs)
		self.time = np.arange(self.xzero, self.xzero + self.xincr * (self.nrpt - 1), self.xincr)

	def GetMaximum(self):
		pos_window_low = int(round((self.windowLow-self.xzero)/self.xincr))
		pos_window_high = int(round((self.windowHigh-self.xzero)/self.xincr))
		# ipdb.set_trace(context=5)
		pos_max = np.argmax(self.volts[pos_window_low:pos_window_high])
		self.peak_val = float(self.volts[pos_window_low:pos_window_high][pos_max])
		self.peak_pos = float(self.time[pos_window_low:pos_window_high][pos_max])

	def ReadBothData(self):
		self.inst.write('HEADer 0')
		self.WaitForSRQ()
		self.bindata = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.peak_val2 = float(self.inst.ask('MEASUrement:IMMed:VALue?'))

	def ValidateMeasuredData(self):
		if self.doWaves:
			return True
		return self.peak_val > 0

	def DoMeasurement(self):
		if self.doWaves:
			self.ReadData()
			self.FormatData()
			self.GetMaximum()
		else:
			self.ReadMeasurementFromDevice()
		if self.verb:
			if self.doWaves:
				print 't = {t:.4}; v = {v:.4}'.format(t=self.peak_pos, v=self.peak_val)
			else:
				print 'v = {v:.4}'.format(v=self.peak_val)

	def WriteCSVRow(self, iteration):
		iterArray = np.full(len(self.volts), iteration)
		event = np.array((iterArray, self.time, self.volts), 'f').transpose()
		self.fileWavesWriter.writerows(event)

	def TakeMeasurements(self):
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
			valid = False
			while not valid:
				self.DoMeasurement()
				valid = self.ValidateMeasuredData()
			self.WriteCSVRow(iter)
			self.peak_values = np.append(self.peak_values, self.peak_val)
			self.peak_times = np.append(self.peak_times, self.peak_pos)
			bar.update(iter+1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1-t0), 's'

	def TakeBothMeasurements(self):
		if self.doCont:
			print 'Can\'t do in continuous mode. Run without -c'
			return
		t0 = time.time()
		for iter in xrange(self.meas):
			self.ReadBothData()
			self.FormatData()
			self.GetMaximum()
			self.peak_values_waves = np.append(self.peak_values_waves, self.peak_val)
			self.peak_values_measu = np.append(self.peak_values_measu, self.peak_val2)
		t1 = time.time()
		print 'Total time:', str(t1-t0), 's'

	def SetOutputFilesNames(self):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))

		def AddSuffix(string1):
			string1 += '_Pos' if self.bias > 0 else '_Neg'
			string1 += '_{b}V'.format(b=abs(int(self.bias)))
			string1 += '_Waves' if self.doWaves else '_Measure'
			string1 += '_Continuous' if self.doCont else '_SingleEvts'
			string1 += '.csv'
			return string1

		self.outString1 = '{dir}/Runs/data_{f}'.format(dir=self.outdir, f=self.filename)
		self.outString1 = AddSuffix(self.outString1)
		self.outString2 = '{dir}/Runs/point_{f}'.format(dir=self.outdir, f=self.filename)
		self.outString2 = AddSuffix(self.outString2)
		self.outString3 = '{dir}/Runs/waves_{f}'.format(dir=self.outdir, f=self.filename)
		self.outString3 = AddSuffix(self.outString3)

	def SaveMeasurements(self):
		np.savetxt(self.outString1, self.peak_values, delimiter=';')
		filep = open(self.outString2, 'w')
		filep.write('{m:.4}'.format(m=self.peak_values.mean(0, 'f')))
		filep.write('{s:.4}'.format(s=self.peak_values.std(0, 'f')))
		filep.close()

	def DoTriggerSweep(self):
		self.inst.write('SELect:CH{s} OFF'.format(s=self.sigCh))
		# Trigger
		# Aqcuisition
		self.inst.write('ACQuire:MODe HIRes')  # SAMple
		self.inst.write('ACQuire:NUMEnv 1')
		if not self.doCont:
			self.inst.write('ACQuire:STOPAfter SEQuence')
		else:
			self.inst.write('ACQuire:STOPAfter RUNSTOP')
		self.inst.write('ACQuire:STATE ON')
		# Measurement
		self.inst.write('MEASUrement:GATing OFF')
		# Cursors
		self.inst.write('CURSor:FUNCtion OFF')

		initialValue = float(raw_input('Enter the initial value to do the trigger sweep: '))
		finalValue = float(raw_input('Enter the final value to do the trigger sweep: '))
		stepValue = float(raw_input('Enter the step for the trigger sweep: '))
		timeCollect = float(raw_input('Enter the time in seconds for each measurement: '))

		trigThs = np.arange(initialValue, finalValue, stepValue)
		trigThs = trigThs if finalValue in trigThs else np.append(trigThs, finalValue)

		widgets = [
			progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			# ' ', progressbar.AdaptativeTransferSpeed()
		]
		bar = progressbar.ProgressBar(widgets=widgets, maxval=int(trigThs.size * timeCollect)+1)
		bar.start()
		numTrigs = np.zeros(trigThs.size, 'i')
		for iter in xrange(trigThs.size):
			self.inst.write('TRIGger:A:LEVel:CH{t} {trVal}'.format(t=self.trigCh, trVal=trigThs[iter]))
			t0 = time.time()
			t1 = time.time()
			while t1 - t0 < timeCollect:
				t1 = time.time()
				self.inst.write('HEADer 0')
				self.WaitForSRQ()
				numTrigs[iter] += 1
				self.inst.write('HEADer 1')
				t2 = time.time()
				if t2 - t1 > 5:
					t0 += (t2 - t1)
				t1 = t2
				bar.update(int(iter * timeCollect + t1 - t0))
		bar.finish()
		outf = open('{dir}/Runs/trigger_ctrl_{c}V_low_{l:4<}V_high_{h:4<}V.csv'.format(dir=self.outdir, c=self.trigVal,
																				 l=abs(initialValue), h=abs(finalValue))
					,'wb')
		trigWriter = csv.writer(outf, delimiter=',')
		outArray = np.array((trigThs, numTrigs)).T
		trigWriter.writerows(outArray)
		print 'The results are:'
		for i in xrange(trigThs.size):
			print '{t}\t: {num}'.format(t=trigThs[i], num=numTrigs[i])
		outf.close()

	def Quit(self):
		self.inst.before_close()
		self.inst.close()
		# sys.exit(0)
		exit()

	def TakeRandomTriggers(self):
		widgets = [
			progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			# ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.inst.write('TRIGger:A:MODe AUTO')
		bar = progressbar.ProgressBar(widgets=widgets, maxval=self.meas)
		bar.start()
		t0 = time.time()
		for iter in xrange(self.meas):
			valid = False
			while not valid:
				self.DoRandomMeasurement()
				valid = self.ValidateMeasuredData()
			self.WriteCSVRow(iter)
			self.peak_values = np.append(self.peak_values, self.peak_val)
			self.peak_times = np.append(self.peak_times, self.peak_pos)
			bar.update(iter + 1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1 - t0), 's'

	def DoRandomMeasurement(self):
		if self.doWaves:
			self.RandomReadData()
			self.FormatData()
			self.GetMaximum()
		else:
			self.ReadMeasurementFromDevice()
		if self.verb:
			if self.doWaves:
				print 't = {t:.4}; v = {v:.4}'.format(t=self.peak_pos, v=self.peak_val)
			else:
				print 'v = {v:.4}'.format(v=self.peak_val)

	def RandomReadData(self):
		self.inst.write('HEADer 0')
		if not self.doCont:
			self.WaitForRandomSRQ()
		self.bindata = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.inst.write('HEADer 1')

	def WaitForRandomSRQ(self):
		self.inst.write('ACQuire:STATE OFF')
		self.inst.write('ACQuire:STOPAfter RUnsTOP')
		# self.inst.write('DESE 1')  # enable device event status enable register
		# self.inst.write('*ESE 1')  # enable event status enable register
		# self.inst.write('*SRE 32')  # Generate a service request from the oscilloscope
		self.inst.write('ACQuire:STATE ON')
		# self.inst.write('*OPC')  # wait for interruption
		time.sleep(abs(np.random.normal(0.5, 0.2)))
		# while True:
		# 	ev = int(self.inst.ask('*ESR?'))  # event in queue
		# 	if int(self.inst.ask('EVENT?')) == 402:  # message code for completion
		# 		break
		return

	# TODO ONLY ONE WAVE FORM AT A TIME. DO !!!!!!

	def SetOutputFormatThree(self):
		self.points = 10000
		self.inst.write('DLC')
		self.inst.write('*CLS')
		# Trigger
		self.inst.write('SELect:CH1 ON')
		self.inst.write('SELect:CH2 ON')
		self.inst.write('SELect:CH3 ON')
		self.inst.write('TRIGger:A:TYPe EDG')
		self.inst.write('TRIGger:A:MODe NORMal')
		self.inst.write('TRIGger:A:EDGE:COUP DC')
		self.inst.write('TRIGger:A:EDGE:SLOpe FALL')
		self.inst.write('TRIGger:A:EDGE:SOUrce CH{t}'.format(t=self.trigCh))
		self.inst.write('TRIGger:A:LEVel:CH1 {trVal}'.format(t=self.trigCh, trVal=self.trigVal))
		# Horizontal
		self.inst.write('HORizontal:DELay:MODe OFF')
		self.inst.write('HORizontal:POSition 10')  # Sets the horizontal position trigger to 10%
		self.inst.write('HORizontal:RECOrdlength {p}'.format(p=self.points))
		self.inst.write('HORizontal:SCAle 1E-6')
		# Vertical
		self.inst.write('CH1:TERmination 50')
		self.inst.write('CH2:TERmination 50')
		self.inst.write('CH1:POSition 4.5')
		self.inst.write('CH2:POSition 4.5')
		self.inst.write('CH1:COUPling DC')
		self.inst.write('CH2:COUPling DC')
		self.inst.write('CH1:SCAle 20E-3') # 50E-3 # 10E-3 for attenuated calibrations
		self.inst.write('CH2:SCAle 10E-3') # 50E-3 # 10E-3 for attenuated calibrations
		self.inst.write('CH1:YUNits "V"')
		self.inst.write('CH2:YUNits "V"')
		self.inst.write('CH1:LABel "Trigger1"')
		self.inst.write('CH2:LABel "Trigger2"')
		self.inst.write('CH3:TERmination 50')
		if self.bias < 0:
			print 'Positive waves'
			self.inst.write('CH3:POSition -4')
		else:
			print 'Negative waves'
			self.inst.write('CH3:POSition 4')
		self.inst.write('CH3:COUPling DC')
		self.inst.write('CH3:SCAle 60E-3')  # 60E-3# 27.5E-3 #50E-3 #2E-3 for attenuated calibrations
		self.inst.write('CH3:YUNits "V"')
		self.inst.write('CH3:LABel "Signal"')
		# if self.bias >= 0:
		# 	self.inst.write('CH3:INVert ON')
		# 	print 'Inverting signals to get positive waves'
		# else:
		# 	self.inst.write('CH3:INVert OFF'.format(s=self.sigCh))
		# 	print 'Not inverting signals. They are already positive'
		# Aqcuisition
		self.inst.write('ACQuire:MODe HIRes') # SAMple
		self.inst.write('ACQuire:NUMEnv 1')
		if not self.doCont:
			self.inst.write('ACQuire:STOPAfter SEQuence')
		else:
			self.inst.write('ACQuire:STOPAfter RUNSTOP')
		self.inst.write('ACQuire:STATE ON')
		self.nrpt = np.zeros(3, 'I')
		self.xzero = np.zeros(3, 'f8')
		self.xincr = np.zeros(3, 'd')
		self.xunit = np.empty(3, 'U')
		self.yzero = np.zeros(3, 'f8')
		self.ymult = np.zeros(3, 'f8')
		self.yoffs = np.zeros(3, 'f8')
		self.yunit = np.empty(3, 'U')

	def SetupOutputChannel(self, ch=3):
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

	def TakeThreeWaves(self):
		self.SetOutputFormatThree()
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
			valid = False
			while not valid:
				self.DoMeasurementThree()
				valid = self.ValidateMeasuredData()
			self.WriteCSVRowThree(iter)
			self.peak_values = np.append(self.peak_values, self.peak_val)
			self.peak_times = np.append(self.peak_times, self.peak_pos)
			bar.update(iter+1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1-t0), 's'

	def DoMeasurementThree(self):
		self.ReadThreeData()
		self.FormatThreeData()

	def ReadThreeData(self):
		self.inst.write('HEADer 0')
		if not self.doCont:
			self.WaitForSRQ()
		self.SetupOutputChannel(1)
		self.bindata1 = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.SetupOutputChannel(2)
		self.bindata2 = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.SetupOutputChannel(3)
		self.bindata3 = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
		self.inst.write('HEADer 1')
	
	def FormatThreeData(self):
		self.volts1 = self.yzero[0] + self.ymult[0] * (self.bindata1 - self.yoffs[0])
		self.volts2 = self.yzero[1] + self.ymult[1] * (self.bindata2 - self.yoffs[1])
		self.volts3 = self.yzero[2] + self.ymult[2] * (self.bindata3 - self.yoffs[2])
		self.time1 = np.linspace(self.xzero[0], self.xzero[0] + self.xincr[0] * (self.nrpt - 1), self.nrpt)
		self.time2 = np.linspace(self.xzero[1], self.xzero[1] + self.xincr[1] * (self.nrpt - 1), self.nrpt)
		self.time3 = np.linspace(self.xzero[2], self.xzero[2] + self.xincr[2] * (self.nrpt - 1), self.nrpt)

	def WriteCSVRowThree(self, iteration):
		iterArray = np.full(len(self.volts1), iteration)
		# ipdb.set_trace()
		event = np.array((iterArray, self.time1, self.volts1, self.volts2, self.volts3), 'f').transpose()
		self.fileWavesWriter.writerows(event)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-i', '--ip', dest='ip', default='169.254.2.226', type='string', help='IP address of the instrument. e.g. 169.254.2.226')
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string', help='Relative path to output directory where the Runs folder is e.g. ./')
	parser.add_option('-f', '--filename', dest='filename', default='', type='string', help='Filename to save the data. e.g. output_2017_06')
	parser.add_option('-s', '--signal', dest='sig', default=2, type='int', help='The channel connected to the signal. e.g. 2')
	parser.add_option('-t', '--trigger', dest='trig', default=1, type='int', help='The channel connected to the trigger (e.g. scintillator) e.g. 1')
	parser.add_option('-l', '--trigval', dest='trigval', default=-0.08, type='float', help='The value of the trigger e.g. -0.04')
	parser.add_option('-m', '--measurements', dest='meas', default=100, type='int', help='The number of measurements to make e.g. 100')
	parser.add_option('-w', '--waveforms', dest='waves', default=False, action='store_true', help='Toggles full waveforms calculation instead of measurement acquisition')
	parser.add_option('-c', '--continuous', dest='cont', default=False, action='store_true', help='Toggles continuous data sampling instead of successive single measurements')
	parser.add_option('-p', '--pulser', dest='pulser', default=False, action='store_true', help='Toggles trigger from pulser')
	parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true', help='Toggles automatic analysis')
	parser.add_option('-b', '--bias', dest='bias', default=-400, type='int', help='The HV bias used on the sample e.g. -400')
	parser.add_option('-g', '--trigsweep', dest='trigsweep', default=False, action='store_true', help='Toggles trigger threshold analysis')
	parser.add_option('-r', '--randomtrig', dest='ran', default=False, help='Toggles random Triggers taken each 0.5sec with a std of 0.2 sec', action='store_true')
	parser.add_option('-d', '--doubletrig', dest='dtrig', default=False, help='Toggles double trigger analysis', action='store_true')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	ip = str(options.ip)
	outdir = str(options.outdir)
	filename = str(options.filename)
	sig = int(options.sig)
	trig = int(options.trig)
	trigval = float(options.trigval)
	meas = int(options.meas)
	waves = bool(options.waves)
	cont = bool(options.cont)
	pulser = bool(options.pulser)
	automatic = bool(options.automatic)
	bias = int(options.bias)
	trigsweep = bool(options.trigsweep)
	ran = bool(options.ran)
	dtrig = bool(options.dtrig)
	verb = bool(options.verb)
	if dtrig:
		filename += '_ThreeWaves'
	z = TektronixReadout(ip, outdir, filename, sig, trig, trigval, meas, waves, cont, pulser, bias, verb)
	if trigsweep:
		z.DoTriggerSweep()
	elif automatic:
		z.SetOutputFilesNames()
		z.fileWaves = open(z.outString3, 'wb')
		z.fileWavesWriter = csv.writer(z.fileWaves, delimiter=',')
		if ran:
			z.TakeRandomTriggers()
		elif dtrig:
			z.TakeThreeWaves()
		else:
			z.TakeMeasurements()
		print 'Results:', str(z.peak_values.mean(0, 'f')), '+/-', str(z.peak_values.std(0, 'f')), 'V. For ', str(z.meas), 'measurements :)'
		z.fileWaves.close()
		z.SaveMeasurements()
		z.inst.before_close()
		z.inst.close()
		sys.stdout.write('\a\a\a')
		sys.stdout.flush()
