#!/usr/bin/env python
import visa
import numpy as np
from struct import unpack
import time, os
from optparse import OptionParser
import progressbar
# from TektronixReadout import TektronixReadout

class DataAcquisition:
	def __init__(self, tektronixReadout=None, verbose=False):
		self.verb = verbose
		self.tektronixRO = tektronixReadout
		self.inst = None if self.tektronixRO == None else self.tektronixRO.inst
		self.doCont = None if self.tektronixRO == None else self.tektronixRO.doCont
		self.doWaves = None if self.tektronixRO == None else self.tektronixRO.doWaves

		self.bindata = None
		self.volts, self.time, self.peak_val, self.peak_pos = None, None, None, None
		self.peak_values = np.empty(0, 'f')
		self.peak_times = np.empty(0, 'f')
		# self.peak_val2 = None
		# self.peak_values_waves = np.empty(0,'f')
		# self.peak_values_measu = np.empty(0,'f')

	def SetInstrument(self, ipn):
		self.tektronixRO = TektronixReadout(verb=self.verb)
		self.inst = self.tektronixRO.inst
		self.doCont = self.tektronixRO.doCont
		self.doWaves = self.tektronixRO.doWaves

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
		pos_max = np.argmax(self.volts)
		self.peak_val = float(self.volts[pos_max])
		self.peak_pos = float(self.time[pos_max])

	# def ReadBothData(self):
	# 	self.inst.write('HEADer 0')
	# 	self.WaitForSRQ()
	# 	self.bindata = self.inst.query_binary_values('CURVE?', datatype='H', is_big_endian=False, container=np.array)
	# 	self.peak_val2 = float(self.inst.ask('MEASUrement:IMMed:VALue?'))

	def ValidateMeasuredData(self):
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
			self.peak_values = np.append(self.peak_values, self.peak_val)
			self.peak_times = np.append(self.peak_times, self.peak_pos)
			bar.update(iter+1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1-t0), 's'

	# def TakeBothMeasurements(self):
	# 	if self.doCont:
	# 		print 'Can\'t do in continuous mode. Run without -c'
	# 		return
	# 	t0 = time.time()
	# 	for iter in xrange(self.meas):
	# 		self.ReadBothData()
	# 		self.FormatData()
	# 		self.GetMaximum()
	# 		self.peak_values_waves = np.append(self.peak_values_waves, self.peak_val)
	# 		self.peak_values_measu = np.append(self.peak_values_measu, self.peak_val2)
	# 	t1 = time.time()
	# 	print 'Total time:', str(t1-t0), 's'

	def SaveMeasurements(self):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))
		string1 = '{dir}/Runs/data_{f}'.format(dir=self.outdir, f=self.filename)
		string1 += '_Waves' if self.doWaves else '_Measure'
		string1 += '_Continuous' if self.doCont else '_SingleEvts'
		string1 += '.csv'
		np.savetxt(string1, self.peak_values, delimiter=';')
		string2 = '{dir}/Runs/point_{f}'.format(dir=self.outdir, f=self.filename)
		string2 += '_Waves' if self.doWaves else '_Measure'
		string2 += '_Continuous' if self.doCont else '_SingleEvts'
		string2 += '.csv'
		filep = open(string2, 'w')
		filep.write('{m:.4}'.format(m=self.peak_values.mean(0, 'f')))
		filep.write('{s:.4}'.format(s=self.peak_values.std(0, 'f')))
		filep.close()

	def Quit(self):
		self.inst.before_close()
		self.inst.close()
		# sys.exit(0)
		exit()


if __name__ == '__main__':
	z = DataAcquisition(None, True)
