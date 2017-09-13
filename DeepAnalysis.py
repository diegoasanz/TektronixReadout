#!/usr/bin/env python
import visa
import numpy as np
from struct import unpack
import time, os
from optparse import OptionParser
import progressbar
# from TektronixReadout import TektronixReadout

class DeepAnalysis:
	def __init__(self, tektronixReadout=None):
		self.tektRO = tektronixReadout
		self.verb = self.tektRO.verb
		self.inst = None if self.tektRO is None else self.tektRO.inst
		self.doCont = None if self.tektRO is None else self.tektRO.doCont
		self.doWaves = None if self.tektRO is None else self.tektRO.doWaves
		self.meas = 0 if self.tektRO is None else int(self.tektRO.meas)
		self.nrpt = 0 if self.tektRO is None else self.tektRO.nrpt
		self.xzero = 0 if self.tektRO is None else self.tektRO.xzero
		self.xincr = 0 if self.tektRO is None else self.tektRO.xincr
		self.xunit = '' if self.tektRO is None else self.tektRO.xunit
		self.yzero = 0 if self.tektRO is None else self.tektRO.yzero
		self.ymult = 0 if self.tektRO is None else self.tektRO.ymult
		self.yoffs = 0 if self.tektRO is None else self.tektRO.yoffs
		self.yunit = '' if self.tektRO is None else self.tektRO.yunit
		self.ped = None if self.tektRO is None else self.tektRO.ped
		self.prev_meas = None
		self.ped_points = 90

		self.bindata = None
		self.volts, self.time, self.peak_val, self.peak_pos = None, None, None, None
		self.peak_values = np.empty(0, 'f')
		self.signals = np.empty(0, 'f')
		self.peak_times = np.empty(0, 'f')
		self.ped_values = np.empty(0, 'f')
		self.ped_th = 1

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

	def ValidateMeasuredData(self):
		if self.ped is not None:
			return self.peak_val > self.ped['sigma'] * self.ped_th
		else:
			if self.prev_meas is None:
				self.prev_meas = self.volts
			else:
				if not (self.prev_meas == self.volts).all():
					self.prev_meas = self.volts
					return True
				else:
					return False

	def DoMeasurement(self):
		if self.doWaves:
			self.ReadData()
			self.FormatData()
			if self.ped is not None:
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
			if self.ped is not None:
				self.peak_values = np.append(self.peak_values, self.peak_val)
				self.peak_times = np.append(self.peak_times, self.peak_pos)
			else:
				self.ped_values = np.append(self.ped_values, self.volts)
			bar.update(iter+1)
		t1 = time.time()
		bar.finish()
		print 'Total time:', str(t1-t0), 's'

if __name__ == '__main__':
	z = DeepAnalysis(None)
