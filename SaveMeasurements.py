#!/usr/bin/env python
import visa
import numpy as np
from struct import unpack
import time, os
from optparse import OptionParser
import progressbar
# from TektronixReadout import TektronixReadout

class SaveMeasurements:
	def __init__(self, tektronixReadout=None):
		self.tektRO = tektronixReadout
		self.verb = False if self.tektRO is None else self.tektRO.verb
		self.outdir = None if self.tektR0 is None else self.tektRO.outdir
		self.filename = None if self.tektRO is None else self.tektRO.filename
		self.doCont = None if self.tektRO is None else self.tektRO.doCont
		self.doWaves = None if self.tektRO is None else self.tektRO.doWaves
		self.bindata = None
		self.volts, self.time, self.peak_val, self.peak_pos = None, None, None, None
		self.peak_values = np.empty(0, 'f')
		self.peak_times = np.empty(0, 'f')
		# self.peak_val2 = None
		# self.peak_values_waves = np.empty(0,'f')
		# self.peak_values_measu = np.empty(0,'f')

	def SaveMeasurements(self, values):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))
		string1 = '{dir}/Runs/data_{f}'.format(dir=self.outdir, f=self.filename)
		string1 += '_Waves' if self.doWaves else '_Measure'
		string1 += '_Continuous' if self.doCont else '_SingleEvts'
		string1 += '.csv'
		np.savetxt(string1, values, delimiter=';')
		string2 = '{dir}/Runs/point_{f}'.format(dir=self.outdir, f=self.filename)
		string2 += '_Waves' if self.doWaves else '_Measure'
		string2 += '_Continuous' if self.doCont else '_SingleEvts'
		string2 += '.csv'
		filep = open(string2, 'w')
		filep.write('{m:.4}'.format(m=values.mean(axis=0, dtype='f')))
		filep.write('{s:.4}'.format(s=values.std(axis=0, dtype='f', ddof=1)))
		filep.close()

if __name__ == '__main__':
	z = SaveMeasurements(None)
