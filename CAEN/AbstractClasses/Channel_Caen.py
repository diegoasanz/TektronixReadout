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


class Channel_Caen:
	def __init__(self, ch, ch_type, verb=False):
		self.ch = ch
		self.type = ch_type
		self.verb = verb

		self.base_line_u_adcs = 0
		self.dc_offset_percent = 0
		self.thr_counts = 0
		self.edge = -1
		self.name = str(ch_type)

	def Set_Channel(self, settings):
		if self.type == 'signal':
			self.edge = -int(settings.bias/abs(settings.bias))
			self.Calculate_DC_Offset_Percentage(settings)
		elif self.type == 'trigger':
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.trig_base_line, settings.sigRes)
			self.thr_counts = settings.trig_thr_counts
			self.Calculate_DC_Offset_Percentage(settings)
			self.edge = -1
		elif self.type == 'veto':
			self.base_line_u_adcs = self.Calculate_Universal_ADCs(settings.ac_base_line, settings.sigRes)
			self.thr_counts = settings.ac_thr_counts
			self.Calculate_DC_Offset_Percentage(settings)
			self.edge = -1

	def Calculate_DC_Offset_Percentage(self, settings):
		if self.type == 'signal':
			self.dc_offset_percent = 45 if settings.bias < 0 else -45
		else:
			self.dc_offset_percent = round(100 * np.divide(3 * self.thr_counts + self.base_line_u_adcs, 2.0**settings.dig_bits - 1.0, dtype='f8') - 50)

	def Calculate_Universal_ADCs(self, value_volts, sig_res):
		return np.divide(value_volts, sig_res, dtype='f8')

	def Correct_Base_Line(self, mean_volts, settings):
		self.base_line_u_adcs = self.Calculate_Universal_ADCs(mean_volts, settings.sigRes)
		self.Calculate_DC_Offset_Percentage(settings)

	def ADC_to_Volts(self, adcs, sigres):
		return np.multiply(sigres, np.add(adcs, np.multiply(2**14 - 1, self.dc_offset_percent/100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')

if __name__ == '__main__':
	print 'bla'