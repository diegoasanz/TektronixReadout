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
		self.sigma_adcs = 10
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
			limit = self.base_line_u_adcs + 10 * self.sigma_adcs
			if limit < 0:
				self.dc_offset_percent = -50
			else:
				self.dc_offset_percent = int(round(100 * (limit/float(2**settings.dig_bits - 1) - 0.5)))

	def Calculate_Universal_ADCs(self, value_volts, sig_res):
		return np.divide(value_volts, sig_res, dtype='f8')

	# def Correct_Base_Line(self, mean_volts, sigma_counts, settings):
	# 	self.base_line_u_adcs = self.Calculate_Universal_ADCs(mean_volts, settings.sigRes)
	# 	self.sigma_adcs = sigma_counts
	# 	self.Calculate_DC_Offset_Percentage(settings)

	def Correct_Base_Line2(self, mean_adc, sigma_adc, settings):
		self.sigma_adcs = sigma_adc
		variable = mean_adc + 50 * sigma_adc
		self.base_line_u_adcs = self.Calculate_Universal_ADCs(self.ADC_to_Volts(mean_adc, settings.sigRes, settings.dig_bits), settings.sigRes)
		if variable > (2**settings.dig_bits - 1) * (0.5 + self.dc_offset_percent / 100.0):
			self.dc_offset_percent += int(round(100.0 * (variable + 1 - 2.0**settings.dig_bits) / (2.0**settings.dig_bits - 1.0)))
		else:
			self.dc_offset_percent = -50

	def Correct_Threshold(self, sigma):
		if self.type == 'trigger':
			self.thr_counts = int(round(max(10*self.sigma_adcs, self.thr_counts)))
		elif self.type == 'veto':
			self.thr_counts = int(round(max(4*self.sigma_adcs, self.thr_counts)))

	def ADC_to_Volts(self, adcs, sigres, nbits=14):
		return np.multiply(sigres, np.add(adcs, np.multiply(2**nbits - 1, self.dc_offset_percent/100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')

if __name__ == '__main__':
	print 'bla'