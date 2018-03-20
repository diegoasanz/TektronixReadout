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
from copy import deepcopy
import glob


# from DataAcquisition import DataAcquisition

class Utils:
	def __init__(self):
		self.bar = None

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

def FindRedundantEvents(settings):
	channels = [settings.sigCh, settings.trigCh, settings.acCh] if settings.ac_enable else [settings.sigCh, settings.trigCh]
	types = {settings.sigCh: 'signal', settings.trigCh: 'trigger', settings.acCh: 'veto'} if settings.ac_enable else {settings.sigCh: 'signal', settings.trigCh: 'trigger'}
	for ch in channels:
		print '\n', types[ch], ':\n'
		filename = '{d}/Runs/{f}_{t}.dat'.format(d=settings.outdir, f=settings.filename, t=types[ch])
		f0 = open(filename, 'rb')
		events = int(np.floor(os.path.getsize(filename) / settings.struct_len))
		for ev0 in xrange(events):
			f0.seek(ev0 * settings.struct_len)
			dataev0 = f0.read(settings.struct_len)
			dataev0 = struct.Struct(settings.struct_fmt).unpack_from(dataev0)
			dataev0 = np.array(dataev0, 'H')
			for ev1 in xrange(ev0+1, events):
				f0.seek(ev1 * settings.struct_len)
				dataev1 = f0.read(settings.struct_len)
				dataev1 = struct.Struct(settings.struct_fmt).unpack_from(dataev1)
				dataev1 = np.array(dataev1, 'H')

				if (dataev0 == dataev1).all():
					print '', ev0, '\t', ev1
		f0.close()

def IsInt(i):
	try:
		int(i)
		return True
	except ValueError:
		return False

def IsFloat(f):
	try:
		float(f)
		return True
	except ValueError:
		return False

def PlotHisto2DLimits(xmin=-0.48e-6, xmax=4.68e-6, ymin=-0.5, ymax=0.1, xres=2e-9, yres=2.15/(2.0**14-1)):
	xbins = int(round((xmax - xmin)/xres))
	ybins = int(round((ymax - ymin)/yres))
	print '({nx},{minx},{maxx},{ny},{miny},{maxy})'.format(nx=xbins, minx=xmin, maxx=xmax, ny=ybins, miny=ymin,maxy=ymax)

if __name__ == '__main__':
	print 'blaaaa'


