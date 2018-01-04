#!/usr/bin/env python
import csv
import numpy as np
from struct import unpack
import time, os
from optparse import OptionParser
import progressbar
from ROOT import TProfile2D, TTree, gStyle, TFile, vector
import ipdb
# from TektronixReadout import TektronixReadout

class WaveAnalysis:
	def __init__(self, outputDir, inputDir, fileName, bias, verbose):
		self.outDir, self.indir, self.fileName, self.bias, self.verb = outputDir, inputDir, fileName, bias, verbose
		self.tree = TTree('tree', 'tree')
		t0 = time.time()
		self.tree.ReadFile('{ind}/{f}'.format(ind=self.indir, f=self.fileName), 'event:time:voltage')
		t0 = time.time() - t0
		print 'Time reading file:', t0, 'seconds'
		self.entries = self.tree.GetEntries()
		self.tree.GetEntry(self.entries - 1)
		lastEvent = self.tree.event
		self.tree.GetEntry(0)
		firstEvent = self.tree.event
		self.events = lastEvent - firstEvent + 1
		self.ptsWave = int(self.entries/(self.events))
		self.tree.SetEstimate(self.entries)
		leng = self.tree.Draw('time:voltage', '', 'goff')
		tempTime = self.tree.GetVal(0)
		tempVolt = self.tree.GetVal(1)
		arrayTime = [[tempTime[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries/self.ptsWave))]
		arrayVolt = [[tempVolt[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries/self.ptsWave))]
		self.timeVect = np.array(arrayTime, 'f4')
		self.voltVect = np.array(arrayVolt, 'f4')
		self.timeBra = np.zeros(self.ptsWave, 'f4')
		self.voltBra = np.zeros(self.ptsWave, 'f4')
		self.event = np.zeros(1, 'I')
		rawtreename = 'raw_tree_pos_{b}V'.format(b=self.bias) if self.bias >= 0 else 'raw_tree_neg_{b}V'.format(b=self.bias)
		self.fileRaw = TFile('{r}.root'.format(r=rawtreename), 'RECREATE')
		self.treeRaw = TTree('{n}'.format(n=rawtreename), '{n}'.format(n=rawtreename))
		self.CreateRawRootFile()
		self.signalTLow, self.signalTHigh = 0.5e-6, 5e-6
		self.pedestalTPos = 0

	def CreateRawRootFile(self):
		self.treeRaw.Branch('event', self.event, 'event/i')
		self.treeRaw.Branch('time', self.timeBra, 'time[{s}]/F'.format(s=self.ptsWave))
		self.treeRaw.Branch('voltage', self.voltBra, 'voltage[{s}]/F'.format(s=self.ptsWave))
		t1 = time.time()
		for i in xrange(int(self.events)):
			self.event.fill(i)
			np.putmask(self.timeBra, 1 - np.zeros(self.ptsWave, '?'), self.timeVect[i])
			np.putmask(self.voltBra, 1 - np.zeros(self.ptsWave, '?'), self.voltVect[i])
			numFil = self.treeRaw.Fill()
		t2 = time.time() - t1
		self.fileRaw.Write()
		self.fileRaw.Close()
		print 'Time creating raw Tree', t2, 'seconds'

	def DefineSignalWindow(self):
		doChange = raw_input('Type 1 to change signal window; 0 to leave it as it is (current: [{l}, {h}])? '.format(l=self.signalTLow, h=self.signalTHigh))
		try:
			doChange = bool(int(doChange))
		except Exception:
			print 'The value entered cannot be understood. Leaving signal window as it is'
			doChange = False
		if doChange:
			cont = False
			while not cont:
				low = raw_input('Enter the time in seconds after the trigger for the lower time window (should be positive): ')
				try:
					low = float(low)
				except Exception:
					print 'The value is not correct. Try again'
				if low < 0:
					print 'The value is negative. Try again'
				elif low >= self.timeVect[0, -1]:
					print 'The value should be less than the acquisition maximum ({tm}). Try again'.format(tm=self.timeVect[0, -1])
				else:
					cont = True
					self.signalTLow = low
			cont = False
			while not cont:
				high = raw_input('Enter the time in seconds after the trigger for the higher time window (should be positive): ')
				try:
					high = float(high)
				except Exception:
					print 'The value is not correct. Try again'
				if high < self.signalTLow:
					print 'The value should be greater than the lower time window {l}'.format(l=self.signalTLow)
				elif high > self.timeVect[0, -1]:
					print 'The value should be at most the acquisition maximum ({tm}). Try again'.format(tm=self.timeVect[0, -1])
				else:
					cont = True
					self.signalTHigh = high
			print 'The signal window is now: [{l}, {h}] seconds'.format(l=self.signalTLow, h=self.signalTHigh)

	def Quit(self):
		# sys.exit(0)
		exit()


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string',
					  help='Relative path to output directory e.g. ./')
	parser.add_option('-i', '--indir', dest='indir', default='./', type='string',
					  help='Relative path to input directory e.g. ./')
	parser.add_option('-f', '--filename', dest='filename', default='', type='string',
					  help='Filename to save the data. e.g. output_2017_06')
	parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true',
					  help='Toggles automatic analysis')
	parser.add_option('-b', '--bias', dest='bias', default=-400, type='int',
					  help='The HV bias used on the sample e.g. -400')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	outdir = str(options.outdir)
	indir = str(options.indir)
	filename = str(options.filename)
	automatic = bool(options.automatic)
	bias = int(options.bias)
	verb = bool(options.verb)
	z = WaveAnalysis(outdir, indir, filename, bias, verb)
	if automatic:
		pass