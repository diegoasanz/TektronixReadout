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
from Utils import *
from collections import OrderedDict
# from DataAcquisition import DataAcquisition


class CCD_Analysis:
	def __init__(self, infile='', outdir='', bias=0, hasCal=True, verbose=False):
		self.outDir, self.inputFile, self.bias, self.hasCal, self.verb = outdir, infile, bias, hasCal, verbose
		self.treeName = self.inputFile.split('.')[0]
		self.fileRaw, self.treeRaw = None, None
		self.ptsWave, self.event, self.events = 0, np.zeros(1, 'I'), 0
		self.eventVect = np.empty(0, 'f8')
		self.timeVect, self.signalWaveVect, self.triggerWaveVect, self.vetoWaveVect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')

		self.ped, self.pedSigma, self.sigAndPed, self.sigAndPedSigma = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.vetoedEvent, self.badShape, self.badPedestal = np.empty(0, '?'), np.empty(0, np.dtype('int8')), np.empty(0, '?')
		self.voltageHV, self.currentHV = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.timeHV = np.empty(0, 'f8')

		self.signalWaveMeanVect, self.signalWaveSigmaVect = None, None

		self.pedestalIntegrationTime = 0.4e-6
		self.pedestalTEndPos = -20e-9
		self.peakTime = 2.126e-6 if self.bias >= 0 else 2.120e-6
		self.peakForward = self.pedestalIntegrationTime/2.0
		self.peakBackward = self.pedestalIntegrationTime/2.0

		self.branches1D = ['event', 'vetoedEvent', 'badShape', 'badPedestal', 'voltageHV', 'currentHV', 'timeHV']
		self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV': 'float64'}
		self.branchesWaves = ['time', 'voltageSignal', 'voltageTrigger', 'voltageVeto']
		self.branchesWavesType = {'time': 'float64', 'voltageSignal': 'float64', 'voltageTrigger': 'float64', 'voltageVeto': 'float64'}
		self.branchesAll = self.branches1D + self.branchesWaves

		self.dicBraVect1D = OrderedDict()
		self.dicBraVectWaves = OrderedDict()

		self.hasBranch = {}

	def OpenROOTFile(self, mode='READ'):
		if not os.path.isdir(self.outDir):
			print 'Directory:', self.outDir, '; does not exist. Exiting!'
			exit()
		if not os.path.isfile('{o}/{f}'.format(o=self.outDir, f=self.fileRaw)):
			print 'File:', self.fileRaw, '; does not exist in:', self.outDir, '. Exiting!'
			exit()
		if self.fileRaw:
			if self.fileRaw.IsOpen():
				if self.fileRaw.GetOption().lower() != mode.lower():
					if self.fileRaw.ReOpen(mode) == -1:
						print 'Could not reopen file:', self.fileRaw, 'in mode:', mode, '. Exiting!'
						exit()
				return
			else:
				self.fileRaw = None
		self.fileRaw = ro.TFile('{o}/{f}'.format(o=self.outDir, f=self.fileRaw), mode)

	def LoadTree(self):
		if not self.fileRaw.IsOpen():
			print 'First OpenROOTFile before Loading Tree'
			return
		self.treeRaw = self.fileRaw.Get(self.treeName)
		self.hasBranch = {branch: self.TreeHasBranch(branch) for branch in self.branchesAll}
		self.UpdateBranchesLists()
		self.IsTimeHVTimeStamp()

	def UpdateBranchesLists(self):
		for branch in self.branches1D[:]:
			if not self.hasBranch[branch]:
				self.branches1D.remove(branch)
		for branch in self.branchesWaves[:]:
			if not self.hasBranch[branch]:
				self.branchesWaves.remove(branch)

	def IsTimeHVTimeStamp(self):
		if self.hasBranch['timeHV']:
			if self.treeRaw.GetLeaf('timeHV').GetTypeName() != 'TDatime':
				self.branches1D = ['timeHV.AsDouble()' if branch == 'timeHV' else branch for branch in self.branches1D]
				self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.AsDouble()': 'float64'}
			else:
				self.branches1D = ['timeHV.Convert()' if branch == 'timeHV' else branch for branch in self.branches1D]
				self.branches1DType = {'event': 'uint32', 'vetoedEvent': 'bool', 'badShape': 'int8', 'badPedestal': 'bool', 'voltageHV': 'float32', 'currentHV': 'float32', 'timeHV.Convert()': 'uint32'}

	def TreeHasBranch(self, branch):
		if self.treeRaw.GetLeaf(branch):
			return True
		else:
			return False

	def LoadVectorsFromTree(self):
		self.ptsWave = self.treeRaw.GetLeaf('time').GetLen()
		self.events = self.treeRaw.GetEntries()

		leng = self.treeRaw.Draw(':'.join(self.branches1D), '', 'goff para', self.events, 0)
		if leng == -1:
			print 'Error, could not load the branches: {b}. Try again :('.format(b=':'.join(self.branches1D))
			return
		while leng >= self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(leng)
			leng = self.treeRaw.Draw(':'.join(self.branches1D), '', 'goff para', self.events, 0)
		self.events = leng
		for pos, branch in enumerate(self.branches1D):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = self.treeRaw.GetVal(pos)
			self.dicBraVect1D[branch] = np.array([temp[ev] for ev in xrange(self.events)], dtype=np.dtype(self.branches1DType[branch]))
			if self.verb: print 'Done'

		leng = self.treeRaw.Draw(':'.join(self.branchesWaves), '', 'goff para', self.events, 0)
		if leng == -1:
			print 'Error, could not load the branches {b}. Try again :('.format(b=':'.join(self.branchesWaves))
			return
		while leng >= self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(leng)
			leng = self.treeRaw.Draw(':'.join(self.branchesWaves), '', 'goff para', self.events, 0)
		for pos, branch in enumerate(self.branchesWaves):
			if self.verb: print 'Vectorising branch:', branch, '...', ; sys.stdout.flush()
			temp = self.treeRaw.GetVal(pos)
			self.dicBraVectWaves[branch] = np.array([[temp[ev * self.ptsWave + pt] for pt in xrange(self.ptsWave)] for ev in xrange(self.events)], dtype=np.dtype(self.branchesWavesType[branch]))
			if self.verb: print 'Done'

	def ExplicitVectorsFromDictionary(self):
		if self.hasBranch['voltageSignal']:
			self.signalWaveVect = self.dicBraVectWaves['voltageSignal']
		if self.hasBranch['voltageTrigger']:
			self.triggerWaveVect = self.dicBraVectWaves['voltageTrigger']
		if self.hasBranch['voltageVeto']:
			self.vetoWaveVect = self.dicBraVectWaves['voltageVeto']
		if self.hasBranch['time']:
			self.timeVect = self.dicBraVectWaves['time']
		if self.hasBranch['event']:
			self.eventVect = self.dicBraVectWaves['event']
		if self.hasBranch['vetoedEvent']:
			self.vetoedEvent = self.dicBraVectWaves['vetoedEvent']
		if self.hasBranch['badShape']:
			self.badShape = self.dicBraVectWaves['badShape']
		if self.hasBranch['badPedestal']:
			self.badPedestal = self.dicBraVectWaves['badPedestal']
		if self.hasBranch['voltageHV']:
			self.voltageHV = self.dicBraVectWaves['voltageHV']
		if self.hasBranch['currentHV']:
			self.currentHV = self.dicBraVectWaves['currentHV']
		if self.hasBranch['timeHV']:
			key = 'timeHV.Convert()' if 'timeHV.Convert()' in self.branches1D else 'timeHV.AsDouble()'
			self.timeHV = self.dicBraVectWaves[key]

	def ExtractMeanOfWaveforms(self):
		self.signalWaveMeanVect = self.signalWaveVect.mean(axis=0)
		self.signalWaveSigmaVect = self.signalWaveVect.std(axis=0)

	def ReadInputFile(self):
		parser = ConfigParser()
		if self.infile != 'None':
			if os.path.isfile(self.infile):
				print 'Reading input file: {f} ...'.format(f=self.infile)
				parser.read(self.infile)
				if parser.has_section('OPTILINK'):
					if parser.has_option('OPTILINK', 'link'):
						self.optlink = parser.getint('OPTILINK', 'link')
					if parser.has_option('OPTILINK', 'node'):
						self.node = parser.getint('OPTILINK', 'node')
					if parser.has_option('OPTILINK', 'vme_base_address'):
						self.vme_b_addr = parser.getint('OPTILINK', 'vme_base_address')
					if parser.has_option('OPTILINK', 'wavedump_path'):
						self.wd_path = parser.get('OPTILINK', 'wavedump_path')
				if parser.has_section('RUN'):
					if parser.has_option('RUN', 'signal_channel') and self.sigCh == -1:
						self.sigCh = parser.getint('RUN', 'signal_channel')
					if parser.has_option('RUN', 'trig_channel') and self.trigCh == -1:
						self.trigCh = parser.getint('RUN', 'trig_channel')
					if parser.has_option('RUN', 'trig_val'):
						self.trigVal = parser.getfloat('RUN', 'trig_val')
					if parser.has_option('RUN', 'time'):
						self.points = int(np.ceil(parser.getfloat('RUN', 'time') * 1e-6 / self.time_res))
					if parser.has_option('RUN', 'post_trigger_percent'):
						self.post_trig_percent = parser.getint('RUN', 'post_trigger_percent')
				if parser.has_section('OUTPUT'):
					if parser.has_option('OUTPUT', 'dir') and self.outdir == 'None':
						self.outdir = parser.get('OUTPUT', 'dir')
					if parser.has_option('OUTPUT', 'prefix'):
						self.prefix = parser.get('OUTPUT', 'prefix')
					if parser.has_option('OUTPUT', 'suffix'):
						self.suffix = parser.get('OUTPUT', 'suffix')
			else:
				print 'Input file {f} does not exist. Loading default values...'.format(f=self.infile)
				self.LoadDefaults()
		else:
			self.LoadDefaults()

	def LoadDefaults(self):
		self.optlink = 1
		self.node = 0
		self.vme_b_addr = 32100000
		self.sigCh = 3 if self.sigCh == -1 else self.sigCh
		self.trigCh = 7 if self.trigCh == -1 else self.trigCh
		self.trigVal = 0.7138
		self.points = 5000
		self.post_trig_percent = 90

	def SetOutputFilesNames(self):
		if not os.path.isdir('{dir}/Runs'.format(dir=self.outdir)):
			os.makedirs('{dir}/Runs'.format(dir=self.outdir))

		def AddSuffix(string1):
			string1 += '_Input' if self.calv else '_Output'
			string1 += '_Pos' if self.bias >= 0 else '_Neg'
			string1 += '_{b}mV'.format(b=abs(1000 * self.bias))
			string1 += self.suffix
			return string1

		self.outString3 = '{p}_{f}'.format(dir=self.outdir, p=self.prefix, f=self.filename)
		self.outString3 = AddSuffix(self.outString3)
		self.filename = self.outString3

	def TakeTwoWaves(self):
		t0 = time.time()
		self.SetupDigitiser()
		print 'Starting getting data using wavedump...'
		p = subp.Popen(['wavedump', '{d}/WaveDumpConfig_CCD_cal.txt'.format(d=self.outdir)], bufsize=-1, stdin=subp.PIPE)
		t1 = time.time()
		# while p.poll() is None:
		self.Delay(4)
		p.stdin.write('c')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('W')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('P')
		p.stdin.flush()
		self.Delay(1.5)
		p.stdin.write('s')
		p.stdin.flush()
		while p.poll() is None:
			continue
		t0 = time.time() - t0
		print 'Total time saving {m} events:'.format(m=self.meas), t0, 'seconds'
		self.CreateRootFile()
		self.MoveBinaryFiles()

	def Delay(self, ti=1.0):
		t0 = time.time()
		while time.time() - t0 < ti:
			continue
		return

	def SetupDigitiser(self):
		print 'Creating digitiser CAEN V1730D configuration file... ', ; sys.stdout.flush()
		rfile = open('{d}/WaveDumpConfig_CCD_cal.txt'.format(d=self.outdir), 'w')
		rfile.write('[COMMON]')
		rfile.write('\n\n# open the digitezer')
		rfile.write('\nOPEN PCI {ol} {n} {ba}'.format(ol=int(self.optlink), n=int(self.node), ba=int(self.vme_b_addr)))
		rfile.write('\n\n# GNUPLOT path, normally /usr/bin/')
		rfile.write('\nGNUPLOT_PATH\t"/usr/bin/"')
		rfile.write('\n\n# output format can be BINARY or ASCII')
		rfile.write('\nOUTPUT_FILE_FORMAT\tBINARY')
		rfile.write('\n\n# if OUTPUT_FILE_HEADER is YES, the structure of the event has to be changed in self.struct_fmt to include the header')
		rfile.write('\nOUTPUT_FILE_HEADER NO')
		rfile.write('\n\n# specify the amount of samples to save. This defines the event window')
		rfile.write('\nRECORD_LENGTH\t{p}'.format(p=int(self.points)))
		rfile.write('\n\n# number of events to save in the file')
		rfile.write('\nMAX_NUM_EVENTS\t{n}'.format(n=int(self.meas)))
		rfile.write('\n\nTEST_PATTERN\tNO')
		rfile.write('\n\nENABLE_DES_MODE\tNO')
		rfile.write('\n\n# use external trigger. Options are: DISABLED, ACQUISITION_ONLY, ACQUISITION_AND_TRGOUT')
		rfile.write('\nEXTERNAL_TRIGGER\tDISABLED')
		rfile.write('\n\n# specify maximum number of events to read out in one Block Transfer. Must be between 1 and 1023')
		rfile.write('\nMAX_NUM_EVENTS_BLT\t100')
		rfile.write('\n\n# the percentage of the amount of data stored after the trigger in the event window. Has an offset of ~1.6... only accepts integers')
		rfile.write('\nPOST_TRIGGER\t{pt}'.format(pt=int(round(self.post_trig_percent*0.9996 - 1.6384))))
		rfile.write('\n\n# number of events that have to be ready before readout when the IRQ is asserted. 0 means run continuously. 1023 is the maximum')
		rfile.write('\nUSE_INTERRUPT\t0')
		rfile.write('\n\n# type of the fornt panel LEMO connectors: NIM, TTL')
		rfile.write('\n\nFPIO_LEVEL\tNIM')
		rfile.write('\n\nSKIP_STARTUP_CALIBRATION\tNO')
		rfile.write('\n\nCHANNEL_TRIGGER\tDISABLED')

		sig_polarity = 'POSITIVE' if self.bias >= 0 else 'NEGATIVE'

		rfile.write('\n\n# configuration for each channel [0] to [7]')
		for ch in xrange(16):
			rfile.write('\n\n[{ch}]'.format(ch=ch))
			if ch == self.sigCh or ch == self.trigCh:
				rfile.write('\nENABLE_INPUT\tYES')
			else:
				rfile.write('\nENABLE_INPUT\tNO')
			if ch == self.sigCh:
				rfile.write('\nPULSE_POLARITY\t{sp}'.format(sp=sig_polarity))
				rfile.write('\nDC_OFFSET\t{o}'.format(o=self.sig_offset))
				rfile.write('\nCHANNEL_TRIGGER\tDISABLED')
			if ch == self.trigCh:
				rfile.write('\nPULSE_POLARITY\tPOSITIVE')
				rfile.write('\nDC_OFFSET\t{o}'.format(o=self.trig_offset))
				rfile.write('\nCHANNEL_TRIGGER\tACQUISITION_ONLY')
				rfile.write('\nTRIGGER_THRESHOLD\t{th}'.format(th=int(round(np.divide(self.trigVal+1-self.trig_offset/50.0, self.sigRes)))))
		rfile.write('\n')
		rfile.close()
		print 'Done'

	def CreateRootFile(self):
		print 'Start creating root file'
		t0 = time.time()
		self.rawFile = ro.TFile('{d}/Runs/{r}.root'.format(d=self.outdir, r=self.filename), 'RECREATE')
		self.treeRaw = ro.TTree(self.filename, self.filename)
		fs = open('wave{s}.dat'.format(s=self.sigCh), 'rb')
		ft = open('wave{t}.dat'.format(t=self.trigCh), 'rb')
		eventBra = np.zeros(1, 'I')
		voltBra = np.zeros(self.points, 'f8')
		timeBra = np.zeros(self.points, 'f8')
		self.treeRaw.Branch('event', eventBra, 'event/i')
		self.treeRaw.Branch('time', timeBra, 'time[{s}]/D'.format(s=self.points))
		self.treeRaw.Branch('voltageSignal', voltBra, 'voltageSignal[{s}]/D'.format(s=self.points))
		self.CreateProgressBar(self.meas)
		self.bar.start()
		for ev in xrange(self.meas):
			fs.seek(ev * self.struct_len)
			ft.seek(ev * self.struct_len)
			datas = fs.read(self.struct_len)
			datat = ft.read(self.struct_len)
			if not datas or not datat:
				print 'No event in files... exiting'
				exit()
			s = struct.Struct(self.struct_fmt).unpack_from(datas)
			signalADCs = np.array(s, 'H')
			t = struct.Struct(self.struct_fmt).unpack_from(datat)
			triggADCs = np.array(t, 'H')
			signalVolts = np.array(np.multiply(signalADCs, self.sigRes) + self.sig_offset / 50.0 - 1, 'f8')
			triggVolts = np.array(np.multiply(triggADCs, self.sigRes) + self.trig_offset / 50.0 - 1, 'f8')
			left, right = np.double(triggVolts[0]), np.double(triggVolts[-1])
			mid = np.double((left + right) / 2.0)
			distFromMid = np.array(np.abs(triggVolts - mid), 'f8')
			midPos = distFromMid.argmin()
			timeVect = np.linspace(-midPos * self.time_res, self.time_res * (self.points - 1 - midPos), self.points, dtype='f8')
			eventBra.fill(ev)
			np.putmask(timeBra, 1 - np.zeros(self.points, '?'), timeVect)
			np.putmask(voltBra, 1 - np.zeros(self.points, '?'), signalVolts)
			numFil = self.treeRaw.Fill()
			self.bar.update(ev + 1)
		self.bar.finish()
		self.rawFile.Write()
		self.rawFile.Close()
		fs.close()
		ft.close()
		t0 = time.time() - t0
		print 'Time creating root tree:', t0, 'seconds'

	def MoveBinaryFiles(self):
		print 'Moving binary files... ', ; sys.stdout.flush()
		shutil.move('wave{chs}.dat'.format(chs=self.sigCh), '{d}/Runs/{f}_signal.dat'.format(d=self.outdir, f=self.filename))
		shutil.move('wave{cht}.dat'.format(cht=self.trigCh), '{d}/Runs/{f}_trigger.dat'.format(d=self.outdir, f=self.filename))
		print 'Done'

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


if __name__ == '__main__':
	print 'blaaaa'
