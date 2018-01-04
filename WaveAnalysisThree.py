#!/usr/bin/env python
import csv
import numpy as np
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
from ROOT import TProfile2D, TTree, gStyle, TFile, vector, gDirectory, TCut, TBranch, TH2F, TCanvas, TH2D, TH1D, TGraph, TGraphErrors
from scipy.optimize import curve_fit
import ipdb
# from TektronixReadout import TektronixReadout

class WaveAnalysisThree:
	def __init__(self, outputDir, inputfile, bias, cal, verbose):
		self.outDir, self.inputFile, self.bias, self.cal, self.verb = outputDir, inputfile, bias, cal, verbose
		self.fileRaw, self.treeRaw = None, None
		self.ptsWave, self.event, self.events = 0, np.zeros(1, 'I'), 0
		self.ped = np.zeros(1, 'f')
		self.rawtreename = 'raw_tree_' if not self.cal else 'raw_tree_cal_'
		self.rawtreename += 'pos_' if self.bias >= 0 else 'neg_'
		self.rawtreename += '{b}V_3_waves'.format(b=abs(self.bias)) if not self.cal else '{b}mV_2_waves'.format(b=abs(self.bias * 1000))
		self.boolConverted = self.ConvertFile()
		self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'READ')
		self.treeRaw = self.fileRaw.Get(self.rawtreename)
		self.treeRaw.GetEntry(0)
		self.ptsWave = self.treeRaw.GetLeaf('time').GetLen() if not self.boolConverted else self.ptsWave
		self.events = self.events if self.boolConverted else self.treeRaw.GetEntries()
		self.timeVect, self.eventVect, self.acceptEvents = self.treeRaw.time, None, np.array(1-np.zeros(self.events, '?'), '?')
		self.timeVect = np.array([self.timeVect[i] for i in xrange(self.ptsWave)], 'f')
		self.treeRaw.Draw('event', '', 'goff')
		self.eventVect = self.treeRaw.GetVal(0)
		self.eventVect = np.array([self.eventVect[i] for i in xrange(self.events)], 'I')
		self.signalVect = np.zeros(self.ptsWave, 'f')
		self.signalTLow, self.signalTHigh = 0.5e-6, 5e-6
		self.pedestalIntegrationTime = 0.98e-6
		self.pedestalTEndPos = -20e-9
		self.peakTime = 2.097e-6 if self.bias >= 0 else 2.079e-6
		self.peakForward = self.pedestalIntegrationTime/2
		self.peakBackward = self.pedestalIntegrationTime/2
		self.ped, self.pedSigma = np.zeros(1, 'f'), np.zeros(1, 'f')
		self.sigAndPed, self.sigAndPedSigma = np.zeros(1, 'f'), np.zeros(1, 'f')
		self.pedVect = np.empty(0, 'f')
		self.sigAndPedVect = np.empty(0, 'f')
		self.sigVect = np.empty(0, 'f')
		self.signal, self.signalSigma = np.zeros(1, 'f'), np.zeros(1, 'f')
		self.meanWaveGraph = None
		self.voltWaves = None
		self.voltWavesMean = None
		self.voltWavesMeanErrors = None
		if self.boolConverted:
			if self.cal:
				# TODO: Don't consider if during pedestal abs(signal - mean) >= 4 sigma
				self.measCalIntegrationTime = 0.5e-6
				self.riseTime = 20e-9
			self.CalculateAcceptedEvents()
			self.CalculatePedestalIntegral()
			if self.cal:
				self.measCalV = self.bias
				self.CalculateCalibrationVoltage()
				self.CalculateSignalIntegral(self.peakBackward, self.peakForward)
		# self.fileRaw.Close()

	def ConvertFile(self):
		self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'NEW')
		if self.fileRaw.IsOpen():
			print 'The file has not been created. Converting csv file to ROOT file...'
			self.csvtree = TTree('cvsTree', 'cvsTree')
			t0 = time.time()
			if not self.cal:
				self.csvtree.ReadFile(self.inputFile, 'event:time:voltage1:voltage2:voltage3')
			else:
				self.csvtree.ReadFile(self.inputFile, 'event:time:voltage1:voltage2')
			t0 = time.time() - t0
			print 'Time reading file:', t0, 'seconds'
			self.entries = self.csvtree.GetEntries()
			self.csvtree.GetEntry(self.entries - 1)
			lastEvent = self.csvtree.event
			self.csvtree.GetEntry(0)
			firstEvent = self.csvtree.event
			self.events = int(lastEvent - firstEvent + 1)
			self.ptsWave = int(self.entries / (self.events))
			self.csvtree.SetEstimate(self.entries)
			leng = self.csvtree.Draw('time:voltage3:voltage1:voltage2', '', 'goff') if not self.cal else self.csvtree.Draw('time:voltage2:voltage1', '', 'goff')
			tempTime = self.csvtree.GetVal(0)
			tempVolt1 = self.csvtree.GetVal(1)
			if not self.cal:
				tempVolt2 = self.csvtree.GetVal(2)
				tempVolt3 = self.csvtree.GetVal(3)
			arrayTime = [[tempTime[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries / self.ptsWave))]
			arrayVolt1 = [[tempVolt1[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries / self.ptsWave))]
			if not self.cal:
				arrayVolt2 = [[tempVolt2[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries / self.ptsWave))]
				arrayVolt3 = [[tempVolt3[ev * self.ptsWave + t] for t in xrange(self.ptsWave)] for ev in xrange(int(self.entries / self.ptsWave))]
			self.timeVect = np.array(arrayTime, 'f4')
			self.voltVect1 = np.array(arrayVolt1, 'f4')
			if not self.cal:
				self.voltVect2 = np.array(arrayVolt2, 'f4')
				self.voltVect3 = np.array(arrayVolt3, 'f4')
			self.timeBra = np.zeros(self.ptsWave, 'f4')
			self.voltBra1 = np.zeros(self.ptsWave, 'f4')
			if not self.cal:
				self.voltBra2 = np.zeros(self.ptsWave, 'f4')
				self.voltBra3 = np.zeros(self.ptsWave, 'f4')
			self.fileRaw.Write()
			self.fileRaw.Close()
			self.event = np.zeros(1, 'I')
			self.CreateRawRootFile()
			return True
		else:
			print 'The file already exists. Loading existing file...'
			del self.fileRaw
			del self.treeRaw
			return False

	def CalculateAcceptedEvents(self):
		self.CalculateGoodPedestalCut()

	def CalculateGoodPedestalCut(self):
		print 'Selecting good events based on pedestal'
		timeCut1, timeCut2 = TCut('time>{vl}'.format(vl=(-self.pedestalTEndPos-self.pedestalIntegrationTime))), TCut('time<{vh}'.format(vh=(-self.pedestalTEndPos)))
		timeCuts = timeCut1+timeCut2
		self.treeRaw.SetBranchStatus('*', 0)
		self.treeRaw.SetBranchStatus('voltageSignal', 1)
		self.treeRaw.SetBranchStatus('time', 1)
		num = self.treeRaw.Draw('voltageSignal>>htemp', timeCuts, 'goff')
		htemp = gDirectory.Get('htemp')
		mean = htemp.GetMean()
		sigm = htemp.GetStdDev()
		meansq = mean ** 2 + sigm * (num - 1)/num
		for i in xrange(5):
			for ev in self.eventVect:
				if self.acceptEvents[ev]:
					numev = self.treeRaw.Draw('voltageSignal>>htemp', timeCuts, 'goff', 1, ev)
					htemp = gDirectory.Get('htemp')
					evMean = htemp.GetMean()
					evSigm = htemp.GetStdDev()
					evMeansq = evMean ** 2 + evSigm * (numev - 1)/numev
					if abs(evMean - mean) > 4 * sigm:
						self.eventVect[ev] = False
						mean = (mean * num - evMean * numev) / (num - numev)
						sigm = (num * meansq - numev * evMeansq - ((num * mean - numev * evMean) ** 2)/(num-numev))/(num-numev-1)
						num = num - numev
						meansq = mean ** 2 + sigm * (num - 1)/num

	def CalculatePedestalIntegral(self, saveFile=True):
		print 'Calculating pedestals...'
		timeCut1, timeCut2 = TCut('time>{vl}'.format(vl=(self.pedestalTEndPos-self.pedestalIntegrationTime))), TCut('time<{vh}'.format(vh=(self.pedestalTEndPos)))
		timeCuts = timeCut1+timeCut2
		self.treeRaw.SetBranchStatus('*', 0)
		self.treeRaw.SetBranchStatus('voltageSignal', 1)
		self.treeRaw.SetBranchStatus('time', 1)
		for ev in self.eventVect:
			if self.acceptEvents[ev]:
				numev = self.treeRaw.Draw('voltageSignal>>htemp', timeCuts, 'goff', 1, ev)
				htemp = gDirectory.Get('htemp')
				evMean = htemp.GetMean()
				evSigm = htemp.GetStdDev()
				evMeansq = evMean ** 2 + evSigm * (numev - 1) / numev
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), evMean)
			else:
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), -100)
		self.fileRaw.Close()
		if saveFile:
			self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'UPDATE')
			self.treeRaw = self.fileRaw.Get(self.rawtreename)
			pedBranch = self.treeRaw.Branch('pedestal', self.ped, 'pedestal/F')
			# self.treeRaw.SetBranchStatus('*', 0)
			self.treeRaw.SetBranchStatus('pedestal', 1)
			for ev in self.eventVect:
				self.treeRaw.GetEntry(ev)
				self.ped.fill(self.pedVect[ev])
				pedBranch.Fill()
			self.treeRaw.Write()
			self.treeRaw.Delete()
			self.fileRaw.Close()
		else:
			temp = (np.delete(self.pedVect, np.where(self.pedVect == -100)))
			return (temp.mean(), temp.std())

	def CalculateSignalIntegral(self, back, forth, saveTree=True):
		print 'Calculating signals...'
		self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'READ')
		self.treeRaw = self.fileRaw.Get(self.rawtreename)
		timeCut1, timeCut2 = TCut('time>{vl}'.format(vl=self.peakTime-back)), TCut('time<{vh}'.format(vh=self.peakTime+forth))
		timeCuts = timeCut1+timeCut2
		self.treeRaw.SetBranchStatus('*', 0)
		self.treeRaw.SetBranchStatus('voltageSignal', 1)
		self.treeRaw.SetBranchStatus('time', 1)
		for ev in self.eventVect:
			if self.acceptEvents[ev]:
				numev = self.treeRaw.Draw('voltageSignal>>htemp2', timeCuts, 'goff', 1, ev)
				htemp = gDirectory.Get('htemp2')
				evMean = htemp.GetMean()
				evSigm = htemp.GetStdDev()
				evMeanMinusPed = evMean - self.pedVect[ev]
				self.sigAndPedVect = np.insert(self.sigAndPedVect, len(self.sigAndPedVect), evMean)
				self.sigVect = np.insert(self.sigVect, len(self.sigVect), evMeanMinusPed)
			else:
				self.sigAndPedVect = np.insert(self.sigAndPedVect, len(self.sigAndPedVect), -100)
				self.sigVect = np.insert(self.sigVect, len(self.sigVect), -100)
		self.fileRaw.Close()
		if saveTree:
			self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'UPDATE')
			self.treeRaw = self.fileRaw.Get(self.rawtreename)
			sigPedBranch = self.treeRaw.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
			sigBranch = self.treeRaw.Branch('signal', self.signal, 'signal/F')
			for ev in self.eventVect:
				self.treeRaw.GetEntry(ev)
				self.sigAndPed.fill(self.sigAndPedVect[ev])
				self.signal.fill(self.sigVect[ev])
				sigPedBranch.Fill()
				sigBranch.Fill()
			self.treeRaw.Write()
			self.treeRaw.Delete()
			self.fileRaw.Close()
		else:
			return (self.sigVect.mean(), self.sigVect.std())

	def SNRCalculation(self, left, right):
		self.AnalysisAllWaves(False, left, right)
		print '{l} left {r} right '.format(l=left, r=right)

	def SNRMap(self):
		if self.voltWaves is None:
			self.GetTimeVoltageWaveformsVectors()
		if self.voltWavesMean is None:
			self.ExtractMeanOfWaveforms()
		self.fileRaw.Close()
		self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'READ')
		self.treeRaw = self.fileRaw.Get(self.rawtreename)
		name = 'snrMap_{v}mV'.format(v=self.bias*1000)
		nameline = 'snrline_{v}mV'.format(v=self.bias*1000)
		filebla = TFile('{n}.root'.format(n=name), 'RECREATE')
		# blacanv = TCanvas('blac', 'blac', 1)
		self.fileRaw.cd()
		minbl = -0.1e-6
		maxbl = 1e-6
		binsl = int((maxbl-minbl)/0.01e-6 + 1)
		minbr = -0.1e-6
		maxbr = 1e-6
		binsr = int((maxbr-minbr)/0.01e-6 + 1)
		left = np.linspace(minbl, maxbl, binsl, dtype='f8')
		right = np.linspace(minbr, maxbr, binsr, dtype='f8')
		self.snrMap = TH2D(name, name, binsl, minbl-(maxbl-minbl)/(2*(binsl-1)), maxbl+(maxbl-minbl)/(2*(binsl-1)), binsr, minbr-(maxbr-minbr)/(2*(binsr-1)), maxbr+(maxbr-minbr)/(2*(binsr-1)))
		lineValues = [2, 5, 10, 20, 30, 40 , 50 , 60 , 70 , 80 , 90, 95, 98]
		self.snrLines = {val: TH1D('{nl}_{v}'.format(nl=nameline, v=val), '{nl}_{v}'.format(nl=nameline, v=val), binsl, minbl-(maxbl-minbl)/(2*(binsl-1)), maxbl+(maxbl-minbl)/(2*(binsl-1))) for val in lineValues}
		for i in xrange(len(left)):
			for j in xrange(len(right)):
				self.pedestalIntegrationTime = left[i]+right[j]
				self.pedestalTEndPos = -20e-9
				if 1e-6 + self.pedestalTEndPos >= self.pedestalIntegrationTime > 1e-10:
					self.SNRCalculation(left[i], right[j])
					if self.pedSigma[0] != 0:
						self.snrMap.SetBinContent(i+1, j+1, abs(float(self.signal[0])/float(self.pedSigma[0])))
						for val in lineValues:
							if abs(self.pedestalIntegrationTime - val*1e-8) <= 5e-9:
								self.snrLines[val].SetBinContent(i+1, abs(float(self.signal[0])/float(self.pedSigma[0])))
					else:
						self.snrMap.SetBinContent(i+1, j+1, 0)
				else:
					self.snrMap.SetBinContent(i + 1, j + 1, 0)
		filebla.cd()
		self.snrMap.Write()
		for val in lineValues:
			self.snrLines[val].Write()
		filebla.Write()
		filebla.Close()

	def CalculateCalibrationVoltage(self):
		pass

	def CreateRawRootFile(self):
		self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'UPDATE')
		self.treeRaw = TTree(self.rawtreename, self.rawtreename)
		self.treeRaw.Branch('event', self.event, 'event/i')
		self.treeRaw.Branch('time', self.timeBra, 'time[{s}]/F'.format(s=self.ptsWave))
		self.treeRaw.Branch('voltageSignal', self.voltBra1, 'voltageSignal[{s}]/F'.format(s=self.ptsWave))
		if not self.cal:
			self.treeRaw.Branch('voltageTrigger1', self.voltBra2, 'voltageTrigger1[{s}]/F'.format(s=self.ptsWave))
			self.treeRaw.Branch('voltageTrigger2', self.voltBra3, 'voltageTrigger2[{s}]/F'.format(s=self.ptsWave))
		t1 = time.time()
		for i in xrange(int(self.events)):
			self.event.fill(i)
			np.putmask(self.timeBra, 1 - np.zeros(self.ptsWave, '?'), self.timeVect[i])
			np.putmask(self.voltBra1, 1 - np.zeros(self.ptsWave, '?'), self.voltVect1[i])
			if not self.cal:
				np.putmask(self.voltBra2, 1 - np.zeros(self.ptsWave, '?'), self.voltVect2[i])
				np.putmask(self.voltBra3, 1 - np.zeros(self.ptsWave, '?'), self.voltVect3[i])
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

	def GetTimeVoltageWaveformsVectors(self):
		self.GetTimeVector()
		self.GetVoltageWaveforms()

	def GetTimeVector(self):
		self.treeRaw.SetBranchStatus('*', 0)
		self.treeRaw.SetBranchStatus('time', 1)
		self.pointsWave = self.treeRaw.Draw('time', '', 'goff', 1, 0)
		if self.pointsWave > self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(self.pointsWave)
			self.pointsWave = self.treeRaw.Draw('time', '', 'goff', 1)
		self.timeVect = self.treeRaw.GetVal(0)
		self.timeVect = np.array([self.timeVect[t] for t in xrange(self.pointsWave)], 'f8')

	def GetVoltageWaveforms(self):
		self.treeRaw.SetBranchStatus('*', 0)
		self.treeRaw.SetBranchStatus('voltageSignal', 1)
		numelem = self.treeRaw.Draw('voltageSignal', '', 'goff')
		if numelem > self.treeRaw.GetEstimate():
			self.treeRaw.SetEstimate(numelem)
			numelem = self.treeRaw.Draw('voltageSignal', '', 'goff')
		self.voltWaves = self.treeRaw.GetVal(0)
		self.voltWaves = np.array([[self.voltWaves[ev * self.pointsWave + t] for t in xrange(self.pointsWave)] for ev in xrange(self.events)], 'f8')

	def ExtractMeanOfWaveforms(self):
		if self.voltWaves is None:
			self.GetTimeVoltageWaveformsVectors()
		self.voltWavesMean = self.voltWaves.mean(axis=0)
		self.voltWavesMeanErrors = self.voltWaves.std(axis=0)

	def CreateMeanWaveformGraph(self):
		self.meanWaveGraph = TGraphErrors(self.pointsWave, self.timeVect, self.voltWavesMean, np.zeros(self.pointsWave, 'f8'), self.voltWavesMeanErrors)

	def FindPedestalPositions(self, left, right):
		self.pedestalIntegrationTime = left + right
		self.pedestalTimeIndices = np.argwhere(np.bitwise_and(self.pedestalTEndPos - self.pedestalIntegrationTime <= self.timeVect, self.timeVect <= self.pedestalTEndPos)).flatten()

	def FindSignalPosition(self, left, right):
		self.pedestalIntegrationTime = left + right
		self.signalTimeIndices = np.argwhere(abs(self.timeVect - self.peakTime + (right - left)/2.0) <= self.pedestalIntegrationTime/2.0).flatten()

	def FindSignalPositionReal(self, left, right):
		self.pedestalIntegrationTime = left + right
		self.signalTimeIndicesReal = np.argwhere(abs(self.timeVect - self.peakTimeReal + (right - left)/2.0) <= self.pedestalIntegrationTime/2.0).flatten()

	def FindRealPeakPosition(self):
		if self.meanWaveGraph is None:
			self.CreateMeanWaveformGraph()
		mpos = self.voltWavesMean.argmin() if self.bias >= 0 else self.voltWavesMean.argmax()
		xmin, xmax = self.timeVect[mpos] - 200e-9, self.timeVect[mpos] + 200e-9
		fit = self.meanWaveGraph.Fit('pol2', 'QMN0FS', '', xmin, xmax)
		# ipdb.set_trace()
		c, b, a = fit.Parameter(0), fit.Parameter(1), fit.Parameter(2)
		return -b/(2*a)

	def AnalysisMeanWaves(self):
		if self.voltWaves is None:
			self.ExtractMeanOfWaveforms()
		self.peakBackward, self.peakForward = 0.469e-6, 0.511e-6
		self.peakTime = 2.127e-6 if self.bias >= 0 else 2.1195e-6
		self.peakTimeReal = self.FindRealPeakPosition()
		self.FindPedestalPositions(self.peakBackward, self.peakForward)
		self.FindSignalPosition(self.peakBackward, self.peakForward)
		self.FindSignalPositionReal(self.peakBackward, self.peakForward)
		self.pedestalMeanWaveAverage = self.voltWavesMean[self.pedestalTimeIndices].mean()
		self.signalAndPedestalMeanWaveAverage = self.voltWavesMean[self.signalTimeIndices].mean()
		self.signalMeanWaveAverage = self.signalAndPedestalMeanWaveAverage - self.pedestalMeanWaveAverage
		self.signalAndPedestalMeanWaveAverageReal = self.voltWavesMean[self.signalTimeIndicesReal].mean()
		self.signalMeanWaveAverageReal = self.signalAndPedestalMeanWaveAverageReal - self.pedestalMeanWaveAverage

	def CalculateSignalsAllWaves(self, saveFile=True):
		print 'Calculating pedestals and signals...'
		for ev in self.eventVect:
			if self.acceptEvents[ev]:
				evPedMean = self.voltWaves[ev, self.pedestalTimeIndices].mean() if self.pedestalTimeIndices.size >= 1 else -10
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), evPedMean)
				evMeasMean = self.voltWaves[ev, self.signalTimeIndices].mean() if self.signalTimeIndices.size >= 1 else -10
				evMeasMeanReal = self.voltWaves[ev, self.signalTimeIndicesReal].mean() if self.signalTimeIndicesReal.size >= 1 else -10
				evSignalMean = evMeasMean - evPedMean
				evSignalMeanReal = evMeasMeanReal - evPedMean
				self.sigAndPedVect = np.insert(self.sigAndPedVect, len(self.sigAndPedVect), evMeasMean)
				self.sigAndPedVectReal = np.insert(self.sigAndPedVectReal, len(self.sigAndPedVectReal), evMeasMeanReal)
				self.sigVect = np.insert(self.sigVect, len(self.sigVect), evSignalMean)
				self.sigVectReal = np.insert(self.sigVectReal, len(self.sigVectReal), evSignalMeanReal)
			else:
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), -10)
				self.sigAndPedVect = np.insert(self.sigAndPedVect, len(self.sigAndPedVect), -10)
				self.sigAndPedVectReal = np.insert(self.sigAndPedVectReal, len(self.sigAndPedVectReal), -10)
				self.sigVect = np.insert(self.sigVect, len(self.sigVect), -10)
				self.sigVectReal = np.insert(self.sigVectReal, len(self.sigVectReal), -10)
		if saveFile:
			if self.fileRaw.IsOpen():
				self.fileRaw.Close()
			self.fileRaw = TFile('{o}/{r}.root'.format(o=self.outDir, r=self.rawtreename), 'UPDATE')
			self.treeRaw = self.fileRaw.Get(self.rawtreename)
			pedBranch = self.treeRaw.Branch('pedestal', self.ped, 'pedestal/F')
			sigPedBranch = self.treeRaw.Branch('signalAndPedestal', self.sigAndPed, 'signalAndPedestal/F')
			sigPedRealBranch = self.treeRaw.Branch('signalAndPedestalReal', self.sigAndPedReal, 'signalAndPedestalReal/F')
			sigBranch = self.treeRaw.Branch('signal', self.signal, 'signal/F')
			sigRealBranch = self.treeRaw.Branch('signalReal', self.signalReal, 'signalReal/F')
			self.treeRaw.SetBranchStatus('pedestal', 1)
			self.treeRaw.SetBranchStatus('signalAndPedestal', 1)
			self.treeRaw.SetBranchStatus('signalAndPedestalReal', 1)
			self.treeRaw.SetBranchStatus('signal', 1)
			self.treeRaw.SetBranchStatus('signalReal', 1)
			for ev in self.eventVect:
				self.treeRaw.GetEntry(ev)
				self.ped.fill(self.pedVect[ev])
				self.sigAndPed.fill(self.sigAndPedVect[ev])
				self.sigAndPedReal.fill(self.sigAndPedVectReal[ev])
				self.signal.fill(self.sigVect[ev])
				self.signalReal.fill(self.sigVectReal[ev])
				pedBranch.Fill()
				sigPedBranch.Fill()
				sigPedRealBranch.Fill()
				sigBranch.Fill()
				sigRealBranch.Fill()
			self.treeRaw.Write()
			self.treeRaw.Delete()
			self.fileRaw.Close()
		else:
			tempPed = np.delete(self.pedVect, np.where(self.pedVect == -10))
			tempSAP = np.delete(self.sigAndPedVect, np.where(self.sigAndPedVect == -10))
			tempSAPR = np.delete(self.sigAndPedVectReal, np.where(self.sigAndPedVectReal == -10))
			tempS = np.delete(self.sigVect, np.where(self.sigVect == -10))
			tempSR = np.delete(self.sigVectReal, np.where(self.sigVectReal == -10))
			if tempPed.size > 1 and tempSAP.size > 0 and tempS.size > 0 and tempSAPR.size > 0 and tempSR.size > 0:
				self.signal.fill(tempS.mean())
				self.pedSigma.fill(tempPed.std())
			else:
				self.signal.fill(0)
				self.pedSigma.fill(0)

	def AnalysisAllWaves(self, saveFile=True, left=0.47e-6, right=0.51e-6):
		if self.voltWaves is None:
			self.GetTimeVoltageWaveformsVectors()
		if self.voltWavesMean is None:
			self.ExtractMeanOfWaveforms()

		self.pedVect = np.empty(0, 'f8')
		self.sigAndPedVect, self.sigAndPedVectReal = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.sigVect, self.sigVectReal = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.ped, self.pedSigma = np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.sigAndPed, self.sigAndPedSigma = np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.sigAndPedReal = np.zeros(1, 'f8')
		self.signal, self.signalSigma = np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.signalReal = np.zeros(1, 'f8')

		self.peakBackward, self.peakForward = left, right
		self.peakTime = 2.1270e-6 if self.bias >= 0 else 2.1185e-6
		self.peakTimeReal = self.FindRealPeakPosition()
		self.FindPedestalPositions(self.peakBackward, self.peakForward)
		self.FindSignalPosition(self.peakBackward, self.peakForward)
		self.FindSignalPositionReal(self.peakBackward, self.peakForward)
		self.CalculateSignalsAllWaves(saveFile)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string',
					  help='Relative path to output directory e.g. ./')
	parser.add_option('-i', '--infile', dest='infile', type='string',
					  help='path to file to convert and analyse')
	parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true',
					  help='Toggles automatic analysis')
	parser.add_option('-b', '--bias', dest='bias', default=-400, type='float',
					  help='The HV bias used on the sample e.g. -400')
	parser.add_option('-c', '--calib', dest='calib', default=False, action='store_true',
					  help='Enables the the two wave analysis for callibration. Default: False')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	outdir = str(options.outdir)
	infile = str(options.infile)
	automatic = bool(options.automatic)
	bias = float(options.bias)
	cal = bool(options.calib)
	verb = bool(options.verb)
	z = WaveAnalysisThree(outdir, infile, bias, cal, verb)
	if automatic:
		pass
