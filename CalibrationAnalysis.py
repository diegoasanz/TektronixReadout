#!/usr/bin/env python
import csv
import numpy as np
import ROOT as ro
import scipy.constants as sc
from struct import unpack
import time, os, sys
from optparse import OptionParser
import progressbar
import ConfigParser
import pymouse as pym
import pykeyboard as pyk

# from ROOT import TProfile2D, TTree, gStyle, TFile, vector, gDirectory, TCut, TBranch, TH2F, TCanvas, TH2D, TH1D, TGraph, TGraphErrors
#from scipy.optimize import curve_fit
import ipdb
# from TektronixReadout import TektronixReadout

class CalibrationAnalysis:
	def __init__(self, outputDir, inputfile, bias, verbose):
		ro.gErrorIgnoreLevel = ro.kError
		self.time0 = time.time()
		self.outDir, self.inputFile, self.vcal, self.verb = outputDir, inputfile, bias, verbose
		self.ptsWaveOutput, self.eventOutput, self.eventsOutput = 0, np.zeros(1, 'I'), 0
		self.ptsWaveInput, self.eventInput, self.eventsInput = 0, np.zeros(1, 'I'), 0
		self.eventVectOutput, self.timeVectOutput, self.wavesVectOutput, self.wavesMeanVectOutput, self.wavesMeanErrorsVectOutput = None, None, None, None, None
		self.eventVectInput, self.timeVectInput, self.wavesVectInput, self.wavesMeanVectInput, self.wavesMeanErrorsVectInput = None, None, None, None, None
		self.pedVect, self.pedSigmaVect, self.pedSignalVect, self.pedSignalRealVect, self.signalVect, self.signalRealVect = None, None, None, None, None, None
		self.pedInput1Vect, self.pedInput2Vect, self.pedInput3Vect, self.pedInput4Vect, self.pedInputSigma1Vect, self.pedInputSigma2Vect, self.pedInputSigma3Vect, self.pedInputSigma4Vect = None, None, None, None, None, None, None, None
		self.pedCal1Vect, self.pedCal2Vect, self.pedCal3Vect, self.pedCal4Vect, self.calVoltsReal1Vect, self.calVoltsReal2Vect, self.calVoltsReal3Vect, self.calVoltsReal4Vect = None, None, None, None, None, None, None, None
		self.ped, self.pedSigma, self.pedSignal, self.pedSignalReal, self.signal, self.signalReal = {}, {}, {}, {}, {}, {}
		self.pedInput1, self.pedInput2, self.pedInput3, self.pedInput4, self.pedInputSigma1, self.pedInputSigma2, self.pedInputSigma3, self.pedInputSigma4 = {}, {}, {}, {}, {}, {}, {}, {}
		self.pedInput1Pos, self.pedInput2Pos, self.pedInput3Pos, self.pedInput4Pos = None, None, None, None
		self.pedCal1Pos, self.pedCal2Pos, self.pedCal3Pos, self.pedCal4Pos = None, None, None, None
		self.pedCal1, self.pedCal2, self.pedCal3, self.pedCal4, self.calVoltsReal1, self.calVoltsReal2, self.calVoltsReal3, self.calVoltsReal4 = {}, {}, {}, {}, {}, {}, {}, {}
		self.calVoltsRealSigma1, self.calVoltsRealSigma2, self.calVoltsRealSigma3, self.calVoltsRealSigma4 = {}, {}, {}, {}
		self.meanInputWaveGraph, self.meanOutputWaveGraph, self.vCalSignalGraph, self.chargeSignalGraph = None, None, None, None
		self.vCalSignalFit, self.chargeSignalFit = None, None
		self.bar = None
		self.timePeakReal, self.posHalfCal = None, None
		self.fileCal = None

		self.boolInputConverted, self.boolInputHasTreeCSV, self.boolInputHasVectors, self.boolInputHasScalars, self.boolInputIsAnalysed = {}, {}, {}, {}, {}
		self.boolOutputConverted, self.boolOutputHasTreeCSV, self.boolOutputHasVectors, self.boolOutputHasScalars, self.boolOutputIsAnalysed = {}, {}, {}, {}, {}
		self.fileOutputRaw, self.treeRawOutput, self.treeOutputCSV = None, None, None
		self.fileInputRaw, self.treeRawInput, self.treeInputCSV = None, None, None

		self.peakBack, self.peakForth = 469e-9, 511e-9
		self.averagingTime = self.peakBack + self.peakForth
		self.pedestalTEndPos = -20e-9
		self.vcal_syst_factor = 0.04

		self.calDist1, self.calDist2, self.calDist3, self.calDist4 = 10e-9, 20e-9, 50e-9, 1e-6

		self.calVolts, self.calVoltsReal, self.runInputFiles, self.runOutputFiles, self.runFilesDir = None, {}, {}, {}, None
		self.inputSuffix = self.outputSuffix = self.inputPrefix = self.outputPrefix = ''
		self.factInputData = self.factOutputData = 1
		self.ReadInputFile()

		self.timeBraOutput, self.voltBraOutput = None, None
		self.timeBraInput, self.voltBraInput = None, None
		self.rawInputTreeNames, self.rawOutputTreeNames = {}, {}
		self.pedestalTimeIndices, self.signalTimeIndices, self.signalTimeIndicesReal = None, None, None
		self.acceptEventsOutput, self.acceptEventsInput = None, None
		for vcal in self.calVolts:
			self.vcal = vcal
			print 'Loading data for vcal: {v}mV...'.format(v=1000*self.vcal)
			self.timePeak = 2.1270e-6 if self.vcal >= 0 else 2.1185e-6
			self.boolInputConverted[self.vcal], self.boolInputHasVectors[self.vcal], self.boolInputHasTreeCSV[self.vcal], self.boolInputHasScalars[self.vcal], self.boolInputIsAnalysed[self.vcal] = False, False, False, False, False
			self.boolOutputConverted[self.vcal], self.boolOutputHasVectors[self.vcal], self.boolOutputHasTreeCSV[self.vcal], self.boolOutputHasScalars[self.vcal], self.boolOutputIsAnalysed[self.vcal] = False, False, False, False, False
			self.ped[self.vcal], self.pedSigma[self.vcal], self.pedSignal[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.pedSignalReal[self.vcal], self.signal[self.vcal], self.signalReal[self.vcal] = np.zeros(1,'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.pedInput1[self.vcal], self.pedInput2[self.vcal], self.pedInput3[self.vcal], self.pedInput4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.pedInputSigma1[self.vcal], self.pedInputSigma2[self.vcal], self.pedInputSigma3[self.vcal], self.pedInputSigma4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.pedCal1[self.vcal], self.pedCal2[self.vcal], self.pedCal3[self.vcal], self.pedCal4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.calVoltsReal1[self.vcal], self.calVoltsReal2[self.vcal], self.calVoltsReal3[self.vcal], self.calVoltsReal4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.calVoltsRealSigma1[self.vcal], self.calVoltsRealSigma2[self.vcal], self.calVoltsRealSigma3[self.vcal], self.calVoltsRealSigma4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
			self.rawInputTreeNames[self.vcal] = 'raw_in_tree_cal_pos_{b}mV_waves'.format(b=abs(self.vcal * 1000)) if self.vcal >= 0 else 'raw_in_tree_cal_neg_{b}mV_waves'.format(b=abs(self.vcal * 1000))
			self.rawOutputTreeNames[self.vcal] = 'raw_out_tree_cal_pos_{b}mV_waves'.format(b=abs(self.vcal * 1000)) if self.vcal >= 0 else 'raw_out_tree_cal_neg_{b}mV_waves'.format(b=abs(self.vcal * 1000))
			self.ConvertFiles()
		self.time1 = time.time()
		tempt = self.time1 - self.time0
		print 'Time taken to load and prepare all files:', tempt, 'seconds'

#   TODO: MAKE IT DO IT FROM A TEXT FILE OR A JSON!!!!
	def ReadInputFile(self):

		parser = ConfigParser.ConfigParser()
		parser.read(self.inputFile)

		self.calVolts = np.array([-2, -1.589, -1.262, -1.002, -0.7962, -0.6325, -0.5024, -0.3991, -0.3170, -0.2518, -0.2000, -0.1589, -0.1262, -0.1002, -0.07962, -0.06325, -0.05024, -0.03991, -0.0317,
		                 -0.02518, -0.0200, -0.01589, -0.01262, -0.01002, -0.007962, -0.006325, -0.005024, 0.005024, 0.006325, 0.007962, 0.01002, 0.01262, 0.01589, 0.0200, 0.02518, 0.0317, 0.03991, 0.05024,
		                 0.06325, 0.07962, 0.1002, 0.1262, 0.1589, 0.2000, 0.2518, 0.3170, 0.3991, 0.5024, 0.6325, 0.7962, 1.002, 1.262, 1.589, 2], 'f8')

		if parser.has_section('CALIBRATES'):
			if parser.has_option('CALIBRATES', 'voltages'):
				calVoltsString = parser.get('CALIBRATES', 'voltages')
				calVoltsString = str.split(calVoltsString, ',')
				self.calVolts = np.array([np.double(calVoltsString[i]) for i in xrange(len(calVoltsString))], 'f8')

		self.runFilesDir = self.outDir
		self.inputPrefix = 'waves_Calibration_Input_'
		self.inputSuffix = 'mV.csv'

		self.outputPrefix = 'waves_Calibration_Output_'
		self.outputSuffix = 'mV.csv'

		if parser.has_section('INPUT'):
			if parser.has_option('INPUT', 'inputFilesDir'):
				self.runFilesDir = parser.get('INPUT', 'inputFilesDir')
			if parser.has_option('INPUT', 'prefix'):
				self.inputPrefix = parser.get('INPUT', 'prefix')
			if parser.has_option('INPUT', 'suffix'):
				self.inputSuffix = parser.get('INPUT', 'suffix')
			if parser.has_option('INPUT', 'factorAnalysis'):
				self.factInputData = parser.getfloat('INPUT', 'factorAnalysis')

		if parser.has_section('OUTPUT'):
			if parser.has_option('OUTPUT', 'prefix'):
				self.outputPrefix = parser.get('OUTPUT', 'prefix')
			if parser.has_option('OUTPUT', 'suffix'):
				self.outputSuffix = parser.get('OUTPUT', 'suffix')
			if parser.has_option('OUTPUT', 'factorAnalysis'):
				self.factOutputData = parser.getfloat('OUTPUT', 'factorAnalysis')

		for val in self.calVolts:
			if val >= 0:
				self.runInputFiles[val] = '{p}Pos_{v}{s}'.format(p=self.inputPrefix, v=abs(1000*val), s=self.inputSuffix)
				self.runOutputFiles[val] = '{p}Pos_{v}{s}'.format(p=self.outputPrefix, v=abs(1000*val), s=self.outputSuffix)
			else:
				self.runInputFiles[val] = '{p}Neg_{v}{s}'.format(p=self.inputPrefix, v=abs(1000*val), s=self.inputSuffix)
				self.runOutputFiles[val] = '{p}Neg_{v}{s}'.format(p=self.outputPrefix, v=abs(1000*val), s=self.outputSuffix)

		# self.runFiles = {val: '{p}{v:.1f}{s}'.format(p=self.inputPrefix, v=val*1000, s=self.inputSuffix) if abs(val) >= 0.1 else '{p}{v:.1f}{s}'.format(p=self.inputPrefix, v=val, s=self.inputSuffix) for val in self.calVolts}
		#
		#   500
		# 	self.runFiles = {-2: 'waves_Calibration_Neg_2000.0mV.csv',
		# 	                 -1.589: 'waves_Calibration_Neg_1589.0mV.csv',
		# 	                 -1.262: 'waves_Calibration_Neg_1262.0mV.csv',
		# 	                 -1.002: 'waves_Calibration_Neg_1002.0mV.csv',
		# 	                 -0.7962: 'waves_Calibration_Neg_796.2mV.csv',
		# 	                 -0.6325: 'waves_Calibration_Neg_632.5mV.csv',
		# 	                 -0.5024: 'waves_Calibration_Neg_502.4mV.csv',
		# 	                 -0.3991: 'waves_Calibration_Neg_399.1mV.csv',
		# 	                 -0.3170: 'waves_Calibration_Neg_317.0mV.csv',
		# 	                 -0.2518: 'waves_Calibration_Neg_251.8mV.csv',
		#   1000
		# 	                 -0.2000: 'waves_Calibration_Neg_200.0mV.csv',
		# 	                 -0.1589: 'waves_Calibration_Neg_158.9mV.csv',
		# 	                 -0.1262: 'waves_Calibration_Neg_126.2mV.csv',
		# 	                 -0.1002: 'waves_Calibration_Neg_100.2mV.csv',
		# 	                 -0.07962: 'waves_Calibration_Neg_79.62mV.csv',
		#   2000
		# 	                 -0.06325: 'waves_Calibration_Neg_63.25mV.csv',
		# 	                 -0.05024: 'waves_Calibration_Neg_50.24mV.csv',
		# 	                 -0.03991: 'waves_Calibration_Neg_39.91mV.csv',
		# 	                 -0.0317: 'waves_Calibration_Neg_31.7mV.csv',
		# 	                 -0.02518: 'waves_Calibration_Neg_25.18mV.csv',
		#   3000
		# 	                 -0.0200: 'waves_Calibration_Neg_20.0mV.csv',
		# 	                 -0.01589: 'waves_Calibration_Neg_15.89mV.csv',
		# 	                 -0.01262: 'waves_Calibration_Neg_12.62mV.csv',
		# 	                 -0.01002: 'waves_Calibration_Neg_10.02mV.csv',
		# 	                 -0.007962: 'waves_Calibration_Neg_7.962mV.csv',
		#   4000
		# 	                 -0.006325: 'waves_Calibration_Neg_6.325mV.csv',
		# 	                 -0.005024: 'waves_Calibration_Neg_5.024mV.csv',
		# 	                 0.005024: 'waves_Calibration_Pos_5.024mV.csv',
		# 	                 0.006325: 'waves_Calibration_Pos_6.325mV.csv',
		#   3000
		# 	                 0.007962: 'waves_Calibration_Pos_7.962mV.csv',
		# 	                 0.01002: 'waves_Calibration_Pos_10.02mV.csv',
		# 	                 0.01262: 'waves_Calibration_Pos_12.62mV.csv',
		# 	                 0.01589: 'waves_Calibration_Pos_15.89mV.csv',
		# 	                 0.0200: 'waves_Calibration_Pos_20.0mV.csv',
		#   2000                                                            TODO AQUI ESTAMOS
		# 	                 0.02518: 'waves_Calibration_Pos_25.18mV.csv',
		# 	                 0.0317: 'waves_Calibration_Pos_31.7mV.csv',
		# 	                 0.03991: 'waves_Calibration_Pos_39.91mV.csv',
		# 	                 0.05024: 'waves_Calibration_Pos_50.24mV.csv',  TODO quitar atenuador 20db
		# 	                 0.06325: 'waves_Calibration_Pos_63.25mV.csv',
		#   1000
		# 	                 0.07962: 'waves_Calibration_Pos_79.62mV.csv',
		# 	                 0.1002: 'waves_Calibration_Pos_100.2mV.csv',
		# 	                 0.1262: 'waves_Calibration_Pos_126.2mV.csv',
		# 	                 0.1589: 'waves_Calibration_Pos_158.9mV.csv',
		# 	                 0.2000: 'waves_Calibration_Pos_200.0mV.csv',
		#   500
		# 	                 0.2518: 'waves_Calibration_Pos_251.8mV.csv',
		# 	                 0.3170: 'waves_Calibration_Pos_317.0mV.csv',
		# 	                 0.3991: 'waves_Calibration_Pos_399.1mV.csv',
		# 	                 0.5024: 'waves_Calibration_Pos_502.4mV.csv',
		# 	                 0.6325: 'waves_Calibration_Pos_632.5mV.csv',
		# 	                 0.7962: 'waves_Calibration_Pos_796.2mV.csv',
		# 	                 1.002: 'waves_Calibration_Pos_1002.0mV.csv',
		# 	                 1.262: 'waves_Calibration_Pos_1262.0mV.csv',
		# 	                 1.589: 'waves_Calibration_Pos_1589.0mV.csv',
		# 	                 2: 'waves_Calibration_Pos_2000.0mV.csv'}

		self.calVoltsReal = {-2: -1.990, -1.589: -1.599, -1.262: -1.285, -1.002: -1.016, -0.7962: -8.101E-01, -0.6325: -6.393E-01, -0.5024: -5.078E-01, -0.3991: -4.021E-01, -0.3170: -3.232E-01,
		                     -0.2518: -2.570E-01, -0.2000: -2.031E-01, -0.1589: -1.609E-01, -0.1262: -1.297E-01, -0.1002: -1.024E-01, -0.07962: -8.142E-02, -0.06325: -6.491E-02, -0.05024: -5.122E-02,
		                     -0.03991: -4.142E-02, -0.0317: -3.290E-02, -0.02518: -2.608E-02, -0.0200: -2.065E-02, -0.01589: -1.643E-02, -0.01262: -1.293E-02, -0.01002: -1.025E-02,
		                     -0.007962: -8.155E-03, -0.006325: -6.510E-03, -0.005024: -5.152E-03, 0.005024: 4.883E-03, 0.006325: 6.238E-03, 0.007962: 7.883E-03, 0.01002: 1.020E-02, 0.01262: 1.284E-02,
		                     0.01589: 1.622E-02, 0.0200: 2.053E-02, 0.02518: 2.590E-02, 0.0317: 3.282E-02, 0.03991: 4.126E-02, 0.05024: 5.027E-02, 0.06325: 6.376E-02, 0.07962: 8.078E-02,
		                     0.1002: 1.009E-01, 0.1262: 1.272E-01, 0.1589: 1.599E-01, 0.2000: 2.017E-01, 0.2518: 2.530E-01, 0.3170: 3.191E-01, 0.3991: 4.020E-01, 0.5024: 5.106E-01, 0.6325: 6.331E-01,
		                     0.7962: 8.081E-01, 1.002: 1.012, 1.262: 1.280, 1.589: 1.599, 2: 1.996}

		# self.CalculateChargeFromVcal()
		
	def CalculateChargeFromVcal(self, realCal=4):
		def CalculateFactor():
			e_charge = sc.codata.value('elementary charge')
			r5a, r5b, r5c = 56.0, 467.0, 4.505  # r5a, r5b and r5c are input resistance, and the voltage divisor resistors respectively
			c9, cp = 1.8e-12, 0.2e-12  # c9 and cp are the input capacitance similar to sensor and the parasitic capacitance estimated by Ulf respectively
			tr5a = tr5b = tr5c = 0.012  # tolerances of the resistors r5a, r5b and r5c respectively
			tc9, tcp = 0.02, 0.2  # tolerances of the c9 capacitor and an estimate of the tolerance of Ulf's estimation respectively
			# Returns a tuple with the value of the factor and its uncertainty
			return (r5c * (c9 + cp)/(r5c + r5b)) / e_charge, (r5c * np.sqrt(r5b**2 * (c9 + cp)**2 * (tr5b**2 + tr5c**2)/(r5b + r5c)**2 + (c9 * tc9)**2 + (cp * tcp)**2) / (r5b + r5c)) / e_charge

		(fact, fact_uncert) = CalculateFactor()
		vcal_sys_percent = self.vcal_syst_factor  # systematic uncertainty due to raise time of the pulser

		if realCal == 1:
			self.calCharge = {vcal: self.calVoltsReal1[vcal]*fact for vcal in self.calVolts}
			self.calChargeSigma = {vcal: np.sqrt((self.calVoltsReal1[vcal]*fact_uncert)**2 + fact**2 * (self.calVoltsRealSigma1[vcal]**2 + (self.calVoltsReal1[vcal] * vcal_sys_percent)**2)) for vcal in self.calVolts}
		elif realCal == 2:
			self.calCharge = {vcal: self.calVoltsReal2[vcal]*fact for vcal in self.calVolts}
			self.calChargeSigma = {vcal: np.sqrt((self.calVoltsReal2[vcal]*fact_uncert)**2 + fact**2 * (self.calVoltsRealSigma2[vcal]**2 + (self.calVoltsReal2[vcal] * vcal_sys_percent)**2)) for vcal in self.calVolts}
		elif realCal == 3:
			self.calCharge = {vcal: self.calVoltsReal3[vcal]*fact for vcal in self.calVolts}
			self.calChargeSigma = {vcal: np.sqrt((self.calVoltsReal3[vcal]*fact_uncert)**2 + fact**2 * (self.calVoltsRealSigma3[vcal]**2 + (self.calVoltsReal3[vcal] * vcal_sys_percent)**2)) for vcal in self.calVolts}
		else:
			self.calCharge = {vcal: self.calVoltsReal4[vcal]*fact for vcal in self.calVolts}
			self.calChargeSigma = {vcal: np.sqrt((self.calVoltsReal4[vcal]*fact_uncert)**2 + fact**2 * (self.calVoltsRealSigma4[vcal]**2 + (self.calVoltsReal4[vcal] * vcal_sys_percent)**2)) for vcal in self.calVolts}

	def OpenInputFile(self, mode='READ'):
		if self.fileInputRaw is not None:
			if self.fileInputRaw.IsOpen:
				self.fileInputRaw.Close()
		self.fileInputRaw = ro.TFile('{o}/{r}.root'.format(o=self.runFilesDir, r=self.rawInputTreeNames[self.vcal]), mode)

	def OpenOutputFile(self, mode='READ'):
		if self.fileOutputRaw is not None:
			if self.fileOutputRaw.IsOpen:
				self.fileOutputRaw.Close()
		self.fileOutputRaw = ro.TFile('{o}/{r}.root'.format(o=self.runFilesDir, r=self.rawOutputTreeNames[self.vcal]), mode)

	def LoadInputTreeCSV(self):
		t0 = time.time()
		self.treeInputCSV = self.fileInputRaw.Get('treeInputCSV_pos_{b}mV'.format(b=abs(self.vcal*1000))) if self.vcal >= 0 else self.fileInputRaw.Get('treeInputCSV_neg_{b}mV'.format(b=abs(self.vcal*1000)))
		if not self.treeInputCSV:
			print 'The file does not contain the CSV tree... exiting. Delete the damaged file and run again'
			exit()
		t0 = time.time() - t0
		print 'Time loading CSV file:', t0, 'seconds'
		self.LoadInputWaveVectors(csvLoaded=True)
		self.boolInputHasTreeCSV[self.vcal] = True

	def LoadOutputTreeCSV(self):
		t0 = time.time()
		self.treeOutputCSV = self.fileOutputRaw.Get('treeOutputCSV_pos_{b}mV'.format(b=abs(self.vcal*1000))) if self.vcal >= 0 else self.fileOutputRaw.Get('treeOutputCSV_neg_{b}mV'.format(b=abs(self.vcal*1000)))
		if not self.treeOutputCSV:
			print 'The file does not contain the CSV tree... exiting. Delete the damaged file and run again'
			exit()
		t0 = time.time() - t0
		print 'Time loading CSV file:', t0, 'seconds'
		self.LoadOutputWaveVectors(csvLoaded=True)
		self.boolOutputHasTreeCSV[self.vcal] = True

	def CheckTreeRawInput(self):
		t0 = time.time()
		self.treeRawInput = self.fileInputRaw.Get(self.rawInputTreeNames[self.vcal])
		if not self.treeRawInput:
			print 'The raw tree has not been created. Creating...'
			self.CreateInputRawRootFile()
		else:
			t0 = time.time() - t0
			print 'Time loading raw tree:', t0, 'seconds'
			self.boolInputHasVectors[self.vcal] = True

	def CheckTreeRawOutput(self):
		t0 = time.time()
		self.treeRawOutput = self.fileOutputRaw.Get(self.rawOutputTreeNames[self.vcal])
		if not self.treeRawOutput:
			print 'The raw tree has not been created. Creating...'
			self.CreateOutputRawRootFile()
		else:
			t0 = time.time() - t0
			print 'Time loading raw tree:', t0, 'seconds'
			self.boolOutputHasVectors[self.vcal] = True

	def LoadInputTreeRaw(self, doLoadVectors=True):
		t0 = time.time()
		self.treeRawInput = self.fileInputRaw.Get(self.rawInputTreeNames[self.vcal])
		if not self.treeRawInput:
			print 'The raw tree has not been created. Creating...'
			self.CreateInputRawRootFile()
		else:
			t0 = time.time() - t0
			print 'Time loading raw tree:', t0, 'seconds'
			if doLoadVectors:
				self.LoadInputWaveVectorsFromROOT()
			self.boolInputHasVectors[self.vcal] = True
			self.CalculateAcceptedEventsInput()
			
	def LoadOutputTreeRaw(self, doLoadVectors=True):
		t0 = time.time()
		self.treeRawOutput = self.fileOutputRaw.Get(self.rawOutputTreeNames[self.vcal])
		if not self.treeRawOutput:
			print 'The raw tree has not been created. Creating...'
			self.CreateOutputRawRootFile()
		else:
			t0 = time.time() - t0
			print 'Time loading raw tree:', t0, 'seconds'
			if doLoadVectors:
				self.LoadOutputWaveVectorsFromROOT()
			self.boolOutputHasVectors[self.vcal] = True
			self.CalculateAcceptedEventsOutput()

	def ConvertFiles(self):
		self.ConvertInputFile()
		self.ConvertOutputFile()

	def ConvertInputFile(self):
		self.OpenInputFile(mode='NEW')
		if self.fileInputRaw.IsOpen():
			print 'The file has not been created. Converting csv file to ROOT file...'
			nameCSV = 'treeInputCSV_pos_{b}mV'.format(b=abs(self.vcal*1000)) if self.vcal >= 0 else 'treeInputCSV_neg_{b}mV'.format(b=abs(self.vcal*1000))
			self.treeInputCSV = ro.TTree(nameCSV, nameCSV)
			t0 = time.time()
			self.treeInputCSV.ReadFile('{d}/{f}'.format(d=self.runFilesDir, f=self.runInputFiles[self.vcal]), 'event/I:time/D:voltage1/F:voltage2/D')
			t0 = time.time() - t0
			self.boolInputHasTreeCSV[self.vcal] = True
			print 'Time reading csv file:', t0, 'seconds'
			# self.LoadInputWaveVectors(csvLoaded=True)
			self.fileInputRaw.Write()
			self.fileInputRaw.Close()
			self.CreateInputRawRootFile()
			self.boolInputConverted[self.vcal] = True
		else:
			print 'The file already exists. Loading existing file...'
			self.OpenInputFile(mode='UPDATE')
			self.CheckInputFileForScalars()
			if not self.boolInputHasScalars[self.vcal]:
				self.CheckTreeRawInput()

	def ConvertOutputFile(self):
		self.OpenOutputFile(mode='NEW')
		if self.fileOutputRaw.IsOpen():
			print 'The file has not been created. Converting csv file to ROOT file...'
			nameCSV = 'treeOutputCSV_pos_{b}mV'.format(b=abs(self.vcal*1000)) if self.vcal >= 0 else 'treeOutputCSV_neg_{b}mV'.format(b=abs(self.vcal*1000))
			self.treeOutputCSV = ro.TTree(nameCSV, nameCSV)
			t0 = time.time()
			self.treeOutputCSV.ReadFile('{d}/{f}'.format(d=self.runFilesDir, f=self.runOutputFiles[self.vcal]), 'event/I:time/D:voltage1/F:voltage2/D')
			t0 = time.time() - t0
			self.boolOutputHasTreeCSV[self.vcal] = True
			print 'Time reading csv file:', t0, 'seconds'
			# self.LoadOutputWaveVectors(csvLoaded=True)
			self.fileOutputRaw.Write()
			self.fileOutputRaw.Close()
			self.CreateOutputRawRootFile()
			self.boolOutputConverted[self.vcal] = True
		else:
			print 'The file already exists. Loading existing file...'
			self.OpenOutputFile(mode='UPDATE')
			self.CheckOutputFileForScalars()
			if not self.boolOutputHasScalars[self.vcal]:
				self.CheckTreeRawOutput()

	def LoadOutputWaveVectors(self, csvLoaded=True):
		if not csvLoaded:
			self.LoadOutputTreeCSV()
		t0 = time.time()
		self.treeOutputCSV.SetBranchStatus('*', 0)
		self.treeOutputCSV.SetBranchStatus('event', 1)
		entries = self.treeOutputCSV.GetEntries()
		self.treeOutputCSV.GetEntry(entries - 1)
		lastEvent = self.treeOutputCSV.event
		self.treeOutputCSV.GetEntry(0)
		firstEvent = self.treeOutputCSV.event
		self.eventsOutput = int(lastEvent - firstEvent + 1)
		self.ptsWaveOutput = int(entries / self.eventsOutput)
		self.treeOutputCSV.SetEstimate(entries)
		self.treeOutputCSV.SetBranchStatus('*', 0)
		self.treeOutputCSV.SetBranchStatus('time', 1)
		self.treeOutputCSV.SetBranchStatus('voltage2', 1)
		leng = self.treeOutputCSV.Draw('time:voltage2', '', 'goff')
		if leng > self.treeOutputCSV.GetEstimate():
			self.treeOutputCSV.SetEstimate(leng)
			leng = self.treeOutputCSV.Draw('time:voltage2', '', 'goff')
		tempTime = self.treeOutputCSV.GetVal(0)
		tempVolt = self.treeOutputCSV.GetVal(1)
		self.timeVectOutput = np.array([[tempTime[ev * self.ptsWaveOutput + t] for t in xrange(self.ptsWaveOutput)] for ev in xrange(self.eventsOutput)], 'f8')
		self.timeVectOutput = np.array(self.timeVectOutput[0], 'f8')
		self.wavesVectOutput = np.array([[tempVolt[ev * self.ptsWaveOutput + t] for t in xrange(self.ptsWaveOutput)] for ev in xrange(self.eventsOutput)], 'f8')
		self.wavesMeanVectOutput = self.wavesVectOutput.mean(axis=0)
		self.wavesMeanErrorsVectOutput = self.wavesVectOutput.std(axis=0)
		self.eventVectOutput = np.arange(0, self.eventsOutput, dtype='i')
		t0 = time.time() - t0
		self.treeOutputCSV.SetBranchStatus('*', 1)
		print 'Time loading vectors from CSV file:', t0, 'seconds'

	def LoadInputWaveVectors(self, csvLoaded=True):
		if not csvLoaded:
			self.LoadInputTreeCSV()
		t0 = time.time()
		self.treeInputCSV.SetBranchStatus('*', 0)
		self.treeInputCSV.SetBranchStatus('event', 1)
		entries = self.treeInputCSV.GetEntries()
		self.treeInputCSV.GetEntry(entries - 1)
		lastEvent = self.treeInputCSV.event
		self.treeInputCSV.GetEntry(0)
		firstEvent = self.treeInputCSV.event
		self.eventsInput = int(lastEvent - firstEvent + 1)
		self.ptsWaveInput = int(entries / self.eventsInput)
		self.treeInputCSV.SetEstimate(entries)
		self.treeInputCSV.SetBranchStatus('*', 0)
		self.treeInputCSV.SetBranchStatus('time', 1)
		self.treeInputCSV.SetBranchStatus('voltage2', 1)
		leng = self.treeInputCSV.Draw('time:voltage2', '', 'goff')
		if leng > self.treeInputCSV.GetEstimate():
			self.treeInputCSV.SetEstimate(leng)
			leng = self.treeInputCSV.Draw('time:voltage2', '', 'goff')
		tempTime = self.treeInputCSV.GetVal(0)
		tempVolt = self.treeInputCSV.GetVal(1)
		self.timeVectInput = np.array([[tempTime[ev * self.ptsWaveInput + t] for t in xrange(self.ptsWaveInput)] for ev in xrange(self.eventsInput)], 'f8')
		self.timeVectInput = np.array(self.timeVectInput[0], 'f8')
		self.wavesVectInput = np.array([[tempVolt[ev * self.ptsWaveInput + t] for t in xrange(self.ptsWaveInput)] for ev in xrange(self.eventsInput)], 'f8')
		self.wavesMeanVectInput = self.wavesVectInput.mean(axis=0)
		self.wavesMeanErrorsVectInput = self.wavesVectInput.std(axis=0)
		self.eventVectInput = np.arange(0, self.eventsInput, dtype='i')
		t0 = time.time() - t0
		self.treeInputCSV.SetBranchStatus('*', 1)
		print 'Time loading vectors from CSV file:', t0, 'seconds'

	def LoadInputWaveVectorsFromROOT(self):
		t0 = time.time()
		self.treeRawInput.SetBranchStatus('*', 0)
		self.treeRawInput.SetBranchStatus('event', 1)
		entries = self.treeRawInput.GetEntries()
		self.treeRawInput.GetEntry(entries - 1)
		lastEvent = self.treeRawInput.event
		self.treeRawInput.GetEntry(0)
		firstEvent = self.treeRawInput.event
		self.eventsInput = int(lastEvent - firstEvent + 1)
		self.ptsWaveInput = int(entries / self.eventsInput)
		self.treeRawInput.SetEstimate(1000000)
		self.treeRawInput.SetBranchStatus('*', 0)
		self.treeRawInput.SetBranchStatus('time', 1)
		self.ptsWaveInput = self.treeRawInput.Draw('time', '', 'goff', 1)
		tempTime = self.treeRawInput.GetVal(0)
		self.timeVectInput = np.array([tempTime[t] for t in xrange(self.ptsWaveInput)], 'f8')
		self.treeRawInput.SetBranchStatus('*', 0)
		self.treeRawInput.SetBranchStatus('voltageSignal', 1)
		leng = self.treeRawInput.Draw('voltageSignal', '', 'goff')
		if leng > self.treeRawInput.GetEstimate():
			self.treeRawInput.SetEstimate(leng)
			leng = self.treeRawInput.Draw('voltageSignal', '', 'goff')
		tempVolt = self.treeRawInput.GetVal(0)
		self.wavesVectInput = np.array([[tempVolt[ev * self.ptsWaveInput + t] for t in xrange(self.ptsWaveInput)] for ev in xrange(self.eventsInput)], 'f8')
		self.wavesMeanVectInput = self.wavesVectInput.mean(axis=0)
		self.wavesMeanErrorsVectInput = self.wavesVectInput.std(axis=0)
		self.eventVectInput = np.arange(0, self.eventsInput, dtype='i')
		t0 = time.time() - t0
		self.treeRawInput.SetBranchStatus('*', 1)
		print 'Time loading vectors from ROOT file:', t0, 'seconds'

	def LoadOutputWaveVectorsFromROOT(self):
		t0 = time.time()
		self.treeRawOutput.SetBranchStatus('*', 0)
		self.treeRawOutput.SetBranchStatus('event', 1)
		entries = self.treeRawOutput.GetEntries()
		self.treeRawOutput.GetEntry(entries - 1)
		lastEvent = self.treeRawOutput.event
		self.treeRawOutput.GetEntry(0)
		firstEvent = self.treeRawOutput.event
		self.eventsOutput = int(lastEvent - firstEvent + 1)
		self.ptsWaveOutput = int(entries / self.eventsOutput)
		self.treeRawOutput.SetEstimate(1000000)
		self.treeRawOutput.SetBranchStatus('*', 0)
		self.treeRawOutput.SetBranchStatus('time', 1)
		self.ptsWaveOutput = self.treeRawOutput.Draw('time', '', 'goff', 1)
		tempTime = self.treeRawOutput.GetVal(0)
		self.timeVectOutput = np.array([tempTime[t] for t in xrange(self.ptsWaveOutput)], 'f8')
		self.treeRawOutput.SetBranchStatus('*', 0)
		self.treeRawOutput.SetBranchStatus('voltageSignal', 1)
		leng = self.treeRawOutput.Draw('voltageSignal', '', 'goff')
		if leng > self.treeRawOutput.GetEstimate():
			self.treeRawOutput.SetEstimate(leng)
			leng = self.treeRawOutput.Draw('voltageSignal', '', 'goff')
		tempVolt = self.treeRawOutput.GetVal(0)
		self.wavesVectOutput = np.array([[tempVolt[ev * self.ptsWaveOutput + t] for t in xrange(self.ptsWaveOutput)] for ev in xrange(self.eventsOutput)], 'f8')
		self.wavesMeanVectOutput = self.wavesVectOutput.mean(axis=0)
		self.wavesMeanErrorsVectOutput = self.wavesVectOutput.std(axis=0)
		self.eventVectOutput = np.arange(0, self.eventsOutput, dtype='i')
		t0 = time.time() - t0
		self.treeRawOutput.SetBranchStatus('*', 1)
		print 'Time loading vectors from ROOT file:', t0, 'seconds'

	def CreateOutputRawRootFile(self):
		t0 = time.time()
		self.eventOutput = np.zeros(1, 'I')
		if self.fileOutputRaw.IsOpen():
			self.fileOutputRaw.Close()
		self.OpenOutputFile(mode='UPDATE')
		self.LoadOutputTreeCSV()
		self.timeBraOutput = np.zeros(self.ptsWaveOutput, 'f8')
		self.voltBraOutput = np.zeros(self.ptsWaveOutput, 'f8')
		self.treeRawOutput = ro.TTree(self.rawOutputTreeNames[self.vcal], self.rawOutputTreeNames[self.vcal])
		self.treeRawOutput.Branch('event', self.eventOutput, 'event/i')
		self.treeRawOutput.Branch('time', self.timeBraOutput, 'time[{s}]/D'.format(s=self.ptsWaveOutput))
		self.treeRawOutput.Branch('voltageSignal', self.voltBraOutput, 'voltageSignal[{s}]/D'.format(s=self.ptsWaveOutput))
		self.CreateProgressBar(self.eventsOutput)
		self.bar.start()
		for i in xrange(self.eventsOutput):
			self.eventOutput.fill(i)
			np.putmask(self.timeBraOutput, 1 - np.zeros(self.ptsWaveOutput, '?'), self.timeVectOutput)
			np.putmask(self.voltBraOutput, 1 - np.zeros(self.ptsWaveOutput, '?'), self.wavesVectOutput[i])
			numFil = self.treeRawOutput.Fill()
			self.bar.update(i+1)
		self.bar.finish()
		self.fileOutputRaw.Write()
		self.fileOutputRaw.Close()
		t0 = time.time() - t0
		print 'Time creating raw tree:', t0, 'seconds'
		self.boolOutputHasVectors[self.vcal] = True
		self.CalculateAcceptedEventsOutput()

	def CreateInputRawRootFile(self):
		t0 = time.time()
		self.eventInput = np.zeros(1, 'I')
		if self.fileInputRaw.IsOpen():
			self.fileInputRaw.Close()
		self.OpenInputFile(mode='UPDATE')
		self.LoadInputTreeCSV()
		self.timeBraInput = np.zeros(self.ptsWaveInput, 'f8')
		self.voltBraInput = np.zeros(self.ptsWaveInput, 'f8')
		self.treeRawInput = ro.TTree(self.rawInputTreeNames[self.vcal], self.rawInputTreeNames[self.vcal])
		self.treeRawInput.Branch('event', self.eventInput, 'event/i')
		self.treeRawInput.Branch('time', self.timeBraInput, 'time[{s}]/D'.format(s=self.ptsWaveInput))
		self.treeRawInput.Branch('voltageSignal', self.voltBraInput, 'voltageSignal[{s}]/D'.format(s=self.ptsWaveInput))
		self.CreateProgressBar(self.eventsInput)
		self.bar.start()
		for i in xrange(self.eventsInput):
			self.eventInput.fill(i)
			np.putmask(self.timeBraInput, 1 - np.zeros(self.ptsWaveInput, '?'), self.timeVectInput)
			np.putmask(self.voltBraInput, 1 - np.zeros(self.ptsWaveInput, '?'), self.wavesVectInput[i])
			numFil = self.treeRawInput.Fill()
			self.bar.update(i+1)
		self.bar.finish()
		self.fileInputRaw.Write()
		self.fileInputRaw.Close()
		t0 = time.time() - t0
		print 'Time creating raw tree:', t0, 'seconds'
		self.boolInputHasVectors[self.vcal] = True
		self.CalculateAcceptedEventsInput()

	def LoadExistingFile(self):
		self.LoadOutputTreeRaw()

	def SNRMap(self):
		print 'Creating SNR Map on the current voltage in self.vcal'
		t0 = time.time()
		if self.fileOutputRaw.IsOpen():
			self.fileOutputRaw.Close()
		self.OpenOutputFile('READ')
		self.LoadOutputTreeRaw()
		name = 'snr_pos_{v}mV'.format(v=self.vcal*1000) if self.vcal >= 0 else 'snr_neg_{v}mV'.format(v=abs(self.vcal*1000))
		nameline = 'snrline_pos_{v}mV'.format(v=self.vcal*1000) if self.vcal >= 0 else 'snrline_neg_{v}mV'.format(v=abs(self.vcal*1000))
		fileSNR = ro.TFile('{n}.root'.format(n=name), 'RECREATE')
		self.fileOutputRaw.cd()
		minbl = -0.1e-6
		maxbl = 1e-6
		binsl = int((maxbl-minbl)/0.01e-6 + 1)
		minbr = -0.1e-6
		maxbr = 1e-6
		binsr = int((maxbr-minbr)/0.01e-6 + 1)
		left = np.linspace(minbl, maxbl, binsl, dtype='f8')
		right = np.linspace(minbr, maxbr, binsr, dtype='f8')
		snrMap = ro.TH2D(name+'_map', name+'_map', binsl, minbl-(maxbl-minbl)/(2*(binsl-1)), maxbl+(maxbl-minbl)/(2*(binsl-1)), binsr, minbr-(maxbr-minbr)/(2*(binsr-1)), maxbr+(maxbr-minbr)/(2*(binsr-1)))
		lineValues = [2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98]
		self.snrLines = {val: ro.TH1D('{nl}_{v}'.format(nl=nameline, v=val), '{nl}_{v}'.format(nl=nameline, v=val), binsl, minbl-(maxbl-minbl)/(2*(binsl-1)), maxbl+(maxbl-minbl)/(2*(binsl-1))) for val in lineValues}
		self.CreateProgressBar(len(left)*len(right))
		self.bar.start()
		for i in xrange(len(left)):
			for j in xrange(len(right)):
				self.averagingTime = left[i]+right[j]
				if 1e-6 + self.pedestalTEndPos >= self.averagingTime > 1e-10:
					self.SNRCalculation(left[i], right[j])
					if self.pedSigma[self.vcal] != 0:
						snrMap.SetBinContent(i+1, j+1, abs(float(self.signal[self.vcal])/float(self.pedSigma[self.vcal])))
						for val in lineValues:
							if abs(self.averagingTime - val*1e-8) <= 5e-9:
								self.snrLines[val].SetBinContent(i+1, abs(float(self.signal[self.vcal])/float(self.pedSigma[self.vcal])))
					else:
						snrMap.SetBinContent(i+1, j+1, 0)
				else:
					snrMap.SetBinContent(i + 1, j + 1, 0)
				self.bar.update(i*len(right) + j + 1)
		self.bar.finish()
		fileSNR.cd()
		snrMap.Write()
		for val in lineValues:
			self.snrLines[val].Write()
		fileSNR.Write()
		fileSNR.Close()
		t0 = time.time() - t0
		print 'Time creating the SNR analysis and plots:', t0, 'seconds'
		self.fileOutputRaw.cd()
		self.fileOutputRaw.Close()

	def SNRCalculation(self, left, right):
		self.AnalysisAllOutputWaves(False, left, right)

	def AnalysisAllInputWaves(self, saveFile=True):
		self.pedInput1Vect, self.pedInput2Vect, self.pedInput3Vect, self.pedInput4Vect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')
		self.pedInputSigma1Vect, self.pedInputSigma2Vect, self.pedInputSigma3Vect, self.pedInputSigma4Vect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')
		self.pedCal1Vect, self.pedCal2Vect, self.pedCal3Vect, self.pedCal4Vect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')
		self.calVoltsReal1Vect, self.calVoltsReal2Vect, self.calVoltsReal3Vect, self.calVoltsReal4Vect = np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8'), np.empty(0, 'f8')

		self.pedInput1[self.vcal], self.pedInput2[self.vcal], self.pedInput3[self.vcal], self.pedInput4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.pedInputSigma1[self.vcal], self.pedInputSigma2[self.vcal], self.pedInputSigma3[self.vcal], self.pedInputSigma4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.pedCal1[self.vcal], self.pedCal2[self.vcal], self.pedCal3[self.vcal], self.pedCal4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.calVoltsReal1[self.vcal], self.calVoltsReal2[self.vcal], self.calVoltsReal3[self.vcal], self.calVoltsReal4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.calVoltsRealSigma1[self.vcal], self.calVoltsRealSigma2[self.vcal], self.calVoltsRealSigma3[self.vcal], self.calVoltsRealSigma4[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8'), np.zeros(1, 'f8')

		self.posHalfCal = self.FindTimeHalfCalPosition()
		(self.pedInput1Pos, self.pedInput2Pos, self.pedInput3Pos, self.pedInput4Pos) = self.FindInputPedestalPositions()
		(self.pedCal1Pos, self.pedCal2Pos, self.pedCal3Pos, self.pedCal4Pos) = self.FindInputPedCalPositions()
		self.CalculateRealCalibratesAllWaves(saveFile)
		self.boolInputHasScalars[self.vcal] = True
		
	def FindTimeHalfCalPosition(self):
		if self.meanInputWaveGraph is None:
			self.CreateMeanInputWaveformGraph()
		left, right = np.double(self.wavesMeanVectInput[0]), np.double(self.wavesMeanVectInput[-1])
		mid = np.double((left+right)/2.0)
		distFromMid = np.array(np.abs(self.wavesMeanVectInput - mid), 'f8')
		midPos = distFromMid.argmin()
		return midPos
	
	def FindInputPedestalPositions(self):
		pos1 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] + self.calDist1), 'f8').argmin()
		pos2 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] + self.calDist2), 'f8').argmin()
		pos3 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] + self.calDist3), 'f8').argmin()
		pos4 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] + self.calDist4), 'f8').argmin()
		return pos1, pos2, pos3, pos4

	def FindInputPedCalPositions(self):
		pos1 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] - self.calDist1), 'f8').argmin()
		pos2 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] - self.calDist2), 'f8').argmin()
		pos3 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] - self.calDist3), 'f8').argmin()
		pos4 = np.array(np.abs(self.timeVectInput - self.timeVectInput[self.posHalfCal] - self.calDist4), 'f8').argmin()
		return pos1, pos2, pos3, pos4

	def AnalysisAllOutputWaves(self, saveFile=True, left=0.469e-6, right=0.511e-6):
		self.pedVect = np.empty(0, 'f8')
		self.pedSigmaVect = np.empty(0, 'f8')
		self.pedSignalVect = np.empty(0, 'f8')
		self.pedSignalRealVect = np.empty(0, 'f8')
		self.signalVect, self.signalRealVect = np.empty(0, 'f8'), np.empty(0, 'f8')
		self.ped[self.vcal], self.pedSigma[self.vcal] = np.zeros(1, 'f8'), np.zeros(1, 'f8')
		self.pedSignal[self.vcal] = np.zeros(1, 'f8')
		self.pedSignalReal[self.vcal] = np.zeros(1, 'f8')
		self.signal[self.vcal] = np.zeros(1, 'f8')
		self.signalReal[self.vcal] = np.zeros(1, 'f8')

		self.peakBack, self.peakForth = left, right
		self.timePeakReal = self.FindRealPeakPosition()
		self.FindOutputPedestalPositions(self.peakBack, self.peakForth)
		self.FindSignalPosition(self.peakBack, self.peakForth)
		self.FindSignalPositionReal(self.peakBack, self.peakForth)
		self.CalculateSignalsAllWaves(saveFile)
		self.boolOutputHasScalars[self.vcal] = True

	def FindRealPeakPosition(self):
		if self.meanOutputWaveGraph is None:
			self.CreateMeanOutputWaveformGraph()
		mpos = self.wavesMeanVectOutput.argmin() if self.vcal >= 0 else self.wavesMeanVectOutput.argmax()
		xmin, xmax = self.timeVectOutput[mpos] - 200e-9, self.timeVectOutput[mpos] + 200e-9
		fit = self.meanOutputWaveGraph.Fit('pol2', 'QMN0FS', '', xmin, xmax)
		c, b, a = fit.Parameter(0), fit.Parameter(1), fit.Parameter(2)
		return -b/(2.0*a)

	def CreateMeanInputWaveformGraph(self):
		if self.fileInputRaw.GetOption() == 'READ':
			self.fileInputRaw.ReOpen('UPDATE')
		self.meanInputWaveGraph = ro.TGraphErrors(self.ptsWaveInput, self.timeVectInput, self.wavesMeanVectInput, np.zeros(self.ptsWaveInput, 'f8'), self.wavesMeanErrorsVectInput)
		namegraph = 'Mean_input_wave_vcal_pos_{v}mV'.format(v=1000*self.vcal) if self.vcal >= 0 else 'Mean_input_wave_vcal_neg_{v}mV'.format(v=abs(1000*self.vcal))
		self.meanInputWaveGraph.SetNameTitle(namegraph, namegraph)
		self.meanInputWaveGraph.GetXaxis().SetTitle('time/s')
		self.meanInputWaveGraph.GetYaxis().SetTitle('voltage/V')
		self.fileInputRaw.cd()
		self.meanInputWaveGraph.Write()
		self.fileInputRaw.Write()

	def CreateMeanOutputWaveformGraph(self):
		if self.fileOutputRaw.GetOption() == 'READ':
			self.fileOutputRaw.ReOpen('UPDATE')
		self.meanOutputWaveGraph = ro.TGraphErrors(self.ptsWaveOutput, self.timeVectOutput, self.wavesMeanVectOutput, np.zeros(self.ptsWaveOutput, 'f8'), self.wavesMeanErrorsVectOutput)
		namegraph = 'Mean_output_wave_vcal_pos_{v}mV'.format(v=1000*self.vcal) if self.vcal >= 0 else 'Mean_output_wave_vcal_neg_{v}mV'.format(v=abs(1000*self.vcal))
		self.meanOutputWaveGraph.SetNameTitle(namegraph, namegraph)
		self.meanOutputWaveGraph.GetXaxis().SetTitle('time/s')
		self.meanOutputWaveGraph.GetYaxis().SetTitle('voltage/V')
		self.fileOutputRaw.cd()
		self.meanOutputWaveGraph.Write()
		self.fileOutputRaw.Write()

	def FindOutputPedestalPositions(self, left, right):
		self.averagingTime = left + right
		self.pedestalTimeIndices = np.argwhere(np.bitwise_and(self.pedestalTEndPos - self.averagingTime <= self.timeVectOutput, self.timeVectOutput <= self.pedestalTEndPos)).flatten()

	def FindSignalPosition(self, left, right):
		self.averagingTime = left + right
		self.signalTimeIndices = np.argwhere(abs(self.timeVectOutput - self.timePeak + (right - left)/2.0) <= self.averagingTime/2.0).flatten()

	def FindSignalPositionReal(self, left, right):
		self.averagingTime = left + right
		self.signalTimeIndicesReal = np.argwhere(abs(self.timeVectOutput - self.timePeakReal + (right - left)/2.0) <= self.averagingTime/2.0).flatten()
		
	def CalculateRealCalibratesAllWaves(self, saveFile=True):
		if saveFile:
			t0 = time.time()
			print 'Calculating real calibration input signals...'
			self.CreateProgressBar(self.eventsInput)
			self.bar.start()
		for ev in self.eventVectInput:
			if self.acceptEventsInput[ev]:
				evPed1 = self.wavesVectInput[ev, self.pedInput1Pos]
				evPed2 = self.wavesVectInput[ev, self.pedInput2Pos]
				evPed3 = self.wavesVectInput[ev, self.pedInput3Pos]
				evPed4 = self.wavesVectInput[ev, self.pedInput4Pos]
				self.pedInput1Vect = np.insert(self.pedInput1Vect, len(self.pedInput1Vect), evPed1)
				self.pedInput2Vect = np.insert(self.pedInput2Vect, len(self.pedInput2Vect), evPed2)
				self.pedInput3Vect = np.insert(self.pedInput3Vect, len(self.pedInput3Vect), evPed3)
				self.pedInput4Vect = np.insert(self.pedInput4Vect, len(self.pedInput4Vect), evPed4)
				evPedCal1 = self.wavesVectInput[ev, self.pedCal1Pos]
				evPedCal2 = self.wavesVectInput[ev, self.pedCal2Pos]
				evPedCal3 = self.wavesVectInput[ev, self.pedCal3Pos]
				evPedCal4 = self.wavesVectInput[ev, self.pedCal4Pos]
				self.pedCal1Vect = np.insert(self.pedCal1Vect, len(self.pedCal1Vect), evPedCal1)
				self.pedCal2Vect = np.insert(self.pedCal2Vect, len(self.pedCal2Vect), evPedCal2)
				self.pedCal3Vect = np.insert(self.pedCal3Vect, len(self.pedCal3Vect), evPedCal3)
				self.pedCal4Vect = np.insert(self.pedCal4Vect, len(self.pedCal4Vect), evPedCal4)
				evCal1 = evPedCal1 - evPed1
				evCal2 = evPedCal2 - evPed2
				evCal3 = evPedCal3 - evPed3
				evCal4 = evPedCal4 - evPed4
				self.calVoltsReal1Vect = np.insert(self.calVoltsReal1Vect, len(self.calVoltsReal1Vect), evCal1)
				self.calVoltsReal2Vect = np.insert(self.calVoltsReal2Vect, len(self.calVoltsReal2Vect), evCal2)
				self.calVoltsReal3Vect = np.insert(self.calVoltsReal3Vect, len(self.calVoltsReal3Vect), evCal3)
				self.calVoltsReal4Vect = np.insert(self.calVoltsReal4Vect, len(self.calVoltsReal4Vect), evCal4)
			else:
				self.pedInput1Vect = np.insert(self.pedInput1Vect, len(self.pedInput1Vect), -10)
				self.pedInput2Vect = np.insert(self.pedInput2Vect, len(self.pedInput2Vect), -10)
				self.pedInput3Vect = np.insert(self.pedInput3Vect, len(self.pedInput3Vect), -10)
				self.pedInput4Vect = np.insert(self.pedInput4Vect, len(self.pedInput4Vect), -10)
				self.pedCal1Vect = np.insert(self.pedCal1Vect, len(self.pedCal1Vect), -10)
				self.pedCal2Vect = np.insert(self.pedCal2Vect, len(self.pedCal2Vect), -10)
				self.pedCal3Vect = np.insert(self.pedCal3Vect, len(self.pedCal3Vect), -10)
				self.pedCal4Vect = np.insert(self.pedCal4Vect, len(self.pedCal4Vect), -10)
				self.calVoltsReal1Vect = np.insert(self.calVoltsReal1Vect, len(self.calVoltsReal1Vect), -10)
				self.calVoltsReal2Vect = np.insert(self.calVoltsReal2Vect, len(self.calVoltsReal2Vect), -10)
				self.calVoltsReal3Vect = np.insert(self.calVoltsReal3Vect, len(self.calVoltsReal3Vect), -10)
				self.calVoltsReal4Vect = np.insert(self.calVoltsReal4Vect, len(self.calVoltsReal4Vect), -10)
			if saveFile:
				self.bar.update(ev + 1)
		if saveFile:
			self.bar.finish()
			t0 = time.time() - t0
			print 'Time Calculating all real calibration input voltages:', t0, 'seconds'
			print 'Saving real calibration input voltages in file...'
			t0 = time.time()
			if self.fileInputRaw.IsOpen():
				self.fileInputRaw.Close()
			self.OpenInputFile('UPDATE')
			self.LoadInputTreeRaw()
			ped1Branch = self.treeRawInput.Branch('pedestal1', self.pedInput1[self.vcal], 'pedestal1/D')
			ped2Branch = self.treeRawInput.Branch('pedestal2', self.pedInput2[self.vcal], 'pedestal2/D')
			ped3Branch = self.treeRawInput.Branch('pedestal3', self.pedInput3[self.vcal], 'pedestal3/D')
			ped4Branch = self.treeRawInput.Branch('pedestal4', self.pedInput4[self.vcal], 'pedestal4/D')
			pedCal1Branch = self.treeRawInput.Branch('pedestalAndCalV1', self.pedCal1[self.vcal], 'pedestalAndCalV1/D')
			pedCal2Branch = self.treeRawInput.Branch('pedestalAndCalV2', self.pedCal2[self.vcal], 'pedestalAndCalV2/D')
			pedCal3Branch = self.treeRawInput.Branch('pedestalAndCalV3', self.pedCal3[self.vcal], 'pedestalAndCalV3/D')
			pedCal4Branch = self.treeRawInput.Branch('pedestalAndCalV4', self.pedCal4[self.vcal], 'pedestalAndCalV4/D')
			calV1Branch = self.treeRawInput.Branch('calVolt1', self.calVoltsReal1[self.vcal], 'calVolt1/D')
			calV2Branch = self.treeRawInput.Branch('calVolt2', self.calVoltsReal2[self.vcal], 'calVolt2/D')
			calV3Branch = self.treeRawInput.Branch('calVolt3', self.calVoltsReal3[self.vcal], 'calVolt3/D')
			calV4Branch = self.treeRawInput.Branch('calVolt4', self.calVoltsReal4[self.vcal], 'calVolt4/D')
			self.treeRawInput.SetBranchStatus('pedestal1', 1)
			self.treeRawInput.SetBranchStatus('pedestal2', 1)
			self.treeRawInput.SetBranchStatus('pedestal3', 1)
			self.treeRawInput.SetBranchStatus('pedestal4', 1)
			self.treeRawInput.SetBranchStatus('pedestalAndCalV1', 1)
			self.treeRawInput.SetBranchStatus('pedestalAndCalV2', 1)
			self.treeRawInput.SetBranchStatus('pedestalAndCalV3', 1)
			self.treeRawInput.SetBranchStatus('pedestalAndCalV4', 1)
			self.treeRawInput.SetBranchStatus('calVolt1', 1)
			self.treeRawInput.SetBranchStatus('calVolt2', 1)
			self.treeRawInput.SetBranchStatus('calVolt3', 1)
			self.treeRawInput.SetBranchStatus('calVolt4', 1)
			self.CreateProgressBar(self.eventsInput)
			self.bar.start()
			for ev in self.eventVectInput:
				self.treeRawInput.GetEntry(ev)
				self.pedInput1[self.vcal].fill(self.pedInput1Vect[ev])
				self.pedInput2[self.vcal].fill(self.pedInput2Vect[ev])
				self.pedInput3[self.vcal].fill(self.pedInput3Vect[ev])
				self.pedInput4[self.vcal].fill(self.pedInput4Vect[ev])
				self.pedCal1[self.vcal].fill(self.pedCal1Vect[ev])
				self.pedCal2[self.vcal].fill(self.pedCal2Vect[ev])
				self.pedCal3[self.vcal].fill(self.pedCal3Vect[ev])
				self.pedCal4[self.vcal].fill(self.pedCal4Vect[ev])
				self.calVoltsReal1[self.vcal].fill(self.calVoltsReal1Vect[ev])
				self.calVoltsReal2[self.vcal].fill(self.calVoltsReal2Vect[ev])
				self.calVoltsReal3[self.vcal].fill(self.calVoltsReal3Vect[ev])
				self.calVoltsReal4[self.vcal].fill(self.calVoltsReal4Vect[ev])
				ped1Branch.Fill()
				ped2Branch.Fill()
				ped3Branch.Fill()
				ped4Branch.Fill()
				pedCal1Branch.Fill()
				pedCal2Branch.Fill()
				pedCal3Branch.Fill()
				pedCal4Branch.Fill()
				calV1Branch.Fill()
				calV2Branch.Fill()
				calV3Branch.Fill()
				calV4Branch.Fill()
				self.bar.update(ev + 1)
			self.bar.finish()
			self.treeRawInput.Write()
			# self.treeRawInput.Delete()
			self.fileInputRaw.Close()
			t0 = time.time() - t0
			print 'Time saving real calibration input voltages branches:', t0, 'seconds'
		tempPed1 = np.delete(self.pedInput1Vect, np.where(self.pedInput1Vect == -10))
		tempPed2 = np.delete(self.pedInput2Vect, np.where(self.pedInput2Vect == -10))
		tempPed3 = np.delete(self.pedInput3Vect, np.where(self.pedInput3Vect == -10))
		tempPed4 = np.delete(self.pedInput4Vect, np.where(self.pedInput4Vect == -10))
		tempPedCal1 = np.delete(self.pedCal1Vect, np.where(self.pedCal1Vect == -10))
		tempPedCal2 = np.delete(self.pedCal2Vect, np.where(self.pedCal2Vect == -10))
		tempPedCal3 = np.delete(self.pedCal3Vect, np.where(self.pedCal3Vect == -10))
		tempPedCal4 = np.delete(self.pedCal4Vect, np.where(self.pedCal4Vect == -10))
		tempCal1 = np.delete(self.calVoltsReal1Vect, np.where(self.calVoltsReal1Vect == -10))
		tempCal2 = np.delete(self.calVoltsReal2Vect, np.where(self.calVoltsReal2Vect == -10))
		tempCal3 = np.delete(self.calVoltsReal3Vect, np.where(self.calVoltsReal3Vect == -10))
		tempCal4 = np.delete(self.calVoltsReal4Vect, np.where(self.calVoltsReal4Vect == -10))
		self.pedInput1[self.vcal].fill(tempPed1.mean())
		self.pedInputSigma1[self.vcal].fill(tempPed1.std())
		self.pedCal1[self.vcal].fill(tempPedCal1.mean())
		self.calVoltsReal1[self.vcal].fill(tempCal1.mean())
		self.calVoltsRealSigma1[self.vcal].fill(tempCal1.std())
		self.pedInput2[self.vcal].fill(tempPed2.mean())
		self.pedInputSigma2[self.vcal].fill(tempPed2.std())
		self.pedCal2[self.vcal].fill(tempPedCal2.mean())
		self.calVoltsReal2[self.vcal].fill(tempCal2.mean())
		self.calVoltsRealSigma2[self.vcal].fill(tempCal2.std())
		self.pedInput3[self.vcal].fill(tempPed3.mean())
		self.pedInputSigma3[self.vcal].fill(tempPed3.std())
		self.pedCal3[self.vcal].fill(tempPedCal3.mean())
		self.calVoltsReal3[self.vcal].fill(tempCal3.mean())
		self.calVoltsRealSigma3[self.vcal].fill(tempCal3.std())
		self.pedInput4[self.vcal].fill(tempPed4.mean())
		self.pedInputSigma4[self.vcal].fill(tempPed4.std())
		self.pedCal4[self.vcal].fill(tempPedCal4.mean())
		self.calVoltsReal4[self.vcal].fill(tempCal4.mean())
		self.calVoltsRealSigma4[self.vcal].fill(tempCal4.std())

	def CalculateSignalsAllWaves(self, saveFile=True):
		if saveFile:
			t0 = time.time()
			print 'Calculating pedestals and signals...'
			self.CreateProgressBar(self.eventsOutput)
			self.bar.start()
		for ev in self.eventVectOutput:
			if self.acceptEventsOutput[ev]:
				evPedMean = self.wavesVectOutput[ev, self.pedestalTimeIndices].mean() if self.pedestalTimeIndices.size >= 1 else -10
				evPedSigma = self.wavesVectOutput[ev, self.pedestalTimeIndices].std() if self.pedestalTimeIndices.size >= 2 else -10
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), evPedMean)
				self.pedSigmaVect = np.insert(self.pedSigmaVect, len(self.pedSigmaVect), evPedSigma)
				evMeasMean = self.wavesVectOutput[ev, self.signalTimeIndices].mean() if self.signalTimeIndices.size >= 1 else -10
				evMeasMeanReal = self.wavesVectOutput[ev, self.signalTimeIndicesReal].mean() if self.signalTimeIndicesReal.size >= 1 else -10
				evSignalMean = evMeasMean - evPedMean
				evSignalMeanReal = evMeasMeanReal - evPedMean
				self.pedSignalVect = np.insert(self.pedSignalVect, len(self.pedSignalVect), evMeasMean)
				self.pedSignalRealVect = np.insert(self.pedSignalRealVect, len(self.pedSignalRealVect), evMeasMeanReal)
				self.signalVect = np.insert(self.signalVect, len(self.signalVect), evSignalMean)
				self.signalRealVect = np.insert(self.signalRealVect, len(self.signalRealVect), evSignalMeanReal)
			else:
				self.pedVect = np.insert(self.pedVect, len(self.pedVect), -10)
				self.pedSigmaVect = np.insert(self.pedSigmaVect, len(self.pedSigmaVect), -10)
				self.pedSignalVect = np.insert(self.pedSignalVect, len(self.pedSignalVect), -10)
				self.pedSignalRealVect = np.insert(self.pedSignalRealVect, len(self.pedSignalRealVect), -10)
				self.signalVect = np.insert(self.signalVect, len(self.signalVect), -10)
				self.signalRealVect = np.insert(self.signalRealVect, len(self.signalRealVect), -10)
			if saveFile:
				self.bar.update(ev + 1)
		if saveFile:
			self.bar.finish()
			t0 = time.time() - t0
			print 'Time Calculating all pedestals and signals:', t0, 'seconds'
			print 'Saving pedestals and signals in file...'
			t0 = time.time()
			if self.fileOutputRaw.IsOpen():
				self.fileOutputRaw.Close()
			self.OpenOutputFile('UPDATE')
			self.LoadOutputTreeRaw()
			pedBranch = self.treeRawOutput.Branch('pedestal', self.ped[self.vcal], 'pedestal/D')
			pedSigmaBranch = self.treeRawOutput.Branch('pedestalSigma', self.pedSigma[self.vcal], 'pedestalSigma/D')
			sigPedBranch = self.treeRawOutput.Branch('pedestalAndSignal', self.pedSignal[self.vcal], 'pedestalAndSignal/D')
			sigPedRealBranch = self.treeRawOutput.Branch('pedestalAndSignalReal', self.pedSignalReal[self.vcal], 'pedestalAndSignalReal/D')
			sigBranch = self.treeRawOutput.Branch('signal', self.signal[self.vcal], 'signal/D')
			sigRealBranch = self.treeRawOutput.Branch('signalReal', self.signalReal[self.vcal], 'signalReal/D')
			self.treeRawOutput.SetBranchStatus('pedestal', 1)
			self.treeRawOutput.SetBranchStatus('pedestalSigma', 1)
			self.treeRawOutput.SetBranchStatus('pedestalAndSignal', 1)
			self.treeRawOutput.SetBranchStatus('pedestalAndSignalReal', 1)
			self.treeRawOutput.SetBranchStatus('signal', 1)
			self.treeRawOutput.SetBranchStatus('signalReal', 1)
			self.CreateProgressBar(self.eventsOutput)
			self.bar.start()
			for ev in self.eventVectOutput:
				self.treeRawOutput.GetEntry(ev)
				self.ped[self.vcal].fill(self.pedVect[ev])
				self.pedSigma[self.vcal].fill(self.pedSigmaVect[ev])
				self.pedSignal[self.vcal].fill(self.pedSignalVect[ev])
				self.pedSignalReal[self.vcal].fill(self.pedSignalRealVect[ev])
				self.signal[self.vcal].fill(self.signalVect[ev])
				self.signalReal[self.vcal].fill(self.signalRealVect[ev])
				pedBranch.Fill()
				pedSigmaBranch.Fill()
				sigPedBranch.Fill()
				sigPedRealBranch.Fill()
				sigBranch.Fill()
				sigRealBranch.Fill()
				self.bar.update(ev + 1)
			self.bar.finish()
			self.treeRawOutput.Write()
			# self.treeRawOutput.Delete()
			self.fileOutputRaw.Close()
			t0 = time.time() - t0
			print 'Time saving pedestal and signal branches:', t0, 'seconds'
		tempPed = np.delete(self.pedVect, np.where(self.pedVect == -10))
		tempSAP = np.delete(self.pedSignalVect, np.where(self.pedSignalVect == -10))
		tempSAPR = np.delete(self.pedSignalRealVect, np.where(self.pedSignalRealVect == -10))
		tempS = np.delete(self.signalVect, np.where(self.signalVect == -10))
		tempSR = np.delete(self.signalRealVect, np.where(self.signalRealVect == -10))
		if tempPed.size > 1 and tempSAP.size > 0 and tempS.size > 0 and tempSAPR.size > 0 and tempSR.size > 0:
			self.ped[self.vcal].fill(tempPed.mean())
			self.pedSigma[self.vcal].fill(tempPed.std())
			self.pedSignal[self.vcal].fill(tempSAP.mean())
			self.pedSignalReal[self.vcal].fill(tempSAPR.mean())
			self.signal[self.vcal].fill(tempS.mean())
			self.signalReal[self.vcal].fill(tempSR.mean())
		else:
			self.ped[self.vcal].fill(0)
			self.pedSigma[self.vcal].fill(0)
			self.pedSignal[self.vcal].fill(0)
			self.pedSignalReal[self.vcal].fill(0)
			self.signal[self.vcal].fill(0)
			self.signalReal[self.vcal].fill(0)

	def AnalysisMeanWaves(self):
		self.timePeakReal = self.FindRealPeakPosition()
		self.FindOutputPedestalPositions(self.peakBack, self.peakForth)
		self.FindSignalPosition(self.peakBack, self.peakForth)
		self.FindSignalPositionReal(self.peakBack, self.peakForth)
		self.pedestalMeanWaveAverage = self.wavesMeanVectOutput[self.pedestalTimeIndices].mean()
		self.signalAndPedestalMeanWaveAverage = self.wavesMeanVectOutput[self.signalTimeIndices].mean()
		self.signalMeanWaveAverage = self.signalAndPedestalMeanWaveAverage - self.pedestalMeanWaveAverage
		self.signalAndPedestalMeanWaveAverageReal = self.wavesMeanVectOutput[self.signalTimeIndicesReal].mean()
		self.signalMeanWaveAverageReal = self.signalAndPedestalMeanWaveAverageReal - self.pedestalMeanWaveAverage

	def AnalyseAllCalibrationFiles(self, realCal=4):
		t0 = time.time()
		for vcal in self.calVolts:
			self.vcal = vcal
			self.AnalysisCalibrationFiles(realCal)
		t0 = time.time() - t0
		print 'Time to get all the values to do the regression:', t0, 'seconds'

	def AnalysisCalibrationFiles(self, realCal=4):
		self.AnalysisCalibrationInputFile(realCal)
		self.AnalysisCalibrationOutputFile()

	def AnalysisCalibrationInputFile(self, realCal=4):
		t0 = time.time()
		print 'Analysing {vc}mV calibration signal...'.format(vc=self.vcal*1000)
		self.OpenInputFile('UPDATE')
		if not self.boolInputHasScalars[self.vcal]:
			self.LoadInputTreeRaw()
			self.AnalysisAllInputWaves(True)
		else:
			self.LoadInputScalars()
		self.CalculateChargeFromVcal(realCal)
		t0 = time.time() - t0
		print 'Time getting mean calibration input voltage and noise for {vc}mV calibration: {t} seconds'.format(vc=self.vcal, t=t0)

	def AnalysisCalibrationOutputFile(self):
		t0 = time.time()
		print 'Analysing {vc}mV calibration signal...'.format(vc=self.vcal*1000)
		self.OpenOutputFile('UPDATE')
		if not self.boolOutputHasScalars[self.vcal]:
			self.LoadOutputTreeRaw()
			self.AnalysisAllOutputWaves(True, self.peakBack, self.peakForth)
		else:
			self.LoadOutputScalars()
		t0 = time.time() - t0
		print 'Time getting mean signal and noise for {vc}mV calibration: {t} seconds'.format(vc=self.vcal, t=t0)

	def CheckInputFileForScalars(self):
		t0 = time.time()
		print 'Checking if input file has been analysed...'
		if not self.fileInputRaw.IsOpen():
			self.OpenInputFile('UPDATE')
		self.LoadInputTreeRaw(False)
		if not self.treeRawInput.GetBranch('pedestal1'):
			print 'The input raw tree does not have the scalar branches'
		else:
			self.boolInputHasScalars[self.vcal] = True
			print 'The input raw tree has the scalar branches'
		t0 = time.time() - t0
		print 'Time checking for scalar branches in the input tree:', t0, 'seconds'

	def CheckOutputFileForScalars(self):
		t0 = time.time()
		print 'Checking if output file has been analysed...'
		if not self.fileOutputRaw.IsOpen():
			self.OpenOutputFile('UPDATE')
		self.LoadOutputTreeRaw(False)
		if not self.treeRawOutput.GetBranch('pedestal'):
			print 'The output raw tree does not have the scalar branches'
		else:
			self.boolOutputHasScalars[self.vcal] = True
			print 'The output raw tree has the scalar branches'
		t0 = time.time() - t0
		print 'Time checking for scalar branches in the output tree:', t0, 'seconds'

	def LoadOutputScalars(self):
		self.LoadOutputTreeRaw(False)
		self.treeRawOutput.SetBranchStatus('*', 0)
		self.treeRawOutput.SetBranchStatus('pedestal', 1)
		self.treeRawOutput.SetBranchStatus('pedestalSigma', 1)
		self.treeRawOutput.SetBranchStatus('pedestalAndSignal', 1)
		self.treeRawOutput.SetBranchStatus('pedestalAndSignalReal', 1)
		self.treeRawOutput.SetBranchStatus('signal', 1)
		self.treeRawOutput.SetBranchStatus('signalReal', 1)
		leng = self.treeRawOutput.Draw('pedestal:pedestalSigma:pedestalAndSignal', '', 'goff')
		if leng > self.treeRawOutput.GetEstimate():
			self.treeRawOutput.SetEstimate(leng)
			leng = self.treeRawOutput.Draw('pedestal:pedestalSigma:pedestalAndSignal', '', 'goff')
		tempPed = self.treeRawOutput.GetVal(0)
		tempPedSigma = self.treeRawOutput.GetVal(1)
		tempPedSignal = self.treeRawOutput.GetVal(2)
		tempPed = np.array([tempPed[ev] for ev in xrange(leng)], 'f8')
		tempPedSigma = np.array([tempPedSigma[ev] for ev in xrange(leng)], 'f8')
		tempPedSignal = np.array([tempPedSignal[ev] for ev in xrange(leng)], 'f8')
		leng = self.treeRawOutput.Draw('pedestalAndSignalReal:signal:signalReal', '', 'goff')
		if leng > self.treeRawOutput.GetEstimate():
			self.treeRawOutput.SetEstimate(leng)
			leng = self.treeRawOutput.Draw('pedestalAndSignalReal:signal:signalReal', '', 'goff')
		tempPedSignalReal = self.treeRawOutput.GetVal(0)
		tempSignal = self.treeRawOutput.GetVal(1)
		tempSignalReal = self.treeRawOutput.GetVal(2)
		tempPedSignalReal = np.array([tempPedSignalReal[ev] for ev in xrange(leng)], 'f8')
		tempSignal = np.array([tempSignal[ev] for ev in xrange(leng)], 'f8')
		tempSignalReal = np.array([tempSignalReal[ev] for ev in xrange(leng)], 'f8')

		tempPed = np.delete(tempPed, np.where(np.bitwise_or((tempPed == -10), (tempPed == -10.0))))
		tempPedSigma = np.delete(tempPedSigma, np.where(np.bitwise_or((tempPedSigma == -10), (tempPedSigma == -10.0))))
		tempPedSignal = np.delete(tempPedSignal, np.where(np.bitwise_or((tempPedSignal == -10), (tempPedSignal == -10.0))))
		tempPedSignalReal = np.delete(tempPedSignalReal, np.where(np.bitwise_or((tempPedSignalReal == -10), (tempPedSignalReal == -10.0))))
		tempSignal = np.delete(tempSignal, np.where(np.bitwise_or((tempSignal == -10), (tempSignal == -10.0))))
		tempSignalReal = np.delete(tempSignalReal, np.where(np.bitwise_or((tempSignalReal == -10), (tempSignalReal == -10.0))))

		if tempPed.size > 1 and tempPedSignal.size > 0 and tempSignal.size > 0 and tempPedSignalReal.size > 0 and tempSignalReal.size > 0:
			self.ped[self.vcal].fill(tempPed.mean())
			self.pedSigma[self.vcal].fill(tempPedSigma.std())
			self.pedSignal[self.vcal].fill(tempPedSignal.mean())
			self.pedSignalReal[self.vcal].fill(tempPedSignalReal.mean())
			self.signal[self.vcal].fill(tempSignal.mean())
			self.signalReal[self.vcal].fill(tempSignalReal.mean())
		else:
			self.ped[self.vcal].fill(0)
			self.pedSigma[self.vcal].fill(0)
			self.pedSignal[self.vcal].fill(0)
			self.pedSignalReal[self.vcal].fill(0)
			self.signal[self.vcal].fill(0)
			self.signalReal[self.vcal].fill(0)

	def LoadInputScalars(self):
		self.LoadInputTreeRaw(False)
		self.treeRawInput.SetBranchStatus('*', 0)
		self.treeRawInput.SetBranchStatus('pedestal1', 1)
		self.treeRawInput.SetBranchStatus('pedestalAndCalV1', 1)
		self.treeRawInput.SetBranchStatus('calVolt1', 1)
		self.treeRawInput.SetBranchStatus('pedestal2', 1)
		self.treeRawInput.SetBranchStatus('pedestalAndCalV2', 1)
		self.treeRawInput.SetBranchStatus('calVolt2', 1)
		self.treeRawInput.SetBranchStatus('pedestal3', 1)
		self.treeRawInput.SetBranchStatus('pedestalAndCalV3', 1)
		self.treeRawInput.SetBranchStatus('calVolt3', 1)
		self.treeRawInput.SetBranchStatus('pedestal4', 1)
		self.treeRawInput.SetBranchStatus('pedestalAndCalV4', 1)
		self.treeRawInput.SetBranchStatus('calVolt4', 1)
		
		leng = self.treeRawInput.Draw('pedestal1:pedestalAndCalV1:calVolt1', '', 'goff')
		if leng > self.treeRawInput.GetEstimate():
			self.treeRawInput.SetEstimate(leng)
			leng = self.treeRawInput.Draw('pedestal1:pedestalAndCalV1:calVolt1', '', 'goff')
		tempPed1 = self.treeRawInput.GetVal(0)
		tempPedCal1 = self.treeRawInput.GetVal(1)
		tempCal1 = self.treeRawInput.GetVal(2)
		tempPed1 = np.array([tempPed1[ev] for ev in xrange(leng)], 'f8')
		tempPedCal1 = np.array([tempPedCal1[ev] for ev in xrange(leng)], 'f8')
		tempCal1 = np.array([tempCal1[ev] for ev in xrange(leng)], 'f8')
		leng = self.treeRawInput.Draw('pedestal2:pedestalAndCalV2:calVolt2', '', 'goff')
		if leng > self.treeRawInput.GetEstimate():
			self.treeRawInput.SetEstimate(leng)
			leng = self.treeRawInput.Draw('pedestal2:pedestalAndCalV2:calVolt2', '', 'goff')
		tempPed2 = self.treeRawInput.GetVal(0)
		tempPedCal2 = self.treeRawInput.GetVal(1)
		tempCal2 = self.treeRawInput.GetVal(2)
		tempPed2 = np.array([tempPed2[ev] for ev in xrange(leng)], 'f8')
		tempPedCal2 = np.array([tempPedCal2[ev] for ev in xrange(leng)], 'f8')
		tempCal2 = np.array([tempCal2[ev] for ev in xrange(leng)], 'f8')
		leng = self.treeRawInput.Draw('pedestal3:pedestalAndCalV3:calVolt3', '', 'goff')
		if leng > self.treeRawInput.GetEstimate():
			self.treeRawInput.SetEstimate(leng)
			leng = self.treeRawInput.Draw('pedestal3:pedestalAndCalV3:calVolt3', '', 'goff')
		tempPed3 = self.treeRawInput.GetVal(0)
		tempPedCal3 = self.treeRawInput.GetVal(1)
		tempCal3 = self.treeRawInput.GetVal(2)
		tempPed3 = np.array([tempPed3[ev] for ev in xrange(leng)], 'f8')
		tempPedCal3 = np.array([tempPedCal3[ev] for ev in xrange(leng)], 'f8')
		tempCal3 = np.array([tempCal3[ev] for ev in xrange(leng)], 'f8')
		leng = self.treeRawInput.Draw('pedestal4:pedestalAndCalV4:calVolt4', '', 'goff')
		if leng > self.treeRawInput.GetEstimate():
			self.treeRawInput.SetEstimate(leng)
			leng = self.treeRawInput.Draw('pedestal4:pedestalAndCalV4:calVolt4', '', 'goff')
		tempPed4 = self.treeRawInput.GetVal(0)
		tempPedCal4 = self.treeRawInput.GetVal(1)
		tempCal4 = self.treeRawInput.GetVal(2)
		tempPed4 = np.array([tempPed4[ev] for ev in xrange(leng)], 'f8')
		tempPedCal4 = np.array([tempPedCal4[ev] for ev in xrange(leng)], 'f8')
		tempCal4 = np.array([tempCal4[ev] for ev in xrange(leng)], 'f8')

		tempPed1 = np.delete(tempPed1, np.where(np.bitwise_or((tempPed1 == -10), (tempPed1 == -10.0))))
		tempPedCal1 = np.delete(tempPedCal1, np.where(np.bitwise_or((tempPedCal1 == -10), (tempPedCal1 == -10.0))))
		tempCal1 = np.delete(tempCal1, np.where(np.bitwise_or((tempCal1 == -10), (tempCal1 == -10.0))))
		tempPed2 = np.delete(tempPed2, np.where(np.bitwise_or((tempPed2 == -10), (tempPed2 == -10.0))))
		tempPedCal2 = np.delete(tempPedCal2, np.where(np.bitwise_or((tempPedCal2 == -10), (tempPedCal2 == -10.0))))
		tempCal2 = np.delete(tempCal2, np.where(np.bitwise_or((tempCal2 == -10), (tempCal2 == -10.0))))
		tempPed3 = np.delete(tempPed3, np.where(np.bitwise_or((tempPed3 == -10), (tempPed3 == -10.0))))
		tempPedCal3 = np.delete(tempPedCal3, np.where(np.bitwise_or((tempPedCal3 == -10), (tempPedCal3 == -10.0))))
		tempCal3 = np.delete(tempCal3, np.where(np.bitwise_or((tempCal3 == -10), (tempCal3 == -10.0))))
		tempPed4 = np.delete(tempPed4, np.where(np.bitwise_or((tempPed4 == -10), (tempPed4 == -10.0))))
		tempPedCal4 = np.delete(tempPedCal4, np.where(np.bitwise_or((tempPedCal4 == -10), (tempPedCal4 == -10.0))))
		tempCal4 = np.delete(tempCal4, np.where(np.bitwise_or((tempCal4 == -10), (tempCal4 == -10.0))))

		self.pedInput1[self.vcal].fill(tempPed1.mean())
		self.pedInputSigma1[self.vcal].fill(tempPed1.std())
		self.pedCal1[self.vcal].fill(tempPedCal1.mean())
		self.calVoltsReal1[self.vcal].fill(tempCal1.mean())
		self.calVoltsRealSigma1[self.vcal].fill(tempCal1.std())
		self.pedInput2[self.vcal].fill(tempPed2.mean())
		self.pedInputSigma2[self.vcal].fill(tempPed2.std())
		self.pedCal2[self.vcal].fill(tempPedCal2.mean())
		self.calVoltsReal2[self.vcal].fill(tempCal2.mean())
		self.calVoltsRealSigma2[self.vcal].fill(tempCal2.std())
		self.pedInput3[self.vcal].fill(tempPed3.mean())
		self.pedInputSigma3[self.vcal].fill(tempPed3.std())
		self.pedCal3[self.vcal].fill(tempPedCal3.mean())
		self.calVoltsReal3[self.vcal].fill(tempCal3.mean())
		self.calVoltsRealSigma3[self.vcal].fill(tempCal3.std())
		self.pedInput4[self.vcal].fill(tempPed4.mean())
		self.pedInputSigma4[self.vcal].fill(tempPed4.std())
		self.pedCal4[self.vcal].fill(tempPedCal4.mean())
		self.calVoltsReal4[self.vcal].fill(tempCal4.mean())
		self.calVoltsRealSigma4[self.vcal].fill(tempCal4.std())

		# TODO: use options output real or not real, and input 1, 2, 3 or 4 for considering different cases due to our ignorance...
	def SaveCalibrations(self, suffix='', realCal=4):
		self.fileCal = ro.TFile('calibrations_{s}.root'.format(s=suffix), 'RECREATE')
		self.CreateVcalVSignalRegression(realCal)
		self.CreateChargeVSignalRegression(realCal)
		self.fileCal.Close()

	def CreateVcalVSignalRegression(self, realCal=4):
		if self.fileCal.GetOption() == 'READ':
			self.fileCal.ReOpen('UPDATE')
		syst = self.vcal_syst_factor
		if realCal == 1:
			xaxis = np.array([self.calVoltsReal1[vcal] for vcal in self.calVolts], 'f8')
			xaxisSigma = np.array([np.sqrt((self.calVoltsReal1[vcal] * syst)**2 + self.calVoltsRealSigma1[vcal]**2) for vcal in self.calVolts], 'f8')
		elif realCal == 2:
			xaxis = np.array([self.calVoltsReal2[vcal] for vcal in self.calVolts], 'f8')
			xaxisSigma = np.array([np.sqrt((self.calVoltsReal2[vcal] * syst)**2 + self.calVoltsRealSigma2[vcal]**2) for vcal in self.calVolts], 'f8')
		elif realCal == 3:
			xaxis = np.array([self.calVoltsReal3[vcal] for vcal in self.calVolts], 'f8')
			xaxisSigma = np.array([np.sqrt((self.calVoltsReal3[vcal] * syst)**2 + self.calVoltsRealSigma3[vcal]**2) for vcal in self.calVolts], 'f8')
		else:
			xaxis = np.array([self.calVoltsReal4[vcal] for vcal in self.calVolts], 'f8')
			xaxisSigma = np.array([np.sqrt((self.calVoltsReal4[vcal] * syst)**2 + self.calVoltsRealSigma4[vcal]**2) for vcal in self.calVolts], 'f8')
		yaxis = np.array([self.signal[vcal] for vcal in self.calVolts], 'f8')
		yaxisSigma = np.array([self.pedSigma[vcal] for vcal in self.calVolts], 'f8')
		self.vCalSignalGraph = ro.TGraphErrors(len(self.calVolts), xaxis, yaxis, xaxisSigma, yaxisSigma)
		namegraph = 'Signal_vs_Vcal'
		self.vCalSignalGraph.SetNameTitle(namegraph, namegraph)
		self.vCalSignalGraph.GetXaxis().SetTitle('Vcal/V')
		self.vCalSignalGraph.GetYaxis().SetTitle('Signal/V')
		self.vCalSignalFit = self.vCalSignalGraph.Fit('pol1', 'QMFS')
		self.fileCal.cd()
		self.vCalSignalGraph.Write()
		self.vCalSignalFit.Write()
		self.fileCal.Write()

	def CreateChargeVSignalRegression(self, realCal=4):
		if self.fileCal.GetOption() == 'READ':
			self.fileCal.ReOpen('UPDATE')
		self.CalculateChargeFromVcal(realCal)
		xaxis = np.array([self.signal[vcal] for vcal in self.calVolts], 'f8')
		xaxisSigma = np.array([self.pedSigma[vcal] for vcal in self.calVolts], 'f8')
		yaxis = np.array([self.calCharge[vcal] for vcal in self.calVolts], 'f8')
		yaxisSigma = np.array([self.calChargeSigma[vcal] for vcal in self.calVolts], 'f8')
		self.chargeSignalGraph = ro.TGraphErrors(len(self.calVolts), xaxis, yaxis, xaxisSigma, yaxisSigma)
		namegraph = 'Charge_vs_Signal'
		self.chargeSignalGraph.SetNameTitle(namegraph, namegraph)
		self.chargeSignalGraph.GetXaxis().SetTitle('Signal/V')
		self.chargeSignalGraph.GetYaxis().SetTitle('Charge/e')
		self.chargeSignalFit = self.chargeSignalGraph.Fit('pol1', 'QMFS')
		self.fileCal.cd()
		self.chargeSignalFit.Write()
		self.chargeSignalGraph.Write()
		self.fileCal.Write()

	#   TODO: Redo this method
	def CalculateAcceptedEventsOutput(self):
		self.acceptEventsOutput = 1 - np.zeros(self.eventsOutput, '?')
	# 	self.CalculateGoodPedestalCut()
	#
	# def CalculateGoodPedestalCut(self):
	# 	print 'Selecting good events based on pedestal'
	# 	timeCut1, timeCut2 = TCut('time>{vl}'.format(vl=(-self.pedestalTEndPos-self.averagingTime))), TCut('time<{vh}'.format(vh=(-self.pedestalTEndPos)))
	# 	timeCuts = timeCut1+timeCut2
	# 	self.treeRawOutput.SetBranchStatus('*', 0)
	# 	self.treeRawOutput.SetBranchStatus('voltageSignal', 1)
	# 	self.treeRawOutput.SetBranchStatus('time', 1)
	# 	num = self.treeRawOutput.Draw('voltageSignal>>htemp', timeCuts, 'goff')
	# 	htemp = gDirectory.Get('htemp')
	# 	mean = htemp.GetMean()
	# 	sigm = htemp.GetStdDev()
	# 	meansq = mean ** 2 + sigm * (num - 1)/num
	# 	for i in xrange(5):
	# 		for ev in self.eventVectOutput:
	# 			if self.acceptEventsOutput[ev]:
	# 				numev = self.treeRawOutput.Draw('voltageSignal>>htemp', timeCuts, 'goff', 1, ev)
	# 				htemp = gDirectory.Get('htemp')
	# 				evMean = htemp.GetMean()
	# 				evSigm = htemp.GetStdDev()
	# 				evMeansq = evMean ** 2 + evSigm * (numev - 1)/numev
	# 				if abs(evMean - mean) > 4 * sigm:
	# 					self.eventVectOutput[ev] = False
	# 					mean = (mean * num - evMean * numev) / (num - numev)
	# 					sigm = (num * meansq - numev * evMeansq - ((num * mean - numev * evMean) ** 2)/(num-numev))/(num-numev-1)
	# 					num = num - numev
	# 					meansq = mean ** 2 + sigm * (num - 1)/num

	def CalculateAcceptedEventsInput(self):
		self.acceptEventsInput = 1 - np.zeros(self.eventsInput, '?')
	# 	self.CalculateGoodPedestalCut()
	#
	# def CalculateGoodPedestalCut(self):
	# 	print 'Selecting good events based on pedestal'
	# 	timeCut1, timeCut2 = TCut('time>{vl}'.format(vl=(-self.pedestalTEndPos-self.averagingTime))), TCut('time<{vh}'.format(vh=(-self.pedestalTEndPos)))
	# 	timeCuts = timeCut1+timeCut2
	# 	self.treeRawInput.SetBranchStatus('*', 0)
	# 	self.treeRawInput.SetBranchStatus('voltageSignal', 1)
	# 	self.treeRawInput.SetBranchStatus('time', 1)
	# 	num = self.treeRawInput.Draw('voltageSignal>>htemp', timeCuts, 'goff')
	# 	htemp = gDirectory.Get('htemp')
	# 	mean = htemp.GetMean()
	# 	sigm = htemp.GetStdDev()
	# 	meansq = mean ** 2 + sigm * (num - 1)/num
	# 	for i in xrange(5):
	# 		for ev in self.eventVectInput:
	# 			if self.acceptEventsInput[ev]:
	# 				numev = self.treeRawInput.Draw('voltageSignal>>htemp', timeCuts, 'goff', 1, ev)
	# 				htemp = gDirectory.Get('htemp')
	# 				evMean = htemp.GetMean()
	# 				evSigm = htemp.GetStdDev()
	# 				evMeansq = evMean ** 2 + evSigm * (numev - 1)/numev
	# 				if abs(evMean - mean) > 4 * sigm:
	# 					self.eventVectInput[ev] = False
	# 					mean = (mean * num - evMean * numev) / (num - numev)
	# 					sigm = (num * meansq - numev * evMeansq - ((num * mean - numev * evMean) ** 2)/(num-numev))/(num-numev-1)
	# 					num = num - numev
	# 					meansq = mean ** 2 + sigm * (num - 1)/num

	def Quit(self):
		# sys.exit(0)
		exit()

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
	parser = OptionParser()
	parser.add_option('-o', '--outdir', dest='outdir', default='./', type='string',
					  help='Relative path to output directory e.g. ./')
	parser.add_option('-i', '--infile', dest='infile', type='string',
					  help='path to file to convert and analyse')
	parser.add_option('-a', '--automatic', dest='automatic', default=False, action='store_true',
					  help='Toggles automatic analysis')
	parser.add_option('-b', '--bias', dest='bias', default=-400, type='float',
					  help='The HV bias used on the sample e.g. -400')
	parser.add_option('-v', '--verbose', dest='verb', default=False, help='Toggles verbose', action='store_true')
	(options, args) = parser.parse_args()
	outdir = str(options.outdir)
	infile = str(options.infile)
	automatic = bool(options.automatic)
	bias = float(options.bias)
	verb = bool(options.verb)
	z = CalibrationAnalysis(outdir, infile, bias, verb)
	if automatic:
		z.AnalyseAllCalibrationFiles()
