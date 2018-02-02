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

def CreateProgressBar(maxVal=1, bar=None):
	widgets = [
		'Processed: ', progressbar.Counter(),
		' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
		' ', progressbar.Bar(marker='>'),
		' ', progressbar.Timer(),
		' ', progressbar.ETA()
		# ' ', progressbar.AdaptativeETA(),
		#  ' ', progressbar.AdaptativeTransferSpeed()
	]
	bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)
	return bar


def LookFor0Time(trigVolts):


if __name__ == '__main__':

	run_dir_location = str(sys.argv[1])
	filename = str(sys.argv[2])
	signal_ch = int(sys.argv[3])
	trigger_ch = int(sys.argv[4])
	anti_co_ch = int(sys.argv[5])
	points = int(sys.argv[6])
	num_events = int(sys.argv[7])
	struct_len = int(sys.argv[8])
	struct_fmt = str(sys.argv[9])
	adc_res = np.double(sys.argv[10])
	sig_offset = float(sys.argv[11])
	trig_offset = float(sys.argv[12])
	anti_co_offset = float(sys.argv[13])
	time_res = np.double(sys.argv[14])
	simultaneous_conversion = bool(sys.argv[15] != '0')

	if simultaneous_conversion:
		print 'Start creating root file simultaneously with data taking'
	else:
		print 'Start creating root file'
	t0 = time.time()
	raw_file = ro.TFile('{d}/Runs/{r}.root'.format(d=run_dir_location, r=filename), 'RECREATE')
	raw_tree = ro.TTree(filename, filename)
	fs = open('raw_wave{s}.dat'.format(s=signal_ch), 'rb')
	ft = open('raw_wave{t}.dat'.format(t=trigger_ch), 'rb')
	if anti_co_ch != -1:
		fa = open('raw_wave{a}.dat'.format(a=anti_co_ch), 'rb')
		vetoBra = np.zeros(points, 'f8')
	eventBra = np.zeros(1, 'I')
	voltBra = np.zeros(points, 'f8')
	trigBra = np.zeros(points, 'f8')
	timeBra = np.zeros(points, 'f8')
	raw_tree.Branch('event', eventBra, 'event/i')
	raw_tree.Branch('time', timeBra, 'time[{s}]/D'.format(s=points))
	raw_tree.Branch('voltageSignal', voltBra, 'voltageSignal[{s}]/D'.format(s=points))
	raw_tree.Branch('voltageTrigger', trigBra, 'voltageTrigger[{s}]/D'.format(s=points))
	raw_tree.Branch('voltageVeto', vetoBra, 'voltageVeto[{s}]/D'.format(s=points))

	if not simultaneous_conversion:
		bar = None
		bar = CreateProgressBar(num_events, bar)
		bar.start()

	for ev in xrange(num_events):
		t1 = time.time()

		wait_for_data = True
		while wait_for_data:
			if time.time() - t1 > 60:
				print 'No data has been saved in file for event {ev}... exiting!'.format(ev=ev)
				exit()
			fs.seek(0, 2)
			signal_written_events = int(round(fs.tell() / struct_len))
			trigger_written_events = int(round(ft.tell() / struct_len))
			if anti_co_ch != -1:
				anti_co_written_events = int(round(fa.tell() / struct_len))
			if signal_written_events + trigger_written_events <= 2 * ev:
				if anti_co_ch != -1:
					if anti_co_written_events <= ev:
						pass
				else:
					pass
			else:
				if anti_co_ch != -1:
					if anti_co_written_events > ev:
						wait_for_data = False

		fs.seek(ev * struct_len)
		ft.seek(ev * struct_len)
		datas = fs.read(struct_len)
		datat = ft.read(struct_len)
		if anti_co_ch != -1:
			fa.seek(ev * struct_len)
			dataa = fa.read(struct_len)
		if not datas or not datat:
			print 'No event in signal or trigger files... exiting'
			exit()
		elif anti_co_ch != -1:
			if not dataa:
				print 'No event in veto trigger file... exiting'

		s = struct.Struct(struct_fmt).unpack_from(datas)
		signalADCs = np.array(s, 'H')
		t = struct.Struct(struct_fmt).unpack_from(datat)
		triggADCs = np.array(t, 'H')
		signalVolts = np.array(np.multiply(signalADCs, adc_res) + sig_offset / 50.0 - 1, 'f8')
		triggVolts = np.array(np.multiply(triggADCs, adc_res) + trig_offset / 50.0 - 1, 'f8')
		left, right = np.double(triggVolts[0]), np.double(triggVolts[-1])
		mid = np.double((left + right) / 2.0)
		distFromMid = np.array(np.abs(triggVolts - mid), 'f8')
		midPos = distFromMid.argmin()
		timeVect = np.linspace(-midPos * self.time_res, self.time_res * (points - 1 - midPos), points,
		                       dtype='f8')
		eventBra.fill(ev)
		np.putmask(timeBra, 1 - np.zeros(points, '?'), timeVect)
		np.putmask(voltBra, 1 - np.zeros(points, '?'), signalVolts)
		numFil = raw_tree.Fill()
		self.bar.update(ev + 1)
	self.bar.finish()
	raw_file.Write()
	raw_file.Close()
	fs.close()
	ft.close()
	t0 = time.time() - t0
	print 'Time creating root tree:', t0, 'seconds'
