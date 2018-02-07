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


def LookForTime0(trigVolts, points, percent_post, time_res, trigVal):
	guess_pos = int(round(points * (100.0 - percent_post)/100.0))
	array_points = np.arange(points, dtype=np.dtype('uint32'))
	condition_trigg = np.array(np.abs(array_points - guess_pos) <= int(round(0.5e-6/time_res)), dtype='?')
	condition_no_trigg = np.array(1 - condition_trigg, dtype='?')
	mean = np.extract(condition_no_trigg, trigVolts).mean()
	sigma = np.extract(condition_no_trigg, trigVolts).std()
	temp_trig_volts = deepcopy(trigVolts)
	np.putmask(temp_trig_volts, condition_no_trigg, 100)
	volt_min_pos = temp_trig_volts.argmin()
	condition_trigg = np.bitwise_and(condition_trigg, np.array(array_points <= volt_min_pos))
	np.putmask(temp_trig_volts, np.bitwise_not(condition_trigg), 100)
	position0 = np.abs(temp_trig_volts - trigVal).argmin()
	return position0

def ADC_to_Volts(adcs, sigres, nbits, offset_p):
	return np.multiply(sigres, np.add(adcs, np.multiply(2 ** nbits - 1, offset_p / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')

if __name__ == '__main__':
	run_dir_location = str(sys.argv[1])
	working_dir_location = str(sys.argv[2])
	filename = str(sys.argv[3])
	signal_ch = int(sys.argv[4])
	trigger_ch = int(sys.argv[5])
	anti_co_ch = int(sys.argv[6])
	points = int(sys.argv[7])
	num_events = int(sys.argv[8])
	struct_len = int(sys.argv[9])
	struct_fmt = str(sys.argv[10])
	adc_res = np.double(sys.argv[11])
	sig_offset = float(sys.argv[12])
	trig_offset = float(sys.argv[13])
	anti_co_offset = float(sys.argv[14])
	time_res = np.double(sys.argv[15])
	post_trig_percent = float(sys.argv[16])
	trig_value = np.double(sys.argv[17])
	dig_bits = int(sys.argv[18])
	time_recalib = float(sys.argv[19])
	simultaneous_conversion = bool(sys.argv[20] != '0')

	if simultaneous_conversion:
		print 'Start creating root file simultaneously with data taking'
	else:
		print str(sys.argv)
		print 'Start creating root file'
	t0 = time.time()
	raw_file = ro.TFile('{wd}/{d}/Runs/{r}.root'.format(wd=working_dir_location, d=run_dir_location, r=filename), 'RECREATE')
	raw_tree = ro.TTree(filename, filename)
	if anti_co_ch != -1:
		vetoBra = np.zeros(points, 'f8')
	eventBra = np.zeros(1, 'I')
	voltBra = np.zeros(points, 'f8')
	trigBra = np.zeros(points, 'f8')
	timeBra = np.zeros(points, 'f8')
	raw_tree.Branch('event', eventBra, 'event/i')
	raw_tree.Branch('time', timeBra, 'time[{s}]/D'.format(s=points))
	raw_tree.Branch('voltageSignal', voltBra, 'voltageSignal[{s}]/D'.format(s=points))
	raw_tree.Branch('voltageTrigger', trigBra, 'voltageTrigger[{s}]/D'.format(s=points))
	if anti_co_ch != -1:
		raw_tree.Branch('voltageVeto', vetoBra, 'voltageVeto[{s}]/D'.format(s=points))
	if not simultaneous_conversion:
		bar = None
		bar = CreateProgressBar(num_events, bar)
		bar.start()

	signal_written_events = int(round(os.path.getsize('{d}/raw_wave{s}.dat'.format(d=working_dir_location, s=signal_ch)) / struct_len))
	trigger_written_events = int(round(os.path.getsize('{d}/raw_wave{t}.dat'.format(d=working_dir_location, t=trigger_ch)) / struct_len))
	fs = open('{wd}/raw_wave{s}.dat'.format(wd=working_dir_location, s=signal_ch), 'rb')
	ft = open('{wd}/raw_wave{t}.dat'.format(wd=working_dir_location, t=trigger_ch), 'rb')
	if anti_co_ch != -1:
		anti_co_written_events = int(round(os.path.getsize('{d}/raw_wave{a}.dat'.format(d=working_dir_location, a=anti_co_ch)) / struct_len))
		fa = open('{wd}/raw_wave{a}.dat'.format(wd=working_dir_location, a=anti_co_ch), 'rb')

	for ev in xrange(num_events):
		t1 = time.time()
		wait_for_data = True if (signal_written_events <= ev or trigger_written_events <= ev) else False
		if anti_co_ch != -1:
			wait_for_data = bool(wait_for_data or (anti_co_written_events <= ev))
		while wait_for_data:
			if simultaneous_conversion:
				if time.time() - t1 - 10 > time_recalib:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=(time_recalib + 10))
					exit()
				if not fs.closed:
					fs.close()
				if not ft.closed:
					ft.close()
				if anti_co_ch != -1:
					if not fa.closed:
						fa.close()
				signal_written_events = int(round(os.path.getsize('{d}/raw_wave{s}.dat'.format(d=working_dir_location, s=signal_ch)) / struct_len))
				trigger_written_events = int(round(os.path.getsize('{d}/raw_wave{t}.dat'.format(d=working_dir_location, t=trigger_ch)) / struct_len))
				wait_for_data = True if (signal_written_events <= ev or trigger_written_events <= ev) else False
				if anti_co_ch != -1:
					anti_co_written_events = int(round(os.path.getsize('{d}/raw_wave{a}.dat'.format(d=working_dir_location, a=anti_co_ch)) / struct_len))
					wait_for_data = bool(wait_for_data or (anti_co_written_events <= ev))
				if not wait_for_data:
					fs = open('{wd}/raw_wave{s}.dat'.format(wd=working_dir_location, s=signal_ch), 'rb')
					ft = open('{wd}/raw_wave{t}.dat'.format(wd=working_dir_location, t=trigger_ch), 'rb')
					if anti_co_ch != -1:
						fa = open('{wd}/raw_wave{a}.dat'.format(wd=working_dir_location, a=anti_co_ch), 'rb')
			else:
				print 'The data is corrupted... exiting'
				exit()

		fs.seek(ev * struct_len, 0)
		datas = fs.read(struct_len)
		ft.seek(ev * struct_len, 0)
		datat = ft.read(struct_len)
		if anti_co_ch != -1:
			fa.seek(ev * struct_len, 0)
			dataa = fa.read(struct_len)
		if not datas or not datat:
			print 'No event in signal or trigger files... exiting'
			exit()
		elif anti_co_ch != -1:
			if not dataa:
				print 'No event in veto trigger file... exiting'

		s = struct.Struct(struct_fmt).unpack_from(datas)
		signalADCs = np.array(s, 'H')
		signalVolts = ADC_to_Volts(signalADCs, adc_res, dig_bits, sig_offset)
		t = struct.Struct(struct_fmt).unpack_from(datat)
		triggADCs = np.array(t, 'H')
		triggVolts = ADC_to_Volts(triggADCs, adc_res, dig_bits, trig_offset)
		trigPos = LookForTime0(trigVolts=triggVolts, points=points, percent_post=post_trig_percent, time_res=time_res, trigVal=trig_value)
		timeVect = np.linspace(-trigPos * time_res, time_res * (points - 1 - trigPos), points, dtype='f8')
		if anti_co_ch != -1:
			ac = struct.Struct(struct_fmt).unpack_from(dataa)
			vetoADCs = np.array(ac, 'H')
			vetoVolts = ADC_to_Volts(vetoADCs, adc_res, dig_bits, anti_co_offset)

		eventBra.fill(ev)
		np.putmask(timeBra, np.bitwise_not(np.zeros(points, '?')), timeVect)
		np.putmask(voltBra, np.bitwise_not(np.zeros(points, '?')), signalVolts)
		np.putmask(trigBra, np.bitwise_not(np.zeros(points, '?')), triggVolts)
		if anti_co_ch != -1:
			np.putmask(vetoBra, np.bitwise_not(np.zeros(points, '?')), vetoVolts)
		numFil = raw_tree.Fill()
		if not simultaneous_conversion:
			bar.update(ev + 1)
	if not simultaneous_conversion:
		bar.finish()
	raw_file.Write()
	raw_file.Close()
	fs.close()
	ft.close()
	if anti_co_ch != -1:
		fa.close()
	t0 = time.time() - t0
	print 'Time creating root tree:', t0, 'seconds'
	exit()
