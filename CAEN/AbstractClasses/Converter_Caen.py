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
	# ipdb.set_trace(context=7)
	guess_pos = int(round(points * (100.0 - percent_post)/100.0))
	array_points = np.arange(points, dtype=np.dtype('int32'))
	condition_trigg = np.array(np.abs(array_points - guess_pos) <= int(round(0.1e-6/time_res)), dtype='?')
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

def IsEventVetoed(vetoADC, points, trigPos, time_res, vetoVal=47):
	array_points = np.arange(points, dtype=np.dtype('int32'))
	window_around_trigg = 50e-9
	condition_base_line = np.array(np.abs(array_points - trigPos) > int(round(window_around_trigg/float(time_res))), dtype='?')
	condition_search = np.array(1 - condition_base_line, dtype='?')
	# condition_no_search = np.array(1 - condition_search, dtype='?')
	mean = np.extract(condition_base_line, vetoADC).mean()
	sigma = np.extract(condition_base_line, vetoADC).std()
	vetoValNew = 4 * sigma if vetoVal < 0.9 * 4 * vetoVal else vetoVal
	veto_event = bool((np.extract(condition_search, vetoADC) - mean + vetoValNew).min() <= 0)

	return veto_event

def ADC_to_Volts(adcs, sigres, nbits, offset_p):
	return np.multiply(sigres, np.add(adcs, np.multiply(2 ** nbits - 1, offset_p / 100.0 - 0.5, dtype='f8'), dtype='f8'), dtype='f8')

if __name__ == '__main__':
	run_dir_location = str(sys.argv[1])  # output directory location
	working_dir_location = str(sys.argv[2])  # current directory location, where the raw_wave files are
	filename = str(sys.argv[3])  # file name for the files
	signal_ch = int(sys.argv[4])  # caen channel for the ccd signal
	trigger_ch = int(sys.argv[5])  # caen channel for the trigger signal
	anti_co_ch = int(sys.argv[6])  # caen channel for the veto scintillator signal
	points = int(sys.argv[7])  # caen number of points in each event
	num_events = int(sys.argv[8])  # number of events to convert
	struct_len = int(sys.argv[9])  # structure length (obtained from struct_fmt and struct.calcsize(self.struct_fmt)
	struct_fmt = str(sys.argv[10])  # structure format of the data per event i.e. '@{p}H'.format(p=self.points)
	adc_res = np.double(sys.argv[11])  # adc resolution of the digitizer
	sig_offset = float(sys.argv[12])  # percentage offset for ccd signal (value within [-50, 50]) look at wavedump config file
	trig_offset = float(sys.argv[13])  # percentage offset for trigger signal (value within [-50, 50]) look at wavedump config file
	anti_co_offset = float(sys.argv[14])  # percentage offset for veto signal (value within [-50, 50]) look at wavedump config file
	time_res = np.double(sys.argv[15])  # time resolution of digitizer. i.e. 2e-9
	post_trig_percent = float(sys.argv[16])  # percentage of the acquired data after trigger (look at wave dump config file)
	trig_value = np.double(sys.argv[17])  # counts below base line for trigger
	veto_value = int(sys.argv[18])  # counts below base line on the veto signal for vetoing
	dig_bits = int(sys.argv[19])  # number of bits of the ADC i.e. 14
	simultaneous_conversion = bool(sys.argv[20] != '0')  # whether or not to do the simultaneous conversion while taking data

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
		vetoedBra = np.zeros(1, '?')
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
		raw_tree.Branch('vetoedEvent', vetoedBra, 'vetoedEvent/O')
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
				if time.time() - t1 > 60:
					print 'No data has been saved in file for event {ev} in the past {t} seconds... exiting!'.format(ev=ev, t=60)
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
			vetoed_event = IsEventVetoed(vetoADCs, points, trigPos, time_res)

		eventBra.fill(ev)
		np.putmask(timeBra, np.bitwise_not(np.zeros(points, '?')), timeVect)
		np.putmask(voltBra, np.bitwise_not(np.zeros(points, '?')), signalVolts)
		np.putmask(trigBra, np.bitwise_not(np.zeros(points, '?')), triggVolts)
		if anti_co_ch != -1:
			np.putmask(vetoBra, np.bitwise_not(np.zeros(points, '?')), vetoVolts)
			vetoedBra.fill(vetoed_event)
		numFil = raw_tree.Fill()
		bar.update(ev + 1)
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
