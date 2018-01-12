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
		self.edge = -1
		self.name = str(ch)

if __name__ == '__main__':
	print 'bla'