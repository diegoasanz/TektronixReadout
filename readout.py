import visa
import numpy as np
from struct import unpack
import time #just for timing

if __name__ == '__main__':
	nmax = 50

	try:
		rm = visa.ResourceManager('@py')
		inst = rm.open_resource('TCPIP0::192.168.1.13::inst0::INSTR')
		print('Connected to: ', inst.ask('*idn?'))
	except:
		print('Check the IP address!')

	#tell scope what we want
	inst.write('DATA:SOU CH1')
	inst.write('DATA:WIDTH 10')
	inst.write('DATA:ENC RPB')

	#get the scalings
	ymult = float(inst.ask('WFMPRE:YMULT?'))
	yzero = float(inst.ask('WFMPRE:YZERO?'))
	yoff = float(inst.ask('WFMPRE:YOFF?'))
	xincr = float(inst.ask('WFMPRE:XINCR?'))

	#wf readout loop
	t0 = time.time()
	count = 0
	while(count<nmax):
		inst.write('CURVE?')
		data = inst.read_raw()
		headerlen = 2 + int(data[1])
		header = data[:headerlen]
		ADC_wave = data[headerlen:-1]

		ADC_wave = np.array(unpack('%sB' % len(ADC_wave),ADC_wave))
		Volts = (ADC_wave - yoff) * ymult  + yzero
		Time = np.arange(0, xincr * len(Volts), xincr)		
		count += 1
		print(count)
	t1 = time.time()

	print('Readout rate: ' + str(round(nmax/(t1-t0),3)) + ' Hz')
