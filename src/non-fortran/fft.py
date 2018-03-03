#! /usr/bin/env python
from __future__ import print_function
import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt
from scipy.signal import blackman
import sys
import os.path
import math

def joinwords(*words):
    joined = ''.join(words)
    return joined

def message(*objs):
    print(joinwords(*objs), file=sys.stderr)

def out_message(*objs):
    print(joinwords(*objs), file=sys.stdout)

def error(*objs):
    message("ERROR: ", *objs)
    sys.exit(0)

def warning(*objs):
    message("WARNING: ", *objs)

def read_warning():
    warning('In input file "',masterin,'" at line ',str(nl),'. Parameter name "', p, '" not recognised.')

def find_nearest(array_in, value):
    idx = (np.abs(np.asarray(array_in)-value)).argmin()
    return idx, array_in[idx]

narg = len(sys.argv)-1

if narg > 0:
    masterin = str(sys.argv[1])
else:
    masterin = 'fft.in'

if not os.path.isfile(masterin):
    fid = open(masterin, 'w')
    print('# Input file for fft.py',file=fid)
    print('',file=fid)
    print('perform       fft # max / integrate / average / fft',file=fid)
    print('infile        test.out',file=fid)
    print('field_file     test.field.out',file=fid)
    print('outfile       test.out.ft',file=fid)
    print('readcol       6',file=fid)
    print('field_col     2',file=fid)
    print('delay         50',file=fid)
    print('from          1900',file=fid)
    print('to            2000',file=fid)
    print('#units_time    fs',file=fid)
    print('#damp_factor   0.1 # Lorentzian damping in eV',file=fid)
    print('#yambo_delta   # Correct field to match YPP output',file=fid)
    print('#yambo_timestep 0.1e-2  # Timestep used in yambo_rt (in fs)',file=fid)
    print('#EELS          # Output 1/FFT also', file=fid)
    fid.close()
    os.system('vi '+masterin)
    exit()
#    error('Input file "',masterin,'" does not exist. Temporary fft.in written.')

#Defaults
readcol = 2
outfile = ''
field_file = ''
start_from = 0.0
go_to = np.inf
delay = 0.0
div_by_field = 0
perform = 'fft'
time_par_type = 'au'
fs_to_au = 1.0
damp_factor = 0.0
damp_function = 1.0
start_id = 0  ################### to shift delta over for yambo
speed_of_light = 137.03599911
yambo_delta = 0
yambo_timestep = 1
calc_eels = 0

nl = 0
with open(masterin) as foo:
    for line in foo:
        nl = nl + 1
        columns = line.split()
        if len(columns) > 0:
            p = columns[0]
            try: a = columns[1]
            except: a = 0
            if p[0] == '#':
                continue
            if p == str('infile'):
                infile = a
                if not os.path.isfile(infile):
                    error('In input file "',masterin,'" at line ',str(nl),'. File "', infile, '" does not exist.')

            elif p == 'field_file':
                field_file = a
                if not os.path.isfile(field_file):
                    error('In input file "',masterin,'" at line ',str(nl),'.  File "', field_file, '" does not exist.')
            elif p == str('outfile'):
                outfile = a
            elif p == 'readcol':
                readcol = int(a)
            elif p == 'from':
                start_from = float(a)
            elif p == 'to':
                go_to = float(a)
            elif p == 'delay':
                delay = float(a)
            elif p == 'field_col':
                field_col = int(a)
                div_by_field = 1
            elif p == 'perform':
                perform = a
            elif p == 'units_time':
                time_par_type = a
            elif p == 'damp_factor':
                damp_factor = float(a)
            elif p == 'yambo_delta':
                yambo_delta = 1
                start_id = 1
            elif p == 'yambo_timestep':
                yambo_timestep = float(a)
            elif p == 'EELS':
                calc_eels = 1
            else:
                read_warning()
foo.close()
#infile = 'ft_00000100.out'
#outfile = 'test.dat'
#readcol = 2

freq_au_to_ev = 27.2113834

x = []
y = []
if div_by_field == 1: f = [] 

if time_par_type == 'fs':
    fs_to_au = 41.341373336561361

delay      = delay * fs_to_au
start_from = start_from * fs_to_au
go_to      = go_to * fs_to_au
damp_factor = damp_factor / freq_au_to_ev

nl = 0
with open(infile) as foo:
    for line in foo:
        li = line.strip()
        if li.startswith("#"):
            continue
        nl = nl+1
        columns = line.split()
        if nl == 1:
            try:
                xtmp = float(columns[0])*fs_to_au
                name = joinwords('col_',str(readcol))
            except:
                name = str(columns[readcol - 1])
                continue
        if len(columns) > 0:
            try:
                xtmp = float(columns[0])*fs_to_au
            except:
                if nl == 2:
                    error('Error in reading input file "',infile,'".')
                else:
                    message('Input file "',infile,'" read successfully.')
                    break
            if np.abs(xtmp) >= start_from and np.abs(xtmp) <= go_to :
                x.append(xtmp)
                y.append(float(columns[readcol - 1]))
                if (div_by_field == 1 and field_file==''): 
                    f.append(float(columns[field_col - 1]))
foo.close()

if not field_file == '':
    with open(field_file) as foo:
        for line in foo:
            li = line.strip()
            if li.startswith("#"):
                continue
            nl = nl+1
            columns = line.split()
            if nl == 1:
                try:
                    xtmp = float(columns[0])*fs_to_au
                    name = joinwords('col_',str(field_col))
                except:
                    name = str(columns[field_col - 1])
                    continue
            if len(columns) > 0:
                try:
                    xtmp = float(columns[0])*fs_to_au
                except:
                    if nl == 2:
                        error('Error in reading input file "',infile,'".')
                    else:
                        message('Input file "',infile,'" read successfully.')
                        break
                if np.abs(xtmp) >= start_from and np.abs(xtmp) <= go_to :
                    if (div_by_field == 1): 
                        f.append(float(columns[field_col - 1]))
foo.close()

T = x[1] - x[0]

delayid, delayapp = find_nearest(x, delay)

N = len(y)

if perform == 'integrate':
    area = T*sum(y[0:len(y)-1])
    message('Area between ',str(x[0]),' and ',str(x[len(x)-1]), ' is:')
    out_message(str(area))
elif perform == 'average':
    av = sum(y)/float(N)
    message('Average value between ',str(x[0]),' and ',str(x[len(x)-1]), ' is:')
    out_message(str(av))
elif perform == 'max':
    maxm = max(y)
    message('Maximum value between ',str(x[0]),' and ',str(x[len(x)-1]), 'is:')
    out_message(str(maxm))
elif perform == 'fft':

    if outfile == '':
        outfile = joinwords(infile,'.',name,'.fft')
    fid = open(outfile, 'w')


    if not damp_factor == 0.0:
        for i in range(N):
            y[i] = y[i] * np.exp(-damp_factor * x[i])

    if not delay == 0.0:
    ##    ytmp = y[N-delayid:N]
    ##    ytmp = ytmp + y[0:N-delayid]
        ytmp = y[delayid:N]
        ytmp = ytmp + y[0:delayid]

        if div_by_field == 1:
            ftmp = f[delayid:N]
            ftmp = ftmp + f[0:delayid]
    #    ytmp = y[delayid:N]
    #    ytmp = ytmp + delayid*[0]
    #    plt.plot(x,ytmp, 'r')
    #    plt.plot(x,y, 'k')
    #    plt.show()
    else:
        ytmp = y
        if (div_by_field == 1): ftmp = f

    #yf = fft(y)
    if yambo_delta and div_by_field: ftmp = [ftmp[i]*fs_to_au*yambo_timestep/T/5.14220652e11 for i in range(len(ftmp))]

    yf = fft(ytmp[start_id:len(ytmp)])
    if (div_by_field == 1): ff = fft(ftmp[start_id:len(ftmp)])

    N = len(yf)

#    for i in range(len(ytmp)): print(i,ytmp[i].imag, ytmp[i].real)
#    for i in range(len(yf)): print(i,yf[i].imag, yf[i].real)

    if div_by_field:
        eps = [yf[i]/ff[i] for i in range(len(yf))]
    else:
        eps = yf

    if yambo_delta and div_by_field: eps = [1. + eps[i]*4.*np.pi for i in range(len(eps))]

    if calc_eels: eels = [1./eps[i] for i in range(len(eps))]

    #w = blackman(N)
    #ywf = fft(y*w)

    xf = np.linspace(0.0, 1.0/(2.0*T ), math.floor(N/2))
#    for i in range(N/2):
#        yf[i] = yf[i]/np.exp(-1j*xf[i]*0.001*fs_to_au)/5.33802520488724E-11*4.*np.pi
    xf = 2.0*np.pi * freq_au_to_ev * xf

    xf = [2.*np.pi*freq_au_to_ev*float(n)/T/float(N) for n in range(N)]


    #### NB Note that we print out the conjugate ####
    if calc_eels:
        print("%-22s %-22s %-22s %-22s %-22s %-22s" % ('Freq.', 'FFT.real', 'FFT.imag', 'FFT.abs', 'EELS.real', 'EELS.imag'), file=fid) 
    else:
        print("%-22s %-22s %-22s %-22s" % ('Freq.', 'FFT.real', 'FFT.imag', 'FFT.abs'), file=fid)

    if (div_by_field == 1):
        for i in range(math.floor(N/2)):
            if calc_eels:
                print("%22.13G %22.13G %22.13G %22.13G %22.13G %22.13G" % (xf[i], eps[i].real, -eps[i].imag, np.abs(eps[i]), eels[i].real, eels[i].imag), file=fid)
            else:
                print("%22.13G %22.13G %22.13G %22.13G" % (xf[i], eps[i].real, -eps[i].imag, np.abs(eps[i])), file=fid)
    else:
        for i in range(math.floor(N/2)):
            print("%22.13G %22.13G %22.13G %22.13G" % (xf[i], 2.0/float(N)*yf[i].real, -2.0/float(N)*yf[i].imag, np.abs(2.0/float(N)*yf[i])), file=fid)

    if not delay == 0.0: message('Time Delay: ', str(delayapp), '.')
    message('FFT of column ',str(readcol), ' successfully output to "',outfile,'".')
    message('Resolution: ',str(xf[1]-xf[0]),' eV.')
    message('Max Energy: ',str(xf[math.floor(N/2)-1]),' eV.')
    fid.close()
else:
    error('Perform program: "',perform,'" not recognised.')
