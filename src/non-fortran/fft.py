#! /usr/bin/env python
from __future__ import print_function
import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt
from scipy.signal import blackman
import sys
import os.path

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
    error('Input file "',masterin,'" does not exist.')

#Defaults
readcol = 2
outfile = ''
start_from = 0.0
go_to = np.inf
delay = 0.0
div_by_field = 0
perform = 'fft'

nl = 0
with open(masterin) as foo:
    for line in foo:
        nl = nl + 1
        columns = line.split()
        if len(columns) > 0:
            p = columns[0]
            a = columns[1]
            if p[0] == '#':
                continue
            if p == str('infile'):
                infile = a
                if not os.path.isfile(infile):
                    error('In input file "',masterin,'" at line ',str(nl),'. File "', infile, '" does not exist.')

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
            else:
                read_warning()

#infile = 'ft_00000100.out'
#outfile = 'test.dat'
#readcol = 2

freq_au_to_ev = 27.211

x = []
y = []
if div_by_field == 1: f = [] 

nl = 0
with open(infile) as foo:
    for line in foo:
        nl = nl+1
        columns = line.split()
        if nl == 1:
            try:
                xtmp = float(columns[0])
                name = joinwords('col_',str(readcol))
            except:
                name = str(columns[readcol - 1])
                continue
        if len(columns) > 0:
            try:
                xtmp = float(columns[0])
            except:
                if nl == 2:
                    error('Error in reading input file "',infile,'".')
                else:
                    message('Input file "',infile,'" read successfully.')
                    break
            if np.abs(xtmp) >= start_from and np.abs(xtmp) <= go_to :
                x.append(xtmp)
                y.append(float(columns[readcol - 1]))
                if (div_by_field == 1): f.append(float(columns[field_col - 1]))

if outfile == '':
    outfile = joinwords(infile,'.',name,'.fft')
fid = open(outfile, 'w')

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
elif perform == 'fft':

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
    yf = fft(ytmp)
    if (div_by_field == 1): ff = fft(ftmp)

    #w = blackman(N)
    #ywf = fft(y*w)

    xf = np.linspace(0.0, 1.0/(2.0*T ), N/2)
    xf = 2.0*np.pi * freq_au_to_ev * xf

    #### NB Note that we print out the conjugate ####
    print("%-22s %-22s %-22s %-22s" % ('Freq.', 'FFT.real', 'FFT.imag', 'FFT.abs'), file=fid)
    if (div_by_field == 1):
        for i in range(N/2):
            print("%22.13G %22.13G %22.13G %22.13G" % (xf[i], (yf[i]/ff[i]).real, -(yf[i]/ff[i]).imag, np.abs(yf[i]/ff[i])), file=fid)
    else:
        for i in range(N/2):
            print("%22.13G %22.13G %22.13G %22.13G" % (xf[i], 2.0/N*yf[i].real, -2.0/N*yf[i].imag, np.abs(2.0/N*yf[i])), file=fid)

    if not delay == 0.0: message('Time Delay: ', str(delayapp), '.')
    message('FFT of column ',str(readcol), ' successfully output to "',outfile,'".')
    message('Resolution: ',str(xf[1]-xf[0]),' eV.')
    message('Max Energy: ',str(xf[N/2-1]),' eV.')
else:
    error('Perform program: "',perform,'" not recognised.')
