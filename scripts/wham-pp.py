#!/usr/bin/env python3

import os
import sys
import getopt
import numpy as np

N_BIN_HIST = 0
OMIT = 1
SKIP = 1
PLOT = 0
PLUMED = 0


def usage():
    print('wham_pp.py -h[HELP] -o[omit_logline] -b[hist_bins] -p[PLOT] -s[skip] -d[PLUMED]')


opts, args = getopt.getopt(sys.argv[1:], 'hpb:o:s:d')
for op, value in opts:
    if op == '-b':
        N_BIN_HIST = int(value)
    elif op == '-p':
        PLOT = 1
    elif op == '-d':
        PLUMED = 1
    elif op == '-o':
        OMIT = int(value)
    elif op == '-s':
        SKIP = int(value)
    elif op == '-h':
        usage()
        sys.exit()


def split_str_digit_alpha(string):
    digit = ''
    alpha = ''
    elements = []
    for i in string:
        if i.isdigit():
            digit += i
            if alpha != '':
                elements.append(alpha)
                alpha = ''
        else:
            alpha += i
            if digit != '':
                elements.append(int(digit))
                digit = ''
    if digit != '':
        elements.append(int(digit))
    elif alpha != '':
        elements.append(alpha)
    return elements


r0_list = []
k_list = []
for i in sorted(os.listdir(os.getcwd()), key=split_str_digit_alpha):
    if os.path.isdir(i) and '-' in i and i[0].isdigit() and i[-1].isdigit():
        r0_list.append(i.split('-')[0])
        k_list.append(i.split('-')[1])

outf_meta = open('wham-meta', 'w')
for r0, k in zip(r0_list, k_list):
    outfile = 'wham-' + r0 + '-' + k + '.txt'
    outf = open(outfile, 'w')
    countline = 0
    if not PLUMED:
        cycle = 0
        start = False
        inpf = open(r0 + '-' + k + '/log')
        for line in inpf:
            if line.startswith('Loop') and cycle == 2:
                break
            if line.startswith('Step'):
                cycle += 1
                if cycle == 2:
                    start = True
            if start:
                countline += 1
                if countline > OMIT + 1 and countline % SKIP == 0:
                    str = line.strip().split()
                    time = float(str[0])
                    r = float(str[6])
                    outf.write('%.1f %.6f\n' % (time, r))
    else:
        inpf = open(r0 + '-' + k + '/DIS.txt')
        for line in inpf:
            countline += 1
            if countline > OMIT + 1 and countline % SKIP == 0:
                str = line.strip().split()
                time = float(str[0])
                r = float(str[1])
                outf.write('%.1f %.6f\n' % (time, r))
    inpf.close()
    outf.close()
    outf_meta.write('%s %s %s\n' % (outfile, r0, k))
    outf_meta.close()

hist_min = min(map(float, r0_list))
hist_max = max(map(float, r0_list))

if PLOT:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig = plt.gcf()
    fig.set_size_inches(16, 8)
    plt.grid(True)
    for i in sorted(os.listdir(os.getcwd()), key=split_str_digit_alpha):
        if i.startswith('wham-') and i.endswith('.txt'):
            array = np.genfromtxt(i)
            data = array[:, 1]
            if N_BIN_HIST:
                hist, bins = np.histogram(data, bins=N_BIN_HIST)
            else:
                hist, bins = np.histogram(data)
            center = (bins[:-1] + bins[1:]) / 2.
            plt.plot(center, hist, linewidth=2)
    plt.savefig('wham-hist.png')

# os.system('wham %f %f %s %s 298.15 0 wham-meta wham-pmf' %(hist_min, hist_max, bins, tol))
