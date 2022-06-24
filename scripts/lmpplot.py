#!/usr/bin/env python3
# coding=utf-8

import os, sys, re, math
import numpy as np

OMITSTEP = -1
ENDSTEP = -1
if len(sys.argv) >= 3 and sys.argv[2].isdigit():
    OMITSTEP = int(sys.argv[2])
if len(sys.argv) >= 4 and sys.argv[3].isdigit():
    ENDSTEP = int(sys.argv[3])

def read_log(log_file):
    inf=open(log_file).read()
    nMin = len(re.findall('\nMinimization', inf))
    nRun = len(re.findall('\nStep', inf)) - nMin
    print('%i Min, %i Run' %(nMin, nRun))

    run = -nMin
    types: [str]
    data: [[float]]
    START = False
    for line in open(log_file):
        line = line.strip()
        if line.startswith('Step'):
            run += 1
            if run > 0:
                START = True
            if run == 1:
                types = line.strip().split()
                data=[[] for i in types]
            continue
        if run > 0 and line.startswith('Loop'):
            START = False
            continue
        if run > 0 and START:
            try:
                step = int(line.strip().split()[0])
            except:
                continue
            if step < OMITSTEP:
                continue
            if ENDSTEP > 0 and step > ENDSTEP:
                break
            for i in range(0, len(types)):
                data[i].append(float(line.strip().split()[i]))
    return types, data

def average_of_blocks(l, nblock=5):
    ave_block = []
    bsize = int(math.ceil(len(l)/nblock))
    for i in range(nblock):
        block = l[i*bsize:(i+1)*bsize]
        ave_block.append(np.mean(block))
    return ave_block

def block_average(l, nblock=5):
    ave_block = average_of_blocks(l, nblock)
    return np.mean(ave_block), np.std(ave_block, ddof=1)/math.sqrt(nblock)

def plot_data(types, data):
    print('Select the property to plot, or input any letter to quit:')
    option = 'File: %s, Steps: %i-%i, Samples: %i\n' %(sys.argv[1], data[0][0], data[0][-1], len(data[0]))
    for i in range(1, len(types)):
        ave, stderr = block_average(data[i])
        inv_blocks = [1000 / ave for ave in average_of_blocks(data[i])]
        inv_ave = np.mean(inv_blocks)
        inv_stderr = np.std(inv_blocks, ddof=1) / math.sqrt(len(inv_blocks))

        option += '%6i: %14s %10.4g  %10.4g  1E3/ %10.4g  %10.4g\n' %(i, types[i], ave, stderr, inv_ave, inv_stderr)
        # option += '%6i: %14s %10.4g  %10.4g\n' %(i, types[i], ave, stderr)
    print(option, end='')
    while True:
        plottype = input()
        if not plottype.isdigit():
            sys.exit()
        plottype = int(plottype)

        if plottype < 1 or plottype >= len(data):
            print('not valid')
        else:
            if not 'plt' in dir():
                import matplotlib.pyplot as plt
            plt.plot(data[0], data[plottype])
            plt.xlabel('Step')
            plt.ylabel(types[plottype])
            plt.show()

if __name__ == '__main__':
    types, data = read_log(sys.argv[1])
    plot_data(types, data)

