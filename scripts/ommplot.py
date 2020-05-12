#!/usr/bin/env python3
# coding=utf-8

import argparse
import sys
import math
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='Data file')
parser.add_argument('-b', '--begin', default=-1, type=float, help='Begin from this step')
parser.add_argument('-e', '--end', default=-1, type=float, help='End at this step')
parser.add_argument('--noplot', default=False, action='store_true', help='Do not plot')
parser.add_argument('--reciprocal', default=False, action='store_true', help='Calculate reciprocal')
args = parser.parse_args()

OMITSTEP = args.begin
ENDSTEP = args.end


def read_log(log_file):
    types: [str] = []
    data: [[float]] = []
    START = False
    for line in open(log_file):
        if line.startswith('#"'):
            types = [s.strip('"') for s in line.strip('#').strip().split('\t')]
            data = [[] for _ in types]
            START = True
            continue
        if START:
            words = line.strip().split()
            try:
                step = float(words[0])
            except:
                continue
            if step < OMITSTEP:
                continue
            if ENDSTEP > 0 and step > ENDSTEP:
                break
            for i in range(len(types)):
                data[i].append(float(words[i]))
    return types, data


def average_of_blocks(l, nblock=5):
    ave_block = []
    var_block = []
    bsize = int(math.ceil(len(l) / nblock))
    for i in range(nblock):
        block = l[i * bsize:(i + 1) * bsize]
        ave_block.append(np.mean(block))
        var_block.append(np.var(block))
    return ave_block, var_block


def block_average(l, nblock=5):
    ave_block, var_block = average_of_blocks(l, nblock)
    return np.mean(ave_block), np.std(ave_block, ddof=1) / math.sqrt(nblock), \
           np.mean(var_block), np.std(var_block, ddof=1) / math.sqrt(nblock)


def plot_data(types, data):
    option = 'File: %s, Steps: %i-%i, Samples: %i\n' % (
        args.input, data[0][0], data[0][-1], len(data[0]))
    for i in range(1, len(types)):
        ave, stderr, var, var_stderr = block_average(data[i])
        if not args.reciprocal:
            option += '%6i: %14s %10.4g %10.4g %10.4g %10.4g\n' % (
            i, types[i], ave, stderr, var, var_stderr)
        else:
            inv_blocks = [1000 / ave for ave in average_of_blocks(data[i])]
            inv_ave = np.mean(inv_blocks)
            inv_stderr = np.std(inv_blocks, ddof=1) / math.sqrt(len(inv_blocks))
            option += '%6i: %14s %10.4g %10.4g 1E3/ %10.4g %10.4g\n' % (
                i, types[i], ave, stderr, inv_ave, inv_stderr)

    print(option, end='')

    if args.noplot:
        return

    print('Select the property to plot, or input any letter to quit:')
    while True:
        plottype = input()
        if not plottype.isdigit():
            sys.exit()
        plottype = int(plottype)

        if plottype < 1 or plottype >= len(data):
            print('not valid')
        else:
            import matplotlib.pyplot as plt
            plt.plot(data[0], data[plottype])
            plt.xlabel(data[0])
            plt.ylabel(types[plottype])
            plt.show()


if __name__ == '__main__':
    types, data = read_log(args.input)
    plot_data(types, data)
