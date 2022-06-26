#!/usr/bin/env python3
# coding=utf-8

import argparse
import sys
import math
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('input', type=str, help='data file')
parser.add_argument('-b', '--begin', default=-1, type=float, help='begin from this step')
parser.add_argument('-e', '--end', default=-1, type=float, help='end at this step')
parser.add_argument('-c', '--converge', default=False, action='store_true', help='detect convergence')
parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot data')
parser.add_argument('-r', '--reciprocal', default=False, action='store_true', help='calculate reciprocal')
parser.add_argument('-d', '--dist', default=False, action='store_true', help='plot distribution')
args = parser.parse_args()


def read_log(log_file):
    labels: [str] = []
    data_list: [[float]] = []
    _START = False
    for line in open(log_file):
        if line.startswith('#"'):
            if not _START:
                labels = [s.strip('"') for s in line.strip('#').strip().split('\t')]
                data_list = [[] for _ in labels]
                _START = True
            continue
        if _START:
            words = line.strip().split()
            try:
                step = float(words[0])
            except:
                continue
            if step < args.begin:
                continue
            if args.end > 0 and step > args.end:
                break
            for i in range(len(labels)):
                data_list[i].append(float(words[i]))
    return labels, data_list


def detect_converge(data_list, when_list):
    import pandas as pd
    from mstk.analyzer.series import is_converged
    for i in range(1, len(data_list)):
        series = pd.Series(data_list[i], index=list(range(len(data_list[0]))))
        converged, when = is_converged(series, frac_min=0)
        when_list[i] = when


def average_of_blocks(l, n_block=5):
    ave_block = []
    var_block = []
    bsize = int(math.ceil(len(l) / n_block))
    for i in range(n_block):
        block = l[i * bsize:(i + 1) * bsize]
        ave_block.append(np.mean(block))
        var_block.append(np.var(block))
    return np.array(ave_block), np.array(var_block)


def block_average(l, n_block=5):
    ave_block, var_block = average_of_blocks(l, n_block)
    return np.mean(ave_block), np.std(ave_block, ddof=1) / math.sqrt(n_block), \
           np.mean(var_block), np.std(var_block, ddof=1) / math.sqrt(n_block)


def show_data(types, data_list, when_list):
    option = 'File: %s, Steps: %i-%i, Samples: %i\n' % (
        args.input, data_list[0][0], data_list[0][-1], len(data_list[0]))
    for i in range(1, len(types)):
        when = when_list[i]
        data = data_list[i][when:]
        ave, stderr, var_block, var_stderr = block_average(data)
        var = np.var(data)
        if not args.reciprocal:
            option += '%6i: %14s %10.5g %6.2g %10.4g %6.2g %8.4g\n' % (
                i, types[i], ave, stderr, var, var_stderr, data_list[0][when])
        else:
            ave_block, var_block = average_of_blocks(data)
            inv_blocks = 1000 / ave_block
            inv_ave = inv_blocks.mean()
            inv_stderr = inv_blocks.std(ddof=1) / math.sqrt(len(inv_blocks))
            option += '%6i: %14s %10.5g %6.2g 1E3/ %10.4g %6.2g %8.4g\n' % (
                i, types[i], ave, stderr, inv_ave, inv_stderr, data_list[0][when])

    print(option, end='')


def plot_data(types, data, when_list):
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
            when = when_list[plottype]
            if not args.dist:
                fig, (ax1) = plt.subplots(1, 1, figsize=(5, 4))
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
                ax2.hist(data[plottype][when:], density=True, bins=40, color='C1')
                ax2.set_xlabel(types[plottype])
                ax2.set_ylabel('Probability')
            ax1.plot(data[0][:when], data[plottype][:when])
            ax1.plot(data[0][when:], data[plottype][when:])
            ax1.set_xlabel(types[0])
            ax1.set_ylabel(types[plottype])
            fig.tight_layout()
            plt.savefig(types[plottype] + '.png')
            plt.show()


if __name__ == '__main__':
    labels, data = read_log(args.input)
    when_list = [0] * len(data)
    if args.converge:
        detect_converge(data, when_list)
    show_data(labels, data, when_list)
    if args.plot:
        plot_data(labels, data, when_list)
