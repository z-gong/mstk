#!/usr/bin/env python3
# coding=utf-8

import argparse
import sys
import re
from mstk.analyzer.series import block_average


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', type=str, help='log file of OpenMM or LAMMPS')
    parser.add_argument('-b', '--begin', default=-1, type=float, help='begin from this sample')
    parser.add_argument('-e', '--end', default=-1, type=float, help='end at this sample')
    parser.add_argument('-c', '--converge', default=False, action='store_true', help='detect convergence')
    parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot data')
    parser.add_argument('-d', '--dist', default=False, action='store_true', help='plot distribution')
    parser.add_argument('--type', type=str, default='openmm', choices=['openmm', 'lammps'])

    return parser.parse_args()


def read_openmm_log(log_file):
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
            if len(words) != len(labels):
                continue
            try:
                values = list(map(float, words))
            except ValueError:
                continue
            step = values[0]
            if step < args.begin:
                continue
            if args.end > 0 and step > args.end:
                break
            for i in range(len(labels)):
                data_list[i].append(values[i])

    return labels, data_list


def read_lammps_log(log_file):
    with open(log_file) as f:
        content = f.read()
    nMin = len(re.findall('\nMinimization', content))
    nRun = len(re.findall('\nStep', content)) - nMin
    print('%i Min, %i Run' % (nMin, nRun))

    labels: [str] = []
    data_list: [[float]] = [[]]
    _START = False
    i_cycle = -nMin
    for line in content.splitlines():
        if line.startswith('Step'):
            i_cycle += 1
            if i_cycle > 0:
                _START = True
            if i_cycle == 1:
                labels = line.strip().split()
                data_list = [[] for i in labels]
            continue
        if i_cycle > 0 and line.startswith('Loop'):
            _START = False
            continue
        if i_cycle > 0 and _START:
            try:
                step = int(line.strip().split()[0])
            except:
                continue
            if step < args.begin:
                continue
            if 0 < args.end < step:
                break
            for i in range(len(labels)):
                data_list[i].append(float(line.strip().split()[i]))
    return labels, data_list


def detect_converge(data_list, when_list):
    import pandas as pd
    from mstk.analyzer.series import is_converged
    for i in range(1, len(data_list)):
        series = pd.Series(data_list[i], index=list(range(len(data_list[0]))))
        converged, when = is_converged(series, frac_min=0)
        when_list[i] = when


def show_data(types, data_list, when_list):
    string = 'File: %s, Steps: %i-%i, Samples: %i\n' % (
        args.input, data_list[0][0], data_list[0][-1], len(data_list[0]))
    string += ' %4s %12s %10s %8s %8s %8s\n' % ('ID', 'LABEL', 'MEAN', 'STDERR', 'STDEV', 'WHEN')
    for i in range(1, len(types)):
        when = when_list[i]
        data = data_list[i][when:]
        (ave, err_ave), (std, err_std) = block_average(data)
        string += ' %4i %12s %10.5g %8.2g %8.2g %8.4g\n' % (
            i, types[i], ave, err_ave, std, data_list[0][when])

    print(string, end='')


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
    args = parse_args()
    if args.type == 'openmm':
        labels, data = read_openmm_log(args.input)
    elif args.type == 'lammps':
        labels, data = read_lammps_log(args.input)
    else:
        raise Exception('Invalid type')

    when_list = [0] * len(data)
    if args.converge:
        detect_converge(data, when_list)
    show_data(labels, data, when_list)
    if args.plot:
        plot_data(labels, data, when_list)
