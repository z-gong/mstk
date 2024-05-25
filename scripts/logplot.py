#!/usr/bin/env python3

import argparse
import os
import re
import pandas as pd
from mstk.analyzer.series import block_average, is_converged
from mstk.analyzer.fitting import polyfit, polyval
from mstk.utils import print_data_to_file


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='XVG file or log file of OpenMM or LAMMPS')
    parser.add_argument('-b', '--begin', default=-1, type=float, help='begin from this data')
    parser.add_argument('-e', '--end', default=-1, type=float, help='end at this data')
    parser.add_argument('--skip', default=1, type=int, help='read every N rows')
    parser.add_argument('-o', '--output', type=str, help='output data to a file')
    parser.add_argument('-c', '--converge', default=False, action='store_true', help='detect convergence')
    parser.add_argument('-p', '--plot', default=False, action='store_true', help='plot data')
    parser.add_argument('--type', type=str, choices=['openmm', 'lammps', 'xvg'])
    parser.add_argument('--fit', action='store_true', help='perform linear fitting')
    parser.add_argument('--dist', action='store_true', help='plot distribution')

    return parser.parse_args()


class Analyzer:
    def __init__(self, log_file, file_type=None, begin=-1, end=-1, skip=1):
        self.labels: [str] = []
        self.data_list: [[float]] = []
        self.when_list: [int] = []
        self.fit_coeff: [(float)] = []
        self.log_file = log_file

        _ext_filetype = {
            'xvg': 'xvg',
            'out': 'openmm',
            'txt': 'openmm',
            'log': 'lammps'
        }
        ext = os.path.splitext(log_file)[1].lstrip('.').lower()
        file_type = file_type or _ext_filetype.get(ext, ext)
        if file_type == 'xvg':
            self.read_xvg(begin, end, skip)
        elif file_type == 'openmm':
            self.read_openmm_log(begin, end, skip)
        elif file_type == 'lammps':
            self.read_lammps_log(begin, end, skip)
        else:
            raise Exception('unknown log file type: ' + file_type)

        # in case times in all frames are zero
        if all(i == 0 for i in self.data_list[0]):
            self.data_list[0] = [i for i in range(len(self.data_list[0]))]

    @staticmethod
    def str2float(string):
        try:
            val = float(string)
        except ValueError:
            val = 0
        return val

    def read_openmm_log(self, begin, end, skip):
        n_row = -1
        _START = False
        for line in open(self.log_file):
            if line.startswith('#"'):
                if not _START:
                    self.labels = [s.strip('"') for s in line.strip('#').strip().split('\t')]
                    self.data_list = [[] for _ in self.labels]
                    self.when_list = [0 for _ in self.labels]
                    _START = True
                continue
            if _START:
                words = line.strip().split()
                if len(words) != len(self.labels):
                    continue
                try:
                    values = list(map(Analyzer.str2float, words))
                except ValueError:
                    continue
                step = values[0]
                if step < begin:
                    continue
                if 0 < end < step:
                    break
                n_row += 1
                if n_row % skip != 0:
                    continue
                for i in range(len(self.labels)):
                    self.data_list[i].append(values[i])

    def read_lammps_log(self, begin, end, skip):
        with open(self.log_file) as f:
            content = f.read()
        nMin = len(re.findall('\nMinimization', content))
        nRun = len(re.findall('\nStep', content)) - nMin
        print('%i Min, %i Run' % (nMin, nRun))

        n_row = -1
        _START = False
        i_cycle = -nMin
        for line in content.splitlines():
            if line.startswith('Step'):
                i_cycle += 1
                if i_cycle > 0:
                    _START = True
                if i_cycle == 1:
                    self.labels = line.strip().split()
                    self.data_list = [[] for _ in self.labels]
                    self.when_list = [0 for _ in self.labels]
                continue
            if i_cycle > 0 and line.startswith('Loop'):
                _START = False
                continue
            if i_cycle > 0 and _START:
                try:
                    step = int(line.strip().split()[0])
                except:
                    continue
                if step < begin:
                    continue
                if 0 < end < step:
                    break
                n_row += 1
                if n_row % skip != 0:
                    continue
                for i in range(len(self.labels)):
                    self.data_list[i].append(float(line.strip().split()[i]))

    def read_xvg(self, begin, end, skip):
        n_row = -1
        for line in open(self.log_file):
            if line.startswith('#'):
                continue
            if line.startswith('@'):
                words = line.strip().split()
                if words[1] in ('xaxis', 'yaxis'):
                    self.labels.append(words[3].strip('"').split()[0])
                    self.data_list.append([])
                    self.when_list.append(0)
                continue
            words = line.strip().split()
            step = float(words[0])
            if step < begin:
                continue
            if 0 < end < step:
                break
            n_row += 1
            if n_row % skip != 0:
                continue
            for i in range(len(self.data_list)):
                self.data_list[i].append(float(words[i]))

    def detect_converge(self):
        for i in range(1, len(self.data_list)):
            series = pd.Series(self.data_list[i], index=list(range(len(self.data_list[0]))))
            converged, when = is_converged(series, frac_min=0)
            self.when_list[i] = when

    def fit(self):
        self.fit_coeff.append((0.0, 1.0))
        for i in range(1, len(self.data_list)):
            when = self.when_list[i]
            coeff, score = polyfit(self.data_list[0][when:], self.data_list[i][when:], degree=1)
            self.fit_coeff.append(coeff)

    def print_summary(self):
        string = 'File: %s, Steps: %i-%i, Samples: %i\n' % (
            self.log_file, self.data_list[0][0], self.data_list[0][-1], len(self.data_list[0]))
        string += ' %2s %12s %11s %8s %8s %8s %9s %9s\n' % (
            'ID', 'LABEL', 'MEAN', 'STDERR', 'STDEV', 'WHEN', 'INTERCEPT', 'SLOPE')
        for i in range(1, len(self.labels)):
            when = self.when_list[i]
            data = self.data_list[i][when:]
            if self.fit_coeff:
                intercept, slope = self.fit_coeff[i]
            else:
                intercept, slope = 0, 0
            (ave, err_ave), (std, err_std) = block_average(data)
            string += ' %2i %12s %11.5g %8.2g %8.2g %8.4g %9.3g %9.3g\n' % (
                i, self.labels[i], ave, err_ave, std, self.data_list[0][when], intercept, slope)

        print(string, end='')

    def plot_data(self, plot_distribution=False):
        import matplotlib
        import matplotlib.pyplot as plt

        matplotlib.use('Agg')
        matplotlib.rcParams.update({'font.size': 15})

        for idx in range(1, len(self.data_list)):
            when = self.when_list[idx]
            if not plot_distribution:
                fig, ax1 = plt.subplots()
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2)
                ax2.hist(self.data_list[idx][when:], density=True, bins=40, color='C1')
                ax2.set_xlabel(self.labels[idx])
                ax2.set_ylabel('Probability')
            ax1.plot(self.data_list[0][:when], self.data_list[idx][:when])
            x_list, y_list = self.data_list[0][when:], self.data_list[idx][when:]
            ax1.plot(x_list, y_list)
            if self.fit_coeff:
                y_pred_list = [polyval(x, self.fit_coeff[idx]) for x in x_list]
                ax1.plot(x_list, y_pred_list, '--')
            ax1.set_xlabel(self.labels[0])
            ax1.set_ylabel(self.labels[idx])
            fig.tight_layout()
            plt.savefig(self.labels[idx] + '.png')
            plt.show()


if __name__ == '__main__':
    args = parse_args()

    analyzer = Analyzer(args.input, args.type, args.begin, args.end, args.skip)
    if args.converge:
        analyzer.detect_converge()
    if args.fit:
        analyzer.fit()
    analyzer.print_summary()
    if args.plot:
        analyzer.plot_data(plot_distribution=args.dist)
    if args.output:
        print_data_to_file(dict(zip(analyzer.labels, analyzer.data_list)), file=args.output)
