#!/usr/bin/env python3

import os
import argparse
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--begin', type=float, default=0)
    parser.add_argument('-p', '--plot', action='store_true', help='plot the histogram')
    parser.add_argument('--skip', type=int, default=1)
    parser.add_argument('--nbin', type=int, help='number of bins for histogram')
    parser.add_argument('--type', type=str, default='openmm', choices=['openmm', 'lammps', 'plumed'])

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    BEGIN = args.begin
    SKIP = args.skip

    r0_k_list = []
    files = filter(os.path.isdir, os.listdir(os.getcwd()))
    for i in files:
        if os.path.isdir(i) and '-' in i and \
                (i[0].isdigit() or i[0].startswith('_')) and i[-1].isdigit():
            r0 = i.split('-')[0]
            k = i.split('-')[1]
            r0_k_list.append((r0, k))
    r0_k_list.sort(key=lambda x: float(x[0].replace('_', '-')))

    outf_meta = open('wham-meta', 'w')
    for r0, k in r0_k_list:
        outfile = 'wham-' + r0 + '-' + k + '.txt'
        outf = open(outfile, 'w')
        if args.type == 'openmm':
            _START = False
            with open(r0 + '-' + k + '/_job.out') as f:
                lines = f.read().splitlines()
            for line in lines:
                if line.startswith('#"'):
                    _START = True
                    continue
                if not _START:
                    continue
                try:
                    words = line.strip().split()
                    step = int(words[0])
                except:
                    continue
                if step >= BEGIN and step % SKIP == 0:
                    str = line.strip().split()
                    time = float(words[0])
                    r = float(words[-1])
                    outf.write('%.1f %.6f\n' % (time, r))
        elif args.type == 'lammps':
            _CYCLE = 0
            _START = False
            countline = 0
            with open(r0 + '-' + k + '/log') as f:
                lines = f.read().splitlines()
            for line in lines:
                if line.startswith('Loop') and _CYCLE == 2:
                    break
                if line.startswith('Step'):
                    _CYCLE += 1
                    if _CYCLE == 2:
                        _START = True
                if _START:
                    countline += 1
                    if countline > BEGIN + 1 and countline % SKIP == 0:
                        str = line.strip().split()
                        time = float(str[0])
                        r = float(str[6])
                        outf.write('%.1f %.6f\n' % (time, r))
        elif args.type == 'plumed':
            countline = 0
            with open(r0 + '-' + k + '/DIS.txt') as f:
                lines = f.read().splitlines()
            for line in lines:
                countline += 1
                if countline > BEGIN + 1 and countline % SKIP == 0:
                    str = line.strip().split()
                    time = float(str[0])
                    r = float(str[1])
                    outf.write('%.1f %.6f\n' % (time, r))
        outf.close()
        outf_meta.write('%s %s %s\n' % (outfile, r0.replace('_', '-'), k))  # replace the leading _ with -
    outf_meta.close()

    if args.plot:
        import matplotlib

        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig = plt.gcf()
        fig.set_size_inches(16, 8)
        plt.grid(True)
        files = filter(lambda x: x.startswith('wham-') and x.endswith('.txt'), os.listdir(os.getcwd()))
        for r0, k in r0_k_list:
            array = np.genfromtxt('wham-%s-%s.txt' % (r0, k))
            data = array[:, 1]
            if args.nbin:
                hist, bins = np.histogram(data, bins=args.nbin)
            else:
                hist, bins = np.histogram(data)
            center = (bins[:-1] + bins[1:]) / 2.
            plt.plot(center, hist, linewidth=2)
        plt.savefig('wham-hist.png')

    # os.system('wham %f %f %s %s 298.15 0 wham-meta wham-pmf' %(hist_min, hist_max, bins, tol))
