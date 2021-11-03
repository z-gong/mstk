#!/usr/bin/env python3

import sys
import argparse
import multiprocessing
import subprocess as sp
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--nxvg', type=int, required=True, help='Number of xvg files')
parser.add_argument('--time', type=float, required=True, help='Time length for calculating ACF')
parser.add_argument('--nproc', type=int, default=1, help='Number of parallel processes')
parser.add_argument('--nomp', type=int, default=1, help='Number of OpenMP threads')
parser.add_argument('--gro', type=str, help='Gro file used to get volume of box')

args = parser.parse_args()

n_xvg = args.nxvg
corr_time = args.time
n_proc = args.nproc
n_omp = args.nomp

vol = 1.0
if args.gro is not None:
    with open(args.gro) as f:
        lines = f.read().splitlines()
        words = lines[-1].strip().split()
        box = list(map(float, words))
        vol = box[0] * box[1] * box[2]


def calc_acf(xvg, acf, vis, corr_time):
    cmd = 'OMP_NUM_THREADS=%i vis-gk %s %s %s %i' % (n_omp, xvg, acf, vis, corr_time)
    sp.check_call(cmd, shell=True)

def wrapper(args):
    return calc_acf(*args)

args_list = []
for i in range(n_xvg):
    xvg = 'energy%i.xvg' % i
    acf = 'acf%i.txt' % i
    vis = 'vis%i.txt' % i
    args_list.append([xvg, acf, vis, corr_time])

if n_proc > 0:
    print('Calculating ACF...')
    with multiprocessing.Pool(n_proc) as pool:
        pool.map(wrapper, args_list)

df_vis = pd.DataFrame()
for i in range(n_xvg):
    vis = 'vis%i.txt' % i
    df = pd.read_table(vis, index_col=0, header=None, names=['vis'])
    df_vis['vis%i' % i] = df.vis

ave_series = df_vis.mean(axis=1) * vol
std_series = df_vis.std(axis=1) * vol

with open('GK-ave.txt', 'w')  as f:
    for i in range(len(ave_series)):
        f.write('%f\t%f\t%f\n' % (ave_series.index[i], ave_series.iloc[i], std_series.iloc[i]))
