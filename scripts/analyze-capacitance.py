#!/usr/bin/env python3

import argparse
import pandas as pd

from mstools.constant import *
from mstools.analyzer.fitting import polyfit, polyval_derivative

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+', type=str, help='files of voltage profiles')
parser.add_argument('--bulk', type=float, default=2.0, help='length of bulk region (in nm)')
parser.add_argument('--cathode', type=float, default=0.0, help='z coordinate of cathode (in nm)')
parser.add_argument('--anode', type=float, default=6.5, help='z coordinate of anode (in nm)')

args = parser.parse_args()

V_Q_electrode = {}
V_Q_cell = {}
V_zc = 0

for inp in args.input:
    df = pd.read_csv(inp, sep='\t', index_col=0)
    z_array = df.index.array
    dz = z_array[1] - z_array[0]
    z_mid = z_array.mean()

    V_bulk = df.V[z_mid - args.bulk / 2: z_mid + args.bulk / 2].mean()
    V_cathode = df.V.iloc[df.index.get_loc(args.cathode - dz, method='nearest')] - V_bulk
    V_anode = df.V.iloc[df.index.get_loc(args.anode + dz, method='nearest')] - V_bulk
    V_drop = V_cathode - V_anode

    q_cathode = df.rho_q[args.cathode] * dz / NANO ** 2 * ELEMENTARY_CHARGE / MILLI  # mC/m^2
    q_anode = df.rho_q[args.anode] * dz / NANO ** 2 * ELEMENTARY_CHARGE / MILLI  # mC/m^2

    V_Q_electrode[V_cathode] = q_cathode
    # ignore anode if no voltage drop
    if q_cathode == 0:
        V_zc = V_cathode
    else:
        V_Q_electrode[V_anode] = q_anode

    V_Q_cell[V_drop] = q_cathode

V_list, Q_list = zip(*sorted(V_Q_electrode.items()))
coeff4, score4 = polyfit(V_list, Q_list, 4)
coeff5, score5 = polyfit(V_list, Q_list, 5)
print('Zero charge voltage: %.4f V' % V_zc)
print('RSQ for 4th and 5th polynomial fitting Cdiff: %.4f %.4f' % (score4, score5))

print('%10s %10s %10s %10s %10s %10s' % ('electrode', 'V', 'Q', 'Cint', 'Cdiff-4th', 'Cdiff-5th'))
for V, Q in sorted(V_Q_electrode.items()):
    Cdiff4 = polyval_derivative(V, coeff4)[1]
    Cdiff5 = polyval_derivative(V, coeff5)[1]
    Cint = Q / (V - V_zc) if V - V_zc != 0 else 0
    print('%10s %10.4f %10.4f %10.4f %10.4f %10.4f' % ('', V, Q, Cint, Cdiff4, Cdiff5))
print('%10s %10s %10s %10s' % ('cell', 'V', 'Q', 'Cint'))
for V, Q in sorted(V_Q_cell.items()):
    Cint = Q / V if V != 0 else 0
    print('%10s %10.4f %10.4f %10.4f' % ('', V, Q, Q / V))
