#!/usr/bin/env python3

import sys
import argparse

from mstools.forcefield import ForceField

parser = argparse.ArgumentParser()
parser.add_argument('input', nargs='+', type=str, help='force field files')
parser.add_argument('-o', '--output', required=True, type=str, help='output file')
args = parser.parse_args()

ff = ForceField.open(*args.input)
print('\nForce field info:\n'
      '%6i atom types\n'
      '%6i self vdW terms\n'
      '%6i pairwise vdW terms\n'
      '%6i charge increment terms\n'
      '%6i bond terms\n'
      '%6i angle terms\n'
      '%6i dihedral terms\n'
      '%6i improper terms\n'
      '%6i polarizable terms' % (
          len(ff.atom_types), len(ff.vdw_terms), len(ff.pairwise_vdw_terms),
          len(ff.charge_increment_terms), len(ff.bond_terms), len(ff.angle_terms),
          len(ff.dihedral_terms), len(ff.improper_terms), len(ff.polarizable_terms)
      ))

ff.write(args.output)
