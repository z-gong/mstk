#!/usr/bin/env python3

import os
import pytest
from mstools.topology import Topology, Molecule

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    pass


def test_write():
    zmat = Topology.open(cwd + '/files/im11.zmat')
    zmat.write(cwd + '/files/zmat-out.pdb')
    mol = Molecule.from_smiles('C[n+]1cn(cc1)CCCC[B-](F)(F)F')
    top = Topology()
    top.init_from_molecules([mol])
    top.write(cwd + '/files/smi-out.pdb')
