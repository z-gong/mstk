#!/usr/bin/env python3

import os
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import Ppf, Padua, Zfp
from mstools.simsys import System

import simtk.openmm as mm
from simtk.openmm import app
from simtk.unit import kilocalorie_per_mole as kcal_mol

cwd = os.path.dirname(os.path.abspath(__file__))


ff = Ppf(cwd + '/files/10-benzene.ppf')
top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
top.assign_charge_from_ff(ff)
system = System(top, ff)
system.export_gmx()
