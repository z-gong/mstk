#!/usr/bin/env python3

import os
import sys
import pytest
from mstools.topology import Topology, UnitCell
from mstools.forcefield import ForceField
from mstools.simsys import System

cwd = os.path.dirname(os.path.abspath(__file__))


def test_gmx():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_gromacs(gro_out=cwd + '/files/conf.gro',
                          top_out=cwd + '/files/topol.top',
                          mdp_out=cwd + '/files/grompp.mdp')


def test_gmx_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_gromacs(gro_out=cwd + '/files/conf-drude.gro',
                          top_out=cwd + '/files/topol-drude.top',
                          mdp_out=None)


def test_lmp_drude():
    ff = ForceField.open(cwd + '/../forcefield/files/CLP.ff', cwd + '/../forcefield/files/CLPol-alpha.ff')
    top = Topology.open(cwd + '/files/5-Im21-BF4-drude.lmp')
    top.generate_angle_dihedral_improper()
    top.generate_drude_particles(ff)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_lammps(data_out=cwd + '/files/data-drude.lmp',
                         in_out=cwd + '/files/in-drude.lmp')


def test_charmm():
    ff = ForceField.open(cwd + '/files/10-benzene.ppf')
    top = Topology.open(cwd + '/files/10-benzene.lmp', improper_center=3)
    top.assign_charge_from_ff(ff)
    system = System(top, ff)
    system.export_charmm(pdb_out=cwd + '/files/conf.pdb',
                         psf_out=cwd + '/files/psf.psf',
                         prm_out=cwd + '/files/ff.prm')
