#!/usr/bin/env python3

import os
import tempfile
import filecmp
import pytest
from mstools.topology import Topology, Molecule, Atom, TIP4PSite
from mstools.forcefield import ForceField

cwd = os.path.dirname(os.path.abspath(__file__))


def test_tip4p():
    top = Topology.open(os.path.join(cwd, 'files/TIP3P.zmat'))
    mol = top.molecules[0]

    ff = ForceField.open(os.path.join(cwd, '../forcefield/files/TIP4P.zfp'))
    mol.generate_virtual_sites(ff)
    mol.assign_charge_from_ff(ff)

    avsite = top.atoms[-1]
    assert avsite.name == 'VS1'
    assert avsite.type == 'MW'
    assert pytest.approx(avsite.position, abs=1E-6) == [0.0091824, -0.011861, 0.]

    assert pytest.approx([atom.charge for atom in top.atoms], abs=1E-6) == [0, 0.52, 0.52, -1.04]
