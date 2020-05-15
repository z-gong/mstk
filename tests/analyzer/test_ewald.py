#!/usr/bin/env python3

import pytest
from mstools.analyzer.ewald import EwaldSum

import os

cwd = os.path.dirname(os.path.abspath(__file__))


def test_energy():
    ewald = EwaldSum.create_test(100, seed=0, cutoff=1.2)
    assert pytest.approx(ewald.alpha, abs=1E-4) == 2.1902
    assert ewald.kmax_x == 7
    assert ewald.kmax_y == 9
    assert ewald.kmax_z == 11

    energy, forces = ewald.calc_energy_forces(ewald=True)
    assert pytest.approx(energy, rel=1E-6) == -1089.0515
    assert pytest.approx(forces[-1], rel=1E-4) == [-1.6694e+02, 5.3097e+02, -1.4420e+02]

    energy, forces = ewald.calc_energy_forces(ewald=False)
    assert pytest.approx(energy, rel=1E-6) == -1011.3025
    assert pytest.approx(forces[-1], rel=1E-4) == [-1.5860e+02, 4.9182e+02, -7.5005e+01]

    ewald = EwaldSum.create_test(100, seed=0, cutoff=1.4)
    assert pytest.approx(ewald.alpha, abs=1E-4) == 1.8773
    assert ewald.kmax_x == 5
    assert ewald.kmax_y == 7
    assert ewald.kmax_z == 9

    energy, forces = ewald.calc_energy_forces(ewald=True)
    assert pytest.approx(energy, rel=1E-6) == -1089.1090
    assert pytest.approx(forces[-1], rel=1E-4) == [-1.6675e+02, 5.3089e+02, -1.4421e+02]

    energy, forces = ewald.calc_energy_forces(ewald=False)
    assert pytest.approx(energy, rel=1E-6) == -1011.3025
    assert pytest.approx(forces[-1], rel=1E-4) == [-1.5860e+02, 4.9182e+02, -7.5005e+01]
