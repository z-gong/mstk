#!/usr/bin/env python3

import os
import pytest
from mstk.topology import UnitCell

cwd = os.path.dirname(os.path.abspath(__file__))


def test_cell():
    cell = UnitCell()
    assert cell.is_rectangular == True
    assert cell.volume == 0.0
    assert pytest.approx(cell.size, abs=1E-6) == [0, 0, 0]
    assert pytest.approx(cell.lengths, abs=1E-6) == [0, 0, 0]
    assert pytest.approx(cell.angles, abs=1E-6) == [90, 90, 90]

    cell = UnitCell([0, 0, 0])
    assert cell.is_rectangular == True
    assert cell.volume == 0.0
    assert pytest.approx(cell.size, abs=1E-6) == [0, 0, 0]
    assert pytest.approx(cell.lengths, abs=1E-6) == [0, 0, 0]
    assert pytest.approx(cell.angles, abs=1E-6) == [90, 90, 90]

    cell = UnitCell([1, 2, 3])
    assert cell.is_rectangular == True
    assert cell.volume == 6.0
    assert pytest.approx(cell.size, abs=1E-6) == [1, 2, 3]
    assert pytest.approx(cell.lengths, abs=1E-6) == [1, 2, 3]
    assert pytest.approx(cell.angles, abs=1E-6) == [90, 90, 90]

    cell = UnitCell([[1, 0, 0], [0, 2, 0], [0, 0, 3]])
    assert cell.is_rectangular == True
    assert cell.volume == 6.0
    assert pytest.approx(cell.size, abs=1E-6) == [1, 2, 3]
    assert pytest.approx(cell.lengths, abs=1E-6) == [1, 2, 3]
    assert pytest.approx(cell.angles, abs=1E-6) == [90, 90, 90]

    cell = UnitCell([[1, 0, 0], [0.3, 2, 0], [0.3, 0.5, 3]])
    assert cell.is_rectangular == False
    assert cell.volume == 6.0
    assert pytest.approx(cell.size, abs=1E-6) == [1, 2, 3]
    assert pytest.approx(cell.lengths, abs=1E-6) == [1, 2.022375, 3.056141]
    assert pytest.approx(cell.angles, abs=1E-6) == [79.842392, 84.366600, 81.469230]
