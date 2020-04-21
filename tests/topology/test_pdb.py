#!/usr/bin/env python3

import os
import pytest
from mstools.topology import Topology

cwd = os.path.dirname(os.path.abspath(__file__))


def test_read():
    pass


def test_write():
    zmat = Topology.open(cwd + '/files/im11.zmat')
    zmat.write(cwd + '/files/zmat-out.pdb')
