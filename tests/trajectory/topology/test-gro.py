#!/usr/bin/env python3

import numpy as np
from mstools.trajectory.gro import GRO

gro = GRO('test-gro.gro')
assert gro.n_atom == 20680
assert gro.n_frame == 2

frame = gro.read_frame(0)
assert frame.has_velocity == False
assert np.all(frame.positions[0] == np.array([-0.121, 0.423, -0.771]))

frame = gro.read_frame(1)
assert frame.xhi == 3.170
assert frame.yhi == 3.296
assert frame.zhi == 13.900
assert np.all(frame.positions[-1] == np.array([2.821, 0.480, 11.564]))
