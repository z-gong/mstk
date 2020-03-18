#!/usr/bin/env python3

import sys
from mstools.topology import Topology
from mstools.trajectory import Trajectory

if sys.argv[2] == sys.argv[1]:
    raise Exception('output filename should be different from input')

psf = Topology.open(sys.argv[1])
gro = Trajectory.open('conf.gro')
frame = gro.read_frame(0)
psf.set_positions(frame.positions)

psf_out = Topology.open(sys.argv[2], 'w')
psf_out.is_drude = psf.is_drude
psf_out.init_from_molecules(psf.molecules, deepcopy=True)
for atom in psf_out.atoms:
    if atom.type=='SMo' and -0.1 < atom.position[2] < 0.1:
        print(atom, atom.position, atom.charge)
        atom.charge += 0.054 / 100 * float(sys.argv[3])
for atom in psf_out.atoms:
    if atom.type=='SMo' and 6.4 < atom.position[2] < 6.6:
        print(atom, atom.position, atom.charge)
        atom.charge -= 0.054 / 100 * float(sys.argv[3])

psf_out.write()
psf_out.close()

