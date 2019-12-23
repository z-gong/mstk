import sys
import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})

from mstools.utils import histogram
from mstools.trajectory.topology import Atom, Molecule, Topology
from mstools.trajectory.lammpsdata import LammpsData
from mstools.trajectory.lammpstrj import LammpsTrj

top = LammpsData('c8c1im-dca/data.lmp')
print(f'{top.n_atom} atoms, {top.n_molecule} molecules')

trj = LammpsTrj('c8c1im-dca/5V.lammpstrj')
print(f'{trj.n_atom} atoms, {trj.n_frame} frames')


def get_atoms(mol: Molecule, names: [str]):
    return [list(filter(lambda x: x.name == name, mol.atoms))[0] for name in names]


def get_com_position(positions, atoms):
    atom_masses = np.array([atom.mass for atom in atoms])
    com_mass = sum(atom_masses)
    atom_ids = [atom.id for atom in atoms]
    atom_positions = positions[atom_ids]
    com_position = np.sum(atom_positions * np.transpose(np.array([atom_masses] * 3)), axis=0) / com_mass
    return com_position


ring_atoms_list = []
dca_atoms_list = []

frame = trj.read_frame(0)
positions = frame.positions
for mol in top.molecules:
    if mol.name == 'c8c1im+':
        ring_atoms = get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'C9', 'C11', 'C14'])
        ring_com = get_com_position(positions, ring_atoms)
        if ring_com[2] > 64:
            ring_atoms_list.append(ring_atoms)

    if mol.name == 'dca-':
        dca_atoms = get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'N9'])
        dca_com = get_com_position(positions, dca_atoms)
        if dca_com[2] < 5:
            dca_atoms_list.append(dca_atoms)

t_list = []
ring_com_z_list = [[] for i in ring_atoms_list]
dca_com_z_list = [[] for i in dca_atoms_list]

for i in range(0, trj.n_frame):
    frame = trj.read_frame(i)
    sys.stdout.write('\r    step %i' % frame.step)

    t_list.append(frame.step / 1000) # ps
    for i, ring_atoms in enumerate(ring_atoms_list):
        ring_com = get_com_position(frame.positions, ring_atoms)
        ring_com_z_list[i].append(ring_com[2])
    for i, dca_atoms in enumerate(dca_atoms_list):
        dca_com = get_com_position(frame.positions, dca_atoms)
        dca_com_z_list[i].append(dca_com[2])

print('')

fig = plt.figure()
ax = fig.add_subplot('111')
for z_list in ring_com_z_list:
    ax.plot(t_list, z_list)
fig.savefig('5V-diffusion-ring.png')
plt.cla()

for z_list in dca_com_z_list:
    ax.plot(t_list, z_list)
fig.savefig('5V-diffusion-dca.png')
plt.cla()
