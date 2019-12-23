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

trj = LammpsTrj('c8c1im-dca/15V.lammpstrj')
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


bins = np.linspace(0, 70, 701)
dz = bins[1] - bins[0]
z_array = (bins[1:] + bins[-1]) / 2
charge_array = np.array([0.] * len(z_array))
z_ring_list = []
z_ring_com_list = []
z_tail_list = []
z_dca_list = []
for i in range(30, trj.n_frame):
    frame = trj.read_frame(i)
    sys.stdout.write('\r    step %i' % frame.step)
    positions = frame.positions
    for mol in top.molecules:
        if mol.name == 'c8c1im+':
            ring_atoms = get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'C9', 'C11', 'C14'])
            z_ring_list += [positions[atom.id][2] for atom in ring_atoms]
            ring_com = get_com_position(positions, ring_atoms)
            z_ring_com_list.append(ring_com[2])
            tail_atoms = get_atoms(mol, ['C21', 'C25', 'C29', 'C33', 'C37', 'C41', 'C45'])
            z_tail_list += [positions[atom.id][2] for atom in tail_atoms]
        if mol.name == 'dca-':
            dca_atoms = get_atoms(mol, ['N1', 'C3', 'N5', 'C7', 'N9'])
            z_dca_list += [positions[atom.id][2] for atom in dca_atoms]

    for atom in top.atoms:
        if atom.molecule.name not in ['c8c1im+', 'dca-']:
            continue
        z = frame.positions[atom.id][2]
        q = frame.charges[atom.id]
        z_idx = math.floor((z - bins[0]) / dz)
        charge_array[z_idx] += q

print('')

fig = plt.figure()
ax = fig.add_subplot('111')
x, y = histogram(z_ring_list, bins=bins)
ax.plot(x, y / frame.box[0] / frame.box[1] / dz / trj.n_frame, label='c8c1im+ ring')
x, y = histogram(z_tail_list, bins=bins)
ax.plot(x, y / frame.box[0] / frame.box[1] / dz / trj.n_frame, label='c8c1im+ tail')
x, y = histogram(z_dca_list, bins=bins)
ax.plot(x, y / frame.box[0] / frame.box[1] / dz / trj.n_frame, label='dca-')
ax.legend()
fig.savefig('15V-dist.png')
plt.cla()

x, y = histogram(z_ring_com_list, bins=bins)
ax.plot(x, y / frame.box[0] / frame.box[1] / dz / trj.n_frame, label='c8c1im+ ring COM')
fig.savefig('15V-dist-com.png')
plt.cla()

ax.plot(z_array, charge_array / frame.box[0] / frame.box[1] / dz / trj.n_frame, label='c8c1im+ ring COM')
fig.savefig('15V-charge.png')
plt.cla()
