import numpy as np
import pandas as pd
from io import StringIO
from . import Trajectory, Frame
from .topology import Atom, Molecule, Topology
from ..forcefield import Element


class LammpsData(Topology):
    def __init__(self, file):
        super(LammpsData, self).__init__()
        self._file = open(file)
        self._file.readline()
        self._file.readline()
        self.n_atom = int(self._file.readline().split()[0])
        self.n_bond = int(self._file.readline().split()[0])
        self.n_angle = int(self._file.readline().split()[0])
        self.n_dihedral = int(self._file.readline().split()[0])
        self.n_atom_type = int(self._file.readline().split()[0])

        self.atoms = [Atom() for i in range(self.n_atom)]
        self._type_masses = [0.] * (self.n_atom_type + 1)
        self._type_names = [''] * (self.n_atom_type + 1)
        self._molecule_dict = {}  # because the length of molecules is unknown. So use dict instead of list

        while True:
            line = self._file.readline()
            if line == '':
                break
            if line.startswith('Masses'):
                self.parse_masses()
            if line.startswith('Atoms'):
                self.parse_atoms()
            if line.startswith('Bonds'):
                self.parse_bonds()
                break

    def parse_masses(self):
        self._file.readline()
        for i in range(self.n_atom_type):
            words = self._file.readline().split()
            type_id = int(words[0])
            self._type_masses[type_id] = float(words[1])
            if words[2] == '#':
                if words[-1] == 'DP':
                    self._type_names[type_id] = 'DP' # Drude particles
                else:
                    self._type_names[type_id] = words[3]
            else:
                self._type_names[type_id] = str(type_id)

    def parse_atoms(self):
        self._file.readline()
        for i in range(self.n_atom):
            words = self._file.readline().split()
            atom_id = int(words[0])
            mol_id = int(words[1])
            type_id = int(words[2])
            charge = float(words[3])

            self.n_molecule = max(self.n_molecule, mol_id)
            if mol_id not in self._molecule_dict.keys():
                mol_name = words[9] if words[7] == '#' else 'UNK'
                mol = Molecule(mol_name)
                mol.id = mol_id - 1  # mol.id starts from 0
                self._molecule_dict[mol_id] = mol
            else:
                mol = self._molecule_dict[mol_id]

            atom = self.atoms[atom_id - 1]
            atom.id = atom_id - 1  # atom.id starts from 0
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = self._type_masses[type_id]
            atom.type = self._type_names[type_id]
            atom.element = Element.guess_from_atom_type(atom.type).symbol

        self.molecules = [self._molecule_dict[i + 1] for i in range(self.n_molecule)]
        for mol in self.molecules:
            for i, atom in enumerate(mol.atoms):
                atom.name = atom.element + str(i + 1)  # atomic symbol + index inside mol starting from 1

    def parse_bonds(self):
        pass


class LammpsTrj(Trajectory):
    '''
    Read step, box and atomic positions (and charges optionally) from dump file of LAMMPS
    Because the topology information are detailed in data file, the mol, type, element in dump file will be ignored
    '''

    def __init__(self, trj_file, mode='r'):
        super().__init__()
        self._file = open(trj_file, mode)
        self._get_info()

    def _get_info(self):
        '''
        Read the number of atoms and record the offset of lines and frames,
        so that we can read arbitrary frame later
        '''
        try:
            self._file.readline()
            self._file.readline()
            self._file.readline()
            self.n_atom = int(self._file.readline())
        except:
            print('Invalid lammpstrj file')
            raise
        self._file.seek(0)

        # read in the file once and build a list of line offsets
        self._line_offset = []
        offset = 0
        for line in self._file:
            self._line_offset.append(offset)
            offset += len(line)
        # the last element is the length of whole file
        self._line_offset.append(offset)
        self._file.seek(0)

        # build a list of frame offsets
        self.n_frame = len(self._line_offset) // (9 + self.n_atom)
        self._frame_offset = []
        for i in range(self.n_frame + 1):
            line_start = (9 + self.n_atom) * i
            self._frame_offset.append(self._line_offset[line_start])

    def read_frame(self, i_frame):
        return self.read_frames([i_frame])[0]

    def read_frames(self, i_frames: [int]) -> [Frame]:
        frames = []
        for i in i_frames:
            # skip to frame i and read only this frame
            self._file.seek(self._frame_offset[i])
            string = self._file.read(self._frame_offset[i + 1] - self._frame_offset[i])
            frames.append(self.read_frame_from_string(string))

        return frames

    def read_frame_from_string(self, string: str):
        frame = Frame(self.n_atom)
        lines = string.splitlines()
        frame.step = int(lines[1])
        frame.xlo, frame.xhi = tuple(map(float, lines[5].split()))
        frame.ylo, frame.yhi = tuple(map(float, lines[6].split()))
        frame.zlo, frame.zhi = tuple(map(float, lines[7].split()))
        frame.box = np.array([frame.xhi - frame.xlo, frame.yhi - frame.ylo, frame.zhi - frame.zlo])
        tokens = lines[8].split()
        title = tokens[2:]

        data = '\n'.join(lines[9:])
        df = pd.read_csv(StringIO(data), header=None, index_col=None, names=title, sep='\s+')
        wrapped = False
        if 'x' in df.columns:
            wrapped = True
        elif 'xs' in df.columns:
            df['x'] = df.xs * (frame.xhi - frame.xlo) + frame.xlo
            df['y'] = df.ys * (frame.yhi - frame.ylo) + frame.ylo
            df['z'] = df.zs * (frame.zhi - frame.zlo) + frame.zlo
            wrapped = True
        elif 'xu' in df.columns:
            df['x'] = df.xu
            df['y'] = df.yu
            df['z'] = df.zu
        elif 'xsu' in df.columns:
            df['x'] = df.xsu * (frame.xhi - frame.xlo) + frame.xlo
            df['y'] = df.ysu * (frame.yhi - frame.ylo) + frame.ylo
            df['z'] = df.zsu * (frame.zhi - frame.zlo) + frame.zlo
        if wrapped:
            if 'ix' not in df.columns:
                print('warning: image flag not found for wrapped positions')
            else:
                df.x += df['ix'] * frame.box[0]  # ix is revered words for pandas, so use df['ix'] instead of df.ix
                df.y += df['iy'] * frame.box[1]
                df.z += df['iz'] * frame.box[2]

        frame.has_charge = 'q' in df.columns
        for row in df.itertuples():
            frame.positions[row.id - 1] = np.array([row.x, row.y, row.z])
            if frame.has_charge:
                frame.charges[row.id - 1] = float(row.q)

        return frame

    def close(self):
        self._file.close()
