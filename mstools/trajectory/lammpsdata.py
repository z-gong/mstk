from .topology import Atom, Molecule, Topology, atomic_symbol


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
        self.type_masses = [0.] * (self.n_atom_type + 1)
        self.type_names = [''] * (self.n_atom_type + 1)
        self.n_molecule = 0
        self.molecule_dict = {}  # because i don't know the length of molecules. So use dict instead of list
        self.molecules: [Molecule]

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
            self.type_masses[type_id] = float(words[1])
            if words[2] == '#':
                self.type_names[type_id] = '_'.join(words[3:])
            else:
                self.type_names[type_id] = str(type_id)

    def parse_atoms(self):
        self._file.readline()
        for i in range(self.n_atom):
            words = self._file.readline().split()
            atom_id = int(words[0])
            mol_id = int(words[1])
            type_id = int(words[2])
            charge = float(words[3])

            self.n_molecule = max(self.n_molecule, mol_id)
            if mol_id not in self.molecule_dict.keys():
                mol_name = words[9] if words[7] == '#' else 'UNK'
                mol = Molecule(mol_name)
                mol.id = mol_id - 1  # mol.id starts from 0
                self.molecule_dict[mol_id] = mol
            else:
                mol = self.molecule_dict[mol_id]

            atom = self.atoms[atom_id - 1]
            atom.id = atom_id - 1  # atom.id starts from 0
            mol.add_atom(atom)
            atom.charge = charge
            atom.mass = self.type_masses[type_id]
            atom.type = self.type_names[type_id]

        self.molecules = [self.molecule_dict[i + 1] for i in range(self.n_molecule)]
        for mol in self.molecules:
            for i, atom in enumerate(mol.atoms):
                if atom.name == 'UNK':
                    atom.name = atomic_symbol(atom.type) + str(i + 1) # atomic symbol + index inside mol starting from 1

    def parse_bonds(self):
        pass
