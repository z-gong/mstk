from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology


class Smi():
    '''
    Generate Topology from SMI file.

    Each line is a molecule represented by SMILES string.
    Molecule name can optionally be put after the SMILES string, separated by spaces.
    Therefore, the moelcule name itself cannot contain space.
    Lines started with # will be treated as comments.
    The positions will be generated using RDKit.

    Parameters
    ----------
    file : str
    kwargs : dict
        Ignored

    Attributes
    ----------
    topology : Topology

    Examples
    --------
    >>> smi = Smi('input.zmat')
    >>> topology = smi.topology
    '''

    def __init__(self, file, **kwargs):
        self.topology = Topology()
        self._parse(file)

    def _parse(self, file):
        molecules = []
        with open(file) as f:
            for line in f:
                line = line.strip()
                if line == '' or line.startswith('#'):
                    continue
                mol = Molecule.from_smiles(line)
                molecules.append(mol)

        self.topology.update_molecules(molecules)


Topology.register_format('.smi', Smi)
