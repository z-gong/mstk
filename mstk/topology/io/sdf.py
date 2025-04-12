from rdkit import Chem
from mstk.topology.atom import Atom
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology
from mstk.chem.element import Element


class Sdf:
    '''
    Generate Topology from SDF file.

    RDKit is used to parse SDF file. Atom symbol, positions, formal charges and bonds are parsed.
    Data fields are ignored.

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
    >>> sdf = Sdf('input.sdf')
    >>> topology = sdf.topology
    '''

    def __init__(self, file, **kwargs):
        suppl = Chem.SDMolSupplier(file, removeHs=False)
        mols = [Molecule.from_rdmol(rdmol) for rdmol in suppl]

        self.topology = Topology(mols)


Topology.register_format('.sdf', Sdf)
