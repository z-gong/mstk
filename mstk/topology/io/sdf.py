from rdkit import Chem
from mstk.topology.molecule import Molecule
from mstk.topology.topology import Topology


class Sdf:
    '''
    Generate Topology from SDF file.

    RDKit is used to parse SDF file. Atom symbol, positions, formal charges and bonds are parsed.
    Data fields are ignored.

    SDF file can cantain multiple entries. Each entry is parsed as a separate molecule.
    The generated topology contains all molecules (entries) in the SDF file.

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

    @staticmethod
    def save_to(top, file, **kwargs):
        '''
        Save topology into an SDF file.

        Each molecule in the topology is saved as a separate entry in the SDF file.

        Parameters
        ----------
        top : Topology
        file : str
        '''
        writer = Chem.SDWriter(file)
        writer.SetForceV3000(True)
        writer.SetKekulize(True)
        for mol in top.molecules:
            rdmol = mol.rdmol
            if mol.has_position:
                rdmol.RemoveAllConformers()
                rdmol.AddConformer(Chem.Conformer(mol.n_atom))
                rdmol.GetConformer().SetPositions(mol.positions * 10)  # nm-> A
            writer.write(rdmol)
        writer.close()


Topology.register_format('.sdf', Sdf)
