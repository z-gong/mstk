import os
import tempfile
from .typer import Typer
from ..errors import *
from ...wrapper import DFF

try:
    import openbabel as ob
    import pybel
except:
    OPENBABEL_IMPORTED = False
else:
    OPENBABEL_IMPORTED = True


class DffTyper(Typer):
    '''
    DffTyper calls pre-installed DFF to assign the atom types based on the chemical environment.

    A type definition file is required by DFF for calculating the atom type.
    The native definition files shipped with DFF can be used by providing the full path of one of those files.

    Parameters
    ----------
    dff : DFF
    file : str
        Path of the DFF type definition file
    '''

    def __init__(self, dff, file):
        if not OPENBABEL_IMPORTED:
            raise Exception('Cannot import openbabel')

        self.dff = dff
        self.file = file

    def type_molecule(self, molecule):
        '''
        Assign atom types in the molecule with rules in type definition file.

        The :attr:`~mstools.topology.Atom.type` attribute of all atoms in the molecule will be updated.

        DffTyper use DFF to analyzing the molecular structure, which requires correct bond orders information.
        It expects :attr:`~mstools.topology.Molecule.obmol` attribute to be available in the molecule.
        Usually it means the molecule should be initialized from SMILES or pybel Molecule
        with :func:`~mstools.topology.Molecule.from_smiles` or :func:`~mstools.topology.Molecule.from_pybel`

        If :attr:`~mstools.topology.Molecule.obmol` attribute is None, an Exception will be raised.

        Parameters
        ----------
        molecule : Molecule
        '''
        if molecule.obmol is None:
            raise TypingNotSupportedError('obmol attribute not found for %s' % str(molecule))

        tmpdir = tempfile.mkdtemp()
        mol2 = os.path.join(tmpdir, 'mol.mol2')
        msd = os.path.join(tmpdir, 'mol.msd')
        convert = os.path.join(tmpdir, 'dff_convert')
        setfc = os.path.join(tmpdir, 'dff_setfc')
        typing = os.path.join(tmpdir, 'dff_typing')

        m = pybel.Molecule(molecule.obmol)
        m.write('mol2', mol2)

        self.dff.convert_model_to_msd(mol2, msd, convert)
        self.dff.set_formal_charge([msd], setfc)
        self.dff.typing([msd], self.file, typing)

        from ...topology import Topology
        top = Topology.open(msd)
        for i, atom in enumerate(molecule.atoms):
            atom.type = top.atoms[i].type
