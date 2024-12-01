import os
from mstk import DIR_MSTK
from mstk.topology import Bond
from .smarts_typer import SmartsTyper


class GaffTyper(SmartsTyper):
    '''
    `GaffTyper` is specifically designed for GAFF atom type assignment.
    By default, it uses the type definition file `gaff.smt`

    Compared with `SmartsTyper`, this class will handle conjugated atom types according to the GAFF convention,
    e.g. cc/cd, ce/cf etc...

    The official tool for GAFF - antechamber has several issues:
    - Mysterious aromatic assignment
    - A lost of parameters missing. `parmchk2` gives unreasonable guess for torsion parameters, e.g. biphenyl

    Parameters
    ----------
    file : str or file-like object, optional
        Type definition file.

    Notes
    -----
    * SMARTS is parsed by using RDKit package. Make sure it is installed.
    * In type definition file, empty lines are ignored, and comments should start with ##.

    '''

    def __init__(self, file=None):
        file = file or os.path.join(DIR_MSTK, 'data', 'forcefield', 'gaff.smt')
        super().__init__(file)

    def _type_molecule(self, molecule):
        '''
        Assign atom types in the molecule with rules in type definition file.

        The :attr:`~mstk.topology.Atom.type` attribute of all atoms in the molecule will be updated.

        GaffTyper use RDKit to do SMARTS matching, therefore the orders must be set for all the bonds in the molecule.
        Usually it means the molecule be initialized from SMILES.
        with :func:`~mstk.topology.Molecule.from_smiles` or :func:`~mstk.topology.Molecule.from_rdmol`

        If an atom can not match any type by the predefined SMARTS patterns, an Exception will be raised.

        Parameters
        ----------
        molecule : Molecule
        '''
        super()._type_molecule(molecule)

        # find the chain of conjugated atom types, add chose one from the two types for each atom
        conj_atoms = [atom for atom in molecule.atoms if '|' in atom.type]
        if not conj_atoms:
            return

        conj1_partners = {atom: [] for atom in conj_atoms}
        conj2_partners = {atom: [] for atom in conj_atoms}
        for atom in conj_atoms:
            for bond in atom.bonds:
                partner = bond.atom2 if bond.atom1 is atom else bond.atom1
                if partner not in conj_atoms:
                    continue
                if bond.order == Bond.Order.DOUBLE or bond.is_aromatic:
                    conj2_partners[atom].append(partner)
                else:
                    conj1_partners[atom].append(partner)

        flag = {atom: False for atom in conj_atoms}

        def set_conj_type_partners(atom):
            type_suffix = atom.type[-1]
            for partner in conj1_partners[atom]:
                if flag[partner]:
                    continue
                candidates = partner.type.split('|')
                try:
                    partner.type = next(typ for typ in candidates if typ.endswith(type_suffix))
                except StopIteration:
                    partner.type = candidates[0]
                flag[partner] = True
                set_conj_type_partners(partner)
            for partner in conj2_partners[atom]:
                if flag[partner]:
                    continue
                candidates = partner.type.split('|')
                try:
                    partner.type = next(typ for typ in candidates if not typ.endswith(type_suffix))
                except StopIteration:
                    partner.type = candidates[0]
                flag[partner] = True
                set_conj_type_partners(partner)

        for atom in conj_atoms:
            if flag[atom]:
                continue
            atom.type = atom.type.split('|')[0]
            flag[atom] = True
            set_conj_type_partners(atom)
