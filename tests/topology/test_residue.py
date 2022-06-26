#!/usr/bin/env python3

import pytest
import copy
from mstk.topology import Topology, Molecule


def test_residue():
    im41 = Molecule.from_smiles('C[n+]1cn(cc1)CCCC Im41')
    atom1, atom2, atom3 = im41.atoms[0], im41.atoms[1], im41.atoms[-1]
    assert im41.n_residue == 1

    residue = im41.residues[0]
    assert residue.id == -1
    assert residue.id_in_mol == 0
    assert residue.name == 'Im41'
    assert residue.n_atom == im41.n_atom

    assert atom1.residue.name == 'Im41'
    assert atom2.residue.name == 'Im41'
    assert atom3.residue.name == 'Im41'

    im41.add_residue('IM', im41.atoms[1:6])
    assert im41.n_residue == 2
    res1, res2 = im41.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'Im41'
    assert res1.n_atom == im41.n_atom - 5

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'IM'
    assert res2.n_atom == 5

    assert atom1.residue.name == 'Im41'
    assert atom2.residue.name == 'IM'
    assert atom3.residue.name == 'Im41'

    im41.add_residue('C', im41.atoms[0:1])
    im41.add_residue('TAIL', im41.atoms[6:])
    assert im41.n_residue == 3
    res1, res2, res3 = im41.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'IM'
    assert res1.n_atom == 5

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'C'
    assert res2.n_atom == 1

    assert res3.id == -1
    assert res3.id_in_mol == 2
    assert res3.name == 'TAIL'
    assert res3.n_atom == im41.n_atom - 6

    assert atom1.residue.name == 'C'
    assert atom2.residue.name == 'IM'
    assert atom3.residue.name == 'TAIL'

    im41.remove_residue(res2)
    assert im41.n_residue == 3
    res1, res2, res3 = im41.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'Im41'
    assert res1.n_atom == 1

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'IM'
    assert res2.n_atom == 5

    assert res3.id == -1
    assert res3.id_in_mol == 2
    assert res3.name == 'TAIL'
    assert res3.n_atom == im41.n_atom - 6

    assert atom1.residue.name == 'Im41'
    assert atom2.residue.name == 'IM'
    assert atom3.residue.name == 'TAIL'

    im41 = copy.deepcopy(im41)
    assert im41.n_residue == 3
    res1, res2, res3 = im41.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'Im41'
    assert res1.n_atom == 1

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'IM'
    assert res2.n_atom == 5

    assert res3.id == -1
    assert res3.id_in_mol == 2
    assert res3.name == 'TAIL'
    assert res3.n_atom == im41.n_atom - 6

    assert atom1.residue.name == 'Im41'
    assert atom2.residue.name == 'IM'
    assert atom3.residue.name == 'TAIL'


def test_top_residue():
    ben = Molecule.from_smiles('c1ccccc1')
    sol = Molecule.from_smiles('O    SOL')
    top = Topology([ben, sol])
    assert top.n_residue == 2
    res1, res2 = top.residues
    assert res1.id == 0
    assert res1.id_in_mol == 0
    assert res1.name == 'C6H6'
    assert res1.n_atom == 12

    assert res2.id == 1
    assert res2.id_in_mol == 0
    assert res2.name == 'SOL'
    assert res2.n_atom == 3

    merged = Molecule.merge([ben, sol])
    assert merged.n_residue == 2
    res1, res2 = merged.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'C6H6'
    assert res1.n_atom == 12

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'SOL'
    assert res2.n_atom == 3

    ben.add_residue('C', ben.atoms[:6])
    merged = Molecule.merge([ben, sol])
    assert merged.n_residue == 3
    res1, res2, res3 = merged.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'C6H6'
    assert res1.n_atom == 6

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'C'
    assert res2.n_atom == 6

    assert res3.id == -1
    assert res3.id_in_mol == 2
    assert res3.name == 'SOL'
    assert res3.n_atom == 3

    ben, sol = merged.split()
    assert ben.n_residue == 2
    res1, res2 = ben.residues
    assert res1.id == -1
    assert res1.id_in_mol == 0
    assert res1.name == 'C'
    assert res1.n_atom == 6

    assert res2.id == -1
    assert res2.id_in_mol == 1
    assert res2.name == 'C6H6'
    assert res2.n_atom == 6

    assert sol.n_residue == 1
    res3 = sol.residues[0]
    assert res3.id == -1
    assert res3.id_in_mol == 0
    assert res3.name == 'SOL'
    assert res3.n_atom == 3
