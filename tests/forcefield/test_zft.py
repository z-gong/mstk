#!/usr/bin/env python3

from mstools.forcefield import Padua, Ppf, Zfp, ZftTyper
from mstools.forcefield.ffterm import *
from mstools.topology import *

import os

cwd = os.path.dirname(os.path.abspath(__file__))

im41 = Molecule.from_smiles('C[n+]1cn(cc1)CCCCCC')
dca = Molecule.from_smiles('N#C[N-]C#N')
im2eben = Molecule.from_smiles('C[n+]1cn(cc1)CCc1ccccc1')
top = Topology([im41, dca, im2eben])

typer = ZftTyper(cwd + '/files/CLP-define.zft')
typer.type(top)

top.write(cwd + '/files/zft.xyz')
