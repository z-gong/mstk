#!/usr/bin/env python3

from mstools.forcefield import PaduaFFSet, PpfFFSet, ZfpFFSet
from mstools.forcefield.ffterm import *

import os

cwd = os.path.dirname(os.path.abspath(__file__))
params = PaduaFFSet(cwd + '/files/clp.ff')
params_il = PpfFFSet(cwd + '/files/TEAM_IL.ppf')

ZfpFFSet.write(params, cwd + '/files/out-clp.zfp')
ZfpFFSet.write(params_il, cwd + '/files/out-TEAM_IL.zfp')

zfp = ZfpFFSet(cwd + '/files/out-TEAM_IL.zfp')
term = DrudeTerm('c_4', 1.1, 2.6)
zfp.drude_terms[term.name] = term
ZfpFFSet.write(zfp, cwd + '/files/out-TEAM_IL-2.zfp')
