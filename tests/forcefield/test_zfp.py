#!/usr/bin/env python3

from mstools.forcefield import PaduaFFSet, PpfFFSet, ZfpFFSet

import os

cwd = os.path.dirname(os.path.abspath(__file__))
params = PaduaFFSet(cwd + '/files/clp.ff')
params_il = PpfFFSet(cwd + '/files/TEAM_IL.ppf')

ZfpFFSet.write(params, cwd + '/files/out-clp.zfp')
ZfpFFSet.write(params_il, cwd + '/files/out-TEAM_IL.zfp')
