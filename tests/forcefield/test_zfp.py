#!/usr/bin/env python3

from mstools.forcefield import PaduaFFSet, PpfFFSet, ZfpFFSet
from mstools.forcefield.ffterm import *

import os
cwd = os.path.dirname(os.path.abspath(__file__))

clp = PaduaFFSet(cwd + '/files/clp.ff', cwd + '/files/clp-alpha.ff')
ZfpFFSet.save_to(clp, cwd + '/files/out-clp.zfp')

il = PpfFFSet(cwd + '/files/TEAM_IL.ppf')
ZfpFFSet.save_to(il, cwd + '/files/out-TEAM_IL.zfp')
