#!/usr/bin/env python3

from mstools.forcefield import Padua, Ppf, Zfp
from mstools.forcefield.ffterm import *

import os
cwd = os.path.dirname(os.path.abspath(__file__))

clp = Padua(cwd + '/files/clp.ff', cwd + '/files/clp-alpha.ff')
Zfp.save_to(clp, cwd + '/files/out-clp.zfp')

il = Ppf(cwd + '/files/TEAM_IL.ppf')
Zfp.save_to(il, cwd + '/files/out-TEAM_IL.zfp')
