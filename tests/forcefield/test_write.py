#!/usr/bin/env python3

from mstools.forcefield import PaduaFFSet, PpfFFSet

import os
cwd = os.path.dirname(os.path.abspath(__file__))
params = PaduaFFSet(cwd + '/files/clp.ff', cwd+'/files/clp-alpha.ff')

PpfFFSet.write(params, cwd + '/files/out-clp.ppf')
