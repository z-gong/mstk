#!/usr/bin/env python3

from mstools.forcefield import FFToolParameterSet, PpfParameterSet

import os
cwd = os.path.dirname(os.path.abspath(__file__))
params = FFToolParameterSet(cwd + '/files/clp.ff')

PpfParameterSet.write(params, cwd + '/files/out-clp.ppf')
