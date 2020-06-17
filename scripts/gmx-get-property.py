#!/usr/bin/env python3
# coding=utf-8

import sys
import subprocess
from subprocess import Popen, PIPE

cmd = 'gmx -quiet -nobackup energy -f %s -b %s' % (sys.argv[1], sys.argv[3])

prop_str = sys.argv[2]
sp = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=PIPE)
sp_out = sp.communicate(input=prop_str.encode())[0]

for line in sp_out.decode().splitlines():
    if line.lower().startswith(prop_str.lower()):
        words = line.split()
        prop = words[0]
        ave = float(words[1])
        stderr = float(words[2])
        if prop_str.lower().startswith('#sur'):
            ave /= 20
            stderr /= 20
        if prop_str.lower().startswith('dens'):
            ave /= 1000
            stderr /= 1000
        print(prop, ave, stderr)
