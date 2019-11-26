#!/usr/bin/env python3

import os
import sys
import shutil

sys.path.append('..')

from mstools.wrapper.gmx import GMX


def test_modify_lj96():
    # os.remove('new.top')
    # os.remove('new.itp')
    # os.remove('new.mdp')
    # os.remove('new.xvg')
    shutil.copy('150-C6.itp', 'new.itp')
    shutil.copy('150-C6.top', 'new.top')
    shutil.copy('150-C6.mdp', 'new.mdp')
    GMX.modify_lj96(['new.itp'],['new.top'], ['new.mdp'], ['new.xvg'])


if __name__ == '__main__':
    test_modify_lj96()
