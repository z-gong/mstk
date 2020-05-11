import os
import random
import subprocess
import sys
from subprocess import PIPE

from mstools.errors import PackmolError


class Packmol:
    '''
    wrappers for Packmol
    '''
    pass

    def __init__(self, packmol_bin):
        self.PACKMOL_BIN = packmol_bin
        self.numbers: [int]
        self.size: [float]

    def build_box(self, files: [str], numbers: [int], output: str,
                  size: [float] = None, length: float = None,
                  slab=None, slab_multiple=False,
                  tolerance: float = 1.8, seed: int = None,
                  inp_file='build.inp', silent=False) -> [int]:
        '''
        Build box directly from files
        '''

        if len(files) == 0:
            raise PackmolError('No files provided')
        if len(files) != len(numbers):
            raise PackmolError('Invalid numbers')
        self.numbers = numbers

        extensions = {filename.split('.')[-1].lower() for filename in files}
        if len(extensions) > 1:
            raise PackmolError('All file types should be the same')
        filetype = extensions.pop()

        if size is not None:
            if len(size) != 3:
                raise PackmolError('Invalid box size')
            else:
                self.size = size
        elif length is not None:
            self.size = [length, length, length]
        else:
            raise PackmolError('Box size needed')

        seed = seed or random.randint(1e7, 1e8)

        inp = (
            'filetype {filetype}\n'
            'tolerance {tolerance}\n'
            'output {output}\n'
            'seed {seed}\n'.format(filetype=filetype, tolerance=tolerance, output=output, seed=seed)
        )

        # liquid-gas interface
        if slab is not None:
            box_liq = '0 0 0 %f %f %f' % (self.size[0], self.size[1], slab)
            box_gas = '0 0 %f %f %f %f' % (slab, self.size[0], self.size[1], self.size[2])
            for i, filename in enumerate(files):
                # put 1/50 molecules in gas phase. Do not put too many in case of nucleation in gas phase
                n_gas = numbers[i] // 50
                n_liq = numbers[i] - n_gas
                inp += (
                    'structure {filename}\n'
                    'number {n_liq}\n'
                    'inside box {box_liq}\n'
                    'end structure\n'
                    'structure {filename}\n'
                    'number {n_gas}\n'
                    'inside box {box_gas}\n'
                    'end structure\n'.format(filename=filename, n_liq=n_liq, n_gas=n_gas,
                                             box_liq=box_liq, box_gas=box_gas)
                )

        else:
            for i, filename in enumerate(files):
                number = numbers[i]
                # slab model for multiple components
                if slab_multiple:
                    lz_per_slab = self.size[3] / len(numbers)
                    box = '0 0 %f %f %f %f' % (i * lz_per_slab, self.size[0], self.size[1], (i + 1) * lz_per_slab)
                else:
                    box = '0 0 0 %f %f %f' % tuple(self.size)

                inp += (
                    'structure {filename}\n'
                    'number {number}\n'
                    'inside box {box}\n'
                    'end structure\n'.format(filename=filename, number=number, box=box)
                )

        with open(inp_file, 'w') as f:
            f.write(inp)

        # TODO subprocess PIPE not work for Packmol new version, do not know why
        if silent:
            if os.name == 'nt':
                os.system(self.PACKMOL_BIN + ' < %s > nul' % inp_file)
            else:
                os.system(self.PACKMOL_BIN + ' < %s > /dev/null' % inp_file)
        else:
            os.system(self.PACKMOL_BIN + ' < %s' % inp_file)

            # (stdout, stderr) = (PIPE, PIPE) if silent else (None, None)
            # sp = subprocess.Popen([self.PACKMOL_BIN], stdin=PIPE, stdout=stdout, stderr=stderr)
            # sp.communicate(input=inp.encode())
