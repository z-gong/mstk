import os
import random
import subprocess
import sys
from subprocess import PIPE

from mstk.errors import PackmolError


class Packmol:
    '''
    Wrapper for Packmol for building initial coordinates for simulation

    Parameters
    ----------
    packmol_bin : str
        Path of the packmol binary executable
    '''

    def __init__(self, packmol_bin):
        self.PACKMOL_BIN = packmol_bin

    def build_box(self, files, numbers, output, size=None, length=None, slab=None, slab_multiple=False,
                  tolerance=0.2, seed=None, inp_file='build.inp', silent=False):
        '''
        Build a box directly from molecule files.

        Currently, rectangular box of following models are supported:

        * Homogeneous bulk liquid or gas model
        * Liquid-vapor interface model
        * Multi-component liquid-liquid interface model

        The input molecule files can be in XYZ or PDB format.
        The format for all input files and output files should be the same.

        Parameters
        ----------
        files : list of str
            List of input XYZ or PDB files
        numbers : list of int
            Numbers to be packed of each molecules
        output : str
            File name of the ouput XYZ or PDB file
        size : list of float, or None
            Size of the packed rectangular box, in unit of nanometer
            If not set, argument `length` should be provided.
        length : float or None
            Length of the packed cubic box, in unit of nanometer
            If not set, argument `size` should be provided.
        slab : float or None
            If a vapor-liquid interface model is to be packed, it gives the z-coordinates of the interface, in unit of nanometer
            It conflicts with argument `slab_multiple`.
        slab_multiple : bool or None
            Set it to True if a liquid-liquid interface model is to be packed.
            It conflicts with argument `slab`.
        tolerance : float
            The minimum distance between atoms belongs to different molecules, in unit of nanometer
        seed : int or None
            The seed for randomizing positions.
            If not provided, a randomly generated seed will be used.
        inp_file : str
            The input file for running packmol to be written.
        silent : bool
            If set to True, packmol will run silently without writing output on the screen.
            The error message will still be written on the screen.
        '''
        Packmol.gen_inp(files, numbers, output, size, length, slab, slab_multiple, tolerance, seed, inp_file)

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

    @staticmethod
    def gen_inp(files, numbers, output, size=None, length=None, slab=None, slab_multiple=False,
                tolerance=0.2, seed=None, inp_file='build.inp'):
        '''
        Generate input file for running packmol.

        It is the same as :func:`build_box` except that packmol is not called after input file generated.
        See :func:`build_box` for the details of parameters.
        '''
        if len(files) == 0:
            raise PackmolError('No files provided')
        if len(files) != len(numbers):
            raise PackmolError('Invalid numbers')

        extensions = {filename.split('.')[-1].lower() for filename in files}
        if len(extensions) > 1:
            raise PackmolError('All file types should be the same')
        filetype = extensions.pop()

        if size is not None:
            if len(size) != 3:
                raise PackmolError('Invalid box size')
            else:
                size = size[:]
        elif length is not None:
            size = [length, length, length]
        else:
            raise PackmolError('Box size needed')

        seed = seed or random.randint(1e7, 1e8)

        inp = (
            'filetype {filetype}\n'
            'tolerance {tolerance}\n'
            'output {output}\n'
            'seed {seed}\n'.format(filetype=filetype, tolerance=tolerance * 10, output=output, seed=seed)
        )

        # liquid-gas interface
        if slab is not None:
            box_liq = '0 0 0 %f %f %f' % (size[0] * 10, size[1] * 10, slab * 10)  # nm -> A
            box_gas = '0 0 %f %f %f %f' % (slab * 10, size[0] * 10, size[1] * 10, size[2] * 10)  # nm -> A
            for i, filename in enumerate(files):
                # put 1/50 molecules in gas phase. Do not put too many in case of nucleation in gas phase
                if numbers[i] == 0:
                    continue
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
                if numbers[i] == 0:
                    continue
                # slab model for multiple components
                if slab_multiple:
                    lz_per_slab = size[3] / len(numbers)
                    box = '0 0 %f %f %f %f' % (
                        i * lz_per_slab * 10, size[0] * 10, size[1] * 10, (i + 1) * lz_per_slab * 10)
                else:
                    box = '0 0 0 %f %f %f' % tuple(element * 10 for element in size)

                inp += (
                    'structure {filename}\n'
                    'number {number}\n'
                    'inside box {box}\n'
                    'end structure\n'.format(filename=filename, number=numbers[i], box=box)
                )

        with open(inp_file, 'w') as f:
            f.write(inp)
