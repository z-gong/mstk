import random
import subprocess


class Packmol:
    '''
    wrappers for Packmol
    '''
    pass

    @staticmethod
    def build_box(files: [str], numbers: [int], output: str,
                  size: [float] = None, length: float = None,
                  tolerance: float = None,
                  seed: int = None,
                  save_input: bool = False):
        '''
        Build box directly from files

        :param files:
        :param numbers:
        :param output:
        :param size:
        :param length:
        :param tolerance:
        :param seed:
        :param save_input:
        :return:
        '''
        if len(files) == 0:
            raise Exception('no files provided')
        if len(files) != len(numbers):
            raise Exception('invalid numbers')

        extensions = {filename.strip().split('.')[-1] for filename in files}
        if len(extensions) > 1:
            raise Exception('all file types should be the same')
        filetype = extensions.pop()

        if size != None:
            if len(size) != 3:
                raise Exception('Invalid box size')
            else:
                box = size
        elif length != None:
            box = [length, length, length]
        else:
            raise Exception('box size needed')

        tolerance = tolerance or 2.0

        seed = seed or random.randint(1e7, 1e8)

        inp = '''filetype {filetype}
tolerance {tolerance}
output {output}
seed {seed}
'''.format(filetype=filetype, tolerance=tolerance, output=output, seed=seed)

        for i, filename in enumerate(files):
            number = numbers[i]
            inp += '''
structure {filename}
  number {number}
  inside box 0 0 0 {box_size}
end structure
'''.format(filename=filename, number=number, box_size=' '.join(map(str, box)))

        if save_input:
            with open('build.inp', 'w') as f:
                f.write(inp)

        sp = subprocess.Popen(['packmol'], stdin=subprocess.PIPE)
        sp.communicate(input=inp.encode())

