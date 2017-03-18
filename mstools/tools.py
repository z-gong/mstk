import os

def count_atoms(filename):
    '''
    count atom numbers in PDB file
    '''
    if not os.path.exists(filename):
        raise Exception('file not exist')
    filetype = filename.split('.')[-1].lower()
    if filetype == 'pdb':
        with open(filename) as f:
            return f.read().count('ATOM')
    else:
        return 1

def gcd(numbers):
    '''
    calculte the greatest common divisor
    '''
    minimal = min(numbers)
    for i in range(minimal, 1, -1):
        flag = True
        for number in numbers:
            if number % i != 0:
                flag = False
        if flag:
            return i
    return 1

