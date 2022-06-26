'''
A set of commonly used small functions.
'''

import os
import subprocess
import random
import string
import numpy as np


def greatest_common_divisor(numbers):
    '''
    Calculate the greatest common divisor.

    Parameters
    ----------
    numbers : list of int

    Returns
    divisor : int
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


def random_string(length=8):
    '''
    Generate a random string in specified length. The string contains only upper and lower case ascii letters.


    Parameters
    ----------
    length : int

    Returns
    -------
    string : str
    '''
    return ''.join(random.sample(string.ascii_letters, length))


def cd_or_create_and_cd(dir):
    '''
    Go to the target directory. If not exist, create this directory and go to it.

    Parameters
    ----------
    dir : str
    '''
    if not os.path.exists(dir):
        try:
            os.makedirs(dir)
        except:
            raise Exception('Cannot create directory: %s' % dir)

    try:
        os.chdir(dir)
    except:
        raise Exception('Cannot read directory: %s' % dir)


def estimate_density_from_formula(f) -> float:
    # in unit of g/mL
    from mstk.chem.formula import Formula
    formula = Formula(f)
    string = formula.to_str()
    density = {
        'H2': 0.07,
        'He': 0.15,
    }
    if string in density.keys():
        return density.get(string)

    nAtoms = formula.n_heavy + formula.n_h
    nC = formula.atoms.get('C', 0)
    nH = formula.atoms.get('H', 0)
    nO = formula.atoms.get('O', 0)
    nN = formula.atoms.get('N', 0)
    nS = formula.atoms.get('S', 0)
    nF = formula.atoms.get('F', 0)
    nCl = formula.atoms.get('Cl', 0)
    nBr = formula.atoms.get('Br', 0)
    nI = formula.atoms.get('I', 0)
    nOther = nAtoms - nC - nH - nO - nN - nS - nF - nCl - nBr - nI
    return (1.175 * nC + 0.572 * nH + 1.774 * nO + 1.133 * nN + 2.184 * nS
            + 1.416 * nF + 2.199 * nCl + 5.558 * nBr + 7.460 * nI
            + 0.911 * nOther) / nAtoms


def n_diff_lines(f1: str, f2: str):
    with open(f1) as f:
        l1 = f.readlines()
    with open(f2) as f:
        l2 = f.readlines()

    n1 = len(l1)
    n2 = len(l2)

    n_diff = 0
    for i in range(min(n1, n2)):
        if l1[i].strip() != l2[i].strip():
            n_diff += 1
    n_diff += abs(n1 - n2)
    return n_diff


def get_last_line(filename):
    # TODO implement windows version
    cmd = 'tail -n 1 %s' % filename
    try:
        out = subprocess.check_output(cmd.split()).decode()
    except:
        raise Exception('Cannot open file: %s' % filename)

    try:
        string = out.splitlines()[-1]
    except:
        string = ''
    return string


def histogram(data, bins, normed=False):
    y, _x = np.histogram(data, bins=bins, density=normed)
    x = (_x[1:] + _x[:-1]) / 2
    return x, y


def print_data_to_file(name_column_dict, file):
    if len(name_column_dict) == 0:
        raise Exception('Columns are empty')

    with open(file, 'w') as f:
        f.write('#')
        for name in name_column_dict.keys():
            f.write('"%s"\t' % name)
        f.write('\n')
        for i in range(len(list(name_column_dict.values())[0])):
            for column in name_column_dict.values():
                f.write('%f\t' % column[i])
            f.write('\n')


def align_mpl_axis(ax1, ax2):
    ylims1 = list(ax1.axes.get_ylim())
    ylims1[0] = min(0, ylims1[0])
    ylims1[1] = max(0, ylims1[1])
    ylims2 = list(ax2.axes.get_ylim())
    ylims2[0] = min(0, ylims2[0])
    ylims2[1] = max(0, ylims2[1])

    if ylims1[0] == 0:
        if ylims2[0] == 0:
            return
        elif ylims2[1] == 0:
            ylims1[0] = -ylims1[1]
            ylims2[1] = -ylims2[0]
        else:
            ylims1[0] = ylims1[1] / ylims2[1] * ylims2[0]
    elif ylims1[1] == 0:
        if ylims2[0] == 0:
            ylims1[1] = -ylims1[0]
            ylims2[0] = -ylims2[1]
        elif ylims2[1] == 0:
            return
        else:
            ylims1[1] = ylims1[0] / ylims2[0] * ylims2[1]
    else:
        if ylims2[0] == 0:
            ylims2[0] = ylims2[1] / ylims1[1] * ylims1[0]
        elif ylims2[1] == 0:
            ylims2[1] = ylims2[0] / ylims1[0] * ylims1[1]
        else:
            ratio1 = -ylims1[1] / ylims1[0]
            ratio2 = -ylims2[1] / ylims2[0]
            if ratio1 > ratio2:
                ylims2[1] = -ylims2[0] * ratio1
            else:
                ylims1[1] = -ylims1[0] * ratio2

    ax1.set(ylim=ylims1)
    ax2.set(ylim=ylims2)
