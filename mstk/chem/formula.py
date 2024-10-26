class Formula:
    '''
    Parse elements, numbers and charge from a non-standardized chemical formula.

    e.g. H4C3(COH2)2 will be parsed as five carbon, eight hydrogen and two oxygen.

    Parameters
    ----------
    mol_str : str

    Attributes
    ----------
    atoms : dict, [str, int]
        Dict of the symbol and number of each chemical element.
    charge : int
        Net charge of the molecule.
    '''

    def __init__(self, mol_str):
        chars = Formula._extract_chars(mol_str)
        chars = Formula._expand_chars(chars)
        atoms, charge = Formula._count_atoms(chars)
        self.atoms = Formula._hill_order(atoms)
        self.charge = charge

    @staticmethod
    def _extract_chars(formula):
        chars = []
        i = 0
        while i < len(formula):
            c = formula[i]
            if c.isupper() or c == '+' or c == '-':
                chars.append(c)
                i += 1
                while i < len(formula):
                    if formula[i].islower():
                        chars[-1] += formula[i]
                        i += 1
                    else:
                        break
            elif c.isdigit():
                chars.append(c)
                i += 1
                while i < len(formula):
                    if formula[i].isdigit():
                        chars[-1] += formula[i]
                        i += 1
                    else:
                        break
            elif c == '(' or c == ')':
                chars.append(c)
                i += 1
            else:
                raise Exception('Invalid character: %s' % c)

        if chars[0].isdigit():
            raise Exception('Invalid formula')
        if chars.count('(') != chars.count(')'):
            raise Exception('Unmatched brackets')
        return chars

    @staticmethod
    def _expand_chars(chars):
        chars = chars[:]
        temp_expanded_chars = []
        temp_chars = []
        while '(' in chars:
            _start = False
            i = 0
            while i < len(chars):
                c = chars[i]
                if c == '(':
                    _start = True
                    temp_expanded_chars = chars[:i]
                    temp_chars = []
                    i += 1
                elif c == ')':
                    if not _start:
                        raise Exception('Unmatched brackets')
                    if i + 1 < len(chars) and chars[i + 1].isdigit():
                        n = int(chars[i + 1])
                        temp_expanded_chars += temp_chars * n
                        if i + 2 < len(chars):
                            temp_expanded_chars += chars[i + 2:]
                    else:
                        n = 1
                        temp_expanded_chars += temp_chars * n
                        if i + 1 < len(chars):
                            temp_expanded_chars += chars[i + 1:]
                    chars = temp_expanded_chars
                    break
                else:
                    if _start:
                        temp_chars.append(c)
                    i += 1

        return chars

    @staticmethod
    def _count_atoms(chars):
        counts = {}
        charge = 0
        for i in range(len(chars)):
            c = chars[i]
            if c.isdigit():
                counts[chars[i - 1]] += int(c) - 1
            else:
                if c not in counts:
                    counts[c] = 1
                else:
                    counts[c] += 1

        for k in list(counts):
            if k.startswith('+') or k.startswith('-'):
                if len(k) > 1:
                    raise Exception('Invalid formula')
                if k == '+':
                    charge += counts.pop(k)
                elif k == '-':
                    charge -= counts.pop(k)

        return counts, charge

    @staticmethod
    def _hill_order(atoms):
        count_C = atoms.pop('C') if 'C' in atoms else 0
        count_H = atoms.pop('H') if 'H' in atoms else 0
        symbols = list(atoms.keys())
        symbols.sort()
        atoms_sorted = {}
        if count_C != 0:
            atoms_sorted['C'] = count_C
        if count_H != 0:
            atoms_sorted['H'] = count_H
        for s in symbols:
            atoms_sorted[s] = atoms[s]
        return atoms_sorted

    def to_str(self):
        '''
        Return the standardized formula in hill order.

        Returns
        -------
        formula_str : str
        '''
        string = ''
        for symbol, num in self.atoms.items():
            string += symbol
            if num != 1:
                string += str(num)
        if self.charge > 0:
            string += '+'
        if self.charge < 0:
            string += '-'
        if abs(self.charge) > 1:
            string += str(abs(self.charge))

        return string

    @property
    def n_heavy(self):
        '''
        Number of heavy atoms (non-hydrogen atoms) in this formula.

        Returns
        -------
        n_heavy : int
        '''
        n = 0
        for k, v in self.atoms.items():
            if k != 'H':
                n += v
        return n

    @property
    def n_h(self):
        '''
        Number of hydrogen atoms in this formula.

        Returns
        -------
        n_h : int
        '''
        n = 0
        for k, v in self.atoms.items():
            if k == 'H':
                n += v
        return n
