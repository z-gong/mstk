class Formula():
    '''
    Parse elements and numbers from a non-standardized chemical formula.

    e.g. H4C3(COH2)2 will be parsed as five carbon, eight hydrogen and two oxygen.

    Parameters
    ----------
    mol_str : str
    
    Attributes
    ----------
    atomlist : list of tuple
        Each element of this list is a tuple represents the symbol and number of each chemical element.
    atomdict : dict, [str, int]
        Dict of the symbol and number of each chemical element.
    '''
    def __init__(self, mol_str=None):
        self.atomlist = []  # list of (atom, number)
        self.atomdict = {}
        if mol_str is not None:
            self._load(mol_str)

    def _load(self, mol_str):
        token_list = self._get_token(mol_str)
        self._calculate(token_list)
        self._count()
        self.atomlist = self._sort_hill()

    @staticmethod
    def read(mol_str):
        '''
        An alias of the constructor for compatibility concern.

        Parameters
        ----------
        mol_str : str

        Returns
        -------
        formula : Formula
        '''
        mol = Formula()
        token_list = mol._get_token(mol_str)
        mol._calculate(token_list)
        mol._count()
        mol.atomlist = mol._sort_hill()
        return mol

    def _get_token(self, mol_str):
        tmp = ''
        tmp_num = ''
        token_list = []

        for i in mol_str:
            if i in ['(', ')']:
                if tmp:
                    token_list.append(tmp)
                    tmp = ''
                if tmp_num:
                    token_list.append(tmp_num)
                    tmp_num = ''
                elif token_list and token_list[-1] == ')':
                    token_list.append('1')

                token_list.append(i)
            elif i.isdigit():
                if tmp:
                    token_list.append(tmp)
                    tmp = ''
                tmp_num += i

            else:
                if tmp_num:
                    token_list.append(tmp_num)
                    tmp_num = ''
                elif token_list and token_list[-1] == ')':
                    token_list.append('1')
                if i.isupper():
                    if tmp:
                        token_list.append(tmp)
                    tmp = i
                else:
                    tmp += i
        if tmp:
            token_list.append(tmp)
        if tmp_num:
            token_list.append(tmp_num)
        elif token_list and token_list[-1] == ')':
            token_list.append('1')

        return token_list

    def _calculate(self, token_list):
        for token in token_list:
            if token == '(':
                self.atomlist.append(('(', 0))
            elif token == ')':
                tmp_list = []
                while self.atomlist[-1][0] != '(':
                    tmp_list.append(self.atomlist.pop())
                self.atomlist.pop()
                self.atomlist.append((tmp_list, 1))
            elif token.isdigit():

                if isinstance(self.atomlist[-1][0], list):
                    tmp_list = self.atomlist.pop()[0]
                    for atom, cnt0 in tmp_list:
                        self.atomlist.append((atom, cnt0 * int(token)))
                else:
                    atom, cnt0 = self.atomlist.pop()
                    self.atomlist.append((atom, cnt0 * int(token)))
            else:
                self.atomlist.append((token, 1))

    def to_str(self):
        '''
        Return the standardized formula in hill order.

        Returns
        -------
        formula_str : str
        '''
        return ''.join([name + Formula._to_num(num) for name, num in self.atomlist])

    def _count(self):
        for name, num in self.atomlist:
            if name not in self.atomdict:
                self.atomdict[name] = num
            else:
                self.atomdict[name] += num

    def _sort_hill(self):
        outlist = []
        C_cnt = None
        H_cnt = None
        if 'C' in self.atomdict:
            outlist.append(('C', self.atomdict['C']))
            C_cnt = self.atomdict['C']
            del self.atomdict['C']
            if 'H' in self.atomdict:
                outlist.append(('H', self.atomdict['H']))
                H_cnt = self.atomdict['H']
                del self.atomdict['H']

        outlist += list(sorted(self.atomdict.items(), key=lambda x: x[0]))
        if C_cnt:
            self.atomdict['C'] = C_cnt
        if H_cnt:
            self.atomdict['H'] = H_cnt
        return outlist

    @staticmethod
    def _to_num(num):
        if num == 1:
            return ''
        else:
            return str(num)

    @property
    def n_heavy(self):
        '''
        Number of heavy atoms (non-hydrogen atoms) in this formula.

        Returns
        -------
        n_heavy : int
        '''
        n = 0
        for k, v in self.atomdict.items():
            if k != 'H':
                n += v
        return n

    @property
    def n_heavy_atom(self) -> int:
        '''
        Alias of n_heavy for compatibility concern.

        Returns
        -------
        n_heavy : int
        '''
        return self.n_heavy

    @property
    def n_h(self):
        '''
        Number of hydrogen atoms in this formula.

        Returns
        -------
        n_h : int
        '''
        n = 0
        for k, v in self.atomdict.items():
            if k == 'H':
                n += v
        return n
