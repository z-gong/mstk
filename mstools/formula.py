class Formula:
    def __init__(self):
        self.atomlist = []  # list of (atom, number)
        self.atomdict = {}

    @staticmethod
    def read(mol_str):
        mol = Formula()
        token_list = mol.get_token(mol_str)
        mol.calculate(token_list)
        mol.count()
        mol.atomlist = mol.sort()
        return mol

    def get_token(self, mol_str):
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

    def calculate(self, token_list):
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
        return ''.join([name + Formula.to_num(num) for name, num in self.atomlist])

    def count(self):
        for name, num in self.atomlist:
            if name not in self.atomdict:
                self.atomdict[name] = num
            else:
                self.atomdict[name] += num

    def sort(self):
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
    def to_num(num):
        if num == 1:
            return ''
        else:
            return str(num)

    @property
    def n_heavy_atom(self) -> int:
        n = 0
        for element in self.atomdict.keys():
            if element != 'H':
                n += 1
        return n

    @property
    def n_h(self) -> int:
        n = 0
        for element in self.atomdict.keys():
            if element == 'H':
                n += 1
        return n

    def estimate_density(self) -> float:
        formula = self.to_str()
        density = {
            'H2': 0.07,
            'He': 0.15,
        }
        if formula in density.keys():
            return density.get(formula)

        # estimate density based on the ratio of hydrogen atoms and heavy atoms
        try:
            density = 1.8 / (self.n_h / self.n_heavy_atom)
        except:
            density = 0.8

        return min(density, 0.8)
