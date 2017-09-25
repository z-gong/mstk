from collections import OrderedDict
from typing import Dict


class Parameter():
    def __init__(self, string):
        string = string.strip()
        if string.endswith('*'):
            self.fixed = True
            self.value = float(string[:-1])
        else:
            self.fixed = False
            self.value = float(string)

    def __str__(self):
        if self.fixed:
            return str(self.value) + '*'
        else:
            return str(self.value)


class FFTerm():
    def __init__(self, term, key, value, comment):
        self.term = term.strip()
        self.key = key.strip()
        self.value = value.strip()
        self.comment = comment.strip()

    def __str__(self):
        return '%s: %s: %s: %s' % (self.term, self.key, self.value, self.comment)


class BINC(FFTerm):
    def __init__(self, key, value, comment=''):
        super().__init__('BINC', key, value, comment)
        atom_types = key.strip().split(',')
        self.atom1 = atom_types[0].strip()
        self.atom2 = atom_types[1].strip()
        self.binc = Parameter(value)

    def __str__(self):
        return '%s: %s, %s: %s: %s' % (self.term, self.atom1, self.atom2, self.binc, self.comment)


class LJ(FFTerm):
    def __init__(self, key, value, comment=''):
        super().__init__('N12_6', key, value, comment)
        self.atom = key.strip()
        self.r0 = Parameter(value.split(',')[0])
        self.e0 = Parameter(value.split(',')[1])

    def __str__(self):
        return '%s: %s: %s, %s: %s' % (self.term, self.atom, self.r0, self.e0, self.comment)


class BHARM(FFTerm):
    def __init__(self, key, value, comment=''):
        super().__init__('BHARM', key, value, comment)
        atom_types = key.strip().split(',')
        self.atom1 = atom_types[0].strip()
        self.atom2 = atom_types[1].strip()
        self.b0 = Parameter(value.split(',')[0])
        self.kb = Parameter(value.split(',')[1])

    def __str__(self):
        return '%s: %s, %s: %s, %s: %s' % (self.term, self.atom1, self.atom2, self.b0, self.kb, self.comment)


class AHARM(FFTerm):
    def __init__(self, key, value, comment=''):
        super().__init__('AHARM', key, value, comment)
        atom_types = key.strip().split(',')
        self.atom1 = atom_types[0].strip()
        self.atom2 = atom_types[1].strip()
        self.atom3 = atom_types[2].strip()
        self.a0 = Parameter(value.split(',')[0])
        self.ka = Parameter(value.split(',')[1])

    def __str__(self):
        return '%s: %s, %s, %s: %s, %s: %s' % (
            self.term, self.atom1, self.atom2, self.atom3, self.a0, self.ka, self.comment)


class TCOSP(FFTerm):
    def __init__(self, key, value, comment=''):
        super().__init__('TCOSP', key, value, comment)
        atom_types = key.strip().split(',')
        self.atom1 = atom_types[0].strip()
        self.atom2 = atom_types[1].strip()
        self.atom3 = atom_types[2].strip()
        self.atom4 = atom_types[3].strip()

        self.k1 = None
        self.k2 = None
        self.k3 = None
        para_list = value.split(',')
        para_list = [para.strip() for para in para_list]
        n_multi = len(para_list) // 3
        for i in range(n_multi):
            k = Parameter(para_list[3 * i + 1])
            if para_list[3 * i + 2].startswith('1'):
                self.k1 = k
            elif para_list[3 * i + 2].startswith('2'):
                self.k2 = k
            elif para_list[3 * i + 2].startswith('3'):
                self.k3 = k

    def __str__(self):
        return '%s: %s, %s, %s, %s: %s: %s' % (self.term, self.atom1, self.atom2, self.atom3, self.atom4,
                                               self.para_str, self.comment)

    @property
    def para_str(self):
        para_list = []
        if self.k3 is not None:
            para_list.append('0*, %s, 3*' % self.k3)
        if self.k1 is not None:
            para_list.append('0*, %s, 1*' % self.k1)
        if self.k2 is not None:
            para_list.append('180*, %s, 2*' % self.k2)
        return ', '.join(para_list)

    def freeze(self):
        if self.k1 is not None:
            self.k1.fixed = True
        if self.k2 is not None:
            self.k2.fixed = True
        if self.k3 is not None:
            self.k3.fixed = True


def get_ppf_term_from_line(line):
    words = line.strip().split(':')
    if len(words) >= 4:
        comment = words[3]
    else:
        comment = ''

    if line.startswith('BINC'):
        term = BINC(words[1], words[2], comment)
    elif line.startswith('N12_6'):
        term = LJ(words[1], words[2], comment)
    elif line.startswith('BHARM'):
        term = BHARM(words[1], words[2], comment)
    elif line.startswith('AHARM'):
        term = AHARM(words[1], words[2], comment)
    elif line.startswith('TCOSP'):
        term = TCOSP(words[1], words[2], comment)
    else:
        term = FFTerm(words[0], words[1], words[2], comment)
    return term


class PPF():
    def __init__(self, ppf_file=None, string=None):
        lines = []
        if ppf_file is not None:
            with open(ppf_file) as f:
                lines = f.read().splitlines()
        if string is not None:
            lines = string.splitlines()

        self.pre_lines = []
        self.terms = []

        PPF_START = False
        for line in lines:
            line = line.strip()
            if line == '#DFF:PPF':
                PPF_START = True
            if line.startswith('#') or not PPF_START:
                self.pre_lines.append(line)
            else:
                self.terms.append(get_ppf_term_from_line(line))

    def __str__(self):
        string = '\n'.join(self.pre_lines)
        for term in self.terms:
            string += '\n' + str(term)
        return string

    def write(self, ppf_out):
        with open(ppf_out, 'w') as f:
            f.write(str(self))

    def get_adj_nb_paras(self) -> OrderedDict:
        adj_paras = OrderedDict()
        for term in self.terms:
            if term.term == 'N12_6':
                if not term.r0.fixed:
                    adj_paras[term.atom + '_r0'] = term.r0.value
                if not term.e0.fixed:
                    adj_paras[term.atom + '_e0'] = term.e0.value

            elif term.term == 'BINC':
                if not term.binc.fixed:
                    adj_paras['%s_%s_bi' % (term.atom1, term.atom2)] = term.binc.value

        return adj_paras

    def set_nb_paras(self, new_paras: Dict, delta=False):
        for term in self.terms:
            if term.term == 'N12_6':
                key = term.atom + '_e0'
                if key in new_paras.keys():
                    term.e0.value = new_paras[key]
                    term.e0.fixed = False
                key = term.atom + '_r0'
                if key in new_paras.keys():
                    term.r0.value = new_paras[key]
                    term.r0.fixed = False
            elif term.term == 'BINC':
                key = '%s_%s_bi' % (term.atom1, term.atom2)
                if key in new_paras.keys():
                    term.binc.value = new_paras[key]
                    term.binc.fixed = False

        ### temperature dependent
        if delta:
            for term in self.terms:
                if term.term == 'N12_6':
                    ### scale r0
                    key = term.atom[:3] + '_dr'  # c_4_dr, h_1_dr
                    if key in new_paras.keys():
                        term.r0.value *= (1 + new_paras[key])
                    else:
                        key = 'all_dr'  # all_dr
                        if key in new_paras.keys():
                            term.r0.value *= (1 + new_paras[key])

                    ### scale e0
                    key = term.atom[:3] + '_de'  # c_4_de, h_1_de
                    if key in new_paras.keys():
                        term.e0.value *= (1 + new_paras[key])
                    else:
                        key = 'all_de'  # all_de
                        if key in new_paras.keys():
                            term.e0.value *= (1 + new_paras[key])

                    ### scale C6
                    key = term.atom[:3] + '_dl'  # c_4_dl, h_1_dl
                    if key in new_paras.keys():
                        term.r0.value /= (1 + new_paras[key]) ** (1 / 6)
                        term.e0.value *= (1 + new_paras[key]) ** 2
                    else:
                        key = 'all_dl'
                        if key in new_paras.keys():
                            term.r0.value /= (1 + new_paras[key]) ** (1 / 6)
                            term.e0.value *= (1 + new_paras[key]) ** 2

    def freeze_torsions(self):
        for term in self.terms:
            if term.term == 'TCOSP':
                term.freeze()

    def relax_torsion(self, torsion_key):
        for term in self.terms:
            if term.term == 'TCOSP':
                key_words = torsion_key.split(',')
                key_words = [w.strip() for w in key_words]
                term_key_words = [term.atom1, term.atom2, term.atom3, term.atom4]
                if term_key_words == key_words or term_key_words == list(reversed(key_words)):
                    # term.k1 = Parameter('0')
                    # term.k2 = Parameter('0')
                    # term.k3 = Parameter('0')
                    if term.k1 is not None:
                        term.k1.fixed = False
                    if term.k2 is not None:
                        term.k2.fixed = False
                    if term.k3 is not None:
                        term.k3.fixed = False

    def modify_torsion(self, torsion_key, n, delta):
        for term in self.terms:
            if term.term == 'TCOSP':
                key_words = torsion_key.split(',')
                key_words = [w.strip() for w in key_words]
                term_key_words = [term.atom1, term.atom2, term.atom3, term.atom4]
                if term_key_words == key_words or term_key_words == list(reversed(key_words)):
                    if n == 1:
                        term.k1.value += delta
                    elif n == 2:
                        term.k2.value += delta
                    elif n == 3:
                        term.k3.value += delta

    def fit_torsion(self, dff_root, qmd, msd, restraint, torsion_key, dfi_name='fit_torsion'):
        import os
        import random
        from .dff import DFF

        ### relax only one torsion
        self.freeze_torsions()
        self.relax_torsion(torsion_key)

        ### Backup adj_nb_paras
        adj_nb_paras = self.get_adj_nb_paras()

        dff = DFF(dff_root=dff_root)
        ppf_tmp = 'tmp-%i.ppf' % random.randint(1E5, 1E6)
        ppf_out_tmp = 'tmp-%i.ppf' % random.randint(1E5, 1E6)
        self.write(ppf_tmp)
        try:
            # should set charge for MSD first
            dff.set_charge([msd], ppf_tmp)
            dff.fit_torsion(qmd, msd, ppf_tmp, ppf_out_tmp, restraint, dfi_name=dfi_name)
        except:
            raise
        else:
            self.__init__(ppf_out_tmp)
            ### restore adj_nb_paras
            self.set_nb_paras(adj_nb_paras)
            os.remove(ppf_tmp)
            os.remove(ppf_out_tmp)
            os.remove(dfi_name + '.dfi')
            os.remove(dfi_name + '.dfo')
            os.remove('set_charge.dfi')
            os.remove('set_charge.dfo')

def delta_ppf(ppf_file, ppf_out, T, paras_delta=None):
    if paras_delta is None:
        paras_delta = {
            'h_1_dr': -0.01 * (T - 298) / 100,
            'h_1_de': 0.056 * (T - 298) / 100,
            'c_4_dr': -0.01 * (T - 298) / 100,
            'c_4_de': 0.056 * (T - 298) / 100,
            'c_3_dr': -0.004 * (T - 298) / 100,
            'c_3_de': 0.005 * (T - 298) / 100,
        }
    ppf = PPF(ppf_file)
    ppf.set_nb_paras(paras_delta, delta=True)
    ppf.write(ppf_out)
