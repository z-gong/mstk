import itertools
from mstk.forcefield.forcefield import ForceField
from mstk.forcefield.ffterm import *


class Zff:
    '''
    Parse ForceField from ZFF file.

    ZFF is the default plain text force field format in `mstk`.

    In ZFF file, there is no leading 1/2 for harmonic energy terms.
    The length is in unit of nm, and angle is in unit of degree.
    The energies for harmonic bond, harmonic angle and vdW are in unit of kJ/mol/nm^2, kJ/mol/rad^2 and kJ/mol.
    The two values for LJ interaction are epsilon and sigma.

    A force field term describing one topology element should only appear once.
    If there are duplicated terms, an Exception will be raised.

    There are three types of comments:
    1. Lines start with '#' are completely ignored by the parser.
    2. Lines start with '*' are comments for the force field. They will be stored in the ForceField object as a list of str.
       Each such line will be an element of comments.
    3. '#' in the middle of line is separator for parameters and comments. These comments will be stored in the FFTerm object as a list of str.
       The comments in one line can contain several elements, separated by ';'

    Parameters
    ----------
    file : str

    Attributes
    ----------
    forcefield : ForceField
    '''

    def __init__(self, file):
        self.forcefield = ForceField()
        self._parse(file)

    def _parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        ff = self.forcefield
        _duplicated = []
        for line in lines:
            if line.startswith('*'):
                ff.comments.append(line.lstrip('*').strip())
                continue
            parts = line.split('#')
            str_term = parts[0].strip()
            if str_term == '':
                continue
            if str_term.startswith('Setting'):
                words = str_term.split()
                key, str_val = words[1:]
                if key not in ForceField.SETTING_ATTRS:
                    raise Exception('Invalid setting: ' + key)
                setattr(ff, key, ForceField.SETTING_ATTRS[key](str_val))
                continue

            try:
                term = FFTermFactory.create_from_zff(str_term)
            except:
                raise Exception('Invalid force field line: ' + str_term)

            if len(parts) > 1:
                for comment in parts[1].split(';'):
                    c = comment.strip()
                    if c:
                        term.comments.append(c)

            try:
                ff.add_term(term)
            except:
                _duplicated.append(term)

        if _duplicated:
            msg = f'{len(_duplicated)} duplicated terms from {file}'
            for term in _duplicated:
                msg += f'\n        {term}'
            logger.error(msg)
            raise Exception('Duplicated terms in ZFF file')

    @staticmethod
    def save_to(ff, file):
        '''
        Save ForceField to a ZFF file

        Parameters
        ----------
        ff : ForceField
        file : str
        '''
        line = ''
        if ff.comments:
            for comment in ff.comments:
                line += f'* {comment}\n'
            line += '\n'

        for attr in ForceField.SETTING_ATTRS:
            line += '%-16s %-16s %s\n' % ('Setting', attr, getattr(ff, attr))
        line += '\n'

        _header_printed = set()
        for term in itertools.chain(ff.atom_types.values(),
                                    ff.virtual_site_terms.values(),
                                    ff.bci_terms.values(),
                                    ff.vdw_terms.values(),
                                    ff.pairwise_vdw_terms.values(),
                                    ff.bond_terms.values(),
                                    ff.angle_terms.values(),
                                    ff.dihedral_terms.values(),
                                    ff.improper_terms.values(),
                                    ff.polar_terms.values()):
            alias = type(term).get_alias()
            if alias not in _header_printed:
                line += term.to_zff_header() + '\n'
                _header_printed.add(alias)
            line += term.to_zff()
            if term.comments:
                line += '  # ' + '; '.join(term.comments)
            line += '\n'

        with open(file, 'wb') as f:
            f.write(line.encode())


ForceField.register_format('.zff', Zff)
