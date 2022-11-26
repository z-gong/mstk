import itertools
from mstk import logger
from mstk.forcefield.forcefield import ForceField
from mstk.forcefield.ffterm import *


class Zff():
    '''
    Parse ForceField from ZFF file.

    ZFF is the default plain text force field format in `mstk`.
    Comparing to ZFP, ZFF is more human-friendly but less program-friendly.

    Several files can be read at one time.
    Only one ForceField object will be created and it will contain force field terms from all of the files.
    If there are duplicated terms, only the first term will be kept.
    The force field settings will be read from the first file.

    Parameters
    ----------
    files : list of str

    Attributes
    ----------
    forcefield : ForceField
    '''

    def __init__(self, *files):
        self.forcefield = ForceField()
        self._setting_read = False
        for file in files:
            self._parse(file)

    def _parse(self, file):
        with open(file) as f:
            lines = f.read().splitlines()

        ff = self.forcefield
        _duplicated = []
        for line in lines:
            parts = line.split('#')
            str_term = parts[0].strip()
            if str_term == '':
                continue
            if str_term.startswith('Setting'):
                if self._setting_read:
                    continue

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
                term.comments.append(parts[1].strip())

            try:
                ff.add_term(term)
            except:
                _duplicated.append(term)

        if _duplicated != []:
            msg = 'Following duplicated terms are omitted:'
            for term in _duplicated:
                msg += ' ' + str(term)
            logger.warning(msg)

        self._setting_read = True

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
        for attr in ForceField.SETTING_ATTRS:
            line += '%-16s %-16s %s\n' % ('Setting', attr, getattr(ff, attr))
        line += '\n'

        _header_printed = set()
        for term in itertools.chain(ff.atom_types.values(),
                                    ff.virtual_site_terms.values(),
                                    ff.qinc_terms.values(),
                                    ff.vdw_terms.values(),
                                    ff.pairwise_vdw_terms.values(),
                                    ff.bond_terms.values(),
                                    ff.angle_terms.values(),
                                    ff.dihedral_terms.values(),
                                    ff.improper_terms.values(),
                                    ff.polarizable_terms.values()):
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


ForceField.registor_format('.zff', Zff)
