import xml.etree.ElementTree as ET
from mstk.forcefield.forcefield import ForceField
from mstk.forcefield.ffterm import *
from mstk import logger


class Zfp:
    '''
    Load ForceField from ZFP file.

    XML language is used by ZFP file to serialize the setting and all the FFTerms in a ForceField object.

    In ZFP file, there is no leading 1/2 for harmonic energy terms.
    The length is in unit of nm, and angle is in unit of degree.
    The energies for harmonic bond, harmonic angle and vdW are in unit of kJ/mol/nm^2, kJ/mol/rad^2 and kJ/mol.

    A force field term describing one topology element should only appear once.
    If there are duplicated terms, an Exception will be raised.
    e.g. HarmonicAngleTerm('c_4', 'c_4', 'h_1') and SDKAngleTerm('c_4', 'c_4', 'h_1') are considered as duplicated,
    because both of them are describing the angle type ('c_4', 'c_4', 'h_1').
    LJ126Term('c_4', 'h_1') and MieTerm('c_4', 'h_1') are also duplicated,
    because both of them are describing the vdW interactions between atom types 'c_4' and 'h_1'.

    ZFP format is deprecated. New features are not supported by ZFP format
    1. Comments for force field and ffterm

    Parameters
    ----------
    file : str

    Attributes
    ----------
    forcefield : ForceField
    '''

    def __init__(self, file):
        logger.warning('ZFP format is deprecated. Consider converting it to a ZFF file with ffconv.py')
        self.forcefield = ForceField()
        self._parse(file)

    def _parse(self, file):
        try:
            tree = ET.ElementTree(file=file)
        except:
            raise Exception('Invalid ZFP file')

        root = tree.getroot()
        if root is None:
            raise Exception('Empty ZFP file')

        ff = self.forcefield

        node = root.find('Setting')
        ff.vdw_cutoff = float(node.attrib['vdw_cutoff'])
        ff.vdw_long_range = node.attrib['vdw_long_range']
        ff.lj_mixing_rule = node.attrib['lj_mixing_rule']
        ff.scale_14_vdw = float(node.attrib['scale_14_vdw'])
        ff.scale_14_coulomb = float(node.attrib['scale_14_coulomb'])

        tags = {
            'AtomTypes'           : ff.atom_types,
            'VirtualSiteTerms'    : ff.virtual_site_terms,
            'ChargeIncrementTerms': ff.bci_terms,
            'VdwTerms'            : ff.vdw_terms,
            'PairwiseVdwTerms'    : ff.pairwise_vdw_terms,
            'BondTerms'           : ff.bond_terms,
            'AngleTerms'          : ff.angle_terms,
            'DihedralTerms'       : ff.dihedral_terms,
            'ImproperTerms'       : ff.improper_terms,
            'PolarTerms'          : ff.polar_terms,
        }
        _duplicated = []
        for tag, d in tags.items():
            node = root.find(tag)
            if node is None:
                continue
            for element in node:
                try:
                    term = FFTermFactory.create_from_zfp(element.tag, element.attrib)
                except:
                    raise Exception('Invalid tag or attributes: %s, %s' %
                                    (element.tag, str(element.attrib)))
                if term.name in d.keys():
                    _duplicated.append(term)
                else:
                    d[term.name] = term

        if _duplicated:
            msg = f'{len(_duplicated)} duplicated terms from {file}'
            for term in _duplicated:
                msg += f'\n        {term}'
            logger.error(msg)
            raise Exception('Duplicated terms in ZFP file')

    @staticmethod
    def save_to(ff, file):
        '''
        Save ForceField to a ZFP file

        Parameters
        ----------
        ff : ForceField
        file : str
        '''
        raise NotImplementedError('Writing to ZFP file has been deprecated')


ForceField.register_format('.zfp', Zfp)
