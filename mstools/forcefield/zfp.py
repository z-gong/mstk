import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom
from .forcefield import ForceField
from .ffterm import *


class Zfp():
    '''
    Generate ForceField from ZFP file.

    ZFP is the default format for storing force field in mstools.
    '''
    @staticmethod
    def read(*files):
        '''
        Parse a ZFP file

        Parameters
        ----------
        files : list of str

        Returns
        -------
        ff : ForceField
        '''
        ff = ForceField()
        for file in files:
            Zfp._parse(ff, file)

        return ff

    @staticmethod
    def _parse(ff, file):
        try:
            tree = ET.ElementTree(file=file)
        except:
            raise Exception('Invalid ZFP file')

        root = tree.getroot()
        if root is None:
            raise Exception('Empty ZFP file')

        node = root.find('Setting')
        ff.vdw_cutoff = float(node.attrib['vdw_cutoff'])
        ff.vdw_long_range = node.attrib['vdw_long_range']
        ff.lj_mixing_rule = node.attrib['lj_mixing_rule']
        ff.scale_14_vdw = float(node.attrib['scale_14_vdw'])
        ff.scale_14_coulomb = float(node.attrib['scale_14_coulomb'])

        tags = {
            'AtomTypes'           : ff.atom_types,
            'ChargeIncrementTerms': ff.charge_increment_terms,
            'VdwTerms'            : ff.vdw_terms,
            'PairwiseVdwTerms'    : ff.pairwise_vdw_terms,
            'BondTerms'           : ff.bond_terms,
            'AngleTerms'          : ff.angle_terms,
            'DihedralTerms'       : ff.dihedral_terms,
            'ImproperTerms'       : ff.improper_terms,
            'PolarizableTerms'    : ff.polarizable_terms,
        }
        for tag, d in tags.items():
            node = root.find(tag)
            if node is None:
                continue
            for element in node:
                try:
                    cls = getattr(sys.modules['mstools.forcefield.ffterm'], element.tag)
                except:
                    raise Exception('Invalid tag %s' % element.tag)
                term = cls.from_zfp(element.attrib)
                if term.name in d.keys():
                    raise Exception('Duplicated term: %s' % str(term))
                d[term.name] = term

    @staticmethod
    def save_to(ff, file):
        '''
        Save ForceField to a ZFP file

        Parameters
        ----------
        ff : ForceField
        file : str
        '''
        root = ET.Element('ForceFieldTerms')

        attrib = {
            'vdw_cutoff'      : str(ff.vdw_cutoff),
            'vdw_long_range'  : ff.vdw_long_range,
            'lj_mixing_rule'  : ff.lj_mixing_rule,
            'scale_14_vdw'    : str(ff.scale_14_vdw),
            'scale_14_coulomb': str(ff.scale_14_coulomb),
        }
        node = ET.SubElement(root, 'Setting', attrib=attrib)

        tags = {
            'AtomTypes'           : ff.atom_types,
            'ChargeIncrementTerms': ff.charge_increment_terms,
            'VdwTerms'            : ff.vdw_terms,
            'PairwiseVdwTerms'    : ff.pairwise_vdw_terms,
            'BondTerms'           : ff.bond_terms,
            'AngleTerms'          : ff.angle_terms,
            'DihedralTerms'       : ff.dihedral_terms,
            'ImproperTerms'       : ff.improper_terms,
            'PolarizableTerms'    : ff.polarizable_terms,
        }
        for tag, d in tags.items():
            node = ET.SubElement(root, tag)
            for term in d.values():
                ET.SubElement(node, term.__class__.__name__, attrib=term.to_zfp())

        str_xml = minidom.parseString(ET.tostring(root)).toprettyxml(indent='  ')
        with open(file, 'w') as f:
            f.write(str_xml)
