import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom
from .forcefield import ForceField
from .ffterm import *


class Zfp(ForceField):
    def __init__(self, *files):
        super().__init__()

        for file in files:
            self._parse(file)

    def _parse(self, file):
        try:
            tree = ET.ElementTree(file=file)
        except:
            raise Exception('Invalid ZFP file')

        root = tree.getroot()
        if root is None:
            raise Exception('Empty ZFP file')

        node = root.find('Setting')
        self.vdw_cutoff = float(node.attrib['vdw_cutoff'])
        self.vdw_long_range = node.attrib['vdw_long_range']
        self.lj_mixing_rule = node.attrib['lj_mixing_rule']
        self.scale_14_vdw = float(node.attrib['scale_14_vdw'])
        self.scale_14_coulomb = float(node.attrib['scale_14_coulomb'])

        tags = {
            'AtomTypes'           : self.atom_types,
            'ChargeIncrementTerms': self.charge_increment_terms,
            'VdwTerms'            : self.vdw_terms,
            'PairwiseVdwTerms'    : self.pairwise_vdw_terms,
            'BondTerms'           : self.bond_terms,
            'AngleTerms'          : self.angle_terms,
            'DihedralTerms'       : self.dihedral_terms,
            'ImproperTerms'       : self.improper_terms,
            'PolarizableTerms'    : self.polarizable_terms,
        }
        for tag, d in tags.items():
            node = root.find(tag)
            if node is None:
                continue
            for element in node:
                try:
                    cls = getattr(sys.modules['mstools.forcefield.ffterm'], element.tag)
                except:
                    raise Exception('Invalid term under node %s: %s' % (node, element.tag))
                term = cls.from_zfp(element.attrib)
                if term.name in d.keys():
                    raise Exception('Duplicated term: %s' % str(term))
                d[term.name] = term

    @staticmethod
    def save_to(params: ForceField, file):
        root = ET.Element('ForceFieldTerms')

        attrib = {
            'vdw_cutoff'      : str(params.vdw_cutoff),
            'vdw_long_range'  : params.vdw_long_range,
            'lj_mixing_rule'  : params.lj_mixing_rule,
            'scale_14_vdw'    : str(params.scale_14_vdw),
            'scale_14_coulomb': str(params.scale_14_coulomb),
        }
        node = ET.SubElement(root, 'Setting', attrib=attrib)

        tags = {
            'AtomTypes'           : params.atom_types,
            'ChargeIncrementTerms': params.charge_increment_terms,
            'VdwTerms'            : params.vdw_terms,
            'PairwiseVdwTerms'    : params.pairwise_vdw_terms,
            'BondTerms'           : params.bond_terms,
            'AngleTerms'          : params.angle_terms,
            'DihedralTerms'       : params.dihedral_terms,
            'ImproperTerms'       : params.improper_terms,
            'PolarizableTerms'    : params.polarizable_terms,
        }
        for tag, d in tags.items():
            node = ET.SubElement(root, tag)
            for term in d.values():
                ET.SubElement(node, term.__class__.__name__, attrib=term.to_zfp())

        str_xml = minidom.parseString(ET.tostring(root)).toprettyxml(indent='  ')
        with open(file, 'w') as f:
            f.write(str_xml)
