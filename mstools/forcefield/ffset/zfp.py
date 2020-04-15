import sys
import xml.etree.ElementTree as ET
from xml.dom import minidom
from .ffset import FFSet
from ..ffterm import *


class ZfpFFSet(FFSet):
    def __init__(self, file):
        super().__init__()
        self.lj_mixing_rule = self.LJ_MIXING_LB
        self.scale_14_vdw = 0.5
        self.scale_14_coulomb = 1.0 / 1.2

        self.parse(file)

    def parse(self, file):
        try:
            tree = ET.ElementTree(file=file)
        except:
            raise Exception('Invalid ZFP file')

        root = tree.getroot()
        if root is None:
            raise Exception('Empty ZFP file')

        tags = {
            'AtomTypes'           : self.atom_types,
            'ChargeIncrementTerms': self.charge_increment_terms,
            'VdwTerms'            : self.vdw_terms,
            'PairwiseVdwTerms'    : self.pairwise_vdw_terms,
            'BondTerms'           : self.bond_terms,
            'AngleTerms'          : self.angle_terms,
            'DihedralTerms'       : self.dihedral_terms,
            'ImproperTerms'       : self.improper_terms,
            'DrudeTerms'          : self.drude_terms,
        }
        for tag, d in tags.items():
            node = root.find(tag)
            if node is None:
                continue
            for element in node:
                cls = getattr(sys.modules['mstools.forcefield.ffterm'], element.tag)
                at = cls.from_zfp(element.attrib)
                d[at.name] = at

    @staticmethod
    def write(params: FFSet, file):
        root = ET.Element('ForceFieldTerms')
        tags = {
            'AtomTypes'           : params.atom_types,
            'ChargeIncrementTerms': params.charge_increment_terms,
            'VdwTerms'            : params.vdw_terms,
            'PairwiseVdwTerms'    : params.pairwise_vdw_terms,
            'BondTerms'           : params.bond_terms,
            'AngleTerms'          : params.angle_terms,
            'DihedralTerms'       : params.dihedral_terms,
            'ImproperTerms'       : params.improper_terms,
            'DrudeTerms'          : params.drude_terms,
        }
        for tag, d in tags.items():
            node = ET.SubElement(root, tag)
            for term in d.values():
                ET.SubElement(node, term.__class__.__name__, attrib=term.to_zfp())

        str_xml = minidom.parseString(ET.tostring(root)).toprettyxml(indent='  ')
        with open(file, 'w') as f:
            f.write(str_xml)
