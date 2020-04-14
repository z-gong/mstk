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

    @staticmethod
    def write(params: FFSet, file):
        root = ET.Element('ForceFieldTerms')
        node_atype = ET.SubElement(root, 'AtomTypes')
        node_increment = ET.SubElement(root, 'ChargeIncrementTerms')
        node_vdw = ET.SubElement(root, 'VdwTerms')
        node_pair = ET.SubElement(root, 'PairwiseVdwTerms')
        node_bond = ET.SubElement(root, 'BondTerms')
        node_angle = ET.SubElement(root, 'AngleTerms')
        node_dihedral = ET.SubElement(root, 'DihedralTerms')
        node_improper = ET.SubElement(root, 'ImproperTerms')

        for term in params.atom_types.values():
            ET.SubElement(node_atype, term.identifier, attrib=term.to_zfp())
        for term in params.charge_increment_terms.values():
            ET.SubElement(node_increment, term.identifier, attrib=term.to_zfp())
        for term in params.vdw_terms.values():
            ET.SubElement(node_vdw, term.identifier, attrib=term.to_zfp())
        for term in params.pairwise_vdw_terms.values():
            ET.SubElement(node_pair, term.identifier, attrib=term.to_zfp())
        for term in params.bond_terms.values():
            ET.SubElement(node_bond, term.identifier, attrib=term.to_zfp())
        for term in params.angle_terms.values():
            ET.SubElement(node_angle, term.identifier, attrib=term.to_zfp())
        for term in params.dihedral_terms.values():
            ET.SubElement(node_dihedral, term.identifier, attrib=term.to_zfp())
        for term in params.improper_terms.values():
            ET.SubElement(node_improper, term.identifier, attrib=term.to_zfp())

        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(indent='  ')
        with open(file, 'w') as f:
            f.write(xmlstr)
