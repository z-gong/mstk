import sys
import math
from simtk.openmm.app.gromacsgrofile import GromacsGroFile
from simtk.unit import nanometer, picosecond, norm, is_quantity

class GroFile(GromacsGroFile):
    @staticmethod
    def writeFile(topology, time, positions, vectors, file):
        GroFile.writeHeader(time, file)
        GroFile.writeModel(topology, positions, file)
        GroFile.writeFooter(vectors, file)

    @staticmethod
    def writeHeader(time, file=sys.stdout):
        """Write out the header for a PDB file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the molecular system being written
        file : file=stdout
            A file to write the file to
        """
        print("written by openmm t = %.3f ps" % time.value_in_unit(picosecond), file=file)

    @staticmethod
    def writeModel(topology, positions, file=sys.stdout):
        """Write out a model to a PDB file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the model to write
        positions : list
            The list of atomic positions to write
        file : file=stdout
            A file to write the model to
        modelIndex : int=None
            If not None, the model will be surrounded by MODEL/ENDMDL records
            with this index
        keepIds : bool=False
            If True, keep the residue and chain IDs specified in the Topology
            rather than generating new ones.  Warning: It is up to the caller to
            make sure these are valid IDs that satisfy the requirements of the
            PDB format.  No guarantees are made about what will happen if they
            are not, and the output file could be invalid.
        extraParticleIdentifier : string=' '
            String to write in the element column of the ATOM records for atoms whose element is None (extra particles)
        """

        if len(list(topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometer)
        if any(math.isnan(norm(pos)) for pos in positions):
            raise ValueError('Particle position is NaN')
        if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError('Particle position is infinite')
        print('%i' % len(positions), file=file)
        atomIndex = 1
        for (chainIndex, chain) in enumerate(topology.chains()):
            residues = list(chain.residues())
            for (resIndex, res) in enumerate(residues):
                resName = res.name[:4]
                resId = res.id
                for atom in res.atoms():
                    atomName = atom.name[:5]
                    coords = positions[atomIndex-1]
                    line = '%5i%5s%5s%5i%8.3f%8.3f%8.3f' %(
                        int(resId), resName, atomName, atomIndex,
                        coords[0], coords[1], coords[2])
                    print(line, file=file)
                    atomIndex += 1

    @staticmethod
    def writeFooter(periodicBoxVectors, file=sys.stdout):
        """Write out the footer for a PDB file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the molecular system being written
        file : file=stdout
            A file to write the file to
        """
        vectors = periodicBoxVectors.value_in_unit(nanometer)
        xx, xy, xz = vectors[0]
        yx, yy, yz = vectors[1]
        zx, zy, zz = vectors[2]
        print("%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f" % (
            xx, yy, zz, xy, xz, yx, yz, zx, zy), file=file)
