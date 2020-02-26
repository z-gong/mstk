import sys
from simtk.openmm.app.gromacsgrofile import GromacsGroFile
from simtk.unit import nanometer, picosecond, norm, is_quantity

class GroFile(GromacsGroFile):
    @staticmethod
    def writeFile(topology, time, positions, vectors, file, subset):
        GroFile.writeHeader(time, file)
        GroFile.writeModel(topology, positions, file, subset)
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
    def writeModel(topology, positions, file=sys.stdout, subset=None):
        """Write out a model to a PDB file.

        Parameters
        ----------
        topology : Topology
            The Topology defining the model to write
        positions : list
            The list of atomic positions to write
        file : file=stdout
            A file to write the model to
        subset : list(int)=None
            If not None, only the selected atoms will be written
        """

        atoms = list(topology.atoms())
        if len(atoms) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometer)

        if subset is None:
            subset = list(range(len(atoms)))

        print('%i' % len(subset), file=file)
        for ii, i in enumerate(subset):
            atom = atoms[i]
            residue = atom.residue
            coords = positions[i]
            line = '%5i%5s%5s%5i%8.3f%8.3f%8.3f' % (
                int(residue.id), residue.name, atom.name, int(atom.id),
                coords[0], coords[1], coords[2])
            print(line, file=file)


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
        print(' %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f' %(
            xx, yy, zz, xy, xz, yx, yz, zx, zy), file=file)
