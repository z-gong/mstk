import math

from ..utils import create_mol_from_smiles, estimate_density_from_formula
from ..wrapper import Packmol, DFF


class Simulation():
    '''
    Base class of predefined simulation protocols.

    Simulation defines protocol for running specific type of simulation to get desired properties.
    e.g. NPT simulation with GROMACS to get cohesive energy, heat capacity and density,
    NPT_PPM simulation with GROMACS to get viscosity,
    QM calculation with Gaussian to get intramolecular heat capacity, etc...

    Parameters
    ----------
    packmol : Packmol
    dff : DFF
    jobmanager : subclass of JobManager
    packmol_bin : str, optional
        Binary of packmol. Deprecated.
    dff_root : str, optional
        Root directory of DFF. Deprecated.
    dff_db : str, optional
        Default database of DFF. Deprecated.
    dff_table : str, optional
        Default table of DFF. Deprecated.
    '''

    def __init__(self, packmol=None, dff=None, packmol_bin=None, dff_root=None, dff_db=None,
                 dff_table=None, jobmanager=None):
        if packmol is not None:
            self.packmol = packmol
        elif packmol_bin is not None:
            self.packmol = Packmol(packmol_bin=packmol_bin)
        if dff is not None:
            self.dff = dff
        elif dff_root is not None:
            self.dff = DFF(dff_root=dff_root, default_db=dff_db, default_table=dff_table)
        self.jobmanager = jobmanager
        self.procedure = None

        # minimum numbers of atoms and molecules for running this kind of simulation
        # if n_mol_list or n_atoms is specified in set_system(), n_atom_default and n_mol_default will be ignored
        # otherwise, both n_atom_default and n_mol_default should be satisfied
        self.n_atom_default: int = None
        self.n_mol_default: int = None

        self.n_mol_list: [int]
        self.msd = 'init.msd'
        self.pdb = 'init.pdb'

        # self._single_msd and self._single_pdb are used for exporting force field files
        self._single_msd = '_single.msd'
        self._single_pdb = '_single.pdb'

    def set_system(self, smiles_list, n_mol_list=None, n_atoms=None, n_mol_ratio=None,
                   length=None, density=None, name_list=None):
        '''
        Specify the components and size of the simulation system.

        The species are defined as SMILES string. The system can be mixtures.
        e.g. for pure hexane, smiles_list=['CCCCCC'];
        for the mixture of benzene and water, smiles_list=['c1ccccc1','O']
        This can used for calculating the phase equilibrium of mixtures or the solvation free energy.

        n_mol_list, n_atoms, n_mol_ratio are used to specify the number of molecules.
        n_mol_list defines the molecular number of different species.
        e.g. a system contains 200 hexane molecules, smiles_list=['CCCCCC'], n_mol_list=[200]
        a mixture system contains 1 benzene and 1000 water molecules, smiles_list=['c1ccccc1','O'], n_mol_list=[1, 1000]

        If n_mol_list is defined, n_atoms and n_mol_ratio will be ignored.
        If n_mol_list is None, then n_atoms and n_mol_ratio are used to determine the molecular number of species.
        n_atoms is the *approximate* number of atoms (including hydrogen) the system contains.
        n_mol_ratio is the ratio of number of molecules of species. If not specified, all species will have the same number of molecules
        e.g. for water system, smiles_list=['O'], n_atoms=100, the generated system will contains 34 water molecules
        for mixture system, the n_mol_ratio should be specified, which provides the ratio of number of species.
        e.g. for the mixture of benzene and water, smiles_list=['c1ccccc1','O'], n_atoms=100, n_mol_ratio=[1, 2],
        the generated system will contain 6 benzene and 12 water molecules. (the number of atoms will be 108)

        If n_mol_list is None and n_atoms is None, then the components will be determined by self.n_atom_default and self.n_mol_default
        self.n_atom_default gives the minimum number of atoms in the system
        self.n_mol_default gives the minimum number of molecules in the system
        both of these two criterion should be satisfied

        length and density determine the size of the simulation box.
        If length is defined, a cubic box with specified length will be generated, and density will be ignored.
        If length is None, then density is used to determine the size of the cubic box
        If length is None and density is None, then the density of the systems will be guessed based on the components.

        name_list specifies the residue name of species. It will be used as the molecular names in GROMACS top file.
        Be careful that name should contains at most three characters and starts with Alphabet letters.
        e.g. for the mixture of benzene and water, smiles_list=['c1ccccc1','O'], name_list=['BEN', 'WAT']
        If name_list is None, the names of species will be ['MO0', 'MO1', 'MO2'...]
        Normally, it is not recommended to define molecular names. The default names are more general.

        Parameters
        ----------
        smiles_list : list of str
            List of SMILES of the molecules to be simulated
        n_mol_list : list of int, optional
            List of numbers of molecules
        n_atoms : int, optional
            Approximate number of atoms
        n_mol_ratio : list of int, optional
            The ratio between numbers of molecules for mixture system
        length : float
            The length of cubic simulation box
        density : float
            The density of cubic simulation box
        name_list : list of str
            The name of molecules. They will be used in mol2, msd and GROMACS top file
        '''
        if type(smiles_list) != list:
            raise Exception('smiles_list should be list')
        self.smiles_list = smiles_list[:]
        self.pdb_list = []
        self.mol2_list = []
        n_components = len(smiles_list)
        n_atom_list = []  # number of atoms of each molecule
        molwt_list = []  # molecule weight of each molecule
        density_list = []  # estimated density of each molecule
        for i, smiles in enumerate(smiles_list):
            pdb = 'mol-%i.pdb' % i
            mol2 = 'mol-%i.mol2' % i
            if name_list is not None:
                resname = name_list[i]
            else:
                resname = 'MO%i' % i
            py_mol = create_mol_from_smiles(smiles, pdb_out=pdb, mol2_out=mol2, resname=resname)
            self.pdb_list.append(pdb)
            self.mol2_list.append(mol2)
            n_atom_list.append(len(py_mol.atoms))
            molwt_list.append(py_mol.molwt)
            density_list.append(estimate_density_from_formula(
                py_mol.formula) * 0.9)  # * 0.9, build box will be faster

        if n_mol_list is not None:
            self.n_mol_list = n_mol_list[:]
        else:
            if n_mol_ratio is None:
                n_mol_ratio = [1] * n_components
            n_atom_all = sum([n_atom_list[i] * n for i, n in enumerate(n_mol_ratio)])
            if n_atoms is None:
                n_atoms = n_atom_all * math.ceil(self.n_mol_default / sum(n_mol_ratio))
                n_atoms = max(n_atoms, self.n_atom_default)
            self.n_mol_list = [math.ceil(n_atoms / n_atom_all) * n for n in n_mol_ratio]

        mass = sum([molwt_list[i] * self.n_mol_list[i] for i in range(n_components)])

        if length is not None:
            self.length = length
            self.box = [length, length, length]
            self.vol = self.length ** 3
        else:
            if density is None:
                density = sum(
                    [density_list[i] * self.n_mol_list[i] for i in range(n_components)]) / sum(
                    self.n_mol_list)
            self.vol = 10 / 6.022 * mass / density
            self.length = self.vol ** (1 / 3)  # assume cubic box
            self.box = [self.length, self.length, self.length]

    def build(self):
        '''
        Build the initial configuration.

        This method should be implemented by subclasses.
        '''
        pass

    def prepare(self):
        '''
        Generate input files for simulation engines, and return the commands for running the simulation.

        This method should be implemented by subclasses.

        Returns
        -------
        commands : list of str
            The commands for running the simulation.
        '''
        pass

    def run(self):
        '''
        Submit the job to PBS scheduler
        '''
        self.jobmanager.submit()

    def check_finished(self):
        '''
        Check if the job is successfully finished.

        This method should be implemented by subclasses.
        '''
        pass

    def analyze(self):
        '''
        Analyze simulation data and obtain raw results like energy, density...

        This method should be implemented by subclasses.
        '''
        pass

    def clean(self):
        '''
        Delete intermediate simulation data.

        This method should be implemented by subclasses.
        '''
        pass

    def post_process(self, **kwargs):
        '''
        Post process the raw simulation results,
        e.g. fit the density and energy at different temperatures and pressures.

        This method should be implemented by subclasses.
        '''
        pass

    def get_post_data(self, **kwargs):
        '''
        Derive other results from the post-processed data,
        e.g. heat capacity, thermal expansion, critical point...

        This method should be implemented by subclasses.
        '''
        pass
