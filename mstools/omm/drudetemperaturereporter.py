from __future__ import absolute_import

import simtk.openmm as mm
from simtk.openmm.app import Simulation, Topology
from simtk import unit
from simtk.unit import kelvin, kilojoule_per_mole as kJ_mol
import numpy as np

class DrudeTemperatureReporter(object):
    """
    DrudeTemperatureReporter reports the temperature of Drude simulation
    The temperatures for three degrees of freedom are reported
    -- molecular center of mass, internal atomic and Drude temperature
    It's better to set the reportInterval larger than 10000 to avoid significant performance penalty
    """

    def __init__(self, file, reportInterval):
        """Create a DrudeTemperatureReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        """
        self._reportInterval = reportInterval
        self._out = open(file, 'w')
        self._hasInitialized = False

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, True, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        system: mm.System = simulation.system
        if not self._hasInitialized:
            topology: Topology = simulation.topology
            self.n_atom = system.getNumParticles()
            self.n_residue = topology.getNumResidues()
            self.residue_atoms = np.zeros(self.n_atom, dtype=int)  # record which residue the atoms are in
            self.mass_residues = np.zeros(self.n_residue) # record the mass of residues
            self.residues = [[] for i in range(self.n_residue)] # record the index of atoms in each residues
            for i, residue in enumerate(list(topology.residues())):
                for atom in residue.atoms():
                    self.residues[i].append(atom.index)
                    self.residue_atoms[atom.index] = i
                    self.mass_residues[i] += system.getParticleMass(atom.index).value_in_unit(unit.dalton)

            self.dof_com = np.count_nonzero(self.mass_residues) * 3
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                self.dof_com -= 3

            self.dof_atom = self.dof_drude = 0
            for i in range(self.n_atom):
                if system.getParticleMass(i) > 0 * unit.dalton:
                    self.dof_atom += 3
            self.dof_atom -= (self.dof_com + system.getNumConstraints())

            drude_set = set()
            self.pair_set = set()
            force = [f for f in system.getForces() if type(f) == mm.DrudeForce][0]
            self.dof_atom -= 3 * force.getNumParticles()
            self.dof_drude += 3 * force.getNumParticles()
            for i in range(force.getNumParticles()):
                i_drude, i_core = force.getParticleParameters(i)[:2]
                drude_set.add(i_drude)
                self.pair_set.add((i_drude, i_core))
            self.drude_array = np.array(list(drude_set))
            self.atom_array = np.array(list(set(range(self.n_atom)) - drude_set))

            self._hasInitialized = True
            print('#"Step"\t"T_COM"\t"T_Atom"\t"T_Drude"\t"KE_COM"\t"KE_Atom"\t"KE_Drude"', file=self._out)

        velocities = state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond)
        masses = np.array([system.getParticleMass(i).value_in_unit(unit.dalton) for i in range(self.n_atom)])

        vel_res = np.zeros([self.n_residue, 3])
        for i, atoms in enumerate(self.residues):
            if self.mass_residues[i] == 0:
                continue
            mv = masses[atoms][:, np.newaxis] * velocities[atoms]
            vel_res[i] = np.sum(mv, axis=0) / self.mass_residues[i]
        mvv_com = self.mass_residues * np.sum(vel_res ** 2, axis=1)
        ke_com = mvv_com.sum() / 2 * (unit.nanometer / unit.picosecond) ** 2 * unit.dalton
        t_com = (2 * ke_com / (self.dof_com * unit.MOLAR_GAS_CONSTANT_R))

        velocities -= np.array([vel_res[self.residue_atoms[i]] for i in range(self.n_atom)])
        for i_drude, i_core in self.pair_set:
            v_drude = velocities[i_drude]
            v_core = velocities[i_core]
            m_drude = masses[i_drude]
            m_core = masses[i_core]
            m_com = m_drude + m_core
            m_rel = m_drude * m_core / m_com
            v_com = (m_drude * v_drude + m_core * v_core) / m_com
            v_rel = v_drude - v_core
            velocities[i_drude] = v_rel
            velocities[i_core] = v_com
            masses[i_drude] = m_rel
            masses[i_core] = m_com
        mvv = masses * np.sum(velocities ** 2, axis=1)
        ke = mvv[self.atom_array].sum() / 2 * (unit.nanometer / unit.picosecond) ** 2 * unit.dalton
        ke_drude = mvv[self.drude_array].sum() / 2 * (unit.nanometer / unit.picosecond) ** 2 * unit.dalton
        t = (2 * ke / (self.dof_atom * unit.MOLAR_GAS_CONSTANT_R))
        t_drude = (2 * ke_drude / (self.dof_drude * unit.MOLAR_GAS_CONSTANT_R))
        print(simulation.currentStep,
              t_com.value_in_unit(kelvin), t.value_in_unit(kelvin), t_drude.value_in_unit(kelvin),
              ke_com.value_in_unit(kJ_mol), ke.value_in_unit(kJ_mol), ke_drude.value_in_unit(kJ_mol),
              sep='\t', file=self._out)


        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()
