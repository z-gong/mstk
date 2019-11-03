"""
pdbreporter.py: Outputs simulation trajectories in PDB format

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2012 Stanford University and the Authors.
Authors: Peter Eastman
Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

import simtk.openmm as mm
from simtk import unit

class DrudeTemperatureReporter(object):
    """PDBReporter outputs a series of frames from a Simulation to a PDB file.

    To use it, create a PDBReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=None):
        """Create a PDBReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        self._out = open(file, 'w')
        self._topology = None
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
        return (steps, False, True, False, False, self._enforcePeriodicBox)

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
            self.dof_atom = self.dof_drude = 0
            for i in range(system.getNumParticles()):
                if system.getParticleMass(i) > 0*unit.dalton:
                    self.dof_atom += 3
            self.dof_atom -= system.getNumConstraints()
            if any(type(system.getForce(i)) == mm.CMMotionRemover for i in range(system.getNumForces())):
                self.dof_atom -= 3

            self.drude_core_set = set()
            self.pair_set = set()
            for idx in range(system.getNumForces()):
                force: mm.DrudeForce = system.getForce(idx)
                if type(force) != mm.DrudeForce:
                    continue
                self.dof_atom -= 3 * force.getNumParticles()
                self.dof_drude += 3 * force.getNumParticles()
                for i in range(force.getNumParticles()):
                    i_drude, i_core = force.getParticleParameters(i)[:2]
                    self.drude_core_set.add(i_drude)
                    self.drude_core_set.add(i_core)
                    self.pair_set.add((i_drude, i_core))
            self._hasInitialized = True
            print('#Step\tT_Atom\tT_Drude\tKE_Atom\tKE_Drude', file=self._out)

        ke_com = 0 * unit.kilojoule_per_mole
        ke_rel = 0 * unit.kilojoule_per_mole
        velocities = state.getVelocities()
        for i_drude, i_core in self.pair_set:
            v_drude = velocities[i_drude]
            v_core = velocities[i_core]
            m_drude = system.getParticleMass(i_drude)
            m_core = system.getParticleMass(i_core)
            m_com = m_drude + m_core
            m_rel = m_drude * m_core / m_com
            v_com = (m_drude * v_drude + m_core * v_core) / m_com
            v_rel = v_drude - v_core
            ke_com += m_com * (v_com[0] * v_com[0] + v_com[1] * v_com[1] + v_com[2] * v_com[2]) / 2
            ke_rel += m_rel * (v_rel[0] * v_rel[0] + v_rel[1] * v_rel[1] + v_rel[2] * v_rel[2]) / 2
        for ia in range(system.getNumParticles()):
            if ia in self.drude_core_set:
                continue
            v = velocities[ia]
            m = system.getParticleMass(ia)
            ke_com += m * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / 2

        t_com = (2 * ke_com / (self.dof_atom * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
        t_rel = (2 * ke_rel / (self.dof_drude * unit.MOLAR_GAS_CONSTANT_R)).value_in_unit(unit.kelvin)
        print(simulation.currentStep, t_com, t_rel, ke_com.value_in_unit(unit.kilojoule_per_mole),
              ke_rel.value_in_unit(unit.kilojoule_per_mole), sep='\t', file=self._out)


        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()
