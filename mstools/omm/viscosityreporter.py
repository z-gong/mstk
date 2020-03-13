from simtk.openmm import app
from simtk import unit


class ViscosityReporter(object):
    """DrudeTemperatureReporter report the real and Drude temperature of Drude simulation
    """

    def __init__(self, file, reportInterval):
        """Create a PDBReporter.

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
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation: app.Simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        if not self._hasInitialized:
            self._hasInitialized = True
            print('#"Step"\t"Acceleration (nm/ps^2)"\t"Velocity at z=0 (nm/ps)"\t"1/Viscosity (1/Pa.s)"', file=self._out)

        try:
            acceleration = simulation.integrator.getCosAcceleration().value_in_unit(unit.nanometer / unit.picosecond**2)
        except:
            print('This integrator does not support periodic perturbation', file=self._out)
        else:
            vMax, invVis = simulation.integrator.getViscosity()
            vMax = vMax.value_in_unit(unit.nanometer / unit.picosecond)
            invVis = invVis.value_in_unit((unit.pascal*unit.second)**-1)
            print(simulation.currentStep, acceleration, vMax, invVis, sep='\t', file=self._out)

        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()
