import math
from .grofile import GroFile

class GroReporter(object):
    """GroReporter outputs a series of frames from a Simulation to a GRO file.

    To use it, create a PDBReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=False, subset=None, reportVelocity=False, append=False):
        """Create a GroReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int, str
            The interval (in time steps) at which to write frames
            if set to "logfreq", then write trajectory at [10, 20, 30, ..., 90, 100, 200, ..., 900, 1000, 2000, ...]
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        subset : list(int)=None
            If not None, only the selected atoms will be written
        """
        if not (type(reportInterval) is int or reportInterval=='logfreq'):
            raise ValueError('reportInterval should be a integer or "logfreq"')
        self._reportInterval = reportInterval
        self._enforcePeriodicBox = enforcePeriodicBox
        if append:
            self._out = open(file, 'a')
        else:
            self._out = open(file, 'w')
        self._reportVelocity = reportVelocity

        if subset is None:
            self._subset = None
        else:
            self._subset = subset[:]

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
        if type(self._reportInterval) is int:
            steps = self._reportInterval - simulation.currentStep % self._reportInterval
        else:
            if simulation.currentStep < 10:
                _base = 10
            else:
                _base = 10 ** math.floor(math.log10(simulation.currentStep))
            steps = _base - simulation.currentStep % _base

        return (steps, True, self._reportVelocity, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """
        time = state.getTime()
        positions = state.getPositions(asNumpy=True)
        velocities = state.getVelocities(asNumpy=True) if self._reportVelocity else None
        vectors = state.getPeriodicBoxVectors()
        GroFile.writeFile(simulation.topology, time, positions, vectors, self._out, self._subset, velocities)

        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        self._out.close()
