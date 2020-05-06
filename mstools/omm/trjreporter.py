from simtk import unit
from mstools.topology import Topology, Molecule, Atom
from mstools.trajectory import Frame, Trajectory


class TrjReporter(object):
    """TrjReporter outputs a series of frames from a Simulation to a Trajectory file supported by mstools.

    To use it, create a TrjReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, enforcePeriodicBox=False, subset=None, append=False):
        """Create a DCDReporter.

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
        subset : list(int)=None
            If not None, only the selected atoms will be written
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        if subset is None:
            self._subset = None
        else:
            self._subset = subset[:]
        if append:
            mode = 'a'
        else:
            mode = 'w'
        self._trj = Trajectory.open(file, mode)

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
        return (steps, True, False, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if not self._initialized:
            self._topology = Topology()
            self._topology.init_from_omm_topology(simulation.topology)
            self._frame = Frame(self._topology.n_atom)
            self._initialized = True

        self._frame.step = simulation.currentStep
        self._frame.time = state.getTime().value_in_unit(unit.picosecond)
        self._frame.positions = state.getPositions(asNumpy=True)
        self._frame.cell.set_box(state.getPeriodicBoxVectors(asNumpy=True))

        self._trj.write_frame(self._frame, self._topology, subset=self._subset)

    def __del__(self):
        self._trj.close()
