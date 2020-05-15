import os
import shutil


class CheckpointReporter(object):
    """CheckpointReporter saves periodic checkpoints of a simulation.
    The checkpoints will overwrite one another -- only the last checkpoint
    will be saved in the file.

    To use it, create a CheckpointReporter, then add it to the Simulation's
    list of reporters. To load a checkpoint file and continue a simulation,
    use the following recipe:

    >>> with open('checkput.chk', 'rb') as f:
    >>>     simulation.context.loadCheckpoint(f.read())

    Notes:
    A checkpoint contains not only publicly visible data such as the particle
    positions and velocities, but also internal data such as the states of
    random number generators.  Ideally, loading a checkpoint should restore the
    Context to an identical state to when it was written, such that continuing
    the simulation will produce an identical trajectory.  This is not strictly
    guaranteed to be true, however, and should not be relied on.  For most
    purposes, however, the internal state should be close enough to be
    reasonably considered equivalent.

    A checkpoint contains data that is highly specific to the Context from
    which it was created. It depends on the details of the System, the
    Platform being used, and the hardware and software of the computer it was
    created on.  If you try to load it on a computer with different hardware,
    or for a System that is different in any way, loading is likely to fail.
    Checkpoints created with different versions of OpenMM are also often
    incompatible.  If a checkpoint cannot be loaded, that is signaled by
    throwing an exception.

    """

    def __init__(self, file, reportInterval):
        """Create a CheckpointReporter.

        Parameters
        ----------
        file : string
            The file to write to. Any current contents will be overwritten.
        reportInterval : int
            The interval (in time steps) at which to write checkpoints.
        """

        self._reportInterval = reportInterval
        self._file = file

        if type(file) is not str:
            raise Exception('file should be str')

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        # Do a safe save.

        with open(self._file, 'w+b', 0) as out:
            out.write(simulation.context.createCheckpoint())

        filename = self._file + '_%i' % simulation.currentStep
        file_prev_prev = self._file + '_%i' % (simulation.currentStep - 2 * self._reportInterval)
        if os.path.exists(file_prev_prev):
            os.remove(file_prev_prev)
        shutil.copy(self._file, filename)
