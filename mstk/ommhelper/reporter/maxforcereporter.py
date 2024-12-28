import numpy as np


class MaxForceReporter:
    '''
    MaxForceReporter outputs the atom experiencing the largest force

    Parameters
    ----------
    file : string or file
        The file to write to, specified as a file name or file object
    reportInterval : int
        The interval (in time steps) at which to write frames
    append : bool
        If set to True, will append to file
    '''

    def __init__(self, file, reportInterval, append=False):
        self._reportInterval = reportInterval
        self._openedFile = isinstance(file, str)
        if self._openedFile:
            if append:
                self._out = open(file, 'a')
            else:
                self._out = open(file, 'w')
        else:
            self._out = file
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

        return (steps, False, False, True, False, False)

    def report(self, simulation, state):
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
            print('#"Step"\t"Atom id"\t"Force (kJ/mol/nm)"', file=self._out)

        forces = state.getForces(asNumpy=True)._value
        fsq = np.sum(forces * forces, axis=1)
        imax = np.argmax(fsq)
        fmax = fsq[imax] ** 0.5
        self._out.write(f'{simulation.currentStep}\t{imax}\t{fmax:.3f}\n')

        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()

    def __del__(self):
        if self._openedFile:
            self._out.close()
