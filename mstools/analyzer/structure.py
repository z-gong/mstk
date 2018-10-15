import numpy as np
from pandas import Series
from .ptre import Ptre


def is_crystal() -> bool:
    pass


def check_vle_density(series: Series) -> (bool, bool, [float]):
    '''
    Check in a density series is interface

    :param series:  density series. The index is the z axis
    :return: is interface: bool
             center is gas phase: bool
             location of changing nodes [float]
    '''
    interval = series.index[1] - series.index[0]

    result, peaks = Ptre.test_data(series, interval)  # the peaks starts from 0, not original index
    peaks.sort(key=lambda x: x[0])

    if result != Ptre._PATTERN_A:
        return False, None, None

    nodes = [p[0] + series.index[0] for p in peaks]
    return True, peaks[0][1] > 0, nodes


def N_vaporize_condense(phases) -> (int, int):
    '''
    Check the number of vaporize and condensation

    :param list:  ['l', 'i', 'g'], liquid, interface, gas
    :return: (N_vapor, N_condense)
    '''
    N_vapor = N_condense = 0
    phases = [i for i in phases if i != 'i']
    while len(phases) > 1:
        if phases[0] == 'l' and phases[1] == 'g':
            N_vapor += 1
        elif phases[0] == 'g' and phases[1] == 'l':
            N_condense += 1
        phases.pop(0)

    return N_vapor, N_condense


def check_interface(series: Series, debug=False) -> (bool, list):
    interval = series.index[1] - series.index[0]

    result, peaks = Ptre.test_data(series, interval, debug=debug)

    if debug:
        try:
            import pylab
        except:
            print('matplotlib not found, cannot debug')
        else:
            t = np.array(series.index) - series.index[0]
            pylab.plot(t, series)
            if result in [Ptre._PATTERN_A, Ptre._PATTERN_D]:
                pylab.vlines([p[0] for p in peaks], 0, 1, colors='red')
            pylab.show()

            for i in range(0, len(peaks)):
                if peaks[i][1] < 0:
                    state = 'INCREASING'
                else:
                    state = 'DECREASING'
                print(state, peaks[i][0])

    if result == Ptre._PATTERN_A:
        return True, peaks
    else:
        return False, peaks


def angular_momentum(com, mass_list, corr_list, vel_list, vel_com):
    ### TODO Not sure whether it is correct
    Lx = Ly = Lz = 0
    for i, corr in enumerate(corr_list):
        dx = corr[0] - com[0]
        dy = corr[1] - com[1]
        dz = corr[2] - com[2]
        vx = vel_list[i][0] - vel_com[0]
        vy = vel_list[i][1] - vel_com[1]
        vz = vel_list[i][2] - vel_com[2]
        mass = mass_list[i]
        Lx += mass * (dy * vz - dz * vy)
        Ly += mass * (dz * vx - dx * vz)
        Lz += mass * (dx * vy - dy * vx)
    return Lx, Ly, Lz


def velocity_com(mass_list, vel_list):
    Vx = Vy = Vz = 0
    for i, vel in enumerate(vel_list):
        vx = vel_list[i][0]
        vy = vel_list[i][1]
        vz = vel_list[i][2]
        mass = mass_list[i]

        Vx += mass * vx
        Vy += mass * vy
        Vz += mass * vz
    total_mass = sum(mass_list)
    Vx /= total_mass
    Vy /= total_mass
    Vz /= total_mass

    return Vx, Vy, Vz
