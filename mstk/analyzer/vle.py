import numpy as np
from pandas import Series
from .canny import canny1d

__all__ = [
    'check_vle_density',
    'check_interface',
    'N_vaporize_condense'
]


class Ptre():
    _PATTERN_A = 'a'
    _PATTERN_B = 'b'
    _PATTERN_C = 'c'
    _PATTERN_D = 'd'

    @staticmethod
    def test_data(data, interval, nms_t=11, low_t=0.3, hpw_t=1.0, debug=False):
        """ Test data and entry of other tests.
            Arguments:
                nms_t, low_t is written in canny.py
        """
        peaks = canny1d(interval, data, nms_t, low_t, debug=debug)
        ave = sum([abs(p[1]) for p in peaks]) / len(peaks)

        if debug:
            print('ave peak:', ave)

        if len(peaks) == 0:
            result = 'b'
        elif len(peaks) == 1:
            # if peaks[0][2] * 2 * hpw_t < len(data):
            #     result = 'a'
            # else:
            #     result = 'c'
            result = 'c'
        elif len(peaks) == 2:
            if not Ptre.test_updown(peaks):
                result = 'Error: Two peaks with same direction.'
            if Ptre.test_hpw(peaks, len(data) * interval, hpw_t, debug):
                result = 'a'
            else:
                result = 'c'
        else:
            if ave > 0.017:  # this value change with gauss_sigma. be careful when modify either
                result = 'd'
            else:
                result = 'b'

        return result, peaks

    # half peak width test
    @staticmethod
    def test_hpw(peaks, length, threshold, debug=False):

        for i in range(0, len(peaks) - 1):
            hpw_ratio = (peaks[i + 1][0] - peaks[i][0]) / (peaks[i][2] + peaks[i + 1][2])
            if debug:
                print('hpw ratio:', i, hpw_ratio)
            if hpw_ratio < threshold:
                return False
        hpw_ratio = (peaks[0][0] + length - peaks[-1][0]) / (peaks[-1][2] + peaks[0][2])
        if debug:
            print('hpw ratio:', -1, hpw_ratio)
        if hpw_ratio < threshold:
            return False
        return True

    @staticmethod
    def test_updown(peaks):
        for i in range(len(peaks) - 1):
            if peaks[i][0] * peaks[i + 1][0] > 0:
                return False
        return True


def check_vle_density(series):
    '''
    Check whether or not a density profile represents a vapor-liquid interface with canny edge method

    Parameters
    ----------
    series : Series
        The density series along z axis. The index of the series is z coordinates

    Returns
    -------
    is_interface : bool
        Whether or not there is a vapor-liquid interface
    is_center_gas : bool
        Whether or not the center region is gas phase
    nodes : list of float
        The location of interfaces in z axis

    '''
    interval = series.index[1] - series.index[0]

    result, peaks = Ptre.test_data(series, interval)  # the peaks starts from 0, not original index
    peaks._sort_hill(key=lambda x: x[0])

    if result != Ptre._PATTERN_A:
        return False, None, None

    nodes = [p[0] + series.index[0] for p in peaks]
    return True, peaks[0][1] > 0, nodes


def N_vaporize_condense(phases):
    '''
    Check how many times one molecules have vaporized or condensed

    Parameters
    ----------
    phases : list of str
        The time evolution of the phase one molecule stayed in.
        The elements in this list can only be 'l', 'i' or 'g', which means liquid, interface and gas phase.

    Returns
    -------
    N_vapor: int
        Times of vaporization (transit from liquid to gas phase)
    N_condense: int
        Times of condensation (transit from gas to liquid phase)

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
