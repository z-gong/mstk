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
