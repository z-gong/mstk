import numpy as np
from pandas import Series
from .ptre import Ptre


def is_crystal() -> bool:
    pass


def is_interface(series: Series, debug=False) -> bool:
    interval = series.index[1] - series.index[0]

    result, peaks = Ptre.test_data(series, interval, debug=debug)
    peaks = [(0, -peaks[0][1], None)] + peaks + [(len(series) * interval, -peaks[-1][1], None)]

    if debug:
        try:
            import pylab
        except:
            print('matplotlib not found, cannot debug')
        else:
            t = np.linspace(0, interval * len(series), len(series))
            pylab.plot(t, series)
            if result in [Ptre._PATTERN_A, Ptre._PATTERN_D]:
                pylab.vlines([p[0] for p in peaks], 0, 1, colors='red')
            pylab.show()

            for i in range(0, len(peaks) - 1):
                print(peaks[i][1] < 0, '%d-%d' % (peaks[i][0], peaks[i + 1][0]))

    if result == Ptre._PATTERN_A:
        return True
    else:
        return False
