from mstools.analyzer import canny


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
        peaks = canny.canny1d(interval, data, nms_t, low_t, debug=debug)
        ave = sum([abs(p[1]) for p in peaks]) / len(peaks)

        if debug:
            print('ave peak:', ave)

        if len(peaks) == 0:
            result = 'b'
        elif len(peaks) == 1:
            if peaks[0][2] * 2 * hpw_t < len(data):
                result = 'a'
            else:
                result = 'c'
        elif len(peaks) == 2:
            if not Ptre.test_updown(peaks):
                result = 'Error: Two peaks with same direction.'
            if Ptre.test_hpw(peaks, hpw_t, debug):
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
    def test_hpw(peaks, threshold, debug=False):

        for i in range(0, len(peaks) - 1):
            hpw_ratio = (peaks[i + 1][0] - peaks[i][0]) / (peaks[i][2] + peaks[i + 1][2])
            if debug:
                print('hpw ratio:', i, hpw_ratio)
            if hpw_ratio < threshold:
                return False
        return True

    @staticmethod
    def test_updown(peaks):
        for i in range(len(peaks) - 1):
            if peaks[i][0] * peaks[i + 1][0] > 0:
                return False
        return True
