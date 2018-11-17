def polyfit(x: [float], y: [float], degree: int, weight: [float] = None):
    """
    least squares polynomial fit
    when degree ==n, y = c0 + c1 * x + c2 * x**2 + ... + cn * x**n

    :param x:
    :param y:
    :param degree:
    :param weight:
    :return: coeff: array, np.array([c0, c1, ... , cn])
             score: int
    """
    from sklearn import linear_model
    from sklearn.preprocessing import PolynomialFeatures

    skx = list(zip(x))
    skv = list(y)

    poly = PolynomialFeatures(degree)
    skx_ = poly.fit_transform(skx)
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(skx_, skv, sample_weight=weight)
    return clf.coef_, clf.score(skx_, skv)


def polyval(x, coeff):
    from numpy.polynomial.polynomial import polyval as np_polyval
    return np_polyval(x, coeff)


def polyval_derivative(x: float, coeff: [float]) -> (float, float):
    '''
    when degree == 2, len(coeff) = 3,
        y = c0 + c1 * x + c2 * xx
        dy/dx = c1 + 2*c2 * x
    when degree == 3, len(coeff) = 4,
        y = c0 + c1 * x + c2 * xx + c3 * xxx
        dy/dx = c1 + 2*c2 * x + 3*c3 * xx

    :param x:
    :param coeff: [c0, c1, c2, ...]
    :return: y: float
             dy/dx: float
    '''
    from numpy.polynomial.polynomial import polyval as np_polyval
    y = np_polyval(x, coeff)

    degree = len(coeff) - 1
    dydx = 0
    for i in range(degree):
        dydx += (i + 1) * coeff[i + 1] * x ** i

    return y, dydx


def polyfit_2d(x: [float], y: [float], z: [float], degree: int, weight: [float] = None):
    """
    least squares polynomial fit
    when degree == 3, z = c0 + c1 * x + c2 * y + c3 * xx + c4 * xy + c5 * yy + c6 * xxx + c7 * xxy + c8 * xyy + c9 * yyy
    when degree == 4, z = ...

    :param x:
    :param y:
    :param z:
    :param degree:
    :param weight:
    :return: coeff: array, np.array([c0, c1, c2, ...])
             score: int
    """
    from sklearn import linear_model
    from sklearn.preprocessing import PolynomialFeatures

    skx = list(zip(x, y))
    skv = list(z)

    poly = PolynomialFeatures(degree)
    skx_ = poly.fit_transform(skx)
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(skx_, skv, sample_weight=weight)
    return clf.coef_, clf.score(skx_, skv)


def polyval_derivative_2d(x: float, y: float, degree: int, coeff: [float]) -> (float, float, float):
    '''
    when degree == 3, z = c0 + c1 * x + c2 * y + c3 * xx + c4 * xy + c5 * yy + c6 * xxx + c7 * xxy + c8 * xyy + c9 * yyy
    when degree == 4, z = ...

    :param x:
    :param y:
    :param degree:
    :param coeff: [c0, c1, c2, ...]
    :return: z: float
             dz/dx: float
             dz/dy: float
    '''
    if degree == 3:
        k0, k1, k2, k3, k4, k5, k6, k7, k8, k9 = coeff
        z = k0 \
            + k1 * x + k2 * y \
            + k3 * x * x + k4 * x * y + k5 * y * y \
            + k6 * x * x * x + k7 * x * x * y + k8 * x * y * y + k9 * y * y * y
        dzdx = k1 + 2 * k3 * x + k4 * y + 3 * k6 * x * x + 2 * k7 * y * x + k8 * y * y
        dzdy = k2 + k4 * x + 2 * k5 * y + k7 * x * x + 2 * k8 * x * y + 3 * k9 * y * y

    elif degree == 4:
        k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14 = coeff
        z = k0 \
            + k1 * x + k2 * y \
            + k3 * x * x + k4 * x * y + k5 * y * y \
            + k6 * x * x * x + k7 * x * x * y + k8 * x * y * y + k9 * y * y * y \
            + k10 * x ** 4 + k11 * x ** 3 * y + k12 * x ** 2 * y ** 2 + k13 * x * y ** 3 + k14 * y ** 4
        dzdx = k1 + 2 * k3 * x + k4 * y + 3 * k6 * x * x + 2 * k7 * y * x + k8 * y * y \
               + 4 * k10 * x ** 3 + 3 * k11 * y * x ** 2 + 2 * k12 * y ** 2 * x + k13 * y ** 3
        dzdy = k2 + k4 * x + 2 * k5 * y + k7 * x * x + 2 * k8 * x * y + 3 * k9 * y * y \
               + k11 * x ** 3 + 2 * k12 * x ** 2 * y + 3 * k13 * x * y ** 2 + 4 * k14 * y ** 3

    else:
        raise Exception('Available degree: 3, 4')

    return z, dzdx, dzdy


def curve_fit_rsq(func, x_list, y_list, guess=None, bounds=None) -> ((float), float):
    import numpy as np
    from scipy.optimize import curve_fit

    x_array = np.array(x_list)
    y_array = np.array(y_list)

    popt, pcov = curve_fit(func, x_array, y_array, guess, bounds=bounds)
    ss_tot = ((y_array - y_array.mean()) ** 2).sum()
    predict = np.array([func(x, *popt) for x in x_array])
    ss_res = ((y_array - predict) ** 2).sum()
    rsq = 1 - ss_res / ss_tot

    return tuple(popt), rsq


def logistic(x, A1, A2, x0, p):
    return (A1 - A2) / (1 + (x / x0) ** p) + A2


def logistic_derivative(x, A1, A2, x0, p):
    y = logistic(x, A1, A2, x0, p)
    dydx = -(A1 - A2) / (1 + (x / x0) ** p) ** 2 * p * (x / x0) ** p / x
    return y, dydx


def fit_logistic(x_list: [float], y_list: [float], guess: [float] = None, bounds=None) -> ((float), float):
    import numpy as np

    guess = guess or [y_list[0], 2 * y_list[-1] - y_list[0], x_list[-1], 1.0]
    bounds = bounds or ([-np.inf, -np.inf, -np.inf, 0], np.inf)

    return curve_fit_rsq(logistic, x_list, y_list, guess, bounds)


def fit_vle_tanh(x_list: [float], d_list: [float], guess: [float] = None, bounds=None) -> ((float), float):
    import numpy as np

    def func(x, c, A, r, s):
        return c + A * np.tanh((x - r) / s)

    guess = guess or [0, 1, 0, 1]
    bounds = bounds or (-np.inf, np.inf)

    return curve_fit_rsq(func, x_list, d_list, guess, bounds)


def vle_dminus(T, Tc, B):
    return B * (1 - T / Tc) ** 0.325


def fit_vle_dminus(T_list: [float], dminus_list: [float], guess=None, bounds=None) -> ((float), float):
    '''
    Fit critical temperature using VLE density
    dliq - dgas = B(1-T/Tc)**0.325

    :param x: [float], temperatures
    :param y: [float], dliq - dgas
    :return: ((Tc, B), score)
    '''
    import numpy as np

    guess = guess or [max(T_list) / 0.8, 1.0]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(vle_dminus, T_list, dminus_list, guess, bounds)


def vle_dplus(T, Dc, A, Tc):
    return 2 * (Dc + A * (1 - T / Tc))


def fit_vle_dplus(T_list: [float], dplus_list: [float], Tc, guess=None, bounds=None) -> ((float), float):
    '''
    Fit critical density using VLE density and critical temperature
    dliq + dgas = 2(Dc+A(1-T/Tc))

    :param x:  [float], temperatures
    :param y: [float], dliq + dgas
    :param Tc: float, critical temperature
    :return: ((Dc, A), score)
    '''
    import numpy as np

    guess = guess or [0.3, 1.0]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(lambda T, Dc, A: vle_dplus(T, Dc, A, Tc), T_list, dplus_list, guess, bounds)


def vle_st(T, A, n, Tc):
    return A * (1 - T / Tc) ** n


def fit_vle_st(T_list, st_list, Tc, guess=None, bounds=None):
    import numpy as np

    guess = guess or [50, 1.22]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(lambda T, A, n: vle_st(T, A, n, Tc), T_list, st_list, guess, bounds)


def vle_log10pvap(T, A, B):
    return A - B / T


def vle_pvap(T, A, B):
    return 10 ** vle_log10pvap(T, A, B)


def fit_vle_pvap(T_list, pvap_list, guess=None, bounds=None):
    import numpy as np

    y_array = np.log10(np.array(pvap_list))

    guess = guess or [10, 3000]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(vle_log10pvap, T_list, y_array, guess, bounds)
