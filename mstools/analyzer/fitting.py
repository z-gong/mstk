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


def fit_tanh(x_list: [float], d_list: [float], guess: [float] = None) -> (float, float, float):
    import numpy as np
    from scipy.optimize import curve_fit

    def func(x, c, A, r, s):
        return c + A * np.tanh((x - r) / s)

    x_array = np.array(x_list)
    y_array = np.array(d_list)

    guess = guess or [0, 1, 0, 1]

    popt, pcov = curve_fit(func, x_array, y_array, guess)
    return tuple(popt)


def vle_d_minus(T, Tc, B):
    return B * (1 - T / Tc) ** 0.325


def fit_vle_d_minus(T_list: [float], d_minus_list: [float], guess=None):
    '''
    Fit critical temperature using VLE density
    d_liq - d_gas = B(1-T/Tc)**0.325

    :param x: [float], temperatures
    :param y: [float], dliq - dgas
    :return: float, critical temeprature
    '''
    import numpy as np
    from scipy.optimize import curve_fit

    T_array = np.array(T_list)
    y_array = np.array(d_minus_list)

    guess = guess or [500.0, 1.0]

    popt, pcov = curve_fit(vle_d_minus, T_array, y_array, guess)

    return popt


def vle_d_plus(T, Dc, A, Tc):
    return 2 * (Dc + A * (1 - T / Tc))


def fit_vle_d_plus(T_list: [float], d_plus_list: [float], Tc, guess=None):
    '''
    Fit critical density using VLE density and critical temperature
    d_liq + d_gas = 2(d_critical+A(1-T/Tc))

    :param x:  [float], temperatures
    :param y: [float], dliq + dgas
    :param Tc: float, critical temperature
    :return:  float, critical density
    '''
    import numpy as np
    from scipy.optimize import curve_fit

    T_array = np.array(T_list)
    y_array = np.array(d_plus_list)

    guess = guess or [0.3, 1.0]

    popt, pcov = curve_fit(lambda T, Dc, A: vle_d_plus(T, Dc, A, Tc), T_array, y_array, guess)

    return popt


def vle_st(T, A, n, Tc):
    return A * (1 - T / Tc) ** n


def fit_vle_st(T_list, st_list, Tc, guess=None):
    import numpy as np
    from scipy.optimize import curve_fit

    T_array = np.array(T_list)
    y_array = np.array(st_list)

    guess = guess or [50, 1.22]

    popt, pcov = curve_fit(lambda T, A, n: vle_st(T, A, n, Tc), T_array, y_array, guess)

    return popt



def vle_log10pvap(T, A, B):
    return A - B / T

def vle_pvap(T, A, B):
    return 10 ** vle_log10pvap(T, A, B)

def fit_vle_pvap(T_list, pvap_list, guess=None):
    import numpy as np
    from scipy.optimize import curve_fit

    T_array = np.array(T_list)
    y_array = np.log10(np.array(pvap_list))

    guess = guess or [10, 3000]

    popt, pcov = curve_fit(vle_log10pvap, T_array, y_array, guess)

    return popt
