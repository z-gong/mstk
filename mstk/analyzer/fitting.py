import numpy as np


def polyfit(x, y, degree, weight=None):
    '''
    Least square n-th order polynomial fitting:
    y = c0 + c1 * x + c2 * x**2 + ... + cn * x**n

    Parameters
    ----------
    x : list of float
    y : list of float
    degree : int
    weight : list of float, optional

    Returns
    -------
    coeff : ndarray
    rsq : float
    '''
    from sklearn import linear_model
    from sklearn.preprocessing import PolynomialFeatures

    skx = list(zip(x))
    skv = list(y)

    skx = np.array(skx).astype(np.float64)
    skv = np.array(skv).astype(np.float64)

    poly = PolynomialFeatures(degree)
    skx_ = poly.fit_transform(skx)
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(skx_, skv, sample_weight=weight)
    return clf.coef_, clf.score(skx_, skv)


def polyval(x, coeff):
    '''
    Evaluate the n-th order polynomial result:
    y = c0 + c1 * x + c2 * x**2 + ... + cn * x**n

    Parameters
    ----------
    x : float
    coeff : list of float

    Returns
    -------
    y : float
    '''
    from numpy.polynomial.polynomial import polyval as np_polyval
    return np_polyval(x, coeff)


def polyval_derivative(x, coeff):
    '''
    Evaluate the n-th order polynomial result and derivative:
    y = c0 + c1 * x + c2 * x**2 + ... + cn * x**n;
    dy/dx = c1 + 2*c2 * x + ... + n*cn * x**(n-1)

    Parameters
    ----------
    x : float
    coeff : list of float

    Returns
    -------
    y : float
    dy/dx : float
    '''
    from numpy.polynomial.polynomial import polyval as np_polyval
    y = np_polyval(x, coeff)

    degree = len(coeff) - 1
    dydx = 0
    for i in range(degree):
        dydx += (i + 1) * coeff[i + 1] * x ** i

    return y, dydx


def polyfit_2d(x, y, z, degree, weight=None):
    '''
    Least square n-th order 2-D polynomial fitting.
    e.g. the 3rd polynomial fitting:
    z = c0 + c1 * x + c2 * y + c3 * xx + c4 * xy + c5 * yy + c6 * xxx + c7 * xxy + c8 * xyy + c9 * yyy

    Parameters
    ----------
    x : list of float
    y : list of float
    z : list of float
    degree : int
    weight : list of float, optional

    Returns
    -------
    coeff : ndarray
    rsq : float
    '''
    from sklearn import linear_model
    from sklearn.preprocessing import PolynomialFeatures

    skx = list(zip(x, y))
    skv = list(z)

    skx = np.array(skx).astype(np.float64)
    skv = np.array(skv).astype(np.float64)

    poly = PolynomialFeatures(degree)
    skx_ = poly.fit_transform(skx)
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(skx_, skv, sample_weight=weight)
    return clf.coef_, clf.score(skx_, skv)


def polyval_derivative_2d(x, y, degree, coeff):
    '''
    Evaluate the n-th order 2-D polynomial result and derivative. Only 3-rd and 4-th order are supported.
    e.g. the 3rd polynomial fitting:
    z = c0 + c1 * x + c2 * y + c3 * xx + c4 * xy + c5 * yy + c6 * xxx + c7 * xxy + c8 * xyy + c9 * yyy

    Parameters
    ----------
    x : float
    y : float
    degree : [3, 4]
    coeff : list of float

    Returns
    -------
    z : float
    dz/dx : float
    dz/dy : float
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


def curve_fit_rsq(func, x_list, y_list, guess=None, bounds=None):
    '''
    Least square curve fitting

    Parameters
    ----------
    func : function
    x_list : list of float
    y_list : list of float
    guess : list of float, optional
    bounds : list of float, optional

    Returns
    -------
    coeff: tuple of float
    rsq : float
    '''
    from scipy.optimize import curve_fit

    x_array = np.array(x_list)
    y_array = np.array(y_list)

    popt, pcov = curve_fit(func, x_array, y_array, guess, bounds=bounds)
    ss_tot = ((y_array - y_array.mean()) ** 2).sum()
    predict = np.array([func(x, *popt) for x in x_array])
    ss_res = ((y_array - predict) ** 2).sum()
    rsq = 1 - ss_res / ss_tot

    return tuple(popt), rsq


def fit_opls_dihedral(x_list, y_list, guess=None, bounds=None):
    def func(x, k1, k2, k3):
        return k1 * (1 + np.cos(x)) + k2 * (1 - np.cos(2 * x)) + k3 * (1 + np.cos(3 * x))

    guess = guess or [0.0, 0.0, 0.0]
    bounds = bounds or (-np.inf, np.inf)

    return curve_fit_rsq(func, x_list, y_list, guess, bounds)


def tg_hyperbola(T, T0, d0, a, b, c):
    H = (T - T0) / 2 + np.sqrt((T - T0) ** 2 / 4 + np.exp(c))
    return d0 - a * (T - T0) - b * H


def fit_tg_hyperbola(T_list, d_list, guess=None, bounds=None):
    guess = guess or [400.0, 1.0, 0.0, 0.0, 0.0]
    bounds = bounds or ([0.0, 0.0, 0.0, 0.0, -np.inf], [1000.0, 2.0, np.inf, np.inf, np.inf])

    return curve_fit_rsq(tg_hyperbola, T_list, d_list, guess, bounds)


def logistic(x, A1, A2, x0, p):
    return (A1 - A2) / (1 + (x / x0) ** p) + A2


def logistic_derivative(x, A1, A2, x0, p):
    y = logistic(x, A1, A2, x0, p)
    dydx = -(A1 - A2) / (1 + (x / x0) ** p) ** 2 * p * (x / x0) ** p / x
    return y, dydx


def fit_logistic(x_list, y_list, guess=None, bounds=None):
    guess = guess or [y_list[0], 2 * y_list[-1] - y_list[0], x_list[-1], 1.0]
    bounds = bounds or ([-np.inf, -np.inf, -np.inf, 0], np.inf)

    return curve_fit_rsq(logistic, x_list, y_list, guess, bounds)


def fit_vle_tanh(x_list, d_list, guess=None, bounds=None):
    def func(x, c, A, r, s):
        return c + A * np.tanh((x - r) / s)

    guess = guess or [0, 1, 0, 1]
    bounds = bounds or (-np.inf, np.inf)

    return curve_fit_rsq(func, x_list, d_list, guess, bounds)


def vle_dminus(T, Tc, B):
    return B * (1 - T / Tc) ** 0.325


def fit_vle_dminus(T_list, dminus_list, guess=None, bounds=None):
    '''
    Fit critical temperature using VLE density:
    dliq - dgas = B(1-T/Tc)**0.325

    Parameters
    ----------
    T_list : list of float
    dminus_list : list of float
    guess : list of float, optional
    bounds : list of float, optional

    Returns
    -------
    coeff : tuple of float
        Tc and B
    rsq : float
    '''
    guess = guess or [max(T_list) / 0.8, 1.0]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(vle_dminus, T_list, dminus_list, guess, bounds)


def vle_dplus(T, Dc, A, Tc):
    return 2 * (Dc + A * (1 - T / Tc))


def fit_vle_dplus(T_list, dplus_list, Tc, guess=None, bounds=None):
    '''
    Fit critical density using VLE density and critical temperature:
    dliq + dgas = 2(Dc+A(1-T/Tc))

    Parameters
    ----------
    T_list : list of float
    dplus_list : list of float
    Tc : float
    guess : list of float, optional
    bounds : list of float, optional

    Returns
    coeff : tuple of float
        Tc and B
    rsq : float
    -------
    '''
    guess = guess or [0.3, 1.0]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(lambda T, Dc, A: vle_dplus(T, Dc, A, Tc), T_list, dplus_list, guess, bounds)


def vle_st(T, A, n, Tc):
    return A * (1 - T / Tc) ** n


def fit_vle_st(T_list, st_list, Tc, guess=None, bounds=None):
    guess = guess or [50, 1.22]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(lambda T, A, n: vle_st(T, A, n, Tc), T_list, st_list, guess, bounds)


def vle_log10pvap(T, A, B):
    return A - B / T


def vle_pvap(T, A, B):
    return 10 ** vle_log10pvap(T, A, B)


def fit_vle_pvap(T_list, pvap_list, guess=None, bounds=None):
    y_array = np.log10(np.array(pvap_list))

    guess = guess or [10, 3000]
    bounds = bounds or (0, np.inf)

    return curve_fit_rsq(vle_log10pvap, T_list, y_array, guess, bounds)
