def polyfit_2d(x, y, z, degree):
    from sklearn import linear_model
    from sklearn.preprocessing import PolynomialFeatures

    skx = list(zip(x, y))
    skv = list(z)

    poly = PolynomialFeatures(degree)
    skx_ = poly.fit_transform(skx)
    clf = linear_model.LinearRegression(fit_intercept=False)
    clf.fit(skx_, skv)
    return clf.coef_, clf.score(skx_, skv)


def polyval_derivative_2d(x, y, degree, coeff) -> (float, float, float):
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
