import random, math

from pandas import Series
import numpy as np


def is_converged(series: Series, tolerance=0.9, debug=False) -> (bool, float):
    n_block = 5
    start = series.index[0]
    end = series.index[-1]
    span = end - start
    block_size = span / n_block

    blocks = []
    dists = []
    for i in range(n_block):
        block = series.loc[start + i * block_size: end + i * block_size]
        blocks.append(block)
        mu = np.mean(block)
        sigma = np.std(block)
        dists.append([mu, sigma])

    minimal = dists[0][0] - dists[0][1] * 3.8  # 3.8 means 0.9999 confidential
    maximal = dists[0][0] + dists[0][1] * 3.8
    height = max([0.3989 / dist[1] for dist in dists])  # height of normal distribution
    for dist in dists[1:]:
        minimal = min(minimal, dist[0] - dist[1] * 3.8)
        maximal = max(maximal, dist[0] + dist[1] * 3.8)

    mc_points = int(1E5)
    random_x_list = np.ndarray(shape=(mc_points,), dtype=float)
    random_y_list = np.ndarray(shape=(mc_points,), dtype=float)
    gauss_y_list = np.ndarray(shape=(mc_points, n_block), dtype=float)
    for n in range(mc_points):
        x = random.random() * (maximal - minimal) + minimal
        y = random.random() * height
        random_x_list[n] = x
        random_y_list[n] = y
        for i in range(n_block):
            gauss_y_list[n][i] = gauss(dists[i][0], dists[i][1], x)

    if debug:
        import pylab
        from scipy.stats import norm
        pylab.xlim([minimal, maximal])
        pylab.ylim([0, height])
        pylab.scatter(random_x_list, random_y_list, s=1)
        for i in range(n_block):
            dist = dists[i]
            x_series = np.linspace(minimal, maximal, 100)
            pylab.plot(x_series, norm.pdf(x_series, dist[0], dist[1]), linewidth=1 + i * 0.5)
        pylab.show()

    for i in range(n_block - 2):
        overlap = 0.0
        for n in range(mc_points):
            if random_y_list[n] < min(gauss_y_list[n][i:]):
                overlap += 1
        overlap /= mc_points
        overlap *= ((maximal - minimal) * height)
        if debug:
            print(i, overlap)
        if overlap >= tolerance:
            converged_from = start + block_size * i
            return True, converged_from

    else:
        return False, 0


def gauss(mu, sigma, x):
    h = 0.3989 / sigma
    exponent = -(x - mu) ** 2 / 2 / sigma ** 2
    return h * np.exp(exponent)


def block_average(series, n_block=5) -> [float, float]:
    '''
    Get block average and standard error
    '''
    n_points = len(series)
    if n_block is None:
        n_block = 5

    block_size = n_points // n_block
    blocks = []
    for n in range(n_block - 1):
        blocks.append(series.iloc[block_size * n:block_size * (n + 1)])
    blocks.append(series.iloc[block_size * (n_block - 1):])

    block_aves = [np.mean(b) for b in blocks]
    return [np.mean(series), np.std(block_aves, ddof=1) / math.sqrt(n_block)]
