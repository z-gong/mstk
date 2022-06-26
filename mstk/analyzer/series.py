import math

import numpy as np
from pandas import Series


def block_average(series, n_block=5):
    '''
    Calculate the block average and standard error

    Parameters
    ----------
    series : Series
        Time series
    n_block : int
        Number of blocks to use

    Returns
    -------
    ave : float
        Block average
    stderr : float
        Standard error

    '''
    block_aves = average_of_blocks(series, n_block)
    ave, stderr = np.mean(block_aves), np.std(block_aves, ddof=1) / math.sqrt(n_block)
    stderr = float('%.1e' % stderr)  # 2 effective number for stderr
    return ave, stderr


def average_of_blocks(series, n_block=5):
    '''
    Split data to several blocks and return the average of each block

    Parameters
    ----------
    series : Series
        Time series
    n_block : int
        Number of blocks

    Returns
    -------
    aves : list of float
        The average of each block
    '''
    n_points = len(series)
    block_size = n_points // n_block
    blocks = []
    for n in range(n_block - 1):
        blocks.append(series.iloc[block_size * n:block_size * (n + 1)])
    blocks.append(series.iloc[block_size * (n_block - 1):])
    block_aves = [np.mean(b) for b in blocks]
    return block_aves


def is_converged(series: Series, frac_min=0.5):
    '''
    Determine whether a time series has converged or not

    Parameters
    ----------
    series : Series
        Time series
    frac_min : float
        Consider this time series is converged only if the fraction of converged parts relative to the full series
        is larger than this threshold

    Returns
    -------
    converged : bool
        Converged or not
    when : float
        From when this times series converged
    '''
    from pymbar import timeseries

    n_points = len(series)
    array = np.array(series)
    t0, g, Neff_max = timeseries.detectEquilibration(array, nskip=max(1, n_points // 100))
    if t0 > n_points * (1 - frac_min):
        return False, series.index[t0]
    return True, series.index[t0]


def efficiency_with_block_size(data):
    '''
    Calculate the statistical efficiency of block average with different block size.
    It can be used to determine the maximum number of block size.

    Parameters
    ----------
    data : list of float

    Returns
    -------
    block_size : list of int
        List of block size
    efficiency : list of float
        List of statistic efficiency with different block size
    '''
    array = np.array(data)
    n_points = len(data)
    bsize_list = []
    n_block_list = []
    s_list = []

    for bsize in range(1, int(math.sqrt(n_points))):
        n_block = int(n_points / bsize)
        bsize_list.append(bsize)
        n_block_list.append(n_block)

    for n_block in range(int(math.sqrt(n_points)), 4, -1):
        bsize = int(n_points / n_block)
        bsize_list.append(bsize)
        n_block_list.append(n_block)

    for i, bsize in enumerate(bsize_list):
        n_block = n_block_list[i]
        blocks = np.array_split(array, n_block)
        ave_blocks = [np.mean(block) for block in blocks]
        std_ave_blocks = np.std(ave_blocks, ddof=1)
        s = bsize * std_ave_blocks ** 2 / np.std(array) ** 2
        s_list.append(s)

    return bsize_list, s_list


def mean_and_uncertainty(series: Series, inefficiency=None) -> (float, float):
    from pymbar import timeseries

    ave = np.mean(series)
    array = np.array(series)
    if inefficiency == None:
        inefficiency = timeseries.statisticalInefficiency(array)
    return ave, np.std(array, ddof=1) / math.sqrt(len(array) / inefficiency)

