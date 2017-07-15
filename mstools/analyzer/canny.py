""" One dimension canny algorithm, for detecting edges
"""
# Last modified: 2017.2

import math
import numpy as np
def canny1d(interval, y, nms_width, peak_threshold=0.1, gauss_sigma=1, gauss_width=5, grad_width=3, debug=False):
    """ 1-dimension canny algorithm. 
     Arguments:
        interval: time interval;
        y: numpy array for raw data;
        gauss_sigma: Gaussian filter sigma, default is 0.35, determined by test, no modify recommended;
        guass_width: Gaussian filter width, default is 5;
        nms_width: Non-maximum supression width. It is the scan width of window in order to delete maximum.
            Modify this only when two much near peaks are found;
        peak_threshold: Peaks with absolute value lower than it will be deleted.
            Not recommend to set larger than 0.3;
        grad_width: Gradient operator width. Higher will make the curve smoother. Less than 5 is recommended.
    Returns:
        A list of tuple (t, value) for peaks;
    """
    if debug:
        try:
            import matplotlib.pyplot as plt
        except:
            print('matplotlib not found, cannot debug')
            debug = False

    y_oper = np.concatenate((y, y)) # stack twice for boundary results
    t = np.linspace(0,interval*len(y),len(y))
    t_oper = np.concatenate((t, t+len(y)*interval))

    if debug:
        plt.figure(1)
        plt.plot(t_oper, y_oper)
    smooth_y = uniform(np.convolve(y_oper, get_gaussian_1d(gauss_sigma, gauss_width), mode='same'))
    if debug:
        plt.plot(t_oper, smooth_y)
    grad_y = np.convolve(smooth_y, get_grad_operator(grad_width), mode='same')
    if debug:
        plt.figure(2)
        plt.plot(t_oper, grad_y)
    grad_y_nms = non_max_suppress(np.abs(grad_y), nms_width)
    if debug:
        plt.figure(2)
        plt.plot(t_oper, grad_y_nms)
    
    max_peak = np.max(grad_y_nms[len(y)//2:len(y)*3//2])
    if debug:
        print('max peak', max_peak)

    raw_peak_idx = np.where(grad_y_nms > peak_threshold*max_peak)[0]
    raw_peak_idx = raw_peak_idx[(raw_peak_idx >= len(y)/2)*(raw_peak_idx < len(y)*1.5)]

    hpws = get_halfwidth(np.abs(grad_y), raw_peak_idx.tolist())

    peaks = []

    t_oper = np.concatenate((t, t))

    for idx, hpw in zip(raw_peak_idx, hpws):
        peaks.append((t_oper[idx], grad_y[idx], hpw*interval))

# used in debug
    if debug:
        plt.show()
        print('raw peaks are (total' , len(peaks), '):')
        print(peaks)
    
    return sorted(peaks, key=lambda x: x[0])

def get_gaussian_1d(sigma, width):
    """ return 1d gauss blur kernel. width must be odd."""
    half_width = (width-1)//2
    kernel = np.zeros(width)
    for i in range(width):
        kernel[i] = (i-half_width)**2
    return np.exp(-kernel/(2*sigma*sigma))/(2*math.pi*sigma*sigma)


def get_grad_operator(width):
    """ return 1d gradient operator. width must be odd
    """
    half_width = (width-1)//2
    grad_operator = np.zeros(width)
    for i in range(0, half_width):
        grad_operator[i] = - i - 1

    for i in range(half_width + 1, width):
        grad_operator[i] = width - i
    return grad_operator/(width - 1)


def non_max_suppress(data, width):
    """ set non-local-maximum to 0
    """
    half_width = (width+1)//2
    localmax = np.zeros_like(data)
    for i in range(len(data)):
        if data[i] == np.max(data[max(i-half_width,0):min(i+half_width,len(data))]):
            localmax[i] = data[i]
    return localmax

def get_halfwidth(data, peaks):
    hpw = []
    peaks_tmp = [0] + peaks + [len(data)]
    for i in range(1,len(peaks_tmp)-1):
        segment = data[peaks_tmp[i-1]:peaks_tmp[i+1]]
        hpw.append(len(segment[segment > data[i]/2]))
    return hpw

def uniform(array):
    """ Set max to 1, min to 0.
    """
    return (array - np.min(array))/(np.max(array) - np.min(array))

