#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
There are some funtions for post processing.
"""

import numpy as np
import sys, os
import shutil
import argparse
from numpy.core.fromnumeric import transpose


def check_path(path_check):
    path_root = os.path.abspath(path_check)
    if (os.access(path_root, os.F_OK) == False
            or os.path.isdir(path_root) == False):
        raise TypeError('Not accessible path.')
    return path_root


def is_digit(string):
    try:
        float(string)
    except:
        return False
    else:
        return True


def str2bool(string):
    if isinstance(string, bool):
        return string
    if string.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif string.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def first(iterable, condition=lambda x: True):
    try:
        return next(i for i, x in enumerate(iterable) if condition(x))
    except StopIteration:
        return None


def item_file_remove(item_path):
    for item in item_path:
        try:
            os.remove(item)
        except FileNotFoundError as e:
            pass
        else:
            pass


def interpolated_intercepts(x, y1, y2):
    """Find the intercepts of two curves, given by the same x data"""
    def intercept(point1, point2, point3, point4):
        """
        find the intersection between two lines
        the first line is defined by the line between point1 and point2
        the first line is defined by the line between point3 and point4
        each point is an (x,y) tuple.
        Returns: the intercept, in (x,y) format
        """
        def line(p1, p2):
            A = (p1[1] - p2[1])
            B = (p2[0] - p1[0])
            C = (p1[0] * p2[1] - p2[0] * p1[1])
            return A, B, -C

        def intersection(L1, L2):
            D = L1[0] * L2[1] - L1[1] * L2[0]
            Dx = L1[2] * L2[1] - L1[1] * L2[2]
            Dy = L1[0] * L2[2] - L1[2] * L2[0]
            if D == 0:
                raise ZeroDivisionError
            x = Dx / D
            y = Dy / D
            return x, y

        L1 = line([point1[0], point1[1]], [point2[0], point2[1]])
        L2 = line([point3[0], point3[1]], [point4[0], point4[1]])
        R = intersection(L1, L2)
        return R

    idxs = np.argwhere(np.diff(np.sign(y1 - y2)) != 0)
    xcs = []
    ycs = []
    for idx in idxs:
        xc, yc = intercept((x[idx], y1[idx]), ((x[idx + 1], y1[idx + 1])),
                           ((x[idx], y2[idx])), ((x[idx + 1], y2[idx + 1])))
        xcs.append(xc)
        ycs.append(yc)
    return np.array(xcs), np.array(ycs)


def obtain_sorted_timeseries(case_path):
    try:
        time_series = list(filter(is_digit, os.listdir(case_path)))
        time_series = list(
            filter(os.path.isdir,
                   [os.path.join(case_path, item) for item in time_series]))
        time_series = [os.path.split(item)[1] for item in time_series]
        time_series.sort(key=lambda item: eval(item))
    except:
        print("Failed to get time series in {0}".format(case_path))
        raise RuntimeError
    return time_series


def obtain_sorted_itemseries(case_path):
    # Obtain both the item and file suffix
    itemFileLists = os.listdir(case_path)
    itemSet = set()
    itemSuffixDict = dict()
    for itemFile in itemFileLists:
        if len(itemFile.split('_', 1)) != 2:
            continue
        itemNow = itemFile.split('_', 1)[0]
        if itemNow not in itemSuffixDict:
            itemSuffixDict.update({itemNow: set()})
        itemSuffixDict[itemNow].add(itemFile.split('_', 1)[1])
    return itemSuffixDict


def obtain_sorted_surfaceitemseries(case_path, suffix='raw'):
    # Obtain both the item and file suffix
    itemFileLists = os.listdir(case_path)
    itemSet = set()
    itemSuffixDict = dict()
    for itemFile in itemFileLists:
        if os.path.splitext(itemFile)[1] == '.' + suffix:
            # Check to see from which item
            if len(os.path.splitext(itemFile)[0].split('_')) < 2:
                continue
            itemNow = os.path.splitext(itemFile)[0].split('_')[-1]
            if itemNow not in itemSuffixDict:
                itemSuffixDict.update({itemNow: set()})
            itemSuffixDict[itemNow].add(
                os.path.splitext(itemFile)[0][0:-len(itemNow) - 1])
    return itemSuffixDict


def moving_derivative(vec, halfwindowlength, order, dt):
    n = np.shape(vec)[0]
    if n != np.size(dt):
        print("The size of vec is not equal to dt.")
        raise RuntimeError
    if n < (halfwindowlength * 2 + 1):
        print("The halfwindowwidth is to large compared to vec.")
        raise RuntimeError
    dvec = np.zeros_like(vec)
    for i in range(halfwindowlength, n - halfwindowlength):
        p = np.polyfit(dt[i - halfwindowlength:i + halfwindowlength],
                       vec[i - halfwindowlength:i + halfwindowlength], order)
        pD = np.array([np.polyder(i, m=1) for i in p.transpose()]).transpose()
        dvec[i] = np.polyval(pD, dt[i])
    # For patch at each end
    p = np.polyfit(dt[0:halfwindowlength], vec[0:halfwindowlength], order)
    pD = np.array([np.polyder(i, m=1) for i in p.transpose()]).transpose()
    for i in range(halfwindowlength):
        dvec[i] = np.polyval(pD, dt[i])

    p = np.polyfit(dt[n - halfwindowlength:n], vec[n - halfwindowlength:n],
                   order)
    pD = np.array([np.polyder(i, m=1) for i in p.transpose()]).transpose()
    for i in range(halfwindowlength, 0, -1):
        dvec[n - i] = np.polyval(pD, dt[n - i])
    return dvec


def computedYdx(x, y, averageWidth=32):
    dYdx = np.zeros_like(x)
    for indexS in range(averageWidth, x.shape[0] - averageWidth - 1, 1):
        dataX = x[indexS - averageWidth:indexS + averageWidth + 1]
        dataY = y[indexS - averageWidth:indexS + averageWidth + 1]
        dYdx[indexS] = (np.mean(dataX * dataY) -
                        np.mean(dataY) * np.mean(dataX)) / (np.mean(
                            dataX * dataX) - np.mean(dataX) * np.mean(dataX))
    return dYdx


def savgol_filter_nonuniform(y, halfwindowlength, polynom, x):
    """
    Applies a Savitzky-Golay filter to y with non-uniform spacing as defined in x

    This is based on https://dsp.stackexchange.com/questions/1676/savitzky-golay-smoothing-filter-for-not-equally-spaced-data
    The borders are interpolated like scipy.signal.savgol_filter would do

    Parameters
    ----------
    y : array_like
        List of floats representing the y values.
    window : int
        Half window length of datapoints. Must be smaller than x
    polynom : int
        The order of polynom used. Must be smaller than the window size
    x : array_like
        List of floats representing the x values of the data. x.size == y.size.

    Returns
    -------
    np.array of float
        The smoothed y values
    """
    n = np.shape(y)[0]
    if y.ndim == 1:
        y = np.reshape(y, (n, 1))
    if y.ndim > 2:
        raise ValueError('The data is not properly prepared.')
    if n != np.size(x):
        print("The size of vec is not equal to dt.")
        raise ValueError('"x" and "y" must be of the same size.')
    if n < (halfwindowlength * 2 + 1):
        print("The halfwindowwidth is to large compared to vec.")
        raise ValueError('The halfwindowwidth is to large compared to vec.')
    if type(polynom) is not int:
        raise TypeError('"polynom" must be an integer')
    window = halfwindowlength * 2 + 1
    if polynom >= window:
        raise ValueError('"polynom" must be less than "window"')
    polynom += 1

    # Initialize variables
    A = np.empty((window, polynom))  # Matrix
    tA = np.empty((polynom, window))  # Transposed matrix
    t = np.empty(window)  # Local x variables
    y_smoothed = np.zeros_like(y)  # np.full(len(y), np.nan)

    # Compute the matrix coefficients
    for i in range(halfwindowlength, n - halfwindowlength, 1):
        # Center a window of x values on x[i]
        for j in range(0, window, 1):
            t[j] = x[i + j - halfwindowlength] - x[i]
        # Create the initial matrix A and its transposed form tA
        for j in range(0, window, 1):
            r = 1.0
            for k in range(0, polynom, 1):
                A[j, k] = r
                tA[k, j] = r
                r *= t[j]
        # Invert the product of the two matrices
        tAA = np.linalg.inv(np.matmul(tA, A))
        # Calculate the pseudoinverse of the design matrix
        coeffs = np.matmul(tAA, tA)

        # Calculate c0 which is also the y value for y[i]
        y_smoothed[i] *= 0
        for j in range(0, window, 1):
            y_smoothed[i] += coeffs[0, j] * y[i + j - halfwindowlength]

        # If at the end or beginning, interpolate the results at the left border
        if i == halfwindowlength:
            first_coeffs = np.zeros((polynom, y.shape[1]))
            for j in range(0, window, 1):
                for k in range(polynom):
                    first_coeffs[k] += coeffs[k, j] * y[j]
        elif i == len(x) - halfwindowlength - 1:
            last_coeffs = np.zeros((polynom, y.shape[1]))
            for j in range(0, window, 1):
                for k in range(polynom):
                    last_coeffs[k] += coeffs[k, j] * y[len(y) - window + j]

    # Interpolate the result at the left border
    for i in range(0, halfwindowlength, 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += first_coeffs[j] * x_i
            x_i *= x[i] - x[halfwindowlength]

    # Interpolate the result at the right border
    for i in range(len(x) - halfwindowlength, len(x), 1):
        y_smoothed[i] = 0
        x_i = 1
        for j in range(0, polynom, 1):
            y_smoothed[i] += last_coeffs[j] * x_i
            x_i *= x[i] - x[-halfwindowlength - 1]
    return y_smoothed
