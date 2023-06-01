# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:09:12 2022

@author: richy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import math
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy.fft as fft

#%%
#Smoothing function since python doesn't seem to have a built in one
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    import numpy as np
    from math import factorial
    
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

#Check that two arrays are the same length and shorten the longer one if not
def LengthCheck(array,X,xcm):
    if array.shape[1] != np.array([X]).shape[1]: #Sometimes the number of indices between minWN and maxWN varies by 1. This fixes that
        if len(array) < np.array([X]).shape[1]:
            X = X[:-1]
            xcm = xcm[:-1]
        else:
            array = array[:,:-1]
    return array,X,xcm

#Baseline correction; subtracting a 3rd degree polynomial fit instead of just the mean because that's what Eric did
def Baseline(pos,intfgm):
    pos_corr = pos-np.min(pos)
    p = np.poly1d(np.polyfit(pos_corr,intfgm,3))
    intfgm_corr = intfgm - p(pos_corr)
    return pos_corr,intfgm_corr

#Apodization function; mostly finds important parameters (center burst location and FWHM), then feeds them to another function so we can switch if desired
def cdf(pos,HalfMaxLoc,Stdv):
    return (1/2)*(1+erf((pos-HalfMaxLoc+Stdv)/(np.sqrt(2)*Stdv)))
def ExpDecay(pos,MaxLoc,expTC):
    return np.exp(-(pos-MaxLoc)/expTC)

def Apodization(pos_corr,intfgm_corr):
    #This part just finds the middle of the center burts and the HWHM (called Stdv)
    INTfgm = np.abs(intfgm_corr)
    mask = [i for i in np.arange(1,len(INTfgm)-2,1) if INTfgm[i] > INTfgm[i+1] and INTfgm[i] > INTfgm[i-1]] #Find all local maxima
    INT_smooth = savitzky_golay(INTfgm[mask],11,3)
    POS_smooth = pos_corr[mask]
    centerBurst = [i for i in np.arange(0,len(INT_smooth)-1,1) if INT_smooth[i] > np.max(INT_smooth)/2] #Find all points greater than half the maximum
    centerI = [i for i in np.arange(0,len(INT_smooth)-1,1) if INT_smooth[i]==np.max(INT_smooth)] #Find the index of the maximum
    HalfMaxLoc = POS_smooth[centerBurst[0]]
    MaxLoc = POS_smooth[centerI]
    Stdv = MaxLoc-HalfMaxLoc
    expTC = (np.max(pos_corr)-np.min(pos_corr))/4
    #Now we can use the location and width of the center burst to define the apodization function
    return cdf(pos_corr,HalfMaxLoc,Stdv)*ExpDecay(pos_corr,MaxLoc,expTC)

#%%
#This is the main function; it implements the apodization function and does the fft
def FourierTransform(pos_corr,intfgm_corr,zpad):
    INTfgm = Apodization(pos_corr,intfgm_corr)*intfgm_corr
    #Zero-pad interferogram
    targetL = 2**(math.ceil(math.log2(len(intfgm_corr)))+zpad)
    INTpad = np.append(INTfgm,np.zeros(targetL-len(INTfgm)))
    #Shift interferogram so maximum is at the edge (makes phase unwrapped)
    [centerI] = [i for i in np.arange(0,len(INTpad)-1,1) if INTpad[i]==np.max(INTpad)]
    INTpad_shift = np.roll(INTpad,-centerI)
    #Take the Fourier transform
    Dpos = (np.max(pos_corr)-np.min(pos_corr))/(len(pos_corr)*10**4) #Delta mirror position
    FT_sig = fft.fftshift(fft.fft(INTpad_shift))
    FT_freq = fft.fftshift(fft.fftfreq(FT_sig.size,Dpos*2)) #Multiply by two to ignore the negative frequencies
    return FT_sig,FT_freq