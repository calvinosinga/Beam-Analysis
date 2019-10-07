"""

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy import optimize
import itertools

# def avg_data(data,length):
#     """
    
#     """
#     np.arange()
    


def visibilities(imag, real):
    """
    Takes the imaginary and real gains and returns the corresponding visibility numpy array.
    The array can have any dimensions, as long as imag and real have the same dimensions.
    """ 
    return np.array(list(map(sqrtsumofsq, imag,real)))

def sqrtsumofsq(x,y):
    return np.sqrt(x**2+y**2)
    
def fit_gaus(data, time):
    """
    Fits a gaussian to the given data. Time is the equivalent time values for the gains (can be degrees,
    seconds, minutes, seconds, whichever). Returns the fit
    """
    fit=np.zeros(4)
    try:
        params, _ = optimize.curve_fit(gaussian,time, data)
    except RuntimeError:
        fit= np.array([0,0,1,0])
        
    

# def fwhm(data):

# def check_calibration():

# def map():

# def stat():

def seconds_to_degrees(seconds,declination, observ, newcenter=0):
    """
    Converts the time values of seconds and converts them to degrees, and centers the
    gaussian at zero. Observ should contain the visibilities, and can be a 1D or 2D numpy
    array. The method will find the maximum and use that as the center, or take the average
    index of the maximum if observ is 2D, and use that as the center. The user may also specify
    the center, if that is desired.
    """
    degrees = np.zeros(len(seconds))
    for s in range(len(seconds)):
        degrees[s]=(seconds[s]/24.0/3600.0*360.0*np.cos(float(declination)))
    dim = observ.shape
    maxlist = np.zeros(dim[0])
    if len(dim) !=1:
        for f in range(len(observ)):
            maxlist[f]=np.where(observ[f] == np.max(observ[f]))[0][0]
    else:
        maxlist[0] = np.where(observ[f]== np.max(observ[f]))[0][0]
    avg =int(np.sum(maxlist)/float(len(maxlist)))
    if newcenter==0:
        newcenter = degrees[avg]

    newdegrees = []
    for d in degrees:
        newdegrees.append(d-newcenter)
    return newdegrees
    
def get_column(data, index):
    """
    From the given data, formatted like the original text files with 512 columns representing frequencies
    and 3600 rows representing time, will find the visibilities for a frequency specified by the index (0-511).
    """

    amp = []
    for second in data:
        amp.append(second[index])
    return amp

def gaussian(x,a,b,c,d):
    """
    Returns the functional form of the gaussian, with the parameters in this arrangement: a*exp(-((x-b)/c)**2)+d
    """
    return a*np.exp((-((x-b)/c)**2)/2.0)+d

def line(x,m,b):
    """
    Returns the functional form for a line.
    """
    return m*x+b