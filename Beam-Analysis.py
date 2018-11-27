import numpy as np
import matplotlib.pyplot as plt
from io import open
from scipy import optimize
from scipy import stats
from scipy import integrate
import os

"""
This file automates the process of analyzing a beam. It fits a Gaussian to the data, plots and saves the visibility
data with the fit, plots and fits a line to the full width half maximum versus frequency plots, categorizes the baseline
into "good" and "bad" baselines, ...
"""
class Baseline:
    def __init__(self):
        """
        Constructs a baseline object that will contain information on its characteristics.
        """

def main(identifier, timestamps, frequencies, chnno_to_feeds):
    """
    Main runs the program that does the beam analysis. identifier is a String that is the unique identifier for data
    taken on a specific source. (in the case of the M1 2018 data it was 3srcNP). timestamps is a list of hours for the
    data, as a list of Strings. (for M1, the timestamps are 20180101234415 to 20180102014415). Frequencies is a list of
    the frequencies for which there is data, as floats. chnno_to_feeds is a list with two indices, the first being
    the chnno, with the second being another list of length two that has the pair of feed numbers.
    """
    # First step is to create the gaussian and fwhm plots and save them to text files.
    cwd = setup()
    available_baselines = [] #keeps track of which baselines have files for vis vals
    for bl in range(1,528):
        print(bl) #helps keep track of how far program has run
        try:
            ampl = get(cwd, identifier, timestamps, bl)
        except FileNotFoundError:
            #the baseline either wasn't found, or wasn't supposed to be analyzed
        else: #I think the else block only runs if except doesn't
            available_baselines.append(bl)
            #this averages the visibilities in bins of 20 frequencies.
            avg = avg_frequencies(ampl)
            median = get_median_frequencies(frequencies)
            os.mkdir(cwd+) #TODO: create directory to store this baselines info
            # TODO: need method to eliminate calibration
            fit = fit_data(avg)







    # Second step is to analyze the data for characteristics and then make the maps for them.

def fit_data(data):
    """
    Takes the data, and returns both the gaussian, and the fitted points for the gaussian.
    """

    initial = [4.0,3550.0,450,40]
    fits =[] #keeps the data points for the fits
    gaussians = [] #keeps the parameters for the gaussians
    for fr in avg:
        result = optimize.curve_fit(gaussian, times, fr, p0=initial)


def get_median_frequencies(frequencies):
    median = []
    for x in range(6,len(frequencies)-20,20): #test if stop is correct
        median.append(frequencies[x])
    return median

def avg_frequencies(data, freq_list):
    """
    Averages together 20 frequency data in order to reduce noise. Starts on 6th frequency and ends on the 506th (of the 512
    frequencies). Returns 2D array, first index is time, second is frequency (of the ones averaged)
    """
    freq = []
    total = []
    avg = []
    for x in range(6,len(freq_list)-20,20): #check this later; possible cause of errors as stop may not be early enough.
        total = get_frequency(data, x)
        for i in range(1,20):
            freq = get_frequency(data, x+i)
            for j in range(len(total)):
                total[j] += freq[j]
        avg.append(total)
    #divides all items in avg by 20
    avg = [[item/20.0 for item in subl] for subl in avg]
    return avg

def get_frequency(data, index):
    """
    From the given data, formatted like the original text files with 512 columns representing frequencies
    and 3600 rows representing time, will find the data for a frequency specified by the index (0-511).
    """

    amp = []
    for second in data:
        amp.append(second[index])
    return amp

def setup():
    """
    This method sets up the current working directory so that the output is organized as intended.
    """
    basepath = os.getcwd()
    # this makes the directories FWHM, Gaussians, Maps
    os.mkdir(basepath + '/FWHM')
    os.mkdir(basepath+'/Gaussians')
    os.mkdir(basepath + '/Maps')
    return basepath


def get(cwd, idf, timestamps, chnno):
    """
    Returns the amplitudes for all the hours of the specified baseline.
    """
    imag = []
    real = []

    for t in range(len(timestamps)-1):
        i= get_data_from_file(cwd, idf, timestamps[t], timestamps[t+1], str(chnno), 'i')
        r= get_data_from_file(cwd, idf, timestamps[t], timestamps[t+1], str(chnno), 'r')
        # if file isn't found, the method should end here
        imag.extend(i)
        real.extend(r)
    amplitude = []
    for row in range(len(imag)):
        amp_for_one_second = []
        for col in range(len(imag[row])):
            amp = np.sqrt(imag[row][col]**2 + real[row][col]**2)
            amp_for_one_second.append(amp)
        amplitude.append(amp_for_one_second)
    return amplitude



def get_data_from_file(cwd, idf, start, stop, chnno, i_or_r ):
    """
    This is a helper method for get(...) that returns the values in one specific file, which is specified by the parameters of the
    method. Baseline is a string or integer in the ChnNo form (see Santanu's ChnNo->Baseline file). hour is an integer
    that tells which hour's file (since the data is split into hour text files) to access (0 or 1 in the case of the M1 signal),
    since the M1 signal only appears in two hours there won't be an array for accessing any of the other hours.
    i_or_r tells whether to access the imaginary or real values (a string). Must be an 'I' or an 'R'.
    Returns data as a double array, first index is the second(time) the second is the frequency.
    """

    path=cwd+'/visibilities/'+i_or_r.capitalize()+'_'+idf+'_'+start+'_'+stop+'.txt'
    f = open(path, 'r')
    # if file isn't found, this method should end here
    lines = list(f)
    data = []
    for line in lines:
        one_line_as_float = []
        one_line_as_string = line.split()
        for point in one_line_as_string:
            one_line_as_float.append(float(point))
        data.append(one_line_as_float)
    return data


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
