import numpy as np
import matplotlib.pyplot as plt
from io import open
from scipy import optimize
from scipy import stats
import os

"""
This file automates the process of analyzing a beam. It fits a Gaussian to the data, plots and saves the visibility
data with the fit, plots and fits a line to the full width half maximum versus frequency plots, categorizes the baseline
into "good" and "bad" baselines, ...
"""

def main(cwd, identifier, timestamps, frequencies, chnno_to_feeds, remove_seconds):
    """
    Main runs the program that does the beam analysis. identifier is a String that is the unique identifier for data
    taken on a specific source. (in the case of the M1 2018 data it was 3srcNP). timestamps is a list of hours for the
    data, as a list of Strings. (for M1, the timestamps are 20180101234415 to 20180102014415). Frequencies is a list of
    the frequencies for which there is data, as floats. chnno_to_feeds is a list with two indices, the first being
    the chnno, with the second being another list of length two that has the pair of feed numbers. remove_seconds is
    a list of the seconds that should be removed, in order to isolate the beam (remove calibration and extraneous data)
    """
    # First step is to create the gaussian and fwhm plots and save them to text files.
    setup(cwd)
    error_stream = ''
    available_baselines = [] #keeps track of which baselines have files for vis vals
    fwhmfile = open(cwd+'/FWHM/fwhm_analysis.txt', 'w') # here I write the fwhm pvals and analysis
    outstring = '' #this keeps track of the fwhm analysis information
    for bl in range(1,528):
        print(bl) #helps keep track of how far program has run
        try:
            ampl = get(cwd, identifier, timestamps, bl, remove_seconds)
        except IOError:
            #the baseline either wasn't found, or wasn't supposed to be analyzed
            error_stream += 'Baseline '+ str(bl) + ' was not found \n'
        else: #I think the else block only runs if except doesn't
            available_baselines.append(bl)
            #trying to remove calibrations
            temp = []
            for pt in range(len(ampl)):
                if pt not in remove_seconds:
                    temp.append(ampl[pt])
            ampl = temp
            #this averages the visibilities in bins of 20 frequencies.
            avg = avg_frequencies(ampl,frequencies)

            median = get_median_frequencies(frequencies)
            if not os.path.exists(cwd+'/Gaussians/'+str(bl)):
                os.mkdir(cwd+'/Gaussians/'+str(bl))
            fits,gaussians,errors,times = fit_data(avg, remove_seconds,median)
            outfile = open(cwd+'/Gaussians/'+str(bl)+'/gaussians'+str(bl), 'w')
            #this is the file to write the baselines gaussian data in.
            #TODO: finish the explanation of what the numbers mean
            outfile.write(unicode('a*np.exp((-((x-b)/c)**2)/2.0)+d... the numbers are a,b,c,d,' + '\n'))
            for g in gaussians:
                outfile.write(unicode(str(g[0])+' '+str(g[1])+' '+str(g[2])+' '+str(g[3])+'\n'))
            plot_gaussians(times, avg, fits, cwd+'/Gaussians/'+str(bl), bl,median)
            outstring = outstring + plot_fwhm(cwd, bl, gaussians, median)

            #TODO: make error text file
            #TODO: make fwhm plots
    # here I write the fwhm information
    fwhmfile.write(unicode(outstring))
    # Second step is to analyze the data for characteristics and then make the maps for them.

def plot_fwhm(cwd, bl, gaussians, freqindex):
    """
    Plots the fwhm values for a baseline in seconds, saving it to the FWHM folder.
    """
    #TODO: convert from seconds to degrees
    filename = cwd+'/FWHM/'+str(bl)
    #I write the output file for the chi squared info outside of this method
    fwhms = []
    for frequency in gaussians:
        fwhms.append(frequency[2]*2* np.sqrt(2 * np.log1p(2)))
    plt.clf()
    plt.scatter(freqindex, fwhms)
    params, _ = optimize.curve_fit(line, freqindex, fwhms)
    fitline = []
    for fr in freqindex:
        fitline.append(line(fr,params[0], params[1]))
    plt.plot(freqindex, fitline, label = 'y='+ str(params[0]) + 'x+' + str(params[1]))
    chisq, pval = stats.chisquare(fwhms, fitline)
    plt.legend()
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Time (seconds)')
    plt.savefig(filename)
    outstring = str(bl)+' '+str(params[0])+' '+str(params[1])+' '+str(chisq)+' '+str(pval)+'\n'
    return outstring


def plot_gaussians(times, msrd, fit, filepath, baseline, median_frequencies):
    median=median_frequencies
    #something going wrong with the indices
    for x in range(len(msrd)):
        plt.plot(times,msrd[x], label='Measured Signal')
        plt.xlabel('Time (seconds)')
        plt.ylabel('Visibility')
        if fit[x]!=[0,0,0,0]:
            plt.plot(times, fit[x], label='Fitted Data')
        plt.legend()
        plt.savefig(filepath+'/'+str(baseline)+'_'+str(median[x])+'MHz.png')
        plt.gcf().clear()
        plt.clf()

def fit_data(data, remove_seconds, median_frequencies):
    """
    Takes the data, and returns both the gaussian, and the fitted points for the gaussian.
    """
    bound = [[0,0,0,-np.inf],[np.inf, 100000.0,100000.0,np.inf]] #no initial points works well
    fits =[] #keeps the data points for the fits
    gaussians = [] #keeps the parameters for the gaussians
    times = []
    errors = [] #if optimize could not find a fit, this list will keep track of itself.
    for i in range(len(data[1])):
        times.append(i)
    n=0
    median=median_frequencies
    for fr in data:#i think data
        try:
            params,_ = optimize.curve_fit(gaussian, times, fr,bounds=bound)
        except RuntimeError:
            errors.append(median[n])
            #any flat lines would indicate an issue
            params = [0,0,1,0] #the one prevents a divide by zero error
        finally:
            gaussians.append(params)
            fitted_data=[]
            for t in times:
                fitted_data.append(gaussian(t,params[0],params[1],params[2],params[3]))
            fits.append(fitted_data)
            n+=1
    return fits,gaussians,errors,times

def get_median_frequencies(frequencies):
    median = []
    for x in range(6,len(frequencies)-20,20): #test if stop is correct
        median.append(int(frequencies[x]))
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

def setup(basepath):
    """
    This method sets up the current working directory so that the output is organized as intended.
    """

    # this makes the directories FWHM, Gaussians, Maps
    if not os.path.exists(basepath+'/FWHM'):
        os.mkdir(basepath + '/FWHM')
        os.mkdir(basepath+'/Gaussians')
        os.mkdir(basepath + '/Maps')
    return basepath


def get(cwd, idf, timestamps, chnno,remove_seconds):
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

def baseline_to_feeds(cwd):
    """
    Returns two lists with corresponding Chnno -> feed no. pairs.
    """
    Chnfeed = open(cwd+'/dish_bl_order.txt', 'r')
    lines=list(Chnfeed)
    feeds = []
    Chnnos = []
    for val in range(0, len(lines), 2):
        lines[val]=str(lines[val][1:-3])
        lines[val] = lines[val].split()
        for i in range(len(lines[val])):
            lines[val][i]=int(lines[val][i])
        feeds.append(lines[val])
    for val in range(1, len(lines), 2):
        lines[val]=int(str(lines[val][:-2]))
        Chnnos.append(lines[val])
    return feeds,Chnnos

def get_data_from_file(cwd, idf, start, stop, chnno, i_or_r):
    """
    This is a helper method for get(...) that returns the values in one specific file, which is specified by the parameters of the
    method. Baseline is a string or integer in the ChnNo form (see Santanu's ChnNo->Baseline file). hour is an integer
    that tells which hour's file (since the data is split into hour text files) to access (0 or 1 in the case of the M1 signal),
    since the M1 signal only appears in two hours there won't be an array for accessing any of the other hours.
    i_or_r tells whether to access the imaginary or real values (a string). Must be an 'I' or an 'R'.
    Returns data as a double array, first index is the second(time) the second is the frequency.
    """

    path=cwd+'/visibilities/'+i_or_r.capitalize()+'_'+ str(chnno)+'_'+idf+'_'+start+'_'+stop+'.txt'
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
