import numpy as np
import matplotlib.pyplot as plt
from io import open
from scipy import optimize
from scipy import stats
import os
from  matplotlib import colors as mcolors
from PIL import Image


#TODO: add fwhm plots to maps
# for presentation present possible issues, suggestions on fixing them.
# make map for dead feeds
# make map for weak signals
"""
This file automates the process of analyzing a beam. It fits a Gaussian to the data, plots and saves the visibility
data with the fit, plots and fits a line to the full width half maximum versus frequency plots, categorizes the baseline
into "good" and "bad" baselines, ...
"""

def main(cwd, identifier, timestamps, frequencies, calibration_min, declination,remove_seconds=[]):
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
    deadlist = [] # keeps track of dead baselines
    error_stream = ''
    fwhmfile = open(cwd+'/FWHM/fwhm_analysis.txt', 'w') # here I write the fwhm pvals and analysis
    outstring = '' #this keeps track of the fwhm analysis information
    pvals = [] #stores the pvalues of the fwhm plots
    feeds, chnnos = baseline_to_feeds(cwd)
    chnno_to_feeds=[feeds,chnnos]
    for_calibration = [] #used to find the calibration times
    for bl in range(1,529):
        print(bl) #helps keep track of how far program has run
        try:
            ampl, imag, real = get(cwd, identifier, timestamps, bl)
        except IOError:
            #the baseline either wasn't found, or wasn't supposed to be analyzed
            error_stream += 'Baseline '+ str(bl) + ' was not found \n'
        else: #I think the else block only runs if except doesn't
            #trying to remove calibrations
            temp = []
            tempi = []
            tempr = []
            temps = []
            seconds = range(len(ampl))
            if not for_calibration: #checks if for calibration is empty
                for i in range(len(ampl)):
                    for_calibration.append(ampl[i][int(len(ampl[i])/2)])
                remove_seconds.extend(find_calibrations(for_calibration, calibration_min)) #will this be an issue in the future??
                plt.plot(for_calibration)
                plt.savefig(cwd+'/InfoStream/raw_data_for_calibration_.png')
                plt.clf()
            
            for pt in range(len(ampl)):
                if pt not in remove_seconds:
                    temp.append(ampl[pt])
                    tempi.append(imag[pt])
                    tempr.append(real[pt])
                    temps.append(seconds[pt]/24.0/3600.0*360.0*np.cos(float(declination)))
            ampl = temp
            imag = tempi
            real = tempr
            degrees= temps
            #this averages the visibilities in bins of 20 frequencies.
            avg = avg_frequencies(ampl,frequencies)
            avgi = avg_frequencies(imag, frequencies)
            avgr = avg_frequencies(real, frequencies)
            median = get_median_frequencies(frequencies)
            degrees=center(avg,degrees) #centers the degrees at the peak
            if not os.path.exists(cwd+'/Gaussians/'+str(bl)):
                os.mkdir(cwd+'/Gaussians/'+str(bl))
            fits,gaussians,errors = fit_data(avg,median,degrees)
            for e in errors:
                error_stream+='Python could not get fit for baseline' +str(bl)+', frequency:' +str(e)
            outfile = open(cwd+'/Gaussians/'+str(bl)+'/gaussians'+str(bl), 'w')
            #this is the file to write the baselines gaussian data in.
            outfile.write(unicode('a*np.exp((-((x-b)/c)**2)/2.0)+d... the numbers are a,b,c,d,' + '\n'))
            for g in gaussians:
                outfile.write(unicode(str(g[0])+' '+str(g[1])+' '+str(g[2])+' '+str(g[3])+'\n'))
            plot_gaussians(degrees, avg, fits, cwd+'/Gaussians/'+str(bl), bl,median)
            plot_ir(cwd+'/Gaussians/'+str(bl), degrees, avgi, avgr, median)
            fwhminfo,fwhmpval = plot_fwhm(cwd, bl, gaussians, median)
            outstring = outstring + fwhminfo
            pvals.append([bl,fwhmpval]) #I add the bl number so I know what the pval corresponds to 
            if np.sum(avg) > 1.0:
                deadlist.append(bl)
 
    # here I write the fwhm information in a text file
    fwhmfile.write(unicode(outstring))
    # Second step is to analyze the data for characteristics and then make the maps for them.
    bllist = []
    for pv in pvals:
        if pv[1] < .90: 
            bllist.append(pv[0])
    feedlist = []
    deadfeeds = []
    for chn in range(len(chnno_to_feeds[0])):
        if int(chnno_to_feeds[1][chn]) in bllist:
            feedlist.append(chnno_to_feeds[0][chn])
        if int(chnno_to_feeds[1][chn]) in deadlist:
            deadfeeds.append(chnno_to_feeds[0][chn])
    make_map(feedlist,cwd, 'pvalbelow90')
    make_map(deadfeeds, cwd, 'dead_baselines')

    efile = open(cwd+'/InfoStream/errorstream.txt', 'w')
    efile.write(unicode(error_stream))
   
def center(signal,degrees):
    maxlist = []
    for freq in signal:
        curmax = 0
        for pt in range(len(freq)):
            if freq[pt] > curmax:
                curmax=freq[pt]
        maxlist.append(freq.index(curmax))
            
    avg =int(np.sum(maxlist)/float(len(maxlist)))
    newcenter = degrees[avg]
    newdegrees = []
    for d in degrees:
        newdegrees.append(d-newcenter)
    return newdegrees

def find_calibrations(data, calibration_min):
    seconds = []
    for s in data:
        if s> calibration_min:
            seconds.append(data.index(s))
    return seconds

def make_map(feedlist, cwd, plotname):
    """
    This method takes the list of feeds given and plots and makes a map displaying the geometry of those
    baselines, as well as the polarizations. cwd is the directory where the plot should go or is the folder
    in which Maps is located. Plotname is the name of the plot that will be used for the file name.
    """
    colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
    colors = colors.keys()
    colorindex = 0
    ax = plt.subplot(111, projection = 'polar')
    ax.set_theta_zero_location('N')
    theta = [-40,0,40,80,120,160,200,240,280,-30,30,90,150,210,270,330]
    radius= [2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,0]
    theta = np.radians(theta)
    ax.scatter(theta, radius, marker = 'o', s=1, c='k')
    base_size = 10
    for pair in feedlist:
        if pair[0]%2==0:
            mark1 = 'o'
            size1= base_size/2
        else:
            mark1 = 'x'
            size1 = base_size
        if pair[1]%2==0:
            mark2='o'
            size2 = base_size/2
        else:
            mark2='x'
            size2=base_size
        # if the same polarizations, should make thicker dashed line
        if mark1 != mark2:
            linewidth = 3
            linestyle = '--'
        else:
            linewidth=1
            linestyle='-'
        n1 = int((pair[0]-1)/2)
        n2 = int((pair[1]-1)/2)
        th = [theta[n1],theta[n2]]
        rad = [radius[n1], radius[n2]]
        ax.plot(th,rad, color= colors[colorindex%len(colors)], linewidth=linewidth, linestyle=linestyle)
        if th[0]==th[1] and rad[0]==rad[1]:

            ax.plot(th[0],rad[0], marker=mark1, markersize=size1,color=colors[colorindex%len(colors)], label = str(pair))
            ax.plot(th[1],rad[1], marker=mark2, markersize=size2,color=colors[colorindex%len(colors)])
        else:
            ax.plot(th[0],rad[0], color=colors[colorindex%len(colors)], label = str(pair))
            ax.plot(th[1],rad[1], color=colors[colorindex%len(colors)])
        colorindex+=1
    if os.path.isdir(cwd+'/Maps'):
        plt.savefig(cwd+'/Maps/'+plotname+'.png')
        images = [Image.open(cwd+'/dishmap.png'), Image.open(cwd+'/Maps/'+plotname+'.png')]
    else:
        plt.savefig(cwd+'/'+plotname+'.png')
        images = [Image.open(cwd+'/dishmap.png'), Image.open(cwd+'/Maps/'+plotname+'.png')]

    widths, heights = zip(*(i.size for i in images))
    totalwidth = np.sum(widths)
    maxheight = np.max(heights)
    new_im = Image.new('RGB', (totalwidth, maxheight))
    x_offset = 0
    for im in images:
        new_im.paste(im, (x_offset,0))
        x_offset += im.size[0]

    if os.path.isdir(cwd+'/Maps'):
        new_im.save(cwd+'/Maps/'+plotname+'.png')
    else:
        new_im.save(cwd+'/'+plotname+'.png')

        

def plot_ir(filepath, times, imag, real, median_frequencies):
    median= median_frequencies
    for i in range(len(imag)):
        plt.clf()
        plt.plot(times, imag[i], label='imaginary')
        plt.plot(times, real[i], label='real')
        plt.legend()
        plt.savefig(filepath+'/'+str(median[i])+'ir.png')
    plt.clf()

def plot_fwhm(cwd, bl, gaussians, freqindex):
    """
    Plots the fwhm values for a baseline in degrees, saving it to the FWHM folder. 
    freqindex is a list of the frequencies that will be plotted on the x-axis.
    Returns the p-values of the chi-squared analysis of the fits for the plots,
    and returns outstring; a string containing the parameters for the fwhm fit, and the
    chi-squared and p-values that can be written to a file, and pvals; a . 
    """
    filename = cwd+'/FWHM/'+str(bl)
    #I write the output file for the chi squared info outside of this method
    fwhms = []
    #take out the corresponding x values
    index = 0
    for frequency in gaussians:
        if frequency[2] != 1:
            fwhms.append(frequency[2]*2* np.sqrt(2 * np.log1p(2)))
        else:
            freqindex.remove(freqindex[index])

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
    plt.ylabel('Time (degrees)')
    plt.savefig(filename)
    
    plt.clf()
    outstring = str(bl)+' '+str(params[0])+' '+str(params[1])+' '+str(chisq)+' '+str(pval)+'\n'
    return outstring, pval

def plot_gaussians(times, msrd, fit, filepath, baseline, median_frequencies):
    """
    Plots the data and the fitted gaussians for each frequency on a baseline. Times is the x-values for the plots,
    in degrees. msrd is the actual data, fit is the fitted data. filepath, baseline, median_frequencies are params
    that establish a naming convention for the plots.
    """
    median=median_frequencies
    for x in range(len(msrd)):
        plt.plot(times,msrd[x], label='Measured Signal')
        plt.xlabel('Time (degrees)')
        plt.ylabel('Visibility')
        if fit[x]!=[0,0,1,0]: #if the fit is this, there was an error trying to find it, so we don't plot it
            plt.plot(times, fit[x], label='Fitted Data')
        plt.legend()
        plt.savefig(filepath+'/'+str(baseline)+'_'+str(median[x])+'MHz.png')
        plt.gcf().clear()
        plt.clf()

def fit_data(data, median_frequencies, degrees):
    """
    Takes the data for a single baseline, and returns both the gaussian's parameters, 
    and the fitted points for the gaussian. Parameters: data [list, 2x]; the data we fit the gaussian to,
    is the data for one baseline. median_frequencies [list ,1x]; used to name the plots properly. degrees[list,1x];
    the x-values that we plot the gaussians to. returns fits [list, 2x], the fitted data points. gaussians [list, 2x];
    the parameters a,b,c,d (see method gaussian(...)) for each frequency, errors [list, 1x]; a list of the
    frequencies for which python could not find a fit.
    """
    bound = [[0,0.0,0,-100.0],[1000.0, 500.0,1000.0,100.0]] #no initial points works well
    fits =[] #keeps the data points for the fits
    gaussians = [] #keeps the parameters for the gaussians
    errors = [] #if optimize could not find a fit, will add info to errors
    n=0
    median=median_frequencies
    for fr in data:
        try:
            params,_ = optimize.curve_fit(gaussian, degrees, fr,bounds=bound)
        except RuntimeError:
            errors.append(median[n])
            #any flat lines would indicate an issue with finding a fit
            params = [0,0,1,0] #the one prevents a divide by zero error
        finally:
            gaussians.append(params)
            fitted_data=[]
            for t in degrees:
                fitted_data.append(gaussian(t,params[0],params[1],params[2],params[3]))
            fits.append(fitted_data)
            n+=1
    return fits,gaussians,errors

def get_median_frequencies(frequencies):
    """
    After averaging the each frequencies' visibilities in bins of 20, this method returns the median frequency
    (as an integer) for each of those bins, when given the frequecies that had data collected for them. 
    """
    median = []
    for x in range(6,len(frequencies)-20,20):
        median.append(int(frequencies[x]))
    return median

def avg_frequencies(data, freq_list):
    """
    Averages together 20 frequency data in order to reduce noise. data is the visibility amplitudes, freq_list
    is the list of frequencies that we have data for.
    """
    freq = []
    total = []
    avg = []
    for x in range(6,len(freq_list)-20,20): 
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
    and 3600 rows representing time, will find the visibilities for a frequency specified by the index (0-511).
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
    if not os.path.exists(basepath+'/Gaussians'):
        os.mkdir(basepath+'/Gaussians')
    if not os.path.exists(basepath+'/Maps'):
        os.mkdir(basepath + '/Maps')
    if not os.path.exists(basepath+'/InfoStream'):
        os.mkdir(basepath+'/InfoStream')

def get(cwd, idf, timestamps, chnno):
    """
    Returns the amplitudes for all the hours of the specified baseline. cwd is the current working directory,
    idf is the identifier for a particular signal, timestamps are the hour timestamps that are used in the
    filename, and chnno is the baseline number.
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

    return amplitude, imag, real

def baseline_to_feeds(cwd):
    """
    Returns two lists with corresponding Chnno -> feed no. pairs. cwd is the current working directory,
    where the text file dish_bl_order should be found, in accordance to the organization in github.
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
    This is a helper method for get(...) that returns the values in one specific file, which is to be combined
    later with the corresponding part. All of the parameters are used to specify the file that we are trying
    to access. cwd is the current working directory, idf is the identifier for the signal, start is the starting
    hour, stop is the stopping hour, chnno is the baseline number, i_or_r specifies the imaginary or real parts.
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
