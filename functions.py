# Structure:
# data cleaning function
# frequencies
# resconstruct the light curve
# lightcurve
# lightcurve
# magnitude_distance

# imports
import numpy as np
import nfft
import matplotlib.pyplot as plt

# cleanup function

def Load_cleanup(raw_data):
    '''
    Will cleanup the data file in order for better processing
    Includes: Outlier removal, get time to start from 0 (its in julian dates)
    and time normalization, and center the intensities around the mean
    
    Inputs: raw data file from OGLE
    
    Returns: cleaned up data file
    '''

    # the format of the data is col 1: time, col 2: intensity
    raw_data = np.loadtxt(raw_data)
    
    # get the times
    times = raw_data[:,0] - np.min(raw_data[:,0])   # sld be the first cols minus the starting time taking

    # first we want to do outlier removal for the intensity
    intensity = raw_data[:,1]
    # get quantiles
    q1 = np.quantile(intensity,.25)
    q3 = np.quantile(intensity,.75)
    # get iqr range
    iqr = q3 - q1
    # get outlier bounds
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr

    # now we filter
    times = times[(intensity > lower) & (intensity < upper)]
    intensity = intensity[(intensity > lower) & (intensity < upper)]

    # normalize times to [-0.5, 0.5]
    t_norm = (times - np.min(times)) / (np.max(times) - np.min(times)) - 0.5
    # normalize intensities around 0
    i_norm = intensity - np.mean(intensity)


    # and now we return the cleaned up data as an array
    return np.vstack((t_norm,i_norm,times,intensity)).T

# now 
def frequencies(data):
    '''
    Will take the clean data file in order to determine the max frequencies in it.
    Includes an nonuniform fourier transform and the max frequency obtained from it.
    
    Inputs: the cleaned data file from OGLE (cleaned using cleanup function)
    
    Returns: the most probably frequency, the most probable period, the frequency spectrum
    '''
    # get the times and intensities that have been normalized
    t_norm = data[:,0]
    i_norm = data[:,1]
    times = data[:,2]
    
    # get double the size
    N = t_norm.size * 2
    # and we want to make sure this is even
    N += N % 2 != 0

    # get the list of frequencies from the 
    f_k = nfft.nfft_adjoint(t_norm, i_norm, N)
    freqs = np.arange(-N//2, N//2)

    # we really only care about the frequencies above a power of like 100
    f_k[np.abs(f_k) <= 100] = 0

    plt.plot(freqs, np.abs(f_k))
    plt.savefig('power_graph.png')

    # then we want to get the normalized period
    T_recon_norm = -1 / freqs[np.argmax(np.abs(f_k))]
    # and the unnormalized period
    T_recon = T_recon_norm * (np.max(times) - np.min(times))

    # we'll take the inverse of this to get the unnormlized frequency
    max_f = 1 / T_recon

    # and we return
    return max_f, T_recon


A = Load_cleanup('OGLE-BLG-CEP-001.dat')
f, t = frequencies(A)
print(t)
    







    




