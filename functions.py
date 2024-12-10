import numpy as np
import nfft
import matplotlib.pyplot as plt



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

def frequencies(data, T_range=(1, 100), filter_threshold=0.2, plot_frequencies=False):
    '''
    Calculates a filtered (normalized) frequency spectrum and the most likely Cepheid period

    Implements a non-uniform FFT to calculate the frequency spectrum. Can optionally create a plot of the frequency spectrum
    
    Args:
        data (numpy.ndarray): 2d array of data from the cleanup function, with the following columns:
            0: normalized times
            1: normalized intensities
            2: raw times
            3: raw intensities
        T_range (tuple): tuple of length 2 containing the minimum and maximum periods (in days) to search for. Defaults to (1, 100)
        filter_threshold (float): float between 0 and 1 that determines the threshold of the filter. All frequency components
                                  with strengths less than the maximum strength * filter_threshold are discarded. Defaults to 0.2
        plot_frequencies (bool): Determines whether to plot the frequency spectrum. Defaults to False

    Returns:
        T (float): the most likely Cepheid period, calculated from the maximum frequency component
        f_k (numpy.ndarray): the filtered and normalized frequency spectrum
        N (int): the number of calculated frequency components
    '''
    # get the times and intensities that have been normalized
    t_norm = data[:,0]
    i_norm = data[:,1]
    # get the range of unnormalized times
    time_range = np.max(data[:,2]) - np.min(data[:,2])
    
    # calculate the minimum and maximum normalized frequencies
    f_norm_min = 1 / (T_range[1] / time_range)
    f_norm_max = 1 / (T_range[0] / time_range)

    # the number of frequencies to calculate (including positive and negative)
    # not interested in frequencies greater than the max
    N = int(f_norm_max.max()) * 2
    # get the frequencies
    k = np.arange(-N//2, N//2)

    # calculate the Fourier coefficients for each frequency 
    f_k = nfft.nfft_adjoint(t_norm, i_norm, N)

    # exclude frequencies less than the minimum frequency
    f_k[np.abs(k) < int(f_norm_min.min())] = 0
    # exclude frequencies with Fourier coefficient magnitudes less than the max coefficient * filter_threshold
    f_k[np.abs(f_k) < np.abs(f_k).max() * filter_threshold] = 0

    # calculate the normalized period from the strongest frequency component
    T_norm = -1 / k[np.argmax(np.abs(f_k))]
    # unormalize the period
    T = T_norm * time_range

    # plot the filtered frequency spectrum
    if plot_frequencies:
        # unnormalize the frequencies
        k_unnorm = k / time_range

        # plot the magnitude of the Fourier coefficients using only the positive frequencies
        plt.plot(k_unnorm[k>0], np.abs(f_k[k>0]) / N)
        plt.title('Filtered frequency spectrum ')
        plt.xlabel('Frequencies [Hz]')
        plt.ylabel('Magnitude of Fourier Coefficients')

        plt.savefig('freq_spect.png', bbox_inches='tight')
        plt.close()

    # return the period, the normalized frequency spectrum, and the number of frequencies
    return T, f_k, N

def reconstruct_lightcurve(data, T, f_k, N):
    '''
    Reconstructs the original light curve from the filtered frequency spectrum and plots the phase-folded light curve
    
    Args:
        data (numpy.ndarray): 2d array of data from the cleanup function, with the following columns:
            0: normalized times
            1: normalized intensities
            2: raw times
            3: raw intensities
        T (float): the most likely Cepheid period
        f_k (numpy.ndarray): the filtered and normalized frequency spectrum
        N (int): the number of calculated frequency components
    '''

    # reconstructs the normalized intensities
    i_norm_recon = nfft.nfft(data[:,0], f_k, N).real / N
    # unnormalizes the intensities
    i_recon = i_norm_recon + np.mean(data[:,3])

    # plots the phase-folded light curve
    plt.plot(data[:,2] % T / T, i_recon, '.')
    plt.title('Phase-folded reconstructed light curve')
    plt.xlabel('Phase')
    plt.ylabel('Reconstructed magnitude')

    plt.savefig('reco_light_curve.png', bbox_inches='tight')
    plt.close()

def mag_dist(P, cleaned_data, m = None):
    '''
    Will take the period and the data and calculate the absolute magnitude 
    and distance to the cepheid star using the Period-Luminosity relation and 
    distance modulus. 
    
    Inputs: the period, the cleaned data file, and optional the apparent magnitude
    
    Returns: the absolute magnitude, the distance in parsecs, the distance in ly's
    '''
    
    # start by calculating the absolute magnitude of the star
    M = -2.76 * (np.log10(P) - 1) - 4.218

    # and now we calculate the apparent magnitude
    # which is just the average of the intensities
    if m == None:
        m = np.mean(cleaned_data[:,3])

    # the distance calculation, gives it in parsecs
    d_parsec = 10 ** ((m - M + 5)/5)

    # we'll also return the light years
    d_ly = d_parsec * 3.26156

    return M, d_parsec, d_ly



if __name__ == '__main__':
    # load data
    A = Load_cleanup('data/OGLE-LMC-CEP-1812.dat')

    # calculate period and frequency spectrum
    T, f_k, N = frequencies(A, plot_frequencies=True)
    # plot the phase-folded light curve
    reconstruct_lightcurve(A, T, f_k, N)

    # calculate the distance to the Cepheid using the period-luminosity relation and the distance modulus
    M, dp, dly = mag_dist(T,A)

    print(f'Period: {T:.3f} days')
    print(f'Distance: {int(dp)} pc')
