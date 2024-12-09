# Structure:
# data cleaning function
# frequencies
# resconstruct the light curve
# lightcurve
# lightcurve
# magnitude_distance

# imports
import numpy as np
from scipy.signal import lombscargle

# cleanup function

def cleanup(raw_data):
    '''
    Will cleanup the data file in order for better processing
    Includes: Outlier removal, get time to start from 0 (its in julian dates), and center the magnitude around the mean
    
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
    q2 = np.quantile(intensity,.50)
    q3 = np.quantile(intensity,.75)
    # get iqr range
    iqr = q3 - q1
    # get outlier bounds
    lower = q1 - 1.5*iqr
    upper = q3 + 1.5*iqr

    # now we filter
    times = times[(intensity > lower) & (intensity < upper)]
    intensity = intensity[(intensity > lower) & (intensity < upper)]

    # and now we center the data around the mean
    intensity = intensity - np.mean(intensity)

    # and now we return the cleaned up data as an array
    return np.vstack((times,intensity)).T

array = cleanup('cephied_001')
print(array[:10,:])





    




