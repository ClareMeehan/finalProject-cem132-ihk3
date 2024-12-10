# finalProject-cem132-ihk3



## Description

This project estimates the periods, light curves, and distances to classical Cepheids, a type of variable star. It cleans the data and prepares it for the Non-Uniform Fast Fourier Transform and then performs this transform in order to get the most probable frequencies of the Cepheid’s pulsation. These frequencies are then used to find the period of the cepheid, its light curve, and the distance to it. 
The structure of the project is as follows. There are 4 functions, Load_cleanup, frequencies, reconstruct_lightcurve, and mag_dist in the file functions.py. The Load_cleanup function loads in the data and cleans it up. The frequencies function performs the NFFT and finds the max frequency and the period. The reconstruct_lightcurve function will calculate the lightcurve of the data based on the calculated period and list of frequencies. The mag_dist function will calculate the absolute magnitude of the cepheid and the distance to it. 

Please note: some changes were made since the code plan was submitted. They are as follows:
- Rather than using the lomb-scargle, we elected to use the NFFT package. This decision was made as there aren’t really inverse fft’s from the lomb-scargle, but a quasi-inverse exists for the NFFT package.
- The reconstruct and lightcurve functions were combined into one, as they are highly related and both would be quite short as separate ones. 
- The apparent magnitude was made into an optional argument as it can be estimated from the given data. 


## Background

Cepheid variable stars are a type of star whose luminosities vary in a periodic way, in other words they pulsate according to a certain frequency. This change in luminosity comes from both a change in radius and temperature, and the result is a star whose brightness is measured to vary over time. Because of this, light curves for cepheid stars are periodic and show a regular brightening and dimming of the star. Cepheids also have the property that their magnitude is directly correlated with the log of their period, and so using the measured frequency of a cepheid, one can determine its absolute magnitude (its magnitude at 10 parsecs). Through comparison with its apparent magnitude (its magnitude as measured on Earth), the distance to the cepheid can be determined, and is currently a rung on the “distance ladder”* of incredible accuracy! So, through analysis of the periodic form of time vs. luminosity data of a cepheid star, the magnitude and distance to it can be determined, along with its brightness shown over its phase. 

*The cosmic distance ladder is a term used for the various ‘rungs’ of determining distance to astronomical objects in cosmology. It roughly goes as follows: earth parallax, star parallax, RR Lyrae variable stars, Cepheid Variable, Type Ia Supernova, Hubble constant. 

The Non-Uniform Fast Fourier Transform (NFFT) is a variant of the Fast Fourier Transform (FFT) for time series data that is sampled at irregular intervals. Ordinarily, the FFT performs a Forward Discrete Fourier Transform on evenly-spaced time series data, returning a list of frequencies of the same length. An Inverse FFT (IFFT) can be performed on these frequencies, reconstructing the original signal in the time domain.

These transforms can be modeled as linear transformations. A vector of time domain data is multiplied by a transformation matrix to create a vector in the frequency domain, and a vector in the frequency domain is multiplied by the adjoint of the transformation matrix (i.e. the complex conjugate) to create a vector in the time domain. In the case of even sampling, the complex conjugate is proportional to the inverse (by a factor of the length of the time and frequency vectors). Therefore, the inverse always exists and is easily computable, allowing for the IFFT to exactly undo the FFT.

However, in the case of the NFFT, the transformation matrix does not obey such nice properties. With the NFFT, the number of times no longer has to equal the number of frequencies, as the spacing between points is no longer constant. Therefore, the transformation matrix is often not square, and even when it is, it may not be invertible. Therefore, a canonical inverse to the NFFT does not exist. However, the adjoint of the NFFT (multiplication by the complex conjugate of the transformation matrix), or ANFFT provides a good approximation for an inverse NFFT. Therefore, an irregularly spaced time domain signal can be approximately reconstructed from a frequency signal using the ANFFT.


## Features

In order to make the light curve and determine the period of the cepheid, its frequency must be determined from the data. In order to do this, the data is cleaned up, a non-uniform fast fourier transform is performed, the data is then reconstructed with the adjoint non-uniform fast fourier transform, and the frequency with the strongest power is determined as the frequency of the cepheid and used in the distance calculations/light curves. 

The cleaning of the data works to prepare it for the NFFT function such that the function will yield the best result. The first step of this is to remove the outliers of the intensities through the Interquartile range rule, ie. taking out data points further than 1.5 times the IQR from the IQR out of the data set. The NFFT implementation here requires the time samples to range between -0.5 and 0.5, so the times are normalized as such. To avoid the presence of a DC frequency component, the intensities are normalized such that they are centered around 0. Later in the function they are un-normalized so as to get the final results with the proper scaling. 

After normalization, the frequency spectrum is calculated using the forward NFFT. This is done by using the `nfft.nfft_adjoint` function. (Note: the nomenclature for which transform is forward and which is adjoint differs between authors. In this package, the nomenclature is opposite what was described above. `nfft.nfft_adjoint` computes frequencies from irregularly spaced time data, and `nfft.nfft` reconstructs a time domain signal from frequencies.) The range of frequencies to calculate are constrained by the Cepheid periods we are searching for, which is 1 to 200 days by default. Additionally, frequencies with Fourier coefficients that are less than the filter_threshold times the maximum Fourier coefficient are removed. This filters the frequency spectrum and removes noise from the signal. We also optionally produce a plot of the magnitudes of the Fourier coefficients of the frequency spectrum. (Only positive frequencies are plotted. The intensity data is real valued, meaning the Fourier coefficients of the negative frequencies are the complex conjugates of the positive frequencies and thus have the same magnitude. The DC coefficient at 0 Hz is necessarily 0, as the intensities are normalized, so it is also excluded.) Finally, we take the frequency with the maximum Fourier coefficient to be the frequency of oscillation for the Cepheid, and use it to calculate the period.

The intensity data is reconstructed from the frequency spectrum using the ANNFT. We do this using `nfft.nfft` (which, as mentioned above, is what we refer to as the adjoint transformation). After unnormalizing the times and reconstructed intensities, we plot the phase-folded light curve. This is done by taking the modulus of the times by the period, so that intensities are plotted against their phase. This produces a clear visualization of one cycle of the Cepheid’s oscillation.

The magnitude-distance calculation of the data uses the Leavitt-Law (aka the Period-Luminosity relation) in order to find the absolute magnitude and distance to the cepheid star. First the absolute magnitude is calculated using the Leavitt-Law for classical cepheids, and the apparent magnitude is found through taking the mean of the intensities recorded. From here, the distance modulus (the apparent magnitude minus the absolute magnitude) is used to calculate the distance to the cepheid in parsecs (which is 10^(distance modulus / 5)). We also return the distance in light years. 


## Usage

This program is designed to estimate the absolute magnitude of, distance to, and period of a cepheid variable star. It plots the frequency spectrum of its magnitude periodicity to identify the primary frequencies, and it also plots the light curve of the data phase folded with the calculated period in order to examine the functional form of its oscillation. 

The required python packages are given in [`requirements.txt`](requirements.txt),

To use this program, go to the OGLE database of variable stars and download the I band photometry for a classical cepheid. (Note that V band photometry will work too, but the distance calculations are the most precise in the I band.) This project is intended to work on fundamental-mode classical cepheid variable stars and provides the best results for stars that have a period in the typical cepheid range and works best when the data reflects a clear periodic form (i.e. little noise and precise measurements). The idea of trying to be in the typical cepheid range specifically applies to cepheids with very small periods (usually << 1 day), as this causes the frequency normalized relative to the total observation time to be very large, resulting in poor NFFT performance. So, for the default options, it is best to use a star with a typical period range (1 day < period < 200 days). 

The frequencies function also takes the filter threshold as an input. By default, this filters out frequencies with strengths less than 0.2 * the maximum strength. However, for signals with few data points and/or small periods, the frequency signal is noisier due to less available data and spectral leakage. Therefore, for these functions, a higher threshold (e.g. 0.4) can result in more effective noise reduction.

In general, the usage pattern is as follows:

In `cepheid_analysis.py`:
- Cleaned_data = Load_cleanup(‘file name’)
- Period, list of frequencies, number of frequencies = frequencies(Cleaned_Data)
- reconstruct_lightcurve(Cleaned_data, Period, list of frequencies, number of frequencies)
- Magnitude, dis_parsecs, dis_lys = mag_dist(Period, Cleaned Data)


## Testing

Our testing of this project was on classical cepheids with clearly defined light curves on the OGLE collection of variable stars (IDs: OGLE-BLG-CEP-001, OGLE-LMC-CEP-1812). We then calculated their periods and lightcurves and then compared the value of the period and shape of the lightcurve to the value and lightcurve on the OGLE collection of variable stars for that star. Unfortunately the OGLE database does not list the distance to the star on the star’s page, so a named star was chosen (the second ID) in order to also be able to compare the calculated distance to the accepted distance to that star. The comparisons are listed below:

- 001: True Period: 2.5975725, Calculated Period: 2.5977603. 
- 1812: True Period: 1.3129039, Calculated Period: 1.3131914, True Distance: roughly 49,600 parsecs, Calculated Distance: 40,000 parsecs (40,100 when using the apparent magnitude from OGLE)

As shown from the testing results, the period calculation is very close to the accepted value for the period of the cepheids. The distance calculation is a tiny bit off, but considering the distance too it is just roughly 49,600 (and less than 10,000 parsecs off in distance is considered pretty close), it does a fine enough job calculating the distance. Therefore, our program is quite close to correct.


## Citations

Iwanek, P., Soszyński, I., Kozłowski, S., Poleski, R., Pietrukowicz, P., Skowron, J., ... & Ratajczak, M. (2022). The OGLE Collection of Variable Stars: Nearly 66,000 Mira Stars in the Milky Way. The Astrophysical Journal Supplement Series, 260(2), 46.

Jens Keiner, Stefan Kunis, and Daniel Potts. 2009. Using NFFT 3---A Software Library for Various Nonequispaced Fast Fourier Transforms. ACM Trans. Math. Softw. 36, 4, Article 19 (August 2009), 30 pages. https://doi.org/10.1145/1555386.1555388
