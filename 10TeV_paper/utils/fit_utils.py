import numpy as np
from scipy.optimize import curve_fit

# Function for creating a Gaussian fit
def gaussian(x, a, mu, sigma):
    return a * np.exp(-0.5 * ((x - mu) / sigma)**2)

def double_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    return gaussian(a1,mu1,sigma1) + gaussian(a2,mu2,sigma2)

# Function for fitting a Gaussian to the data
def fit_gaussian(slice_data, bins, mean = 0, rms = 0.01):
    """
    Fit a Gaussian to the input data slice.

    Parameters:
        slice_data (numpy.ndarray): Input data slice.
        bins (numpy.ndarray, optional): Binning for histogram. Default is np.linspace(-1, 1, 300).
        mean (float, optional): Mean value for the Gaussian fit. Default is 0.
        rms (float, optional): RMS value for the Gaussian fit. Default is 0.01.

    Returns:
        tuple: Tuple containing fit parameters (popt), covariance matrix (pcov), and bin centers.
    """
    if mean is None:
        mean = np.mean(slice_data)
    if rms is None:
        rms = np.sqrt(np.mean(np.square(slice_data - mean)))
    if bins is None:
        bins = np.linspace(np.min(slice_data), np.max(slice_data), int(np.sqrt(len(slice_data))))
    counts, bin_edges = np.histogram(slice_data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, counts, p0=[max(counts), mean, rms])
    return popt, pcov, bin_centers

def fit_double_gaussian(slice_data, bins=np.linspace(-0.001, 0.001, 300), mean1=0, rms1=0.0001, mean2=0, rms2=0.0005):
    """
    Fit a double Gaussian to the input data slice.

    Parameters:
        slice_data (numpy.ndarray): Input data slice.
        bins (numpy.ndarray, optional): Binning for histogram. Default is np.linspace(-1, 1, 300).
        mean1 (float, optional): Mean value for the first Gaussian. Default is 0.
        rms1 (float, optional): RMS value for the first Gaussian.
        mean2 (float, optional): Mean value for the second Gaussian. Default is 0.
        rms2 (float, optional): RMS value for the second Gaussian.

    Returns:
        tuple: Tuple containing fit parameters (popt), covariance matrix (pcov), and bin centers.
    """
    if bins is None:
        bins = np.linspace(-2.5*np.sqrt(np.mean(np.square(slice_data))), -2.5*np.sqrt(np.mean(np.square(slice_data))), int(np.sqrt(len(slice_data))))
    counts, bin_edges = np.histogram(slice_data, bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    p0 = [max(counts), mean1, rms1, max(counts) / 2, mean2, rms2]
    popt, pcov = curve_fit(double_gaussian, bin_centers, counts, p0=p0)
    return popt, pcov, bin_centers

# Function to calculate the mean and RMS of a double Gaussian distribution
def double_gaussian_mean_rms(popt, pcov):
    """
    Calculate a representative mean and RMS for a double Gaussian distribution.

    Parameters:
        popt (list or numpy.ndarray): Optimal parameters from the double Gaussian fit.
        Expected order: amp1, mean1, rms1, amp2, mean2, rms2.
    Returns:
        tuple: A tuple containing the representative mean and RMS.
    """
    amp1, mean1, rms1, amp2, mean2, rms2 = popt
    total_amp = amp1 + amp2

    # Calculate weighted mean
    mean = (amp1 * mean1 + amp2 * mean2) / total_amp

    # Calculate weighted variance
    variance1 = rms1 ** 2
    variance2 = rms2 ** 2
    mean_variance = ((amp1 * (mean1 - mean) ** 2 + amp2 * (mean2 - mean) ** 2) / total_amp)
    weighted_variance = (amp1 * variance1 + amp2 * variance2) / total_amp + mean_variance
    rms1_uncertainty = np.sqrt(pcov[2, 2])  # Uncertainty of rms1
    rms2_uncertainty = np.sqrt(pcov[5, 5])  # Uncertainty of rms2

    # And assuming w1 and w2 are your weights for rms1 and rms2:
    # For example, let's use the amplitudes as weights
    w1 = amp1 / (amp1 + amp2)
    w2 = amp2 / (amp1 + amp2)

    # Then, a simplified weighted uncertainty for the RMS might be calculated as:
    weighted_rms_uncertainty = np.sqrt((w1**2 * rms1_uncertainty**2) + (w2**2 * rms2_uncertainty**2))
    # Calculate representative RMS
    rms = np.sqrt(weighted_variance)
    sem = np.sqrt(weighted_rms_uncertainty)

    return mean, rms, sem