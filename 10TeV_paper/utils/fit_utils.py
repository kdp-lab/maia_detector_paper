import uuid
import numpy as np
import ROOT as rt
from scipy.optimize import curve_fit

def RN():
    return str(uuid.uuid4())

# Function for creating a Gaussian fit
# def gaussian(x, a, mu, sigma):
#     return a * np.exp(-0.5 * ((x - mu) / sigma)**2)

def gaussian(x,p):
    return p[0] * np.exp(-0.5 * np.square((x[0] - p[1]) / p[2]))

def double_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    return gaussian(a1,mu1,sigma1) + gaussian(a2,mu2,sigma2)

def two_sided_gaussian(x,p):
    if(x[0] >= p[1]):
        sigma = p[2] * (1. + p[3])
    else:
        sigma = p[2] * (1. - p[3])
    return p[0] * np.exp(-0.5 * np.square((x[0] - p[1]) / sigma))

# Function for fitting a Gaussian to the data
def fit_gaussian(slice_data, bins, mean = 0, rms = None, mean_bounds=(-0.001,0.001),debug=False):
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

    if bins is None:
        bins = np.linspace(np.min(slice_data), np.max(slice_data), int(np.sqrt(len(slice_data))))

    h = rt.TH1D(RN(),'',len(bins)-1,bins)
    for entry in slice_data:
        h.Fill(entry)

    if mean is None:
        mean = h.GetMean() #np.mean(slice_data)
    if rms is None:
        rms = h.GetRMS() # NOTE: This is the standard deviation, see ROOT docs! | np.std(slice_data) # np.sqrt(np.mean(np.square(slice_data - mean))) # - mean

    # f = rt.TF1('f_{}'.format(RN()),lambda x, p: gaussian(x,p[0],p[1],p[2]),bins[0],bins[-1],3)
    f = rt.TF1('f_{}'.format(RN()),gaussian,bins[0],bins[-1],3)
    fname = f.GetName()
    f.SetParameter(0, 0.9 * h.GetMaximum())
    # print('Set par0 to {:.2f}'.format(f.GetParameter(0)))
    f.SetParameter(1,mean)
    f.SetParameter(2,0.8 * rms)

    f.SetParLimits(0,0.5 * h.GetMaximum(),1.5 * h.GetMaximum())
    f.SetParLimits(1,*mean_bounds)
    f.SetParLimits(2,0.05 * rms,3. * rms)
    # f = rt.TF1(RN(),'gaus',bins[0],bins[-1])

    initial_parameters = np.array([f.GetParameter(x) for x in range(3)])

    fit_result = h.Fit(fname,'RQSI')
    parameters = np.array([f.GetParameter(x) for x in range(3)])
    uncerts    = np.array([f.GetParError(x) for x in range(3)])

    # print fit result pointer here, causes segmentation violation outside (something necessary deleted?)
    if(debug): fit_result.Print()

    return {
        'fit_result_pointer':fit_result,
        'initial_parameters':initial_parameters,
        'parameters':parameters,
        'uncertainties':uncerts,
        'bins':bins,
        'histogram':h
    }

# Function for fitting a Gaussian to the data -- first step fits with fixed mean, then releases for final fit
def fit_gaussian_two_step(slice_data, bins, mean = 0, rms = None,debug=False):
    #print("test 1")

    if bins is None:
        bins = np.linspace(np.min(slice_data), np.max(slice_data), int(np.sqrt(len(slice_data))))
    #print("test 2")
    h = rt.TH1D(RN(),'',len(bins)-1,bins)
    for entry in slice_data:
        h.Fill(entry)

    if mean is None:
        mean = h.GetMean() #np.mean(slice_data)
    if rms is None:
        rms = h.GetRMS() # NOTE: This is the standard deviation, see ROOT docs! | np.std(slice_data) # np.sqrt(np.mean(np.square(slice_data - mean))) # - mean
    #print("test 3")
    # f = rt.TF1('f_{}'.format(RN()),lambda x, p: gaussian(x,p[0],p[1],p[2]),bins[0],bins[-1],3)
    f = rt.TF1('f_{}'.format(RN()),gaussian,bins[0],bins[-1],3)
    fname = f.GetName()
    f.SetParameter(0, 0.9 * h.GetMaximum())
    # print('Set par0 to {:.2f}'.format(f.GetParameter(0)))
    f.FixParameter(1,mean)
    f.SetParameter(2,0.5 * rms)
    #print("test 4")
    f.SetParLimits(0,0.5 * h.GetMaximum(),1.5 * h.GetMaximum())
    #print("test 4.01")
    # f.SetParLimits(1,*mean_bounds)
    f.SetParLimits(2,0.05 * rms,3. * rms)
    #print("test 4.02")
    # f = rt.TF1(RN(),'gaus',bins[0],bins[-1])
    #print("rms:", rms)
    #print("mean:", mean)

    initial_parameters = np.array([f.GetParameter(x) for x in range(3)])
    #print("test 4.1")
    fit_result = h.Fit(fname,'RQSI')
    #print("test 4.2")
    initial_parameters_2 = np.array([f.GetParameter(x) for x in range(3)])
    #print("test 5")
    f.ReleaseParameter(0)
    f.ReleaseParameter(1)
    f.ReleaseParameter(2)

    fit_option = 'RSI'
    if(not debug): fit_option += 'Q'
    #print("test 5.1")
    fit_result = h.Fit(fname,fit_option)
    #print("test 6")
    parameters = np.array([f.GetParameter(x) for x in range(3)])
    uncerts    = np.array([f.GetParError(x) for x in range(3)])

    return {
        'fit_result_pointer':fit_result,
        'initial_parameters':[initial_parameters,initial_parameters_2],
        'parameters':parameters,
        'uncertainties':uncerts,
        'bins':bins,
        'histogram':h
    }

# Function for fitting a two-sided Gaussian to the data.
def fit_gaussian_two_sided(slice_data, bins, mean = 0, rms = None,debug=False):
    #print("slice data gaussian 2 sided:", fit_gaussian_two_sided)
    #print("bins gaussian 2 sided:", bins)
    if bins is None:
        bins = np.linspace(np.min(slice_data), np.max(slice_data), int(np.sqrt(len(slice_data))))

    h = rt.TH1D(RN(),'',len(bins)-1,bins)
    for entry in slice_data:
        h.Fill(entry)

    if mean is None:
        mean = h.GetMean() #np.mean(slice_data)
    if rms is None:
        rms = h.GetRMS() # NOTE: This is the standard deviation, see ROOT docs! | np.std(slice_data) # np.sqrt(np.mean(np.square(slice_data - mean))) # - mean

    # f = rt.TF1('f_{}'.format(RN()),lambda x, p: gaussian(x,p[0],p[1],p[2]),bins[0],bins[-1],3)
    f = rt.TF1('f_{}'.format(RN()),two_sided_gaussian,bins[0],bins[-1],4)
    fname = f.GetName()
    f.SetParameter(0, 0.9 * h.GetMaximum())
    # print('Set par0 to {:.2f}'.format(f.GetParameter(0)))
    f.FixParameter(1,mean)
    f.SetParameter(2,0.5 * rms)
    f.SetParameter(4,0)

    f.SetParLimits(0,0.5 * h.GetMaximum(),1.5 * h.GetMaximum())
    # f.SetParLimits(1,*mean_bounds)
    f.SetParLimits(2,0.05 * rms,3. * rms)
    f.SetParLimits(3,0,0.9)
    # f = rt.TF1(RN(),'gaus',bins[0],bins[-1])

    initial_parameters = np.array([f.GetParameter(x) for x in range(3)])

    fit_result = h.Fit(fname,'RQSI')

    initial_parameters_2 = np.array([f.GetParameter(x) for x in range(3)])

    f.ReleaseParameter(0)
    f.ReleaseParameter(1)
    f.ReleaseParameter(2)

    fit_option = 'RSI'
    if(not debug): fit_option += 'Q'

    fit_result = h.Fit(fname,fit_option)

    parameters = np.array([f.GetParameter(x) for x in range(3)])
    uncerts    = np.array([f.GetParError(x) for x in range(3)])

    return {
        'fit_result_pointer':fit_result,
        'initial_parameters':[initial_parameters,initial_parameters_2],
        'parameters':parameters,
        'uncertainties':uncerts,
        'bins':bins,
        'histogram':h
    }

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
    # For example, let's use the amplitudes as weights # TODO: Code has lots of things that seem AI-written, need to clean-up/remove. -Jan
    w1 = amp1 / (amp1 + amp2)
    w2 = amp2 / (amp1 + amp2)

    # Then, a simplified weighted uncertainty for the RMS might be calculated as:
    weighted_rms_uncertainty = np.sqrt((w1**2 * rms1_uncertainty**2) + (w2**2 * rms2_uncertainty**2))
    # Calculate representative RMS
    rms = np.sqrt(weighted_variance)
    sem = np.sqrt(weighted_rms_uncertainty)

    return mean, rms, sem