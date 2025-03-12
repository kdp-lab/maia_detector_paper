#
# This code is based on the 10TeV_plotting_scripts.ipynb notebook.
# It is definitely in need of some cleanup, I have tried to make
# the minimum modifications to get it working. -Jan T. Offermann

import numpy as np
import awkward as ak
import sys,glob
import matplotlib.pyplot as plt
import mplhep as hep
import argparse as ap
# hep.style.use(hep.style.ROOT) # For now ROOT defaults to CMS

from utils.data_utils import load_data_from_json
from utils.plot_utils import choose_histogram,plot_efficiencies
from utils.calc_utils import calculate_efficiencies
from utils.fit_utils import gaussian, double_gaussian, fit_gaussian, fit_double_gaussian, double_gaussian_mean_rms

# Function to process data and calculate RMS values
def process_data(datax, datay, numbins, bins=None,theta=False):
    """
    Process data to calculate RMS values binned by the provided data.

    Parameters:
        datax (list of numpy.ndarray): List of x-data arrays (Usually theta values).
        datay (list of numpy.ndarray): List of y-data arrays (Either pT or d0 resolution).
        numbins (int): How many bins in theta.
        bins (list of numpy.ndarray, optional): List of binning arrays for histograms.

    Returns:
        list: A list of dictionaries containing processed results for each dataset.
            Each dictionary contains bin centers, rms values, and sem values.
    """
    processed_results = []

    if isinstance(bins, np.ndarray):
        bins = [bins] * len(datay)  # Replicate the array for each dataset

    # Loop over the data
    for j in range(len(datay)):
        # Take flattened and transposed data
        data_flatx = np.array(np.ravel(datax[j])).T
        data_flaty = np.array(np.ravel(datay[j])).T
        x_bins = np.linspace(data_flatx.min(), data_flatx.max(), numbins + 1)
        if theta: x_bins = np.linspace(15,165, numbins + 1)
        rms_values = []
        sem_values = []
        bin_centers = []

        # Loop over the theta bins
        for i in range(numbins):
            # Slice the data based on the theta bins
            slice_data = data_flaty[(data_flatx >= x_bins[i]) & (data_flatx < x_bins[i + 1])]
            try:
                # Fit a Gaussian to the slice data using the fit_gaussian function
                popt, pcov, _ = fit_gaussian(slice_data, bins=bins[j])
                fitted_rms = popt[2]
                sem = (np.sqrt(np.diag(pcov)))[2]
            except:
                try:
                    popt, pcov, _ = fit_double_gaussian(slice_data, bins=bins[j])
                    _, fitted_rms, sem = double_gaussian_mean_rms(popt, pcov)
                    if sem == np.inf:
                        sem = np.std(slice_data, ddof=1) / np.sqrt(2 * (len(slice_data) - 1))
                except:
                    fitted_rms = np.sqrt(np.mean(np.square(slice_data - np.mean(slice_data))))
                    sem = np.std(slice_data, ddof=1) / np.sqrt(2 * (len(slice_data) - 1))

            rms_values.append(np.abs(fitted_rms))
            sem_values.append(sem)
            bin_centers.append((x_bins[i] + x_bins[i + 1]) / 2)

        processed_results.append({
            'bin_centers': np.array(bin_centers),
            'rms_values': np.array(rms_values),
            'sem_values': np.array(sem_values),
            'x_err': (data_flatx.max() - data_flatx.min()) / (2 * numbins)
        })

    return processed_results

# Function to fold the data over theta = 90 because detector is symmetric in theta.
# To be used in the case of low statistics (high pT and BIB data)
def fold_data(data, LC_theta_match, LC_pt_match, bib = True):
    # Filter out data points less than 90
    # if bib == False:
    right_side_data = (data[(ak.flatten((LC_theta_match)[(LC_pt_match)>1000])) > 90])
    left_side_data = (data[(ak.flatten((LC_theta_match)[(LC_pt_match)>1000])) <= 90])
    # else:
    #     right_side_data = (data[(ak.flatten((bib_5TeV['LC_track_theta']) > np.pi/2))])
    #     left_side_data = (data[(ak.flatten((bib_5TeV['LC_track_theta']) <= np.pi/2))])
    # if np.array_equal(data, ak.flatten(bib_5TeV['LC_track_theta'][ak.flatten(bib_5TeV['LC_pt_match']>1000)])):
    if np.array_equal(data, ak.flatten(LC_theta_match[ak.flatten(LC_pt_match>1000)])):
        # Mirror the data onto the left side
        mirrored_data_right = np.pi - (right_side_data)
        mirrored_data_left = np.pi - (left_side_data)
        mirrored_data = np.concatenate((mirrored_data_right, mirrored_data_left))
        #folded_data = np.concatenate((data, mirrored_data_right))
        folded_data = np.concatenate((data, mirrored_data))
    elif np.array_equal(data, ak.flatten(LC_theta_match[np.ravel(LC_pt_match)>1000])):
        mirrored_data_right = 180 - (right_side_data)
        mirrored_data_left = 180 - (left_side_data)
        mirrored_data = np.concatenate((mirrored_data_right, mirrored_data_left))
        #folded_data = np.concatenate((data, mirrored_data_right))
        folded_data = np.concatenate((data, mirrored_data))
    else:
        # Concatenate the mirrored data with the original data
        mirrored_data = np.concatenate((right_side_data, left_side_data))
        #folded_data = np.concatenate((data, right_side_data))
        folded_data = np.concatenate((data, mirrored_data))
    return (folded_data)

def plot_processed_data(processed_results, labels=None, xlabel='', ylabel='', title='', fontsize=20, log=False, xlog=False, ylim=None, xlim=None, savepdf=False, bib=False, outdir=None, savename=''):
    """
    Plot the processed RMS data.

    Parameters:
        processed_results (list): List of dictionaries containing processed data.
        labels (list of str, optional): Labels for the data series.
        xlabel (str, optional): X-axis label.
        ylabel (str, optional): Y-axis label.
        title (str, optional): Plot title.
        fontsize (int, optional): Font size for the plot labels.
        log (bool, optional): Use a logarithmic y-axis scale.
        ylim (tuple, optional): Y-axis limits.
        xlim (tuple, optional): X-axis limits.
        savepdf (bool, optional): Save the plot as a PDF file.
        savename (str, optional): Name of the saved PDF file.
    """
    fig, ax = plt.subplots(figsize=(8,6))

    # Loop over the processed results (dictionary of bin centers, rms values, sem values, and theta width) and create errorbar graphs
    for i, result in enumerate(processed_results):
        label = labels[i] if labels else ''
        ax.errorbar(
            result['bin_centers'],
            result['rms_values'],
            xerr=result['x_err'],
            yerr=result['sem_values'],
            fmt='o',
            markersize=6,
            label=label
        )

    # Custom Muon Collider text
    if len(title) > 0:
        ax.set_title(title, fontsize=fontsize)
    else:
        ax.set_title(r'$\it{MAIA}$ Detector Concept', fontsize=fontsize, loc = 'right')

    plt.text(0.04, 0.92, "Muon Collider", fontweight='bold', style='italic', transform=plt.gca().transAxes)
    label = "Simulation"
    if bib: label += ', with BIB'
    else  : label += ', no BIB'
    plt.text(0.04, 0.85, label, transform=plt.gca().transAxes)

    plt.text(0.04, 0.78, r'Lattice v04', transform=plt.gca().transAxes)
    plt.text(0.04, 0.71, r'$\sqrt{s}$ = 10 TeV', transform=plt.gca().transAxes)

    # handle axes
    ax.set_xlabel(xlabel, loc='right', fontsize=fontsize+5)
    ax.set_ylabel(ylabel, loc='top'  , fontsize=fontsize+5)
    if log : ax.set_yscale('log')
    if xlog: ax.set_xscale('log')
    if ylim is not None: ax.set_ylim(ylim)
    if xlim is not None: ax.set_xlim(xlim)
    if labels:
        ax.legend(fontsize=fontsize-2,loc='upper right')
    ax.tick_params(labelsize=fontsize)
    #ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

    if len(savename)>0:
        plt.savefig(savename+".pdf", format='pdf', bbox_inches='tight')

    plt.show()

def main(args):

    hep.style.use(hep.style.ATLAS)

    # Load all the data.

    # TODO: This needs a major rework! Shouldn't be loading all these different hard-coded files, make this an argument!!! -Jan

    # No BIB data
    data_filepath = glob.glob('**/*v2_noBIB_merged.json', recursive=True)[0] # find the v2_noBIB_merged.json file
    v2_noBIB_merged = load_data_from_json(data_filepath)

    # pt_all_path = glob.glob('**/*v0_noBIB_all.json', recursive=True)[0] # find the v0_noBIB_all.json file
    # pt_all = load_data_from_json(pt_all_path)

    # pt_all_5TeV_path = glob.glob('**/*v0_noBIB_all_5TeV.json', recursive=True)[0] # find the v0_noBIB_all_5TeV.json file
    # pt_all_5TeV = load_data_from_json(pt_all_5TeV_path)

    # BIB data
    # bib_all_path = glob.glob('**/*v0_BIB_all.json', recursive=True)[0] # find the v0_BIB_all.json file
    # bib_all = load_data_from_json(bib_all_path)

    # bib_0_50_path = glob.glob('**/*v0_BIB_0_50.json', recursive=True)[0] # find the v0_BIB_0_50.json file
    # bib_0_50 = load_data_from_json(bib_0_50_path)

    # bib_50_250_path = glob.glob('**/*v0_BIB_50_250.json', recursive=True)[0] # find the v0_BIB_50_250.json file
    # bib_50_250 = load_data_from_json(bib_50_250_path)

    # bib_250_1000_path = glob.glob('**/*v0_BIB_250_1000.json', recursive=True)[0] # find the v0_BIB_250_1000.json file
    # bib_250_1000 = load_data_from_json(bib_250_1000_path)

    # bib_5TeV_path = glob.glob('**/*v0_BIB_5TeV.json', recursive=True)[0] # find the v0_BIB_5TeV.json file
    # bib_5TeV = load_data_from_json(bib_5TeV_path)


    # choose_histogram(bib_all, 'nhits')
    # choose_histogram(bib_all, 'chi2_ndf')
    # choose_histogram(bib_all, 'd0')
    # choose_histogram(bib_all, 'pt')
    # choose_histogram(bib_all, 'eta')

    # Convert theta to degrees
    # bib_track_theta = np.degrees(bib_all['LC_track_theta'])
    # bib_truth_theta = np.degrees(2 * np.arctan(np.exp(-bib_all['mcp_mu_eta'])))
    # LC_track_theta = np.degrees(pt_all['LC_track_theta'])
    # mcp_mu_theta = np.degrees(2 * np.arctan(np.exp(-pt_all['mcp_mu_eta'])))


    # ### No BIB

    # Cell to clean (make cuts) and separate data into Barrel/Endcap

    # In[40]:


    which_data = v2_noBIB_merged
    #which_data=pt_all_5TeV

    # Assign variables to no-BIB data
    LC_pt_match = which_data['LC_pt_match']
    mcp_mu_pt = which_data['mcp_mu_pt']
    LC_nhits = which_data['LC_nhits']
    LC_track_theta = np.degrees(which_data['LC_track_theta'])
    mcp_mu_theta = np.degrees(2 * np.arctan(np.exp(-which_data['mcp_mu_eta'])))

    # Clean track using pT >=1 GeV, d0 <= 0.1 mm, and nhits > 4
    track_clean = (ak.flatten(which_data['LC_track_pt'])>=1) & (ak.flatten(which_data['LC_d0'])<= 0.1) & (ak.flatten(which_data['LC_nhits'])>4)

    # Define the eta transition region from barrel to endcap
    transition_region = 1

    # Separate the data into barrel and endcap
    track_barrel = (np.abs(ak.flatten(which_data['LC_eta_match']))<transition_region)
    truth_barrel = (np.abs(ak.flatten(which_data['mcp_mu_eta']))<transition_region)
    track_endcap = (np.abs(ak.flatten(which_data['LC_eta_match']))>=transition_region)
    truth_endcap = (np.abs(ak.flatten(which_data['mcp_mu_eta']))>=transition_region)

    # Sanity check with some overall efficiencies
    # print(len((LC_pt_match[track_barrel])), len((LC_pt_match[track_endcap])), len(LC_pt_match), len(mcp_mu_pt[truth_barrel]), len((mcp_mu_pt[truth_endcap])), len(mcp_mu_pt))
    print("For no BIB:")
    print("Overall efficiency:", len((LC_pt_match))/len((mcp_mu_pt)))
    print("Barrel Efficiency:", len((LC_pt_match[track_barrel]))/len((mcp_mu_pt[truth_barrel])))
    print("Endcap Efficiency:", len((LC_pt_match[track_endcap]))/len((mcp_mu_pt[truth_endcap])))
    print("Lost Efficiency after cleaning:", len((LC_pt_match))/len((mcp_mu_pt)) - len((LC_pt_match)[track_clean])/len((mcp_mu_pt)))

    # Binned in Theta, comparing before and after cleaning

    results, min_value, max_value = calculate_efficiencies([LC_track_theta,LC_track_theta[track_clean]], [mcp_mu_theta,mcp_mu_theta])
    plot_efficiencies(results, min_value, max_value,
                    xlabel=r"$\theta [\degree]$ ",
                    labels=["Before Cleaning", "After Cleaning"],
                    savename='nobib_eff_vs_theta'
                    )

    # Binned in pT

    custom_bins = [1,2,5,10,20,50,100,200,500,1000,2000,5000]
    results, min_value, max_value = calculate_efficiencies([LC_pt_match[track_barrel],LC_pt_match[track_clean & track_barrel]], [mcp_mu_pt [truth_barrel],mcp_mu_pt[truth_barrel]], custom_bins=custom_bins)
    plot_efficiencies(results, min_value, max_value,
                    xlabel="$p_T$ [GeV]",
                    labels=["Before Cleaning", "After Cleaning"],
                    misctext=r'$40^{\circ}<\theta<140^{\circ}$',
                    savename='nobib_eff_vs_pt_barrel'
                    )
    results, min_value, max_value = calculate_efficiencies([LC_pt_match[track_endcap],LC_pt_match[track_clean & track_endcap]], [mcp_mu_pt [truth_endcap],mcp_mu_pt[truth_endcap]], custom_bins=custom_bins)
    plot_efficiencies(results, min_value, max_value,
                    xlabel="$p_T$ [GeV]",
                    labels=["Before Cleaning", "After Cleaning"],
                    misctext=r'$\theta<40^{\circ}$ or $\theta>140^{\circ}$',
                    savename='nobib_eff_vs_pt_endcap'
                    )

    # # # BIB

    # # In[75]:


    # # Assign variables to BIB data
    # bib_track_pt = bib_all['LC_pt_match']
    # bib_truth_pt = bib_all['mcp_mu_pt']
    # bib_nhits = bib_all['LC_nhits']
    # bib_track_eta = bib_all['LC_eta_match']
    # bib_truth_eta = bib_all['mcp_mu_eta']

    # # Clean track using pT >=1 GeV, d0 <= 0.1 mm, and nhits > 4
    # bib_track_clean = (ak.flatten(bib_track_pt)>=1) & (ak.flatten(bib_all['LC_d0'])<= 0.1) & (ak.flatten(bib_nhits)>4)

    # # Define the eta transition region from barrel to endcap
    # transition_region = 1

    # # Separate the data into barrel and endcap
    # bib_track_barrel = (np.abs(ak.flatten(bib_track_eta))<transition_region)
    # bib_truth_barrel = (np.abs(ak.flatten(bib_truth_eta))<transition_region)
    # bib_track_endcap = (np.abs(ak.flatten(bib_track_eta))>=transition_region)
    # bib_truth_endcap = (np.abs(ak.flatten(bib_truth_eta))>=transition_region)

    # # Sanity check with some overall efficiencies
    # # print(len((LC_pt_match[track_barrel])), len((LC_pt_match[track_endcap])), len(LC_pt_match), len(mcp_mu_pt[truth_barrel]), len((mcp_mu_pt[truth_endcap])), len(mcp_mu_pt))
    # print("For BIB:")
    # print("Total efficiency:", len((bib_track_pt))/len((bib_truth_pt)))
    # print("Barrel efficiency:", len((bib_track_pt[bib_track_barrel]))/len((bib_truth_pt[bib_truth_barrel])))
    # print("Endcap efficiency:", len((bib_track_pt[bib_track_endcap]))/len((bib_truth_pt[bib_truth_endcap])))
    # print("Efficiency lost from cleaning:", len((bib_track_pt))/len((bib_truth_pt)) - len((bib_track_pt)[bib_track_clean])/len((bib_truth_pt)))


    # # Binned in Theta

    # # In[76]:


    # results, min_value, max_value = calculate_efficiencies([bib_track_theta,bib_track_theta[bib_track_clean]], [bib_truth_theta,bib_truth_theta])
    # plot_efficiencies(results, min_value, max_value,
    #                 xlabel=r"$\theta [\degree]$ ",
    #                 labels=["Before Cleaning", "After Cleaning"],
    #                 savename='bib_eff_vs_theta',
    #                 bib=True)


    # # Binned in pT

    # # In[79]:


    # custom_bins = [1,2,5,10,20,50,100,200,500,1000]
    # results, min_value, max_value = calculate_efficiencies([bib_track_pt[bib_track_barrel],bib_track_pt[bib_track_clean & bib_track_barrel]], [bib_truth_pt [bib_truth_barrel],bib_truth_pt[bib_truth_barrel]],custom_bins=custom_bins)
    # plot_efficiencies(results, min_value, max_value,
    #                 xlabel="$p_T$ [GeV]",
    #                 labels=["Before Cleaning", "After Cleaning"],
    #                 savename='bib_eff_vs_pt_barrel',
    #                 misctext=r'$40^{\circ}<\theta<140^{\circ}$',
    #                 bib=True)

    # results, min_value, max_value = calculate_efficiencies([bib_track_pt[bib_track_endcap],bib_track_pt[bib_track_clean & bib_track_endcap]], [bib_truth_pt [bib_truth_endcap],bib_truth_pt[bib_truth_endcap]],custom_bins=custom_bins)
    # plot_efficiencies(results, min_value, max_value,
    #                 xlabel="$p_T$ [GeV]",
    #                 labels=["Before Cleaning", "After Cleaning"],
    #                 savename='bib_eff_vs_pt_endcap',
    #                 misctext=r'$\theta<40^{\circ}$ or $\theta>140^{\circ}$',
    #                 bib=True)


    # # Resolutions Plotting

    # ### Setup before plotting
    #
    # Assign the data to variables

    # In[84]:


    which_data = v2_noBIB_merged
    #which_data = pt_all_5TeV

    LC_nhits = which_data['LC_nhits']
    LC_pixel_nhits = which_data['LC_pixel_nhits']
    LC_d0 = which_data['LC_d0']
    LC_z0 = which_data['LC_z0']
    LC_track_pt = which_data['LC_track_pt']
    LC_track_theta = np.degrees(which_data['LC_track_theta'])
    LC_track_eta = which_data['LC_track_eta']
    LC_pt_res = which_data['LC_pt_res']
    LC_pt_match = which_data['LC_pt_match']
    LC_eta_match = which_data['LC_eta_match']
    LC_theta_match = np.degrees(2 * np.arctan(np.exp(-LC_eta_match)))
    LC_phi_match = which_data['LC_phi_match']
    LC_ndf = which_data['LC_ndf']
    LC_chi2 = which_data['LC_chi2']
    LC_dr = which_data['LC_dr']

    # bib_0_50_theta = np.degrees(2*np.arctan(np.exp(-bib_0_50['LC_eta_match'])))
    # bib_50_250_theta = np.degrees(2*np.arctan(np.exp(-bib_50_250['LC_eta_match'])))
    # bib_250_1000_theta = np.degrees(2*np.arctan(np.exp(-bib_250_1000['LC_eta_match'])))

    ptmask1 = np.ravel(LC_pt_match)<=50
    ptmask2 = (np.ravel(LC_pt_match)>50) & (np.ravel(LC_pt_match)<=250)
    ptmask3 = (np.ravel(LC_pt_match)>250) & (np.ravel(LC_pt_match)<=1000)
    ptmask4 = np.ravel(LC_pt_match)>=1000
    # Split the data into 4 pT ranges (and fold the data for the high pT range)
    theta_all    = [LC_theta_match[ptmask1], LC_theta_match[ptmask2], LC_theta_match[ptmask3], fold_data(ak.flatten(LC_theta_match[ptmask4]), LC_theta_match, LC_pt_match, bib = False)]
    d0_all       = [LC_d0[ptmask1]         , LC_d0[ptmask2]         , LC_d0[ptmask3]         , fold_data(ak.flatten(LC_d0[ptmask4])         , LC_theta_match, LC_pt_match, bib = False)]
    pt_res_all   = [LC_pt_res[ptmask1]     , LC_pt_res[ptmask2]     , LC_pt_res[ptmask3]     , fold_data(ak.flatten(LC_pt_res[ptmask4])     , LC_theta_match, LC_pt_match, bib = False)]
    pt_track_all = [LC_track_pt[ptmask1]   , LC_track_pt[ptmask2]   , LC_track_pt[ptmask3]   , fold_data(ak.flatten(LC_track_pt[ptmask4])   , LC_theta_match, LC_pt_match, bib = False)]
    pt_match_all = [LC_pt_match[ptmask1]   , LC_pt_match[ptmask2]   , LC_pt_match[ptmask3]   , fold_data(ak.flatten(LC_pt_match[ptmask4])   , LC_theta_match, LC_pt_match, bib = False)]
    nhits_all    = [LC_nhits[ptmask1]      , LC_nhits[ptmask2]      , LC_nhits[ptmask3]      , fold_data(ak.flatten(LC_nhits[ptmask4])      , LC_theta_match, LC_pt_match, bib = False)]


    # Assign cuts
    theta_all_masked = []
    d0_all_masked = []
    pt_res_all_masked = []
    pt_res2_all_masked = []
    pt_match_all_masked = []
    for i in range(len(theta_all)):
        x = theta_all[i]
        y = d0_all[i]
        w = pt_match_all[i]
        v = pt_track_all[i]
        z = pt_res_all[i] # (truth-reco)/(truth pT) i
        t = pt_res_all[i] /w # (truth-reco)/(truth pT^2)
        u = nhits_all[i]
        # Create a boolean mask for the condition
        theta_cut = (0 <= x) & (x < 180)
        pt_res_cut = np.abs(z) > 0
        pt_cut = v > 1
        d0_cut = np.abs(y) <= 0.1
        nhits_cut = u > 4
        # Apply the mask to filter the arrays
        mask = theta_cut & pt_res_cut & pt_cut & d0_cut & nhits_cut
        x_masked = x[mask]
        y_masked = y[mask]
        z_masked = z[mask]
        w_masked = w[mask]
        t_masked = t[mask]
        theta_all_masked.append(x_masked)
        d0_all_masked.append(y_masked)
        pt_res_all_masked.append(z_masked)
        pt_res2_all_masked.append(t_masked)
        pt_match_all_masked.append(w_masked)
    # theta_all = ak.concatenate([theta_0_50, theta_50_250, theta_match],axis = 0)
    numpoints = 5
    array1 = np.linspace(-1,1,300)    # pt_0_50_bins
    array2 = np.linspace(-0.5,0.5,300)  # pt_50_250_bins
    array3 = np.linspace(-0.1,0.1,300) # pt_250_1000_bins
    array4 = np.linspace(-0.1,0.1,300)  # pt_1000_5000_bins
    d0_bins = [array1, array2, array3, array4] #

    array1 = np.linspace(-0.1,0.1,300)    # pt_0_50_bins
    array2 = np.linspace(-0.1,0.1,300)  # pt_50_250_bins
    array3 = np.linspace(-0.1,0.1,300) # pt_250_1000_bins
    array4 = np.linspace(-0.1,0.1,300)  # pt_1000_5000_bins
    pt_bins = [array1, array2, array3, array4] #

    array1 = np.linspace(-0.003,0.003,300)    # pt_0_50_bins
    array2 = np.linspace(-0.003,0.003,300)  # pt_50_250_bins
    array3 = np.linspace(-0.002,0.002,300) # pt_250_1000_bins
    array4 = np.linspace(-0.002,0.002,300)  # pt_1000_5000_bins
    pt_2_bins = [array1, array2, array3, array4]
    d0_ylim = (0,0.01)
    pt_ylim = (0,0.01)


    # Code for plotting resolutions for no bib

    # In[94]:


    # d0 resolution versus theta
    processed_data = process_data(
        datax=theta_all_masked,
        datay=d0_all_masked,
        numbins=numpoints,
        bins=d0_bins,
        theta=True
    )
    plot_processed_data(
        processed_results=processed_data,
        labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
        xlabel=r'$\theta [\degree]$',
        ylabel=r'$\sigma(d_0)$ [mm]',
        ylim=(0.001,0.1),
        xlim=(0,180),
        log=True,
        savename="res_d0_v_theta_nobib"
    )

    # pt resolution versus theta
    processed_data = process_data(
        datax=theta_all_masked,
        datay=pt_res_all_masked,
        numbins=numpoints,
        bins=pt_bins,
        theta=True
    )
    plot_processed_data(
        processed_results=processed_data,
        labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
        xlabel=r'$\theta[\degree]$',
        ylabel= r'$\sigma(p_T)/p_T$',
        ylim=(0.001,1.0),
        xlim=(0,180),
        log=True,
        savename="res_pt_v_theta_nobib"

    )
    # pt resolution versus pt
    processed_data = process_data(
        datax=pt_match_all_masked,
        datay=pt_res_all_masked,
        numbins=3,
        bins=pt_bins
    )
    plot_processed_data(
        processed_results=processed_data,
        labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
        xlabel=r'$p_T[GeV]$',
        ylabel= r'$\sigma(p_T)/p_T$',
        # title=r'Single $\mu^{\pm}$ no BIB @ 10TeV',
        ylim=(0,0.1),
        savename="res_pt_v_pt_nobib"
    )
    # pt^2 resolution versus theta
    processed_data = process_data(
        datax=theta_all_masked,
        datay=pt_res2_all_masked,
        numbins=numpoints,
        bins=pt_2_bins,
        theta=True
    )
    plot_processed_data(
        processed_results=processed_data,
        labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
        xlabel=r'$\theta[\degree]$',
        ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
        ylim=(0.00001,0.001),
        log=True,
        xlim=(0,180),
        savename="res_pt2_v_theta_nobib"

    )
    processed_data = process_data(
        datax=pt_match_all_masked,
        datay=pt_res2_all_masked,
        numbins=3,
        bins=pt_2_bins
    )
    plot_processed_data(
        processed_results=processed_data,
        labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
        xlabel=r'$p_T[GeV]$',
        ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
        ylim=(0,0.0005),
        xlog=True,
        savename="res_pt2_v_pt_nobib"

    )


    # In[ ]:





    # In[ ]:





    # ### BIB

    # In[90]:


    # bib_theta_all    = [bib_0_50_theta         , bib_50_250_theta         , bib_250_1000_theta         , np.degrees(fold_data(ak.flatten((bib_5TeV['LC_track_theta'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))))]
    # bib_d0_all       = [bib_0_50['LC_d0']      , bib_50_250['LC_d0']      , bib_250_1000['LC_d0']      , fold_data(ak.flatten(bib_5TeV['LC_d0'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))]
    # bib_pt_res_all   = [bib_0_50['LC_pt_res']  , bib_50_250['LC_pt_res']  , bib_250_1000['LC_pt_res']  , fold_data(ak.flatten(bib_5TeV['LC_pt_res'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))]
    # bib_track_pt_all = [bib_0_50['LC_track_pt'], bib_50_250['LC_track_pt'], bib_250_1000['LC_track_pt'], fold_data(ak.flatten(bib_5TeV['LC_track_pt'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))]
    # bib_pt_match_all = [bib_0_50['LC_pt_match'], bib_50_250['LC_pt_match'], bib_250_1000['LC_pt_match'], fold_data(ak.flatten(bib_5TeV['LC_pt_match'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))]
    # bib_nhits_all    = [bib_0_50['LC_nhits']   , bib_50_250['LC_nhits']   , bib_250_1000['LC_nhits']   , fold_data(ak.flatten(bib_5TeV['LC_nhits'][ak.flatten(bib_5TeV['LC_pt_match']>1000)]))]

    # bib_theta_all_masked = []
    # bib_d0_all_masked = []
    # bib_pt_res_all_masked = []
    # bib_pt_res2_all_masked = []
    # bib_pt_match_all_masked = []
    # for i in range(len(bib_theta_all)):
    #     x = bib_theta_all[i]
    #     y = bib_d0_all[i]
    #     w = bib_pt_match_all[i]
    #     z = bib_pt_res_all[i] # truth-reco)/(truth pT)
    #     t = bib_pt_res_all[i] /w # (truth-reco)/(truth pT^2)
    #     v = bib_track_pt_all[i]
    #     u = bib_nhits_all[i]
    #     # Create a boolean mask for the condition
    #     theta_cut = (0 <= x) & (x < 180)
    #     pt_res_cut = np.abs(z) > 0
    #     pt_cut = v > 1
    #     d0_cut = y <= 0.1
    #     nhits_cut = u > 4
    #     # Apply the mask to filter the arrays
    #     mask = theta_cut & pt_res_cut & pt_cut & d0_cut & nhits_cut
    #     x_masked = x[mask]
    #     y_masked = y[mask]
    #     z_masked = z[mask]
    #     w_masked = w[mask]
    #     t_masked = t[mask]
    #     bib_theta_all_masked.append(x_masked)
    #     bib_d0_all_masked.append(y_masked)
    #     bib_pt_res_all_masked.append(z_masked)
    #     bib_pt_res2_all_masked.append(t_masked)
    #     bib_pt_match_all_masked.append(w_masked)
    # # bib_theta_all = ak.concatenate([theta_0_50, theta_50_250, theta_match],axis = 0)
    # numpoints = 5
    # array1 = np.linspace(-0.9,0.9,200)  # pt_0_50_bins
    # array2 = np.linspace(-0.3,0.3,300)  # pt_50_250_bins
    # array3 = np.linspace(-0.4,0.4,300)  # pt_250_1000_bins
    # array4 = np.linspace(-0.2,0.2,300)  # pt_1000_5000_bins
    # d0_bins = [array1, array2, array3, array4] #

    # array1 = np.linspace(-0.1,0.1,100)  # pt_0_50_bins
    # array2 = np.linspace(-0.1,0.1,300)  # pt_50_250_bins
    # array3 = np.linspace(-0.1,0.1,300)  # pt_250_1000_bins
    # array4 = np.linspace(-0.1,0.1,300)  # pt_1000_5000_bins
    # pt_bins = [array1, array2, array3, array4] #

    # array1 = np.linspace(-0.003,0.003,300)    # pt_0_50_bins
    # array2 = np.linspace(-0.002,0.002,300)  # pt_50_250_bins
    # array3 = np.linspace(-0.002,0.002,300) # pt_250_1000_bins
    # array4 = np.linspace(-0.001,0.001,300)  # pt_1000_5000_bins
    # pt_2_bins = [array1, array2, array3, array4]
    # d0_ylim = (0,0.01)
    # pt_ylim = (0,0.01)


    # # In[92]:


    # #d0 res versus theta
    # processed_data = process_data(
    #     datax=bib_theta_all_masked,
    #     datay=bib_d0_all_masked,
    #     numbins=numpoints,
    #     bins=d0_bins,
    #     theta=True
    # )
    # plot_processed_data(
    #     processed_results=processed_data,
    #     labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
    #     xlabel=r'$\theta [\degree]$',
    #     ylabel=r'$\sigma(d_0)$ [mm]',
    #     ylim=(0.001,0.1),
    #     log=True,
    #     xlim=(0,180),
    #     bib=True,
    #     savename="res_d0_v_theta_bib"
    # )
    # # pt resolution versus theta
    # processed_data = process_data(
    #     datax=bib_theta_all_masked,
    #     datay=bib_pt_res_all_masked,
    #     numbins=numpoints,
    #     bins=pt_bins,
    #     theta=True,
    # )
    # plot_processed_data(
    #     processed_results=processed_data,
    #     labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
    #     xlabel=r'$\theta[\degree]$',
    #     ylabel=r'$\sigma(p_T)/p_T$',
    #     ylim=(0.001,1.0),
    #     xlim=(0,180),
    #     log=True,
    #     bib=True,
    #     savename="res_pt_v_theta_bib"
    # )
    # # pt resolution versus pt
    # processed_data = process_data(
    #     datax=bib_pt_match_all_masked,
    #     datay=bib_pt_res_all_masked,
    #     numbins=3,
    #     bins=(pt_bins),
    # )
    # plot_processed_data(
    #     processed_results=processed_data,
    #     labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
    #     xlabel=r'$p_T$ [GeV]',
    #     ylabel=r'$\sigma(p_T)/p_T$',
    #     ylim=(0,0.1),
    #     xlim=(0,5000),
    #     bib=True,
    #     savename="res_pt_v_pt_bib"
    # )
    # # pt2 resolution versus theta
    # processed_data = process_data(
    #     datax=bib_theta_all_masked,
    #     datay=bib_pt_res2_all_masked,
    #     numbins=numpoints,
    #     bins=pt_2_bins,
    #     theta=True,
    # )
    # plot_processed_data(
    #     processed_results=processed_data,
    #     labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
    #     xlabel=r'$\theta[\degree]$',
    #     ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
    #     ylim=(0.00001,0.001),
    #     log=True,
    #     xlim=(0,180),
    #     bib=True,
    #     savename="res_pt2_v_theta_bib"
    # )
    # # pt2 resolution versus pt
    # processed_data = process_data(
    #     datax=bib_pt_match_all_masked,
    #     datay=bib_pt_res2_all_masked,
    #     numbins=3,
    #     bins=(pt_2_bins)
    # )
    # plot_processed_data(
    #     processed_results=processed_data,
    #     labels=[r'$p_T$ = 0-50 GeV', r'$p_T$ = 50-250 GeV', r'$p_T$ = 250-1000 GeV', r'$p_T$ = 1000-5000 GeV'],
    #     xlabel=r'$p_T$ [GeV]',
    #     ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
    #     ylim=(0,0.0003),
    #     xlim=(1,5000),
    #     xlog=True,
    #     bib=True,
    #     savename="res_pt2_v_pt_bib"
    # )


    # In[ ]:





    # # # Extra

    # # ## Check Individual Gaussians
    # #
    # # If some points have large error bars or are not fitting well, check the individual gaussians

    # # In[533]:


    # def plotrms_slice(datax, datay, x_bins, bins=None, xlim = None, title="", rv = False, sigma5 = False):
    #     """
    #     Plot a 1D histogram and fit a Gaussian to it for a specified x-slice.

    #     Parameters:
    #         datax (numpy.ndarray): x-data array.
    #         datay (numpy.ndarray): y-data array.
    #         x_bins (numpy.ndarray): Binning for the x-slice.
    #         bins (numpy.ndarray, optional): Binning for histogram. Default is None.
    #         xlim (tuple, optional): X-axis limits. Default is None.
    #         title (str, optional): Plot title. Default is an empty string.
    #         rv (bool, optional): Return fitted RMS and mean. Default is False.
    #         sigma5 (bool, optional): Plot sigma-5 lines. Default is False.
    #     Returns:
    #         list: List containing fitted RMS and mean if rv=True and a Gaussian fit was successful; otherwise [0, 0].
    #     """
    #     data_flatx = np.array(np.ravel(datax)).T
    #     data_flaty = np.array(np.ravel(datay)).T
    #     for i in range(len(x_bins)-1):
    #         # Select data points within the specified x slice
    #         slice_data = data_flaty[(data_flatx >= x_bins[i]) & (data_flatx < x_bins[i + 1])]
    #         gaussian_fit = True
    #         double_gauss = False
    #         try:
    #             # Fit a Gaussian to the slice data using the fit_gaussian function
    #             popt, pcov, bin_centers = fit_gaussian(slice_data, bins=bins, mean = None, rms = None)
    #             fitted_mean = popt[1]
    #             fitted_rms = popt[2]
    #             # Print the values
    #             # print(f"Slice {i}:", x_bins[i], 'to', x_bins[i + 1])
    #             print("Mean from Fit:", fitted_mean)
    #             print("Sigma from Fit:", fitted_rms)
    #             print("------------------------------------")
    #         except:
    #                 try:
    #                     popt, pcov, bin_centers = fit_double_gaussian(slice_data, bins=bins)
    #                     fitted_mean, fitted_rms = double_gaussian_mean_rms(popt)
    #                     print(f"Slice {i}:", x_bins[i], 'to', x_bins[i + 1])
    #                     print("Mean from Fit:", fitted_mean)
    #                     print("Sigma from Fit:", fitted_rms)
    #                     print("------------------------------------")
    #                     double_gauss = True
    #                 except:
    #                     print(f"Could not fit Gaussian for Slice {i}.")
    #                     gaussian_fit = False
    #         # Plot the 1D histogram with the Gaussian fit
    #         plt.hist(slice_data, bins, alpha=0.5, label='Data')
    #         if gaussian_fit:
    #             if double_gauss:
    #                 plt.plot(bin_centers, double_gaussian(bin_centers, *popt), 'r--', label='Double Gaussian Fit')
    #                 plt.plot(bin_centers, gaussian(bin_centers, *popt[:3]), 'b-', label='Gaussian Fit 1')
    #                 plt.plot(bin_centers, gaussian(bin_centers, *popt[3:]), 'g-', label='Gaussian Fit 2')
    #             else:
    #                 plt.plot(bin_centers, gaussian(bin_centers, *popt), 'r--', label='Gaussian Fit')
    #             if sigma5 == True:
    #                 plt.axvline(x=fitted_mean - 5*np.abs(fitted_rms), linestyle='dotted', label=f'-5$\sigma$ = {fitted_mean - 5*np.abs(fitted_rms)}')
    #                 plt.axvline(x=fitted_mean + 5*np.abs(fitted_rms), linestyle='dotted', label=f'+5$\sigma$ = {fitted_mean + 5*np.abs(fitted_rms)}')
    #         plt.xlabel('\n'+title, loc = 'right')#+ f' for Slice {i}')
    #         plt.ylabel('Counts', loc = 'top')
    #         if xlim is not None:
    #             plt.xlim(xlim)
    #         #plt.yscale('log')
    #         plt.legend()
    #         plt.show()
    #         if gaussian_fit:
    #             if rv == True:
    #                 return fitted_rms, fitted_mean
    #         else:
    #             return [0,0]


    # # In[779]:


    # # Adjust this binning as needed to get the best gaussian, then use that for the resolution plot above
    # d0_bins = np.linspace(-0.3, 0.3, 100) # Seems to work fine using either linspace or just a number
    # pt_bins =  np.linspace(-0.001,0.001, 300) #np.linspace(-0.0003,0.0003, 100) # Use for LC_pt_res/LC_pt_match; Unless there is a very hard cut on the LC_pt_res (< ~1), use linspace between -1,1, but even then not great in endcaps
    # # pt_bins = np.linspace(-0.5,0.5, 100) # Use for LC_pt_res

    # # In case you want to look at other resolutions
    # eta_bins = np.linspace(-2.7, 2.7, 100)
    # chi2ndf_bins = np.linspace(0, 2, 100)
    # nhits_bins = np.linspace(0, 20, 100)

    # # Titles for the plots
    # d0_title = r'$\Delta d_0$'
    # pt_title = r'$\Delta p_T$'


    # # In[780]:


    # # How many data points? Probably make this same as above for consistency
    # nbins = 5
    # theta_bins = np.linspace(15,165,nbins+1)

    # # Specify which pt range is not working well from: bib_0_50, bib_50_250, bib_250_1000, bib_5TeV
    # which_data = bib_250_1000

    # binning = [min(which_data['mcp_mu_pt']), max(which_data['mcp_mu_pt'])]
    # # which_theta = np.degrees(fold_data(bib_250_1000['LC_track_theta'], bib = True)) # THIS HAS TO BE THE SAME AS which_data, also check for degrees vs radians!!!
    # which_theta = np.degrees((which_data['LC_track_theta']))

    # for i in range(nbins):
    #     print(r'Theta:', theta_bins[i], r'<= theta <', theta_bins[i+1])
    #     theta_bin = ((theta_bins[i] <= which_theta) & (which_theta < theta_bins[i+1]) & [np.abs(LC_pt_resolution[0]) < 1 for LC_pt_resolution in (which_data['LC_pt_res'])])
    #     # print(which_data['LC_pt_res'][theta_cut])
    #     count = 0
    #     for bin in theta_bin:
    #         if bin[0] == True:
    #             count +=1
    #     # print("# of data points (total, pt > 1000):", count, len(ak.flatten(fold_data(which_data['LC_pt_match'], bib = True)[theta_bin[fold_data(which_data['LC_pt_match']>1000, bib = True)]]))) # Not sure how 'count' works but it does so don't worry
    #     # print("# of data points (total):", count, len(ak.flatten((which_data['LC_pt_match'])[theta_bin[(which_data['LC_pt_match']>1000)]]))) # Not sure how 'count' works but it does so don't worry

    #     # /pT resolution
    #     # plotrms_slice(which_data['LC_pt_match'][theta_bin], (which_data['LC_pt_res'])[theta_bin], x_bins = x_bins_250_1000, bins=pt_bins, title=pt_title, rv = False, sigma5 = False)

    #     # /pT^2 resolution
    #     plotrms_slice((which_data['LC_pt_match'])[theta_bin], ((which_data['LC_pt_res']/which_data['LC_pt_match']))[theta_bin], x_bins = binning, bins=pt_bins, xlim = None, title=pt_title+r'/$p_T^2$', rv = False, sigma5 = False)

    #     # d0 resolution
    #     # plotrms_slice(which_data['LC_pt_match'][theta_bin], which_data['LC_d0'][theta_bin], x_bins = [1000,5000], bins=pt_bins, xlim = None, title=d0_title, rv = False, sigma5 = False)
    #     # if i == 5:
    #     #     break


    # # In[544]:


    # from matplotlib.colors import LogNorm
    # import seaborn as sns


    # def makeTestPlot(test_data):
    #     LC_d0       = test_data["LC_d0"]
    #     LC_pt_match = test_data["LC_pt_match"]
    #     LC_pt_track = test_data["LC_track_pt"]
    #     LC_nhits    = test_data["LC_nhits"]
    #     LC_chi2     = test_data["LC_chi2"]
    #     LC_ndof     = test_data["LC_ndf"]
    #     LC_chi2ndf = np.ravel(LC_chi2)/np.ravel(LC_ndof)

    #     #plot 0
    #     fig, ax = plt.subplots(figsize=(6, 4))
    #     x = ak.to_numpy(np.ravel(LC_pt_match))
    #     y = ak.to_numpy(np.ravel(LC_chi2ndf))
    #     plt.hist2d(x,y,bins=(100,10),cmin=1,norm=LogNorm(),)

    #     fig, ax = plt.subplots(figsize=(6, 4))

    #     mask = np.ravel(np.abs(LC_d0)<0.1) & np.ravel(LC_nhits>4) & np.ravel(LC_pt_match>1) & np.ravel(LC_pt_track>1) #& np.ravel(LC_chi2<3)

    #     #plot 1
    #     x = ak.to_numpy(np.ravel(LC_pt_match[mask]))
    #     y = ak.to_numpy( (np.ravel(LC_pt_track[mask])-np.ravel(LC_pt_match[mask]) )/np.ravel(LC_pt_match[mask])/np.ravel(LC_pt_match[mask]) )
    #     hist, xe, ye, img = plt.hist2d(x,y,bins=20,range=[[0,5000],[0,0.1]],cmin=0,norm=LogNorm())

    #     #Compute the profiles
    #     # getting the mean and RMS values of each vertical slice of the 2D distribution
    #     # bin width
    #     xbinw = xe[1]-xe[0]
    #     x_slice_mean, x_slice_rms = [], []
    #     for i,b in enumerate(xe[:-1]):
    #         x_slice_mean.append( y[ (x>xe[i]) & (x<=xe[i+1]) ].mean())
    #         x_slice_rms.append( y[ (x>xe[i]) & (x<=xe[i+1]) ].std())

    #     x_slice_mean = np.array(x_slice_mean)
    #     x_slice_rms = np.array(x_slice_rms)

    #     plt.errorbar(xe[:-1]+ xbinw/2, x_slice_mean, x_slice_rms,fmt='_', ecolor='k', color='k')

    #     plt.colorbar()
    #     plt.show()

    #     ax.set_ylabel(r'track $p_{T}$ [GeV]')
    #     ax.set_xlabel(r'truth $p_{T}$ [GeV]')


    #     #plot 2
    #     fig, ax = plt.subplots(figsize=(6, 4))

    #     mask = np.ravel(np.abs(LC_d0)<0.1) & np.ravel(LC_nhits>4) & np.ravel(LC_pt_match>1)

    #     x = ak.to_numpy(np.ravel(LC_pt_match[mask]))
    #     y = ak.to_numpy(np.ravel(LC_d0[mask]))
    #     hist, xe, ye, img = plt.hist2d(x,y,bins=50,range=[[0,5000],[0,0.2]],cmin=0,norm=LogNorm())

    #     #Compute the profiles
    #     # getting the mean and RMS values of each vertical slice of the 2D distribution
    #     # bin width
    #     xbinw = xe[1]-xe[0]
    #     x_slice_mean, x_slice_rms = [], []
    #     for i,b in enumerate(xe[:-1]):
    #         x_slice_mean.append( y[ (x>xe[i]) & (x<=xe[i+1]) ].mean())
    #         x_slice_rms.append( y[ (x>xe[i]) & (x<=xe[i+1]) ].std())

    #     x_slice_mean = np.array(x_slice_mean)
    #     x_slice_rms = np.array(x_slice_rms)

    #     plt.errorbar(xe[:-1]+ xbinw/2, x_slice_mean, x_slice_rms,fmt='_', ecolor='k', color='k')

    #     plt.colorbar()
    #     plt.show()

    #     ax.set_ylabel(r'track $d_0$ [GeV]')
    #     ax.set_xlabel(r'truth $p_{T}$ [GeV]')


    # makeTestPlot(pt_all_5TeV)
    # #makeTestPlot(bib_all)

if(__name__=='__main__'):
    main(sys.argv)