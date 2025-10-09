import os,uuid
import numpy as np
import awkward as ak
import ROOT as rt

from utils.fit_utils import fit_gaussian, fit_gaussian_two_step, fit_gaussian_two_sided, gaussian, two_sided_gaussian

def RN():
    return str(uuid.uuid4())

def combine_masks(mask_list):
    """
    A bit of a hacky function, but we need to potentially combine
    masks where they were derived from different lists of objects,
    so each is a jagged array but some might have empty entries
    for events where others do not.
    We will simply assume that entries are of length 0 or 1,
    which should be OK for the current usage, and then make
    a flat mask (i.e. single dimension) from this.
    """
    # TODO: Broken
    nevents = len(mask_list[0])
    combined_mask = np.full(nevents,True)

    for i,mask in enumerate(mask_list):
        for j in range(nevents):
            if(len(mask[j]) == 0): continue
            combined_mask[j] *= mask[j,0]
    return combined_mask

def process_data(datax, datay, numbins, x_bins=None, bins=None,theta=False,degrees=False,debug=False,debug_name=None, debug_directory=None, debug_labels=None, debug_xlabel=None):
    """
    Process data to calculate RMS values binned by the provided data.

    Parameters:
        datax (list of numpy.ndarray): List of x-data arrays (Usually theta values).
        datay (list of numpy.ndarray): List of y-data arrays (Either pT or d0 resolution).
        numbins (int): How many bins in theta.
        bins (list of numpy.ndarray, optional): List of binning arrays for histograms.

        Jan: Was this AI-generated? Should rewrite this description, it is unclear.

    Returns:
        list: A list of dictionaries containing processed results for each dataset.
            Each dictionary contains bin centers, rms values, and sem values.
    """
    if(debug and debug_name is None):
        debug_name = 'debug'

    if(debug and debug_directory is None):
        debug_directory = os.getcwd()

    if(debug_directory is not None):
        os.makedirs(debug_directory,exist_ok=True)

    processed_results = []

    if isinstance(bins, np.ndarray):
        bins = [bins] * len(datay)  # Replicate the array for each dataset

    # Loop over the data
    for j in range(len(datay)):
        data_flatx = ak.to_numpy(np.transpose(np.ravel(datax[j])))
        data_flaty = ak.to_numpy(np.transpose(np.ravel(datay[j])))

        # TODO: Deal properly with edge case of len(data_flatx) == 0
        if(len(data_flatx) == 0): continue

        if(x_bins is None):

            x_bins = np.linspace(np.min(data_flatx), np.max(data_flatx), numbins + 1,dtype=float)
            if(theta):
                if(degrees):
                    x_bins = np.linspace(15,165, numbins + 1)
                else:
                    x_bins = np.array([30.*np.pi/180.,
                        40.*np.pi/180.,
                        50.*np.pi/180.,
                        60.*np.pi/180.,
                        70.*np.pi/180.,
                        90.*np.pi/180.,
                        110.*np.pi/180.,
                        120.*np.pi/180.,
                        130.*np.pi/180.,
                        140.*np.pi/180.,
                        150.*np.pi/180.
                        ]
                    )
        numbins = len(x_bins) - 1

        hist = rt.TH1D(RN(),'',numbins,x_bins)

        # Loop over the theta bins # TODO: Theta? I don't think these comments are quite right. -Jan
        for k in range(numbins):
            if(debug):
                print('Bin [{}/{}]'.format(k+1,numbins))

            # Slice the data based on the theta bins
            slice_data = data_flaty[(data_flatx >= x_bins[k]) & (data_flatx < x_bins[k + 1])]

            # try:

            # Gaussian fit
            fit_results = fit_gaussian_two_step(slice_data, bins=bins[j], mean=0,debug=debug)
            assert fit_results['fit_result_pointer'].Status() == 0 # NOTE: If this breaks, we fall back on mean/rms directly from distribution
            popt = fit_results['parameters']
            uncerts = fit_results['uncertainties']
            fitted_rms = popt[2]
            sem = uncerts[2]

            # Two-sided Gaussian fit
            fit_results2 = fit_gaussian_two_sided(slice_data, bins=bins[j], mean=0,debug=debug)
            assert fit_results2['fit_result_pointer'].Status() == 0 # NOTE: If this breaks, we fall back on mean/rms directly from distribution
            popt2 = fit_results['parameters']
            uncerts2 = fit_results['uncertainties']
            fitted_rms2 = popt[2]
            sem2 = uncerts[2]

            use_two_sided = False

            # decide which fit to use, based on which has better reduced chi2
            chi2_2 = fit_results2['fit_result_pointer'].Chi2()/fit_results2['fit_result_pointer'].Ndf()
            chi2 = fit_results['fit_result_pointer'].Chi2()/fit_results['fit_result_pointer'].Ndf()
            if(chi2_2 < chi2):
                fit_results = fit_results2
                popt = popt2
                uncerts = uncerts2
                fitted_rms = fitted_rms2
                sem = sem2
                use_two_sided = True
                if(debug):
                    print('\t\tUsing two-sided Gaussian.')

            if(debug):

                # print('\tNumber of data points: {}'.format(len(slice_data)))
                c = rt.TCanvas('c_{}'.format(RN()),'',800,600)

                dmin = fit_results['bins'][0]
                dmax = fit_results['bins'][-1]

                title = 'j = {}, k = {}'.format(j,k)
                if(debug_labels is not None):
                    debug_label = debug_labels[j].replace('$','')
                    title = '{} | bin {}: {} #in [{:.1f},{:.1f}]'.format(debug_label,k,debug_xlabel,x_bins[k], x_bins[k+1])

                h = rt.TH1D(RN(),title,len(fit_results['bins'])-1,fit_results['bins'])
                for entry in slice_data:
                    h.Fill(entry)
                    c.cd()
                h.Draw('HIST')

                if(use_two_sided):
                    f = rt.TF1('f_{}'.format(RN()),two_sided_gaussian,dmin,dmax,4)
                    for l in range(4):
                        f.SetParameter(l,popt[l])
                else:
                    f = rt.TF1('f_{}'.format(RN()),gaussian,dmin,dmax,3)
                    for l in range(3):
                        f.SetParameter(l,popt[l])
                f.Draw('SAME')
                f.SetLineColor(rt.kRed)
                c.cd()
                f.Draw('SAME')
                f.SetNpx(500)

                if(use_two_sided):
                    f2 = rt.TF1('f2_{}'.format(RN()),two_sided_gaussian,dmin,dmax,4)
                    for l in range(4):
                        f2.SetParameter(l,fit_results['initial_parameters'][0][l]) # NOTE: indexing for fit_gaussian_two_step
                else:
                    f2 = rt.TF1('f2_{}'.format(RN()),gaussian,dmin,dmax,3)
                    for l in range(3):
                        f2.SetParameter(l,fit_results['initial_parameters'][0][l]) # NOTE: indexing for fit_gaussian_two_step
                f2.Draw('SAME')
                f2.SetLineColor(rt.kCyan)
                f2.SetLineStyle(rt.kDotted)
                c.cd()
                f2.Draw('SAME')
                f2.SetNpx(500)

                # write fit parameters on plot
                pave = rt.TPaveText(0.1,0.65,0.3,0.8,'NDC')
                pave.SetTextSize(0.02)
                pave.SetBorderSize(0)
                pave.SetTextFont(102)
                pave.SetFillColorAlpha(rt.kWhite,0.)
                pave.AddText('A = {:.1e} #pm {:.1e}'.format(popt[0],uncerts[0]))
                pave.AddText('#mu = {:.1e} #pm {:.1e}'.format(popt[1],uncerts[1]))
                pave.AddText('#sigma = {:.1e} #pm {:.1e}'.format(popt[2],uncerts[2]))
                if(use_two_sided):
                    pave.AddText('#alpha = {:.1e} #pm {:.1e}'.format(popt[3],uncerts[3]))
                pave.Draw()

                legend = rt.TLegend(0.6,0.6,0.9,0.8)
                legend.SetFillColorAlpha(rt.kWhite,0.)
                legend.SetBorderSize(0)
                legend.AddEntry(f2,'initialized fit','l')
                legend.AddEntry(f,'final fit','l')
                legend.Draw()

                c.Draw()
                c.SaveAs("{}/{}_{}_{}.pdf".format(debug_directory,debug_name,j,k))

            # except:
            #     if(debug):
            #         print('Warning: Failed Gaussian fit.')
            #         print('\t           x_bins = ',x_bins)
            #         print('\t(resolution) bins: min = {:.2e}, max = {:.2e}, nbins = {}'.format(bins[j][0],bins[j][-1],len(bins[j])-1))
            #         print('\tlen(slice_data) = ',len(slice_data))

            #     fitted_rms = np.sqrt(np.mean(np.square(slice_data - np.mean(slice_data))))
            #     sem = np.std(slice_data, ddof=1) / np.sqrt(2 * (len(slice_data) - 1))

            hist.SetBinContent(k+1,np.abs(fitted_rms))
            hist.SetBinError(k+1,sem)
        processed_results.append(hist)
    return processed_results

def calculate_efficiencies(track_data, truth_data, mask_pairs, num_bins=10, min_value=None, max_value=None, custom_bins=None):
    """

    Updated efficiency calculation. A little hard-coded, but I think
    the old version defined efficiency in a somewhat confusing way.
    Here we explicitly compute efficiency before and after track
    cleaning is applied.
    """
    # print('len(track_data) = {}, len(truth_data) = {}, len(track_cleaning) = {}'.format(len(track_data),len(truth_data),len(track_cleaning)))

    # Can assign min and max values to the data, don't have to
    if min_value is None:
        min_value = 0.99 * np.min(np.min(np.ravel(truth_data)))

    if max_value is None:
        max_value = 1.01 * np.max(np.max(np.ravel(truth_data)))

    # Bins based on min and max values
    efficiency_bins = np.linspace(min_value, max_value, num_bins+1,dtype=float)

    # print('efficiency_bins = ',efficiency_bins)
    if custom_bins is not None:
        efficiency_bins = np.array(custom_bins,dtype=float)
        min_value = efficiency_bins[0]
        max_value = efficiency_bins[-1]
        num_bins = len(efficiency_bins) - 1

    results = []
    nevents = len(track_data)
    # print('efficiency_bins = ',efficiency_bins)

    # masks allow us to optionally apply track cleaning as a Boolean array.
    for mask_pair in mask_pairs:
        track_mask, truth_mask = mask_pair
        if(truth_mask is None): # if not masking out via any cuts, just take all events
            truth_mask = np.full(nevents,True)
        if(track_mask is None):
            track_mask = np.full(nevents,True) # NOTE: allowing masks to be different, in case reco has extra cuts (e.g. cleaning)

        numerator   = rt.TH1D(RN(),'',num_bins,efficiency_bins)
        denominator = rt.TH1D(RN(),'',num_bins,efficiency_bins)

        for i in range(nevents): # looping thru awkward arrays
            if(not(truth_mask[i])): continue

            track = track_data[i] # in practice, might be empty -- if we failed to get a track
            muon = truth_data[i][0] # basically assuming this to be of length 1
            # bin_idx = (np.digitize(muon,efficiency_bins) - 1)[0]

            denominator.Fill(muon)

            if(len(track) == 0): continue # will happen if there is no matched track

            if(not(track_mask[i])): continue

            numerator.Fill(muon) # NOTE: filling bin corresponding to the matched *muon* kinematics, not the track kinematics!

        # prevent divide-by-zero errors
        for i in range(denominator.GetNbinsX()):
            if(denominator.GetBinContent(i+1) == 0):
                denominator.SetBinContent(i+1,1)

        efficiency = rt.TEfficiency(numerator,denominator)
        efficiency.SetStatisticOption(rt.TEfficiency.kFNormal) # TODO: Is this right? There are a number of options, frequentist & Bayesian -Jan

        # efficiency = rt.TH1D(numerator)
        # efficiency.Divide(denominator)

        # # NOTE: ROOT might be able to do some of this automatically? Will do by hand, to be safe.
        # for i in range(efficiency.GetNbinsX()):
        #     num = numerator.GetBinContent(i+1)
        #     num_e = numerator.GetBinError(i+1)
        #     denom = denominator.GetBinContent(i+1)
        #     denom_e = denominator.GetBinError(i+1)
        #     eff = efficiency.GetBinContent(i+1)

        #     uncert = 0.
        #     if(num != 0. and denom != 0.):
        #         uncert = eff * np.sqrt( np.square(num_e / num) + np.square(denom_e / denom) )
        #         print('a = {:.2e} +/- {:.2e}'.format(num,num_e))
        #         print('b = {:.2e} +/- {:.2e}'.format(denom,denom_e))
        #         print('\t-> uncertainty = {:.2e}'.format(uncert))
        #     efficiency.SetBinError(i+1,uncert)

        # package things in a list, since this is how old code did it
        results.append(efficiency)
    # results.append((bin_centers, efficiency_clean, uncertainty_clean, bin_widths))
    return results, min_value, max_value