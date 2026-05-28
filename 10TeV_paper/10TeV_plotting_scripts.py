#
# This code is based on the 10TeV_plotting_scripts.ipynb notebook.
# It is definitely in need of some cleanup, I have tried to make
# the minimum modifications to get it working. -Jan T. Offermann

import numpy as np
import ROOT as rt
import awkward as ak
import sys
import matplotlib.pyplot as plt
import mplhep as hep
import argparse as ap
import uproot as ur
# hep.style.use(hep.style.ROOT) # For now ROOT defaults to CMS

from utils.data_utils import DataLoader
from utils.plot_utils import Plotter
from utils.calc_utils import calculate_efficiencies,combine_masks, process_data

# Function to fold the data over theta = 90 because detector is symmetric in theta.
# To be used in the case of low statistics (high pT and BIB data)
def fold_data(data, LC_theta_match, LC_pt_match):
    # TODO: I suspect this function is broken or needs adjustment. Why is LC_pt_match being used internally?
    #       In practice, hasn't the input "data" *already* had the mask applied? Then applying it again will break things.
    #       Also question the use of ak.flatten()... - Jan
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

def CreateKeyMapping(infile):
    """
    Some sort of messy stuff, for handling the keys with which data is saved.
    A bit complex for historical reasons, as we are trying to keep the ability
    to load older datasets where things were saved under different keys.
    """
    key_mapping = {
        'JSON': # the original key names (basically an identity mapping for the code below)
            {
            'LC_pt_match':'LC_pt_match',
            'LC_eta_match':'LC_eta_match',
            'LC_nhits':'LC_nhits',
            'LC_d0':'LC_d0',
            'LC_z0':'LC_z0',
            'LC_chi2':'LC_chi2',
            'LC_ndf':'LC_ndf',
            'LC_pt_res':'LC_pt_res',
            'LC_track_pt':'LC_track_pt',
            'LC_track_theta':'LC_track_theta',
            'mcp_mu_pt':'mcp_mu_pt',
            'mcp_mu_eta':'mcp_mu_eta',

            'fake_pt':'fake_pt',
            'fake_eta':'fake_eta',
            'fake_chi2':'fake_chi2',
            'fake_d0':'fake_d0',
            'fake_nhits':'fake_nhits',
            'fake_ndf':'fake_ndf'
            },
        'ROOT': # mapping the old keys to much more understandable ones!
            {
            'LC_pt_match':'lc-matched_mcp_pt',
            'LC_eta_match':'lc-matched_mcp_eta',
            'LC_nhits':'lc-matched_track_nhits',
            'LC_d0':'lc-matched_track_d0',
            'LC_z0':'lc-matched_track_z0residual',
            'LC_muon_matched_hits':'lc-matched_track_muon_matched_hits',
            'LC_fraction_muon_matched_hits':'lc-matched_track_fraction_muon_matched_hits',
            'LC_qoverp_residual':'lc-matched_track_qoverp_residual',
            'LC_chi2':'lc-matched_track_chi2',
            'LC_ndf':'lc-matched_track_ndf',
            'LC_pt_res':'lc-matched_track_ptres',
            'LC_track_pt':'lc-matched_track_pt',
            'LC_track_theta':'lc-matched_track_theta',
            'LC_track_eta':'lc-matched_track_eta',
            'mcp_mu_pt':'mcp_mu_pt',
            'mcp_mu_eta':'mcp_mu_eta',

            'fake_pt':'fake_track_pt',
            'fake_eta':'fake_track_eta',
            'fake_chi2':'fake_track_chi2',
            'fake_d0':'fake_track_d0',
            'fake_nhits':'fake_track_nhits',
            'fake_ndf':'fake_track_ndf'

            }
    }

    # for historical reasons, ROOT-style keys might have either 'lc-matched' and 'dr-matched'
    # substrings, OR 'lc_matched' and 'dr_matched'. (Latter eventually caused processing issues
    # with some updates to SLCIO-Analyzer).
    # We will open the actual input to check this.
    using_dash = True
    if('.root' in infile):
        f = ur.open(infile)
        t = f['ntuple'] # assuming the treename -- should be OK!
        keys = list(t.keys())
        for key in keys:
            if('lc_matched' in key or 'dr_matched' in key):
                using_dash = False
                break

    if(not using_dash):
        for key,val in key_mapping['ROOT'].items():
            # print(key,val)
            new_val = val.replace('lc-matched','lc_matched').replace('dr-matched','dr_matched')
            key_mapping['ROOT'][key] = new_val
        f.close()
    return key_mapping

def main(args):
    upto1TeV = False # FIXME make this an argument
    #upto10GeV = True 
    hep.style.use(hep.style.ATLAS)
    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptStat(0)

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFile',type=str,required=True)
    parser.add_argument('-is','--inputFile_SiTracks',type=str,required=True) # SiTracks input file
    parser.add_argument('-o','--outputDirectory',type=str,default='output_5000GeV_v8_NoBIB')
    parser.add_argument('-dl','--dataLabel',type=str,default='Simulation, BIB')
    parser.add_argument('-ll','--latticeLabel',type=str,default='EU24 Lattice')
    parser.add_argument('-s','--suffix',type=str,default=None)
    parser.add_argument('-degrees','--degrees',type=int,default=0)

    # allow skipping of the efficiency plots -- mostly useful for debugging, to speed things up
    parser.add_argument('-doEfficiency'  ,   '--doEfficiency', type=int, default=1)
    parser.add_argument('-efficiencyOpts', '--efficiencyOpts', type=int, default=1, help='0 for before/after track cleaning, 1 for before only, 2 for after only')
    parser.add_argument('-upTo10GeVbool', '--upTo10GeVbool', type=bool, default=False, help='False for running over [0, 1] or [0, 5] TeV sample, True for running over (0.5, 10) GeV sample')

    args = vars(parser.parse_args())
    infile = args['inputFile']
    infile_SiTracks = args['inputFile_SiTracks']
    outdir = args['outputDirectory']
    data_label = args['dataLabel']
    lattice_label = args['latticeLabel']
    suffix = args['suffix']
    degrees = args['degrees'] > 0
    do_efficiency = args['doEfficiency'] > 0
    efficiency_opts = args['efficiencyOpts']
    upto10GeV = args['upTo10GeVbool']

    if(suffix is None):
        if('no bib' in data_label.lower()):
            suffix = 'nobib'
        else:
            suffix = 'bib'

    key_mapping = CreateKeyMapping(infile)

    key_mapping_selectedtracks = CreateKeyMapping(infile_SiTracks)

    # Load the data. For ROOT files, uproot will handle things so that memory usage is reasonable.
    data_loader = DataLoader()
    data_loader.SetFilename(infile)
    data_loader.Load()

    data_loader_SiTracks = DataLoader()
    data_loader_SiTracks.SetFilename(infile_SiTracks)
    data_loader_SiTracks.Load()

    plotter = Plotter()
    plotter.SetCOMTev(10)
    plotter.SetDataLabel(data_label)
    plotter.SetLatticeLabel(lattice_label)
    plotter.SetOutputDirectory(outdir)

    #plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_pt'], key_mapping[data_loader.GetMode()]['LC_track_pt'], np.logspace(0, 4, 100), 'Track $p_T$ [GeV]', 'Normalized Count', x_range=(1, 5000), y_scale='log')
    plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_pt'], key_mapping[data_loader.GetMode()]['mcp_mu_pt'], np.linspace(0, 10, 10), 'Track $p_T$ [GeV]', 'Normalized Count', x_range=(0, 10), y_scale='log')
    plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_eta'], key_mapping[data_loader.GetMode()]['LC_track_eta'], np.linspace(-2.5,2.5,25), r'Track $\eta$', 'Normalized Count', y_scale='linear')
    plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_chi2'], key_mapping[data_loader.GetMode()]['LC_chi2'], np.linspace(0,3,30), r'Track $\chi^2/n_{dof}$', 'Normalized Count', x_range=(0, 3), y_scale='linear', custom_data_func=lambda d: (ak.flatten(d[key_mapping[data_loader.GetMode()]['fake_chi2']]) / ak.flatten(d[key_mapping[data_loader.GetMode()]['fake_ndf']]), ak.flatten(d[key_mapping[data_loader.GetMode()]['LC_chi2']]) / ak.flatten(d[key_mapping[data_loader.GetMode()]['LC_ndf']])))
    #plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_chi2'], key_mapping[data_loader.GetMode()]['LC_chi2'], np.logspace(0,6,30), r'Track $\chi^2/n_{dof}$', 'Normalized Count', x_range=(0, 3), x_scale='log', y_scale='log', custom_data_func=lambda d: (ak.flatten(d[key_mapping[data_loader.GetMode()]['fake_chi2']]) / ak.flatten(d[key_mapping[data_loader.GetMode()]['fake_ndf']]), ak.flatten(d[key_mapping[data_loader.GetMode()]['LC_chi2']]) / ak.flatten(d[key_mapping[data_loader.GetMode()]['LC_ndf']])))
    plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_d0'], key_mapping[data_loader.GetMode()]['LC_d0'], np.linspace(-6,6,50), r'Track $d_0$ [mm]', 'Normalized Count', x_range=(-6,6), y_scale='log')
    plotter.plot_fake_and_true_distributions(data_loader, key_mapping[data_loader.GetMode()]['fake_nhits'], key_mapping[data_loader.GetMode()]['LC_nhits'], np.arange(-0.5, 26, 1), r'Track $n_{hits}$', 'Normalized Count', y_scale='linear', custom_data_func=lambda d: (ak.flatten(d[key_mapping[data_loader.GetMode()]['fake_nhits']]), ak.flatten(d[key_mapping[data_loader.GetMode()]['LC_nhits']])))
    # PlotHistogram(data, 'fake_phi', 'mcp_mu_phi', (100), r'Track $\phi$', 'Normalized Count', x_range=(-np.pi, np.pi), y_scale='linear')

    # Assign variables
    LC_pt_match = data_loader[key_mapping[data_loader.GetMode()]['LC_pt_match']]
    LC_SiTrack_pt_match = data_loader_SiTracks[key_mapping_selectedtracks[data_loader.GetMode()]['LC_pt_match']]
    LC_eta_match = data_loader[key_mapping[data_loader.GetMode()]['LC_eta_match']]
    LC_theta_match = 2 * np.arctan(np.exp(-LC_eta_match))
    LC_track_pt = data_loader[key_mapping[data_loader.GetMode()]['LC_track_pt']]
    LC_SiTrack_theta = data_loader_SiTracks[key_mapping_selectedtracks[data_loader.GetMode()]['LC_track_theta']]
    LC_SiTrack_eta = data_loader_SiTracks[key_mapping_selectedtracks[data_loader.GetMode()]['LC_track_eta']]
    LC_track_theta = data_loader[key_mapping[data_loader.GetMode()]['LC_track_theta']]
    LC_track_eta = data_loader[key_mapping[data_loader.GetMode()]['LC_track_eta']]
    LC_track_chi2_red = data_loader[key_mapping[data_loader.GetMode()]['LC_chi2']] / data_loader[key_mapping[data_loader.GetMode()]['LC_ndf']]
    mcp_mu_pt = data_loader[key_mapping[data_loader.GetMode()]['mcp_mu_pt']]
    mcp_mu_eta = data_loader[key_mapping[data_loader.GetMode()]['mcp_mu_eta']]
    mcp_mu_theta = 2 * np.arctan(np.exp(-mcp_mu_eta))
    LC_nhits = data_loader[key_mapping[data_loader.GetMode()]['LC_nhits']]
    LC_d0 = data_loader[key_mapping[data_loader.GetMode()]['LC_d0']]
    LC_z0residual = data_loader[key_mapping[data_loader.GetMode()]['LC_z0']]
    LC_muon_matched_hits = data_loader_SiTracks[key_mapping[data_loader.GetMode()]['LC_muon_matched_hits']]
    LC_fraction_muon_matched_hits = data_loader_SiTracks[key_mapping[data_loader.GetMode()]['LC_fraction_muon_matched_hits']]
    LC_qoverp_residual = data_loader_SiTracks[key_mapping[data_loader.GetMode()]['LC_qoverp_residual']]
    LC_pt_res = data_loader[key_mapping[data_loader.GetMode()]['LC_pt_res']]
    LC_pt_res_SiTracks = data_loader_SiTracks[key_mapping_selectedtracks[data_loader.GetMode()]['LC_pt_res']]
    if(degrees):
        LC_theta_match = np.degrees(LC_theta_match)
        LC_track_theta = np.degrees(LC_track_theta)
        mcp_mu_theta = np.degrees(mcp_mu_theta)

    # Defining cleaning differently here, as an awkward array
    pt_cut_value = 1. # GeV
    d0_cut_value = 0.1 # mm?
    nhits_cut_value = 4
    chi2_red_cut_value = 3.0
    ievt = 0
    print("number of events:", enumerate(LC_nhits))
    for ievt, evt_nhits in enumerate(LC_nhits):
        #print("number of tracks in this event:", len(LC_track_pt[ievt]))
        for itrk, nhits in enumerate(evt_nhits):
            if mcp_mu_pt[ievt][itrk] < 1.0:
                
                print("-------- TRACK DATA --------------")
                print("event:", ievt, "track:", itrk)
                print("LC_fraction_muon_matched_hits:", LC_fraction_muon_matched_hits[ievt][itrk])
                print("LC_track_pt:", LC_track_pt[ievt][itrk])
                print("LC_track_chi2_red:", LC_track_chi2_red[ievt][itrk])
                print("LC_pt_match:", LC_pt_match[ievt][itrk])
                print("mcp_mu_pt:", mcp_mu_pt[ievt][itrk])
                print("LC_track_theta:", LC_track_theta[ievt][itrk])
                print("LC_d0:", LC_d0[ievt][itrk])
                print("LC z0 residual:", LC_z0residual[ievt][itrk])
                print("nhits:", nhits)
        ievt += 1
    print("number of total events:", ievt)
    
    #track_clean = (LC_track_pt>=pt_cut_value) * (LC_nhits>nhits_cut_value) * (LC_d0<= d0_cut_value) * (LC_track_chi2_red < chi2_red_cut_value) 
    #track_clean = (LC_track_pt>=pt_cut_value) * (LC_nhits>nhits_cut_value) * (LC_track_chi2_red < chi2_red_cut_value)
    #track_clean = (LC_track_pt>=pt_cut_value) * (LC_track_chi2_red < chi2_red_cut_value)
    track_clean = 1
    #track_clean = (LC_nhits>nhits_cut_value)
    #track_clean = (LC_track_pt>=pt_cut_value) * (LC_track_chi2_red < chi2_red_cut_value)
    # Define the eta transition region from barrel to endcap
    # Transition region now defined to be 0.6 - beyond this tracks will at least cross through part of VXD endcap
    transition_region = 0.6
    acceptance_edge = 2.44

    # Separate the data into barrel and endcap
    # track_barrel = (np.abs(LC_eta_match)<transition_region)
    
    # now always apply pT acceptance cut requiring that muon pT > 500 MeV, consistent with reco level cut
    truth_barrel = (np.abs(mcp_mu_eta)<transition_region) * (mcp_mu_pt > 0.5) * (np.abs(mcp_mu_eta)<acceptance_edge)
    # track_endcap = (np.abs(LC_eta_match)>=transition_region)
    truth_endcap = (np.abs(mcp_mu_eta)>=transition_region) * (mcp_mu_pt > 0.5) * (np.abs(mcp_mu_eta)<acceptance_edge)
    print("mcp_mu_pt:", mcp_mu_pt)
    #frac = ak.fill_none(LC_fraction_muon_matched_hits, 0.0)
    print(LC_fraction_muon_matched_hits)
    pt_matched_hits_acceptance = (LC_fraction_muon_matched_hits > 0.5)
    pt_acceptance = (mcp_mu_pt > 0.5)
    print("pt_acceptance:", pt_acceptance)

    # FIXME these should be configurable based on the data read in, so they don't need to be manually commented out each time
    #muon_gun_label = [r'Muon particle gun, uniform', r'in $p_T \in (0,1\text{ TeV})$ and $\theta$;']
    if upto10GeV:
        muon_gun_label = [r'Muon particle gun, uniform', r'in $p_T \in (0.5,10\text{ GeV})$ and $\theta$;']
    else: 
        muon_gun_label = [r'Muon particle gun, uniform', r'in $p_T \in (0.5,5000\text{ GeV})$ and $\theta$;']

    # Eta bin edges equivalent to the hardcoded theta bins — used in both efficiency and resolution plots.
    _theta_std_bins   = np.array([10., 40., 70., 90., 110., 140., 170.]) * np.pi / 180.
    _theta_pt_bins    = np.array([10., 30., 50., 60., 70., 80., 90.,
                                  100., 110., 120., 130., 150., 170.]) * np.pi / 180.
    eta_standard_bins = np.sort(-np.log(np.tan(_theta_std_bins / 2.0)))
    eta_pt_bins       = np.sort(-np.log(np.tan(_theta_pt_bins  / 2.0)))
    xlim_eta          = (float(eta_standard_bins[0]), float(eta_standard_bins[-1]))

    if(do_efficiency):
        # Binned in eta (same physics bins as before, just expressed in eta)
        bottom_label = muon_gun_label
        print('Computing reconstruction efficiency as a function of eta.')

        custom_bins = eta_standard_bins if not degrees else None
        xlim = xlim_eta if not degrees else None

        # NOTE: a bit messy, trying to make an "easy" toggle for whether or not before/after cleaning is shown
        mask_pairs = [
            (pt_matched_hits_acceptance, pt_acceptance), # now always apply acceptance cut requiring muon pT > 500 MeV
            (pt_matched_hits_acceptance, pt_acceptance) # no longer apply track cleaning, now done in reconstruction, now always apply acceptance cut requiring muon pT > 500 MeV
        ]
        labels=["Before Cleaning", "After Cleaning"]

        if(efficiency_opts == 1):
            mask_pairs = [mask_pairs[0]]
            labels = None # if before cleaning, don't write anything
            labels=["Loose Tracks","UltraTight Tracks"]
        elif(efficiency_opts == 2):
            mask_pairs = [mask_pairs[1]]
            labels = None # NOTE: also not writing legend if only after cleaning, might want to add misc text for this?
            bottom_label += ['After cleaning']
        results, min_value, max_value = calculate_efficiencies(
            LC_SiTrack_eta if not degrees else LC_SiTrack_theta,
            LC_track_eta   if not degrees else LC_track_theta,
            mcp_mu_eta     if not degrees else mcp_mu_theta,
            mask_pairs=mask_pairs,
            custom_bins=custom_bins
        )
        xlabel = r"True Muon $\eta$"
        if(degrees):
            xlabel = r"Muon $\theta$ [$\degree$]"
        plotter.plot_efficiencies(results, min_value, max_value,
                        xlabel=xlabel,
                        labels=labels,
                        savename='eff_vs_theta_{}'.format(suffix),
                        xlim=xlim,
                        bottom_label=bottom_label,
                        theta_axis=(not degrees)
                        )

        custom_bins_eta = np.linspace(-2.6,2.6,53) 
        results, min_value, max_value = calculate_efficiencies(
            LC_SiTrack_eta,
            LC_track_eta,
            mcp_mu_eta,
            mask_pairs=mask_pairs,
            num_bins=52,
            custom_bins=custom_bins_eta
        )
        xlabel = r"Muon $\eta$ "

        # FIXME xlim and min_value, max_value don't do anything unless 
        xlim = (-2.6,2.6)
        plotter.plot_efficiencies(results, -2.6, 2.6,
                        xlabel=xlabel,
                        labels=labels,
                        savename='eff_vs_eta_fine_{}'.format(suffix),
                        xlim=xlim,
                        bottom_label=bottom_label
                        )

        custom_bins_eta = np.linspace(-2.6,2.6,27) 
        results, min_value, max_value = calculate_efficiencies(
            LC_SiTrack_eta,
            LC_track_eta,
            mcp_mu_eta,
            mask_pairs=mask_pairs,
            num_bins=26,
            custom_bins=custom_bins_eta
        )
        xlabel = r"Muon $\eta$ "

        # FIXME xlim and min_value, max_value don't do anything unless 
        xlim = (-2.6,2.6)
        plotter.plot_efficiencies(results, -2.6, 2.6,
                        xlabel=xlabel,
                        labels=labels,
                        savename='eff_vs_eta_less_fine_{}'.format(suffix),
                        xlim=xlim,
                        bottom_label=bottom_label
                        )

        # Binned in pT, split into barrel and endcap regions
        print('Computing reconstruction efficiency as a function of pT, for barrel region.')
        #custom_bins = [1,2,5,10,20,50,100,200,500,1000,2000,5000]
        #custom_bins = [1,2,5,10,20,50,100,200,500,1000]
        if upto10GeV: 
            custom_bins = [0.5,1,2,3,4,5,6,7,8,9,10]
        else: 
            custom_bins = [1,2,5,10,20,50,100,200,500,1000,2000,5000]
        mask_pairs=[
            (pt_matched_hits_acceptance,truth_barrel), # for muons in the barrel region
            #(combine_masks([track_clean,truth_barrel]),truth_barrel), # for muons in the barrel region, with track cleaning applied to tracks
            ]
        if(efficiency_opts == 1):
            mask_pairs = [mask_pairs[0]]
        elif(efficiency_opts == 2):
            mask_pairs = [mask_pairs[1]]

        results, min_value, max_value = calculate_efficiencies(
            LC_SiTrack_pt_match,
            LC_pt_match,
            mcp_mu_pt,
            mask_pairs=mask_pairs,
            custom_bins=custom_bins
        )

        mask_pairs_endcap=[
            (pt_matched_hits_acceptance,truth_endcap), # for muons in the barrel region
            #(combine_masks([track_clean,truth_endcap]),truth_endcap), # for muons in the barrel region, with track cleaning applied to tracks
            ]

        misctext = r'$|\eta|<0.6$'
        if(degrees):
            misctext = r'$40^{\circ}<\theta<140^{\circ}$'

        plotter.plot_efficiencies(results, min_value, #10,
                        max_value,
                        xlabel="Muon $p_T$ [GeV]",
                        labels=labels,
                        misctext=misctext,
                        savename='eff_vs_pt_barrel_{}'.format(suffix),
                        bottom_label=bottom_label,
                        upTo10GeVBool=upto10GeV
                        )

        print('Computing reconstruction efficiency as a function of pT, for endcap region.')

        results, min_value, max_value = calculate_efficiencies(
            LC_SiTrack_pt_match,
            LC_pt_match,
            mcp_mu_pt,
            mask_pairs=mask_pairs_endcap,
            custom_bins=custom_bins
        )

        misctext = r'$0.6 < |\eta| < 2.44$'
        if(degrees):
            misctext = r'$\theta<40^{\circ}$ or $\theta>140^{\circ}$'
        plotter.plot_efficiencies(results, min_value, 
                        max_value,
                        #10, #max_value,
                        xlabel="Muon $p_T$ [GeV]",
                        labels=labels,
                        misctext=misctext,
                        savename='eff_vs_pt_endcap_{}'.format(suffix),
                        bottom_label=bottom_label,
                        upTo10GeVBool=upto10GeV
                        )

    ptmask1 = LC_pt_match<=50
    ptmask2 = (LC_pt_match>50) & (LC_pt_match<=250)
    ptmask3 = (LC_pt_match>250) & (LC_pt_match<=1000)
    ptmask4 = LC_pt_match>=1000
    # Split the data into 4 pT ranges (and fold the data for the high pT range)
    if(upto1TeV):
        theta_all    = [LC_theta_match[ptmask1], LC_theta_match[ptmask2], LC_theta_match[ptmask3]]
        d0_all       = [LC_d0[ptmask1]         , LC_d0[ptmask2]         , LC_d0[ptmask3]         ]
        z0_residual_all       = [LC_z0residual[ptmask1]         , LC_z0residual[ptmask2]         , LC_z0residual[ptmask3]         ]
        pt_res_all   = [LC_pt_res[ptmask1]     , LC_pt_res[ptmask2]     , LC_pt_res[ptmask3]     ]
        pt_track_all = [LC_track_pt[ptmask1]   , LC_track_pt[ptmask2]   , LC_track_pt[ptmask3]   ]
        pt_match_all = [LC_pt_match[ptmask1]   , LC_pt_match[ptmask2]   , LC_pt_match[ptmask3]   ]
        nhits_all    = [LC_nhits[ptmask1]      , LC_nhits[ptmask2]      , LC_nhits[ptmask3]      ]
    else:
        theta_all    = [LC_theta_match[ptmask1], LC_theta_match[ptmask2], LC_theta_match[ptmask3], fold_data(ak.flatten(LC_theta_match[ptmask4]), LC_theta_match, LC_pt_match)]
        d0_all       = [LC_d0[ptmask1]         , LC_d0[ptmask2]         , LC_d0[ptmask3]         , fold_data(ak.flatten(LC_d0[ptmask4])         , LC_theta_match, LC_pt_match)]
        z0_residual_all       = [LC_z0residual[ptmask1]         , LC_z0residual[ptmask2]         , LC_z0residual[ptmask3]         , fold_data(ak.flatten(LC_z0residual[ptmask4])         , LC_theta_match, LC_pt_match)]
        pt_res_all   = [LC_pt_res[ptmask1]     , LC_pt_res[ptmask2]     , LC_pt_res[ptmask3]     , fold_data(ak.flatten(LC_pt_res[ptmask4])     , LC_theta_match, LC_pt_match)]
        pt_track_all = [LC_track_pt[ptmask1]   , LC_track_pt[ptmask2]   , LC_track_pt[ptmask3]   , fold_data(ak.flatten(LC_track_pt[ptmask4])   , LC_theta_match, LC_pt_match)]
        pt_match_all = [LC_pt_match[ptmask1]   , LC_pt_match[ptmask2]   , LC_pt_match[ptmask3]   , fold_data(ak.flatten(LC_pt_match[ptmask4])   , LC_theta_match, LC_pt_match)]
        nhits_all    = [LC_nhits[ptmask1]      , LC_nhits[ptmask2]      , LC_nhits[ptmask3]      , fold_data(ak.flatten(LC_nhits[ptmask4])      , LC_theta_match, LC_pt_match)]
    
    print("len(z0_residual_all):", len(z0_residual_all))
    print('Checking some min/max stuff.')
    #print('min/max of pt_match_all[0] (  0 -   50 GeV) :',np.min(pt_match_all[0]), np.max(pt_match_all[0]),)
    #print('min/max of pt_match_all[1] ( 50 -  250 GeV) :',np.min(pt_match_all[1]), np.max(pt_match_all[1]),)
    #print('min/max of pt_match_all[2] (250 - 1000 GeV) :',np.min(pt_match_all[2]), np.max(pt_match_all[2]),)

    # print('Sum of intersection of ptmask2 and ptmask3: {}'.format(np.sum(np.logical_and(ptmask2,ptmask3))))

    # print('Printing entries of np.ravel(LC_pt_match), and ptmask2.')
    # for i,entry in enumerate(np.ravel(LC_pt_match)):
    #     print('[{}], {} -> {}'.format(i,entry,ptmask2[i]),end='')
    #     if(ptmask2[i]):
    #         print(' <-------------------')
    #     else:
    #         print()

    # (Re-)apply cleaning, will be used for all resolution plots.
    theta_all_masked = []
    d0_all_masked = []
    z0_residuals_all_masked = []
    pt_res_all_masked = []
    pt_res2_all_masked = []
    pt_match_all_masked = []
    for i in range(len(theta_all)):
        # Create a boolean mask for the condition
        if(degrees):
            theta_cut = (0 <= theta_all[i]) & (theta_all[i] < 180)
        else:
            theta_cut = (0 <= theta_all[i]) & (theta_all[i] < np.pi)
        pt_res_cut = np.abs(pt_res_all[i]) > 0
        pt_cut = pt_track_all[i] > pt_cut_value
        d0_cut = np.abs(d0_all[i]) <= d0_cut_value
        nhits_cut = nhits_all[i] > nhits_cut_value
        mask = pt_res_cut # now only apply pt res cut, since cuts applied in track selector in reconstruction 

        # Apply the mask to filter the arrays
        x_masked = theta_all[i][mask]
        y_masked = d0_all[i][mask]
        z0_residual_masked = z0_residual_all[i][mask]
        z_masked = pt_res_all[i][mask]
        w_masked = pt_match_all[i][mask]
        t_masked = pt_res_all[i][mask]/pt_match_all[i][mask]

        theta_all_masked.append(x_masked)
        d0_all_masked.append(y_masked)
        z0_residuals_all_masked.append(z0_residual_masked)
        pt_res_all_masked.append(z_masked)
        pt_res2_all_masked.append(t_masked)
        pt_match_all_masked.append(w_masked)

    # Convert theta (radians) → eta for resolution-vs-eta plots.
    # theta_all_masked entries are jagged awkward arrays; flatten the same way calc_utils does.
    if degrees:
        eta_all_masked = [-np.log(np.tan(np.radians(ak.to_numpy(np.ravel(x)).astype(float)) / 2.0)) for x in theta_all_masked]
    else:
        eta_all_masked = [-np.log(np.tan(ak.to_numpy(np.ravel(x)).astype(float) / 2.0)) for x in theta_all_masked]
    # eta_standard_bins, eta_pt_bins, xlim_eta already computed above (before efficiency section)

    if(upto1TeV):
        numpoints = 4 # this is number of bins we want 
    else: 
        numpoints = 4 # this is number of bins we want 
    
    array1 = np.linspace(-1,1,300)    # pt_0_50_bins
    array2 = np.linspace(-0.5,0.5,300)  # pt_50_250_bins
    array3 = np.linspace(-0.1,0.1,300) # pt_250_1000_bins
    array4 = np.linspace(-0.1,0.1,300)  # pt_1000_5000_bins

    if(upto1TeV):
        d0_bins = [array1, array2, array3] # for 0, 1 TeV
    else: 
        d0_bins = [array1, array2, array3, array4] # for 0, 5 TeV

    array1 = np.linspace(-5,5,300)    # pt_0_50_bins
    array2 = np.linspace(-5,5,300)  # pt_50_250_bins
    array3 = np.linspace(-5,5,300) # pt_250_1000_bins
    array4 = np.linspace(-5,5,300)  # pt_1000_5000_bins

    if(upto1TeV):
        z0_bins = [array1, array2, array3] # for 0, 1 TeV
    else: 
        z0_bins = [array1, array2, array3, array4] # for 0, 5 TeV
    

    array1 = np.linspace(-0.1,0.1,300)    # pt_0_50_bins
    array2 = np.linspace(-0.1,0.1,300)  # pt_50_250_bins
    array3 = np.linspace(-0.1,0.1,300) # pt_250_1000_bins
    array4 = np.linspace(-0.1,0.1,300)  # pt_1000_5000_bins
    if(upto1TeV):
        pt_bins = [array1, array2, array3] # for 0, 1 TeV
    else: 
        pt_bins = [array1, array2, array3, array4] # for 0, 5 TeV
   

    array1 = np.linspace(-0.003,0.003,300)    # pt_0_50_bins
    array2 = np.linspace(-0.003,0.003,300)  # pt_50_250_bins
    array3 = np.linspace(-0.002,0.002,300) # pt_250_1000_bins
    array4 = np.linspace(-0.002,0.002,300)  # pt_1000_5000_bins
    
    if(upto1TeV):
        pt_2_bins = [array1, array2, array3] # for 0, 1 TeV
    else: 
        pt_2_bins = [array1, array2, array3, array4] # for 0, 5 TeV
    # d0_ylim = (0,0.01)
    # pt_ylim = (0,0.01)

    # Code for plotting resolutions
    if(upto10GeV):
        labels = [r'$p_T$ $\in$ 0.5-10 GeV']
    else:
        labels = [r'$p_T$ $\in$ 0.5-50 GeV', r'$p_T$ $\in$ 50-250 GeV', r'$p_T$ $\in$ 250-1000 GeV', r'$p_T$ $\in$ 1000-5000 GeV']
    
    misctext = muon_gun_label + ['UltraTight Tracks']

    # d0 resolution versus eta
    processed_data = process_data(
        datax=eta_all_masked if not degrees else theta_all_masked,
        datay=d0_all_masked,
        numbins=numpoints,
        bins=d0_bins,
        x_bins=eta_standard_bins if not degrees else None,
        theta=degrees,
        degrees=degrees
    )
    xlabel = r'True Muon $\eta$'
    xlim = xlim_eta
    if degrees:
        xlabel = r'Muon $\theta$ [$\degree$]'
        xlim = (0, 180)
    plotter.plot_processed_data(
        processed_results=processed_data,
        labels=labels,
        xlabel=xlabel,
        ylabel=r'$\sigma(d_0)$ [mm]',
        ylim=(0.001,1.),
        xlim=xlim,
        log=True,
        savename="res_d0_vs_theta_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62),
        theta_axis=(not degrees)
    )

    # z0 resolution versus eta
    print("len(z0_residuals_all_masked):", len(z0_residuals_all_masked))
    processed_data = process_data(
        datax=eta_all_masked if not degrees else theta_all_masked,
        datay=z0_residuals_all_masked,
        numbins=numpoints,
        bins=d0_bins,
        x_bins=eta_standard_bins if not degrees else None,
        theta=degrees,
        degrees=degrees
    )
    plotter.plot_processed_data(
        processed_results=processed_data,
        labels=labels,
        xlabel=xlabel,
        ylabel=r'$\sigma(z_0)$ [mm]',
        ylim=(0.001,1.),
        xlim=xlim,
        log=True,
        savename="res_z0_vs_theta_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62),
        theta_axis=(not degrees)
    )

    # pt resolution versus eta
    processed_data = process_data(
        datax=eta_all_masked if not degrees else theta_all_masked,
        datay=pt_res_all_masked,
        numbins=numpoints,
        bins=pt_bins,
        x_bins=eta_pt_bins if not degrees else None,
        theta=degrees,
        degrees=degrees,
        pTvsTheta=degrees  # pTvsTheta bins only apply in degrees/legacy mode
    )
    xlabel = r'True Muon $\eta$'
    xlim = (float(eta_pt_bins[0]), float(eta_pt_bins[-1]))
    if degrees:
        xlabel = r'Muon $\theta$ [$\degree$]'
        xlim = (0, 180)
    plotter.plot_processed_data(
        processed_results=processed_data,
        labels=labels,
        xlabel=xlabel,
        ylabel= r'$\sigma(p_T)/p_T$',
        ylim=(0.001,10.0),
        xlim=xlim,
        log=True,
        savename="res_pt_vs_theta_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62),
        theta_axis=(not degrees)
    )

    # pt resolution versus pt
    print("len of pt_res_all_masked:", len(pt_res_all_masked))
    if(upto10GeV):
        use2BinsFinalDataSlice = False
    else:
        use2BinsFinalDataSlice = True
    
    print("processing pt vs pt")
    processed_data = process_data(
        datax=pt_match_all_masked,
        datay=pt_res_all_masked,
        numbins=numpoints,
        bins=pt_bins,
        debug=True,
        use2BinsFinalpTSlice=use2BinsFinalDataSlice
    )

    plotter.plot_processed_data(
        processed_results=processed_data,
        labels=labels,
        xlabel=r'Truth muon $p_T [GeV]$',
        ylabel= r'$\sigma(p_T)/p_T$',
        xlog=True,
        log=True,
        # title=r'Single $\mu^{\pm}$ no BIB @ 10TeV',
        ylim=(0,1.0),
        savename="res_pt_vs_pt_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62)
    )
    # pt^2 resolution versus eta
    processed_data = process_data(
        datax=eta_all_masked if not degrees else theta_all_masked,
        datay=pt_res2_all_masked,
        numbins=numpoints,
        bins=pt_2_bins,
        x_bins=eta_standard_bins if not degrees else None,
        theta=degrees,
        degrees=degrees,
        #debug=True,
        debug_directory=outdir + '/debug',
        debug_name='debug_pt2_v_theta',
        debug_xlabel='#eta',
        debug_labels=labels
    )
    xlabel = r'True Muon $\eta$'
    xlim = xlim_eta
    if degrees:
        xlabel = r'Muon $\theta$ [$\degree$]'
        xlim = (0, 180)
    plotter.plot_processed_data(
        processed_results=processed_data,
        labels=labels,
        xlabel=xlabel,
        ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
        ylim=(0.00001,0.6),
        log=True,
        xlim=xlim,
        savename="res_pt2_vs_theta_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62),
        theta_axis=(not degrees)
    )

    processed_data = process_data(
        datax=pt_match_all_masked,
        datay=pt_res2_all_masked,
        numbins=3,
        bins=pt_2_bins,
        # debug=True,
        # debug_directory=outdir,
        # debug_name='debug_res_pt2_vs_pt',
        # debug_labels=['p_{T} #in 0-50 GeV', 'p_{T} #in 50-250 GeV', 'p_{T} #in 250-1000 GeV', 'p_{T} #in 1000-5000 GeV'],
        # debug_xlabel='p_{T}[GeV]'
    )
    plotter.plot_processed_data(
        processed_results=processed_data,
        # labels=[r'$p_T$ $\in$ 0-50 GeV', r'$p_T$ $\in$ 50-250 GeV', r'$p_T$ $\in$ 250-1000 GeV', r'$p_T$ $\in$ 1000-5000 GeV'],
        xlabel=r'Truth muon $p_T [GeV]$',
        ylabel= r'$\sigma(p_T)/p_T^2$ $[GeV^{-1}]$',
        ylim=(0,0.01),
        xlog=True,
        log=True,
        savename="res_pt2_vs_pt_{}".format(suffix),
        misctext=misctext,
        legend_loc=(0.52,0.62)
    )


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
    #     data_flatx_alt = np.array(np.ravel(datax)).T
    #     data_flaty_alt = np.array(np.ravel(datay)).T
    #     for i in range(len(x_bins)-1):
    #         # Select data points within the specified x slice
    #         slice_data = data_flaty_alt[(data_flatx_alt >= x_bins[i]) & (data_flatx_alt < x_bins[i + 1])]
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
    # data_loader = bib_250_1000

    # binning = [min(data_loader['mcp_mu_pt']), max(data_loader['mcp_mu_pt'])]
    # # which_theta = np.degrees(fold_data(bib_250_1000['LC_track_theta'], bib = True)) # THIS HAS TO BE THE SAME AS data_loader, also check for degrees vs radians!!!
    # which_theta = np.degrees((data_loader['LC_track_theta']))

    # for i in range(nbins):
    #     print(r'Theta:', theta_bins[i], r'<= theta <', theta_bins[i+1])
    #     theta_bin = ((theta_bins[i] <= which_theta) & (which_theta < theta_bins[i+1]) & [np.abs(LC_pt_resolution[0]) < 1 for LC_pt_resolution in (data_loader['LC_pt_res'])])
    #     # print(data_loader['LC_pt_res'][theta_cut])
    #     count = 0
    #     for bin in theta_bin:
    #         if bin[0] == True:
    #             count +=1
    #     # print("# of data points (total, pt > 1000):", count, len(ak.flatten(fold_data(data_loader['LC_pt_match'], bib = True)[theta_bin[fold_data(data_loader['LC_pt_match']>1000, bib = True)]]))) # Not sure how 'count' works but it does so don't worry
    #     # print("# of data points (total):", count, len(ak.flatten((data_loader['LC_pt_match'])[theta_bin[(data_loader['LC_pt_match']>1000)]]))) # Not sure how 'count' works but it does so don't worry

    #     # /pT resolution
    #     # plotrms_slice(data_loader['LC_pt_match'][theta_bin], (data_loader['LC_pt_res'])[theta_bin], x_bins = x_bins_250_1000, bins=pt_bins, title=pt_title, rv = False, sigma5 = False)

    #     # /pT^2 resolution
    #     plotrms_slice((data_loader['LC_pt_match'])[theta_bin], ((data_loader['LC_pt_res']/data_loader['LC_pt_match']))[theta_bin], x_bins = binning, bins=pt_bins, xlim = None, title=pt_title+r'/$p_T^2$', rv = False, sigma5 = False)

    #     # d0 resolution
    #     # plotrms_slice(data_loader['LC_pt_match'][theta_bin], data_loader['LC_d0'][theta_bin], x_bins = [1000,5000], bins=pt_bins, xlim = None, title=d0_title, rv = False, sigma5 = False)
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

    # 2D heatmaps: fraction of muon-matched hits (SiTracks) vs. pT resolution and vs. q/p residual
    heatmap_misctext = muon_gun_label + ['Loose Tracks']  # reuse "Selected Tracks" + muon gun label

    plotter.plot_2d_histogram(
        datax=LC_fraction_muon_matched_hits,
        datay=LC_pt_res_SiTracks,
        bins=[25, 50],
        xlabel=r'Fraction of muon-matched hits',
        ylabel=r'$p_T$ resolution $(p_T^{\rm truth} - p_T^{\rm reco})/p_T^{\rm truth}$',
        xlim=(0., 1.),
        ylim=(-10., 10.),
        log_y=True,
        savename='heatmap_fraction_muon_hits_vs_pt_res_{}'.format(suffix),
        misctext=heatmap_misctext,
        log_color=True,
    )

    plotter.plot_2d_histogram(
        datax=LC_fraction_muon_matched_hits,
        datay=LC_qoverp_residual,
        bins=[25, 50],
        xlabel=r'Fraction of muon-matched hits',
        ylabel=r'$(q/p_{\rm truth} - q/p_{\rm reco})^2$ $[\mathrm{GeV}^{-2}]$',
        xlim=(0., 1.),
        square_y=True,
        savename='heatmap_fraction_muon_hits_vs_qoverp_res2_{}'.format(suffix),
        misctext=heatmap_misctext,
        log_color=True,
    )

if(__name__=='__main__'):
    main(sys.argv)