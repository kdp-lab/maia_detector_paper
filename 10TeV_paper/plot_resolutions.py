#
# A modified version f 10TeV_plotting_scripts.py, that just does
# the resolution plots -- for combining BIB and noBIB results.
# Author: Jan T. Offermann
#
import numpy as np
import ROOT as rt
import awkward as ak
import sys
import mplhep as hep
import argparse as ap
import uproot as ur
# hep.style.use(hep.style.ROOT) # For now ROOT defaults to CMS

from utils.data_utils import DataLoader
from utils.plot_utils import Plotter
from utils.calc_utils import calculate_efficiencies,combine_masks, process_data

def MakeDataSubsets(array, masks):
    d = {}
    for key in array.keys():
        d[key] = [array[key][mask[key]] for mask in masks]
    return d

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
            'LC_z0':'lc-matched_track_z0',
            'LC_chi2':'lc-matched_track_chi2',
            'LC_ndf':'lc-matched_track_ndf',
            'LC_pt_res':'lc-matched_track_ptres',
            'LC_track_pt':'lc-matched_track_pt',
            'LC_track_theta':'lc-matched_track_theta',
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

    hep.style.use(hep.style.ATLAS)
    rt.gROOT.SetBatch(True)
    rt.gStyle.SetOptStat(0)

    parser = ap.ArgumentParser()
    parser.add_argument('-i','--inputFile',type=str,required=True,nargs='+')
    parser.add_argument('-k','--keys',type=str,nargs='+',default=[])
    parser.add_argument('-o','--outputDirectory',type=str,default='output')
    parser.add_argument('-dl','--dataLabel',type=str,default='Simulation, BIB')
    parser.add_argument('-ll','--latticeLabel',type=str,default='EU24 Lattice')
    parser.add_argument('-s','--suffix',type=str,default='all')
    parser.add_argument('-degrees','--degrees',type=int,default=0)
    parser.add_argument('-fine_binning','--fine_binning',type=int,default=0)

    parser.add_argument('-d0','--d0',type=int,default=1)
    parser.add_argument('-pt','--pt',type=int,default=1)
    parser.add_argument('-pt2','--pt2',type=int,default=1)

    args = vars(parser.parse_args())
    infiles = args['inputFile']
    keys = args['keys']

    if(len(keys) < len(infiles)):
        print('Warning: Keys not understood.')
        keys = [x.split('/')[-1] for x in infiles]
    elif(len(keys) > len(infiles)):
        print('Warning: Truncating keys.')
        # keys = keys[:len(infiles)]

    outdir = args['outputDirectory']
    data_label = args['dataLabel']
    lattice_label = args['latticeLabel']
    suffix = args['suffix']
    degrees = args['degrees'] > 0

    do_d0  = args['d0' ] > 0
    do_pt  = args['pt' ] > 0
    do_pt2 = args['pt2'] > 0

    fine_binning = args['fine_binning'] > 0

    # Defining cleaning differently here, as an awkward array
    pt_cut_value = 1. # GeV
    d0_cut_value = 0.1 # mm?
    nhits_cut_value = 4

    infiles = {keys[i]:infiles[i] for i in range(len(infiles))}
    key_mapping = {key: CreateKeyMapping(val) for key,val in infiles.items()}

    # Load the data. For ROOT files, uproot will handle things so that memory usage is reasonable.
    data_loader = {key:DataLoader() for key in infiles.keys()}
    for key,loader in data_loader.items():
        loader.SetVerbose(True)
        loader.SetFilename(infiles[key])
        print('Loading data from {}.'.format(infiles[key]))
        loader.Load()


    plotter = Plotter()
    plotter.SetCOMTev(10)
    plotter.SetDataLabel(data_label)
    plotter.SetLatticeLabel(lattice_label)
    plotter.SetOutputDirectory(outdir)

    # Assign variables
    LC_pt_match = {key:loader[key_mapping[key][loader.GetMode()]['LC_pt_match']] for key,loader in data_loader.items()}
    LC_eta_match = {key:loader[key_mapping[key][loader.GetMode()]['LC_eta_match']] for key,loader in data_loader.items()}
    LC_theta_match = {key:2 * np.arctan(np.exp(-LC_eta_match[key])) for key in data_loader.keys()}
    LC_track_pt = {key:loader[key_mapping[key][loader.GetMode()]['LC_track_pt']] for key,loader in data_loader.items()}
    LC_track_theta = {key:loader[key_mapping[key][loader.GetMode()]['LC_track_theta']] for key,loader in data_loader.items()}
    mcp_mu_pt = {key:loader[key_mapping[key][loader.GetMode()]['mcp_mu_pt']] for key,loader in data_loader.items()}
    mcp_mu_eta = {key:loader[key_mapping[key][loader.GetMode()]['mcp_mu_eta']] for key,loader in data_loader.items()}
    mcp_mu_theta = {key:2 * np.arctan(np.exp(-mcp_mu_eta[key])) for key in data_loader.keys()}
    LC_nhits = {key:loader[key_mapping[key][loader.GetMode()]['LC_nhits']] for key,loader in data_loader.items()}
    LC_d0 = {key:loader[key_mapping[key][loader.GetMode()]['LC_d0']] for key,loader in data_loader.items()}
    LC_pt_res = {key:loader[key_mapping[key][loader.GetMode()]['LC_pt_res']] for key,loader in data_loader.items()}

    LC_theta_match_mirrored = {key:2 * np.arctan(np.exp(-np.abs(LC_eta_match[key]))) for key in data_loader.keys()}

    if(degrees):
        LC_theta_match = {key:np.degrees(val) for key,val in LC_theta_match.items()}
        LC_track_theta = {key:np.degrees(val) for key,val in LC_track_theta.items()}
        mcp_mu_theta = {key:np.degrees(val) for key,val in mcp_mu_theta.items()}

    muon_gun_label = [r'Muon particle gun', r'uniform in $p_{T}$ and $\theta$;']

    ptmask1 = {key:LC_pt_match[key]<=50 for key in LC_pt_match.keys()}
    ptmask2 = {key:(LC_pt_match[key]>50) & (LC_pt_match[key]<=250) for key in LC_pt_match.keys()}
    ptmask3 = {key:(LC_pt_match[key]>250) & (LC_pt_match[key]<=1000) for key in LC_pt_match.keys()}
    # ptmask4 = {key:LC_pt_match>=1000 for key in LC_pt_match.keys()}
    # Split the data into 4 pT ranges (and fold the data for the high pT range)
    pt_masks = [ptmask1, ptmask2, ptmask3]

    theta_all = MakeDataSubsets(LC_theta_match,pt_masks)
    d0_all = MakeDataSubsets(LC_d0,pt_masks)
    pt_res_all = MakeDataSubsets(LC_pt_res,pt_masks)
    pt_track_all = MakeDataSubsets(LC_track_pt,pt_masks)
    pt_match_all = MakeDataSubsets(LC_pt_match,pt_masks)
    nhits_all = MakeDataSubsets(LC_nhits,pt_masks)

    theta_mirrored_all = MakeDataSubsets(LC_theta_match_mirrored,pt_masks)

    # (Re-)apply cleaning, will be used for all resolution plots.
    theta_all_masked = {}
    d0_all_masked = {}
    pt_res_all_masked = {}
    pt_res2_all_masked = {}
    pt_match_all_masked = {}

    theta_mirrored_all_masked = {}

    for key in theta_all.keys():

        theta_all_masked[key] = []
        d0_all_masked[key] = []
        pt_res_all_masked[key] = []
        pt_res2_all_masked[key] = []
        pt_match_all_masked[key] = []

        theta_mirrored_all_masked[key] = []

        for i in range(len(theta_all[key])):
            # Create a boolean mask for the condition
            if(degrees):
                theta_cut = (0 <= theta_all[key][i]) & (theta_all[key][i] < 180)
            else:
                theta_cut = (0 <= theta_all[key][i]) & (theta_all[key][i] < np.pi)
            pt_res_cut = np.abs(pt_res_all[key][i]) > 0
            pt_cut = pt_track_all[key][i] > pt_cut_value
            d0_cut = np.abs(d0_all[key][i]) <= d0_cut_value
            nhits_cut = nhits_all[key][i] > nhits_cut_value
            mask = theta_cut & pt_res_cut & pt_cut & d0_cut & nhits_cut

            theta_all_masked[key].append(theta_all[key][i][mask])
            d0_all_masked[key].append(d0_all[key][i][mask])
            pt_res_all_masked[key].append(pt_res_all[key][i][mask])
            pt_res2_all_masked[key].append(pt_res_all[key][i][mask]/pt_match_all[key][i][mask])
            pt_match_all_masked[key].append(pt_match_all[key][i][mask])

            theta_mirrored_all_masked[key].append(theta_mirrored_all[key][i][mask])

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

    array1 = np.linspace(-0.003,0.003,100)    # pt_0_50_bins #300 bins originally (for all)
    array2 = np.linspace(-0.003,0.003,150)  # pt_50_250_bins
    array3 = np.linspace(-0.002,0.002,150) # pt_250_1000_bins
    array4 = np.linspace(-0.002,0.002,150)  # pt_1000_5000_bins
    pt_2_bins = [array1, array2, array3, array4]
    # d0_ylim = (0,0.01)
    # pt_ylim = (0,0.01)

    if(fine_binning):
        x_bins_theta = np.array([30.*np.pi/180.,
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
    else:
        x_bins_theta = np.array([30.*np.pi/180.,
                            50.*np.pi/180.,
                            70.*np.pi/180.,
                            90.*np.pi/180.,
                            110.*np.pi/180.,
                            130.*np.pi/180.,
                            150.*np.pi/180.
                            ]
                        )

    x_bins_theta_mirrored = np.array([30.*np.pi/180.,
                            40.*np.pi/180.,
                            50.*np.pi/180.,
                            60.*np.pi/180.,
                            70.*np.pi/180.,
                            90.*np.pi/180.,
                            ]
                        )

    misctext = muon_gun_label

    labels_template = [r'muon $p_{T}$ $\in$ 0.5-50 GeV, {}', r'muon $p_{T}$ $\in$ 50-250 GeV, {}', r'muon $p_{T}$ $\in$ 250-1000 GeV, {}']
    labels = []
    if(len(keys) == 1):
        labels = [x.replace(', {}','') for x in labels_template]
        if('with bib' in keys[0].lower()):
            misctext = muon_gun_label + ['After cleaning'] # TODO: fragile code

    else:
        for key in data_loader.keys():
            labels += [x.format(key) for x in labels_template]

    # Plot some kinematic distributions
    # TODO: This code is a bit fragile.

    plot_data = []
    for key in theta_all_masked.keys():
        plot_data += theta_all_masked[key]

    plotter.plot_distributions(
        plot_data,
        labels,
        bins=x_bins_theta,
        x_label=r'Muon $\theta$ [rad]',
        y_label='',
        x_range=(0.5,np.pi-0.5),
        y_scale='linear',
        savename='theta_dist_{}'.format(suffix)
    )

    # Code for plotting resolutions
    if(do_d0):
        # d0 resolution versus theta
        processed_data = {
            key:process_data(
                datax=theta_all_masked[key],
                datay=d0_all_masked[key],
                numbins=numpoints,
                bins=d0_bins,
                theta=True,
                degrees=degrees
            ) for key in theta_all_masked.keys()
        }

        # now make a list of processed_data
        # turn things into lists for plot_efficiencies()
        processed_data_list = []
        for key in data_loader.keys():
            processed_data_list += processed_data[key]

        xlabel=r'Truth muon $\theta$ [rad]'
        xlim=(0.5,np.pi-0.5)
        if(degrees):
            xlabel=r'Truth muon $\theta [\degree]$'
            xlim=(0,180)

        plotter.plot_processed_data(
            processed_results=processed_data_list,
            labels=labels,
            xlabel=xlabel,
            ylabel=r'$\sigma(d_0)$ [mm]',
            ylim=(0.001,1.),
            xlim=xlim,
            log=True,
            savename="res_d0_vs_theta_{}".format(suffix),
            misctext=misctext
        )

    # pt resolution versus theta
    if(do_pt):
        processed_data = {
            key:process_data(
                datax=theta_all_masked[key],
                datay=pt_res_all_masked[key],
                numbins=numpoints,
                bins=pt_bins,
                theta=True,
                degrees=degrees
            ) for key in theta_all_masked.keys()
        }
        # now make a list of processed_data
        # turn things into lists for plot_efficiencies()
        processed_data_list = []
        for key in data_loader.keys():
            processed_data_list += processed_data[key]

        xlabel=r'Truth muon $\theta$ [rad]'
        xlim=(0.5,np.pi-0.5)
        if(degrees):
            xlabel=r'Truth muon $\theta [\degree]$'
            xlim=(0,180)
        plotter.plot_processed_data(
            processed_results=processed_data_list,
            labels=labels,
            xlabel=xlabel,
            ylabel= r'$\sigma(p_{T})/p_{T}$',
            ylim=(0.001,1.0),
            xlim=xlim,
            log=True,
            savename="res_pt_vs_theta_{}".format(suffix),
            misctext=misctext
        )

        # pt resolution versus pt
        processed_data = {
            key:process_data(
                datax=pt_match_all_masked[key],
                datay=pt_res_all_masked[key],
                numbins=3,
                bins=pt_bins,
                theta=True,
                degrees=degrees
            ) for key in theta_all_masked.keys()
        }
        # now make a list of processed_data
        # turn things into lists for plot_efficiencies()
        processed_data_list = []
        for key in data_loader.keys():
            processed_data_list += processed_data[key]

        plotter.plot_processed_data(
            processed_results=processed_data_list,
            labels=labels,
            xlabel=r'Truth muon $p_{T} [GeV]$',
            ylabel= r'$\sigma(p_{T})/p_{T}$',
            # title=r'Single $\mu^{\pm}$ no BIB @ 10TeV',
            ylim=(0,0.1),
            savename="res_pt_vs_pt_{}".format(suffix),
            misctext=misctext
        )

    if(do_pt2):
        # pt^2 resolution versus theta
        x_bins = None
        if(not degrees): x_bins = x_bins_theta
        processed_data = {
            key:process_data(
                datax=theta_all_masked[key],
                datay=pt_res2_all_masked[key],
                numbins=numpoints,
                x_bins=x_bins,
                bins=pt_2_bins,
                theta=True,
                degrees=degrees,
                debug=True,
                debug_directory=outdir + '/debug',
                debug_name='debug_pt2_v_theta',
                debug_xlabel='#theta',
                debug_labels=labels
            ) for key in theta_all_masked.keys()
        }
        # now make a list of processed_data
        # turn things into lists for plot_efficiencies()
        processed_data_list = []
        for key in data_loader.keys():
            processed_data_list += processed_data[key]

        xlabel=r'Truth muon $\theta$ [rad]'
        xlim=(0.5,np.pi-0.5)
        if(degrees):
            xlabel=r'Truth muon $\theta [\degree]$'
            xlim=(0,180)
        plotter.plot_processed_data(
            processed_results=processed_data_list,
            labels=labels,
            xlabel=xlabel,
            ylabel= r'$\sigma(p_{T})/p_{T}^2$ $[GeV^{-1}]$',
            ylim=(0.00001,1.),
            log=True,
            xlim=xlim,
            legend_loc=(0.375,0.45),
            # label_block_y_up=0.85,
            savename="res_pt2_vs_theta_{}".format(suffix),
            misctext=misctext,
            # misctext_y_up=0.65
        )


        # pt^2 resolution versus theta_mirrored
        x_bins = x_bins_theta_mirrored
        processed_data = {
            key:process_data(
                datax=theta_mirrored_all_masked[key],
                datay=pt_res2_all_masked[key],
                numbins=numpoints,
                x_bins=x_bins,
                bins=pt_2_bins,
                theta=True,
                degrees=degrees,
                debug=True,
                debug_directory=outdir + '/debug',
                debug_name='debug_pt2_v_theta_mirrored',
                debug_xlabel='#theta_{mirrored}',
                debug_labels=labels
            ) for key in theta_all_masked.keys()
        }
        # now make a list of processed_data
        # turn things into lists for plot_efficiencies()
        processed_data_list = []
        for key in data_loader.keys():
            processed_data_list += processed_data[key]

        xlabel=r'Truth muon $\theta_\text{mirrored}$ [rad]'
        xlim=(0.5,np.pi/2.)
        if(degrees):
            xlabel=r'Truth muon $\theta_\text{mirrored}$ $[\degree]$'
            xlim=(0,180)
        plotter.plot_processed_data(
            processed_results=processed_data_list,
            labels=labels,
            xlabel=xlabel,
            ylabel= r'$\sigma(p_{T})/p_{T}^2$ $[GeV^{-1}]$',
            ylim=(0.00001,1.),
            log=True,
            xlim=xlim,
            legend_loc=(0.375,0.45),
            # label_block_y_up=0.85,
            savename="res_pt2_vs_theta_mirrored_{}".format(suffix),
            misctext=misctext,
            # misctext_y_up=0.65
        )


        # # pt^2 resolution versus pt
        # processed_data = {
        #     key:process_data(
        #         datax=pt_match_all_masked[key],
        #         datay=pt_res2_all_masked[key],
        #         numbins=3,
        #         bins=pt_2_bins,
        #         theta=True,
        #         degrees=degrees
        #     ) for key in theta_all_masked.keys()
        # }
        # # now make a list of processed_data
        # # turn things into lists for plot_efficiencies()
        # processed_data_list = []
        # for key in data_loader.keys():
        #     processed_data_list += processed_data[key]

        # plotter.plot_processed_data(
        #     processed_results=processed_data_list,
        #     # labels=[r'$p_{T}$ $\in$ 0-50 GeV', r'$p_{T}$ $\in$ 50-250 GeV', r'$p_{T}$ $\in$ 250-1000 GeV', r'$p_{T}$ $\in$ 1000-5000 GeV'],
        #     xlabel=r'Truth muon $p_{T} [GeV]$',
        #     ylabel= r'$\sigma(p_{T})/p_{T}^2$ $[GeV^{-1}]$',
        #     ylim=(0,0.01),
        #     xlog=True,
        #     log=True,
        #     savename="res_pt2_vs_pt_{}".format(suffix),
        #     misctext=misctext
        # )

    return


if(__name__=='__main__'):
    main(sys.argv)