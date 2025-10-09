#
# A modified version f 10TeV_plotting_scripts.py, that just does
# the efficiency plots -- for combining BIB and noBIB results.
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
from utils.calc_utils import calculate_efficiencies,combine_masks

# Function to fold the data over theta = 90 because detector is symmetric in theta.
# To be used in the case of low statistics (high pT and BIB data)
def fold_data(data, LC_theta_match, LC_pt_match):
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

    # allow skipping of the efficiency plots -- mostly useful for debugging, to speed things up
    parser.add_argument('-efficiencyOpts', '--efficiencyOpts', type=int, default=0, help='0 for before/after track cleaning, 1 for before only, 2 for after only')

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
    efficiency_opts = args['efficiencyOpts']

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

    if(degrees):
        LC_theta_match = {key:np.degrees(val) for key,val in LC_theta_match.items()}
        LC_track_theta = {key:np.degrees(val) for key,val in LC_track_theta.items()}
        mcp_mu_theta = {key:np.degrees(val) for key,val in mcp_mu_theta.items()}


    track_clean = {
        key:
        (LC_track_pt[key]>=pt_cut_value) * (LC_d0[key]<= d0_cut_value) * (LC_nhits[key]>nhits_cut_value)
        for key in data_loader.keys()
    }

    # Define the eta transition region from barrel to endcap
    transition_region = 1

    # Separate the data into barrel and endcap
    # track_barrel = (np.abs(LC_eta_match)<transition_region)
    truth_barrel = {key:(np.abs(mcp_mu_eta[key])<transition_region) for key in data_loader.keys()}
    # track_endcap = (np.abs(LC_eta_match)>=transition_region)
    truth_endcap = {key:(np.abs(mcp_mu_eta[key])>=transition_region) for key in data_loader.keys()}

    muon_gun_label = [r'Muon particle gun', r'uniform in $p_T$ and $\theta$;']

    # Binned in theta
    bottom_label = muon_gun_label
    print('Computing reconstruction efficiency as a function of theta.')

    custom_bins= [30.*np.pi/180.,
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
    xlim = (custom_bins[0],custom_bins[-1])
    if(degrees):
        custom_bins=None
        xlim = None

    # NOTE: a bit messy, trying to make an "easy" toggle for whether or not before/after cleaning is shown
    mask_pairs = {key:[
        (None,None), # no cleaning,
        (track_clean[key],None) # track cleaning
    ] for key in data_loader.keys()}

    labels={key:["{}: Before Cleaning".format(key), "{}: After Cleaning".format(key)] for key in data_loader.keys()}

    if(efficiency_opts == 1):
        mask_pairs = {key:[mask_pairs[key][0]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][0].split(':')[0]] for key in labels.keys()} # if before cleaning, don't write anything
    elif(efficiency_opts == 2):
        mask_pairs = {key:[mask_pairs[key][1]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][1].split(':')[0]] for key in labels.keys()} # if after cleaning, don't write anything
        bottom_label += ['After cleaning']

    result_dict = {key:calculate_efficiencies(
        LC_track_theta[key],
        mcp_mu_theta[key],
        mask_pairs=mask_pairs[key],
        custom_bins=custom_bins
    ) for key in data_loader.keys()}

    results   = {key:result_dict[key][0] for key in result_dict.keys()}
    min_value = {key:result_dict[key][1] for key in result_dict.keys()}
    max_value = {key:result_dict[key][2] for key in result_dict.keys()}

    # turn things into lists for plot_efficiencies()
    result_list = []
    label_list = []
    for key in data_loader.keys():
        result_list += results[key]
        label_list += labels[key]
    min_value = min_value[list(data_loader.keys())[0]]
    max_value = max_value[list(data_loader.keys())[0]]

    xlabel = r"Muon $\theta$ [rad] "
    if(degrees):
        xlabel = r"Muon $\theta [\degree]$ "
    plotter.plot_efficiencies(result_list, min_value, max_value,
                    xlabel=xlabel,
                    labels=label_list,
                    savename='eff_vs_theta_{}'.format(suffix),
                    xlim=xlim,
                    ylim=(0.95,1.05),
                    bottom_label=bottom_label,
                    label_block_y_up=0.8
                    )

    # print('Produced efficiency vs. theta plots. Number of events used per histogram is:')
    # for i in range(len(result_list)):
    #     print('\t{} : {}'.format(label_list[i],result_list[i].GetTotalHistogram().GetEntries()))


    # Binned in pT, split into barrel and endcap regions
    print('Computing reconstruction efficiency as a function of pT, for barrel region.')
    # custom_bins = [1,2,5,10,20,50,100,200,500,1000,2000,5000]
    custom_bins = [1,2,5,10,20,50,100,200,500,1000]

    mask_pairs = {key:[
        (truth_barrel[key],truth_barrel[key]), # no cleaning,
        (combine_masks([track_clean[key],truth_barrel[key]]),truth_barrel[key]) # track cleaning
    ] for key in data_loader.keys()}

    if(efficiency_opts == 1):
        mask_pairs = {key:[mask_pairs[key][0]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][0].split(':')[0]] for key in labels.keys()} # if before cleaning, don't write anything
    elif(efficiency_opts == 2):
        mask_pairs = {key:[mask_pairs[key][1]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][1].split(':')[0]] for key in labels.keys()} # if after cleaning, don't write anything
        bottom_label += ['After cleaning']

    result_dict = {key:calculate_efficiencies(
        LC_pt_match[key],
        mcp_mu_pt[key],
        mask_pairs=mask_pairs[key],
        custom_bins=custom_bins
    ) for key in data_loader.keys()}

    results   = {key:result_dict[key][0] for key in result_dict.keys()}
    min_value = {key:result_dict[key][1] for key in result_dict.keys()}
    max_value = {key:result_dict[key][2] for key in result_dict.keys()}

    # turn things into lists for plot_efficiencies()
    result_list = []
    label_list = []
    for key in data_loader.keys():
        result_list += results[key]
        label_list += labels[key]
    min_value = min_value[list(data_loader.keys())[0]]
    max_value = max_value[list(data_loader.keys())[0]]

    misctext = r'$|\eta|<1$'
    if(degrees):
        misctext = r'$40^{\circ}<\theta<140^{\circ}$'
    xlabel = "Muon $p_T$ [GeV]"

    plotter.plot_efficiencies(result_list, min_value, max_value,
                    xlabel=xlabel,
                    labels=label_list,
                    savename='eff_vs_pt_barrel_{}'.format(suffix),
                    ylim=(0.95,1.05),
                    bottom_label=bottom_label,
                    misctext=misctext
                    # label_block_y_up=0.8
                    )

    print('Computing reconstruction efficiency as a function of pT, for endcap region.')

    mask_pairs = {key:[
        (truth_endcap[key],truth_endcap[key]), # no cleaning,
        (combine_masks([track_clean[key],truth_endcap[key]]),truth_endcap[key]) # track cleaning
    ] for key in data_loader.keys()}

    if(efficiency_opts == 1):
        mask_pairs = {key:[mask_pairs[key][0]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][0].split(':')[0]] for key in labels.keys()} # if before cleaning, don't write anything
    elif(efficiency_opts == 2):
        mask_pairs = {key:[mask_pairs[key][1]] for key in mask_pairs.keys()}
        labels = {key:[labels[key][1].split(':')[0]] for key in labels.keys()} # if after cleaning, don't write anything
        bottom_label += ['After cleaning']

    result_dict = {key:calculate_efficiencies(
        LC_pt_match[key],
        mcp_mu_pt[key],
        mask_pairs=mask_pairs[key],
        custom_bins=custom_bins
    ) for key in data_loader.keys()}

    results   = {key:result_dict[key][0] for key in result_dict.keys()}
    min_value = {key:result_dict[key][1] for key in result_dict.keys()}
    max_value = {key:result_dict[key][2] for key in result_dict.keys()}

    # turn things into lists for plot_efficiencies()
    result_list = []
    label_list = []
    for key in data_loader.keys():
        result_list += results[key]
        label_list += labels[key]
    min_value = min_value[list(data_loader.keys())[0]]
    max_value = max_value[list(data_loader.keys())[0]]

    misctext = r'$|\eta|>1$'
    if(degrees):
        misctext = r'$\theta<40^{\circ}$ or $\theta>140^{\circ}$'
    xlabel = "Muon $p_T$ [GeV]"

    plotter.plot_efficiencies(result_list, min_value, max_value,
                    xlabel=xlabel,
                    labels=label_list,
                    savename='eff_vs_pt_endcap_{}'.format(suffix),
                    ylim=(0.95,1.05),
                    bottom_label=bottom_label,
                    misctext=misctext
                    # label_block_y_up=0.8
                    )

    return


if(__name__=='__main__'):
    main(sys.argv)