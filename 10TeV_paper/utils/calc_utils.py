import numpy as np
import awkward as ak

def combine_masks(mask_list):
    """
    There is probably a better way to do this, but having issues with broadcasting
    multiplication of jagged boolean awkward arrays, even if dimensions are matched.
    """
    # TODO: Broken
    # nevents = len(mask_list[0])
    mask = ak.copy(mask_list[0])
    for i,entry in enumerate(mask_list):
        if(i==0): continue
        for j in range(len(mask)):
            for k in range(len(mask[j])):
                mask[j,k] = mask[j,k] * entry[j,k]
    return mask

def calculate_efficiences_v2(track_data, truth_data, mask_pairs, num_bins=10, min_value=None, max_value=None, custom_bins=None):
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
    efficiency_bins = np.linspace(min_value, max_value, num_bins+1)

    # print('efficiency_bins = ',efficiency_bins)
    if custom_bins is not None:
        efficiency_bins = np.array(custom_bins)
        min_value = efficiency_bins[0]
        max_value = efficiency_bins[-1]
        num_bins = len(efficiency_bins) - 1

    bin_centers = []
    for i,entry in enumerate(efficiency_bins):
        if(i == 0): continue
        bin_centers.append((entry + efficiency_bins[i-1])/2)

    #bin_width = efficiency_bins[1] - efficiency_bins[0]
    bin_widths = (efficiency_bins[1:] - efficiency_bins[:-1]) / 2

    results = []

    for mask_pair in mask_pairs:

        track_mask, truth_mask = mask_pair

        numerator = np.zeros(num_bins)
        denominator = np.zeros(num_bins)

        nevents = len(track_data)

        if(truth_mask is None):
            truth_mask = np.full(nevents,True)

        if(track_mask is None):
            track_mask = np.full(nevents,True)

        for i in range(nevents): # looping thru awkward arrays

            if(not(truth_mask[i])): continue

            track = track_data[i] # in practice, might be empty -- if we failed to get a track
            muon = truth_data[i]
            bin_idx = np.digitize(muon,efficiency_bins) - 1

            denominator[bin_idx] += 1
            if(len(track) == 0): continue

            if(track_mask[i]):
                numerator[bin_idx] += 1

        efficiency = np.copy(numerator / denominator)

        uncertainty = efficiency * np.sqrt((1 - efficiency) / numerator) # TODO: check this

        # package things in a list, since this is how old code did it
        results.append((bin_centers, efficiency, uncertainty, bin_widths))
    # results.append((bin_centers, efficiency_clean, uncertainty_clean, bin_widths))
    return results, min_value, max_value

# Efficiencies calculated using fraction of truth-matched tracks/truth tracks in each bin
def calculate_efficiencies(track_data_list, truth_data_list, num_bins=10, min_value=None, max_value=None, custom_bins=None):

    # Can assign min and max values to the data, don't have to
    if min_value is None:
        min_value = np.min([np.min(np.ravel(truth_data)) for truth_data in truth_data_list])

    if max_value is None:
        max_value = np.max([np.max(np.ravel(truth_data)) for truth_data in truth_data_list])

    # Bins based on min and max values
    efficiency_bins = np.linspace(min_value, max_value, num_bins+1)

    if custom_bins is not None:
        efficiency_bins = np.array(custom_bins)
        min_value = efficiency_bins[0]
        max_value = efficiency_bins[-1]

    results = []
    #bin_width = efficiency_bins[1] - efficiency_bins[0]
    bin_widths = (efficiency_bins[1:] - efficiency_bins[:-1]) / 2
    # Loop over the track and truth data
    for track_data, truth_data in zip(track_data_list, truth_data_list):

        # print('len(track_data) = {}, len(truth_data) = {}'.format(len(track_data),len(truth_data)))

        efficiencies = []
        errors = []
        bin_centers = []
        # For each bin, assing the min and max values and then calculate the efficiency
        for j in range(len(efficiency_bins) - 1):
            bin_min = efficiency_bins[j]
            bin_max = efficiency_bins[j + 1]
            track_data_in_bin = np.ravel(track_data)[(np.ravel(track_data) >= bin_min) & (np.ravel(track_data) < bin_max)]
            truth_data_in_bin = np.ravel(truth_data)[(np.ravel(truth_data) >= bin_min) & (np.ravel(truth_data) < bin_max)]
            if((len(truth_data_in_bin) != 0) and (len(track_data_in_bin) != 0)):
                efficiency = len(track_data_in_bin) / len(truth_data_in_bin)
                if(efficiency > 1.):
                    print('Warning: Found efficiency > 100% for bin {}'.format(j))
                error = efficiency * np.sqrt((1 - efficiency) / len(track_data_in_bin))
            # If bin is empty, efficiency is 0
            else:
                efficiency = 0
                error = 0
            efficiencies.append(efficiency)
            errors.append(error)
            bin_centers.append((bin_min + bin_max) / 2)
        results.append((bin_centers, efficiencies, errors, bin_widths))

    return results, min_value, max_value
