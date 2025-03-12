import numpy as np

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
        efficiencies = []
        errors = []
        bin_centers = []
        # For each bin, assing the min and max values and then calculate the efficiency
        for j in range(len(efficiency_bins) - 1):
            bin_min = efficiency_bins[j]
            bin_max = efficiency_bins[j + 1]
            track_data_in_bin = np.ravel(track_data)[(np.ravel(track_data) >= bin_min) & (np.ravel(track_data) < bin_max)]
            truth_data_in_bin = np.ravel(truth_data)[(np.ravel(truth_data) >= bin_min) & (np.ravel(truth_data) < bin_max)]
            if len(track_data_in_bin) != 0:
                efficiency = len(track_data_in_bin) / len(truth_data_in_bin)
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
