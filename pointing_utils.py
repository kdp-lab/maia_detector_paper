# Import necessary libraries
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os
import re

def load_data(file_path):
    """
    Load data from a JSON file.
    """
    which_data = ak.from_json(Path(file_path))
    return which_data

# Function to remove duplicate hits
def make_unique(r_xy_vb, r_xyz_vb, hits_layer):
    # Find unique rows in r_xyz_vb
    unique_rows, unique_indices = np.unique(r_xyz_vb, axis=0, return_index=True)

    # Check if the number of unique rows is less than the total number of rows
    if len(unique_rows) < len(r_xyz_vb):
        print("Duplicates found in dataset:")
        
        # Find duplicate indices
        duplicate_indices = np.setdiff1d(np.arange(len(r_xyz_vb)), unique_indices)
        print(len(duplicate_indices))
        # Print duplicate rows
        # for idx in duplicate_indices:
        #     print(idx, r_xyz_vb[idx])
    else:
        print("No duplicates found in dataset.")
    r_xyz_vb = r_xyz_vb[unique_indices]
    r_xy_vb = r_xy_vb[unique_indices]
    vertex_layer = (hits_layer[unique_indices])
    return r_xy_vb, r_xyz_vb, vertex_layer

# Function to create r_xy and r_xyz arrays for the vertex barrel region
def create_r_xy(which_data, pt_cut = 100, d0_cut = 1):
    '''
    Returns [r_xy, z] and [x,y,z] for the vertex barrel region for a given dataset and cuts
    '''
    mask = (which_data['track_pt'] > pt_cut) & (which_data['d0_res'] < d0_cut) & (which_data['hit_detector'] == 1)
    hits_x = np.ravel(which_data['x'][mask])
    hits_y = np.ravel(which_data['y'][mask])
    hits_z = np.ravel(which_data['z'][mask])
    hits_layer = np.ravel(which_data['hit_layer'][mask])
    # print(hits_x)
    r_xy_vb = np.array([[np.sqrt(a**2 + b**2), c] for a, b, c in zip(hits_x, hits_y, hits_z)])
    r_xyz_vb = np.array([[a, b, c] for a, b, c in zip(hits_x, hits_y, hits_z)])
    return make_unique(r_xy_vb, r_xyz_vb, hits_layer)

# Function to find doublet pairs
def find_doublets(r_xy_vb, r_xyz_vb, vertex_layer, num_hits = None, z_tolerance = 1):
    if num_hits is None:
        num_hits = len(r_xy_vb) # 10000 takes < 1 min, 20000 takes ~4 mins

    vertex_layer = vertex_layer[:num_hits]
    r_xyz_vb_reduced = r_xyz_vb[:num_hits]
    r_xy_vb_reduced = r_xy_vb[:num_hits]

    doublet_pairs_xyz = []
    doublet_pairs_rz = []

    for layer in range(0, max(vertex_layer), 2):  
        # Identify hits in the current and next layer
        hits_current_layer = r_xyz_vb_reduced[vertex_layer == layer]
        hits_next_layer = r_xyz_vb_reduced[vertex_layer == layer + 1]
        
        # For each hit in the current layer, find a hit in the next layer within z_tolerance
        for i, hit1 in enumerate(hits_current_layer):
            #print(i, hit1)
            # Calculate the z-distance between hit1 and all hits in the next layer
            z_distances = np.abs(hit1[2] - hits_next_layer[:, 2])
            
            # Find indices of hits in the next layer within the z_tolerance
            within_tolerance = np.where(z_distances < z_tolerance)[0]
            
            # Pair hit1 with all hits within tolerance
            for j in within_tolerance:
                hit2 = hits_next_layer[j]
                doublet_pairs_xyz.append([hit1, hit2])
                # Assuming r_xy_vb_reduced is structured similarly to r_xyz_vb_reduced, index with the same logic
                doublet_pairs_rz.append([r_xy_vb_reduced[vertex_layer == layer][i], r_xy_vb_reduced[vertex_layer == layer + 1][j]])

    # Convert lists to arrays for further processing if necessary
    doublet_pairs_xyz = np.array(doublet_pairs_xyz)
    doublet_pairs_rz = np.array(doublet_pairs_rz)
    print(f"Total pairs found from {num_hits} hits: {len(doublet_pairs_xyz)}")

    return doublet_pairs_xyz, doublet_pairs_rz

# Function to compute the cylindrical coordinates (phi) of a vector
def xy_to_phi(x, y):
    return np.arctan2(y, x)

# Function to compute the longitudinal angle (theta) between a vector and the z-axis
def xyz_to_theta(x, y, z):
    return np.arctan2(np.sqrt(x**2 + y**2), z)

# Function to compute the angular difference between two angles (in radians) with cylindrical symmetry
def angular_difference(angle1, angle2):
    diff = np.abs(angle1 - angle2)
    # Adjust the difference to be within the range [-pi, pi)
    # diff = (diff + np.pi) % (2 * np.pi) - np.pi
    return np.abs(diff)

# Function to compute the pointing angles between pairs of hits
def pointing_angles(doublet_pairs_xyz, max_theta_tolerance = 3e-3, max_phi_tolerance = 35e-3):
    d_theta_d_phi = []
    # Iterate over each pair of hits in doublet_pairs_xyz
    for pair in doublet_pairs_xyz:
        # Extract the Cartesian coordinates of the two hits
        hit1_xyz, hit2_xyz = pair
        
        # Compute the longitudinal angle (theta) and cylindrical coordinate (phi) for each hit
        theta1 = xyz_to_theta(*hit1_xyz)
        theta2 = xyz_to_theta(*hit2_xyz)
        phi1 = xy_to_phi(hit1_xyz[0], hit1_xyz[1])
        phi2 = xy_to_phi(hit2_xyz[0], hit2_xyz[1])
        
        # Compute the angular differences
        d_theta = angular_difference(theta1, theta2)
        d_phi = angular_difference(phi1, phi2)
        
        # Filter based on the maximum tolerated angles
        if d_theta <= max_theta_tolerance and d_phi <= max_phi_tolerance:
            d_theta_d_phi.append(np.array([d_theta, d_phi]))

    d_theta_d_phi = np.array(d_theta_d_phi)

    return d_theta_d_phi


def save2DHistogram(datax, datay, fname, bins=100, weights=None, norm="log", xlim=None, ylim=None):
    data_flatx = np.array(np.ravel(datax)).T
    data_flaty = np.array(np.ravel(datay)).T
    
    if weights is not None:
        values = np.array(np.ravel(weights)).T
    else:
        values = None

    fig, ax = plt.subplots(figsize=(8,6))
    hh = ax.hist2d(data_flatx, data_flaty, bins=bins, weights=values, norm=norm, cmap="viridis")

    pattern = r'(_sim|_digi|_digi_bib|_reco|_reco_bib)\.json$'
    base_filename = os.path.basename(fname)
    #print(base_filename)
    input_filename = re.sub(pattern, '', base_filename)
    #print(input_filename)
    mass_str, lifetime_str = input_filename.split('_')
    #print(mass_str, lifetime_str)
    title = f"Mass: {mass_str} GeV, Lifetime: {lifetime_str} ns"
    
    ax.set_title(title)
    ax.set_xlabel("$\Delta \\theta$ [mrad]", fontsize=16, loc="right")
    ax.set_ylabel(r"$\Delta \phi$ [mrad]", fontsize=16, loc="top")
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    
    cbar = fig.colorbar(hh[3])
    #cbar.ax.tick_params(labelsize=tick_font_size)
    cbar.set_label("")#, fontsize=6)

    output_dir = r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_graphs"
    output_path = os.path.join(output_dir, f"{input_filename}.png")
    
    print(f"Saving histogram to {output_path}")
    plt.savefig(output_path)
    plt.close(fig)  # Close the figure to free up memory

