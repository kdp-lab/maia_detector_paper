# Import necessary libraries
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import json
from pathlib import Path
import os
import re

def load_data(file_path):
    """
    Load data from a JSON file.
    """
    which_data = ak.from_json(Path(file_path))
    return which_data

def xyz_to_r_xy_vb(which_data):
    r_xy = np.array([np.array((np.sqrt(a**2 + b**2),c)) for a,b,c in zip(ak.flatten(which_data['x']),ak.flatten(which_data['y']),ak.flatten(which_data['z']))])
    r_xy_vb = r_xy[(ak.flatten(which_data['hit_detector']) == 1)] # | (ak.flatten(which_data['hit_detector']) == 2)]
    r_xyz = np.array([np.array((a,b,c)) for a,b,c in zip(ak.flatten(which_data['x']),ak.flatten(which_data['y']),ak.flatten(which_data['z']))])
    r_xyz_vb = r_xyz[(ak.flatten(which_data['hit_detector']) == 1)] # | (ak.flatten(which_data['hit_detector']) == 2)]
    # Find unique rows in r_xyz_vb
    unique_rows, unique_indices = np.unique(r_xyz_vb, axis=0, return_index=True)
    # Check if the number of unique rows is less than the total number of rows
    
    r_xyz_vb = r_xyz_vb[unique_indices]
    r_xy_vb = r_xy_vb[unique_indices]
    vertex_layer = ak.flatten(which_data['hit_layer'][(which_data['hit_detector'] == 1)])[unique_indices]
    if len(unique_rows) < len(r_xyz_vb):
        duplicate_indices = np.setdiff1d(np.arange(len(r_xyz_vb)), unique_indices)
        print("Duplicates found in dataset:")
        print(len(duplicate_indices))
    else:
        print("No duplicates found in dataset.")
    return r_xy_vb, r_xyz_vb, vertex_layer

def find_hit_pairs(vertex_layer, r_xyz_vb, r_xy_vb=None, num_hits = None, z_tolerance=1):
    """
    Finds pairs of hits within a specified z tolerance between consecutive layers.

    Parameters:
    - which_data: Dataset containing hits information.
    - r_xyz_vb: Numpy array of xyz coordinates of hits.
    - r_xy_vb: Numpy array of xy coordinates of hits.
    - num_hits: Number of hits to consider for pairing.
    - z_tolerance: Z-axis tolerance for considering hits as a pair (default 1mm).

    Returns:
    - doublet_pairs_xyz: Array of hit pairs in xyz coordinates.
    - doublet_pairs_rz: Array of hit pairs in r(z) coordinates.
    """

    doublet_pairs_xyz = []
    #doublet_pairs_rz = []

    if num_hits is None: num_hits = len(r_xyz_vb)

    r_xyz_vb_reduced = r_xyz_vb[:num_hits]
    #r_xy_vb_reduced = r_xy_vb[:num_hits]
    vertex_layer = vertex_layer[:num_hits]

    for layer in range(0, max(vertex_layer), 2):  
        hits_current_layer = r_xyz_vb_reduced[vertex_layer == layer]
        hits_next_layer = r_xyz_vb_reduced[vertex_layer == layer + 1]
        
        for i, hit1 in enumerate(hits_current_layer):
            z_distances = np.abs(hit1[2] - hits_next_layer[:, 2])
            within_tolerance = np.where(z_distances < z_tolerance)[0]
            
            for j in within_tolerance:
                hit2 = hits_next_layer[j]
                doublet_pairs_xyz.append([hit1, hit2])
                #doublet_pairs_rz.append([r_xy_vb_reduced[vertex_layer == layer][i], r_xy_vb_reduced[vertex_layer == layer + 1][j]])

    doublet_pairs_xyz = np.array(doublet_pairs_xyz)
    #doublet_pairs_rz = np.array(doublet_pairs_rz)
    print(f"Total pairs found from {num_hits} hits: {len(doublet_pairs_xyz)}")
    return doublet_pairs_xyz#, doublet_pairs_rz

def xy_to_phi(x, y):
    """
    Convert Cartesian coordinates (x, y) to cylindrical coordinate (phi).
    """
    return np.arctan2(y, x)

def xyz_to_theta(x, y, z):
    """
    Compute the longitudinal angle (theta) between a vector and the z-axis.
    """
    return np.arctan2(np.sqrt(x**2 + y**2), z)

def angular_difference(angle1, angle2):
    """
    Compute the angular difference between two angles with cylindrical symmetry.
    """
    diff = np.abs(angle1 - angle2)
    return np.abs(diff)

import numpy as np

def filter_hits_by_angle(doublet_pairs_xyz, max_theta_tolerance=3e-3, max_phi_tolerance=3e-3):
    """
    Filters pairs of hits based on angular differences within specified tolerances.

    Parameters:
    - doublet_pairs_xyz: Array of hit pairs in xyz coordinates.
    - max_theta_tolerance: Maximum tolerated longitudinal angle difference (in radians).
    - max_phi_tolerance: Maximum tolerated cylindrical coordinate (phi) difference (in radians).

    Returns:
    - d_theta_d_phi: Array of angular differences (theta, phi) for hit pairs within tolerance.
    """
    d_theta_d_phi = []

    for pair in doublet_pairs_xyz:
        hit1_xyz, hit2_xyz = pair
        theta1 = xyz_to_theta(*hit1_xyz)
        theta2 = xyz_to_theta(*hit2_xyz)
        phi1 = xy_to_phi(hit1_xyz[0], hit1_xyz[1])
        phi2 = xy_to_phi(hit2_xyz[0], hit2_xyz[1])

        d_theta = angular_difference(theta1, theta2)
        d_phi = angular_difference(phi1, phi2)

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

