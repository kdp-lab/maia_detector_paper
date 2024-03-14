import pointing_utils as pu
import numpy as np
from pathlib import Path
import multiprocessing as mp
import os
import re

# To run in powershell: 
# Start-Process -NoNewWindow -FilePath "C:\Users\leoro\anaconda3\python.exe" -ArgumentList "`"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_analysis.py`"" -RedirectStandardOutput "C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_log.log" -RedirectStandardError "C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_log_error.log"
# Start-Process -NoNewWindow 
# -FilePath "C:\Users\leoro\anaconda3\python.exe" 
# -ArgumentList "`"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_analysis.py`"" 
# -RedirectStandardOutput "C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_log.log" 
# -RedirectStandardError "C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\pointing_log_error.log"

def process_dataset(data_file_path, num_hits=None, z_tolerance=1):
    """
    Process a single dataset.
    """
    # Load the data
    which_data = pu.load_data(data_file_path)
    print(f"Processed data {data_file_path}")

    # Convert loaded data to r, xy, vb, and extract vertex layers
    r_xy_vb, r_xyz_vb, vertex_layer = pu.xyz_to_r_xy_vb(which_data)
    print(f"Extracted r_xyz_vb from {data_file_path}")

    # Find hit pairs with a specific z tolerance
    doublet_pairs_xyz = pu.find_hit_pairs(vertex_layer, r_xyz_vb, num_hits=num_hits, z_tolerance=z_tolerance) # z_tolerance in mm
    print(f"Found pairs in {data_file_path}")

    # Filter hit pairs based on angular differences within specified tolerances
    d_theta_d_phi = pu.filter_hits_by_angle(doublet_pairs_xyz)
    print(f"Found {len(doublet_pairs_xyz)} filtered pairs in {data_file_path}")

    pattern = r'(_sim|_digi|_digi_bib|_reco|_reco_bib)\.json$'
    base_filename = os.path.basename(data_file_path)
    input_filename = re.sub(pattern, '', base_filename)
    # Extract directory from data_file_path
    directory = r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data"
    # Construct the new filename for saving
    if '_bib' in str(data_file_path):
        save_filename = os.path.join(directory, f"{input_filename}_d_theta_d_phi_bib.npy")
    else:
        save_filename = os.path.join(directory, f"{input_filename}_d_theta_d_phi.npy")
    # Save d_theta_d_phi array
    np.save(save_filename, d_theta_d_phi)
    print(f"Saved {save_filename}")
    # pu.save2DHistogram(d_theta_d_phi[:,0]*1000, d_theta_d_phi[:,1]*1000, fname = data_file_path, bins=np.linspace(0,1,200))

def main():
    # List of paths to your JSON data files
    data_files = [
        Path(r"X:\local\d1\lrozanov\mucoll-tutorial-2023\reco_Hbb\134_10_reco.json"),
        Path(r"X:\local\d1\lrozanov\mucoll-tutorial-2023\reco_Hbb\134_0.1_reco.json"),
        Path(r"X:\local\d1\lrozanov\mucoll-tutorial-2023\reco_Hbb\134_1_reco.json")
        ,Path(r"X:\local\d1\lrozanov\mucoll-tutorial-2023\reco_Hbb\300_10_reco.json")
        # Add more datasets as needed
        ]

    # Determine the number of processes to use
    num_processes = min(len(data_files), 4)
    num_hits = None
    z_tolerance = 1

    # Create a pool of workers to process the datasets in parallel
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(process_dataset, [(data_file, num_hits, z_tolerance) for data_file in data_files]) 

if __name__ == "__main__":
    main()

