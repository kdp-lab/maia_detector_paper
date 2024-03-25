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

def process_dataset(data_file_path, **kwargs):
    """
    Process a single dataset.
    """
    # Load the data
    which_data = pu.load_data(data_file_path)
    print(f"Processed data {data_file_path}")

    # Unpack kwargs for create_r_xy specifically or adjust parameters accordingly
    r_xy_vb, r_xyz_vb, vertex_layer = pu.create_r_xy(
        which_data, 
        pt_cut=kwargs.get('pt_cut', 100),  # Provide default values if not specified
        d0_cut=kwargs.get('d0_cut', 1)
    )
    # Unpack kwargs for find_doublets
    doublet_pairs_xyz, doublet_pairs_rz = pu.find_doublets(
        r_xy_vb, 
        r_xyz_vb, 
        vertex_layer, 
        num_hits=kwargs.get('num_hits', 10000),  # Provide default values if not specified
        z_tolerance=kwargs.get('z_tolerance', 1)
    )
    # Now call the revised d_theta_d_phi function with correct parameters
    d_theta_d_phi = pu.pointing_angles( 
        doublet_pairs_xyz, 
        max_theta_tolerance=kwargs.get('max_theta_tolerance', 3e-3), 
        max_phi_tolerance=kwargs.get('max_phi_tolerance', 35e-3)
    )
    print(f"The hit survival rate: {len(d_theta_d_phi)/len(r_xy_vb)}")

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

kwargs = {
            'num_hits': None,
            'z_tolerance': 1,
            'max_theta_tolerance': 3e-3,
            'max_phi_tolerance': 35e-3,
            'pt_cut': 100,
            'd0_cut': 1
        }

def wrapper(data_file):
    # Call process_dataset with **kwargs
    process_dataset(data_file, **kwargs)

def main():
    # List of paths to your JSON data files
    data_files = [
        Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_10_reco.json"),
        Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_0.1_reco.json"),
        Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_1_reco.json"),
        Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\300_10_reco.json")
        # Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_0.1_reco_bib.json"),
        # Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_1_reco_bib.json"),
        # Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\300_10_reco_bib.json"),
        # Path(r"C:\Users\leoro\OneDrive\Muon_Collider\Psets\mucol_plotting_scripts-main\macros\local_scripts\LLPs_data\134_0.1_reco_bib.json"),
        ]

    # Determine the number of processes to use
    num_processes = min(len(data_files), 4)

    # Create a pool of workers to process the datasets in parallel
    with mp.Pool(processes=num_processes) as pool:
        pool.starmap(wrapper, [(data_file,) for data_file in data_files])

if __name__ == "__main__":
    main()

