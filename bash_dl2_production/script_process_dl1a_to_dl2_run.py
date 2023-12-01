import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from datetime import datetime
import pickle, json, sys, os, glob
import pandas as pd
import subprocess

import utils

from traitlets.config.loader import Config
from astropy.coordinates     import SkyCoord
from lstchain.io.config      import get_standard_config
from ctapipe.io              import read_table


# Source specifications
source_name = "crab"

# Root path of this script
root = "/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/bash_dl2_production/"
# Path to store the configuration file we are going to use
config_file = root + "objects/standard_config.json"
# Data main directory
root_data = root + f"../../data/cosmic_ray_data_correction/{source_name}/"

# STANDARD paths ---------
# DL1 data root
dl1_root = "/fefs/aswg/data/real/DL1/*/v0.*/tailcut84/"
# RFs root main directory
rfs_root = "/fefs/aswg/data/models/AllSky/20230901_v0.10.4_allsky_base_prod/"
# MCs dl2 main directory
mcs_root = "/fefs/aswg/data/mc/DL2/AllSky/20230901_v0.10.4_allsky_base_prod/TestingDataset/"


# directories for the data
dir_dl1b = root_data + "dl1b/"
dir_dl2  = root_data + "dl2/"
    

def run(string):
    
    run_numbers = [int(string.split("-")[1])]
    subrun_sel  = int(string.split("-")[2])
    scales      = np.array(string.split("-")[-1].split(",")).astype(float)

    # Creating the directories in case they don't exist
    for path in [os.path.dirname(config_file), dir_dl1b, dir_dl2]:
        if not os.path.exists(path):
            os.makedirs(os.path.join(path), exist_ok=True)

    config_dict = get_standard_config()
    # print(config_dict)

    #-------------------
    # Changes in the configuration should be done here

    # We select the heuristic flatfield option in the standard configuration
    config_dict["source_config"]["LSTEventSource"]["use_flatfield_heuristic"] = True

    #-------------------

    with open(config_file, 'w') as json_file:
        json.dump(config_dict, json_file)



    # Getting coordinates of source
    source_coords = SkyCoord.from_name(source_name)

    dict_source = {
        "name"   : source_name,
        "coords" : source_coords,
        "ra"     : source_coords.ra.deg  * u.deg, # ra in degrees
        "dec"    : source_coords.dec.deg * u.deg, # dec in degrees
    }

    # We create a empty dictionary to store all the information needed inside
    DICT = {}
    for run in run_numbers:
        DICT[run] = {
            "run_num" : run,
            "errors"  : "", # log of errors trough the analysis
        }

    DICT = utils.add_dl1_paths_to_dict(DICT, dl1_root)
    DICT = utils.add_dl1_paths_to_dict(DICT, dl1_root, dchecking=True)

    for run in run_numbers:

        tab = read_table(DICT[run]["dchecks"]["runwise"], "/dl1datacheck/cosmics")

        # reading the variables
        _zd,     _az       = 90 - np.rad2deg(np.array(tab["mean_alt_tel"])), np.rad2deg(np.array(tab["mean_az_tel"]))
        _t_start, _t_elapsed = tab["dragon_time"][0][0],                       np.array(tab["elapsed_time"])

        DICT[run]["time"] = {
            "tstart"   : _t_start,            # datetime object
            "telapsed" : np.sum(_t_elapsed),  # s
            "srunwise" : {
                "telapsed" : _t_elapsed,      # s      
            },
        }
        DICT[run]["pointing"] = {
            "zd" : np.mean(_zd),  # deg
            "az" : np.mean(_az),  # deg
            "srunwise" : {
                "zd" : _zd,       # deg
                "az" : _az,       # deg
            },
        }

    # then we also select the RFs and MC files looking at the nodes available
    DICT, dict_nodes = utils.add_mc_and_rfs_nodes(DICT, rfs_root, mcs_root, dict_source)



    for scale in scales:
        for ir, run in enumerate(DICT.keys()):

            # Creating the directories in case they don't exist
            for path in [dir_dl1b + f"{run:05}", dir_dl2 + f"{run:05}"]:
                if not os.path.exists(path):
                    os.makedirs(os.path.join(path), exist_ok=True)



            sruns = [int(path.split(".")[-2]) for path in DICT[run]["dl1a"]["srunwise"]]
            DICT[run]["dl1b"] = {"srunwise" : []}

            for i, srun in enumerate(sruns):
                
                if srun == subrun_sel:

                    input_fname  = DICT[run]["dl1a"]["srunwise"][i]
                    output_fname = dir_dl1b + f"{run:05}/" + f"dl1_LST-1.Run{run:05}.{srun:04}_s{scale:.4f}.h5"

                    print(f"\nComputing dl1b Run {run:5} Subrun {srun:04} - {i/len(sruns)*100:3.1f}% sruns {ir+1}/{len(DICT.keys())} runs")
                    print(f"--> {output_fname}\n")

                    command = f"lstchain_dl1ab -f {input_fname} -o {output_fname} -c {config_file} --no-image --light-scaling {scale}"
                    
                    # with open('tmp.sh', 'w') as f:
                    #     f.write("#! /bin/bash\n\n")
                    #     f.write()
                    subprocess.run(command, shell=True)
        
                    DICT[run]["dl1b"]["srunwise"] = output_fname

            DICT[run]["dl2"] = {"srunwise" : []}

            for i, srun in enumerate(sruns):
                
                
                if srun == subrun_sel:

                    input_fname  = DICT[run]["dl1b"]["srunwise"]
                    output_fname = dir_dl1b + f"{run:05}/" + input_fname.split("/")[-1].replace("dl1", "dl2", 1)
                    rf_node      = DICT[run]["simulations"]["rf"]

                    print(f"\nComputing dl1b Run {run:5} Subrun {srun:04} - {i/len(sruns)*100:3.1f}% sruns {ir+1}/{len(DICT.keys())} runs")
                    print(f"--> {output_fname}\n")
                    
                    command = f"lstchain_dl1_to_dl2 -f {input_fname} -p {rf_node} -o {dir_dl2} -c {config_file}"
                    subprocess.run(command, shell=True)

                    DICT[run]["dl2"]["srunwise"] = output_fname
                    
                    
if __name__ == "__main__":
    input_string = sys.argv[1]
    run(input_string)