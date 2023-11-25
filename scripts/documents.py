import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle, sys, os, json

def load_gammapy_analysis_configuration(path="/fefs/aswg/workspace/juan.jimenez/lst1_systematics/analysis_first_corrections/config/config_gammapy_analysis.json", Print=True):

    with open(path, "r") as json_file:
        dict = json.load(json_file)

    target_name   = dict["target_name"]
    n_off_regions = dict["n_off_regions"]
    _e_reco = dict["e_reco"]
    _e_true = dict["e_true"]

    if Print == True:
        print("Gammapy analysis configuration:\n")
        display(dict)
        
    return target_name, n_off_regions, _e_reco, _e_true