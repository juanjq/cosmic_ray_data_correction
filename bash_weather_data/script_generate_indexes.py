import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from datetime import datetime, timedelta
import pickle, json, sys, os, glob
import pandas as pd
pd.set_option("display.max_columns", None)

# location of the scripts
sys.path.insert(0, os.path.join("/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/scripts"))
import auxiliar as aux
import geometry as geom

dcheck_root = "/fefs/aswg/workspace/abelardo.moralejo/data/datachecks/night_wise/DL1_datacheck_"
ws_database = "/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/analysis_first_corrections/objects/WS2003-22_short.h5"
dir_objects = "/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/analysis_weather/objects"
results_path = "/fefs/aswg/workspace/juan.jimenez/cosmic_ray_data_correction/bash_weather_data/RESULTS.txt"

def run(string):

    init = int(string.split(",")[0])
    end  = int(string.split(",")[1])
    
    with open(dir_objects + "/data_dict.pkl", 'rb') as f:
        dict_dcheck = pickle.load(f)

    df_ws = pd.read_hdf(ws_database)
    dates = dict_dcheck["time"]

    maxdate = np.max(dates)
    mindate = np.min(dates)

    # Assuming df_ws.index is already a NumPy array
    dates_ws = np.array([datetime.fromisoformat(str(d).split(".")[0]) for d in df_ws.index])

    # Combine date filtering in NumPy
    mask = (dates_ws > mindate) & (dates_ws < maxdate)
    dates_ws = dates_ws[mask]
    df_ws = df_ws[mask]

    print("Starting to iterate")
    
    indexes = []
    ref     = []
    for j, date in enumerate(dates):
        
        if j >= init and j <= end:
            
            if j % 200 == 0:
                print(f"{j}/{6000}")

            ref.append(j)
            indexes.append(np.argmin(np.abs(dates[j] - dates_ws)))   

    print("Writting...")
    with open(results_path, "a") as f:
        for i, r in zip(indexes, ref):
            f.write(f"{i},{r}\n")
                    
if __name__ == "__main__":
    input1 = sys.argv[1]
    run(input1)