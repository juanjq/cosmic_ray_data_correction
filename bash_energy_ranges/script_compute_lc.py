import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import pickle, sys, os, json
import astropy.units as u
import pandas as pd
pd.set_option("display.max_columns", None)

from astropy.coordinates import SkyCoord
from matplotlib.dates    import DayLocator, MonthLocator, DateFormatter
from regions             import PointSkyRegion
from astropy.time        import Time
from scipy.stats         import chi2

from gammapy.modeling.models import create_crab_spectral_model, SkyModel, LogParabolaSpectralModel
from gammapy.estimators      import FluxPointsEstimator, LightCurveEstimator, FluxPoints
from gammapy.modeling        import Fit
from gammapy.datasets        import Datasets, SpectrumDataset
from gammapy.makers          import SpectrumDatasetMaker, WobbleRegionsFinder, ReflectedRegionsBackgroundMaker, SafeMaskMaker
from gammapy.maps            import MapAxis, RegionGeom, Map, TimeMapAxis
from gammapy.data            import DataStore

# import scripts
sys.path.insert(0, os.path.join("/fefs/aswg/workspace/juan.jimenez/lst1_systematics/scripts"))
import auxiliar  as aux
import documents as docs

# ============================ #
# dl3 path where dl3 and index files are
dl3_dir = "/fefs/aswg/workspace/juan.jimenez/data/lst1_systematics/dl3"

dicts_dir = "objects/dicts_sed_and_lc/"
# ============================ #


def calculate_chi2_pvalue(table, sys_error=0):
    uncertainty = np.sqrt((sys_error * table["flux"])**2 + table["flux_err"]**2)
    flux = table["flux"]
    mean_flux = (flux/uncertainty**2).sum() / (1/uncertainty**2).sum()
    mean_flux_err = np.sqrt(1/np.sum(1/uncertainty**2))
    print(f"Weighted mean flux: {mean_flux:.3e} +/- {mean_flux_err:.3e} cm-2 s-1")
    
    chi2_value = np.sum((table["flux"] - mean_flux)**2/uncertainty**2)
    ndf = len(table["flux"]) - 1
    pvalue = chi2.sf(x=chi2_value, df=ndf)
    print(f"Chi2: {chi2_value:.1f}, ndf: {ndf}, P-value: {pvalue:.2e}")
    return chi2_value, ndf, pvalue


def weighted_average(table, sys_error=0):
    val = table["flux"]
    uncertainty = np.sqrt((sys_error * table["flux"])**2 + table["flux_err"]**2)
    return (val/uncertainty**2).sum() / (1/uncertainty**2).sum(), np.sqrt(1/np.sum(1/uncertainty**2))



def compute(e_ranges_str):
    
    
    e_lc_min, e_lc_max = np.array(e_ranges_str.split(",")).astype(float)
    
    
    fname_dict = os.path.join(dicts_dir, f"dict_{e_lc_min:.2f}_{int(e_lc_max)}.pkl")
    
    # reading the configuration from the gammapy configuration file
    target_name, n_off_regions, _e_reco, _e_true = docs.load_gammapy_analysis_configuration(Print=False)

    e_reco_min, e_reco_max, e_reco_bin_p_dec = e_lc_min, e_lc_max, _e_reco["bins_p_dec"]
    e_true_min, e_true_max, e_true_bin_p_dec = _e_true["min"], _e_true["max"], _e_true["bins_p_dec"]    
    
    # Opening all the dl3 data in a path
    total_data_store = DataStore.from_dir(dl3_dir)

    # Taking obs ids
    obs_ids = total_data_store.obs_table["OBS_ID"].data
    obs_ids = obs_ids[:]

    # Then we get the observation information from the total data store
    observations = total_data_store.get_observations(
        obs_ids,
        required_irf=["aeff", "edisp", "rad_max"]
    )

    # Defining target position and ON reion
    target_position = SkyCoord.from_name(target_name, frame='icrs')
    on_region = PointSkyRegion(target_position)

    # ============================ #
    # estimated energy axes
    energy_axis = MapAxis.from_energy_bounds(
        e_reco_min, e_reco_max, 
        nbin=e_reco_bin_p_dec, per_decade=True, 
        unit="TeV", name="energy"
    )
    # ============================ #
    # estimated energy axes
    energy_axis_true = MapAxis.from_energy_bounds(
        e_true_min, e_true_max, 
        nbin=e_true_bin_p_dec, per_decade=True, 
        unit="TeV", name="energy_true"
    )
    # ============================ #
    # Energy for the spectrum
    e_fit_min = energy_axis.edges[0].value
    e_fit_max = energy_axis.edges[-1].value
    e_fit_bin_p_dec = e_reco_bin_p_dec

    # Just to have a separate MapAxis for spectral fit energy range
    energy_fit_edges = MapAxis.from_energy_bounds(
        e_fit_min, e_fit_max, 
        nbin=e_fit_bin_p_dec, per_decade=True, 
        unit="TeV"
    ).edges

    # ============================ #
    # Energy for the lightcurve
    e_lc_min = energy_axis.edges[0]
    e_lc_max = energy_axis.edges[-1]

    print("Spectral fit will be done in energy edges:\n", energy_fit_edges)
    print(f"\nLC will be estimated from {e_lc_min:.1f} to {e_lc_max:.1f}")
    

    # geometry defining the ON region and SpectrumDataset based on it
    geom = RegionGeom.create(
        region=on_region, 
        axes=[energy_axis]
    )

    # creating an empty dataset
    dataset_empty = SpectrumDataset.create(
        geom=geom, 
        energy_axis_true=energy_axis_true
    )
    dataset_maker = SpectrumDatasetMaker(
        containment_correction=False,
        selection=["counts", "exposure", "edisp"]
    )

    # tell the background maker to use the WobbleRegionsFinder
    region_finder = WobbleRegionsFinder(n_off_regions=n_off_regions)
    bkg_maker = ReflectedRegionsBackgroundMaker(region_finder=region_finder)
    
    
    # The final object will be stored as a Datasets object
    datasets = Datasets()

    for observation in observations:
        dataset = dataset_maker.run(
            dataset=dataset_empty.copy(name=str(observation.obs_id)),
            observation=observation
        )
        dataset_on_off = bkg_maker.run(
            dataset=dataset, 
            observation=observation
        )
        datasets.append(dataset_on_off) 

    # Stacking all the datasets in one
    stacked_dataset = Datasets(datasets).stack_reduce()    
    
    # defining the model we want to fit and the starting values
    spectral_model = LogParabolaSpectralModel(
        amplitude=1e-12 * u.Unit("cm-2 s-1 TeV-1"),
        alpha=2,
        beta=0.1,
        reference=1 * u.TeV,
    )
    # we will use the crab model in general
    model = SkyModel(
        spectral_model=spectral_model, 
        name="crab"
    )

    # We set the model of all datasets to log parabola
    stacked_dataset.models = model

    # Now we run the fit to extract the parameters of the model
    fit = Fit()
    result = fit.run(datasets=stacked_dataset)
    best_fit_model = model.copy()


    # then extracting the flux points from the data
    fpe = FluxPointsEstimator(
        energy_edges=energy_fit_edges, 
        source=target_name, 
        selection_optional="all"
    )

    # We apply the flux point estiation from the datasets
    flux_points = fpe.run(datasets=stacked_dataset)

    model.parameters["alpha"].frozen = True
    model.parameters["beta"].frozen  = True

    # Create the LC Estimator for each run
    lc_maker_1d = LightCurveEstimator(
        energy_edges=[e_lc_min, e_lc_max], 
        reoptimize=False, # Re-optimizing other free model parameters (not belonging to the source)
        source="crab", 
        selection_optional="all" # Estimates asymmetric errors, upper limits and fit statistic profiles
    )

    # Assigning the fixed parameters model to each dataset
    for data in datasets:
        data.models = model

    lc_runwise = lc_maker_1d.run(datasets)
    lightcurve = lc_runwise.to_table(sed_type="flux", format="lightcurve")
    

    mean_flux, mean_flux_err = weighted_average(lightcurve)

    chi2_val, ndf, pvalue = calculate_chi2_pvalue(lightcurve, sys_error=0.0)   
    
    
    # Start time, duration and central time
    time_min = Time(np.hstack(lightcurve["time_min"]), format='mjd').datetime
    time_max = Time(np.hstack(lightcurve["time_max"]), format='mjd').datetime
    delta_time  = time_max - time_min
    time_center = time_min + delta_time / 2

    # Flux and flux error
    flux_lst1 = np.hstack(lightcurve["flux"])
    flux_stat_err_lst1 = np.hstack(lightcurve["flux_err"])

    # run numbers
    run_num = [int(n) for n in observations.ids]   
    
    
    crab = create_crab_spectral_model("magic_lp")

    crab.amplitude.error = 0.03e-11 * u.Unit("cm-2 s-1 TeV-1")
    crab.alpha.error = 0.01
    crab.beta.error = 0.01/np.log(10)


    flux_crab = crab.integral(e_lc_min, e_lc_max)
    flux_crab_error = flux_crab * 0    
    
    
    dict_total = {

        "dict_model" : best_fit_model.to_dict(), # SkyModel.from_dict(<>)

        "table_sed"  : flux_points.to_table(),   # FluxPoints.from_table(<>)

        "lightcurve" : {

            "run_number" : run_num,
            "t_start"    : time_min,
            "t_stop"     : time_max,
            "timedelta"  : delta_time,
            "flux"       : flux_lst1,
            "e_flux"     : flux_stat_err_lst1,

            "global" : {
                "e_min"               : e_lc_min,
                "e_max"               : e_lc_max,
                "n_off_regions"       : n_off_regions,
                "target_name"         : target_name,
                "crab_reference_flux" : flux_crab,
                "chi2"                : chi2_val,
                "pvalue"              : pvalue
            }
        }
    }


    # Saving the object
    with open(fname_dict, 'wb') as f:
        pickle.dump(dict_total, f, pickle.HIGHEST_PROTOCOL)    
    
    