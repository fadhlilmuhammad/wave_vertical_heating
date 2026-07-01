import xarray as xr 
import numpy as np 
import sys
import glob
import os
import xskillscore as xs
import pandas as pd
import random

from scipy.stats import pearsonr

from datetime import datetime,timedelta


from concurrent.futures import ProcessPoolExecutor
import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta



delta = timedelta(
    days=50,
    seconds=27,
    microseconds=10,
    milliseconds=29000,
    minutes=5,
    hours=8,
    weeks=2
)


np.set_printoptions(threshold=sys.maxsize)

#getting data
wave        = sys.argv[1]
mode        = sys.argv[2]
season      = sys.argv[3]
region_name = sys.argv[4]
bootsamples = 1000


# er MJJASO dry EIO-MC
#getting data
# wave        = 'mjo'
# season      = 'MJJASO'
# mode        = 'dry'
# region_name = 'EIO-MC'
# bootsamples = 2


# Determine control longitude based on wave type
if region_name == 'EIO-MC':
    if ((wave == 'mjo') or (wave == 'kelvin')):
        ctrl_lon = 95
    else:
        ctrl_lon = 133
elif region_name == 'NA-SIO':
    if ((wave == 'mjo') or (wave == 'kelvin')):
        ctrl_lon = 110
    else:
        ctrl_lon = 133
elif region_name == 'PO':
    if ((wave == 'mjo') or (wave == 'kelvin')):
        ctrl_lon = 170
    else:
        ctrl_lon = 250


if region_name == "EIO-MC":
    latN    = 15
    latS    = -15
    lonL    = 60
    lonR    = 150
elif region_name == "NA-SIO":
    latN    = -10
    latS    = -25
    lonL    = 95
    lonR    = 155
elif region_name == "PO":
    latN    = 15
    latS    = -15
    lonL    = 160
    lonR    = 270
else:
    print("region not supported, return back to EIO-MC")
    latN    = 15
    latS    = -15
    lonL    = 60
    lonR    = 150

var     = "rlut"
var2 = "olr"


if wave != 'tdsh' and wave != 'tdnh':
    file = f'local_{mode}_{wave}_phase_{season}_allclim.init_fors2s_{ctrl_lon}.nc'
    waveindex = wave
else:
    waveindex = wave
    file = f'local_{mode}_{waveindex}_phase_{season}_allclim.init_fors2s_{ctrl_lon}.nc'
    wave = 'td'
    
    print(f"{waveindex} change to {wave}")

home_raw="/home/563/fm6730/localrepo/wave_vertical_heating/data/raw"
home_drv="/home/563/fm6730/localrepo/wave_vertical_heating/data/derived"

folder_obs  = f'{home_raw}/observation/{var}'
# folder_s2s  = f'/scratch/v46/fm6730/dataout/access_s2/{var}/{wave}'
folder_s2s  = f'{home_drv}/access-s2_filtered/{var}/{wave}'
folder_init_phase = f'{home_drv}/local_wave_phase/'

# Load phase information
folder_phase = f'{home_drv}/local_wave_phase/'
# ensemble     = ['e01','e02','e03','ensmean_tlag']
ensemble     = ['e01','e02','e03']
print(folder_s2s)

folder_out = f"{home_drv}/acc/{var}/{wave}/{region_name}"
os.makedirs(folder_out, exist_ok=True)

#dictionary for each ensemble
wave_dict    = dict().fromkeys(ensemble, 0)

print(f"Control longitude: {ctrl_lon}\n")

    
    
phase_file = os.path.join(folder_phase, file)

print(f"Loading phase data from: {phase_file}")
if not os.path.exists(phase_file):
    print(f"ERROR: Phase file not found: {phase_file}")
    sys.exit(1)

phases = xr.open_dataset(phase_file)

time_active = pd.to_datetime(phases.time.values)
time_active_str = set(time_active.strftime('%Y%m%d'))

init_sets = [init for init in time_active_str if init != '19810101']
init_sets = sorted(init_sets)  # Sort for consistent ordering

print(f"Processing {len(init_sets)} initializations...\n")


# %%
print(init_sets)

# %%

# Main execution
wave_dict = {}

# Process ensemble members in parallel
with ProcessPoolExecutor(max_workers=4) as executor:
    futures = {
        executor.submit(
            process_ensemble_member, 
            e, 
            folder_s2s, 
            init_sets,
            wave, 
            var, 
            lonL, 
            lonR,
            latN,    # lat_N
            latS,    # lat_S
        ): e for e in ensemble if e != 'ensmean_tlag'
    }
    
    for future in futures:
        result = future.result()
        if result is not None:
            e, datasets_by_init = result
            
            # Concatenate all init dates into 'init' dimension
            ds_concat = xr.concat(datasets_by_init, dim='init')
            wave_dict[e] = ds_concat
            
            print(f"Processed {e}: shape = {ds_concat.sizes}")
            print(f"  init dates: {ds_concat.init.values}")
            print(f"  time_lag: {ds_concat.time_lag.values}")


# %%
with ProcessPoolExecutor(max_workers=4) as executor:
    futures = {
        executor.submit(
            process_truth,
            folder_obs, 
            folder_s2s, 
            init_sets,
            wave, 
            var, 
            lonL, 
            lonR,
            latN,    # lat_N
            latS,    # lat_S
        )
    }
    
    for future in futures:
        result = future.result()
        if result is not None:
            e, datasets_by_init = result
            
            # Concatenate all init dates into 'init' dimension
            ds_concat = xr.concat(datasets_by_init, dim='init')
            wave_dict[e] = ds_concat
            
            print(f"Processed {e}: shape = {ds_concat.sizes}")
            # print(f"  init dates: {ds_concat.init.values}")
            # print(f"  time_lag: {ds_concat.time_lag.values}")



# %%

weights = np.cos(np.deg2rad(wave_dict['true'].lat))
weights.name = "weights"

tmp_ensmean = (wave_dict['e01'] + wave_dict['e02'] + wave_dict['e03'])/3.
wave_dict['ensmean'] = tmp_ensmean.mean(dim='time_lag')
del(tmp_ensmean)

for lag in wave_dict['e01']['time_lag']:
    print(lag)
    wave_dict[f'ensmean_tlag_{lag.values}'] = (wave_dict['e01'].sel(time_lag=lag) + wave_dict['e02'].sel(time_lag=lag) + wave_dict['e03'].sel(time_lag=lag))/3

acc_dict = {}
ensemble_ensmean = ['ensmean', 'ensmean_tlag_-2', 'ensmean_tlag_-1', 'ensmean_tlag_0']
ensemble_expand = ensemble + ensemble_ensmean
print(ensemble_expand)

for e in ensemble_expand:
    acc_dict[e] = calculate_acc(wave_dict[e], wave_dict['true'], weights, latN=float(latN), latS=float(latS)) #OUTPUT IS R-correl

z_acc_dict = acc_dict.copy()
for keys in acc_dict.keys():
    z_acc_dict[keys] = np.arctanh(acc_dict[keys].copy()) #TO Z-transform


for keys in z_acc_dict.keys():
    z_acc_dict[keys] = z_acc_dict[keys].mean(dim='init')


# keys_to_run = [k for k in wave_dict.keys() if k != 'true']
# z_acc_boot = {}
# with ProcessPoolExecutor() as executor:
#     results = executor.map(
#         run_bootstrap_shuffle,
#         keys_to_run,
#         [wave_dict]*len(keys_to_run),
#         [weights]*len(keys_to_run),
#         [wave]*len(keys_to_run),
#         [latN]*len(keys_to_run),
#         [latS]*len(keys_to_run),
#         [bootsamples]*len(keys_to_run),
#     )

# for k, res in zip(keys_to_run, results):
#     z_acc_boot[k] = res

# %%

# z_ci = {}
# for keys in z_acc_boot:
#     z_ci[keys] = z_acc_boot[keys].quantile([0.025, 0.975], dim='n_bootstrap')


# %%

# z_ci_dict = z_ci.copy()

transformed_acc = z_acc_dict.copy()
# transformed_ci = z_ci_dict.copy()
for keys in z_acc_dict.keys():
    transformed_acc[keys] = np.tanh(z_acc_dict[keys].copy())
    transformed_acc[keys].attrs["description"] = "ACC result"
    transformed_acc[keys].attrs["nsamples"] = f"{len(init_sets)}"
    # transformed_acc[keys].to_netcdf(f"acc_{wave}_{mode}_{var}_{keys}.nc")
    
    # transformed_ci[keys] = np.tanh(z_ci_dict[keys].copy())
    # transformed_ci[keys].attrs["description"] = "Bootstrap ACC result"
    # transformed_acc[keys].attrs["nsamples"] = f"{len(init_sets)}"
    # transformed_ci[keys].to_netcdf(f"acc_ci_{wave}_{mode}_{var}_{keys}.nc")d

# %%
transformed_acc

# %%
transformed_acc_list = []
# transformed_ci_list = []
for key in transformed_acc.keys():
    transformed_acc_list.append(transformed_acc[key][wave].rename(key))
    # transformed_ci_list.append(transformed_ci[key].rename(key))

# %%
ds_acc = xr.Dataset({da.name: da for da in transformed_acc_list})
# ds_ci = xr.Dataset({da.name: da for da in transformed_ci_list})

ds_acc.attrs = {
    "description": "Anomaly Correlation Coefficient",
    "wave": str(waveindex),
    "season": season,
    "mode": mode,
    "variable": var,
    "region": region_name,
    "nsamples": f"{len(init_sets)}",
}

# ds_ci.attrs = {
#     "description": "Confidence Interval of ACC",
#     "wave": str(waveindex),
#     "season": season,
#     "mode": mode,
#     "variable": var,
#     "region": region_name,
#     "nsamples": f"{len(init_sets)}",
# }

# %%
ds_acc.to_netcdf(f"{folder_out}/acc_{waveindex}_{season}_{mode}_{var}_{region_name}_withoutcorrection.nc")
ds_ci.to_netcdf(f"{folder_out}/acc_ci_{waveindex}_{season}_{mode}_{var}_{region_name}_withoutcorrection.nc")


# %%



