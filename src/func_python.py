# %%
import xarray as xr 
import numpy as np 
import sys
import glob
import os
import xskillscore as xs
import pandas as pd
import random
import re

from scipy.stats import pearsonr

from datetime import datetime,timedelta


from concurrent.futures import ProcessPoolExecutor
import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta

def extract_date(filename):
    match = re.search(r'\d{8}', filename)
    if not match:
        raise ValueError(f"No date found in {filename}")
    return datetime.strptime(match.group(), "%Y%m%d")

def get_region(region_name,wave):
    
    #CTRL longitude
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
            
    #lat-lon
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
        
    return ctrl_lon,latN,latS,lonL,lonR

def reverse_latitude_if_ascending(ds):
    """
    Reverse latitude values in an xarray Dataset or DataArray if they are in ascending order.
    
    Parameters:
    - ds: xarray.Dataset or xarray.DataArray

    Returns:
    - xarray.Dataset or xarray.DataArray with reversed latitude if needed
    """
    # Check if the input is a Dataset or DataArray and get the latitude variable
    if isinstance(ds, xr.Dataset):
        if 'lat' in ds.coords:
            lat = ds.coords['lat']
        else:
            raise ValueError("Latitude coordinate 'lat' not found in the Dataset.")
    elif isinstance(ds, xr.DataArray):
        if 'lat' in ds.dims:
            lat = ds['lat']
        else:
            raise ValueError("Latitude dimension 'lat' not found in the DataArray.")
    else:
        raise TypeError("Input must be an xarray.Dataset or xarray.DataArray.")

    # Check if latitudes are in ascending order
    if lat.values[0] < lat.values[-1]:
        # print("Latitudes are in ascending order. Reversing...")
        
        # Reverse the latitude and all associated data
        ds = ds.sortby('lat', ascending=False)
    
    return ds


def process_ensemble_member(e, folder_s2s, dates, wave, var, lonL, lonR, lat_N, lat_S):
    """Process a single ensemble member"""
    
    if e == 'ensmean_tlag':
        return None
    
    try:
        # Load all files matching pattern
        file_list = f"{folder_s2s}/{wave}_tlag_da_{var}_*_{e}_remap_padded.nc"
        
        files = [
            f"{folder_s2s}/{wave}_tlag_da_{var}_{date}_{e}_remap_padded.nc"
            for date in dates
        ]

        # print(expected_files)
        # print(files)
        
        datasets_by_init = []
        
        for file_path in files:
            # Extract date from filename (yyyymmdd)
            # Assuming format: wave_tlag_da_var_YYYYMMDD_e_remap_padded.nc
            # print(file_path)
            filename = file_path
            # parts = filename.split('_')
            # date_str = parts[5]  # Extract YYYYMMDD
            date_str = extract_date(filename).strftime("%Y%m%d")
            # print(date_str)
            init_date = datetime.strptime(date_str, '%Y%m%d')
            
            # Load dataset
            ds = xr.open_dataset(str(file_path))
            # print(ds)
            # Reverse latitude if ascending
            ds = reverse_latitude_if_ascending(ds)

            # Set dates
            date_choose = [(datetime.strptime(date_str, "%Y%m%d") + timedelta(days=i)).strftime("%Y%m%d") for i in range(-1, 71)]
            
            # Select spatial bounds
            ds = ds.sel(
                lat=slice(lat_N, lat_S), 
                lon=slice(float(lonL), float(lonR)),
                time=slice(date_choose[0],date_choose[-1]),
            )

            
            # Rewrite time coordinate as integers from -1 to 70 (days from init)
            # Assuming time_lag dimension exists, or rename time to time_lag
            # if 'time_lag' in ds.dims:
            #     time_dim = 'time_lag'
            # else:
            #     time_dim = 'time'
            #     ds = ds.rename({time_dim: 'time_lag'})
            
            # Create integer time coordinates from -1 to 70
            n_times = len(ds['time'])
            time_int = np.arange(-1, n_times - 1)
            ds = ds.assign_coords({'time': time_int})
            
            # Add init date as coordinate
            ds = ds.assign_coords(init=init_date)
            
            datasets_by_init.append(ds)
        
        return (e, datasets_by_init)
    
    except Exception as ex:
        print(f"Error processing {e}: {ex}")
        return None


def process_truth(folder_obs, folder_s2s, dates, wave, var, lonL, lonR, lat_N, lat_S):
    """Process truth"""

    # Load all files matching pattern
    file_list = f"{folder_s2s}/{wave}_tlag_da_{var}_*_e01_remap_padded.nc"
    
    files = [
        f"{folder_s2s}/{wave}_tlag_da_{var}_{date}_e01_remap_padded.nc"
        for date in dates
    ]

    file_obs = f"{var}.{wave}.30_truth_for_S2S.nc"
    
    # print(expected_files)
    # print(files)
    
    datasets_by_init = []
    
    for file_path in files:
        # Extract date from filename (yyyymmdd)
        # Assuming format: wave_tlag_da_var_YYYYMMDD_e_remap_padded.nc
        # print(file_path)
        filename = file_path
        date_str = extract_date(filename).strftime("%Y%m%d")
        # print(date_str)
        init_date = datetime.strptime(date_str, '%Y%m%d')
    
        # Load dataset
        ds = xr.open_dataset(f"{folder_obs}/{file_obs}")
        
        # Reverse latitude if ascending
        ds = reverse_latitude_if_ascending(ds)

        # Set dates
        date_choose = [(datetime.strptime(date_str, "%Y%m%d") + timedelta(days=i)).strftime("%Y%m%d") for i in range(-1, 71)]
        
        # print(date_choose[0])
        # print(date_choose[-1])
        # Select spatial bounds
        ds = ds.sel(
            lat=slice(lat_N, lat_S), 
            lon=slice(float(lonL), float(lonR)),
            time=slice(date_choose[0],date_choose[-1]),
        )
        
        # Create integer time coordinates from -1 to 70
        n_times = len(ds['time'])
        time_int = np.arange(-1, n_times - 1)
        ds = ds.assign_coords({'time': time_int})
        
        # Add init date as coordinate
        ds = ds.assign_coords(init=init_date)
        
        datasets_by_init.append(ds)
        # print(datasets_by_init)
    
    return ("true", datasets_by_init)


def calculate_acc(x1,x2, weights, latN=None, latS=None, ):
    #x1 and x2 is DataArray with 3Dimensions
    # print(x1)
    # print(x2)


    weights_use = weights.copy()
    weights_use = weights.sel(lat=slice(latN,latS))
    weights_use.name = "weights"

    a = x1.sel(lat=slice(latN,latS)).copy()
    b = x2.sel(lat=slice(latN,latS)).copy()

    aweighted = (a*weights_use).copy()
    bweighted = (b*weights_use).copy()

    amean = aweighted.mean(dim=["lat", "lon"])
    bmean = bweighted.mean(dim=["lat", "lon"])

    a_prod_b = weights_use * (a-amean)*(b-bmean)
    abcov = a_prod_b.mean(dim=["lat", "lon"])

    aanom2 = weights_use * (a-amean)**2
    banom2 = weights_use * (b-bmean)**2

    aanom2_sum = aanom2.mean(dim=["lat", "lon"])
    banom2_sum = banom2.mean(dim=["lat", "lon"])

    acc = abcov / (np.sqrt(aanom2_sum)*np.sqrt(banom2_sum))
    
    # acc = xr.corr(a,b,dim=["lat","lon"], weights=weights_use)

    return acc


# %%
def bootstrap_acc(acc, n_bootstrap=1000, dim_use='init'):

    n_sample = len(acc[dim_use])
    time_size = len(acc['time'])
    
    acc_bootstrap = np.zeros((n_bootstrap, time_size))
    # print(f"acc_boot:{acc_bootstrap.shape}")
    
    for i in range(n_bootstrap):
        # print(f"bootstrap-round: {i+1}")
        
        # Resample with replacement along initialization dimension
        boot_idx = np.random.choice(n_sample, n_sample, replace=True)

        boot_chosen = acc.isel(init=boot_idx)
        # print(boot_chosen)
        boot_mean = boot_chosen.mean(dim=dim_use).to_numpy()

        
        # Put across initializations (now bootstrap dimension)
        acc_bootstrap[i,:] = boot_mean

        # if (i + 1) % 100 == 0:
        #     print(f"    Completed {i + 1}/{n_bootstrap} first bootstrap iterations")

    return acc_bootstrap

def shuffle_window_blocks(da, dim="time", window=10, seed=None):
    # RNG
    rng = np.random.default_rng(seed)

    # create storage
    n = da.sizes[dim]
    blocks = []
    
    # Create full blocks
    for start in range(0, n, window):
        end = min(start + window, n)  # <-- residual included in a window too, return the end value as the minimum of n or start+window
        block = da.isel({dim: slice(start, end)}) # create blocks
        blocks.append(block) # append blocks
    
    # Shuffle all blocks (including residual)
    shuffled_order = rng.permutation(len(blocks)) # shuffle blocks index
    shuffled_blocks = [blocks[i] for i in shuffled_order] # shuffle blocks
    
    # Concatenate back
    return xr.concat(shuffled_blocks, dim=dim, join="override")

def shuffle_window_blocks_paired(da1, da2, dim="time", window=10, seed=None):

    rng = np.random.default_rng(seed)

    n = da1.sizes[dim]
    blocks1 = []
    blocks2 = []

    # Create blocks
    for start in range(0, n, window):
        end = min(start + window, n)

        blocks1.append(da1.isel({dim: slice(start, end)}))
        blocks2.append(da2.isel({dim: slice(start, end)}))

    # Same shuffled order for both
    shuffled_order = rng.permutation(len(blocks1))

    shuffled1 = [blocks1[i] for i in shuffled_order]
    shuffled2 = [blocks2[i] for i in shuffled_order]

    # Concatenate
    da1_shuffled = xr.concat(shuffled1, dim=dim, join="override")
    da2_shuffled = xr.concat(shuffled2, dim=dim, join="override")

    return da1_shuffled, da2_shuffled

def get_dim_size(obj, dim_name):
    if dim_name in obj.dims:
        return obj.sizes[dim_name]
    else:
        return None

# def get_dim_size()
#     if dim_name in obj.coords:
#         return obj.coords[dim_name].size
#     else:
#         return None
    # raise ValueError(f"{dim_name} not found in dims or coords")

    
def bootstrap_random_shuffle(mod, obs, weights, wave, tlag=None, latN=None, latS=None, n_bootstrap=1000, dim_use='init'):

    n_sample = len(mod[dim_use])
    time_size = len(mod['time'])
    
    # print(mod)
    tlag_size = get_dim_size(mod, tlag)
    # print(tlag_size)
    
    if tlag_size is None:
        print(f"{dim_use} not found, skipping lag dimension")
        acc_bootstrap = np.zeros((n_bootstrap, time_size))
        # fallback logic here
    else :
        acc_bootstrap = np.zeros((n_bootstrap, time_size, tlag_size))
    # ...

    print(f"acc_boot:{acc_bootstrap.shape}")
    
    
    for i in range(n_bootstrap):
        # print(f"bootstrap-round: {i+1}")
        
        # Resample with replacement along initialization dimension
        boot_idx = np.random.choice(n_sample, n_sample, replace=True)

        print(boot_idx)
        mod_chosen = mod.copy().isel(init=boot_idx)
        obs_chosen = obs.copy().isel(init=boot_idx)
        
        # print(mod_chosen)
        # print(obs_chosen)
        
        # mod_shuffled, obs_shuffled = shuffle_window_blocks_paired(mod_chosen, obs_chosen)
        # mod_shuffled = shuffle_window_blocks(mod_chosen)

        mod_shuffled = mod_chosen.assign_coords(init=np.arange(mod_chosen.sizes[dim_use]))
        obs_shuffled = obs_chosen.assign_coords(init=np.arange(mod_chosen.sizes[dim_use]))
        # mod_shuffled = mod_shuffled.assign_coords(init=np.arange(mod_shuffled.sizes[dim_use])) # reassign 'dim' to make it in order as obs
        # obs_shuffled = obs_shuffled.assign_coords(init=np.arange(obs_shuffled.sizes[dim_use]))
        # print(obs_shuffled)
        
        # print(mod_chosen)
        # print(mod_shuffled)
        
        acc = calculate_acc(mod_shuffled, obs_shuffled, weights, latN, latS)
        # print(boot_chosen)
        # (f"acc:{acc[wave].shape}")
        z_acc = np.arctanh(acc[wave].copy())
        # print(f"z_acc:{z_acc.shape}")
        
        boot_mean = z_acc.mean(dim=dim_use).to_numpy()
        # print(f"z_acc_mean:{boot_mean.shape}")

        # Put across initializations (now bootstrap dimension)
        acc_bootstrap[i,...] = boot_mean

        # if (i + 1) % 100 == 0:
        #     print(f"    Completed {i + 1}/{n_bootstrap} first bootstrap iterations")

    if tlag_size is None:
        print(f"{dim_use} not found, skipping lag dimension")
        print("ACC (Fisher-Z transform), ready to average. Transform to r first for final value")
        ds_acc_z = xr.DataArray(data=acc_bootstrap,
            dims=["n_bootstrap", "time"],
            coords=dict(
                n_bootstrap = (['n_bootstrap'], np.arange(n_bootstrap)),
                time= mod.time,
            ),
            attrs=dict(
                description="ACC (Fisher-Z transform), ready to average. Transform to r first for final value",
            ),
        )
    # fallback logic here
    else :
        print("ACC (Fisher-Z transform), ready to average. Transform to r first for final value")
        ds_acc_z = xr.DataArray(data=acc_bootstrap,
            dims=["n_bootstrap", "time", "time_lag"],
            coords=dict(
                n_bootstrap = (['n_bootstrap'], np.arange(n_bootstrap)),
                time = mod.time,
                time_lag = mod.time_lag,
            ),
            attrs=dict(
                description="ACC (Fisher-Z transform), ready to average. Transform to r first for final value",
            ),
        )

                        
    return ds_acc_z

def run_bootstrap_shuffle(keys, wave_dict, weights, wave, latN, latS, n_bootstrap):
    return bootstrap_random_shuffle(
        wave_dict[keys],
        wave_dict['true'],
        weights,
        wave,
        'time_lag',
        float(latN),
        latS=float(latS),
        n_bootstrap=n_bootstrap,
    )

# %%
# Function to subtract 12 hours if time is exactly 12:00, otherwise keep the same
def adjust_time(time):
    if time.astype('datetime64[m]').astype(str).endswith('12:00'):
        return time - np.timedelta64(12, 'h')
    elif time.astype('datetime64[m]').astype(str).endswith('11:00'):
        return time - np.timedelta64(11, 'h')
    else:
        return time

# %%


def bootstrap_random_shuffle_index(xobs, yobs, xmod, ymod, tlag=None, n_bootstrap=1000, process='bivariate_correlation', dim_use='init'):

    n_sample = len(xmod[dim_use])
    time_size = len(xmod['time'])
    
    # print(mod)
    tlag_size = get_dim_size(xmod, tlag)
    
    if tlag_size is None:
        print(f"tlag not found, skipping lag dimension")
        var_bootstrap = np.zeros((n_bootstrap, time_size))
        # fallback logic here
    else :
        var_bootstrap = np.zeros((n_bootstrap, time_size, tlag_size))
    # ...

    print(f"acc_boot:{var_bootstrap.shape}")
    
    
    for i in range(n_bootstrap):
        # print(f"bootstrap-round: {i+1}")
        
        # Resample with replacement along initialization dimension
        boot_idx = np.random.choice(n_sample, n_sample, replace=True)
        
        # xobs_shuffled, yobs_shuffled, xmod_shuffled, ymod_shuffled = shuffle_window_blocks_paired_index(xobs, yobs, xmod, ymod)
        xmod_shuffled = xmod.isel(init=boot_idx)
        ymod_shuffled = ymod.isel(init=boot_idx)
        # xobs_shuffled, yobs_shuffled, xmod_shuffled, ymod_shuffled = shuffle_window_blocks_paired_index(xobs, yobs, xmod, ymod)
        
        xobs_shuffled = xobs.isel(init=boot_idx)
        yobs_shuffled = yobs.isel(init=boot_idx)
        
        xobs_shuffled = xobs_shuffled.assign_coords(init=np.arange(xobs_shuffled.sizes[dim_use]))
        yobs_shuffled = yobs_shuffled.assign_coords(init=np.arange(yobs_shuffled.sizes[dim_use]))
        xmod_shuffled = xmod_shuffled.assign_coords(init=np.arange(xmod_shuffled.sizes[dim_use]))
        ymod_shuffled = ymod_shuffled.assign_coords(init=np.arange(ymod_shuffled.sizes[dim_use]))
        
        # print(boot_chosen)
        # print("BOOT_XOBS")
        # (f"{xobs_shuffled.shape}")
        
        if process == 'bivariate_correlation':
            bcor = bivariate_correlation(xobs_shuffled,yobs_shuffled,xmod_shuffled,ymod_shuffled)
            # vartmp = np.arctanh(bcor.copy())
            vartmp = bcor
            del(bcor)
        elif process == 'rmse':
            vartmp = rmse(xobs_shuffled,yobs_shuffled,xmod_shuffled,ymod_shuffled)
        elif process == 'amplitude_error': 
            vartmp = amplitude_error(xobs_shuffled,yobs_shuffled,xmod_shuffled,ymod_shuffled)
        elif process == 'phase_error':
            vartmp = phase_error(xobs_shuffled,yobs_shuffled,xmod_shuffled,ymod_shuffled)
        else:
            print('WARNING! PROCESS NOT SUPPORTED. Only "bivariate_correlation", "rmse", "amplitude_error", "phase_error" are supported')
        # print(f"z_acc:{z_acc.shape}")
        
        
        # print(vartmp)
        boot_mean = vartmp.to_numpy()
        # print(f"z_acc_mean:{boot_mean.shape}")

        # Put across initializations (now bootstrap dimension)
        var_bootstrap[i,...] = boot_mean

        # if (i + 1) % 100 == 0:
        #     print(f"    Completed {i + 1}/{n_bootstrap} first bootstrap iterations")

    if tlag_size is None:
        if process == 'bivariate_correlation':
            # print("BCOR (Fisher-Z transform), ready to average. Transform to r first for final value")
            # desc = f"Processed {process}, in Z-transform"
            print("Already in r-correlation")
            desc = f"Processed {process}, in R-correlation"
        else:
            print(f"Process {process} is finished...")
            desc = f"Processed {process}"
            
        ds_acc_z = xr.DataArray(data=var_bootstrap,
            dims=["n_bootstrap", "time"],
            coords=dict(
                n_bootstrap = (['n_bootstrap'], np.arange(n_bootstrap)),
                time= xmod.time,
            ),
            attrs=dict(
                description=desc,
            ),
        )
    # fallback logic here
    else :
        if process == 'bivariate_correlation':
            # print("BCOR (Fisher-Z transform), ready to average. Transform to r first for final value")
            print("Already in r-correlation")
            desc = f"Processed {process}, in R-correlation"
        else:
            print(f"Process {process} is finished...")
            desc = f"Processed {process}"
            
        ds_acc_z = xr.DataArray(data=var_bootstrap,
            dims=["n_bootstrap", "time", "time_lag"],
            coords=dict(
                n_bootstrap = (['n_bootstrap'], np.arange(n_bootstrap)),
                time = xmod.time,
                time_lag = xmod.time_lag,
            ),
            attrs=dict(
                description=desc,
            ),
        )

                        
    return ds_acc_z


def bootstrap_random_shuffle_1D(da, n_bootstrap=400, random_state=10, dim='init'):
    np.random.seed(random_state)
    nsize = da.sizes[dim]

    iboot = np.random.choice(nsize, size=(n_bootstrap, nsize), replace=True)
    boot_samples = []

    for i in range(n_bootstrap):
        boot_sample = da.isel({dim: iboot[i]})
        boot_sample = boot_sample.assign_coords(init=np.arange(0,nsize))
        boot_samples.append(boot_sample)

    # print(boot_samples)

    return xr.concat(boot_samples, dim='bootstrap')
            
            
def phase_error(xobs, yobs, xmod, ymod, dim='init'):
    "xobs is dataarray (init, leadtime) and yobs is (init, leadtime)"
    "xmod is dataarray (init, leadtime) and ymod is (init, leadtime)"

    nom = xobs*ymod - yobs*xmod
    denom = xobs*xmod + yobs*ymod

    tmp = np.rad2deg(np.arctan2(nom, denom))

    err_phase = tmp.mean(dim=dim)

    # print(err_phase.shape)
    return err_phase

def bivariate_correlation(xobs, yobs, xmod, ymod, dim='init'):

    tmp_nom = xobs*xmod + yobs*ymod
    tmp_denom_1 = xobs**2 + yobs**2
    tmp_denom_2 = xmod**2 + ymod**2

    nom = tmp_nom.sum(dim=dim)
    denom_1 = tmp_denom_1.sum(dim=dim)
    denom_2 = tmp_denom_2.sum(dim=dim)

    corr = nom / (np.sqrt(denom_1)*np.sqrt(denom_2))
    # print(corr.shape)
    return corr

def amplitude_error(xobs, yobs, xmod, ymod, dim='init'):
    "xobs is dataarray (init, leadtime) and yobs is (init, leadtime)"
    "xmod is dataarray (init, leadtime) and ymod is (init, leadtime)"

    amp_obs = np.sqrt(xobs**2 + yobs**2)
    amp_mod = np.sqrt(xmod**2 + ymod**2)

    err_amplitude = amp_mod - amp_obs
    # print(err_amplitude.shape)
    return err_amplitude.mean(dim=dim)

def rmse(xobs, yobs, xmod, ymod, dim='init'):
    "xobs is dataarray (init, leadtime) and yobs is (init, leadtime)"
    "xmod is dataarray (init, leadtime) and ymod is (init, leadtime)"

    nom = (xobs - xmod)**2 + (yobs - ymod)**2

    nommean = nom.mean(dim=dim)
    rmse = np.sqrt(nommean)
    # print(rmse.shape)
    return rmse