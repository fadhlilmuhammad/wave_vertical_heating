# %%
import xarray as xr 
import numpy as np

import sys 
import os

from pathlib import Path

# %%
var=sys.argv[1]
folder = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/{var}/era5/")
folder_out = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/{var}/")

# %%
ds_clim = xr.open_mfdataset(os.path.join(folder,f"{var}_era5_oper_pl_merge*.nc"), chunks="auto")

# %%
daily_clim = ds_clim.groupby('time.dayofyear').mean('time')

# %%
ds_clim

# %%
daily_clim.to_netcdf(os.path.join(folder_out,f"{var}_era5_doyclim_1979-2020.nc"))

# %%



