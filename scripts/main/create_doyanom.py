# %%
import xarray as xr 
import numpy as np

import sys 
import os

from pathlib import Path

# %%
var=sys.argv[1]
varfolder=sys.argv[2]
ens=sys.argv[3]
# year=sys.argv[3]

folder = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/{varfolder}/")
folder_out = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_anom/access_s2/{varfolder}/")
folder_in = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/q1_vert_structure/{ens}/")

os.makedirs(folder_out, exist_ok=True)

# %%
ds_clim = xr.open_dataset(os.path.join(folder,f"{varfolder}_era5_doyclim_1979-2020.nc"), chunks="auto")
ds_data = xr.open_dataset()

# %%
da_clim = ds_clim[var]



# %%
daily_clim.to_netcdf(os.path.join(folder_out,f"{var}_era5_doyclim_1979-2020.nc"))

# %%



