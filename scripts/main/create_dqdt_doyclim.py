# %%
import xarray as xr
import numpy as np 
from pathlib import Path
import os 
import sys
import glob

import numpy as np
import os

def make_packing_encoding(da, dtype="int16", n_bits=16):
    """Compute scale_factor/add_offset to pack a float array into int16."""
    vmin = float(da.min().values)
    vmax = float(da.max().values)

    # reserve top value for _FillValue
    n = 2**n_bits
    scale_factor = (vmax - vmin) / (n - 2)
    add_offset = vmin + (n / 2) * scale_factor

    return {
        "dtype": dtype,
        "scale_factor": scale_factor,
        "add_offset": add_offset,
        "_FillValue": -32767,
        "zlib": True,
        "complevel": 4,
    }


# %%
folder = Path("/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/q/era5")

# %%
files = sorted(glob.glob(f"{folder}/*.nc"))

# %%
ds = xr.open_mfdataset(files)

# %%
da = ds.q

# %%
q = da.chunk({"time": -1})
dqdt = q.differentiate("time")


# %%
dqdt

# %%
daily_clim = dqdt.groupby('time.dayofyear').mean('time')

# %%
folder_out = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/dqdt/")

# %%
os.makedirs(folder_out, exist_ok=True)

# %%
# daily_clim.to_netcdf(os.path.join(folder_out,f"dqdt_era5_doyclim_1979-2020.nc"))

encoding = {daily_clim.name: make_packing_encoding(daily_clim)}

daily_clim.to_netcdf(
    os.path.join(folder_out, "dqdt_era5_doyclim_1979-2020.nc"),
    encoding=encoding,
)


# %%



