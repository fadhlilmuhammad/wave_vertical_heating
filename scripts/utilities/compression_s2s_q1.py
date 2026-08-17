
import xarray as xr 
import sys
from pathlib import Path
import os

# init=19810201
# ens="e01"
init=sys.argv[1]
ens=sys.argv[2]

folder = Path(f"/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/q1_vert_structure/{ens}/")
file = os.path.join(folder,f"da_q1q2_yanai_{init}_{ens}_remap.nc")
ds = xr.open_dataset(file)

comp = dict(zlib=True, complevel=7)
encoding = {var: comp for var in ds.data_vars}

output_file = os.path.join(folder,f"da_q1q2_yanai_{init}_{ens}_compressed.nc"
try:
    ds.to_netcdf(output_file, encoding=encoding)

    if output_file.exists():
        os.remove(file)
        print(f"Compression successful. Deleted: {file}")
    else:
        print(f"Output file not found. Kept: {file}")

except Exception as e:
    print(f"Compression failed: {e}")
    print(f"Kept original file: {file}")
    
# ds.to_netcdf(os.path.join(folder,f"da_q1q2_yanai_{init}_{ens}_compressed.nc"), encoding=encoding)

# os.remove(os.path.join(folder,f"da_moisture_budget_{init}_{ens}_remap.nc"))