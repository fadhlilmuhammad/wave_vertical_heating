

load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/utils.ncl"

; init  = toint(getenv("init"))
; ens   = getenv("ens")


folder = "/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/"
start_time = systemfunc("date +%s")
print("Script started at: " + systemfunc("date"))

fili_q =  folder+"/q/q_era5_doyclim_1979-2020.nc"
fili_u =  folder+"/u/u_era5_doyclim_1979-2020.nc"
fili_v =  folder+"/v/v_era5_doyclim_1979-2020.nc"
fili_w =  folder+"/w/w_era5_doyclim_1979-2020.nc"

latN = 40
latS = -latN

fq = addfile(fili_q, "r")
fu = addfile(fili_u, "r")
fv = addfile(fili_v, "r")
fw = addfile(fili_w, "r")


ptop = 100.
pbot = 1000.

g    = 9.81          ; m/s2

q = short2flt(fq->q(:,{pbot:ptop},{latS:latN},:))
u = short2flt(fu->u(:,{pbot:ptop},{latS:latN},:))
v = short2flt(fv->v(:,{pbot:ptop},{latS:latN},:))
w = short2flt(fw->w(:,{pbot:ptop},{latS:latN},:))

lon = u&lon 
lat = u&lat
time = u&dayofyear
ntime           = dimsizes(time)

print("Reading variables...")

printVarSummary(q)
printVarSummary(u)
printVarSummary(v)

ptop = ptop*100
ptop@units = "Pa"
pbot = pbot * 100.0

; lev = q&z1_p_level                      ; hPa
plev = q&level                      ; hPa
plev = plev*100                     ; [100000,...,30000] Pa [kg/(m-s2)]
plev@units = "Pa"

dp 		= dpres_plevel_Wrap(plev, pbot, ptop, 0)
DP    = conform(u,dp,1)

q    = q
q@units = "kg/kg"

ps = DP

dpg  = DP/g
dpg@long_name = "Layer Mass Weighting"
dpg@units     = "kg/m2"      ; dp/g, Pa/(m s-2), reduce to kg m-2


;gridfill if any of the variables contain missing values;
guess_type = 1                ; use zonal means
is_cyclic  = True             ; cyclic [global]
nscan      = 500             ; usually much less than this
epsx       = 1.e-2            ; variable dependent
relc       = 0.6              ; relaxation coefficient
opt        = 0                ; not used

if any(ismissing(u)) then
    print("u contains missing values")
    poisson_grid_fill(u, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("u filled")
end if

if any(ismissing(v)) then
    print("v contains missing values")
    poisson_grid_fill(v, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("v filled")
end if

if any(ismissing(q)) then
    print("q contains missing values")
    poisson_grid_fill(q, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("q filled")
end if

if any(ismissing(w)) then
    print("w contains missing values")
    poisson_grid_fill(w, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("w filled")
end if

print("Calculating advection of moisture...")

rad = get_d2r("float")
re = 6.371e+06
; g=9.81

print("Calculating temporal tendency of moisture...")

  q1 = q
  q1(1:ntime-1,:,:,:) = (/q(0:ntime-2,:,:,:)/)

  dq_dt = (/(q-q1)/86400.0/)
  copy_VarCoords(q, dq_dt)

dq_dt@long_name = "Temporal tendency of moisture"
dq_dt@units = "kg/kg s**-1"

printVarSummary(dq_dt)
printMinMax(dq_dt, False)

delete(q1)

dq_dx = center_finite_diff_n(q, lon * rad, True, 0, 3) / conform(q, re * cos(lat * rad), 2)
dq_dy = center_finite_diff_n(q, lat * rad, False, 0, 2) / re
dq_dp = center_finite_diff_n(q, plev, False, 0, 1) 

adv_q_zonal 		= -u * dq_dx
adv_q_meridional 	= -v * dq_dy
adv_q_vertical 		= -w * dq_dp

; adv_q = adv_q_zonal + adv_q_meridional

; copy_VarCoords(q, adv_q)
; adv_q@long_name = "Advection of moisture"
; adv_q@units = "kg/kg s**-1"

print("Calculating convergence of moisture...")

du_dx = center_finite_diff_n(u, lon * rad, True, 0, 3) / conform(q, re * cos(lat * rad), 2)
dv_dy = center_finite_diff_n(v, lat * rad, False, 0, 2) / re
dw_dp = center_finite_diff_n(w, plev, False, 0, 1)

conv_q_zonal = -q * du_dx
conv_q_meridional = -q * dv_dy
conv_q_vertical = -q * dw_dp

; conv_q = conv_q_zonal + conv_q_meridional

print("Calculating moisture loss (gain) due to the condensational heating")
q2 = dq_dt + (u * dq_dx) + (v * dq_dy) + (w * dq_dp)


;+++++++++++++++++++++++++++++++++++++++++++++++++++++

print("Calculating horizontal convergence of moisture flux...")

; mfc = adv_q + conv_q
; copy_VarCoords(q, mfc)

; mfc@long_name = "Convergence of moisture flux"
; mfc@units = "kg/kg s**-1"

; printVarSummary(mfc)
; printMinMax(mfc, False)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++

; print("Calculating vertical convergence of moisture flux...")

; mfc_vertical = adv_q_vertical + conv_q_vertical
; copy_VarCoords(q, mfc_vertical)

; mfc_vertical@long_name = "Vertical convergence of moisture flux"
; mfc_vertical@units = "kg/kg s**-1"

; printVarSummary(mfc_vertical)
; printMinMax(mfc_vertical, False)

;+++++++++++++++++++++++++++++++++++++++++++++++++++++
; print("Calculating vertical integrals...")

; psfc = max(plev)
; ptop = min(plev)
; dp = dpres_plevel(plev, psfc, ptop, 0)
; printVarSummary(dp)
; print(dp)
    
;+++++++++++++++++++++++++++++++++++++++++++++++++++++

; copy_VarCoords(q, conv_q)
; copy_VarCoords(q, adv_q)
copy_VarCoords(q, dq_dt)
copy_VarCoords(q, adv_q_vertical)
copy_VarCoords(q, adv_q_meridional)
copy_VarCoords(q, adv_q_zonal)
copy_VarCoords(q, conv_q_zonal)
copy_VarCoords(q, conv_q_meridional)
copy_VarCoords(q, conv_q_vertical)
; copy_VarCoords(q, mfc)
; copy_VarCoords(q, mfc_vertical)
copy_VarCoords(q, q2)

coords = addfile("/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/w/w_era5_doyclim_1979-2020.nc", "r")
coords_lat = coords->lat({latS:latN})
coords_lon = coords->lon

fo_q = area_conserve_remap_Wrap (q&lon,q&lat,q, coords_lon,coords_lat({latS:latN}), False)
fo_dq_dt = area_conserve_remap_Wrap (dq_dt&lon,dq_dt&lat,dq_dt, coords_lon,coords_lat({latS:latN}), False) 

; fo_adv_q = area_conserve_remap_Wrap (adv_q&lon,adv_q&lat,adv_q, coords_lon,coords_lat({latS:latN}), False)
fo_adv_q_meridional = area_conserve_remap_Wrap (adv_q_meridional&lon,adv_q_meridional&lat,adv_q_meridional, coords_lon,coords_lat({latS:latN}), False)
fo_adv_q_zonal = area_conserve_remap_Wrap (adv_q_zonal&lon,adv_q_zonal&lat,adv_q_zonal, coords_lon,coords_lat({latS:latN}), False)
fo_adv_q_vertical = area_conserve_remap_Wrap (adv_q_vertical&lon,adv_q_vertical&lat,adv_q_vertical, coords_lon,coords_lat({latS:latN}), False)

; fo_conv_q = area_conserve_remap_Wrap (conv_q&lon,conv_q&lat,conv_q, coords_lon,coords_lat({latS:latN}), False)
fo_conv_q_zonal = area_conserve_remap_Wrap (conv_q_zonal&lon,conv_q_zonal&lat,conv_q_zonal, coords_lon,coords_lat({latS:latN}), False)
fo_conv_q_meridional = area_conserve_remap_Wrap (conv_q_meridional&lon,conv_q_meridional&lat,conv_q_meridional, coords_lon,coords_lat({latS:latN}), False)
fo_conv_q_vertical = area_conserve_remap_Wrap (conv_q_vertical&lon,conv_q_vertical&lat,conv_q_vertical, coords_lon,coords_lat({latS:latN}), False)


; fo_mfc = area_conserve_remap_Wrap (mfc&lon,mfc&lat,mfc, coords_lon,coords_lat({latS:latN}), False)
; fo_mfc_vertical = area_conserve_remap_Wrap (mfc_vertical&lon,mfc_vertical&lat,mfc_vertical, coords_lon,coords_lat({latS:latN}), False)
fo_q2 = area_conserve_remap_Wrap (q2&lon,q2&lat,q2, coords_lon,coords_lat({latS:latN}), False)
printVarSummary(fo_q)
printVarSummary(fo_dq_dt)

; printVarSummary(fo_adv_q)
; printVarSummary(fo_conv_q)
print("Writing output files...")

setfileoption("nc","Format","NetCDF4")

folder1 = "/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/moisture_budget/"
system("mkdir -p "+folder1)

filo1 = folder1+"moisture_budget_era5_doyclim_1979-2020.nc"
system("rm -f "+filo1)

fout = addfile(filo1, "c")

;----------------------------------------------------------
; 1. Dimensions
;----------------------------------------------------------
dimNames = (/"dayofyear","level","lat","lon"/)
dimSizes = (/ 366,        17,     33,   144  /)
dimUnlim = (/ False,      False,  False,False/)
filedimdef(fout, dimNames, dimSizes, dimUnlim)

;----------------------------------------------------------
; 2. Coordinate variables
;----------------------------------------------------------
filevardef(fout, "dayofyear", typeof(fo_dq_dt&dayofyear), "dayofyear")
filevardef(fout, "level",     typeof(fo_dq_dt&level),     "level")
filevardef(fout, "lat",       typeof(fo_dq_dt&lat),       "lat")
filevardef(fout, "lon",       typeof(fo_dq_dt&lon),       "lon")

filevarattdef(fout, "dayofyear", fo_dq_dt&dayofyear)
filevarattdef(fout, "level",     fo_dq_dt&level)
filevarattdef(fout, "lat",       fo_dq_dt&lat)
filevarattdef(fout, "lon",       fo_dq_dt&lon)

fout->dayofyear = (/fo_dq_dt&dayofyear/)
fout->level     = (/fo_dq_dt&level/)
fout->lat       = (/fo_dq_dt&lat/)
fout->lon       = (/fo_dq_dt&lon/)

;----------------------------------------------------------
; 3. Data variables — define + attr each explicitly
;----------------------------------------------------------
; filevardef(fout, "dq_dt", typeof(fo_dq_dt), (/"dayofyear","level","lat","lon"/))
; filevarattdef(fout, "dq_dt", fo_dq_dt)

filevardef(fout, "q", typeof(fo_q), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "q", fo_q)

filevardef(fout, "conv_q_meridional", typeof(fo_conv_q_meridional), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "conv_q_meridional", fo_conv_q_meridional)

filevardef(fout, "conv_q_zonal", typeof(fo_conv_q_zonal), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "conv_q_zonal", fo_conv_q_zonal)

filevardef(fout, "conv_q_vertical", typeof(fo_conv_q_vertical), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "conv_q_vertical", fo_conv_q_vertical)

filevardef(fout, "adv_q_zonal", typeof(fo_adv_q_zonal), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "adv_q_zonal", fo_adv_q_zonal)

filevardef(fout, "adv_q_meridional", typeof(fo_adv_q_meridional), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "adv_q_meridional", fo_adv_q_meridional)

filevardef(fout, "adv_q_vertical", typeof(fo_adv_q_vertical), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "adv_q_vertical", fo_adv_q_vertical)

filevardef(fout, "q2", typeof(fo_q2), (/"dayofyear","level","lat","lon"/))
filevarattdef(fout, "q2", fo_q2)

;----------------------------------------------------------
; 4. Write actual data
;----------------------------------------------------------
; fout->dq_dt              = (/fo_dq_dt/)
fout->q                  = (/fo_q/)
fout->conv_q_meridional  = (/fo_conv_q_meridional/)
fout->conv_q_zonal       = (/fo_conv_q_zonal/)
fout->conv_q_vertical    = (/fo_conv_q_vertical/)
fout->adv_q_zonal        = (/fo_adv_q_zonal/)
fout->adv_q_meridional   = (/fo_adv_q_meridional/)
fout->adv_q_vertical     = (/fo_adv_q_vertical/)
fout->q2                 = (/fo_q2/)

print("checking...")
printVarSummary(fout)

fout@title = "Moisture Budget"
fout@history = "Created by NCL script on " + systemfunc("date")
fout@Conventions = "None"
fout@source_file = fili_q + " and " + fili_u + " and " + fili_v + " and " + fili_w
fout@creation_date = systemfunc("date")
fout@notes = "The dq/dt is on the other files."
fout@institution = "ARC W21C"
fout@source = "ERA5"
fout@Conventions = "CF-1.6"



