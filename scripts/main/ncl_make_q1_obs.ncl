

load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/utils.ncl"

; init = 19980201
; ens = "e01"

; year=2000
year=year

start_time = systemfunc("date +%s")
print("Script started at: " + systemfunc("date"))

folder_in = "/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/era5/daily"
; fili_q = "/g/data/ux62/access-s2/hindcast/raw_model/atmos/hus/daily/"+ens+"/da_hus_"+init+"_"+ens+".nc"
fili_w =  folder_in+"/w/w_era5_oper_pl_merge_1deg_daily_"+year+".nc"
fili_u =  folder_in+"/u/u_era5_oper_pl_merge_1deg_daily_"+year+".nc"
fili_v =  folder_in+"/v/v_era5_oper_pl_merge_1deg_daily_"+year+".nc"
fili_t =  folder_in+"/t/t_era5_oper_pl_merge_1deg_daily_"+year+".nc"
fili_q =  folder_in+"/q/q_era5_oper_pl_merge_1deg_daily_"+year+".nc"


latN = 40
latS = -latN

fq = addfile(fili_q, "r")
fw = addfile(fili_w, "r")
fu = addfile(fili_u, "r")
fv = addfile(fili_v, "r")
ft = addfile(fili_t, "r")


; plev = 400.

ptop = 100.
pbot = 1000.

g    = 9.81          ; m/s2

coord_model = fw->w
printVarSummary(coord_model)

delete(coord_model)

q = short2flt(fq->q(:,{pbot:ptop},{latS:latN},:))
w = short2flt(fw->w(:,{pbot:ptop},{latS:latN},:))
u = short2flt(fu->u(:,{pbot:ptop},{latS:latN},:))
v = short2flt(fv->v(:,{pbot:ptop},{latS:latN},:))
T = short2flt(ft->t(:,{pbot:ptop},{latS:latN},:))

lon = u&lon 
lat = u&lat
time = u&time
print("Reading variables...")

printVarSummary(w)
printVarSummary(u)
printVarSummary(v)
printVarSummary(q)

ntime = dimsizes(time)

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

if any(ismissing(w)) then
    print("w contains missing values")
    poisson_grid_fill(w, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("w filled")
end if


if any(ismissing(q)) then
    print("q contains missing values")
    poisson_grid_fill(q, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("q filled")
end if


if any(ismissing(T)) then
    print("T contains missing values")
    poisson_grid_fill(T, is_cyclic, guess_type, nscan, epsx, relc, opt)
    print("T filled")
end if


rad = get_d2r("float")
re = 6.371e+06
; g=9.81
; Some of this are taken and modified from https://github.com/sandrolubis/Moisture_Budget/blob/main/cal_vint_moist_budget.ncl

print("Calculating temperature tendency...")
T1 = T
T1(1:ntime-1,:,:,:) = (/T(0:ntime-2,:,:,:)/)

dT_dt = (/(T-T1)/86400.0/)
copy_VarCoords(T, dT_dt)

delete(T1)
dT_dt@long_name = "Local time derivative of temperature"
dT_dt@units = "K/s"
printVarSummary(dT_dt)
printMinMax(dT_dt, False)

; print("Calculating temporal tendency of moisture...")

;   q1 = q
;   q1(1:ntime-1,:,:,:) = (/q(0:ntime-2,:,:,:)/)

;   dq_dt = (/(q-q1)/86400.0/)
;   copy_VarCoords(q, dq_dt)

; dq_dt@long_name = "Temporal tendency of moisture"
; dq_dt@units = "kg/kg s**-1"


;****************************************
;---Compute static stability [K/Pa] <==>  or "K-m-s2/kg"
;****************************************

;---Convert hPa -> Pa for function

p_raw = T&level
p       = p_raw*100.                                      ; Pa  [100000,...,10000]
delete(p_raw)
p@units = "Pa"
p!0     = "p"
p&p     =  p                   ; not necessary
printVarSummary(p)
print("---")

ss  = static_stability (p , T, 1, 0) 
printVarSummary(ss)
printMinMax(ss,0)
print("-----")

print("Calculating temperature advection...")
dT_dx = center_finite_diff_n(T, lon * rad, True, 0, 3) / conform(T, re * cos(lat * rad), 2)
dT_dy = center_finite_diff_n(T, lat * rad, False, 0, 2) / re

adv_T_zonal 		= -u * dT_dx
adv_T_meridional 	= -v * dT_dy
copy_VarCoords(T, adv_T_zonal)
adv_T_zonal@long_name = "Zonal advection of temperature"
copy_VarCoords(T, adv_T_meridional)
adv_T_zonal@units = "K/s"
printVarSummary(adv_T_zonal)
printMinMax(adv_T_zonal, False)

adv_T_meridional@long_name = "Meridional advection of temperature"
copy_VarCoords(T, adv_T_meridional)
adv_T_meridional@units = "K/s"
printVarSummary(adv_T_meridional)
printMinMax(adv_T_meridional, False)

adv_T = adv_T_zonal + adv_T_meridional
copy_VarCoords(adv_T_zonal, adv_T)
adv_T@long_name = "Advection of Horizontal temperature"
adv_T@units = "K/s"
printVarSummary(adv_T)
printMinMax(adv_T, False)

;****************************************
;---Net Heat
;****************************************
print("Calculating q1...")
q1 = dT_dt - (w*ss - adv_T)      ; term_1 - term_2
q1@long_name = "q1: Net Heat Flux"
q1@units     = "K/s"
copy_VarCoords(T,q1)
printVarSummary(q1) 
printMinMax(q1,0)
print("-----")

;****************************************
;---Apparent Heat Source
;****************************************

Cp           = 1004.64 
Cp@long_name = "specific heat of dry air at constant pressure"
Cp@units    = "J/(K-kg)"  ; [kg-m2/s2]/(K-kg) => [kg-m2]/[K-kg-s2] => m2/[K-s2]  

Q1_cp  = Cp*q1  ; [J/(K-kg)][K/s] -> [J/(kg-s)] -> [(kg-m2/s2)/(kg-s)) -> [ m2/s3 ] 
            ; [J/(K-kg)][K/s] -> [(J/s)(1/kg)] -> W/kg   ???
Q1_cp@long_name = "Q1: Total Diabatic Heating as the Apparent Heat Source"
Q1_cp@units     = "m2/s3"     
copy_VarCoords(T,Q1_cp)
printVarSummary(Q1_cp) 
printMinMax(Q1_cp,0)
print("-----") 



;****************************************
;---Make q1 to per day
;****************************************

q1 = q1*86400
q1@units     = "K/day"
q1@long_name = "q1: Net Heat Flux per day"
copy_VarCoords(T,q1)
; copy_VarCoords(T,q2)
   

;****************************************
;---Make Q1, Q2 to per day
;****************************************

Q1_cp = Q1_cp*86400.
Q1_cp@units     = "m2/s2 per day"
copy_VarCoords(T,Q1_cp)

coords = addfile("/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/daily_clim/era5/w/w_era5_doyclim_1979-2020.nc", "r")
coords_lat = coords->lat
coords_lon = coords->lon

fo_q1 = area_conserve_remap_Wrap (q1&lon,q1&lat,q1, coords_lon,coords_lat({latS:latN}), False)
; fo_dq_dt = area_conserve_remap_Wrap (dq_dt&lon,dq_dt&lat,dq_dt, coords_lon,coords_lat({latS:latN}), False)
fo_Q1_cp = area_conserve_remap_Wrap (Q1_cp&lon,Q1_cp&lat,Q1_cp, coords_lon,coords_lat({latS:latN}), False)


; fo_Q2_lv = area_conserve_remap_Wrap (Q2_lv&lon,Q2_lv&lat,Q2_lv, coords_lon,coords_lat({latS:latN}), False)


print("Writing output files...")
setfileoption("nc","Format","NETCDF4")
folder1 = "/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/q1/"
system("mkdir -p "+folder1)

filo1 = folder1+"q1_era5_oper_pl_merge_2.5deg_daily_eq_"+year+".nc"
system("rm -f "+filo1)

fout = addfile(filo1, "c")

; sanity check shapes match before writing anything
printVarSummary(fo_q1)
; printVarSummary(q2)
printVarSummary(fo_q1)

dimNames = (/"time","level","lat","lon"/)
dimSizes = (/ -1,    12,          29,   144   /)
dimUnlim = (/True,  False,      False, False/)

filedimdef(fout, dimNames, dimSizes, dimUnlim)
filedimdef(fout, dimNames, dimSizes, dimUnlim)

; coordinate variables first
filevardef(fout, "time", "double", "time")
filevardef(fout, "level", "float", "level")
filevardef(fout, "lat", typeof(fo_q1&lat), "lat")
filevardef(fout, "lon", typeof(fo_q1&lon), "lon")

; data variables
; filevardef(fout, "dq_dt", "float", (/"time","level","lat","lon"/))
filevardef(fout, "q1", "double", (/"time","level","lat","lon"/))
filevardef(fout, "Q1_cp", "double", (/"time","level","lat","lon"/))
; filevardef(fout, "q2", "float", (/"time","level","lat","lon"/))
; filevardef(fout, "qr_eddy", "float", (/"time","level","lat","lon"/))

; attributes
filevarattdef(fout, "time", fo_q1&time)
filevarattdef(fout, "level", fo_q1&level)
filevarattdef(fout, "lat", fo_q1&lat)
filevarattdef(fout, "lon", fo_q1&lon)
filevarattdef(fout, "q1", fo_q1)
filevarattdef(fout, "Q1_cp", fo_Q1_cp)
; filevarattdef(fout, "q2", q2)
; filevarattdef(fout, "qr_eddy", fo_q1)
; filevarattdef(fout, "dq_dt", dq_dt)

; write coords
fout->time       = (/fo_q1&time/)
fout->level      = (/fo_q1&level/)
fout->lat        = (/fo_q1&lat/)
fout->lon        = (/fo_q1&lon/)

; write data
fout->Q1_cp      = (/fo_Q1_cp/)
fout->q1      = (/fo_q1/)
; fout->qr_eddy = (/fo_q1/)
; fout->dq_dt = (/fo_dq_dt/)

print("checking...")
; filevardef(fout, "q2",      typeof(q2),    getvardims(q2))
; filevardef(fout, "qr_eddy", typeof(fo_q1), getvardims(fo_q1))
; fout->Q1_cp = fo_Q1_cp(:,{plev},:,:)
; fout->Q2_lv = fo_Q2_lv(:,{plev},:,:)

fout@title = "Yanai (1973) Q1 terms and dq_dt"
fout@history = "Created by NCL script on " + systemfunc("date")
fout@source_file = fili_u + " and " + fili_v + " and " + fili_w
fout@creation_date = systemfunc("date")
fout@institution = "The University of Melbourne"
fout@data = "ERA5"
; fout@notes = "Including the DQ/DT"