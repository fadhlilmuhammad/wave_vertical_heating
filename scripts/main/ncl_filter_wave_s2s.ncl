load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/s2s_utils.ncl"
load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/kf_filter.ncl"


; ens = ens 
; var = var 
; wave = wave

ens = "e01"
var = "q1"
wave = "mjo"

nint = 277
npad = 365 
nfil = npad*2-nint

latN = 30.
latS = -30

; day_min2 = day_min2
; day_min1 = day_min1 
; day_min0 = day_min0
; day_min2 = 19811230
; day_min1 = 19811231
; day_min0 = 19820101


day_min2 = 19820130
day_min1 = 19810131
day_min0 = 19820201


init = day_min0

year_0 = toint(str_get_cols(tostring(day_min0), 0, 3))
; year_3 = toint(year_0)-3
year_2 = toint(year_0)-2
year_1 = toint(year_0)-1 


; years = (/year_3, year_2, year_1, year_0/)
years = (/year_2, year_1, year_0/)
print(years)

folder = "/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/q1_vert_structure/"+ens+"/"
folder_obs = "/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/"+var+"/"

fname_min2 = "da_q1_yanai_"+day_min2+"_"+ens+"_remap.nc"
fname_min1 = "da_q1_yanai_"+day_min1+"_"+ens+"_remap.nc"
fname_min0 = "da_q1_yanai_"+day_min0+"_"+ens+"_remap.nc"

fname_obs = new(dimsizes(years), "string")

do i = 0, dimsizes(years)-1
    fname_obs(i) = folder_obs + var \
                 + "_era5_oper_pl_merge_2.5deg_daily_eq_" \
                 + years(i) + ".nc"
end do

print(fname_obs)
; fname_obs  = var+"_era5_oper_pl_merge_2.5deg_daily_eq_"+year+".nc"
fname_obsclim = "/home/563/fm6730/localrepo/wave_vertical_heating/data/raw/observation/"+var+"/daily_clim/"+var+"_era5_doyclim_1979-2020.nc"

f_min2 = addfile(folder + fname_min2, "r")
f_min1 = addfile(folder + fname_min1, "r")
f_min0 = addfile(folder + fname_min0, "r")
f_obs   = addfiles(fname_obs, "r")
f_clim  = addfile(fname_obsclim, "r")

x_min2 = f_min2->$var$(:,:,{latN:latS},:)
x_min1 = f_min1->$var$(:,:,{latN:latS},:)
x_min0 = f_min0->$var$(:,:,{latN:latS},:)
xobs   = f_obs[:]->$var$(:,:,{latN:latS},:)
clmDay  = f_clim->$var$(:,:,{latN:latS},:)
smth_clmDay = dble2flt((smthClmDayTLLL(clmDay, 3)))

printVarSummary(x_min2)
x_tlag = time_lag_ens_access_s2_TLLL(x_min2, x_min1, x_min0,nint)

printVarSummary(x_tlag)

ndims = dimsizes(x_tlag)
nens = ndims(0)
ntim = ndims(1)
nlev = ndims(2)
nlat = ndims(3)
nlon = ndims(4)


;***********************************************************
; Read time and create required yyyyddd                    
;***********************************************************
time    = x_tlag&time                          ; time:units = "hours since 1-1-1 00:00:0.0"                               
TIME    = cd_calendar(time, 0)          ; type float 
date    = cd_calendar(time, -2)
year    = toint( TIME(:,0) )
month   = toint( TIME(:,1) )
day     = toint( TIME(:,2) ) 
ddd     = day_of_year(year, month, day) 
yyyyddd = year*1000 + ddd                                  ; needed for input

time_obs = xobs&time 
TIME_obs    = cd_calendar(time_obs, 0)          ; type float 
date_obs  = cd_calendar(time_obs,-2)
year_obs    = toint( TIME_obs(:,0) )
month_obs   = toint( TIME_obs(:,1) )
day_obs    = toint( TIME_obs(:,2) ) 
ddd_obs    = day_of_year(year_obs, month_obs, day_obs) 
yyyyddd_obs = year_obs*1000 + ddd_obs      

init_mod = date(0)
; print(init_mod)
ind_init = toint(ind(init_mod.eq.date_obs))

ndimtlag = dimsizes(x_tlag)
x_tlag_anom = new((/ndimtlag(0), ndimtlag(1), ndimtlag(2), ndimtlag(3),ndimtlag(4)/),"float", x_tlag@_FillValue)
; x_tlag_anom = x_tlag
; printVarSummary(x_tlag_anom)
copy_VarMeta(x_tlag, x_tlag_anom)
;calculate anomaly and detrend

; do i = 0 , dimsizes(ddd)-1 
;     do e = 0, nens-1
;     ; print(doy)
;     d = ddd(i)
;     ; print(d)
;     x_tlag_anom(e,i,:,:,:) = x_tlag(e, i, :,:,:) - smth_clmDay({d},:,:,:)
;     end do 
; end do

printVarSummary(smth_clmDay)


do j = 0, nlev-1 
    do i = 0, nens-1 
        anom_temp = calcDayAnomTLL(x_tlag(i, :, j, :, :), yyyyddd, smth_clmDay(:,j,:,:))
        printVarSummary(anom_temp)
        x_tlag_anom(i,:,j,:,:) = anom_temp

        delete(anom_temp)
    end do
end do


x_anom_detrend = remove_linear_trend_TLLL(x_tlag_anom, nint)
printVarSummary(x_anom_detrend)

; PAD WITH OBSERVATION
print(ind_init)
printVarSummary(x_tlag)

x_obs_anom = xobs
do j = 0, nlev-1
    x_obs_anom(:,j,:,:) =  calcDayAnomTLL(xobs(:, j, :, :), yyyyddd_obs, smth_clmDay(:,j,:,:))
end do

print(date_obs(ind_init))
ind_init_obs = toint(ind(init.eq.date_obs))
; anomaly padded with obs anom
x_tlag_obs_pad = observation_padding_ensemble_TLLL(x_anom_detrend, x_obs_anom,ind_init, npad*2 )

print("half taper ensemble")
xtaper = half_taper_ensemble_TLLL(x_tlag_obs_pad,nfil, 0)

print("xtaper:")
print("  first missing = " + num(ismissing(xtaper)))

print("zero padding ensemble")
xpad_zero = zero_padding_ensemble_TLLL(xtaper, npad*2-nint-1)
printVarSummary((xpad_zero))


print("xpad_zero:")
print("  first missing = " + num(ismissing(xpad_zero)))

printMinMax(xpad_zero, 0)
print("assigning time")
time_padded = xpad_zero&time

; ; GET DATES
; date_pad = cd_calendar(x_tlag_obs_pad&time, -2)
; TIME_pad = cd_calendar(x_tlag_obs_pad&time, 0)
; year_pad    = toint( TIME_pad(:,0) )
; month_pad   = toint( TIME_pad(:,1) )
; day_pad     = toint( TIME_pad(:,2) ) 
; ddd_pad     = day_of_year(year_pad, month_pad, day_pad) 
; yyyyddd_pad = year_pad*1000 + ddd_pad     

;***********************************************************
; Read time and create required yyyyddd                    
;***********************************************************


;***********************************************************
; HIPASS FILTERING                 
;***********************************************************
if wave.eq."mjo".or.wave.eq."er" then 
    print("Will do HF filtering...")
    HF_filter = True 
    nwt = 101
    freq = 100.
else if wave.eq."kelvin".or.wave.eq."td".or.wave.eq."mrg" then 
    HF_filter = True 
    nwt = 101
    freq = 20.
else if wave.eq."LF" then 
    HF_filter = False 
end if
end if
end if

; x_hipass_raw = xpad_zero(time|:, lev|:, lat|:, lon|:, time_lag|:)

if HF_filter then
    print("HF filtering.... freq-"+freq)
    xpass = xpad_zero(:,0:npad*2+nint,:,:,:)
    printMinMax(xpass, 0)

    ; print(xpass(:,0,120,0))
    x_hipass_first = hipass_filtering_daily_ensemble_TLLL(xpass, nwt, freq, nens)
    printVarSummary(x_hipass_first)

    ; print(x_hipass_first(:,0,120,0))
    ; CHange the value to hipassed data
    xpad_zero(:,0:npad*2+nint,:,:,:) = x_hipass_first
    ; printMinMax(, 0)

    ; 
    ; print(x_hipass)

    ; Get rid of the first 50 (missing value not accepted for fft filtering)
    x_hipass = xpad_zero(time|50:,lev|:,lat|:,lon|:,time_lag|:)
    ; Force zero at the end of the tail
    x_hipass(npad*2+nint-50:,:,:,:,:) = 0.
    ; x_hipass(958-50:,:,:,:) = 0.
    if (any(ismissing(x_hipass))) then
        print("x_hipass contains missing values")
    else
        print("x_hipass has NO missing values")
    end if
    ; x_hipass = xpad_zero

    print("After filtering HF...")
    printVarSummary(x_hipass)
    ; print(x_hipass(:,{850},{10},{60},0))

else 
    print("Without HF filtering....")
    x_hipass = xpad_zero
end if

printVarSummary(x_hipass)
x_hipass@units= x_min0@units

print("======================================")
print("Missing value check")

print("x_tlag:")
print("  missing = " + num(ismissing(x_tlag)))

print("x_tlag_anom:")
print("  missing = " + num(ismissing(x_tlag_anom)))

print("x_anom_detrend:")
print("  missing = " + num(ismissing(x_anom_detrend)))

print("x_tlag_obs_pad:")
print("  missing = " + num(ismissing(x_tlag_obs_pad)))

print("xpad_zero:")
print("  missing = " + num(ismissing(xpad_zero)))

print("x_hipass_first:")
print("  missing = " + num(ismissing(x_hipass_first)))

print("x_hipass:")
print("  missing = " + num(ismissing(x_hipass)))
printMinMax(x_hipass, False)
; print("  minmax = " + printMinMax(x_hipass, 0))

print("======================================")

; do it = 0, ntim-1

;     tmp = xpad_zero(it,:,:,:,:)
;     tmp1d = ndtooned(tmp)

;     nmiss = num(ismissing(tmp1d))

;     if (nmiss.gt.0) then

;         miss = ind(ismissing(tmp1d))
;         coords = ind_resolve(miss, dimsizes(tmp))

;         do i = 0, nmiss-1
;             il   = coords(i,0)
;             iy   = coords(i,1)
;             ix   = coords(i,2)
;             ilag = coords(i,3)

;             print("Missing: time=" + it + \
;                   " lev=" + il + \
;                   " lat=" + xpad_zero&lat(iy) + \
;                   " lon=" + xpad_zero&lon(ix) + \
;                   " lag=" + ilag)
;         end do

;     end if

; end do

; system("rm "+"/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/check_hipass.nc")

; fout = addfile("/home/563/fm6730/localrepo/wave_vertical_heating/data/temp/wave_vertical_heating/check_hipass.nc", "c")
; fout->x_hipass = x_hipass
; do it = 0, ntim-1

;     if (any(ismissing(xpad_zero(it,:,:,:,:)))) then
;         print("Missing values at time index = " + it)
;     end if

; end do

;***********************************************************
; WAVE FILTERING                 
;***********************************************************

time_lag_hipass = x_hipass&time_lag
time_hipass = x_hipass&time
lon_hipass  = x_hipass&lon
lat_hipass  = x_hipass&lat
lev_hipass  = x_hipass&lev 
; printVarSummary((x_hipass))
printVarSummary((x_hipass))

filtered = new((/dimsizes(time_hipass), dimsizes(lev_hipass), dimsizes(lat_hipass), dimsizes(lon_hipass), dimsizes(time_lag_hipass)/), float, -999.)

ndims = dimsizes(filtered)

filtered!0 = "time"
filtered!1 = "lev"
filtered!2 = "lat"
filtered!3 = "lon"
filtered!4 = "time_lag"

filtered&time_lag = time_lag_hipass
filtered&time = time_hipass
filtered&lev  = lev_hipass
filtered&lat = lat_hipass
filtered&lon = lon_hipass

filtered@units = x_min0@units
obsPerDay = 1

mis = -999.
x_hipass@_FillValue = mis
filtered@_FillValue = mis


; printVarSum

print("wave filtering")

; wave = waves(i)
wave_type = wave

if (wave.eq."mjo") then
    wavenumber = (/1, 5/)
    period = (/20, 100/)
    depth = (/mis, mis/)
    filtered@long_name = "Madden-Julian Oscillations in Vorticity "
else if (wave.eq."kelvin") then
    wavenumber = (/1, 14/)
    period = (/2.5, 17/)
    depth = (/8, 90/)
    filtered@long_name = "Kelvin Waves in Vorticity "
else if (wave.eq."er") then
    wavenumber = (/-10, -1/)
    period = (/9, 72/)
    depth = (/8, 90/)
    filtered@long_name = "Equatorial Rossby Waves in Vorticity "
else if (wave.eq."mrg") then
    wavenumber = (/-10, -1/)
    period = (/3, 10/)
    depth = (/8, 90/)
    filtered@long_name = "Mixed Rossby-Gravity Waves in Vorticity "
else if (wave.eq."eig") then
    wavenumber = (/1, 14/)
    period = (/1, 5/)
    depth = (/12, 50/)
    wave_type = "ig0"
    filtered@long_name = "Eastward Inertio Gravity Waves in Vorticity "
else if (wave.eq."td") then
    wavenumber = (/-20, -6/)
    ; period = (/2.5, 10/)
    period = (/2.5, 5/)
    depth = (/mis, 90/)
    wave_type = "td-type"
    filtered@long_name = "Tropical Depression-Type Waves in Vorticity "
else if (wave.eq."LF") then
    wavenumber = (/-99, 99/)
    ; period = (/2.5, 10/)
    period = (/120, 9999/)
    depth = (/mis, mis/)
    wave_type = "LF"
    filtered@long_name = "Low Frequency Variabilities in Vorticity "
end if
end if
end if
end if
end if
end if
end if

do l = 0, nlev-1
    print(lev_hipass(l))
    do tl = 0, nens-1
    print(time_lag_hipass(tl))
        do y = 0, (dimsizes(lat_hipass) - 1)
            filtered(time|:, lev|l, lat|y, lon|:, time_lag|tl) = (/kf_filter(x_hipass(time|:, lev|l, lat|y, lon|:, time_lag|tl), obsPerDay, period(0), period(1), wavenumber(0), wavenumber(1), depth(0), depth(1), wave_type)/)
            ; cek = (/kf_filter(x_hipass(time|:, lev|l, lat|y, lon|:, time_lag|tl), obsPerDay, period(0), period(1), wavenumber(0), wavenumber(1), depth(0), depth(1), wave_type)/)
            ; print(cek)
            print((/lat_hipass(y)/))
            ; printVarSummary(filtered)
            ; printMinMax((filtered), 0)
        end do
    end do
end do

; filtered = x_hipass

; printVarSummary(filtered)
; printMinMax(filtered, False)


; print(filtered(0,:,0,120))

rad    = 4.0*atan(1.0)/180.0


filtered_reorder = filtered(time|:,lev|:, {lat|-15:15}, lon|:, time_lag|:)
printVarSummary(filtered_reorder)


clat   = cos(filtered_reorder&lat*rad)
printVarSummary(clat)


opt=2 ;opt>1 means to return the mean only if the number of non-missing values is greater-than or equal to opt.
filtered_vert = dim_avg_wgt_n_Wrap(filtered_reorder, clat, opt, 2)

outname = wave+"_tlag_da_"+var+"_"+init+"_"+ens+"_remap_padded"
folder_out = "/home/563/fm6730/localrepo/wave_vertical_heating/data/derived/filtered/q1_vert_structure/"+wave+"/"

system("mkdir -p "+folder_out)

print("Opening the output file...")
system("rm "+folder_out+wave+"_"+outname+".nc")
g = addfile(folder_out+wave+"_"+outname+".nc", "c")

print("Writing the output file...")

date_filtered = cd_calendar(filtered_vert&time, -2)
idinit = ind(date_filtered.eq.init)
print(idinit)
print(date_filtered(idinit))
if wave.eq."LF" then 
    g->$var$ = filtered_vert(time|:,lev|:,lon|:,time_lag|:)
else
    g->$var$ = filtered_vert(time|idinit:idinit+nint,lev|:,lon|:,time_lag|:)
end if

; print(filtered(time|idinit:idinit+nint,lat|0,lon|0,time_lag|0))
print("Closing the output file...")

; delete(filtered)

check_memory("On-the-fly Data Preparation & Filtering ")
; fout = addfile("check_hipass.nc", "c")
; fout->x_hipass = x_hipass


; plot_2d_latlon(x_hipass, 461, 5, 0, \
;                "hipass Q1 at 500 hPa", \
;                "hipass_500hPa")

; plot_2d_latlon(x_tlag_anom, 0, 4, 0, \
;                "tlag_anom Q1 at 500 hPa", \
;                "check_tlag_anom_500hPa")


; plot_2d_latlon(x_tlag, 0, 4, 0, \
;                "tlag Q1 at 500 hPa", \
;                "check_tlag_500hPa")

