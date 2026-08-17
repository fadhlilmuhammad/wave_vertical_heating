; UTILITY codes
;

load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/wkSpaceTime_mod.ncl"
load "/home/563/fm6730/localrepo/wave_vertical_heating/scripts/utilities/utils.ncl"

procedure check_memory(label:string)
begin
    print("======================================")
    print("Memory check: " + label)
    print(systemfunc("ps -C ncl -o pid,%cpu,%mem,rss,vsz,cmd"))
    ; print(systemfunc("ps -C ncl -o pid,rss,vsz,%mem,cmd"))

    print("======================================")
end

;===========================================================
; Plot one 2-D lat-lon slice from x_hipass
; x_hipass dimensions: (time, lev, lat, lon, time_lag)
;===========================================================
procedure plot_2d_latlon(x, itime, ilev, ilag, title, outfile)

begin

    ; Select one 2-D lat-lon field
    data = x(itime, ilev, :, :, ilag)

    wks = gsn_open_wks("png", outfile)

    res = True
    res@gsnDraw  = False
    res@gsnFrame = False

    res@cnFillOn             = True
    res@cnLinesOn            = False
    res@cnLineLabelsOn      = False
    res@cnFillPalette       = "BlWhRe"

    res@lbLabelBarOn         = True

    res@tiMainString         = title
    res@tiMainFontHeightF    = 0.015

    res@gsnAddCyclic         = False

    plot = gsn_csm_contour_map(wks, data, res)

    draw(plot)
    frame(wks)

    delete(data)
    delete(wks)

end


;BANDPASS FILTERING
undef("bandpass_filtering_daily")
function bandpass_filtering_daily(\
  data[*][*][*]  : numeric,\
  nwgt           : integer, \
  freqlow_in_days   : float, \ 
  freqhi_in_days   : float \ 
)
local None 
begin 
  fca    = 1./freqlow_in_days
  fcb    = 1./freqhi_in_days
  ihp    = 2          
  nsigma = 1.
  wgt    = filwgts_lanczos (nwgt, ihp, fca, fcb, nsigma)  

  filtered = wgt_runave_n_Wrap(data, wgt, 0, 0)
  return(filtered)  
end 

;HIPASS FILTERING
undef("hipass_filtering_daily")
function hipass_filtering_daily(\
  data[*][*][*]  : numeric,\
  nwgt           : integer, \
  freq_in_days   : float \ 
)

local None 
begin 
  fca    = 1./freq_in_days
  ; ; ; print(fcb)
  ihp    = 1          
  nsigma = 1.
  wgt    = filwgts_lanczos (nwgt, ihp, fca, -999., nsigma)  

  filtered = wgt_runave_n_Wrap(data, wgt, 0, 0)
  return(filtered)  
end 

undef("hipass_filtering_daily_ensemble")
function hipass_filtering_daily_ensemble(\
  data[*][*][*][*]  : numeric, \
  nwgt              : integer, \
  freq_in_days      : float,   \
  nens              : integer  \
)
local filtered_ens
begin 
  filtered_ens = data
  copy_VarMeta(data, filtered_ens)
  ; printVarSummary(filtered_ens)

  do tlag = 0, nens-1
  ; ; printVarSummary((data(time_lag|tlag,time|:,lat|:,lon|:), nwgt, freq_in_days))
    filtered_ens(time_lag|tlag,time|:,lat|:,lon|:) = hipass_filtering_daily(data(time_lag|tlag,time|:,lat|:,lon|:), nwgt, freq_in_days)     
  end do

  return(filtered_ens)
end

undef("hipass_filtering_daily_ensemble_TLLL")
function hipass_filtering_daily_ensemble_TLLL(\
  data[*][*][*][*][*]  : numeric, \
  nwgt              : integer, \
  freq_in_days      : float,   \
  nens              : integer  \
)
local filtered_ens
begin 
  filtered_ens = data
  copy_VarMeta(data, filtered_ens)
  ; printVarSummary(filtered_ens)

  nlev = dimsizes(data&lev)

  do tlag = 0, nens-1
    do l = 0, nlev-1
  ; ; printVarSummary((data(time_lag|tlag,time|:,lat|:,lon|:), nwgt, freq_in_days))
      filtered_ens(time_lag|tlag,time|:,lev|l,lat|:,lon|:) = hipass_filtering_daily(data(time_lag|tlag,time|:,lev|l, lat|:,lon|:), nwgt, freq_in_days)     
    end do
  end do

  return(filtered_ens)
end


;LOWPASS FILTERING
undef("lowpass_filtering_daily")
function lowpass_filtering_daily(\
  data[*][*][*]  : numeric,\
  nwgt           : integer, \
  freq_in_days   : float \ 
)

local None 
begin 
  fcb    = 1./freq_in_days
  ihp    = 0          
  nsigma = 1.
  wgt    = filwgts_lanczos (nwgt, ihp, -999., fcb, nsigma)  

  filtered = wgt_runave_n_Wrap(data, wgt, 0, 0)
  return(filtered)  
end 


; data must be 3D (time, lat, lon)
; nint must be integer
; tname is the name of the time dimension, in string

undef ("obs_integration")
function obs_integration (\
  data[*][*][*]  : numeric, \
  nint           : integer, \
  tname          : string   \
  )

local None
begin

  ndims = dimsizes(data)
  ntim  = ndims(0)
  nlat  = ndims(1)
  nlon  = ndims(2)

  time  = data&$tname$

  x     = new((/ntim-nint, nint, nlat, nlon/), float)
  x!0   = "n_init"


  do t = 0, ntim-nint-1
    x(t,:,:,:) = data(t:t+(nint-1),:,:)
  end do

  delete(x&$tname$)
  x!1 = "time"
  x&time = ispan(0,nint-1,1)
  ; printVarSummary(x)

  return(x)
end


;TIME LAGGED ENSEMBLE FOR ACCESS-S2
undef ("time_lag_ens_access_s2")
function time_lag_ens_access_s2 (\
  data_lag2[*][*][*]   : numeric, \
  data_lag1[*][*][*]   : numeric, \
  data_lag0[*][*][*]   : numeric, \
  nint                 : integer \
)

local None
begin
  x_lag2     = data_lag2(2:nint+2-1,:,:)   ;remove the last 2 days
  x_lag1     = data_lag1(1:nint+1-1,:,:)   ;remove the last 1 day
  x_lag0     = data_lag0(0:nint-1,:,:)   ;reference day

  ndims  = dimsizes(x_lag2)
  nlat   = ndims(1)
  nlon   = ndims(2)

  x      = new((/3, nint, nlat, nlon/), float)

  x(0,:,:,:) = x_lag2
  x(1,:,:,:) = x_lag1
  x(2,:,:,:) = x_lag0

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)
end


;TIME LAGGED ENSEMBLE FOR ACCESS-S2
undef ("time_lag_ens_access_s2_TLLL")
function time_lag_ens_access_s2_TLLL (\
  data_lag2[*][*][*][*]   : numeric, \
  data_lag1[*][*][*][*]   : numeric, \
  data_lag0[*][*][*][*]   : numeric, \
  nint                 : integer \
)

local None
begin
  x_lag2     = data_lag2(2:nint+2-1,:,:,:)   ;remove the last 2 days
  x_lag1     = data_lag1(1:nint+1-1,:,:,:)   ;remove the last 1 day
  x_lag0     = data_lag0(0:nint-1,:,:,:)   ;reference day

  ndims  = dimsizes(x_lag2)
  nlev   = ndims(1)
  nlat   = ndims(2)
  nlon   = ndims(3)

  x      = new((/3, nint, nlev, nlat, nlon/), float)

  x(0,:,:,:,:) = x_lag2
  x(1,:,:,:,:) = x_lag1
  x(2,:,:,:,:) = x_lag0

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)
end


;;;;GET THE 4TH HARMONICS

undef ("get_harmonic_cycle")
function get_harmonic_cycle (\
  data[*][*][*]       : numeric,\
  nharm               : integer \
)

local time, TIME, year, month, day, ddd, yyyyddd
begin
   time    = data&time                          ; time:units = "hours since 1-1-1 00:00:0.0"
   TIME    = cd_calendar(time, 0)          ; type float
   year    = toint( TIME(:,0) )
   month   = toint( TIME(:,1) )
   day     = toint( TIME(:,2) )
                                                                 ; check for calendar attribute
   if (isatt(TIME,"calendar")) then           ; default is gregorian
       year@calendar = TIME@calendar
   end if

   ddd     = day_of_year(year, month, day)
   if (isatt(year,"calendar")) then           ; default is gregorian
       ddd@calendar = year@calendar
   end if

   yyyyddd = year*1000 + ddd                                    ; needed for input
   if (isatt(ddd,"calendar")) then           ; default is gregorian
       yyyyddd@calendar = ddd@calendar
   end if

  hClmDay    = clmDayTLL(data, yyyyddd)
  x = smthClmDayTLL(hClmDay, nharm)

return (x)
end


undef ("remove_harmonic_cycle")
function remove_harmonic_cycle(\
  data[*][*][*]       : numeric,\
  data_clm[*][*][*]   : numeric,\
  nharm               : integer \
)

local time, TIME, year, month, day, ddd, yyyyddd
begin
  xclm  = get_harmonic_cycle(data_clm, nharm)

  ; ; printVarSummary(xclm)
  time = data&time
  ; time@calendar = "standard"

  utc_date = cd_calendar(time,0)
  year   = tointeger(utc_date(:,0))    ; Convert to integer for
  month  = tointeger(utc_date(:,1))    ; use s; printi
  day    = tointeger(utc_date(:,2))

  xAnom  = data
  doy    = day_of_year(year, month, day)
  ; ; ; print(doy)

  ; ; printVarSummary(xclm)
  do i = 0, dimsizes(doy)-1
    ; ; ; print(i)
    d = doy(i)
    ; ; printVarSummary(data(i,:,:))
    ; ; printVarSummary(xclm({d},:,:))
    xAnom(i,:,:) = data(i,:,:) - xclm({d},:,:); write loop content
  end do


  ; xAnom = calcDayAnomTLL(data, yyyyddd, xclm)
  xAnom@long_name = "Anomalies from Smooth Daily Climatology harmonic-"+nharm
  copy_VarMeta(data, xAnom)
  ; copy_VarCoords()

  return(xAnom)
end


undef ("remove_harmonic_cycle_ensemble")
function remove_harmonic_cycle_ensemble (\
  data[*][*][*][*]    : numeric,\
  data_clm[*][*][*]   : numeric,\
  nharm               : integer \
)

begin
  ndims = dimsizes(data)
  nens  = ndims(0)

  x     = new((/nens, ndims(1), ndims(2), ndims(3)/), typeof(data))

  do z = 0, nens-1
    x(z,:,:,:) = remove_harmonic_cycle(data(z,:,:,:), data_clm, nharm); write loop content
  end do

  copy_VarMeta(data, x)
  return(x)
end

undef("lead_time_climatology")
function lead_time_climatology(\
  data[*][*][*][*]  : numeric \
)

begin 

  ndims  = dimsizes(data)
  ninit  = ndims(0)
  ntim   = ndims(1)
  nlat   = ndims(2)
  nlon   = ndims(3)

  clim = dim_avg_n_Wrap(data, 0)
  return(clim)

end 


undef("lead_time_climatology_tlag")
function lead_time_climatology_tlag(\
  data[*][*][*][*][*]  : numeric \
)

begin 

  ; printVarSummary(data)
  ndims  = dimsizes(data)
  ninit  = ndims(0)
  ntlag  = ndims(1)
  ntim   = ndims(2)
  nlat   = ndims(3)
  nlon   = ndims(4)

  clim_1 = dim_avg_n_Wrap(data, 0)
  ; printVarSummary(clim_1)
  clim = dim_avg_n_Wrap(clim_1, 0)
  ; printVarSummary(clim)
  return(clim)

end 

undef ("remove_model_mean_bias")
function remove_model_mean_bias(\
  model[*][*][*]      : numeric,\
  nint               : integer, \ ;number of forecast days
  model_mean[*][*][*]    : numeric \
)
begin

  ; ; print("checking dimensions...")
  dimsmodel = dimsizes(model)
  dimsmean  = dimsizes(model_mean)

  if dimsmodel(0).ne.dimmean(0) then
    ; ; print("Error: dimension 0 of model and model_mean must be the same")
  else if dimsmodel(1).ne.dimsmean(1) then 
    ; ; print("Error: dimension 1 of model and model_mean must be the same")
  else if dimsmodel(2).ne.dimsmean(2) then 
    ; ; print("Error: dimension 2 of model and model_mean must be the same")
  end if 
  end if 
  end if

  ; model_mean =  dim_avg_n_Wrap(model(0:ninit-1,:,:), 0)
  ; xanom  = model(0:ninit-1,:,:) - conform_dims(dimsizes(model(0:ninit-1,:,:)), model_mean, (/1,2/))
  xanom  = model(0:nint-1,:,:) - model_mean
  copy_VarMeta(model(0:nint-1,:,:), xanom)
  return(xanom)
end

undef("remove_time_mean")
function  remove_time_mean(\
  model[*][*][*][*]  : numeric, \
  ninit             : integer \
)
begin 

  ndims = dimsizes(model)
  xanom = new((/ndims(0), ninit, ndims(2), ndims(3)/), typeof(model))  
  do e = 0, ndims(0)-1
    ; ; printVarSummary(model(e,0:ninit-1,:,:))
    model_mean =  dim_avg_n_Wrap(model(e,0:ninit-1,:,:), 0)
    xanom(e,:,:,:)  = model(e,0:ninit-1,:,:) - conform_dims(dimsizes(model(e,0:ninit-1,:,:)), model_mean, (/1,2/))
    copy_VarMeta(model(e,0:ninit-1,:,:), xanom(e,:,:,:))

  end do
    xanom!0 = "time_lag"
    xanom&time_lag = (/-2,-1,0/)
    xanom@time_lag = "lagged by 2, 1, and 0 days"

    xanom!1 = "time"
    time = model&time
    xanom&time = time(0:ninit-1)

    xanom!2 = "lat"
    xanom&lat = model&lat

    xanom!3 = "lon"
    xanom&lon = model&lon
    ; printVarSummary(xanom)
    return(xanom)
end 

undef("remove_linear_trend")
function remove_linear_trend(\
  model[*][*][*][*]  : numeric, \
  nint             : integer \
)
begin

  time = model&time
  ; printVarSummary(time)
  ndims = dimsizes(model)
  xanom = new((/ndims(0), nint, ndims(2), ndims(3)/), typeof(model)) 
  ; ; ; print(ndims(0)) 
  do e = 0, ndims(0)-1
    ; ; ; print(e)
    xanom(e,:,:,:) = dtrend_msg_n(time(0:nint-1),model(e,0:nint-1,:,:), True, False,0)
  end do

  xanom!0 = "time_lag"
  xanom&time_lag = (/-2,-1,0/)
  xanom@time_lag = "lagged by 2, 1, and 0 days"

  xanom!1 = "time"
  xanom&time = time(0:nint-1)

  xanom!2 = "lat"
  xanom&lat = model&lat

  xanom!3 = "lon"
  xanom&lon = model&lon
  ; printVarSummary(xanom)
  return(xanom)
end



undef("remove_linear_trend_TLLL")
function remove_linear_trend_TLLL(\
  model[*][*][*][*][*]  : numeric, \
  nint             : integer \
)
begin

  time = model&time
  ; printVarSummary(time)
  ndims = dimsizes(model)
  ; ; ; print(ndims)
  xanom = new((/ndims(0), nint, ndims(2), ndims(3), ndims(4)/), typeof(model)) 
  ; ; ; print(ndims(0)) 
  do e = 0, ndims(0)-1
    do l = 0, ndims(2)-1
    ; ; ; print(e)
      xanom(e,:,l,:,:) = dtrend_msg_n(time(0:nint-1),model(e,0:nint-1,l,:,:), True, False,0)
    end do
  end do

  xanom!0 = "time_lag"
  xanom&time_lag = (/-2,-1,0/)
  xanom@time_lag = "lagged by 2, 1, and 0 days"

  xanom!1 = "time"
  xanom&time = time(0:nint-1)

  xanom!2 = "lev"
  model!2 = "lev"
  xanom&lev = model&lev

  xanom!3 = "lat"
  xanom&lat = model&lat

  xanom!4 = "lon"
  xanom&lon = model&lon
  ; printVarSummary(xanom)
  return(xanom)
end



undef ("remove_model_mean_bias_ensemble")
function remove_model_mean_bias_ensemble(\
  model[*][*][*][*]   : numeric,\
  ninit               : integer, \
  model_mean[*][*][*] : numeric \
)
begin

  ; ; print("checking dimensions...")
  dimsmodel = dimsizes(model)
  dimsmean  = dimsizes(model_mean)


  if dimsmodel(0).ne.dimmean(0) then
    ; ; print("Error: dimension 0 of model and model_mean must be the same")
  else if dimsmodel(1).ne.dimsmean(1) then 
    ; ; print("Error: dimension 1 of model and model_mean must be the same")
  else if dimsmodel(2).ne.dimsmean(2) then 
    ; ; print("Error: dimension 2 of model and model_mean must be the same")
  else if dimsmodel(3).ne.dimsmean(3) then 
    ; ; print("Error: dimension 3 of model and model_mean must be the same")
  end if
  end if 
  end if 
  end if

  ndims = dimsizes(model)

  xanom = new((/ndims(0), ninit, ndims(2), ndims(3)/), typeof(model))    
  ; ; printVarSummary(model_mean)
  ; printVarSummary(xanom)
  do e = 0, ndims(0)-1
    ; ; printVarSummary(model(e,0:ninit-1,:,:))
    ; model_mean =  dim_avg_n_Wrap(model(e,0:ninit-1,:,:), 0)
    ; xanom(e,:,:,:)  = model(e,0:ninit-1,:,:) - conform_dims(dimsizes(model(e,0:ninit-1,:,:)), model_mean, (/1,2/))
    xanom(e,:,:,:) = model(e,0:ninit-1,:,:) - model_mean
    ; copy_VarMeta(model(e,0:ninit-1,:,:), xanom(e,:,:,:))

  end do
    xanom!0 = "time_lag"
    xanom&time_lag = (/-2,-1,0/)
    xanom@time_lag = "lagged by 2, 1, and 0 days"

    xanom!1 = "time"
    time = model&time
    xanom&time = time(0:ninit-1)

    xanom!2 = "lat"
    xanom&lat = model&lat

    xanom!3 = "lon"
    xanom&lon = model&lon
    ; printVarSummary(xanom)
    return(xanom)
end


undef ("observation_padding")
function observation_padding(\
  data[*][*][*]      : numeric,\
  pad[*][*][*]       : numeric,\
  ind_init           : integer,\
  npad               : integer \
)

local time, date, time_add, ndims, ntim,nlat,nlon,x
begin

  data!0 = "time"
  time_f  = data&time
  date  = cd_calendar(time_f, -2)

  ; ; ; print(time)
  ; ; printVarSummary(time)

  time_add = new((/dimsizes(date)+npad/), typeof(date))

  ; ; printVarSummary(time_add)

  time = floattointeger(tofloat(time_f))
  do d = 0, npad-1
    time_add(d) = time(0) - npad + d
  end do

  ; ; printVarSummary(time)

  time_add(npad:) = time
  ; printVarSummary((time_add))
  ; ; ; print(cd_calendar(time_add,-2))

  ndims = dimsizes(data)
  ntim  = ndims(0)
  nlat  = ndims(1)
  nlon  = ndims(2)

  x = new((/ntim+npad, nlat, nlon/), typeof(data))

  ; ; printVarSummary(pad((ind_init - npad):ind_init-1, :, :))

  ; pad = short2flt(pad_in)
  ; ; printVarSummary(pad_in)
  ; ; printVarSummary(data)
  ; ; printVarSummary(pad)
  ; ; print("padding pad...")
  pad_converted = (/tofloat(pad)/)
  x(0:npad-1,:,:) = (/pad_converted((ind_init - npad):ind_init-1, :, :)/)
  ; ; print("padding_data")
  ; printVarSummary(x(npad:,:,:))
  ; printVarSummary(data)

    ; Convert pad array to match x type if needed
  if (typeof(data) .ne. typeof(x)) then
    ; ; print("conversion needed, x and data are different type")
    data_converted = tofloat(data)  ; or todouble() depending on type needed
  else
    ; ; print("x and data are same type")
    data_converted = data
  end if

  x(npad:,:,:)    = (/data_converted/)

  x!0 = "time"
  x!1 = "lat"
  x!2 = "lon"

  time_add@units = time_f@units
  time_add@calendar = time_f@calendar
  time_add@standard_name = time_f@standard_name
  time_add@axis = time_f@axis

  x&time = time_add
  x&lat  = data&lat
  x&lon  = data&lon  

  ; ; print("finish padding")


  ; copy_VarAtts(data&time, x&time)

  ; time_float = tofloat(x&time)
  ; x&time = time_float

  return(x)
end

undef ("observation_padding_ensemble")
function observation_padding_ensemble(\
  data[*][*][*][*]   : numeric,\
  pad[*][*][*]       : numeric,\
  ind_init           : integer,\
  npad               : integer \
)
begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlat  = ndims(2)
  nlon  = ndims(3)
  x = new((/nens, ntim+npad, nlat, nlon/), typeof(data))

  ; ; print(typeof(data))
  pad_converted = pad
  ; delete(pad)
  ; pad = pad_converted

  do e = 0, nens-1
    x(e,:,:,:) = observation_padding(data(e,:,:,:), pad_converted, ind_init, npad)
  end do
  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)
end


undef ("observation_padding_ensemble_TLLL")
function observation_padding_ensemble_TLLL(\
  data[*][*][*][*][*]   : numeric,\
  pad[*][*][*][*]       : numeric,\
  ind_init           : integer,\
  npad               : integer \
)
begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlev  = ndims(2)
  nlat  = ndims(3)
  nlon  = ndims(4)
  x = new((/nens, ntim+npad,  nlev, nlat, nlon/), typeof(data))
  ; printVarSummary(x)

  data!2 = "lev"
  ; ; ; print(typeof(data))
  pad_converted = pad
  ; delete(pad)
  ; pad = pad_converted

  ; ; printVarSummary(data(0,:,0,:,:))
  do e = 0, nens-1
    do lev = 0, nlev-1
      ; ; printVarSummary(observation_padding(data(e,:,lev,:,:), pad_converted(:,lev,:,:), ind_init, npad))
      x(e,:,lev,:,:) = observation_padding(data(e,:,lev,:,:), pad_converted(:,lev,:,:), ind_init, npad)
    end do 
  end do


  x!2 = "lev"
  x&lev = data&lev
  x!0 = "time_lag"

  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)
end


undef ("zero_padding")
function zero_padding(\
  data[*][*][*]     : numeric, \
  npad              : integer
)

local time, date, time_add, ndims, ntim,nlat,nlon
begin
  data!0 = "time"
  time  = data&time
  ; printVarSummary(time)
  date  = cd_calendar(time,-2)

  time_add = new((/dimsizes(date)+npad/), typeof(time))

  ; printVarSummary(time)
  ; printVarSummary(time_add)
  time_add(0:dimsizes(date)-1) = time

  ; ; ; print(time(dimsizes(date)-1)+1)

  ; ; ; print("datae awal padding")
  ; ; ; print(cd_calendar(time_add(dimsizes(date)-1)+1,-2))
  ; ; ; print("datae akhir data asli")
  ; ; ; print(cd_calendar(time(dimsizes(date)-1),-2))
  
  do d = 0, npad-1
    time_add(dimsizes(date)+d) = time(dimsizes(date)-1) + d+1
  end do

  ; ; print("data awal="+cd_calendar(time(0),-2))
  ; ; print("data akhir forecast="+cd_calendar(time(dimsizes(date)-1),-2))
  ; ; print("data awal padding forecast="+cd_calendar(time_add(dimsizes(date)),-2))

  ; ; ; print(time_add)
  ndims = dimsizes(data)
  ntim  = ndims(0)
  nlat  = ndims(1)
  nlon  = ndims(2)

  x = new((/ntim+npad, nlat, nlon/), typeof(data))
  x(0:ntim-1,:,:)  = data
  x(ntim:ntim+npad-1,:,:) = 0

  x!0 = "time"
  x!1 = "lat"
  x!2 = "lon"

  x&time = time_add
  x&lat  = data&lat
  x&lon  = data&lon

  return(x)

end

undef ("zero_padding_ensemble")
function zero_padding_ensemble(\
data[*][*][*][*]  : numeric, \
npad              : integer
)
begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlat  = ndims(2)
  nlon  = ndims(3)

  x    = new((/nens, ntim+npad, nlat, nlon/), float)
  do e = 0, nens-1
     tmp = zero_padding(data(e,:,:,:), npad)
    ;  ; printVarSummary((tmp))
     x(e,:,:,:) = tmp
  end do

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)

end


undef ("zero_padding_ensemble_TLLL")
function zero_padding_ensemble_TLLL(\
data[*][*][*][*][*]  : numeric, \
npad              : integer
)
begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlev  = ndims(2)
  nlat  = ndims(3)
  nlon  = ndims(4)

  x    = new((/nens, ntim+npad, nlev, nlat, nlon/), float)
  do l = 0, nlev-1
    do e = 0, nens-1
      tmp = zero_padding(data(e,:,l,:,:), npad)
      ;  ; printVarSummary((tmp))
      x(e,:,l,:,:) = tmp
    end do
  end do

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  x!2 = "lev" 
  x&lev = data&lev

  return(x)

end


undef ("data_taper_n")
function data_taper_n(\
  data[*][*][*]     : numeric,\
  n                 : integer,\
  days_to_taper     : integer \
)
begin
  ndim  = dimsizes(data)

  nsize = ndim(n)

  p     = tofloat(days_to_taper)/tofloat(nsize)
  ; ; ; print(p)

  x     = taper_n(data, p, 1, n)
  copy_VarMeta(data, x)

  return(x)

end

undef("half_taper")
function half_taper(\
  data[*][*][*]     : numeric,\
  days_to_taper     : integer,\
  force_zero        : integer \
)
begin
  ndim  = dimsizes(data)

  nsize = ndim(0)

  p     = tofloat(days_to_taper)/tofloat(nsize)
  ; ; ; print(p)


  wgt = new((/nsize, ndim(1), ndim(2)/), typeof(data))
  do t = 0, nsize-1
    wgt(t,:,:)   = 0.5*(1 - cos(get_pi("f")*(t+0.5)/nsize))
    if force_zero.eq.1 then
      wgt(t,:,:) = where(wgt(t,:,:).le.0.1, 0, wgt(t,:,:))
    end if
  end do

  ;https://pages.stern.nyu.edu/~churvich/TimeSeries/Handouts/Notes3.pdf
  ; ; ; print(wgt(:,20,100))
  x = data
  x(0:days_to_taper-1,:,:) = data(0:days_to_taper-1,:,:) * wgt(0:days_to_taper-1,:,:)
  ; x     = taper_n(data, p, 1, n)
  copy_VarMeta(data, x)

  return(x)

end

undef ("half_taper_ensemble")
function half_taper_ensemble(\
  data[*][*][*][*]   : numeric,\
  days_to_taper      : integer,\
  force_zero         : integer \
)

begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlat  = ndims(2)
  nlon  = ndims(3)

  x = new((/nens, ntim, nlat, nlon/), typeof(data))
  do e = 0, nens-1
    x(e, :,:,:) = half_taper(data(e,:,:,:), days_to_taper, force_zero)
  end do

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)
end

undef ("half_taper_ensemble_TLLL")
function half_taper_ensemble_TLLL(\
  data[*][*][*][*][*]   : numeric,\
  days_to_taper      : integer,\
  force_zero         : integer \
)

begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlev  = ndims(2)
  nlat  = ndims(3)
  nlon  = ndims(4)

  x = new((/nens, ntim, nlev, nlat, nlon/), typeof(data))
  do e = 0, nens-1
   do l = 0, nlev-1
      x(e, :,l,:,:) = half_taper(data(e,:,l,:,:), days_to_taper, force_zero)
    end do
  end do

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  x!2 = "lev"
  x&lev = data&lev

  return(x)
end


undef ("wkSpaceTime_S2S")
procedure wkSpaceTime_S2S(\
  data[*][*][*]     : numeric,\
  data_clm[*][*][*] : numeric,\
  var               : string, \
  diro              : string, \
  case              : string, \
  latN              : numeric,\
  cpd               : integer,\
  nDayWin           : integer,\
  nDaySkip          : integer,\
  opt \
)

begin
  data_clm!0 = "time"
  time_s2s   = data&time
  time       = data_clm&time

  date_s2s   = cd_calendar(time_s2s, 0)
  date       = cd_calendar(time, 0)
  ; ; ; print(date_s2s)
  mon_s2s    = date_s2s(0,1)

  ; ; ; print(mon_s2s)
  seasname   = name_season(mon_s2s)

  ; ; print("seasname="+seasname)
  if (isatt(opt, "season")) then
    iseason    = choose_season(seasname, date)
    ; ; print(iseason)
    data_back  = data_clm(iseason,:,:)
  
  else
    data_back  = data_clm
  end if

  
  opt_back   = False
  spec_back  = SpaceTime_Background(data_back, diro, case, var             \
              ,latN, cpd, nDayWin, nDaySkip, opt_back  )
  ;
  ; opt@Fig_3a      = log10(/1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 /)
  ; opt@Fig_3b      = log10(/1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 /) ; unequal
  ; opt@pltColorMap = "WhiteBlueGreenYellowRed"
  wkSpaceTime_Seas (data, spec_back,  diro, case, var             \
                ,latN, cpd, nDayWin, nDaySkip, opt  )
end

undef ("wkSpaceTime_obs_S2S")
procedure wkSpaceTime_obs_S2S(\
  data[*][*][*]     : numeric,\
  data_clm[*][*][*] : numeric,\
  seas              : string, \
  var               : string, \
  diro              : string, \
  case              : string, \
  latN              : numeric,\
  cpd               : integer,\
  nDayWin           : integer,\
  nDaySkip          : integer,\
  opt \
)

begin
  ; ; ; print(date_s2s)

  seasname   = seas

  if (isatt(opt, "season")) then
    iseason    = choose_season(seasname, date)
    data_back  = data_clm(iseason,:,:)
  else
    data_back  = data_clm
  end if

  
  opt_back   = False
  spec_back  = SpaceTime_Background(data_back, diro, case, var             \
              ,latN, cpd, nDayWin, nDaySkip, opt_back  )
  ;
  ; opt@Fig_3a      = log10(/1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 /)
  ; opt@Fig_3b      = log10(/1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 /) ; unequal
  ; opt@pltColorMap = "WhiteBlueGreenYellowRed"
  wkSpaceTime_Seas (data, spec_back,  diro, case, var             \
                ,latN, cpd, nDayWin, nDaySkip, opt  )
end


undef ("end_padding_observation")
function end_padding_observation(\
  data[*][*][*]     : numeric, \
  data_obs[*][*][*] : numeric, \
  npad              : integer
)

local time, date, time_add, ndims, ntim,nlat,nlon
begin
  data!0 = "time"
  time  = data&time
  ; printVarSummary(time)
  date  = cd_calendar(time,-2)

  time_add = new((/dimsizes(date)+npad/), typeof(time))

  ; printVarSummary(time)
  ; printVarSummary(time_add)
  time_add(0:dimsizes(date)-1) = time

  do d = 0, npad-1
    time_add(dimsizes(date)+d) = time(dimsizes(date)-1) + d+1
  end do

  ; ; print("data awal="+cd_calendar(time(0),-2))
  ; ; print("data akhir forecast="+cd_calendar(time(dimsizes(date)-1),-2))
  ; ; print("data awal padding forecast="+cd_calendar(time_add(dimsizes(date)),-2))

  ; ; ; print(time_add)
  ndims = dimsizes(data)
  ntim  = ndims(0)
  nlat  = ndims(1)
  nlon  = ndims(2)

  x = new((/ntim+npad, nlat, nlon/), float)
  x(0:ntim-1,:,:)  = data

  time_obs     = data_obs&time 
  date_obs     = cd_calendar(time_obs, -2)
  datestrt_pad = cd_calendar(time_add(dimsizes(date)),-2)
  istrt_pad    = ind(date_obs.eq.datestrt_pad)

  x(ntim:ntim+npad-1,:,:) = data_obs(istrt_pad:istrt_pad+npad-1,:,:)

  x!0 = "time"
  x!1 = "lat"
  x!2 = "lon"

  x&time = time_add
  x&lat  = data&lat
  x&lon  = data&lon

  return(x)

end

undef ("end_padding_observation_ensemble")
function end_padding_observation_ensemble(\
data[*][*][*][*]  : numeric, \
data_obs[*][*][*] : numeric, \
npad              : integer
)
begin
  ndims = dimsizes(data)
  nens  = ndims(0)
  ntim  = ndims(1)
  nlat  = ndims(2)
  nlon  = ndims(3)

  x    = new((/nens, ntim+npad, nlat, nlon/), float)
  do e = 0, nens-1
     tmp = end_padding_observation(data(e,:,:,:), data_obs, npad)
    ;  ; printVarSummary((tmp))
     x(e,:,:,:) = tmp
  end do

  x!0 = "time_lag"
  x&time_lag = (/-2,-1,0/)
  x@time_lag = "lagged by 2, 1, and 0 days"

  return(x)

end


; undef ("observation_padding_taper")
; function observation_padding_taper(\
;   data[*][*][*]      : numeric,\
;   pad[*][*][*]       : numeric,\
;   ind_init           : integer,\
;   npad               : integer,\
;   days_to_taper      : integer \
; )
;
; local time, date, time_add, ndims, ntim,nlat,nlon
; begin
;
;   data!0 = "time"
;   time  = data&time
;   date  = cd_calendar(time,-2)
;
;   time_add = new((/dimsizes(date)+npad/), typeof(time))
;
;   do d = 0, npad-1
;     ; ; ; print(time(d))
;     time_add(d) = time(0) - npad + d
;   end do
;
;   time_add(npad:) = time
;   ; ; ; print(cd_calendar(time_add,-2))
;
;   ndims = dimsizes(data)
;   ntim  = ndims(0)
;   nlat  = ndims(1)
;   nlon  = ndims(2)
;
;   x = new((/ntim+npad, nlat, nlon/), float)
;
;   ; ; printVarSummary(pad((ind_init - npad):(ind_init-npad+days_to_taper)-1, :, :))
;
;   tapered = data_taper_n(pad((ind_init - npad):(ind_init-npad+days_to_taper+1), :, :), 0, days_to_taper )
;   ; ; printVarSummary(tapered)
;   x(0:days_to_taper-1,:,:) = tapered(0:days_to_taper-1,:,:)
;
;   ; printVarSummary(pad(ind_init-days_to_taper:ind_init-1,:,:))
;   ; printVarSummary(x(days_to_taper:npad-1, :, :))
;   x(days_to_taper:npad-1, :, :) = pad((ind_init-days_to_taper):ind_init-1,:,:)
;   x(npad:,:,:)    = data
;
;   x!0 = "time"
;   x!1 = "lat"
;   x!2 = "lon"
;
;   x&time = time_add
;   x&lat  = data&lat
;   x&lon  = data&lon
;
;   return(x)
; end
;
; undef ("observation_padding_taper")
; function observation_padding_taper(\
;   data[*][*][*]      : numeric,\
;   pad[*][*][*]       : numeric,\
;   ind_init           : integer,\
;   npad               : integer,\
;   days_to_taper      : integer \
; )
;
; local time, date, time_add, ndims, ntim,nlat,nlon
; begin
;
;   data!0 = "time"
;   time  = data&time
;   date  = cd_calendar(time,-2)
;
;   time_add = new((/dimsizes(date)+npad/), typeof(time))
;
;   do d = 0, npad-1
;     ; ; ; print(time(d))
;     time_add(d) = time(0) - npad + d
;   end do
;
;   time_add(npad:) = time
;   ; ; ; print(cd_calendar(time_add,-2))
;
;   ndims = dimsizes(data)
;   ntim  = ndims(0)
;   nlat  = ndims(1)
;   nlon  = ndims(2)
;
;   x = new((/ntim+npad, nlat, nlon/), float)
;
;   ; ; printVarSummary(pad((ind_init - npad):(ind_init-npad+days_to_taper)-1, :, :))
;
;   tapered = data_taper_n(pad((ind_init - npad):(ind_init-1), :, :), 0, days_to_taper )
;   ; printVarSummary(tapered)
;   ; x(0:days_to_taper-1,:,:) = tapered(0:days_to_taper-1,:,:)
;   ; ; printVarSummary(x(0:npad-1, :, :))
;   x(0:npad-1, :, :) = tapered(:,:,:)
;   x(npad:,:,:)    = data
;
;   x!0 = "time"
;   x!1 = "lat"
;   x!2 = "lon"
;
;   x&time = time_add
;   x&lat  = data&lat
;   x&lon  = data&lon
;
;   return(x)
; end
