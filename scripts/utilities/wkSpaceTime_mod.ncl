;==================================================================
; ====> Create Wheeler-Kiladis Space-Time  plots.  <====
;==================================================================


;==================================================================
; Note_1: The full logitudinal domain is used.
;         This means that every planetary
;         wavenumber will be represented.
; Note_2: Tapering in time is done to make the variable periodic.
;
; The calculations are also only made for the latitudes
; between '-latBound' and 'latBound'.
;
;********************   REFERENCES  *******************************
; Wheeler, M., G.N. Kiladis
;    Convectively Coupled Equatorial Waves:
;    Analysis of Clouds and Temperature in the
;    Wavenumber-Frequency Domain
;    J. Atmos. Sci., 1999,  56: 374-399.
;---
; Hayashi, Y.
;    A Generalized Method of Resolving Disturbances into
;    Progressive and Retrogressive Waves by Space and
;    Fourier and TimeCross Spectral Analysis
;    J. Meteor. Soc. Japan, 1971, 49: 125-128.
;==================================================================




load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"        ; These three libraries are automatically
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"         ; loaded from NCL V6.2.0 onward.
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"     ; No need for user to explicitly load.

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl" ; Still need to load this library.


undef ("modified_daniell_kernel")
function modified_daniell_kernel(m:integer)
begin
  ; local w

  if m.eq.0 then
    wgt = 1.0
  else
  wgt = new((/m*2+1/), float)
  nweight = dimsizes(wgt)

  wnormal = 1./(2.*m)

  wstart = 0.5*wnormal

  wgt(0) = wstart

  do i = 1, nweight-2
    wgt(i) = wnormal
  end do
  wgt(nweight-1) = wgt(0)
  end if

  return wgt
end


; Smoothing 1-2-1 for coherence
undef ("smooth_121")
function smooth_121(data[*][*]:numeric, dim:integer)
begin

  wgt_first = (/0.75,0.25/)
  wgt_end   = (/0.25,0.75/)

  wgt_middle = (/0.25,0.5,0.25/)

  data_weighted = wgt_runave_n_Wrap(data, wgt_middle, 0, dim)
  ndims = dimsizes(data)
  ndata = ndims(dim)

  if dim.eq.0 then
    data_weighted(0,:) = (data(0,:) * wgt_first(0) + data(1,:) * wgt_first(1))
    data_weighted(ndata-1,:) = (data(ndata-2,:) * wgt_end(0) + data(ndata-1,:) * wgt_end(1))
  else
    data_weighted(:,0) = (data(:,0) * wgt_first(0) + data(:,1) * wgt_first(1))/4.
    data_weighted(:,ndata-1) = (data(:,ndata-2) * wgt_end(0) + data(:,ndata-1,:) * wgt_end(1))
  end if

  return(data_weighted)

end


undef("addHorVertLines_panel")
function addHorVertLines_panel(wks[1]:graphic, plot[1]:graphic \
                         ,x1[*],x2[*],y1[*],y2[*],y3[*],y4[*] )
; freq [y] axis: Add horizontal lines that explicitly
; print time in days. This assumes the units
; of the freq axis are "cpd" [cycles per day]
local gsres, txres
begin
    dom = new(7, graphic)
    gsres = True
    gsres@gsLineDashPattern = 1

    dom(0) = gsn_add_polyline(wks, plot, x1,y1,gsres)
    dom(1) = gsn_add_polyline(wks, plot, x1,y2,gsres)
    dom(2) = gsn_add_polyline(wks, plot, x1,y3,gsres)
    dom(3) = gsn_add_polyline(wks, plot, x2,y4,gsres)

    txres = True
    txres@txJust = "CenterLeft"
    txres@txFontHeightF = 0.013

    dom(4) = gsn_add_text(wks, plot, "3 days" ,-14.7,.345,txres)
    dom(5) = gsn_add_text(wks, plot, "6 days" ,-14.7,.175,txres)
    dom(6) = gsn_add_text(wks, plot, "30 days",-14.7,.045,txres)

    return (dom)
end


undef ("wkSpaceTime_Seas")
procedure wkSpaceTime_Seas (x[*][*][*]:numeric                  \
                      ,xback                               \
                      ,diro[1]:string                      \
                      ,caseName[1]:string                  \
                      ,varName[1]:string                   \
                      ,latBound[1]:numeric                 \
                      ,spd[1]:integer                      \
                      ,nDayWin[1]:integer                  \
                      ,nDaySkip[1]:integer                 \
                      ,opt[1]:logical)

local latN, latS, lonL, lonR, spd, fCrit, tim_taper        \
    , lon_taper, pltType, debug, minwav4smth, maxfrq4plt   \
    , minfrq4plt, maxfrq4plt, minwav4plt, maxwav4plt       \
    , fillVal, nMsg, dimx, ntim, nlat, mlon, nDayTot       \
    , nSampTot, nSampWin, nSampSkip, nWindow, N, dNam, work\
    , rmvMeans, xAS, q, peeAS, nl, ntStrt, ntLast, nw, nt  \
    , ml, psumanti, psumsym, wv, wkdir, caseName, pltFilTit\
    , frqfftwin, wavep1, minfrq, maxfrq, tmFontHgtF, pltTit\
    , tiFontHgtF, lbFontHgtF, txFontHgtF, res, freq        \
    , wavenumber, NWVN, dcres, txres, rlat, Ahe            \
    , nWaveType, nPlanetaryWave, nEquivDepth, Apzwn, Afreq \
    , asym, sym, x1, x2, y1, y2, y3, y4, wks, plot         \
    , psumb, psumsym_nolog, psumanti_nolog,  tt, smthlen   \
    , i, pt8cpd, spec, nCn, nExtra

  begin

  ; if seas.eq."NHSummer" then
  ; else if seas.eq."NHWinter" then
  ; else if seas.eq."NAutumn" then
  ; else if seas.eq. "NHSpring" then
  ; end if
  ; end if
  ; end if
  ; end if
  ;
  if (typeof(x).eq."float" .or. typeof(x).eq."double") then
      print(" ")
  else
      print("wkSpaceTime: input variable must be type float or double")
      print("             input variable type is "+typeof(x))
      exit
  end if

  debug   = False    ; default
  if (opt .and. isatt(opt, "debug") ) then
      debug = opt@debug
  end if

  if (isatt(x,"_FillValue")) then   ; Check for _FillValue .... not allowed
      nMsg = num(ismissing(x))
      if (nMsg.gt.0) then
          print("nMsg="+nMsg+"  User must preprocess to remove _FillValue")
          print("               FFTs do not allow missing values!!       ")
          exit
      end if
      delete(x@_FillValue)     ; avoid warning messages from fft
  end if

  if (debug) then
      printVarSummary( x )
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; x sizes and dimension names
;-------------------------------------------------------------------

  dimx = dimsizes(x)
  ntim = dimx(0)               ; total number of temporal samples
  nlat = dimx(1)
  mlon = dimx(2)

  dNam = getvardims( x )

;-------------------------------------------------------------------
; check to make sure that "x" has full days of data
;-------------------------------------------------------------------

  if ((ntim%spd).ne.0) then
      nExtra = ntim%spd
      print("nExtra="+nExtra+" input array must have complete days only ")
      exit
  end if

;-------------------------------------------------------------------
; Make input arguments into "internal" variables
;-------------------------------------------------------------------
  latN    = latBound
  latS    =-latBound ; make symmetric about the equator

  lonL    = 0        ; -180
  lonR    = 360      ;  180

  fCrit   = 1./nDayWin  ; remove all contributions 'longer'

  tim_taper  = 0.1   ; time taper      [0.1   => 10%]
  lon_taper  = 0.0   ; longitude taper [0.0 for globe; only global supported]

;-------------------------------------------------------------------
; Check for not allowed actions
;-------------------------------------------------------------------

  if (lon_taper.gt.0.0 .or. (lonR-lonL).ne.360.) then
      print("Code does currently allow lon_taper>0 or (lonR-lonL)<360")
      exit
  end if

;-------------------------------------------------------------------
; OPTIONS
;-------------------------------------------------------------------
  pltType = "ps"     ; default
  if (opt .and. isatt(opt, "pltType") ) then
      if (any(opt@pltType.eq.(/"ps", "eps", "x11", "ncgm", "png"/))) then
          pltType = opt@pltType
      end if
  end if

  pltTit   = caseName+"_"+varName
  pltTitle = pltTit+" LOG[Power: "+latBound+"S-"+latBound+"N]"
  if (opt .and. isatt(opt, "pltTitle") ) then
      pltTitle = opt@pltTitle
  end if

  pltFilTit = pltTit
  replaceChars( pltFilTit, " ", "_") ; spaces not allowed unix file names

  pltColorMap = "amwg_blueyellowred"
  if (opt .and. isatt(opt, "pltColorMap") ) then
      pltColorMap = opt@pltColorMap
  end if


;-------------------------------------------------------------------
; Create required temporal sampling variables
;-------------------------------------------------------------------

  nDayTot   = ntim/spd         ; # of days (total) for input variable
  nSampTot  = nDayTot*spd      ; # of samples (total)
  nSampWin  = nDayWin*spd      ; # of samples per temporal window
  nSampSkip = nDaySkip*spd     ; # of samples to skip between window segments
                               ;   neg means overlap
  nWindow   = (nSampTot-nSampWin)/(nSampWin+nSampSkip)  + 1
  N         =  nSampWin        ; convenience [historical]

  if (nDayTot.lt.nDayWin) then
      print("nDayTot="+nDayTot+" is less the nDayWin="+nDayWin)
      print("        This is not allowed !!       ")
      exit
  end if

;-------------------------------------------------------------------
; Remove dominant signals
; (a) Explicitly remove *long term* linear trend
;     For consistency with JET code keep the grid point means.
;     This necessitates that 'dtrend_msg' be used because 'dtrend'
;     always removes the mean(s).
; (b) All variations >= approx 'nDayWin' days if full year available
;-------------------------------------------------------------------

  dNam = getvardims( x )
;;work = x($dNam(1)$|:,$dNam(2)$|:,$dNam(0)$|:)    ; reorder  (lat,lon,time)
;;work = dtrend( work , False )   ; remove mean + overall long term temporal trend
;;work = dtrend_msg(ispan(0,ntim-1,1)*1.0, work, False, False)  ; remove just  trend
;;if (isatt(work,"_FillValue")) then
;;    delete(work@_FillValue)                  ; dtrend_msg adds this att
;;end if
                                                         ; replace with detrended
;;x    = (/ work($dNam(0)$|:,$dNam(1)$|:,$dNam(2)$|:) /) ; values (time,lat,lon)
;;delete(work)
  x    = dtrend_msg_n(ispan(0,ntim-1,1)*1.0, x, False, False, 0)  ; remove just trend

  if (nDayTot.ge.365) then                     ; rmv dominant signals
      rmvMeans = False                         ; original code did not remove
      x = rmvAnnualCycle(x, spd, nDayTot, rmvMeans, fCrit, 0)
  end if

  if (debug) then
      print("===> Post removal of trend and signal <===")
      printVarSummary( x )                     ; (time,lat,lon)
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; Decompose to Symmetric and Asymmetric parts
;-------------------------------------------------------------------
  xAS  = decompose2SymAsym( x , 1 )          ; create Asym and Sym parts
  if (debug) then
      printVarSummary(xAS)                   ; xAS(lat,lon,time) [iret=1]
      printMinMax(xAS, True)
  end if

;-------------------------------------------------------------------
; Because there is the possibility of overlapping *temporal* segments,
; we must use a less efficient approach and detrend/taper
; each window segment as it arises.
;          t0   t1   t2   t3   t4  .................. t(N)
; lon(0):  x00  x01  x02  x03  x04 .................. x0(N)
;     :    :   :   :   :   :                     :
; lon(M):  xM0  xM1  xM2  xM3  xM4 .................. xM(N)
;-------------------------------------------------------------------
; q     - temporary array to hold the 2D complex results
;         for each longitude/time (lon,time) window that is fft'd.
;         This is one instance [realization] of space-time decomposition.
;
; peeAS - symmetric and asymmetric power values in each latitude hemisphere.
;         Add extra lon/time to match JET
;-------------------------------------------------------------------
  q      = new((/2,mlon,nSampWin/)       ,typeof(xAS), "No_FillValue")
  peeAS  = new((/nlat,mlon+1,nSampWin+1/),typeof(xAS), "No_FillValue")
  peeAS  = 0.0                       ; initialize

;-------------------------------------------------------------------
; Operate on the xAS array
; NCL does not have a complex 2D FFT at this time.
; Perform a "poorman's" complex 2D FFT by looping over time and space.
; Loop over all latitudes and then perform summing/averaging
; on the spectral results.
;-------------------------------------------------------------------

  do nl=0,nlat-1

     ntStrt = 0
     ntLast = nSampWin-1
     if (debug) then
         print("==============> nl="+nl+" <==============")
     end if

    do nw=0,nWindow-1                        ; temporal window
       if (debug .and. nl.eq.0) then         ; debug
           print("nw="+nw+"  ntStrt="+ntStrt+"   ntLast="+ntLast)
       end if
       work = dtrend( xAS(:,:,ntStrt:ntLast), False )    ; detrend temporal window
       work = taper ( work, tim_taper, 0)    ; taper in "time" [periodic]
                                             ; work(nlat,mlon,N)
      do nt=0,nSampWin-1                     ; for each time perform complex fft in longitude
         q(:,:,nt) = cfftf( work(nl,:,nt), 0.0, 0)   ; space
      end do
       q = q/mlon                                    ; normalize by # lon samples

      do ml=0,mlon-1                         ; for each lon perform complex fft in time
         q(:,ml,:) = cfftf( q(0,ml,:), q(1,ml,:), 0) ; time
      end do
       q = q/nSampWin                                ; normalize by # time samples

;-------------------------------------------------------------------
; At this point 'q(2,mlon,nSampWin)' contains the
; real and imaginary space-time spectrum for this latitude
; ---
; Use Hayashi method to resolve into Progressive [Eastward]
;     and Retrogressive [Westward] Waves.
;-------------------------------------------------------------------

       pee = resolveWavesHayashi( q, nDayWin, spd )
       peeAS(nl,:,:) = peeAS(nl,:,:) + (pee/nWindow)   ; sum window contribution

       ntStrt = ntLast+nSampSkip+1           ; set index for next temporal window
       ntLast = ntStrt+nSampWin-1
    end do                                   ; nw (windows)

  end do                                     ; nl (lat)

  delete(work)
;-------------------------------------------------------------------
; Add meta data to the Hayashi space-time symmetric and asymmetric power
;-------------------------------------------------------------------

  peeAS!0     = "lat"
  peeAS!1     = "wave"
  peeAS!2     = "freq"
  peeAS&lat   = xAS&$dNam(1)$     ; nominally, xAS&lat
  peeAS&wave  = pee&wave
  peeAS&freq  = pee&freq
  peeAS@long_name =  "symmetric and asymmetric components"

  if (debug) then
      printVarSummary( peeAS )
      printMinMax( peeAS , True)
  end if

  delete(  q  )                              ; no longer needed
  delete( pee )

;-------------------------------------------------------------------
; now that we have the power array for sym and asym: use to
;    1) plot raw power spectrum (some smoothing)
;    2) derive and plot the background spectrum (lots of smoothing)
;    3) derive a denoised spectrum that is raw power/background power
;-------------------------------------------------------------------
; psumanti and psumsym will contain the symmetric and asymmetric power
; summed over latitude
;-------------------------------------------------------------------

  psumanti   = new((/dimsizes(peeAS&wave),dimsizes(peeAS&freq)/),typeof(peeAS), 1e20 )
  psumanti!0 = "wave"
  psumanti!1 = "freq"
  psumanti&wave = peeAS&wave
  psumanti&freq = peeAS&freq

  psumsym       = psumanti

  psumanti@long_name = "Asymmetric summed over lat"
  psumsym@long_name  = "Symmetric  summed over lat"

  if (nlat%2.eq.0)    ; use named dimensions to reorder
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2-1,:,:)   , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2-1)   )
  else
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2+1:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2,:,:)       , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2+1:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2)     )
  end if

;-------------------------------------------------------------------
; since summing over half the array (symmetric,asymmetric) the
; total variance is 2x the half sum
;-------------------------------------------------------------------
  psumanti = 2.0*psumanti
  psumsym  = 2.0*psumsym

;-------------------------------------------------------------------
; set the mean to missing to match original code
;-------------------------------------------------------------------
;
  psumanti(:,{0.0}) = (/ psumanti@_FillValue /)
  psumsym (:,{0.0}) = (/ psumsym@_FillValue /)

  if (debug) then
      printVarSummary( psumanti )         ; (wave,freq)
      printMinMax( psumsym , True)
  end if

;-------------------------------------------------------------------
; Apply smoothing to the spectrum. smooth over limited wave numbers
; Smoothing in frequency only (check if mean should be smoothed not smoothing now)
;--
; Smoothing parameters set these larger than the plotting
; wavenumbers to avoid smoothing artifacts
;-------------------------------------------------------------------
  minwav4smth = -27
  maxwav4smth =  27

  do wv=minwav4smth,maxwav4smth
     wk_smooth121( psumanti({wv},N/2+1:N-1) )
     wk_smooth121( psumsym ({wv},N/2+1:N-1) )
  end do

;-------------------------------------------------------------------
; Log10 scaling
;-------------------------------------------------------------------

  psumanti_nolog = psumanti
  psumsym_nolog  = psumsym

  psumanti       = log10(psumanti)
  psumsym        = log10(psumsym)

  psumanti@long_name = "log10(Asymmetric)"
  psumsym@long_name  = "log10(Symmetric)"

  if (debug) then
      printVarSummary( psumanti )             ; (wave,freq)
      printMinMax( psumanti , True)
      printVarSummary( psumsym )
      printMinMax( psumsym , True)
  end if

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        PLOT CODE follows
; --- set some 'plot variables
; set frequency maximum for plotting
; min(user specified frq, maxfrq in window)
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; The following allow DJS variable naming
; to be used with original plot code.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ;caseName       = case
  wkdir          = diro

  frqfftwin      = peeAS&freq
  frqfftwin&freq = peeAS&freq

  wavep1         = peeAS&wave
  wavep1&wave    = peeAS&wave

  if (debug) then
      printVarSummary( frqfftwin )
      printMinMax( frqfftwin, True )
      printVarSummary( wavep1 )
      printMinMax( wavep1, True )
  end if

;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
;-------------------------------------------------------------------
  minfrq4plt =  0.
  maxfrq4plt =  0.8
  minwav4plt = -15
  maxwav4plt =  15

  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))

  fillVal    = 1e20           ; miscellaneous

;=============================================================
;             Start Common Graphics Resources
;=============================================================
   tmFontHgtF            = 0.015     ; not sure why
   tiFontHgtF            = 0.018
   lbFontHgtF            = 0.015
   txFontHgtF            = 0.013

   res = True
   res@gsnFrame          = False
   res@gsnMaximize       = True
   res@gsnPaperOrientation = "portrait"

   res@gsnLeftString     = "Westward"
   res@gsnRightString    = "Eastward"

  ;res@lbBoxMinorExtentF = 0.18
   res@lbLabelFontHeightF= lbFontHgtF
   res@lbOrientation     = "vertical"

   res@cnFillOn          = True
   if (opt .and. isatt(opt, "cnLinesOn") \
           .and. .not.opt@cnLinesOn) then
       res@cnLinesOn     = False
   else
       res@cnLineThicknessF  = 0.5
   end if

   res@tmYLMode          = "Explicit"
   res@tmYLValues        = fspan(minfrq,maxfrq,9)
   res@tmYLLabels        = fspan(minfrq,maxfrq,9)
   res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)

   res@tmYLLabelFontHeightF = tmFontHgtF
   res@tmXBLabelFontHeightF = tmFontHgtF

   res@tiXAxisString     = "Zonal Wave Number"
   res@tiXAxisFontHeightF= tiFontHgtF

   res@tiYAxisString     = "Frequency (cpd)"
   res@tiYAxisFontHeightF= res@tiXAxisFontHeightF

   if (.not.(pltTitle.eq."" .or. pltTitle.eq." ")) then
       res@tiMainString      = pltTitle
       res@tiMainFontHeightF = tiFontHgtF
   end if
   res@txFontHeightF     = tiFontHgtF

;------------------------------------------------------
; Create a list of variable names that have predefined
; contour intervals.
;------------------------------------------------------
   varCnLevels = (/"FLUT" ,"OLR", "olr","U200","U850"   \
                  ,"PRECT","OMEGA500" /)

   if (any(varCnLevels.eq.varName)) then
       res@cnLevelSelectionMode = "ExplicitLevels"
   else
       res@cnLevelSelectionMode = "AutomaticLevels"
   end if

;-------------------------------
; horizontal dashed lines and text for frequency axis [plot only]
;-------------------------------

   freq       = frqfftwin({freq|minfrq:maxfrq})
   wavenumber = wavep1({wave|minwav4plt:maxwav4plt})
   NWVN       = dimsizes(wavenumber)         ; number of wavenumbers

   x1   = wavenumber
   x1!0 = "wave"
   y1   = new (NWVN,float)
   y1!0 = "freq"
   y1   = 1./3.           ; 3 days
   y2   = y1
   y2   = 1./6.           ; 6 days
   y3   = y1
   y3   = 1./30.          ; 30 days
   x2   = new(25,float)
   x2   = 0.0
   y4   = fspan(0.0,1.0,25)

;---------------------------------------------------------------
; dispersion: curves
;---------------------------------------------------------------

   rlat           = 0.0
   Ahe            = (/50.,25.,12./)
   nWaveType      = 6
   nPlanetaryWave = 50
   nEquivDepth    = dimsizes(Ahe)
   Apzwn          = new((/nWaveType,nEquivDepth,nPlanetaryWave/),"double",fillVal)
   Afreq          = Apzwn
   genDispersionCurves(nWaveType, nEquivDepth, nPlanetaryWave, rlat, Ahe, Afreq, Apzwn )

;---------------------------------------------------------------
; dispersion curve and text plot resources
;---------------------------------------------------------------
   dcres = True
   dcres@gsLineThicknessF  = 2.0
   dcres@gsLineDashPattern = 0

   txres = True
   txres@txJust        = "CenterLeft"
   txres@txPerimOn     = True
   txres@txFontHeightF = 0.013
   txres@txBackgroundFillColor = "Background"

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; plotting params for fig 1; subset for plot
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              ; reorder so freq is "y"
   asym       = psumanti({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   asym!0     = "freq"
   asym&freq  =  freq
   asym!1     = "wave"
   asym&wave  =  wavenumber
   asym@long_name = "Fig_1: log10(Asymmetric)"

   sym        = psumsym({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   sym@long_name  = "Fig_1: log10(Symmetric)"

   if (debug) then
       printVarSummary(asym)
       printMinMax(asym, True)
       printVarSummary(sym)
       printMinMax(sym, True)
   end if

;------------------------------------------------------
; Fig 1: Pre-defined contour levels [15]  for selected variables [non-log10]
;------------------------------------------------------
   nCn = 15

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2     \ ; unequal
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \ ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9, -4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_1")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_1                    ; user specified limits
   end if

;------------------------------------------------------
; Fig 1: ANTI-SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1a: ASYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Asym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)    ; not required

;------------------------------------------------------
; Fig 1: SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1b: SYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Sym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Symmetric"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       if (isatt(opt,"dirNetCDF")) then
           dirNetCDF = opt@dirNetCDF
       else
           dirNetCDF = "./"
       end if
       if (isatt(opt,"filNetCDF")) then
           filNetCDF = opt@filNetCDF
       else
           filNetCDF = "SpaceTime."+varName+".nc"
       end if
       fNam      = dirNetCDF+filNetCDF
       system ("/bin/rm -f "+fNam)

       ncdf      = addfile(fNam, "c")

       ncdf->FIG_1_SYM  =  sym
       ncdf->FIG_1_ASYM = asym
   end if

;-----------------------------------------------------------------------------
; ******  now derive and plot the background spectrum (red noise) ************
; [1] Sum power over all latitude
; [2] Put fill value in mean
; [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
;-----------------------------------------------------------------------------
  ;print("======> BACKGROUND <=====")

   psumb = dim_sum_n_Wrap(peeAS, 0)                  ; sum over all latitudes
 ;;psumb = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|:))  ; sum over all latitudes
   psumb@long_name = "Background Spectrum"

   psumb@_FillValue       = fillVal
   psumb(wave|:,freq|N/2) = fillVal
   psumb@_FillValue       = fillVal

   if (debug) then
       printVarSummary(psumb)             ; (wave,freq)
       printMinMax(psumb, True)
   end if

   do tt = N/2+1,N
      smthlen = maxwav4smth-minwav4smth+1
      if (frqfftwin(tt).lt.0.1) then
        do i = 1,5
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.1.and.frqfftwin(tt).lt.0.2) then
        do i = 1,10
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.2.and.frqfftwin(tt).lt.0.3) then
        do i = 1,20
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if

      if (frqfftwin(tt).ge.0.3) then
        do i = 1,40
          wk_smooth121(psumb(freq|tt,{wave|minwav4smth:maxwav4smth}))
        end do
      end if
   end do

   do nw = minwav4smth,maxwav4smth ; smth frequency up to .8 cycles per day
      pt8cpd  = min((/closest_val(.8,frqfftwin),dimsizes(frqfftwin)-1/))
      smthlen = pt8cpd-(N/2+1)+1
     do i = 1,10
        wk_smooth121( psumb({nw},N/2+1:pt8cpd) )
     end do
   end do

;----------------------------------------------------------------
; [1] Put fill value in mean (again)
; [2] SAVE the background spectrum for plotting Fig 3
; [3] LOGARITHMIC SCALING for plotting the background spectrum
;----------------------------------------------------------------

   psumb(wave|:,freq|N/2) = fillVal
   psumb_nolog = xback
   ; psumb_nolog = psumb

   printVarSummary(xback)
   printVarSummary(psumb)
   ; delete(psumb)

   
   psumb       = log10(xback)            ; LOG10
   ; psumb = log10(psumb)
   psumb@long_name = "log10( background spec )"
   if (debug) then
       print(" ===> Post multiple smoothing by wk_smooth121 <===")
       printVarSummary(psumb_nolog)
       printMinMax(psumb_nolog, True)

       printVarSummary(psumb)
       printMinMax(psumb, True)
   end if

;----------------------------------------------------------------
; set up for plotting  [ subset of frequencies: not necessary]
;----------------------------------------------------------------

   freq      = frqfftwin({freq|minfrq:maxfrq})
   spec      = psumb({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})

   if (debug) then
       printVarSummary(spec)                 ; (freq,wave)
       printMinMax(spec, True)
   end if

;----------------------------------------------------------------
; Fig 2: Predefined explicit contour levels BACKGROUND spectrum  [LOG10]
;----------------------------------------------------------------

   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)       ; allow size to change
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2 \     ; unequal 15
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \   ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9,-4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_2")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_2                   ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig2."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)

   res@gsnCenterString = "Background Power"
   plot = gsn_csm_contour(wks,spec,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)

;*********************************************************************************
; Fig 3a, 3b:  psum_nolog/psumb_nolog  [ratio]
;********************************************************************************

   psumsym_nolog  = psumsym_nolog /psumb_nolog   ; (wave,freq)
   psumanti_nolog = psumanti_nolog/psumb_nolog

   if (debug) then
       printVarSummary(psumanti_nolog)           ; (freq,wave)
       printMinMax(psumanti_nolog, True)
       printVarSummary(psumsym_nolog)            ; (freq,wave)
       printMinMax(psumsym_nolog, True)
   end if

;-----------------------------------------------------------
; ANTI-SYMMETRIC: RATIO ( subset to plot )
;-----------------------------------------------------------

   asym = psumanti_nolog({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   sym  = psumsym_nolog({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})

   if (debug) then
       printVarSummary(asym)                    ; (freq,wave)
       printMinMax(asym, True)
       printVarSummary(sym)                     ; (freq,wave)
       printMinMax(sym, True)
   end if

   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)                     ; allow size to change
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = fspan(0.3, 1.7, nCn)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(0.4, 1.8, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(0.4, 1.8, nCn)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/0.6,0.7 ,0.8,0.9 ,1.0,1.1,1.15,1.2,1.25 \
                       ,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/0.6,0.7,0.8,0.9,1.0,1.1,1.15,1.2,1.25 \
                       ,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if

   if (opt .and. isatt(opt, "Fig_3a")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3a                  ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig3.Asym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric/Background"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

; draw dispersion line

   gsn_polyline(wks,plot,Apzwn(0,0,:),Afreq(0,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,1,:),Afreq(0,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,2,:),Afreq(0,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,0,:),Afreq(1,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,1,:),Afreq(1,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,2,:),Afreq(1,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,0,:),Afreq(2,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,1,:),Afreq(2,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,2,:),Afreq(2,2,:),dcres)

; dispersion labels

   gsn_text(wks,plot,"MRG",-10.0,.15,txres)
   gsn_text(wks,plot,"n=2 IG",-3.0,.58,txres)
   gsn_text(wks,plot,"n=0 EIG",9.5,.50,txres)
   gsn_text(wks,plot,"h=50",-10.0,.78,txres)
   gsn_text(wks,plot,"h=25",-10.0,.63,txres)
   gsn_text(wks,plot,"h=12",-10.0,.51,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)

   frame(wks)
   delete(wks)             ; not required

;------------------------------------------------------------------
; SYMMETRIC
;------------------------------------------------------------------
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3b")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3b                  ; user specified
   end if

   wks  = gsn_open_wks(pltType, wkdir+"/Fig3.Sym."+pltFilTit)
   gsn_define_colormap(wks ,pltColorMap)
   res@gsnCenterString = "Symmetric/Background"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",5.0,.06,txres)
   frame(wks)

   if (debug .or. (opt .and. isatt(opt,"Fig_3_statInfo") \
                       .and. opt@Fig_3_statInfo)       ) then
       statAsymSym( asym, "Fig_3a" )
       statAsymSym(  sym, "Fig_3b" )
   end if

   delete(wks)

   ;------------------------------------------------------------------
   ; ALL
   ;------------------------------------------------------------------

   symasym = sym+asym

   copy_VarMeta(sym,symasym)
   printVarSummary(sym)
   printVarSummary(symasym)
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3c")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3c                  ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig3.symAsym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Sym+Asym/Background"
   plot = gsn_csm_contour(wks,symasym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

; draw dispersion line
   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_polyline(wks,plot,Apzwn(0,0,:),Afreq(0,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,1,:),Afreq(0,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,2,:),Afreq(0,2,:),dcres)

   tdres = dcres
   tdres@gsLineColor = "red"

   ftd1 = 1/2.5
   ftd2 = 1/5.
   gsn_polyline(wks,plot,(/-15,-15/), (/ftd2,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-6,-6/), (/ftd2,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-15,-6/), (/ftd1,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-15,-6/), (/ftd2,ftd2/),tdres)

   ; fmjo1 = 1/20.
   ; fmjo2 = 1/100.
   ; gsn_polyline(wks,plot,(/1,1/), (/fmjo2,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/5,5/), (/fmjo2,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/1,5/), (/fmjo1,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/1,5/), (/fmjo2,fmjo2/),tdres)
   ;
   ; fmrg1 = 1/3.
   ; fmrg2 = 1/10.
   ; kmrg1 = -1
   ; kmrg2 = -10
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg1/), (/fmrg2,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg2,kmrg2/), (/fmrg2,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg2/), (/fmrg1,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg2/), (/fmrg2,fmrg2/),tdres)
   ;
   ;
   ; ferw1 = 1/9.
   ; ferw2 = 1/72.
   ; kerw1 = -1
   ; kerw2 = -10
   ; gsn_polyline(wks,plot,(/kerw1,kerw1/), (/ferw2,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw2,kerw2/), (/ferw2,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw1,kerw2/), (/ferw1,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw1,kerw2/), (/ferw2,ferw2/),tdres)
   ;
   ;
   ; fklv1 = 1/3.
   ; fklv2 = 1/10.
   ; kklv1 = 1
   ; kklv2 = 14
   ; gsn_polyline(wks,plot,(/kklv1,kklv1/), (/fklv2,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv2,kklv2/), (/fklv2,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv1,kklv2/), (/fklv1,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv1,kklv2/), (/fklv2,fklv2/),tdres)



   fer1 = 1/7.
   fer2 = 1/10.
   res_poly = True
   res_poly@gsFillColor = "black"
   gsn_polygon(wks,plot,(/12,15,15,12,12/),(/fer1,fer1,fer2,fer2,fer1/),res_poly)

   ; gsn_polyline(wks,plot,)
   ; gsn_polyline(wks,plot,Apzwn(1,0,:),Afreq(1,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,1,:),Afreq(1,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,2,:),Afreq(1,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,0,:),Afreq(2,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,1,:),Afreq(2,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,2,:),Afreq(2,2,:),dcres)

; dispersion labels

   gsn_text(wks,plot,"MRG",-10.0,.15,txres)
   ; gsn_text(wks,plot,"n=2 IG",-3.0,.58,txres)
   ; gsn_text(wks,plot,"n=0 EIG",9.5,.50,txres)
   ; gsn_text(wks,plot,"h=50",-10.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-10.0,.63,txres)
   ; gsn_text(wks,plot,"h=12",-10.0,.51,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   ; gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   ; gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   ; gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",4.0,.03,txres)

   frame(wks)
   delete(wks)             ; not required

;------------------------------------------------------------------
; SYMMETRIC
;------------------------------------------------------------------
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3b")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3b                  ; user specified
   end if

   wks  = gsn_open_wks(pltType, wkdir+"/Fig3.Sym."+pltFilTit)
   gsn_define_colormap(wks ,pltColorMap)
   res@gsnCenterString = "Symmetric/Background"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",5.0,.06,txres)
   frame(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       ncdf->FIG_3_BACK = spec
       ncdf->FIG_3_SYM  =  sym
       ncdf->FIG_3_ASYM = asym
       ncdf->FIG_3_ALL  = sym+asym
   end if

   if (debug) then
       print("********** FINISHED: Space-Time *****************")
   end if
end

;=============================================================================
undef ("SpaceTime_Background")
function SpaceTime_Background (x[*][*][*]:numeric                  \
                      ,diro[1]:string                      \
                      ,caseName[1]:string                  \
                      ,varName[1]:string                   \
                      ,latBound[1]:numeric                 \
                      ,spd[1]:integer                      \
                      ,nDayWin[1]:integer                  \
                      ,nDaySkip[1]:integer                 \
                      ,opt[1]:logical)

local latN, latS, lonL, lonR, spd, fCrit, tim_taper        \
    , lon_taper, pltType, debug, minwav4smth, maxfrq4plt   \
    , minfrq4plt, maxfrq4plt, minwav4plt, maxwav4plt       \
    , fillVal, nMsg, dimx, ntim, nlat, mlon, nDayTot       \
    , nSampTot, nSampWin, nSampSkip, nWindow, N, dNam, work\
    , rmvMeans, xAS, q, peeAS, nl, ntStrt, ntLast, nw, nt  \
    , ml, psumanti, psumsym, wv, wkdir, caseName, pltFilTit\
    , frqfftwin, wavep1, minfrq, maxfrq, tmFontHgtF, pltTit\
    , tiFontHgtF, lbFontHgtF, txFontHgtF, res, freq        \
    , wavenumber, NWVN, dcres, txres, rlat, Ahe            \
    , nWaveType, nPlanetaryWave, nEquivDepth, Apzwn, Afreq \
    , asym, sym, x1, x2, y1, y2, y3, y4, wks, plot         \
    , psumb, psumsym_nolog, psumanti_nolog,  tt, smthlen   \
    , i, pt8cpd, nCn, nExtra

  begin

  ; if seas.eq."NHSummer" then
  ; else if seas.eq."NHWinter" then
  ; else if seas.eq."NAutumn" then
  ; else if seas.eq. "NHSpring" then
  ; end if
  ; end if
  ; end if
  ; end if
  ;
  if (typeof(x).eq."float" .or. typeof(x).eq."double") then
      print(" ")
  else
      print("wkSpaceTime: input variable must be type float or double")
      print("             input variable type is "+typeof(x))
      exit
  end if

  debug   = False    ; default
  if (opt .and. isatt(opt, "debug") ) then
      debug = opt@debug
  end if

  if (isatt(x,"_FillValue")) then   ; Check for _FillValue .... not allowed
      nMsg = num(ismissing(x))
      if (nMsg.gt.0) then
          print("nMsg="+nMsg+"  User must preprocess to remove _FillValue")
          print("               FFTs do not allow missing values!!       ")
          exit
      end if
      delete(x@_FillValue)     ; avoid warning messages from fft
  end if

  if (debug) then
      printVarSummary( x )
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; x sizes and dimension names
;-------------------------------------------------------------------

  dimx = dimsizes(x)
  ntim = dimx(0)               ; total number of temporal samples
  nlat = dimx(1)
  mlon = dimx(2)

  dNam = getvardims( x )

;-------------------------------------------------------------------
; check to make sure that "x" has full days of data
;-------------------------------------------------------------------

  if ((ntim%spd).ne.0) then
      nExtra = ntim%spd
      print("nExtra="+nExtra+" input array must have complete days only ")
      exit
  end if

;-------------------------------------------------------------------
; Make input arguments into "internal" variables
;-------------------------------------------------------------------
  latN    = latBound
  latS    =-latBound ; make symmetric about the equator

  lonL    = 0        ; -180
  lonR    = 360      ;  180

  fCrit   = 1./nDayWin  ; remove all contributions 'longer'

  tim_taper  = 0.1   ; time taper      [0.1   => 10%]
  lon_taper  = 0.0   ; longitude taper [0.0 for globe; only global supported]

;-------------------------------------------------------------------
; Check for not allowed actions
;-------------------------------------------------------------------

  if (lon_taper.gt.0.0 .or. (lonR-lonL).ne.360.) then
      print("Code does currently allow lon_taper>0 or (lonR-lonL)<360")
      exit
  end if

;-------------------------------------------------------------------
; OPTIONS
;-------------------------------------------------------------------
  pltType = "ps"     ; default
  if (opt .and. isatt(opt, "pltType") ) then
      if (any(opt@pltType.eq.(/"ps", "eps", "x11", "ncgm", "png"/))) then
          pltType = opt@pltType
      end if
  end if

  ; pltTit   = caseName+"_"+varName
  pltTit   = caseName
  ; pltTitle = pltTit+" LOG[Power: "+latBound+"S-"+latBound+"N]"
  pltTitle = pltTit+" LOG[Power: "+latBound+"S-"+latBound+"N]"
  if (opt .and. isatt(opt, "pltTitle") ) then
      pltTitle = opt@pltTitle
  end if

  pltFilTit = pltTit
  replaceChars( pltFilTit, " ", "_") ; spaces not allowed unix file names

  pltColorMap = "amwg_blueyellowred"
  if (opt .and. isatt(opt, "pltColorMap") ) then
      pltColorMap = opt@pltColorMap
  end if


;-------------------------------------------------------------------
; Create required temporal sampling variables
;-------------------------------------------------------------------

  nDayTot   = ntim/spd         ; # of days (total) for input variable
  nSampTot  = nDayTot*spd      ; # of samples (total)
  nSampWin  = nDayWin*spd      ; # of samples per temporal window
  nSampSkip = nDaySkip*spd     ; # of samples to skip between window segments
                               ;   neg means overlap
  nWindow   = (nSampTot-nSampWin)/(nSampWin+nSampSkip)  + 1
  N         =  nSampWin        ; convenience [historical]

  if (nDayTot.lt.nDayWin) then
      print("nDayTot="+nDayTot+" is less the nDayWin="+nDayWin)
      print("        This is not allowed !!       ")
      exit
  end if

;-------------------------------------------------------------------
; Remove dominant signals
; (a) Explicitly remove *long term* linear trend
;     For consistency with JET code keep the grid point means.
;     This necessitates that 'dtrend_msg' be used because 'dtrend'
;     always removes the mean(s).
; (b) All variations >= approx 'nDayWin' days if full year available
;-------------------------------------------------------------------

  dNam = getvardims( x )
;;work = x($dNam(1)$|:,$dNam(2)$|:,$dNam(0)$|:)    ; reorder  (lat,lon,time)
;;work = dtrend( work , False )   ; remove mean + overall long term temporal trend
;;work = dtrend_msg(ispan(0,ntim-1,1)*1.0, work, False, False)  ; remove just  trend
;;if (isatt(work,"_FillValue")) then
;;    delete(work@_FillValue)                  ; dtrend_msg adds this att
;;end if
                                                         ; replace with detrended
;;x    = (/ work($dNam(0)$|:,$dNam(1)$|:,$dNam(2)$|:) /) ; values (time,lat,lon)
;;delete(work)
  x    = dtrend_msg_n(ispan(0,ntim-1,1)*1.0, x, False, False, 0)  ; remove just trend

  if (nDayTot.ge.365) then                     ; rmv dominant signals
      rmvMeans = False                         ; original code did not remove
      x = rmvAnnualCycle(x, spd, nDayTot, rmvMeans, fCrit, 0)
  end if

  if (debug) then
      print("===> Post removal of trend and signal <===")
      printVarSummary( x )                     ; (time,lat,lon)
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; Decompose to Symmetric and Asymmetric parts
;-------------------------------------------------------------------
  xAS  = decompose2SymAsym( x , 1 )          ; create Asym and Sym parts
  if (debug) then
      printVarSummary(xAS)                   ; xAS(lat,lon,time) [iret=1]
      printMinMax(xAS, True)
  end if

;-------------------------------------------------------------------
; Because there is the possibility of overlapping *temporal* segments,
; we must use a less efficient approach and detrend/taper
; each window segment as it arises.
;          t0   t1   t2   t3   t4  .................. t(N)
; lon(0):  x00  x01  x02  x03  x04 .................. x0(N)
;     :    :   :   :   :   :                     :
; lon(M):  xM0  xM1  xM2  xM3  xM4 .................. xM(N)
;-------------------------------------------------------------------
; q     - temporary array to hold the 2D complex results
;         for each longitude/time (lon,time) window that is fft'd.
;         This is one instance [realization] of space-time decomposition.
;
; peeAS - symmetric and asymmetric power values in each latitude hemisphere.
;         Add extra lon/time to match JET
;-------------------------------------------------------------------
  q      = new((/2,mlon,nSampWin/)       ,typeof(xAS), "No_FillValue")
  peeAS  = new((/nlat,mlon+1,nSampWin+1/),typeof(xAS), "No_FillValue")
  peeAS  = 0.0                       ; initialize

;-------------------------------------------------------------------
; Operate on the xAS array
; NCL does not have a complex 2D FFT at this time.
; Perform a "poorman's" complex 2D FFT by looping over time and space.
; Loop over all latitudes and then perform summing/averaging
; on the spectral results.
;-------------------------------------------------------------------

  do nl=0,nlat-1

     ntStrt = 0
     ntLast = nSampWin-1
     if (debug) then
         print("==============> nl="+nl+" <==============")
     end if

    do nw=0,nWindow-1                        ; temporal window
       if (debug .and. nl.eq.0) then         ; debug
           print("nw="+nw+"  ntStrt="+ntStrt+"   ntLast="+ntLast)
       end if
       work = dtrend( xAS(:,:,ntStrt:ntLast), False )    ; detrend temporal window
       work = taper ( work, tim_taper, 0)    ; taper in "time" [periodic]
                                             ; work(nlat,mlon,N)
      do nt=0,nSampWin-1                     ; for each time perform complex fft in longitude
         q(:,:,nt) = cfftf( work(nl,:,nt), 0.0, 0)   ; space
      end do
       q = q/mlon                                    ; normalize by # lon samples

      do ml=0,mlon-1                         ; for each lon perform complex fft in time
         q(:,ml,:) = cfftf( q(0,ml,:), q(1,ml,:), 0) ; time
      end do
       q = q/nSampWin                                ; normalize by # time samples

;-------------------------------------------------------------------
; At this point 'q(2,mlon,nSampWin)' contains the
; real and imaginary space-time spectrum for this latitude
; ---
; Use Hayashi method to resolve into Progressive [Eastward]
;     and Retrogressive [Westward] Waves.
;-------------------------------------------------------------------

       pee = resolveWavesHayashi( q, nDayWin, spd )
       peeAS(nl,:,:) = peeAS(nl,:,:) + (pee/nWindow)   ; sum window contribution

       ntStrt = ntLast+nSampSkip+1           ; set index for next temporal window
       ntLast = ntStrt+nSampWin-1
    end do                                   ; nw (windows)

  end do                                     ; nl (lat)

  delete(work)
;-------------------------------------------------------------------
; Add meta data to the Hayashi space-time symmetric and asymmetric power
;-------------------------------------------------------------------

  peeAS!0     = "lat"
  peeAS!1     = "wave"
  peeAS!2     = "freq"
  peeAS&lat   = xAS&$dNam(1)$     ; nominally, xAS&lat
  peeAS&wave  = pee&wave
  peeAS&freq  = pee&freq
  peeAS@long_name =  "symmetric and asymmetric components"

  if (debug) then
      printVarSummary( peeAS )
      printMinMax( peeAS , True)
  end if

  delete(  q  )                              ; no longer needed
  delete( pee )

;-------------------------------------------------------------------
; now that we have the power array for sym and asym: use to
;    1) plot raw power spectrum (some smoothing)
;    2) derive and plot the background spectrum (lots of smoothing)
;    3) derive a denoised spectrum that is raw power/background power
;-------------------------------------------------------------------
; psumanti and psumsym will contain the symmetric and asymmetric power
; summed over latitude
;-------------------------------------------------------------------

  psumanti   = new((/dimsizes(peeAS&wave),dimsizes(peeAS&freq)/),typeof(peeAS), 1e20 )
  psumanti!0 = "wave"
  psumanti!1 = "freq"
  psumanti&wave = peeAS&wave
  psumanti&freq = peeAS&freq

  psumsym       = psumanti

  psumanti@long_name = "Asymmetric summed over lat"
  psumsym@long_name  = "Symmetric  summed over lat"

  if (nlat%2.eq.0)    ; use named dimensions to reorder
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2-1,:,:)   , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2-1)   )
  else
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2+1:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2,:,:)       , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2+1:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2)     )
  end if

;-------------------------------------------------------------------
; since summing over half the array (symmetric,asymmetric) the
; total variance is 2x the half sum
;-------------------------------------------------------------------
  psumanti = 2.0*psumanti
  psumsym  = 2.0*psumsym

;-------------------------------------------------------------------
; set the mean to missing to match original code
;-------------------------------------------------------------------
;
  psumanti(:,{0.0}) = (/ psumanti@_FillValue /)
  psumsym (:,{0.0}) = (/ psumsym@_FillValue /)

  if (debug) then
      printVarSummary( psumanti )         ; (wave,freq)
      printMinMax( psumsym , True)
  end if

;-------------------------------------------------------------------
; Apply smoothing to the spectrum. smooth over limited wave numbers
; Smoothing in frequency only (check if mean should be smoothed not smoothing now)
;--
; Smoothing parameters set these larger than the plotting
; wavenumbers to avoid smoothing artifacts
;-------------------------------------------------------------------
  minwav4smth = -27
  maxwav4smth =  27

  do wv=minwav4smth,maxwav4smth
     wk_smooth121( psumanti({wv},N/2+1:N-1) )
     wk_smooth121( psumsym ({wv},N/2+1:N-1) )
  end do

;-------------------------------------------------------------------
; Log10 scaling
;-------------------------------------------------------------------

  psumanti_nolog = psumanti
  psumsym_nolog  = psumsym

  psumanti       = log10(psumanti)
  psumsym        = log10(psumsym)

  psumanti@long_name = "log10(Asymmetric)"
  psumsym@long_name  = "log10(Symmetric)"

  if (debug) then
      printVarSummary( psumanti )             ; (wave,freq)
      printMinMax( psumanti , True)
      printVarSummary( psumsym )
      printMinMax( psumsym , True)
  end if

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        PLOT CODE follows
; --- set some 'plot variables
; set frequency maximum for plotting
; min(user specified frq, maxfrq in window)
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; The following allow DJS variable naming
; to be used with original plot code.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ;caseName       = case
  wkdir          = diro

  frqfftwin      = peeAS&freq
  frqfftwin&freq = peeAS&freq

  wavep1         = peeAS&wave
  wavep1&wave    = peeAS&wave

  if (debug) then
      printVarSummary( frqfftwin )
      printMinMax( frqfftwin, True )
      printVarSummary( wavep1 )
      printMinMax( wavep1, True )
  end if

;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
;-------------------------------------------------------------------
  minfrq4plt =  0.
  maxfrq4plt =  0.8
  minwav4plt = -15
  maxwav4plt =  15

  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))

  fillVal    = 1e20           ; miscellaneous

;=============================================================
;             Start Common Graphics Resources
;=============================================================
   tmFontHgtF            = 0.015     ; not sure why
   tiFontHgtF            = 0.018
   lbFontHgtF            = 0.015
   txFontHgtF            = 0.013

   res = True
   res@gsnFrame          = False
   res@gsnMaximize       = True
   res@gsnPaperOrientation = "portrait"

   res@gsnLeftString     = "Westward"
   res@gsnRightString    = "Eastward"

  ;res@lbBoxMinorExtentF = 0.18
   res@lbLabelFontHeightF= lbFontHgtF
   res@lbOrientation     = "vertical"

   res@cnFillOn          = True
   if (opt .and. isatt(opt, "cnLinesOn") \
           .and. .not.opt@cnLinesOn) then
       res@cnLinesOn     = False
   else
       res@cnLineThicknessF  = 0.5
   end if

   res@tmYLMode          = "Explicit"
   res@tmYLValues        = fspan(minfrq,maxfrq,9)
   res@tmYLLabels        = fspan(minfrq,maxfrq,9)
   res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)

   res@tmYLLabelFontHeightF = tmFontHgtF
   res@tmXBLabelFontHeightF = tmFontHgtF

   res@tiXAxisString     = "Zonal Wave Number"
   res@tiXAxisFontHeightF= tiFontHgtF

   res@tiYAxisString     = "Frequency (cpd)"
   res@tiYAxisFontHeightF= res@tiXAxisFontHeightF

   if (.not.(pltTitle.eq."" .or. pltTitle.eq." ")) then
       res@tiMainString      = pltTitle
       res@tiMainFontHeightF = tiFontHgtF
   end if
   res@txFontHeightF     = tiFontHgtF

;------------------------------------------------------
; Create a list of variable names that have predefined
; contour intervals.
;------------------------------------------------------
   varCnLevels = (/"FLUT" ,"OLR", "olr","U200","U850"   \
                  ,"PRECT","OMEGA500" /)

   if (any(varCnLevels.eq.varName)) then
       res@cnLevelSelectionMode = "ExplicitLevels"
   else
       res@cnLevelSelectionMode = "AutomaticLevels"
   end if

;-------------------------------
; horizontal dashed lines and text for frequency axis [plot only]
;-------------------------------

   freq       = frqfftwin({freq|minfrq:maxfrq})
   wavenumber = wavep1({wave|minwav4plt:maxwav4plt})
   NWVN       = dimsizes(wavenumber)         ; number of wavenumbers

   x1   = wavenumber
   x1!0 = "wave"
   y1   = new (NWVN,float)
   y1!0 = "freq"
   y1   = 1./3.           ; 3 days
   y2   = y1
   y2   = 1./6.           ; 6 days
   y3   = y1
   y3   = 1./30.          ; 30 days
   x2   = new(25,float)
   x2   = 0.0
   y4   = fspan(0.0,1.0,25)

;---------------------------------------------------------------
; dispersion: curves
;---------------------------------------------------------------

   rlat           = 0.0
   Ahe            = (/50.,25.,12./)
   nWaveType      = 6
   nPlanetaryWave = 50
   nEquivDepth    = dimsizes(Ahe)
   Apzwn          = new((/nWaveType,nEquivDepth,nPlanetaryWave/),"double",fillVal)
   Afreq          = Apzwn
   genDispersionCurves(nWaveType, nEquivDepth, nPlanetaryWave, rlat, Ahe, Afreq, Apzwn )

;---------------------------------------------------------------
; dispersion curve and text plot resources
;---------------------------------------------------------------
   dcres = True
   dcres@gsLineThicknessF  = 2.0
   dcres@gsLineDashPattern = 0

   txres = True
   txres@txJust        = "CenterLeft"
   txres@txPerimOn     = True
   txres@txFontHeightF = 0.013
   txres@txBackgroundFillColor = "Background"

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; plotting params for fig 1; subset for plot
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              ; reorder so freq is "y"
   asym       = psumanti({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   asym!0     = "freq"
   asym&freq  =  freq
   asym!1     = "wave"
   asym&wave  =  wavenumber
   asym@long_name = "Fig_1: log10(Asymmetric)"

   sym        = psumsym({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   sym@long_name  = "Fig_1: log10(Symmetric)"

   if (debug) then
       printVarSummary(asym)
       printMinMax(asym, True)
       printVarSummary(sym)
       printMinMax(sym, True)
   end if

;------------------------------------------------------
; Fig 1: Pre-defined contour levels [15]  for selected variables [non-log10]
;------------------------------------------------------
   nCn = 15

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2     \ ; unequal
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \ ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9, -4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_1")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_1                    ; user specified limits
   end if

;------------------------------------------------------
; Fig 1: ANTI-SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1a: ASYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Asym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)    ; not required

;------------------------------------------------------
; Fig 1: SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1b: SYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Sym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Symmetric"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       if (isatt(opt,"dirNetCDF")) then
           dirNetCDF = opt@dirNetCDF
       else
           dirNetCDF = "./"
       end if
       if (isatt(opt,"filNetCDF")) then
           filNetCDF = opt@filNetCDF
       else
           filNetCDF = "SpaceTime."+varName+".nc"
       end if
       fNam      = dirNetCDF+filNetCDF
       system ("/bin/rm -f "+fNam)

       ncdf      = addfile(fNam, "c")

       ncdf->FIG_1_SYM  =  sym
       ncdf->FIG_1_ASYM = asym
   end if

;-----------------------------------------------------------------------------
; ******  now derive and plot the background spectrum (red noise) ************
; [1] Sum power over all latitude
; [2] Put fill value in mean
; [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
;-----------------------------------------------------------------------------
  ;print("======> BACKGROUND <=====")

   psumb = dim_sum_n_Wrap(peeAS, 0)                  ; sum over all latitudes
 ;;psumb = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|:))  ; sum over all latitudes
   psumb@long_name = "Background Spectrum"

   psumb@_FillValue       = fillVal
   psumb(wave|:,freq|N/2) = fillVal
   psumb@_FillValue       = fillVal

   if (debug) then
       printVarSummary(psumb)             ; (wave,freq)
       printMinMax(psumb, True)
   end if

   do tt = N/2+1,N
      smthlen = maxwav4smth-minwav4smth+1
      if (frqfftwin(tt).lt.0.1) then
        do i = 1,5
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.1.and.frqfftwin(tt).lt.0.2) then
        do i = 1,10
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.2.and.frqfftwin(tt).lt.0.3) then
        do i = 1,20
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if

      if (frqfftwin(tt).ge.0.3) then
        do i = 1,40
          wk_smooth121(psumb(freq|tt,{wave|minwav4smth:maxwav4smth}))
        end do
      end if
   end do

   do nw = minwav4smth,maxwav4smth ; smth frequency up to .8 cycles per day
      pt8cpd  = min((/closest_val(.8,frqfftwin),dimsizes(frqfftwin)-1/))
      smthlen = pt8cpd-(N/2+1)+1
     do i = 1,10
        wk_smooth121( psumb({nw},N/2+1:pt8cpd) )
     end do
   end do

;----------------------------------------------------------------
; [1] Put fill value in mean (again)
; [2] SAVE the background spectrum for plotting Fig 3
; [3] LOGARITHMIC SCALING for plotting the background spectrum
;----------------------------------------------------------------

   psumb(wave|:,freq|N/2) = fillVal
   psumb_nolog = psumb
   psumb       = log10(psumb)            ; LOG10
   psumb@long_name = "log10( background spec )"
   if (debug) then
       print(" ===> Post multiple smoothing by wk_smooth121 <===")
       printVarSummary(psumb_nolog)
       printMinMax(psumb_nolog, True)

       printVarSummary(psumb)
       printMinMax(psumb, True)
   end if

;----------------------------------------------------------------
; set up for plotting  [ subset of frequencies: not necessary]
;----------------------------------------------------------------

   freq      = frqfftwin({freq|minfrq:maxfrq})
   spec      = psumb({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   return(psumb_nolog)

end


;''


undef ("FR04SpaceTime_Seas")
procedure FR04SpaceTime_Seas (x[*][*][*]:numeric                  \
                      ,xback                               \
                      ,diro[1]:string                      \
                      ,caseName[1]:string                  \
                      ,varName[1]:string                   \
                      ,latBound[1]:numeric                 \
                      ,spd[1]:integer                      \
                      ,nDayWin[1]:integer                  \
                      ,nDaySkip[1]:integer                 \
                      ,mode[1]:string                      \
                      ,opt[1]:logical)

local latN, latS, lonL, lonR, spd, fCrit, tim_taper        \
    , lon_taper, pltType, debug, minwav4smth, maxfrq4plt   \
    , minfrq4plt, maxfrq4plt, minwav4plt, maxwav4plt       \
    , fillVal, nMsg, dimx, ntim, nlat, mlon, nDayTot       \
    , nSampTot, nSampWin, nSampSkip, nWindow, N, dNam, work\
    , rmvMeans, xAS, q, peeAS, nl, ntStrt, ntLast, nw, nt  \
    , ml, psumanti, psumsym, wv, wkdir, caseName, pltFilTit\
    , frqfftwin, wavep1, minfrq, maxfrq, tmFontHgtF, pltTit\
    , tiFontHgtF, lbFontHgtF, txFontHgtF, res, freq        \
    , wavenumber, NWVN, dcres, txres, rlat, Ahe            \
    , nWaveType, nPlanetaryWave, nEquivDepth, Apzwn, Afreq \
    , asym, sym, x1, x2, y1, y2, y3, y4, wks, plot         \
    , psumb, psumsym_nolog, psumanti_nolog,  tt, smthlen   \
    , i, pt8cpd, spec, nCn, nExtra

  begin

; mode "NH" or "SH"


  ; if seas.eq."NHSummer" then
  ; else if seas.eq."NHWinter" then
  ; else if seas.eq."NAutumn" then
  ; else if seas.eq. "NHSpring" then
  ; end if
  ; end if
  ; end if
  ; end if
  ;
  if (typeof(x).eq."float" .or. typeof(x).eq."double") then
      print(" ")
  else
      print("wkSpaceTime: input variable must be type float or double")
      print("             input variable type is "+typeof(x))
      exit
  end if

  debug   = False    ; default
  if (opt .and. isatt(opt, "debug") ) then
      debug = opt@debug
  end if

  if (isatt(x,"_FillValue")) then   ; Check for _FillValue .... not allowed
      nMsg = num(ismissing(x))
      if (nMsg.gt.0) then
          print("nMsg="+nMsg+"  User must preprocess to remove _FillValue")
          print("               FFTs do not allow missing values!!       ")
          exit
      end if
      delete(x@_FillValue)     ; avoid warning messages from fft
  end if

  if (debug) then
      printVarSummary( x )
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; x sizes and dimension names
;-------------------------------------------------------------------

  dimx = dimsizes(x)
  ntim = dimx(0)               ; total number of temporal samples
  nlat = dimx(1)
  mlon = dimx(2)

  dNam = getvardims( x )

;-------------------------------------------------------------------
; check to make sure that "x" has full days of data
;-------------------------------------------------------------------

  if ((ntim%spd).ne.0) then
      nExtra = ntim%spd
      print("nExtra="+nExtra+" input array must have complete days only ")
      exit
  end if

;-------------------------------------------------------------------
; Make input arguments into "internal" variables
;-------------------------------------------------------------------
  latN    = latBound
  latS    =-latBound ; make symmetric about the equator

  lonL    = 0        ; -180
  lonR    = 360      ;  180

  fCrit   = 1./nDayWin  ; remove all contributions 'longer'

  tim_taper  = 0.1   ; time taper      [0.1   => 10%]
  lon_taper  = 0.0   ; longitude taper [0.0 for globe; only global supported]

;-------------------------------------------------------------------
; Check for not allowed actions
;-------------------------------------------------------------------

  if (lon_taper.gt.0.0 .or. (lonR-lonL).ne.360.) then
      print("Code does currently allow lon_taper>0 or (lonR-lonL)<360")
      exit
  end if

;-------------------------------------------------------------------
; OPTIONS
;-------------------------------------------------------------------
  pltType = "ps"     ; default
  if (opt .and. isatt(opt, "pltType") ) then
      if (any(opt@pltType.eq.(/"ps", "eps", "x11", "ncgm", "png"/))) then
          pltType = opt@pltType
      end if
  end if

  pltTit   = caseName+"_"+varName
  pltTitle = pltTit+" LOG[Power: "+latBound+"S-"+latBound+"N]"
  if (opt .and. isatt(opt, "pltTitle") ) then
      pltTitle = opt@pltTitle
  end if

  pltFilTit = pltTit
  replaceChars( pltFilTit, " ", "_") ; spaces not allowed unix file names

  pltColorMap = "amwg_blueyellowred"
  if (opt .and. isatt(opt, "pltColorMap") ) then
      pltColorMap = opt@pltColorMap
  end if


;-------------------------------------------------------------------
; Create required temporal sampling variables
;-------------------------------------------------------------------

  nDayTot   = ntim/spd         ; # of days (total) for input variable
  nSampTot  = nDayTot*spd      ; # of samples (total)
  nSampWin  = nDayWin*spd      ; # of samples per temporal window
  nSampSkip = nDaySkip*spd     ; # of samples to skip between window segments
                               ;   neg means overlap
  nWindow   = (nSampTot-nSampWin)/(nSampWin+nSampSkip)  + 1
  N         =  nSampWin        ; convenience [historical]

  if (nDayTot.lt.nDayWin) then
      print("nDayTot="+nDayTot+" is less the nDayWin="+nDayWin)
      print("        This is not allowed !!       ")
      exit
  end if

;-------------------------------------------------------------------
; Remove dominant signals
; (a) Explicitly remove *long term* linear trend
;     For consistency with JET code keep the grid point means.
;     This necessitates that 'dtrend_msg' be used because 'dtrend'
;     always removes the mean(s).
; (b) All variations >= approx 'nDayWin' days if full year available
;-------------------------------------------------------------------

; mlon = 144
  dNam = getvardims( x )
;;work = x($dNam(1)$|:,$dNam(2)$|:,$dNam(0)$|:)    ; reorder  (lat,lon,time)
;;work = dtrend( work , False )   ; remove mean + overall long term temporal trend
;;work = dtrend_msg(ispan(0,ntim-1,1)*1.0, work, False, False)  ; remove just  trend
;;if (isatt(work,"_FillValue")) then
;;    delete(work@_FillValue)                  ; dtrend_msg adds this att
;;end if
                                                         ; replace with detrended
;;x    = (/ work($dNam(0)$|:,$dNam(1)$|:,$dNam(2)$|:) /) ; values (time,lat,lon)
;;delete(work)
  x    = dtrend_msg_n(ispan(0,ntim-1,1)*1.0, x, False, False, 0)  ; remove just trend

  if (nDayTot.ge.365) then                     ; rmv dominant signals
      rmvMeans = False                         ; original code did not remove
      x = rmvAnnualCycle(x, spd, nDayTot, rmvMeans, fCrit, 0)
  end if

  if (debug) then
      print("===> Post removal of trend and signal <===")
      printVarSummary( x )                     ; (time,lat,lon)
      printMinMax( x, True )
  end if



  ;-------------------------------------------------------------------
  ; Decompose to Symmetric and Asymmetric parts
  ;-------------------------------------------------------------------
  ; xAS  = decompose2SymAsym( x , 1 )          ; create Asym and Sym parts
  ; if (debug) then
  ;     printVarSummary(xAS)                   ; xAS(lat,lon,time) [iret=1]
  ;     printMinMax(xAS, True)
  ; end if


  ; xnow    = xAS({lat|0:latS}, lon|:, time|:)

  if mode.eq."SH" then
    xAS = x
    xnow = xAS({lat|0:latS}, lon|:, time|:)
    delete(xAS)

    xAS   = xnow
    delete(xnow)
    ndims = dimsizes(xAS)
    nlat = ndims(0)

  else if mode.eq."NH" then
    xAS = x
    xnow = xAS({lat|latN:0}, lon|:, time|:)
    delete(xAS)

    xAS   = xnow
    delete(xnow)
    ndims = dimsizes(xAS)
    nlat = ndims(0)
  end if
  end if

  q      = new((/2,mlon,nSampWin/)       ,typeof(xAS), "No_FillValue")
  peeAS  = new((/nlat,mlon+1,nSampWin+1/),typeof(xAS), "No_FillValue")
  peeAS  = 0.0                       ; initialize

;-------------------------------------------------------------------
; Operate on the xAS array
; NCL does not have a complex 2D FFT at this time.
; Perform a "poorman's" complex 2D FFT by looping over time and space.
; Loop over all latitudes and then perform summing/averaging
; on the spectral results.
;-------------------------------------------------------------------


  do nl=0,nlat-1

     ntStrt = 0
     ntLast = nSampWin-1
     if (debug) then
         print("==============> nl="+nl+" <==============")
     end if

    do nw=0,nWindow-1                        ; temporal window
       if (debug .and. nl.eq.0) then         ; debug
           print("nw="+nw+"  ntStrt="+ntStrt+"   ntLast="+ntLast)
       end if
       work = dtrend( xAS(:,:,ntStrt:ntLast), False )    ; detrend temporal window
       work = taper ( work, tim_taper, 0)    ; taper in "time" [periodic]
                                             ; work(nlat,mlon,N)
      do nt=0,nSampWin-1                     ; for each time perform complex fft in longitude
         q(:,:,nt) = cfftf( work(nl,:,nt), 0.0, 0)   ; space
      end do
       q = q/mlon                                    ; normalize by # lon samples

      do ml=0,mlon-1                         ; for each lon perform complex fft in time
         q(:,ml,:) = cfftf( q(0,ml,:), q(1,ml,:), 0) ; time
      end do
       q = q/nSampWin                                ; normalize by # time samples

;-------------------------------------------------------------------
; At this point 'q(2,mlon,nSampWin)' contains the
; real and imaginary space-time spectrum for this latitude
; ---
; Use Hayashi method to resolve into Progressive [Eastward]
;     and Retrogressive [Westward] Waves.
;-------------------------------------------------------------------

       pee = resolveWavesHayashi( q, nDayWin, spd )
       peeAS(nl,:,:) = peeAS(nl,:,:) + (pee/nWindow)   ; sum window contribution

       ntStrt = ntLast+nSampSkip+1           ; set index for next temporal window
       ntLast = ntStrt+nSampWin-1
    end do                                   ; nw (windows)

  end do                                     ; nl (lat)

  delete(work)
;-------------------------------------------------------------------
; Add meta data to the Hayashi space-time symmetric and asymmetric power
;-------------------------------------------------------------------

  peeAS!0     = "lat"
  peeAS!1     = "wave"
  peeAS!2     = "freq"
  peeAS&lat   = xAS&$dNam(1)$     ; nominally, xAS&lat
  peeAS&wave  = pee&wave
  peeAS&freq  = pee&freq
  peeAS@long_name =  "all components"

  if (debug) then
      printVarSummary( peeAS )
      printMinMax( peeAS , True)
  end if

  delete(  q  )                              ; no longer needed
  delete( pee )

;-------------------------------------------------------------------
; now that we have the power array for sym and asym: use to
;    1) plot raw power spectrum (some smoothing)
;    2) derive and plot the background spectrum (lots of smoothing)
;    3) derive a denoised spectrum that is raw power/background power
;-------------------------------------------------------------------
; psumanti and psumsym will contain the symmetric and asymmetric power
; summed over latitude
;-------------------------------------------------------------------

  psum   = new((/dimsizes(peeAS&wave),dimsizes(peeAS&freq)/),typeof(peeAS), 1e20 )
  psum!0 = "wave"
  psum!1 = "freq"
  psum&wave = peeAS&wave
  psum&freq = peeAS&freq

  psum@long_name = "without constraints summed over lat"
  ; psumsym@long_name  = "Symmetric  summed over lat"

  if (nlat%2.eq.0)    ; use named dimensions to reorder
      psum = dim_sum_n_Wrap(peeAS(:,:,:), 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2-1)   )
  else
      psum = dim_sum_n_Wrap(peeAS(:,:,:), 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2+1:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2)     )
  end if


;-------------------------------------------------------------------
; set the mean to missing to match original code
;-------------------------------------------------------------------
;
  psum(:,{0.0}) = (/ psum@_FillValue /)

  if (debug) then
      printVarSummary( psum )         ; (wave,freq)
      printMinMax( psum , True)
  end if

;-------------------------------------------------------------------
; Apply smoothing to the spectrum. smooth over limited wave numbers
; Smoothing in frequency only (check if mean should be smoothed not smoothing now)
;--
; Smoothing parameters set these larger than the plotting
; wavenumbers to avoid smoothing artifacts
;-------------------------------------------------------------------
  minwav4smth = -27
  maxwav4smth =  27

  do wv=minwav4smth,maxwav4smth
     wk_smooth121( psum({wv},N/2+1:N-1) )
  end do

;-------------------------------------------------------------------
; Log10 scaling
;-------------------------------------------------------------------

  psum_nolog = psum

  psum       = log10(psum)

  psum@long_name = "log10"

  if (debug) then
      printVarSummary( psum )             ; (wave,freq)
      printMinMax( psum , True)
  end if

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        PLOT CODE follows
; --- set some 'plot variables
; set frequency maximum for plotting
; min(user specified frq, maxfrq in window)
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; The following allow DJS variable naming
; to be used with original plot code.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ;caseName       = case
  wkdir          = diro

  frqfftwin      = peeAS&freq
  frqfftwin&freq = peeAS&freq

  wavep1         = peeAS&wave
  wavep1&wave    = peeAS&wave

  if (debug) then
      printVarSummary( frqfftwin )
      printMinMax( frqfftwin, True )
      printVarSummary( wavep1 )
      printMinMax( wavep1, True )
  end if

;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
;-------------------------------------------------------------------
  minfrq4plt =  0.
  maxfrq4plt =  0.8
  minwav4plt = -15
  maxwav4plt =  15

  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))

  fillVal    = 1e20           ; miscellaneous

;=============================================================
;             Start Common Graphics Resources
;=============================================================
   tmFontHgtF            = 0.015     ; not sure why
   tiFontHgtF            = 0.018
   lbFontHgtF            = 0.015
   txFontHgtF            = 0.013

   res = True
   res@gsnFrame          = False
   res@gsnMaximize       = True
   res@gsnPaperOrientation = "portrait"

   res@gsnLeftString     = "Westward"
   res@gsnRightString    = "Eastward"

  ;res@lbBoxMinorExtentF = 0.18
   res@lbLabelFontHeightF= lbFontHgtF
   res@lbOrientation     = "vertical"

   res@cnFillOn          = True
   if (opt .and. isatt(opt, "cnLinesOn") \
           .and. .not.opt@cnLinesOn) then
       res@cnLinesOn     = False
   else
       res@cnLineThicknessF  = 0.5
   end if

   res@tmYLMode          = "Explicit"
   res@tmYLValues        = fspan(minfrq,maxfrq,9)
   res@tmYLLabels        = fspan(minfrq,maxfrq,9)
   res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)

   res@tmYLLabelFontHeightF = tmFontHgtF
   res@tmXBLabelFontHeightF = tmFontHgtF

   res@tiXAxisString     = "Zonal Wave Number"
   res@tiXAxisFontHeightF= tiFontHgtF

   res@tiYAxisString     = "Frequency (cpd)"
   res@tiYAxisFontHeightF= res@tiXAxisFontHeightF

   if (.not.(pltTitle.eq."" .or. pltTitle.eq." ")) then
       res@tiMainString      = pltTitle
       res@tiMainFontHeightF = tiFontHgtF
   end if
   res@txFontHeightF     = tiFontHgtF

;------------------------------------------------------
; Create a list of variable names that have predefined
; contour intervals.
;------------------------------------------------------
   varCnLevels = (/"FLUT" ,"OLR", "olr","U200","U850"   \
                  ,"PRECT","OMEGA500" /)

   if (any(varCnLevels.eq.varName)) then
       res@cnLevelSelectionMode = "ExplicitLevels"
   else
       res@cnLevelSelectionMode = "AutomaticLevels"
   end if

;-------------------------------
; horizontal dashed lines and text for frequency axis [plot only]
;-------------------------------

   freq       = frqfftwin({freq|minfrq:maxfrq})
   wavenumber = wavep1({wave|minwav4plt:maxwav4plt})
   NWVN       = dimsizes(wavenumber)         ; number of wavenumbers

   x1   = wavenumber
   x1!0 = "wave"
   y1   = new (NWVN,float)
   y1!0 = "freq"
   y1   = 1./3.           ; 3 days
   y2   = y1
   y2   = 1./6.           ; 6 days
   y3   = y1
   y3   = 1./30.          ; 30 days
   x2   = new(25,float)
   x2   = 0.0
   y4   = fspan(0.0,1.0,25)

;---------------------------------------------------------------
; dispersion: curves
;---------------------------------------------------------------

   rlat           = 0.0
   Ahe            = (/50.,25.,12./)
   nWaveType      = 6
   nPlanetaryWave = 50
   nEquivDepth    = dimsizes(Ahe)
   Apzwn          = new((/nWaveType,nEquivDepth,nPlanetaryWave/),"double",fillVal)
   Afreq          = Apzwn
   genDispersionCurves(nWaveType, nEquivDepth, nPlanetaryWave, rlat, Ahe, Afreq, Apzwn )

;---------------------------------------------------------------
; dispersion curve and text plot resources
;---------------------------------------------------------------
   dcres = True
   dcres@gsLineThicknessF  = 2.0
   dcres@gsLineDashPattern = 0

   txres = True
   txres@txJust        = "CenterLeft"
   txres@txPerimOn     = True
   txres@txFontHeightF = 0.013
   txres@txBackgroundFillColor = "Background"

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; plotting params for fig 1; subset for plot
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              ; reorder so freq is "y"
   asym       = psum({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   asym!0     = "freq"
   asym&freq  =  freq
   asym!1     = "wave"
   asym&wave  =  wavenumber
   asym@long_name = "Fig_1: log10"

   if (debug) then
       printVarSummary(asym)
       printMinMax(asym, True)
   end if

;------------------------------------------------------
; Fig 1: Pre-defined contour levels [15]  for selected variables [non-log10]
;------------------------------------------------------
   nCn = 15

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2     \ ; unequal
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \ ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9, -4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_1")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_1                    ; user specified limits
   end if

;------------------------------------------------------
; Fig 1: ANTI-SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1a: ASYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.All."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)    ; not required

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       if (isatt(opt,"dirNetCDF")) then
           dirNetCDF = opt@dirNetCDF
       else
           dirNetCDF = "./"
       end if
       if (isatt(opt,"filNetCDF")) then
           filNetCDF = opt@filNetCDF
       else
           filNetCDF = "SpaceTime."+varName+".nc"
       end if
       fNam      = dirNetCDF+filNetCDF
       system ("/bin/rm -f "+fNam)

       ncdf      = addfile(fNam, "c")

       ncdf->FIG_1_ALL  =  asym
   end if

;-----------------------------------------------------------------------------
; ******  now derive and plot the background spectrum (red noise) ************
; [1] Sum power over all latitude
; [2] Put fill value in mean
; [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
;-----------------------------------------------------------------------------
  ;print("======> BACKGROUND <=====")

   psumb = dim_sum_n_Wrap(peeAS, 0)                  ; sum over all latitudes
 ;;psumb = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|:))  ; sum over all latitudes
   psumb@long_name = "Background Spectrum"

   psumb@_FillValue       = fillVal
   psumb(wave|:,freq|N/2) = fillVal
   psumb@_FillValue       = fillVal

   if (debug) then
       printVarSummary(psumb)             ; (wave,freq)
       printMinMax(psumb, True)
   end if

   do tt = N/2+1,N
      smthlen = maxwav4smth-minwav4smth+1
      if (frqfftwin(tt).lt.0.1) then
        do i = 1,5
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.1.and.frqfftwin(tt).lt.0.2) then
        do i = 1,10
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.2.and.frqfftwin(tt).lt.0.3) then
        do i = 1,20
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if

      if (frqfftwin(tt).ge.0.3) then
        do i = 1,40
          wk_smooth121(psumb(freq|tt,{wave|minwav4smth:maxwav4smth}))
        end do
      end if
   end do

   do nw = minwav4smth,maxwav4smth ; smth frequency up to .8 cycles per day
      pt8cpd  = min((/closest_val(.8,frqfftwin),dimsizes(frqfftwin)-1/))
      smthlen = pt8cpd-(N/2+1)+1
     do i = 1,10
        wk_smooth121( psumb({nw},N/2+1:pt8cpd) )
     end do
   end do

;----------------------------------------------------------------
; [1] Put fill value in mean (again)
; [2] SAVE the background spectrum for plotting Fig 3
; [3] LOGARITHMIC SCALING for plotting the background spectrum
;----------------------------------------------------------------

   psumb(wave|:,freq|N/2) = fillVal
   psumb_nolog = xback
   psumb       = log10(xback)            ; LOG10
   psumb@long_name = "log10( background spec )"
   if (debug) then
       print(" ===> Post multiple smoothing by wk_smooth121 <===")
       printVarSummary(psumb_nolog)
       printMinMax(psumb_nolog, True)

       printVarSummary(psumb)
       printMinMax(psumb, True)
   end if

;----------------------------------------------------------------
; set up for plotting  [ subset of frequencies: not necessary]
;----------------------------------------------------------------

   freq      = frqfftwin({freq|minfrq:maxfrq})
   spec      = psumb({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})

   if (debug) then
       printVarSummary(spec)                 ; (freq,wave)
       printMinMax(spec, True)
   end if

;----------------------------------------------------------------
; Fig 2: Predefined explicit contour levels BACKGROUND spectrum  [LOG10]
;----------------------------------------------------------------

   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)       ; allow size to change
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2 \     ; unequal 15
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \   ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9,-4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_2")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_2                   ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig2."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)

   res@gsnCenterString = "Background Power"
   plot = gsn_csm_contour(wks,spec,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)

;*********************************************************************************
; Fig 3a, 3b:  psum_nolog/psumb_nolog  [ratio]
;********************************************************************************

   psumsym_nolog  = psum_nolog / psumb_nolog   ; (wave,freq)

   if (debug) then
       printVarSummary(psum_nolog)           ; (freq,wave)
       printMinMax(psum_nolog, True)
       copy_VarMeta(psum_nolog,psumsym_nolog)
   end if


   asym = psumsym_nolog({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   ;------------------------------------------------------------------
   ; ALL
   ;------------------------------------------------------------------

   symasym = asym

   copy_VarMeta(asym,symasym)
   printVarSummary(asym)
   printVarSummary(symasym)
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3c")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3c                  ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig3.symAsym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Sym+Asym/Background"
   plot = gsn_csm_contour(wks,symasym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

   ; draw dispersion line
   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_polyline(wks,plot,Apzwn(0,0,:),Afreq(0,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,1,:),Afreq(0,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,2,:),Afreq(0,2,:),dcres)

   ; gsn_polyline(wks,plot,)
   ; gsn_polyline(wks,plot,Apzwn(1,0,:),Afreq(1,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,1,:),Afreq(1,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,2,:),Afreq(1,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,0,:),Afreq(2,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,1,:),Afreq(2,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,2,:),Afreq(2,2,:),dcres)

   ; dispersion labels

   gsn_text(wks,plot,"MRG",-10.0,.15,txres)
   ; gsn_text(wks,plot,"n=2 IG",-3.0,.58,txres)
   ; gsn_text(wks,plot,"n=0 EIG",9.5,.50,txres)
   ; gsn_text(wks,plot,"h=50",-10.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-10.0,.63,txres)
   ; gsn_text(wks,plot,"h=12",-10.0,.51,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   ; gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   ; gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   ; gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",4.0,.03,txres)

   frame(wks)
   delete(wks)             ; not required

   ; delete(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       ncdf->FIG_3_BACK = spec
       ncdf->FIG_3_ALL  = asym
   end if

   if (debug) then
       print("********** FINISHED: Space-Time *****************")
   end if
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;          CHANGE TO 2DFFT
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef ("wkSpaceTime_hann")
procedure wkSpaceTime_hann (x[*][*][*]:numeric                  \
                      ,diro[1]:string                      \
                      ,caseName[1]:string                  \
                      ,varName[1]:string                   \
                      ,latBound[1]:numeric                 \
                      ,spd[1]:integer                      \
                      ,nDayWin[1]:integer                  \
                      ,nDaySkip[1]:integer                 \
                      ,lon_taper[1]:float                  \
                      ,opt[1]:logical)

local latN, latS, lonL, lonR, spd, fCrit, tim_taper        \
    , lon_taper, pltType, debug, minwav4smth, maxfrq4plt   \
    , minfrq4plt, maxfrq4plt, minwav4plt, maxwav4plt       \
    , fillVal, nMsg, dimx, ntim, nlat, mlon, nDayTot       \
    , nSampTot, nSampWin, nSampSkip, nWindow, N, dNam, work\
    , rmvMeans, xAS, q, peeAS, nl, ntStrt, ntLast, nw, nt  \
    , ml, psumanti, psumsym, wv, wkdir, caseName, pltFilTit\
    , frqfftwin, wavep1, minfrq, maxfrq, tmFontHgtF, pltTit\
    , tiFontHgtF, lbFontHgtF, txFontHgtF, res, freq        \
    , wavenumber, NWVN, dcres, txres, rlat, Ahe            \
    , nWaveType, nPlanetaryWave, nEquivDepth, Apzwn, Afreq \
    , asym, sym, x1, x2, y1, y2, y3, y4, wks, plot         \
    , psumb, psumsym_nolog, psumanti_nolog,  tt, smthlen   \
    , i, pt8cpd, spec, nCn, nExtra

  begin

  ; if seas.eq."NHSummer" then
  ; else if seas.eq."NHWinter" then
  ; else if seas.eq."NAutumn" then
  ; else if seas.eq. "NHSpring" then
  ; end if
  ; end if
  ; end if
  ; end if
  ;
  if (typeof(x).eq."float" .or. typeof(x).eq."double") then
      print(" ")
  else
      print("wkSpaceTime: input variable must be type float or double")
      print("             input variable type is "+typeof(x))
      exit
  end if

  debug   = False    ; default
  if (opt .and. isatt(opt, "debug") ) then
      debug = opt@debug
  end if

  if (isatt(x,"_FillValue")) then   ; Check for _FillValue .... not allowed
      nMsg = num(ismissing(x))
      if (nMsg.gt.0) then
          print("nMsg="+nMsg+"  User must preprocess to remove _FillValue")
          print("               FFTs do not allow missing values!!       ")
          exit
      end if
      delete(x@_FillValue)     ; avoid warning messages from fft
  end if

  if (debug) then
      printVarSummary( x )
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; x sizes and dimension names
;-------------------------------------------------------------------

  dimx = dimsizes(x)
  ntim = dimx(0)               ; total number of temporal samples
  nlat = dimx(1)
  mlon = dimx(2)

  dNam = getvardims( x )

;-------------------------------------------------------------------
; check to make sure that "x" has full days of data
;-------------------------------------------------------------------

  if ((ntim%spd).ne.0) then
      nExtra = ntim%spd
      print("nExtra="+nExtra+" input array must have complete days only ")
      exit
  end if

;-------------------------------------------------------------------
; Make input arguments into "internal" variables
;-------------------------------------------------------------------
  latN    = latBound
  latS    =-latBound ; make symmetric about the equator

  lonL    = 0        ; -180
  lonR    = 360      ;  180

  fCrit   = 1./nDayWin  ; remove all contributions 'longer'

  tim_taper  = 0.1   ; time taper      [0.1   => 10%]
  lon_taper  = lon_taper   ; longitude taper [0.0 for globe; only global supported]

;-------------------------------------------------------------------
; Check for not allowed actions
;-------------------------------------------------------------------

  if (lon_taper.gt.0.0 .or. (lonR-lonL).ne.360.) then
      print("Code does currently allow lon_taper>0 or (lonR-lonL)<360")
      exit
  end if

;-------------------------------------------------------------------
; OPTIONS
;-------------------------------------------------------------------
  pltType = "ps"     ; default
  if (opt .and. isatt(opt, "pltType") ) then
      if (any(opt@pltType.eq.(/"ps", "eps", "x11", "ncgm", "png"/))) then
          pltType = opt@pltType
      end if
  end if

  pltTit   = caseName+"_"+varName
  pltTitle = pltTit+" LOG[Power: "+latBound+"S-"+latBound+"N]"
  if (opt .and. isatt(opt, "pltTitle") ) then
      pltTitle = opt@pltTitle
  end if

  pltFilTit = pltTit
  replaceChars( pltFilTit, " ", "_") ; spaces not allowed unix file names

  pltColorMap = "amwg_blueyellowred"
  if (opt .and. isatt(opt, "pltColorMap") ) then
      pltColorMap = opt@pltColorMap
  end if


;-------------------------------------------------------------------
; Create required temporal sampling variables
;-------------------------------------------------------------------

  nDayTot   = ntim/spd         ; # of days (total) for input variable
  nSampTot  = nDayTot*spd      ; # of samples (total)
  nSampWin  = nDayWin*spd      ; # of samples per temporal window
  nSampSkip = nDaySkip*spd     ; # of samples to skip between window segments
                               ;   neg means overlap
  nWindow   = (nSampTot-nSampWin)/(nSampWin+nSampSkip)  + 1
  N         =  nSampWin        ; convenience [historical]

  if (nDayTot.lt.nDayWin) then
      print("nDayTot="+nDayTot+" is less the nDayWin="+nDayWin)
      print("        This is not allowed !!       ")
      exit
  end if

;-------------------------------------------------------------------
; Remove dominant signals
; (a) Explicitly remove *long term* linear trend
;     For consistency with JET code keep the grid point means.
;     This necessitates that 'dtrend_msg' be used because 'dtrend'
;     always removes the mean(s).
; (b) All variations >= approx 'nDayWin' days if full year available
;-------------------------------------------------------------------

  dNam = getvardims( x )
;;work = x($dNam(1)$|:,$dNam(2)$|:,$dNam(0)$|:)    ; reorder  (lat,lon,time)
;;work = dtrend( work , False )   ; remove mean + overall long term temporal trend
;;work = dtrend_msg(ispan(0,ntim-1,1)*1.0, work, False, False)  ; remove just  trend
;;if (isatt(work,"_FillValue")) then
;;    delete(work@_FillValue)                  ; dtrend_msg adds this att
;;end if
                                                         ; replace with detrended
;;x    = (/ work($dNam(0)$|:,$dNam(1)$|:,$dNam(2)$|:) /) ; values (time,lat,lon)
;;delete(work)
  x    = dtrend_msg_n(ispan(0,ntim-1,1)*1.0, x, False, False, 0)  ; remove just trend

  if (nDayTot.ge.365) then                     ; rmv dominant signals
      rmvMeans = False                         ; original code did not remove
      x = rmvAnnualCycle(x, spd, nDayTot, rmvMeans, fCrit, 0)
  end if

  if (debug) then
      print("===> Post removal of trend and signal <===")
      printVarSummary( x )                     ; (time,lat,lon)
      printMinMax( x, True )
  end if

;-------------------------------------------------------------------
; Decompose to Symmetric and Asymmetric parts
;-------------------------------------------------------------------
  xAS  = decompose2SymAsym( x , 1 )          ; create Asym and Sym parts
  if (debug) then
      printVarSummary(xAS)                   ; xAS(lat,lon,time) [iret=1]
      printMinMax(xAS, True)
  end if

;-------------------------------------------------------------------
; Because there is the possibility of overlapping *temporal* segments,
; we must use a less efficient approach and detrend/taper
; each window segment as it arises.
;          t0   t1   t2   t3   t4  .................. t(N)
; lon(0):  x00  x01  x02  x03  x04 .................. x0(N)
;     :    :   :   :   :   :                     :
; lon(M):  xM0  xM1  xM2  xM3  xM4 .................. xM(N)
;-------------------------------------------------------------------
; q     - temporary array to hold the 2D complex results
;         for each longitude/time (lon,time) window that is fft'd.
;         This is one instance [realization] of space-time decomposition.
;
; peeAS - symmetric and asymmetric power values in each latitude hemisphere.
;         Add extra lon/time to match JET
;-------------------------------------------------------------------
  q      = new((/2,mlon,nSampWin/)       ,typeof(xAS), "No_FillValue")
  peeAS  = new((/nlat,mlon+1,nSampWin+1/),typeof(xAS), "No_FillValue")
  peeAS  = 0.0                       ; initialize

;-------------------------------------------------------------------
; Operate on the xAS array
; NCL does not have a complex 2D FFT at this time.
; Perform a "poorman's" complex 2D FFT by looping over time and space.
; Loop over all latitudes and then perform summing/averaging
; on the spectral results.
;-------------------------------------------------------------------

  do nl=0,nlat-1

     ntStrt = 0
     ntLast = nSampWin-1
     if (debug) then
         print("==============> nl="+nl+" <==============")
     end if

    do nw=0,nWindow-1                        ; temporal window
       if (debug .and. nl.eq.0) then         ; debug
           print("nw="+nw+"  ntStrt="+ntStrt+"   ntLast="+ntLast)
       end if
       work = dtrend( xAS(:,:,ntStrt:ntLast), False )    ; detrend temporal window
       work = taper ( work, tim_taper, 0)    ; taper in "time" [periodic]
                                             ; work(nlat,mlon,N)
      do nt=0,nSampWin-1                     ; for each time perform complex fft in longitude
         q(:,:,nt) = cfftf( work(nl,:,nt), 0.0, 0)   ; space
      end do
       q = q/mlon                                    ; normalize by # lon samples

      do ml=0,mlon-1                         ; for each lon perform complex fft in time
         q(:,ml,:) = cfftf( q(0,ml,:), q(1,ml,:), 0) ; time
      end do
       q = q/nSampWin                                ; normalize by # time samples

;-------------------------------------------------------------------
; At this point 'q(2,mlon,nSampWin)' contains the
; real and imaginary space-time spectrum for this latitude
; ---
; Use Hayashi method to resolve into Progressive [Eastward]
;     and Retrogressive [Westward] Waves.
;-------------------------------------------------------------------

       pee = resolveWavesHayashi( q, nDayWin, spd )
       peeAS(nl,:,:) = peeAS(nl,:,:) + (pee/nWindow)   ; sum window contribution

       ntStrt = ntLast+nSampSkip+1           ; set index for next temporal window
       ntLast = ntStrt+nSampWin-1
    end do                                   ; nw (windows)

  end do                                     ; nl (lat)

  delete(work)
;-------------------------------------------------------------------
; Add meta data to the Hayashi space-time symmetric and asymmetric power
;-------------------------------------------------------------------

  peeAS!0     = "lat"
  peeAS!1     = "wave"
  peeAS!2     = "freq"
  peeAS&lat   = xAS&$dNam(1)$     ; nominally, xAS&lat
  peeAS&wave  = pee&wave
  peeAS&freq  = pee&freq
  peeAS@long_name =  "symmetric and asymmetric components"

  if (debug) then
      printVarSummary( peeAS )
      printMinMax( peeAS , True)
  end if

  delete(  q  )                              ; no longer needed
  delete( pee )

;-------------------------------------------------------------------
; now that we have the power array for sym and asym: use to
;    1) plot raw power spectrum (some smoothing)
;    2) derive and plot the background spectrum (lots of smoothing)
;    3) derive a denoised spectrum that is raw power/background power
;-------------------------------------------------------------------
; psumanti and psumsym will contain the symmetric and asymmetric power
; summed over latitude
;-------------------------------------------------------------------

  psumanti   = new((/dimsizes(peeAS&wave),dimsizes(peeAS&freq)/),typeof(peeAS), 1e20 )
  psumanti!0 = "wave"
  psumanti!1 = "freq"
  psumanti&wave = peeAS&wave
  psumanti&freq = peeAS&freq

  psumsym       = psumanti

  psumanti@long_name = "Asymmetric summed over lat"
  psumsym@long_name  = "Symmetric  summed over lat"

  if (nlat%2.eq.0)    ; use named dimensions to reorder
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2-1,:,:)   , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2-1)   )
  else
      psumanti = dim_sum_n_Wrap(peeAS(nlat/2+1:nlat-1,:,:), 0)
      psumsym  = dim_sum_n_Wrap(peeAS(0:nlat/2,:,:)       , 0)
    ;;psumanti = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|nlat/2+1:nlat-1))
    ;;psumsym  = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|0:nlat/2)     )
  end if

;-------------------------------------------------------------------
; since summing over half the array (symmetric,asymmetric) the
; total variance is 2x the half sum
;-------------------------------------------------------------------
  psumanti = 2.0*psumanti
  psumsym  = 2.0*psumsym

;-------------------------------------------------------------------
; set the mean to missing to match original code
;-------------------------------------------------------------------
;
  psumanti(:,{0.0}) = (/ psumanti@_FillValue /)
  psumsym (:,{0.0}) = (/ psumsym@_FillValue /)

  if (debug) then
      printVarSummary( psumanti )         ; (wave,freq)
      printMinMax( psumsym , True)
  end if

;-------------------------------------------------------------------
; Apply smoothing to the spectrum. smooth over limited wave numbers
; Smoothing in frequency only (check if mean should be smoothed not smoothing now)
;--
; Smoothing parameters set these larger than the plotting
; wavenumbers to avoid smoothing artifacts
;-------------------------------------------------------------------
  minwav4smth = -27
  maxwav4smth =  27

  do wv=minwav4smth,maxwav4smth
     wk_smooth121( psumanti({wv},N/2+1:N-1) )
     wk_smooth121( psumsym ({wv},N/2+1:N-1) )
  end do

;-------------------------------------------------------------------
; Log10 scaling
;-------------------------------------------------------------------

  psumanti_nolog = psumanti
  psumsym_nolog  = psumsym

  psumanti       = log10(psumanti)
  psumsym        = log10(psumsym)

  psumanti@long_name = "log10(Asymmetric)"
  psumsym@long_name  = "log10(Symmetric)"

  if (debug) then
      printVarSummary( psumanti )             ; (wave,freq)
      printMinMax( psumanti , True)
      printVarSummary( psumsym )
      printMinMax( psumsym , True)
  end if

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        PLOT CODE follows
; --- set some 'plot variables
; set frequency maximum for plotting
; min(user specified frq, maxfrq in window)
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; The following allow DJS variable naming
; to be used with original plot code.
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 ;caseName       = case
  wkdir          = diro

  frqfftwin      = peeAS&freq
  frqfftwin&freq = peeAS&freq

  wavep1         = peeAS&wave
  wavep1&wave    = peeAS&wave

  if (debug) then
      printVarSummary( frqfftwin )
      printMinMax( frqfftwin, True )
      printVarSummary( wavep1 )
      printMinMax( wavep1, True )
  end if

;-------------------------------------------------------------------
; plotting parameters freq and wavenumbers to plot
;-------------------------------------------------------------------
  minfrq4plt =  0.
  maxfrq4plt =  0.8
  minwav4plt = -15
  maxwav4plt =  15

  minfrq     = minfrq4plt
  maxfrq     = min((/maxfrq4plt,max(frqfftwin)/))

  fillVal    = 1e20           ; miscellaneous

;=============================================================
;             Start Common Graphics Resources
;=============================================================
   tmFontHgtF            = 0.015     ; not sure why
   tiFontHgtF            = 0.018
   lbFontHgtF            = 0.015
   txFontHgtF            = 0.013

   res = True
   res@gsnFrame          = False
   res@gsnMaximize       = True
   res@gsnPaperOrientation = "portrait"

   res@gsnLeftString     = "Westward"
   res@gsnRightString    = "Eastward"

  ;res@lbBoxMinorExtentF = 0.18
   res@lbLabelFontHeightF= lbFontHgtF
   res@lbOrientation     = "vertical"

   res@cnFillOn          = True
   if (opt .and. isatt(opt, "cnLinesOn") \
           .and. .not.opt@cnLinesOn) then
       res@cnLinesOn     = False
   else
       res@cnLineThicknessF  = 0.5
   end if

   res@tmYLMode          = "Explicit"
   res@tmYLValues        = fspan(minfrq,maxfrq,9)
   res@tmYLLabels        = fspan(minfrq,maxfrq,9)
   res@tmYLMinorValues   = fspan(minfrq,maxfrq,17)

   res@tmYLLabelFontHeightF = tmFontHgtF
   res@tmXBLabelFontHeightF = tmFontHgtF

   res@tiXAxisString     = "Zonal Wave Number"
   res@tiXAxisFontHeightF= tiFontHgtF

   res@tiYAxisString     = "Frequency (cpd)"
   res@tiYAxisFontHeightF= res@tiXAxisFontHeightF

   if (.not.(pltTitle.eq."" .or. pltTitle.eq." ")) then
       res@tiMainString      = pltTitle
       res@tiMainFontHeightF = tiFontHgtF
   end if
   res@txFontHeightF     = tiFontHgtF

;------------------------------------------------------
; Create a list of variable names that have predefined
; contour intervals.
;------------------------------------------------------
   varCnLevels = (/"FLUT" ,"OLR", "olr","U200","U850"   \
                  ,"PRECT","OMEGA500" /)

   if (any(varCnLevels.eq.varName)) then
       res@cnLevelSelectionMode = "ExplicitLevels"
   else
       res@cnLevelSelectionMode = "AutomaticLevels"
   end if

;-------------------------------
; horizontal dashed lines and text for frequency axis [plot only]
;-------------------------------

   freq       = frqfftwin({freq|minfrq:maxfrq})
   wavenumber = wavep1({wave|minwav4plt:maxwav4plt})
   NWVN       = dimsizes(wavenumber)         ; number of wavenumbers

   x1   = wavenumber
   x1!0 = "wave"
   y1   = new (NWVN,float)
   y1!0 = "freq"
   y1   = 1./3.           ; 3 days
   y2   = y1
   y2   = 1./6.           ; 6 days
   y3   = y1
   y3   = 1./30.          ; 30 days
   x2   = new(25,float)
   x2   = 0.0
   y4   = fspan(0.0,1.0,25)

;---------------------------------------------------------------
; dispersion: curves
;---------------------------------------------------------------

   rlat           = 0.0
   Ahe            = (/50.,25.,12./)
   nWaveType      = 6
   nPlanetaryWave = 50
   nEquivDepth    = dimsizes(Ahe)
   Apzwn          = new((/nWaveType,nEquivDepth,nPlanetaryWave/),"double",fillVal)
   Afreq          = Apzwn
   genDispersionCurves(nWaveType, nEquivDepth, nPlanetaryWave, rlat, Ahe, Afreq, Apzwn )

;---------------------------------------------------------------
; dispersion curve and text plot resources
;---------------------------------------------------------------
   dcres = True
   dcres@gsLineThicknessF  = 2.0
   dcres@gsLineDashPattern = 0

   txres = True
   txres@txJust        = "CenterLeft"
   txres@txPerimOn     = True
   txres@txFontHeightF = 0.013
   txres@txBackgroundFillColor = "Background"

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; plotting params for fig 1; subset for plot
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                                              ; reorder so freq is "y"
   asym       = psumanti({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   asym!0     = "freq"
   asym&freq  =  freq
   asym!1     = "wave"
   asym&wave  =  wavenumber
   asym@long_name = "Fig_1: log10(Asymmetric)"

   sym        = psumsym({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   sym@long_name  = "Fig_1: log10(Symmetric)"

   if (debug) then
       printVarSummary(asym)
       printMinMax(asym, True)
       printVarSummary(sym)
       printMinMax(sym, True)
   end if

;------------------------------------------------------
; Fig 1: Pre-defined contour levels [15]  for selected variables [non-log10]
;------------------------------------------------------
   nCn = 15

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2     \ ; unequal
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \ ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9, -4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_1")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_1                    ; user specified limits
   end if

;------------------------------------------------------
; Fig 1: ANTI-SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1a: ASYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Asym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)    ; not required

;------------------------------------------------------
; Fig 1: SYMMETRIC
;------------------------------------------------------
  ;print("======> Fig 1b: SYMMETRIC <=====")

   wks  = gsn_open_wks(pltType,wkdir+"/"+"Fig1.Sym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Symmetric"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)
   delete(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       if (isatt(opt,"dirNetCDF")) then
           dirNetCDF = opt@dirNetCDF
       else
           dirNetCDF = "./"
       end if
       if (isatt(opt,"filNetCDF")) then
           filNetCDF = opt@filNetCDF
       else
           filNetCDF = "SpaceTime."+varName+".nc"
       end if
       fNam      = dirNetCDF+filNetCDF
       system ("/bin/rm -f "+fNam)

       ncdf      = addfile(fNam, "c")

       ncdf->FIG_1_SYM  =  sym
       ncdf->FIG_1_ASYM = asym
   end if

;-----------------------------------------------------------------------------
; ******  now derive and plot the background spectrum (red noise) ************
; [1] Sum power over all latitude
; [2] Put fill value in mean
; [3] Apply smoothing to the spectrum. This smoothing DOES include wavenumber zero.
;-----------------------------------------------------------------------------
  ;print("======> BACKGROUND <=====")

   psumb = dim_sum_n_Wrap(peeAS, 0)                  ; sum over all latitudes
 ;;psumb = dim_sum_Wrap(peeAS(wave|:,freq|:,lat|:))  ; sum over all latitudes
   psumb@long_name = "Background Spectrum"

   psumb@_FillValue       = fillVal
   psumb(wave|:,freq|N/2) = fillVal
   psumb@_FillValue       = fillVal

   if (debug) then
       printVarSummary(psumb)             ; (wave,freq)
       printMinMax(psumb, True)
   end if

   do tt = N/2+1,N
      smthlen = maxwav4smth-minwav4smth+1
      if (frqfftwin(tt).lt.0.1) then
        do i = 1,5
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.1.and.frqfftwin(tt).lt.0.2) then
        do i = 1,10
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if
      if (frqfftwin(tt).ge.0.2.and.frqfftwin(tt).lt.0.3) then
        do i = 1,20
          wk_smooth121( psumb(freq|tt,{wave|minwav4smth:maxwav4smth}) )
        end do
      end if

      if (frqfftwin(tt).ge.0.3) then
        do i = 1,40
          wk_smooth121(psumb(freq|tt,{wave|minwav4smth:maxwav4smth}))
        end do
      end if
   end do

   do nw = minwav4smth,maxwav4smth ; smth frequency up to .8 cycles per day
      pt8cpd  = min((/closest_val(.8,frqfftwin),dimsizes(frqfftwin)-1/))
      smthlen = pt8cpd-(N/2+1)+1
     do i = 1,10
        wk_smooth121( psumb({nw},N/2+1:pt8cpd) )
     end do
   end do

;----------------------------------------------------------------
; [1] Put fill value in mean (again)
; [2] SAVE the background spectrum for plotting Fig 3
; [3] LOGARITHMIC SCALING for plotting the background spectrum
;----------------------------------------------------------------

   psumb(wave|:,freq|N/2) = fillVal
   psumb_nolog = psumb
   psumb       = log10(psumb)            ; LOG10
   psumb@long_name = "log10( background spec )"
   if (debug) then
       print(" ===> Post multiple smoothing by wk_smooth121 <===")
       printVarSummary(psumb_nolog)
       printMinMax(psumb_nolog, True)

       printVarSummary(psumb)
       printMinMax(psumb, True)
   end if

;----------------------------------------------------------------
; set up for plotting  [ subset of frequencies: not necessary]
;----------------------------------------------------------------

   freq      = frqfftwin({freq|minfrq:maxfrq})
   spec      = psumb({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})

   if (debug) then
       printVarSummary(spec)                 ; (freq,wave)
       printMinMax(spec, True)
   end if

;----------------------------------------------------------------
; Fig 2: Predefined explicit contour levels BACKGROUND spectrum  [LOG10]
;----------------------------------------------------------------

   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)       ; allow size to change
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR" .or. varName .eq. "olr") then
       res@cnLevels = (/-1.2,-1.1,-1.0,-0.8,-0.6,-0.4,-0.2 \     ; unequal 15
                       , 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.1,1.2/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/-18.2,-18.0,-17.8,-17.6,-17.5,-17.4,-17.3 \   ; unequal
                       ,-17.2,-17.1,-17.0,-16.9,-16.8,-16.7,-16.6,-16.5/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(-3.3, 0.9, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(-3.25, 0.25, nCn)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = fspan(-5.9,-4.5, nCn)
   end if

   if (opt .and. isatt(opt, "Fig_2")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_2                   ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig2."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)

   res@gsnCenterString = "Background Power"
   plot = gsn_csm_contour(wks,spec,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines
   frame(wks)

;*********************************************************************************
; Fig 3a, 3b:  psum_nolog/psumb_nolog  [ratio]
;********************************************************************************

   psumsym_nolog  = psumsym_nolog /psumb_nolog   ; (wave,freq)
   psumanti_nolog = psumanti_nolog/psumb_nolog

   if (debug) then
       printVarSummary(psumanti_nolog)           ; (freq,wave)
       printMinMax(psumanti_nolog, True)
       printVarSummary(psumsym_nolog)            ; (freq,wave)
       printMinMax(psumsym_nolog, True)
   end if

;-----------------------------------------------------------
; ANTI-SYMMETRIC: RATIO ( subset to plot )
;-----------------------------------------------------------

   asym = psumanti_nolog({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})
   sym  = psumsym_nolog({freq|minfrq:maxfrq},{wave|minwav4plt:maxwav4plt})

   if (debug) then
       printVarSummary(asym)                    ; (freq,wave)
       printMinMax(asym, True)
       printVarSummary(sym)                     ; (freq,wave)
       printMinMax(sym, True)
   end if

   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)                     ; allow size to change
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = fspan(0.3, 1.7, nCn)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = fspan(0.4, 1.8, nCn)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = fspan(0.4, 1.8, nCn)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/0.6,0.7 ,0.8,0.9 ,1.0,1.1,1.15,1.2,1.25 \
                       ,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/0.6,0.7,0.8,0.9,1.0,1.1,1.15,1.2,1.25 \
                       ,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if

   if (opt .and. isatt(opt, "Fig_3a")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3a                  ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig3.Asym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Anti-Symmetric/Background"
   plot = gsn_csm_contour(wks,asym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

; draw dispersion line

   gsn_polyline(wks,plot,Apzwn(0,0,:),Afreq(0,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,1,:),Afreq(0,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,2,:),Afreq(0,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,0,:),Afreq(1,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,1,:),Afreq(1,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(1,2,:),Afreq(1,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,0,:),Afreq(2,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,1,:),Afreq(2,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(2,2,:),Afreq(2,2,:),dcres)

; dispersion labels

   gsn_text(wks,plot,"MRG",-10.0,.15,txres)
   gsn_text(wks,plot,"n=2 IG",-3.0,.58,txres)
   gsn_text(wks,plot,"n=0 EIG",9.5,.50,txres)
   gsn_text(wks,plot,"h=50",-10.0,.78,txres)
   gsn_text(wks,plot,"h=25",-10.0,.63,txres)
   gsn_text(wks,plot,"h=12",-10.0,.51,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)

   frame(wks)
   delete(wks)             ; not required

;------------------------------------------------------------------
; SYMMETRIC
;------------------------------------------------------------------
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3b")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3b                  ; user specified
   end if

   wks  = gsn_open_wks(pltType, wkdir+"/Fig3.Sym."+pltFilTit)
   gsn_define_colormap(wks ,pltColorMap)
   res@gsnCenterString = "Symmetric/Background"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",5.0,.06,txres)
   frame(wks)

   if (debug .or. (opt .and. isatt(opt,"Fig_3_statInfo") \
                       .and. opt@Fig_3_statInfo)       ) then
       statAsymSym( asym, "Fig_3a" )
       statAsymSym(  sym, "Fig_3b" )
   end if

   delete(wks)

   ;------------------------------------------------------------------
   ; ALL
   ;------------------------------------------------------------------

   symasym = sym+asym

   copy_VarMeta(sym,symasym)
   printVarSummary(sym)
   printVarSummary(symasym)
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3c")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3c                  ; user specified
   end if

   wks  = gsn_open_wks(pltType,wkdir+"/Fig3.symAsym."+pltFilTit)
   gsn_define_colormap(wks,pltColorMap)
   res@gsnCenterString = "Sym+Asym/Background"
   plot = gsn_csm_contour(wks,symasym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

; draw dispersion line
   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_polyline(wks,plot,Apzwn(0,0,:),Afreq(0,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,1,:),Afreq(0,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(0,2,:),Afreq(0,2,:),dcres)

   tdres = dcres
   tdres@gsLineColor = "red"

   ftd1 = 1/2.5
   ftd2 = 1/5.
   gsn_polyline(wks,plot,(/-15,-15/), (/ftd2,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-6,-6/), (/ftd2,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-15,-6/), (/ftd1,ftd1/),tdres)
   gsn_polyline(wks,plot,(/-15,-6/), (/ftd2,ftd2/),tdres)

   ; fmjo1 = 1/20.
   ; fmjo2 = 1/100.
   ; gsn_polyline(wks,plot,(/1,1/), (/fmjo2,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/5,5/), (/fmjo2,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/1,5/), (/fmjo1,fmjo1/),tdres)
   ; gsn_polyline(wks,plot,(/1,5/), (/fmjo2,fmjo2/),tdres)
   ;
   ; fmrg1 = 1/3.
   ; fmrg2 = 1/10.
   ; kmrg1 = -1
   ; kmrg2 = -10
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg1/), (/fmrg2,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg2,kmrg2/), (/fmrg2,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg2/), (/fmrg1,fmrg1/),tdres)
   ; gsn_polyline(wks,plot,(/kmrg1,kmrg2/), (/fmrg2,fmrg2/),tdres)
   ;
   ;
   ; ferw1 = 1/9.
   ; ferw2 = 1/72.
   ; kerw1 = -1
   ; kerw2 = -10
   ; gsn_polyline(wks,plot,(/kerw1,kerw1/), (/ferw2,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw2,kerw2/), (/ferw2,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw1,kerw2/), (/ferw1,ferw1/),tdres)
   ; gsn_polyline(wks,plot,(/kerw1,kerw2/), (/ferw2,ferw2/),tdres)
   ;
   ;
   ; fklv1 = 1/3.
   ; fklv2 = 1/10.
   ; kklv1 = 1
   ; kklv2 = 14
   ; gsn_polyline(wks,plot,(/kklv1,kklv1/), (/fklv2,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv2,kklv2/), (/fklv2,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv1,kklv2/), (/fklv1,fklv1/),tdres)
   ; gsn_polyline(wks,plot,(/kklv1,kklv2/), (/fklv2,fklv2/),tdres)



   fer1 = 1/7.
   fer2 = 1/10.
   res_poly = True
   res_poly@gsFillColor = "black"
   gsn_polygon(wks,plot,(/12,15,15,12,12/),(/fer1,fer1,fer2,fer2,fer1/),res_poly)

   ; gsn_polyline(wks,plot,)
   ; gsn_polyline(wks,plot,Apzwn(1,0,:),Afreq(1,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,1,:),Afreq(1,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(1,2,:),Afreq(1,2,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,0,:),Afreq(2,0,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,1,:),Afreq(2,1,:),dcres)
   ; gsn_polyline(wks,plot,Apzwn(2,2,:),Afreq(2,2,:),dcres)

; dispersion labels

   gsn_text(wks,plot,"MRG",-10.0,.15,txres)
   ; gsn_text(wks,plot,"n=2 IG",-3.0,.58,txres)
   ; gsn_text(wks,plot,"n=0 EIG",9.5,.50,txres)
   ; gsn_text(wks,plot,"h=50",-10.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-10.0,.63,txres)
   ; gsn_text(wks,plot,"h=12",-10.0,.51,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   ; gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   ; gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   ; gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   ; gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",4.0,.03,txres)

   frame(wks)
   delete(wks)             ; not required

;------------------------------------------------------------------
; SYMMETRIC
;------------------------------------------------------------------
   if (isatt(res,"cnLevels")) then
       delete(res@cnLevels)
   end if

   if (varName .eq. "FLUT" .or. varName .eq. "OLR".or. varName .eq. "olr") then
       res@cnLevels = (/.3,.4,.5,.6,.7,.8,.9,1.,1.1,1.2,1.4,1.7,2.,2.4,2.8/)
   end if
   if (varName .eq. "U200") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "U850") then
       res@cnLevels = (/.4,.6,.8,1.,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.,2.2,2.4,2.6/)
   end if
   if (varName .eq. "PRECT") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.15,1.2,1.25,1.3,1.35,1.4,1.45,1.5,1.6/)
   end if
   if (varName .eq. "OMEGA500") then
       res@cnLevels = (/.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2./)
   end if

   if (opt .and. isatt(opt, "Fig_3b")) then
       res@cnLevelSelectionMode = "ExplicitLevels"
       res@cnLevels = opt@Fig_3b                  ; user specified
   end if

   wks  = gsn_open_wks(pltType, wkdir+"/Fig3.Sym."+pltFilTit)
   gsn_define_colormap(wks ,pltColorMap)
   res@gsnCenterString = "Symmetric/Background"
   plot = gsn_csm_contour(wks,sym,res)
   addHorVertLines(wks, plot, x1,x2,y1,y2,y3,y4 ) ; add dashed lines

   gsn_polyline(wks,plot,Apzwn(3,0,:),Afreq(3,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,1,:),Afreq(3,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(3,2,:),Afreq(3,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,0,:),Afreq(4,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,1,:),Afreq(4,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(4,2,:),Afreq(4,2,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,0,:),Afreq(5,0,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,1,:),Afreq(5,1,:),dcres)
   gsn_polyline(wks,plot,Apzwn(5,2,:),Afreq(5,2,:),dcres)

   gsn_text(wks,plot,"Kelvin",11.0,.40,txres)
   gsn_text(wks,plot,"n=1 ER",-10.7,.07,txres)
   gsn_text(wks,plot,"n=1 IG",-3.0,.45,txres)
   gsn_text(wks,plot,"h=50",-14.0,.78,txres)
   gsn_text(wks,plot,"h=25",-14.0,.60,txres)
   gsn_text(wks,plot,"h=12",-14.0,.46,txres)
   gsn_text(wks,plot,"TD-type",-10.0,.3,txres)
   gsn_text(wks,plot,"MJO",5.0,.06,txres)
   frame(wks)

;------------------------------------------------------
; Is netCDF option set?
;------------------------------------------------------
   if (opt .and. isatt(opt,"netCDF") .and. opt@netCDF) then
       ncdf->FIG_3_BACK = spec
       ncdf->FIG_3_SYM  =  sym
       ncdf->FIG_3_ASYM = asym
       ncdf->FIG_3_ALL  = sym+asym
   end if

   if (debug) then
       print("********** FINISHED: Space-Time *****************")
   end if
end
