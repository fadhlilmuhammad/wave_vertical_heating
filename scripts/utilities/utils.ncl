

; UTILITY codes
;

;----------------------------------------------------------------------
; Function that retrieves contour levels from a plot and changes the
; the requested contour levels to the requested thicknesses.
;
; This code is somewhat complicated because we are checking to make
; sure we don't override any line thicknesses that might have
; originally been set by the user.
;
; This code is very similar to the set_contour_line_patterns code.
;----------------------------------------------------------------------
undef("set_contour_line_thicknesses")
procedure set_contour_line_thicknesses(plot,levels_to_change[*]:numeric,\
                                       thicknesses[*]:numeric)
local ii, levels, nlevels, line_thicknesses, nchange, nthick, \
      mono_thickness, line_thickness,changed_a_level
begin
;---Retrieve the original line thicknesses (or thickness) used for the plot.
  getvalues plot@contour
    "cnLevels"            : levels
    "cnLineThicknesses"   : line_thicknesses
    "cnMonoLineThickness" : mono_thickness
    "cnLineThicknessF"    : line_thickness
  end getvalues
;----------------------------------------------------------------------
; If the original contour plot used a single dash line thickness
; all contour lines, then make sure we use that again for
; any contour lines that are not being changed. Otherwise,
; we assume the user set an array of line thicknesses, and
; we'll use those.
;----------------------------------------------------------------------

;---Check user set line thicknesses, and apply new thicknesses if needed
  if(mono_thickness) then
    line_thicknesses = line_thickness
  end if

  nlevels = dimsizes(levels)
  nchange = dimsizes(levels_to_change)
  nthick  = dimsizes(thicknesses)

  if(nthick.ne.nchange) then
      print("set_contour_line_thicknesses: error: the contour line thicknesses must be an")
      print("array of the same length as the number of contour levels to change.")
  end if

;---Array to hold index values of contour levels that need to have a thickness applied
  changed_a_level = False
  do n=0,nchange-1
    ii := ind(levels.eq.levels_to_change(n))
    if(ismissing(ii(0))) then
      print("set_contour_line_thicknesses: warning, no contour level equal to " + levels_to_change(n))
    else
      line_thicknesses(ii(0)) = thicknesses(n)
      changed_a_level     = True
    end if
  end do

  if(changed_a_level) then
    setvalues plot@contour
      "cnMonoLineThickness" : False             ; allows an array of line thicknesses to be set
      "cnLineThicknesses"   : line_thicknesses
    end setvalues
  end if
end



;;;MASKING

undef("custom_mask")
function custom_mask(x,lat, lon,maskreg, opt ) ;var in the maskfile
local maskvar, fin
begin

  if maskreg.eq."ocean" then
    nmask = 0
  else if maskreg.eq."land" then
    nmask = 1
  end if
end if

  if opt.eq.False then
    opt@var = "LSMASK"
    opt@maskfile = ncargpath("data") + "/cdf/landsea.nc"
  end if
  ; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
  fin   = addfile(opt@maskfile, "r")
  lsm_x = landsea_mask(fin->$opt@var$,lat,lon)
  xmask = mask(x,lsm_x.eq.0,False)        ; mask out all ocean points from t85

  copy_VarMeta(x,xmask)
  return(xmask)
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; percentiles.ncl
;; Carl Schreck (carl@cicsnc.org)
;; August 2012
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; percent_to_value
;; percent_to_value_2d
;; value_to_percent
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

undef ( "percent_to_value" )
function percent_to_value( \
  i_data           : numeric, \
  i_percentiles[*] : numeric  \
)
local None
begin

  ; print("check1")
  retVal = new( dimsizes(i_percentiles), float )
  ; print(retVal)
  data1d = ndtooned( i_data )

if(all(ismissing(data1d))) then
    ; print("ALL MISSING data...")
    return(data1d@_FillValue)
else
    ; print("check2")
  notMissing = data1d( ind(.not.ismissing(data1d) ) )
  ; print(notMissing)

  qsort(notMissing)

  do p = 0, dimsizes(i_percentiles)-1
    ; print(data1d)
;   pInd = round( i_percentiles(p) * .01 * dimsizes(notMissing) + 0.5, 3 ) -1
;   pInd = where( pInd.ge.dimsizes(notMissing), dimsizes(notMissing)-1, pInd )
    floatInd = i_percentiles(p) * .01 * dimsizes(notMissing) - 0.5
    floorInd = toint( floor(floatInd) )
    floorInd = where( floorInd.lt.0, 0, floorInd )
    ceilInd = toint( ceil(floatInd) )
    ceilInd = where( ceilInd.ge.dimsizes(notMissing), \
        dimsizes(notMissing)-1, ceilInd )
;   print(pInd + " " + dimsizes(notMissing))
    if( ceilInd.eq.floorInd ) then
      retVal(p) = notMissing(floorInd)
    else
      retVal(p) = notMissing(floorInd) * ( ceilInd - floatInd ) \
          + notMissing(ceilInd) * ( floatInd - floorInd )
    end if
  end do

  return(retVal)
end if
end ; percent_to_value


undef("percentile")
function percentile(x:numeric,P:float)
begin
if(P.gt.0) then
P = P/100.
end if
  nall = x
  igood = ind(.not.ismissing(nall))
  nx = nall(igood)
  ; print(xgood)
  ; printVarSummary(nx)

  ; print(nx)
  index = round(P*dimsizes(nx),3)
  print(index)
  print("Index of percentile: "+index+ " of size: "+dimsizes(nx))
  ib = dim_pqsort(nx,2)

  ; print(ib)
  ; print("The "+P*100+" percentile: "+nx(index))
  return(nx(index-1))
end

undef("hann_lon")
function hann_lon(data,c,N) ;data = data, c=central longitude, N=longitude window length
	begin

	N2 = toint(N/2)
	; L  = dimsizes(x)
	
  if (c.lt.(N/2+10)).or.(c.gt.(360-N/2-10)) then
  print("Central longitude is too small to the edge of the coordinates, reordering data")
    datapivot = lonPivot(data, (360.-180)) ; reorder
  else  
  print("Central longitude is ok, no reordering needed")
    datapivot = data
  end if

	x = datapivot&lon
	x!0 = "lon"
	x&lon = x

  ic = ind(x.eq.c)
  ngrid = (x(1)-x(0))
  print("ngrid")
  print(ngrid)
  gridwindow = toint(N2/ngrid)
  print("gridwindow")
  print(gridwindow)

	; print(lonL)
	printVarSummary(x)
	print(x)
  print(ic)
	; xW = x({c-N2:c+N2})
  xW = x(ic-(gridwindow):ic+(gridwindow))
  print(xW)
	L  = dimsizes(xW)

	pi = get_pi("float")
	print(pi)

  print(x(ic-(gridwindow)))
	it = ind(x(ic-(gridwindow)).eq.x)
	print(it)
	w  = new((/dimsizes(x)/), "float")
	w = 0

	k  = new((/L/), "float")
	d2r = get_d2r("float")
	; counter = 1
	; print(L)
		do t = 0, L-1
			; print(t)
			; rad = xW(counter)/(xW(L-1)-1)
			; rad = xW(counter)/(xW(L-1))
			n_N = tofloat(tofloat(t+1)/L)
			; print(L)
			; print(rad)
			; k(t) = sin(pi*n_N)^2
			k(t) = 0.5 * (1 - cos(2*pi*n_N)) ; Hann function

			; counter = counter+1
		end do

	print(k)
	w(it:it+L-1) = k

	copy_VarCoords(x,w)
	w@long_name = "hanning weights"
	; print(w)
	return w
	
	end


undef("hann_taper_longitude")
function hann_taper_longitude(data,c,N, dimlon)
	begin 
	print("tapering....")
	w = hann_lon(data,c,N)
  print(w)

  if (c.lt.(N/2+10)).or.(c.gt.(360-N/2-10)) then
  print("Central longitude is too small to the edge of the coordinates, reordering data")
    data_pivot = lonPivot(data, (360.-180)) ; reorder
  else  
  print("Central longitude is ok, no reordering needed")
    data_pivot = data
  end if
 
	data_tapered = data_pivot * conform_dims(dimsizes(data_pivot), w, dimlon)

	copy_VarMeta(data_pivot, data_tapered)
	data_tapered@description = "Hanning Tapered Data"
	data_tapered@weights = w 
	data_tapered@central_lon = c



	return data_tapered
	
	end 


undef ( "percent_to_value_2d" )
function percent_to_value_2d( \
  i_data[*][*][*]  : numeric, \
  i_percentile[1]  : numeric  \
)
local None
begin

  dims = dimsizes(i_data)
  sorted = i_data
  ip = dim_pqsort_n( sorted, 2, 0 )

  pInd = round( i_percentile * .01 * dims(0) + 0.5, 3 ) -1
  pInd = where( pInd.ge.dims(0), dims(0)-1, pInd )
  retVal = sorted(pInd,:,:)
  copy_VarAtts( i_data, retVal )
  copy_VarCoords( i_data(0,:,:), retVal )

  return(retVal)
end ; percent_to_value_2d

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Choosing season
;; Fadhlil R. Muhammad (fadhlilrizki@student.unimelb.edu.au)
;; May 2023
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; choose_season
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


undef("choose_season")
function choose_season(seas, ymdhms)
local mm
begin
  mm      = floattoint( ymdhms(:,1) )        ; just months


  if seas.eq."SHWinter" .or. seas.eq."JJA" then
    iseason = ind( mm.ge.6 .and. mm.le.8)
  else if seas.eq."SHSpring" .or. seas.eq."SON" then
    iseason = ind( mm.ge.9 .and. mm.le.11)
  else if seas.eq."SHSummer" .or. seas.eq."DJF" then
    iseason = ind( (mm.eq.12).or.(mm.eq.1).or.(mm.eq.2))
  else if seas.eq."SHAutumn" .or. seas.eq."MAM" then
    iseason = ind( mm.ge.3 .and. mm.le.5)
  else if seas.eq."NDJFMA" then
    iseason = ind(mm.ge.11 .or. mm.le.4)
  else if seas.eq."MJJASO" then
    iseason = ind(mm.ge.5 .and. mm.le.10)
  else if seas.eq."All Dates" then
    iseason = ind(mm)
  else if seas.eq."All_dates" then
    iseason = ind(mm)
  else if seas.eq."all" then 
    iseason = ind(mm)
  else if seas.eq."ALL" then 
    iseason = ind(mm)
  else if seas.eq."DJFMAM" then
    iseason = ind((mm.eq.12) .or. ((mm.ge.1).and.(mm.le.5)))

  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if
  end if

return (iseason)
end

;******************************************************************************
; cd_string.ncl
;
; This is based on ut_string.ncl. The only difference between these
; two scripts is that ut_string uses ut_calendar, and cd_string uses
; cd_calendar.
;
; If you make changes to this function, be sure to change ut_string in
; the same fashion, if appropriate.
;
; Carl Schreck (carl@atmos.albany.edu)
; Based on time_axis_labels (originally "tTick") by Arindam Chakraborty.
; Obtained from http://www.ncl.ucar.edu/Support/talk_archives/2007/0763.html
; September 2008
; Modified April 2009 by Carl Schreck (carl@atmos.albany.edu)
;   Replaced "mod" function call on line 279 with modulus operator %
;   Add default input string for ""
; Modified October 2010 by Carl Schreck (carl.schreck@noaa.gov)
;   Add support for "day of year"
; Modified April 2011 by Carl Schreck (carl.schreck@noaa.gov)
;   Now handles missing values without error messages, and removed "ChangeCase"
;******************************************************************************

undef("cd_string")

function cd_string( time[*]:numeric, inFormat:string )
;******************************************************************************
; Purpose: Convert a time into a string
;
; Arguments:
; time: The variable containing the values of time. This variable should
; have an attribute named "units".
;
; inFormat: (OPTIONAL) a string specifying the format of the tick labels.
; This string is parsed as follows: the '%' acts as
; the escape character. The single character after every
; '%' is formatted according to the rule:
; Y => 4-digit year (e.g., 2007).
; y => 2-digit year (e.g., 07).
; C => CAPITAL month abbreviation (e.g., JUN).
; c => Small month abbreviation (e.g., Jun).
; F => CAPITAL full month (e.g., JUNE).
; f => Small full month (e.g., June).
; N => 2-digit month (e.g., 06).
; n => 1 or 2 digit month (e.g., 6 for June, 12 for December).
; D => 2-digit day (e.g., 04).
; d => 1 or 2 digit day (e.g., 4)
; J => 3-digit day-of-year (e.g., 091)
; j => 1, 2, or 3 digit day-of-year (e.g., 4, 91, or 181)
; H => 2-digit hour (e.g., 09).
; h => 1 or 2 digit hour (e.g., 9 or 11).
; M => 2 digit minute (e.g., 08).
; m => 1 or 2 digit minute (e.g., 07 or 56).
; S => 2 digit second (e.g., 02).
; s => 1 or 2 digit second (e.g., 2 or 23).
; Any character at any other place in the format string
; is drawn as is.
; NOTE: a '%' can be drawn using "%%".
; In case Format is absent, a minimal algorithm exists
; which tries to determine the format string depending
; on the length and values of the date-time.
;
; If inFormat is "", it will default to "%H%M UTC %d %c %Y"
; e.g., 1800 UTC 4 Jul 1776
;
; Return value:
;  a string represent the time(s) in the specified format
;
; Exmple 1:
;  timeDouble = cd_inv_calendar( 1776, 07, 04, 18, 00, 00,  "hours since 1800-01-01 00:00", 0 )
;  print( (/ cd_string( timeDouble, "" ) /) )
;  print( (/ cd_string( timeDouble, "%Y-%N-%D %H:%M:%S" ) /) )
;
;  The above will output:
; (0)     1800 UTC 4 Jul 1776
; (0)     1776-07-04 18:00:00
;
;******************************************************************************

local ymdhms0, ymdhms1i, chFormat, months, fmonths, inside, stFormat

begin

  if( inFormat.eq."" ) then
    stFormat = "%H%M UTC %d %c %Y"
  else
    stFormat = inFormat
  end if

  months = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
  fmonths = (/"","January","February","March","April","May","June","July","August","September","October","November","December"/)

  ; check input
  if (.not. isatt(time,"units")) then
    print("cd_string: The variable containing time values must have a 'units' attribute.")
    return
  end if

  ymdhms0 = cd_calendar(time,0)
  nt = dimsizes(time)

  if (typeof(ymdhms0) .eq. "float") then
    ymdhms1i = floattointeger(ymdhms0)
  else
    ymdhms1i = ymdhms0
  end if

  retVal = new(nt,"string")

  chFormat = stringtochar(stFormat)

  do i = 0, dimsizes(time)-1
    retVal(i) = ""
    inside = False
    if( .not.ismissing(time(i)) ) then
      doy = day_of_year( ymdhms1i(i,0), ymdhms1i(i,1), ymdhms1i(i,2) )
    else
      doy = default_fillvalue("integer")
    end if
    j = 0
    do while(j .lt. dimsizes(chFormat))
      if (chFormat(j) .eq. "%") then
        inside = True
        j = j + 1
      end if

      if (.not. inside) then
        retVal(i) = retVal(i)+chFormat(j)
      else
        if (chFormat(j) .eq. "Y") then
          retVal(i) = retVal(i)+sprinti("%04d",ymdhms1i(i,0))
          inside = False
        end if
        if (chFormat(j) .eq. "y") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,0)%100)
          inside = False
        end if
        if (chFormat(j) .eq. "C") then
          if( .not.ismissing(time(i)) ) then
            retVal(i) = retVal(i)+str_upper(months(ymdhms1i(i,1)))
          else
            retVal(i) = retVal(i)+default_fillvalue("integer")
          end if
          inside = False
        end if
        if (chFormat(j) .eq. "c") then
          if( .not.ismissing(time(i)) ) then
            retVal(i) = retVal(i)+months(ymdhms1i(i,1))
          else
            retVal(i) = retVal(i)+default_fillvalue("integer")
          end if
          inside = False
        end if
        if (chFormat(j) .eq. "F") then
          if( .not.ismissing(time(i)) ) then
            retVal(i) = retVal(i)+str_upper(fmonths(ymdhms1i(i,1)))
          else
            retVal(i) = retVal(i)+default_fillvalue("integer")
          end if
          inside = False
        end if
        if (chFormat(j) .eq. "f") then
          if( .not.ismissing(time(i)) ) then
            retVal(i) = retVal(i)+fmonths(ymdhms1i(i,1))
          else
            retVal(i) = retVal(i)+default_fillvalue("integer")
          end if
          inside = False
        end if
        if (chFormat(j) .eq. "N") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,1))
          inside = False
        end if
        if (chFormat(j) .eq. "n") then
          retVal(i) = retVal(i)+sprinti("%d",ymdhms1i(i,1))
          inside = False
        end if
        if (chFormat(j) .eq. "D") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,2))
          inside = False
        end if
        if (chFormat(j) .eq. "d") then
          retVal(i) = retVal(i)+sprinti("%d",ymdhms1i(i,2))
          inside = False
        end if
        if (chFormat(j) .eq. "J") then
          retVal(i) = retVal(i)+sprinti("%03d",doy)
          inside = False
        end if
        if (chFormat(j) .eq. "j") then
          retVal(i) = retVal(i)+sprinti("%d",doy)
          inside = False
        end if
        if (chFormat(j) .eq. "H") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,3))
          inside = False
        end if
        if (chFormat(j) .eq. "h") then
          retVal(i) = retVal(i)+sprinti("%d",ymdhms1i(i,3))
          inside = False
        end if
        if (chFormat(j) .eq. "M") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,4))
          inside = False
        end if
        if (chFormat(j) .eq. "m") then
          retVal(i) = retVal(i)+sprinti("%d",ymdhms1i(i,4))
          inside = False
        end if
        if (chFormat(j) .eq. "S") then
          retVal(i) = retVal(i)+sprinti("%02d",ymdhms1i(i,5))
          inside = False
        end if
        if (chFormat(j) .eq. "s") then
          retVal(i) = retVal(i)+sprinti("%d",ymdhms1i(i,5))
          inside = False
        end if
        if (chFormat(j) .eq. "%") then
          retVal(i) = retVal(i)+"%"
          inside = False
        end if
      end if

      j = j + 1
    end do

  end do

  return( retVal )
end

undef("coast_get")
function coast_get (x[*][*][*]:numeric,face[*]:numeric,ncoast[*]:numeric)

; face=3 is south to north: south is ocean
; face=2 is north to south: north is ocean
; face=1 is east to west: east is ocean
; face=0 is west to east: west is ocean

  begin


  time = x&time
  lons = x&lon
  lats = x&lat

  ntime = dimsizes(time)
  nlats = dimsizes(lats)
  nlons = dimsizes(lons)

  icoast = new((/ntime, nlats,nlons/), "float")

    if (face.eq.2.or.face.eq.3) then

      if face.eq.2 then
        if lats(0).lt.lats(dimsizes(lats)-1) then
          print("Latitude from south to north, changing it to north to south")
          x = x(:,::-1,:)
          lats = x&lat
        else
          print("Latitude already from north to south, go on...")
        end if
      else if face.eq.3 then
        if lats(0).gt.lats(dimsizes(lats)-1) then
          print("Latitude from north to south, changing it to south to north")
          x = x(:,::-1,:)
          lats = x&lat
        else
          print("Latitude already from south to north, go on...")
        end if
      end if
      end if
      ; print(lats)


        do i = 0, dimsizes(lons)-1
          do j = 0, dimsizes(lats)-1
            if any(.not.ismissing(x(:,j,i))) then

                do ninit = 0, ncoast-1
                  icoast(:,j+ninit,i) = 1.
                end do

              break
            end if
          end do
        end do

      outcoast = where(icoast.eq.1, x(:,:,:), x@_FillValue)
      copy_VarMeta(x(:,:,:), outcoast)

      return(outcoast)

    else if (face.eq.0.or.face.eq.1) then

      if face.eq.0 then
        if lons(0).gt.lons(dimsizes(lons)-1) then
          print("Longitude from east to west, changing it to west to east")
          x = x(:,:,:)
          lons = x&lon
        else
          print("Longitude already from west to east, go on...")
        end if
      else if face.eq.1 then
        if lons(0).lt.lons(dimsizes(lons)-1) then
          print("Longitude from west to east, changing it to east to west")
          x = x(:,:,:)
          lons = x&lon
          print(lons)
        else
          print("Longitude already from east to west, go on...")
        end if
      end if
      end if
      ; print(lats)


        do j = 0, dimsizes(lats)-1
          do i = 0, dimsizes(lons)-1
            if any(.not.ismissing(x(:,j,i))) then

              do ninit = 0, ncoast-1
                icoast(:,j,i+ninit) = 1.
              end do

              break
            end if
          end do
        end do

      outcoast = where(icoast.eq.1, x(:,:,:), x@_FillValue)
      copy_VarMeta(x(:,:,:), outcoast)

      return(outcoast)


    end if
    end if




  end



; ----------------------------------------------------------
; Adam Phillips
;
undef("landsea_mask")
function landsea_mask (basemap[*][*]:byte,tlat:numeric,tlon:numeric)
;
;	Given a latitude and longitude array, this function will
;	return a field that contains a land sea mask. It is up to
;	the user to apply the mask to his/her data. The data contains
;	five types of values: 0=ocean, 1=land, 2=lake, 3=small island, 4=ice shelf
;
;	If the user is not happy with the default data file, the user may modify it
;	for their own use. Note! This function expects the input array to be of type byte.
;
;	2D lat/lon fields are accepted. integer/float/doubles are accepted
;	longitude arrays accepted from -180:359.999, if your longitudes
;	are < -180, use the contributed function lonFlip to convert
;	your longitudes to 0->360...
;
;	For all four examples:
;	1.  It is assumed array@_FillValue is set for use in the ncl mask function.
;	2.  The data file must be loaded prior to the function call. The data
;	    file is located in: $NCARG_ROOT/lib/ncarg/data/cdf
;	    To read in the data file:
;	       data_d = ncargpath("data") + "/cdf/landsea.nc"
;    	       f        = addfile(data_d,"r")
;	       bmap  = f->LSMASK
;
;	    The data file may also be downloaded from the ncl web site at:
;	    http://www.ncl.ucar.edu/Document/Functions/list_func_doc.shtml
;	    and read into your script.
;	3.  "bmap" is the name of the data array that contains the original land-sea mask.
;	    (see 2. above)
;
;	Example 1 array is a 3D array w/coord vars (time,lat,lon)
;		  Land,lakes, small islands, and ice shelves will be masked.
;		  lsm = landsea_mask(bmap,array&lat,array&lon)
;	          array = mask(array,conform(array,lsm,(/1,2/)).ge.1,False)
;
;	Example 2 array is a 2D array w/coord vars (lat,lon), land/small islands
;		  will be masked, but not lakes
;		  lsm = landsea_mask(bmap,array&lat,array&lon)
;	          array = mask(array,lsm.eq.1,False)
;		  array = mask(array,lsm.eq.3,False)
;
;	Example 3 array is a 2D array w/2D coordinate variables TLAT/TLON
;		  Land,lakes, small islands, and ice shelves will be masked.
;		  lsm = landsea_mask(bmap,TLAT,TLON)
;		  array = mask(array,lsm.ge.1,False)
;
;	Example 4 single point: lat1 and lon1 represent a single point.
;	          lsm = landsea_mask(bmap,lat1,lon1)
;	          print("Value returned = "+lsm)
;
begin

        rank = dimsizes(dimsizes(tlon))         ; # of dimensions

	if (rank.eq.2) then
	   ulon = ndtooned(tlon)		;convert 2D lat/lon coords to 1D
	   ulat = ndtooned(tlat)
	   if (any(dimsizes(ulat).ne.dimsizes(ulon))) then
	      print("landsea_mask:")
	      print("Fatal error: 2D lat/lon arrays are different sizes, exiting")
	      exit
	   end if
	   rdim = 1
	else
	   ulon = tlon
	   ulat = tlat
	   rdim = 0
	end if

	ulat = ulat + 90			;add 90 to convert to index
	if (min(ulon).lt.0) then
	   ind_lt0 = ind(ulon.lt.0)
	   ulon(ind_lt0) = ulon(ind_lt0)+360	;add 360 to those lons < 0
	   delete(ind_lt0)
	end if

	if (typeof(ulon).eq."integer") then
	   inlon = ulon
	end if
	if (typeof(ulat).eq."integer") then
           inlat = ulat
	end if
	if (typeof(ulon).eq."float") then        ;lop off decimal pt.
	   inlon = floattointeger(ulon)		 ;I want 56.0->56.99999 to refer to data
	end if					 ;pt at 56.5
	if (typeof(ulat).eq."float") then
	   inlat = floattointeger(ulat)
	end if
	if (typeof(ulon).eq."double") then
	   inlon = doubletointeger(ulon)
	end if
	if (typeof(ulat).eq."double") then
	   inlat = doubletointeger(ulat)
	end if

	if (min(inlat).lt.0) then
	   ind_lt0 = ind(inlat.lt.0)		;correction for all indices < 0 and > 179
	   inlat(ind_lt0) = 0
	   delete(ind_lt0)
	end if

	if (max(inlat).gt.179) then
	   ind_gt179 = ind(inlat.gt.179)
	   inlat(ind_gt179) = 179
	   delete(ind_gt179)
	end if

	if (max(inlon).ge.360) then
	   if (max(inlon).ge.361) then
	      print("Fatal error: invalid longitude point > 360.999, outside valid longitude range")
	      print("Exiting")
	      exit
	   end if
	   ind_ge360 = ind(inlon.ge.360.)
	   inlon(ind_ge360) = 0
	   delete(ind_ge360)
	end if

	if (rdim.eq.0) then
	   if (dimsizes(inlat).eq.1.and.dimsizes(inlon).eq.1) then		;user passed in 1 pt
	      newarr = basemap(inlat,inlon)
	   else
              ninlat = dimsizes(inlat)
              ninlon = dimsizes(inlon)
	      newarr = new((/ninlat,ninlon/),"byte")
	      newarr!0   = "lat"
	      newarr&lat = tlat
	      newarr!1   = "lon"
	      newarr&lon = tlon
	      do ss = 0,ninlat-1
	         newarr(ss,: ) = (/ basemap(inlat(ss),inlon(:)) /)
	      end do
	   end if
	end if

	if (rdim.eq.1) then
	   tarr = new((/dimsizes(inlat)/),"byte")

	  ;wcStrt = systemfunc("date")
	   do ss = 0,dimsizes(inlat)-1
	      tarr(ss) = (/ basemap(inlat(ss),inlon(ss)) /)
	   end do
	  ;wallClockElapseTime(wcStrt, "ss do loop", 0)
	   dimt = dimsizes(tlat)
	   newarr = onedtond(tarr,(/dimt(0),dimt(1)/))
	   newarr@lat2d = tlat
	   newarr@lon2d = tlon
	   delete(tarr)
	end if
	return(newarr)
end


undef ("name_season")
function name_season(mon)
begin
  month_abbr = (/"","NDJFMA","NDJFMA","NDJFMA","NDJFMA","MJJASO", "MJJASO","MJJASO","MJJASO","MJJASO","MJJASO","NDJFMA","NDJFMA"/)
  ; print(mon)
  print("season is...."+(month_abbr(toint(mon))))
  return(month_abbr(toint(mon)))

end
