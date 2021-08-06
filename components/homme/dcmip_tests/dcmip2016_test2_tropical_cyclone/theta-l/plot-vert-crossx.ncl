;=======================================================================================
; This NCL code calculates radially-averaged tangential and radial wind components
; as well as T anomaly for DCMIP test case #2 (cyclone)
;=======================================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin

;=======================================================================================
; User options
;=======================================================================================

filename="./movies/dcmip2016_test21.nc"
f = addfile(filename,"r")


out_type="pdf"                ; output format, popular options are x11, png, pdf, eps


;=======================================================================================
; Get grid data from file
;=======================================================================================
print("Loading grid data from file...")

time = f->time
nt = dimsizes(time)-1

lat = f->lat
lon = f->lon
lev = f->lev
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

; Set days of interest
days = (/0, 2, 4, 6, 8, 10/)
ndays = dimsizes(days)
dayIdxs = new(ndays, "integer")
do ii = 0, ndays-1
   dayIdxs(ii) = closest_val(days(ii), time)
end do

; Set latitudes of interest
lats = fspan(15, 30, 5)
nlats = dimsizes(lats)
latIdxs = new(nlats, "integer")
do ii = 0, nlats-1
   latIdxs(ii) = closest_val(lats(ii), lat)
end do

;=======================================================================================
; Import state variable data from file
;=======================================================================================
print("Importing state variables...")
g = 9.807

geo   = f->geo(:,::-1,:,:)
u     = f->u(:,::-1,:,:)
v     = f->v(:,::-1,:,:)
w     = f->w(:,::-1,:,:)    
T     = f->T(:,::-1,:,:)    ; Temperature
pnh   = f->pnh(:,::-1,:,:)  ; Pressure
wvmr  = f->Q(:,::-1,:,:)    ; Water vapor mixing ratio

; Delete long names for plotting
delete(u@long_name)
delete(v@long_name)
;delete(w@long_name) ; w has no netCDF metadata
delete(T@long_name)
delete(pnh@long_name)
delete(wvmr@long_name)

; Interpolate fields from eta coord to z coord
print("Interpolating state variables from geopotential levels to evenly-spaced levels...")
z_eta = geo/g
z     = lev
z     = ispan(500, 29500, 1000)

do ii = 0, ndays-1
   dayIdx = dayIdxs(ii)
   do jj = 0, nlats-1
      latIdx = latIdxs(jj)
      do kk = 0, nlon-1
         u(dayIdx, :, latIdx, kk)    = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      u(dayIdx, :, latIdx, kk),   False, z, 0, 0)
         v(dayIdx, :, latIdx, kk)    = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      v(dayIdx, :, latIdx, kk),   False, z, 0, 0)
         w(dayIdx, :, latIdx, kk)    = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      w(dayIdx, :, latIdx, kk),   False, z, 0, 0)
         T(dayIdx, :, latIdx, kk)    = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      T(dayIdx, :, latIdx, kk),   False, z, 0, 0)
         pnh(dayIdx, :, latIdx, kk)  = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      pnh(dayIdx, :, latIdx, kk), False, z, 0, 0)
         wvmr(dayIdx, :, latIdx, kk) = linint1_n_Wrap(z_eta(dayIdx, :, latIdx, kk), \
                                                      wvmr(dayIdx, :, latIdx, kk), False, z, 0, 0)
     end do
   end do
end do

; Set the longitude to center the plot on, and adjust the state variable arrays
print("Centering state variables on longitude of interest...")
cenlon = 0. 

u    = lonPivot(u, cenlon)
v    = lonPivot(v, cenlon)
w    = lonPivot(w, cenlon)
T    = lonPivot(T, cenlon)
pnh  = lonPivot(pnh, cenlon)
wvmr = lonPivot(wvmr, cenlon)

; Set maximum altitude
maxalt = 10000.
maxaltIdx = closest_val(maxalt, z)

;=======================================================================================
; Set up graphics options common to all plots
;=======================================================================================
wks        = gsn_open_wks(out_type, "vert-crossx")

res                    = True
res@gsnDraw            = False
res@gsnFrame           = False
res@gsnMaximize        = True
res@cnFillOn           = True        ; Turn on contour fill
res@cnLinesOn          = False
res@cnLineLabelsOn     = False
res@cnInfoLabelOn      = False
res@cnConstFEnableFill = True
res@gsnAddCyclic       = False

res@sfYArray           = z
res@trYMaxF            = maxalt

;=======================================================================================
; Draw zonal velocity plots
;=======================================================================================
print("Drawing zonal velocity plots...")

u_contours = new(ndays, graphic)

res_u = res

res_u@cnFillPalette        = "temp_diff_18lev"
res_u@cnSpanFillPalette    = True

res_u@tiYAxisString     = "Altitude (m)"
res_u@gsnRightString    = "(m/s)"

pres = True
pres@gsnFrame = True
pres@gsnPanelYWhiteSpacePercent = 5
pres@gsnPanelBottom = 0.1

; Get bounds for colorbars considering all latitudes of interest
prevcnBnd = 0.d
cnBnd     = 0.d
do ii = 0, nlats-1
   prevcnBnd = cnBnd
   cnBnd     = max( (/prevcnBnd, max(abs(u(:, 0:maxaltIdx, latIdxs(ii), :)))/) )
end do

cnSpc = (2.d*cnBnd)/20.d
res_u@cnLevelSelectionMode = "ManualLevels"
res_u@cnLevelSpacingF      = ceil(cnSpc)
res_u@cnMinLevelValF       = -ceil(cnBnd)
res_u@cnMaxLevelValF       = ceil(cnBnd)

do ii = 0, nlats-1
   pres@txString = "Zonal Wind Velocity - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_u@gsnCenterString = "Day " + tostring(days(jj))
      u_contours(jj) = gsn_csm_contour(wks, u(dayIdxs(jj), :, latIdxs(ii), :), res_u)
   end do

   gsn_panel(wks, u_contours, (/ndays/2,2/), pres)
end do

; Clear stuff that doesn't need to be used any more
delete(u_contours)
delete(res_u)

;=======================================================================================
; Draw meridional velocity plots
;=======================================================================================
print("Drawing meridional velocity plots...")

v_contours = new(ndays, graphic)

res_v = res

res_v@cnFillPalette        = "temp_diff_18lev"
res_v@cnSpanFillPalette    = True

res_v@tiYAxisString     = "Altitude (m)"
res_v@gsnRightString    = "(m/s)"

pres = True
pres@gsnFrame = True
pres@gsnPanelYWhiteSpacePercent = 5
pres@gsnPanelBottom = 0.1

; Get bounds for colorbars considering all latitudes of interest
prevcnBnd = 0.d
cnBnd     = 0.d
do ii = 0, nlats-1
   prevcnBnd = cnBnd
   cnBnd     = max( (/prevcnBnd, max(abs(v(:, 0:maxaltIdx, latIdxs(ii), :)))/) )
end do

cnSpc = (2.d*cnBnd)/20.d
res_v@cnLevelSelectionMode = "ManualLevels"
res_v@cnLevelSpacingF      = ceil(cnSpc)
res_v@cnMinLevelValF       = -ceil(cnBnd)
res_v@cnMaxLevelValF       = ceil(cnBnd)

do ii = 0, nlats-1
   pres@txString = "Meridional Wind Velocity - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_v@gsnCenterString = "Day " + tostring(days(jj))
      v_contours(jj) = gsn_csm_contour(wks, v(dayIdxs(jj), :, latIdxs(ii), :), res_v)
   end do

   gsn_panel(wks, v_contours, (/ndays/2,2/), pres)
end do

; Clear stuff that doesn't need to be used any more
delete(v_contours)
delete(res_v)

;=======================================================================================
; Draw vertical velocity plots
;=======================================================================================
print("Drawing vertical velocity plots...")

w = w*3600.d ; Convert to m/h

w_contours = new(ndays, graphic)

res_w = res

res_w@cnFillPalette        = "temp_diff_18lev"
res_w@cnSpanFillPalette    = True

res_w@tiYAxisString     = "Altitude (m)"
res_w@gsnRightString    = "(m/h)"

pres = True
pres@gsnFrame = True
pres@gsnPanelYWhiteSpacePercent = 5
pres@gsnPanelBottom = 0.1

; Get bounds for colorbars considering all latitudes of interest
prevcnBnd = 0.d
cnBnd     = 0.d
do ii = 0, nlats-1
   prevcnBnd = cnBnd
   cnBnd     = max( (/prevcnBnd, max(abs(w(:, 0:maxaltIdx, latIdxs(ii), :)))/) )
end do

cnSpc = (2.d*cnBnd)/20.d
res_w@cnLevelSelectionMode = "ManualLevels"
res_w@cnLevelSpacingF      = ceil(cnSpc)
res_w@cnMinLevelValF       = -ceil(cnBnd)
res_w@cnMaxLevelValF       = ceil(cnBnd)

do ii = 0, nlats-1
   pres@txString = "Vertical Wind Velocity - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_w@gsnCenterString = "Day " + tostring(days(jj))
      w_contours(jj) = gsn_csm_contour(wks, w(dayIdxs(jj), :, latIdxs(ii), :), res_w)
   end do

   gsn_panel(wks, w_contours, (/ndays/2,2/), pres)
end do

; Clear stuff that doesn't need to be used any more
delete(w_contours)
delete(res_w)

;=======================================================================================
; Draw temperature plots
;=======================================================================================
print("Drawing temperature plots...")

T_contours = new(ndays, graphic)

res_T = res

res_T@cnFillPalette        = "MPL_gnuplot"
res_T@cnSpanFillPalette    = True

res_T@gsnYAxisIrregular2Log = True

res_T@tiYAxisString     = "Altitude (m)"
res_T@gsnRightString    = "(K)"

; Get bounds for colorbars considering all latitudes of interest
cnBndU = 0.d
cnBndL = 10000.d
do ii = 0, nlats-1
   cnBndU = max( (/cnBndU, max(T(:, 0:maxaltIdx, latIdxs(ii), :))/) )
   cnBndL = min( (/cnBndL, min(T(:, 0:maxaltIdx, latIdxs(ii), :))/) )
end do

cnSpc = (cnBndU - cnBndL)/24.d
res_T@cnLevelSelectionMode = "ManualLevels"
res_T@cnLevelSpacingF      = ceil(cnSpc)
res_T@cnMinLevelValF       = floor(cnBndL)
res_T@cnMaxLevelValF       = ceil(cnBndU)

do ii = 0, nlats-1
   pres@txString = "Temperature - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_T@gsnCenterString = "Day " + tostring(days(jj))
      T_contours(jj) = gsn_csm_contour(wks, T(dayIdxs(jj), :, latIdxs(ii), :), res_T)
   end do
   gsn_panel(wks, T_contours, (/ndays/2,2/), pres)
end do


; Clear stuff that doesn't need to be used any more
delete(T_contours)
delete(res_T)

;=======================================================================================
; Draw pressure plots
;=======================================================================================
print("Drawing pressure plots...")

pnh = pnh/100. ; Convert from Pa to hPa

pnh_contours = new(ndays, graphic)

res_pnh = res

res_pnh@cnFillPalette        = "WhBlGrYeRe"
res_pnh@cnSpanFillPalette    = True

res_pnh@gsnYAxisIrregular2Log = True

res_pnh@tiYAxisString     = "Altitude (m)"
res_pnh@gsnRightString    = "(hPa)"

; Get bounds for colorbars considering all latitudes of interest
cnBndU = 0.d
cnBndL = 10000.d
do ii = 0, nlats-1
   cnBndU = max( (/cnBndU, max(pnh(:, 0:maxaltIdx, latIdxs(ii), :))/) )
   cnBndL = min( (/cnBndL, min(pnh(:, 0:maxaltIdx, latIdxs(ii), :))/) )
end do

cnSpc = (cnBndU - cnBndL)/24.d
res_pnh@cnLevelSelectionMode = "ManualLevels"
res_pnh@cnLevelSpacingF      = ceil(cnSpc)
res_pnh@cnMinLevelValF       = floor(cnBndL)
res_pnh@cnMaxLevelValF       = ceil(cnBndU)

do ii = 0, nlats-1
   pres@txString = "Pressure - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_pnh@gsnCenterString = "Day " + tostring(days(jj))
      pnh_contours(jj) = gsn_csm_contour(wks, pnh(dayIdxs(jj), :, latIdxs(ii), :), res_pnh)
   end do

   gsn_panel(wks, pnh_contours, (/ndays/2,2/), pres)
end do

; Clear stuff that doesn't need to be used any more
delete(pnh_contours)
delete(res_pnh)


;=======================================================================================
; Draw water vapor mixing ratio plots
;=======================================================================================
print("Drawing water vapor mixing ratio plots...")


wvmr_contours = new(ndays, graphic)

res_wvmr = res

res_wvmr@cnFillPalette        = "WhViBlGrYeOrRe"
res_wvmr@cnSpanFillPalette    = True

res_wvmr@gsnYAxisIrregular2Log = True

res_wvmr@tiYAxisString     = "Altitude (m)"
res_wvmr@gsnRightString    = "(kg/kg)"

; Set bounds for colorbars
cnBndU = 0.d
cnBndL = 1.d
do ii = 0, nlats-1
   cnBndU = max( (/cnBndU, max(wvmr(:, 0:maxaltIdx, latIdxs(ii), :))/) )
   cnBndL = min( (/cnBndL, min(wvmr(:, 0:maxaltIdx, latIdxs(ii), :))/) )
end do

cnBndU = decimalPlaces(cnBndU, 2, True)
cnBndL = decimalPlaces(cnBndL, 2, False)

cnSpc = decimalPlaces((cnBndU - cnBndL)/16.d, 3, True)
res_wvmr@cnLevelSelectionMode = "ManualLevels"
res_wvmr@cnMaxLevelCount      = 24
res_wvmr@cnLevelSpacingF      = cnSpc
res_wvmr@cnMinLevelValF       = cnBndL
res_wvmr@cnMaxLevelValF       = cnBndU

do ii = 0, nlats-1
   pres@txString = "Water Vapor Mixing Ratio - Latitude: " + sprintf("%5.2f", lats(ii)) + " Degrees"
   do jj = 0, ndays-1
      res_wvmr@gsnCenterString = "Day " + tostring(days(jj))
      wvmr_contours(jj) = gsn_csm_contour(wks, wvmr(dayIdxs(jj), :, latIdxs(ii), :), res_wvmr)
   end do

   gsn_panel(wks, wvmr_contours, (/ndays/2,2/), pres)
end do

; Clear stuff that doesn't need to be used any more
delete(wvmr_contours)
delete(res_wvmr)



end
