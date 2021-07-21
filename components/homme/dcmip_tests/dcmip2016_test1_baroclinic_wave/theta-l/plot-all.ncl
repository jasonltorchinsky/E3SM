;================================================;
;	Example ncl script to produce the set of
;	lat-lon plots for DCMIP-2016 test case 1
; Created by James Kent, Christiane Jablonowski
;	and Paul Ullrich (University of Michigan) for DCMIP-2012 
;
; Modified by Colin Zarzycki for DCMIP-2016
; Modified by Jason Torchinsky to create a single combined PDF.
;
; v1.01 - 06/07/2016 
; v1.02 - 06/08/2016 - Typo in var_choice options
; v1.03 - 07/19/2021 - Combined previous plotting scripts.
;================================================;

; NOTE, if your model does not output T850 OR P at level
; midpoints OR you cannot calculate P from hybrid levels, you
; will need to make model specific additions denoted by
; CALCPHERE below

; PRECL needs to be in m/s for contours to plot correctly

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

;************************************************
begin

  print("Creating concatenated plot PDF...")

;************************************************
; Input parameters
;************************************************

  datafilename = "./movies/dcmip2016_test11.nc"

;************************************************
; Initialize workspace
;************************************************

  wks_type                  = "pdf"
  wks_type@wkOrientation    = "portrait"
  wks = gsn_open_wks(wks_type,"plot_all")

;************************************************
; Read in data common to all plots
;************************************************

  datafile = addfile(datafilename, "r")

  lat   = datafile->lat
  lon   = datafile->lon
  lev   = datafile->lev
  time  = datafile->time

  nlat  = dimsizes(lat)
  nlon  = dimsizes(lon)
  nlev  = dimsizes(lev)
  nt    = dimsizes(time)

  day6   = closest_val( 6, time)
  day9   = closest_val( 9, time)
  day12  = closest_val( 12, time)
  day15  = closest_val( 15, time)

;************************************************
; Create the initial condition plot
;************************************************

  res                     = True
  res@cnFillOn            = True       ; color              [default=False]
  res@cnLinesOn           = True       ; draw contour lines [default=True]
  res@lbLabelAutoStride   = True
  res@vpKeepAspect        = True
  res@vpWidthF            = 0.38       ; user specified shape
  res@vpHeightF           = 0.19
  res@gsnSpreadColors     = True       ; use all colors
  res@gsnSpreadColorStart = 2          ; default=2
  res@gsnSpreadColorEnd   = 23        ; final color of "gui-default"

  res@gsnDraw             = False      ; default=True
  res@gsnFrame            = False      ; default=True

  res_p                   = res        ; pressure plot resources
  res_ps                  = res        ; surface pressure plot resources
  res_u                   = res        ; zonal velocity plot resources
  res_up                  = res        ; zonal velocity pert plot resources
  res_t                   = res        ; temperature plot resources
  res_qv                  = res        ; water vapor plot resources
  res_theta               = res        ; potential temperature plot resources

;************************************************
; Panel plot
;************************************************
  print("Creating initial condition plot...")
  plot = new(6,graphic)
  gsn_merge_colormaps(wks,"gui_default", "BlRe")

  t     = datafile->T(0,::-1,:,:)
  theta = datafile->Th(0,::-1,:,:)
  u     = datafile->u(0,::-1,:,:)
  qv    = datafile->Q(0,::-1,:,:)
  p     = datafile->pnh(0,::-1,:,:)
  geo   = datafile->geo(0,::-1,:,:)
  qv    = qv*1000.0
  p     = p/100.0
  ps    = datafile->ps(0,:,:)

  delete(u@long_name)
  delete(p@long_name)
  delete(qv@long_name)
  delete(theta@long_name)
  delete(t@long_name)

  print("lon(0)="+lon(0))
  if(lon(0) .eq. 0) then
    print("lonFlip")
    t     = lonFlip(t)
    theta = lonFlip(theta)
    u     = lonFlip(u)
    qv    = lonFlip(qv)
    p     = lonFlip(p)
    ps    = lonFlip(ps)
    geo   = lonFlip(geo)
  end if

  up = u
  do i=0,nlon-1
    up(0,:,i) = up(0,:,i) - u(0,:,0)
  end do

  ;_________________________________________________
  ; Interpolate fields from eta coord to z coord

  z_eta = geo/9.81
  z     = lev
  z     = ispan(500,29500,1000)

  print("interpolating from geopotential levels to evenly spaced z")
  theta = linint1_n_Wrap(z_eta,theta,False,z,0,0)
  t     = linint1_n_Wrap(z_eta,t ,False,z,0,0)
  u     = linint1_n_Wrap(z_eta,u ,False,z,0,0)
  qv    = linint1_n_Wrap(z_eta,qv,False,z,0,0)
  p     = linint1_n_Wrap(z_eta,p ,False,z,0,0)

; ---
;  z         = lev * 30.0
  p&lev     = z/1000
  t&lev     = z/1000
  theta&lev = z/1000
  qv&lev    = z/1000
  u&lev     = z/1000

; ---
  res_t@cnLevelSelectionMode= "ManualLevels"
  res_t@cnLevelSpacingF   = 10.
  res_t@cnMinLevelValF    = 180.
  res_t@cnMaxLevelValF    = 310.

  res_t@tiYAxisString     = "Altitude (km)"
  res_t@gsnCenterString   = "Temperature"
  res_t@gsnRightString    = "(K)"
  plot(0) = gsn_csm_contour(wks,t(:,:,0),res_t)

  res_theta@cnLevelSelectionMode= "ExplicitLevels"
  res_theta@cnLevels = (/250., 260., 270., 280., 290., 300., 310., 320., 330., 340., 350., 360., 380., 400., 450., 500., 550., 600., 650., 700./)

  res_theta@tiYAxisString     = "Altitude (km)"
  res_theta@gsnCenterString   = "Potential Temperature"
  res_theta@gsnRightString    = "(K)"
  plot(1) = gsn_csm_contour(wks,theta(:,:,0),res_theta)

; ---

  res_u@cnLevelSelectionMode= "ManualLevels"
  res_u@cnLevelSpacingF   = 2.0
  res_u@cnMinLevelValF    = 2.0
  res_u@cnMaxLevelValF    = 26.0

  res_u@tiYAxisString     = "Altitude (km)"
  res_u@gsnCenterString   = "Zonal velocity"
  res_u@gsnRightString    = "(m/s)"
  plot(2) = gsn_csm_contour(wks,u(:,:,0),res_u)

  res_qv@cnLevelSelectionMode= "ManualLevels"
  res_qv@cnLevelSpacingF   = 2.0
  res_qv@cnMinLevelValF    = 2.0
  res_qv@cnMaxLevelValF    = 18.0

  res_qv@tiYAxisString     = "Altitude (km)"
  res_qv@gsnCenterString   = "Specific Humidity"
  res_qv@gsnRightString    = "(g/kg)"
  plot(3) = gsn_csm_contour(wks,qv(:,:,0),res_qv)

; ---

  res_up@cnLevelSelectionMode= "ManualLevels"
  res_up@cnLevelSpacingF   = 0.1
  res_up@cnMinLevelValF    = 0.1
  res_up@cnMaxLevelValF    = 1.0

  res_up@trXMinF = 0.
  res_up@trXMaxF = 40.
  res_up@trYMinF = 30.
  res_up@trYMaxF = 50.

  res_up@tiYAxisString     = "Latitude"
  res_up@gsnCenterString   = "Zonal wind perturbation"
  res_up@gsnRightString    = "(m/s)"
  plot(4) = gsn_csm_contour(wks,up(0,:,:),res_up)

  res_p@cnLevelSelectionMode= "ExplicitLevels"
  res_p@cnLevels = (/1.e1, 2.e1, 5.e1, 1.e2, 2.e2, 3.e2, 4.e2, 5.e2, 6.e2, 7.e2, 8.e2, 9.e2, 1.e3/)

  res_p@tiYAxisString     = "Altitude (km)"
  res_p@gsnCenterString   = "Pressure"
  res_p@gsnRightString    = "(hPa)"
  plot(5) = gsn_csm_contour(wks,p(:,:,0),res_p) 

 ; res_ps@cnLevelSelectionMode= "ManualLevels"
 ; res_ps@cnLevelSpacingF   = 10.
 ; res_ps@cnMinLevelValF    = 950.
 ; res_ps@cnMaxLevelValF    = 1050.

 ; res_ps@gsnCenterString   = "Surfae Pressure"
 ; res_ps@gsnRightString    = "(hPa)"
 ; plot(5) = gsn_csm_contour(wks,ps(:,:),res_ps) 

; ---

  pres = True
  pres@gsnFrame = False
  pres@gsnPanelYWhiteSpacePercent = 5
  pres@gsnPanelBottom = 0.1

  gsn_panel(wks,plot,(/3,2/),pres)

  frame(wks)

  delete(plot)


;***********************************************
; Temperature plot
;***********************************************
  print("Creating temperature plot...")

  if (isfilevar(datafile, "T850")) then
    temp = datafile->T850
  else
    ; If don't have T at 850hPa, interpolate to it
    tempall  = datafile->T
    temp     = tempall(:,0,:,:)
    P        = datafile->pnh
    plevel   = 85000.0
    temp     = int2p_n(P, tempall, plevel, 2, 1)

    delete(tempall)
    delete(P)
    delete(plevel)
  end if

  temp6    = temp(day6, :,:)
  temp9    = temp(day9, :,:)
  temp12   = temp(day12,:,:)
  temp15   = temp(day15,:,:)

  ; Delete loaded data, don't need it any more
  delete(temp)

  ; Produce the plot
  plot  = new(4, graphic)  ; define plot - need 4 panels

  res1                      = True
  res1@gsnDraw              = False          		; panel plot 
  res1@gsnFrame             = False			; don't draw yet
  res1@cnFillOn             = True
  res1@cnLinesOn            = True
  res1@gsnSpreadColors      = True
  res1@lbLabelAutoStride    = True
  res1@gsnCenterString      = ""
  res1@tiMainString         = ""
  res1@vpWidthF             = 0.38
  res1@vpHeightF            = 0.19
  res1@cnLevelSelectionMode = "ManualLevels"
  res1@cnInfoLabelOn        = False			; don't give each
  res1@cnLineLabelsOn       = False			; panel its own
  res1@lbLabelBarOn         = False			; label bar

  pres                          = True
  pres@gsnPanelLabelBar         = True			; Communal label bar
  pres@gsnPanelLeft             = 0.1
  pres@gsnPanelRight            = 0.9
  pres@pmLabelBarOrthogonalPosF = -0.03
  pres@gsnFrame                 = False
  pres@lbLabelStride            = 1

  res1@sfXArray        = lon            			; uses lon as plot x-axis
  res1@sfYArray        = lat            			; uses lat for y axis
  res1@trYReverse      = False          			; reverses y-axis, false
  res1@tiYAxisString   = ""   	         			; y-axis title
  res1@tiXAxisString   = ""             			; x-axis title

  pltTitle="THETA-L Test 16-1, T850 hPa" 			; Plot title if required
  pres@txString = pltTitle

  ; Change the output type and name
  gsn_define_colormap(wks,"gui_default")  

  res1@cnMaxLevelValF  = 300.0				; max contour color label
  res1@cnMinLevelValF  = 230.0				; min contour color label
  res1@cnLevelSpacingF = 10.0				; contour color spacing			; choose a colormap

  res1@gsnCenterString = "t = 6"
  plot(0) = gsn_csm_contour(wks,temp6(:,:),res1)

  res1@gsnCenterString = "t = 9"
  plot(1) = gsn_csm_contour(wks,temp9(:,:),res1)

  res1@gsnCenterString = "t = 12"
  plot(2) = gsn_csm_contour(wks,temp12(:,:),res1)

  res1@gsnCenterString = "t = 15"
  plot(3) = gsn_csm_contour(wks,temp15(:,:),res1)

  gsn_panel(wks,plot,(/2,2/),pres)		; 2x2 plot

  ; Add latitude and longitude labels

  txres3                = True
  txres3@txAngleF       = 90.
  txres3@txFontHeightF  = 0.02
  gsn_text_ndc(wks,"Latitude",0.08,0.55,txres3)

  txres3@txAngleF       = 0.
  gsn_text_ndc(wks,"Longitude",0.5,0.26,txres3)

  ; Delete unnecessary data
  delete(temp6)
  delete(temp9)
  delete(temp12)
  delete(temp15)

  frame(wks)

;***********************************************
; Large-scale precipitation plot
;***********************************************
  print("Creating large-scale precipitation plot...")

  precl = datafile->precl
  precl = 8.64e7*precl

  precl6  = precl(day6, :,:)
  precl9  = precl(day9, :,:)
  precl12 = precl(day12,:,:)
  precl15 = precl(day15,:,:)

  ; Delete unnecessary data
  delete(precl)

  pltTitle = "THETA-L Test 16-1, Large-scale precipitation"
  pres@txString  = pltTitle
  res1@cnLinesOn = False

  gsn_define_colormap(wks,"gui_default")

  res1@cnMaxLevelValF  = 50.0				; max contour color label
  res1@cnMinLevelValF  = 0.0				; min contour color label
  res1@cnLevelSpacingF = 4.0				; contour color spacing			; choose a colormap

  res1@gsnCenterString = "t = 6"
  plot(0) = gsn_csm_contour(wks,precl6(:,:),res1)

  res1@gsnCenterString = "t = 9"
  plot(1) = gsn_csm_contour(wks,precl9(:,:),res1)

  res1@gsnCenterString = "t = 12"
  plot(2) = gsn_csm_contour(wks,precl12(:,:),res1)

  res1@gsnCenterString = "t = 15"
  plot(3) = gsn_csm_contour(wks,precl15(:,:),res1)

  gsn_panel(wks,plot,(/2,2/),pres)		; 2x2 plot

  ; Add latitude and longitude labels

  txres3                = True
  txres3@txAngleF       = 90.
  txres3@txFontHeightF  = 0.02
  gsn_text_ndc(wks,"Latitude",0.08,0.55,txres3)

  txres3@txAngleF       = 0.
  gsn_text_ndc(wks,"Longitude",0.5,0.26,txres3)

  ; Delete unnecessary data
  delete(precl6)
  delete(precl9)
  delete(precl12)
  delete(precl15)

  frame(wks)

;***********************************************
; Surface Pressure Plot
;***********************************************
  print("Creating surface pressure plot...")

  delete(ps)
  ps = datafile->ps
  ps = ps/100.0 ; Convert to hPa

  ps6  = ps(day6 ,:,:)
  ps9  = ps(day9 ,:,:)
  ps12 = ps(day12,:,:)
  ps15 = ps(day15,:,:)

  delete(ps)

  pltTitle="THETA-L Test 16-1, PS"
  pres@txString  = pltTitle
  res1@cnLinesOn = True 

  gsn_define_colormap(wks,"gui_default")  

  res1@cnMaxLevelValF  = 1020.0				; max contour color label
  res1@cnMinLevelValF  = 920.0				; min contour color label
  res1@cnLevelSpacingF = 10.0				; contour color spacing			; choose a colormap

  res1@gsnCenterString = "t = 6"
  plot(0) = gsn_csm_contour(wks,ps6(:,:),res1)

  res1@gsnCenterString = "t = 9"
  plot(1) = gsn_csm_contour(wks,ps9(:,:),res1)

  res1@gsnCenterString = "t = 12"
  plot(2) = gsn_csm_contour(wks,ps12(:,:),res1)

  res1@gsnCenterString = "t = 15"
  plot(3) = gsn_csm_contour(wks,ps15(:,:),res1)

  gsn_panel(wks,plot,(/2,2/),pres)		; 2x2 plot

  ; Add latitude and longitude labels

  txres3                = True
  txres3@txAngleF       = 90.
  txres3@txFontHeightF  = 0.02
  gsn_text_ndc(wks,"Latitude",0.08,0.55,txres3)

  txres3@txAngleF       = 0.
  gsn_text_ndc(wks,"Longitude",0.5,0.26,txres3)

  ; Delete unnecessary data
  delete(ps6)
  delete(ps9)
  delete(ps12)
  delete(ps15)
  
end
