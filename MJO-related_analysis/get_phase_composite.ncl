;***********************************************************
;
; mjoclivar_16.ncl
;
;***********************************************************
; Generate life cycle composites based upon daily phase space
; If the MJO_INDEX is < 1.0 it is not included
;
; Source: Eun-Pa Lim: Bureau of Meteorology, Australia
; July, 2016
;***********************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

begin
   latS    = -80
   latN    =  80

   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000

   pltSubTitle = "Anomalous: OLR, U850, V850"

   pltDir  = "./"+model+"/"                             ; plot directory
   pltType = "png"                            ; send graphics to PNG file
   pltName = "mjoclivar16_"+expname+"_4mon_lag6_baseonREOF"                      ; yrStrt+"_"+yrLast
  
   diro    =  "./"+model+"/"                            ; output directory   
   diri    =  "./"+model+"/anom/"                            ; input directory   

   filo    = "rlut.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filu    = "u850.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filv    = "w500.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filz    = "zg250.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filp    = "pr.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filt    = "tas.day.anomalies."+yrStrt+"-"+yrLast+".nc"
;************************************************
; create BandPass Filter
;************************************************
  ihp      = 2                             ; bpf=>band pass filter
  nWgt     = 201
  sigma    = 1.0                           ; Lanczos sigma
  fca      = 1./100.
  fcb      = 1./20.
  wgt      = filwgts_lanczos (nWgt, ihp, fca, fcb, sigma )

;***********************************************************
; Find the indicies (subscripts) corresponding to the start/end times
;***********************************************************

   f       = addfile (diri+filu , "r")                         
   TIME    = f->time                          ; days since ...
   YMD     = cd_calendar(TIME, -2)            ; entire (time,6)

   iStrt   = ind(YMD.eq.ymdStrt)              ; index start
   iLast   = ind(YMD.eq.ymdLast)              ; index last 
   delete(TIME)
   delete(YMD )

   time    = f->time(iStrt:iLast)             ; days since ...
   u       = f->u850_anom(iStrt:iLast,{latS:latN},:) 
;***********************************************************
; Read anomalies frpm other fields
;***********************************************************
   f       = addfile (diri+filv , "r")                         
   v       = f->w500_anom(iStrt:iLast,{latS:latN},:)

   f       = addfile (diri+filz , "r")
   z       = f->zg250_anom(iStrt:iLast,{latS:latN},:)

   f       = addfile (diri+filp , "r")
   p      = f->pr_anom(iStrt:iLast,{latS:latN},:)

   f       = addfile (diri+filo , "r")                         
   x       = f->rlut_anom(iStrt:iLast,{latS:latN},:) 

   f       = addfile (diri+filt , "r")
   t       = f->tas_anom(iStrt:iLast,{latS:latN},:)

   dimx    = dimsizes( x )
   ntim    = dimx(0)
   nlat    = dimx(1)
   nlon    = dimx(2)
;************************************************
; Apply the band pass filter to the original anomalies
;************************************************
  x   = wgt_runave_leftdim (x, wgt, 0)
  u   = wgt_runave_leftdim (u, wgt, 0)
  v   = wgt_runave_leftdim (v, wgt, 0)
  z   = wgt_runave_leftdim (z, wgt, 0)
  p   = wgt_runave_leftdim (p, wgt, 0)
  t   = wgt_runave_leftdim (t, wgt, 0)
;***********************************************************
; Open PC components file created in 'mjo_14.ncl'
;***********************************************************
  dirMJO  = "./"+model+"/"                             ; input directory   
  fMJO    = "MJO_PC_"+expname+"_baseonREOF.nc"                ; created in mjo_14.ncl
  f       = addfile (dirMJO+fMJO, "r")

;***********************************************************
; Find the indices corresponding to the start/end times
;***********************************************************
  TIME    = f->time                          ; days since ...
  YMD     = cd_calendar(TIME, -2)            ; entire (time,6)

  iStrt   = ind(YMD.eq.ymdStrt)              ; index start
  iLast   = ind(YMD.eq.ymdLast)              ; index last 
  delete(TIME)
  delete(YMD )
  delete(time)

;***********************************************************
; Read the data for the desired period
;***********************************************************
  pc1     = f->PC1(iStrt:iLast)
  pc2     = f->PC2(iStrt:iLast)
  mjo_indx= f->MJO_INDEX(iStrt:iLast)
  lag = 6
  nnn = dimsizes(pc1)
  mjo_indx(lag:)  = mjo_indx(0:nnn-lag-1)
  mjo_indx(0:lag-1) = 0.
  pc1(lag:)  = pc1(0:nnn-lag-1)
  pc2(lag:)  = pc2(0:nnn-lag-1)
  time    = pc1&time 

  ymdhms  = cd_calendar(time, 0)

  imon    = floattoint( ymdhms(:,1) )   ; convenience
  iday    = floattoint( ymdhms(:,2) )   ; subscripts must be integer

;***********************************************************
; Place each array into an appropriate array
;***********************************************************
  nPhase      = 8
  angBnd      = new( (/2,nPhase/), "float")
  angBnd(0,:) = fspan(  0,315,nPhase)
  angBnd(1,:) = fspan( 45,360,nPhase)

  r2d         = 180./(4.*atan(1.0))
  ang         = atan2(pc2,pc1)*r2d     ; phase space 

  nn          = ind(ang.lt.0)

  ang(nn)     = ang(nn) + 360          ; make 0 to 360

;----------------------------------------------------------
; 0 <= ang < 45 --> MJO Phase 5 (i.e. +ve PC1 & +ve PC2)
;
;  print(pc1(:19)+"  "+pc2(:19)+"  "+ang(:19))
;----------------------------------------------------------

  nDays       = new (nPhase, "integer")
  pLabel      = "P"+ispan(1,nPhase,1)+": "

;------------------------------------------------------------
; PLOTS
;------------------------------------------------------------
  ; pltPath = pltDir+pltName+"_rv_16"
  pltPath = pltDir+pltName
  wks  = gsn_open_wks(pltType,pltPath)
  plot = new(nPhase,graphic)              ; create graphic array

  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
 
  res@mpFillOn             = False        ; turn off map fill
  res@mpMinLatF            = latS         ; zoom in on map
  res@mpMaxLatF            = latN
  res@mpCenterLonF         = 210.
  res@cnFillOn             = True         ; turn on color fill
  res@cnFillPalette        = "ViBlGrWhYeOrRe" ; set color map
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's
  res@gsnScalarContour     = True         ; contour 3rd array 
  res@gsnMajorLatSpacing   = 15
  res@gsnMajorLonSpacing   = 60
  res@tmXBLabelFontHeightF = 0.01
  res@tmYLLabelFontHeightF = 0.01

                                          ; common contours 
 ;mnmxint = nice_mnmxintvl( min(x) , max(x), 16, False)
  res@cnLevelSelectionMode = "ManualLevels"
  res@cnMinLevelValF       =  -40         ; -100; mnmxint(0)
  res@cnMaxLevelValF       =   40         ;   80; mnmxint(1)
  res@cnLevelSpacingF      =    5         ;   20; mnmxint(2)
;print(res)

  res@vcMinDistanceF            = 0.01            ; thin the vector density
  res@vcRefMagnitudeF           = 2.0             ; define vector ref mag
  res@vcRefLengthF              = 0.025           ; define length of vec ref
  res@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
  res@vcRefAnnoArrowLineColor   = "black"         ; change ref vector color
  res@vcRefAnnoArrowUseVecColor = False           ; don't use vec color for ref

                                          ; panel plot only resources
  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@lbLabelFontHeightF  = 0.01
  resP@gsnPanelBottom      = 0.05         ; add some space at bottom
  resP@pmLabelBarWidthF    = 0.8          ; label bar width
  resP@pmLabelBarHeightF   = 0.05
  resP@gsnPanelFigureStringsFontHeightF = 0.0125  ; bit larger than default
 ;resP@pmLabelBarOrthogonalPosF = 0.015   ; move labelbar up a bit

  txres                  = True             
  txres@txFontHeightF    = 0.01
  txid = gsn_create_text(wks, pltSubTitle, txres)

  amres                  = True 
 ;amres@amParallelPosF   =  0.575  
  amres@amOrthogonalPosF =  0.75  
  amres@amJust           = "CenterCenter"
 ;amres@amResizeNotify   = True

;*******************************************
; Loop over each phase
;*******************************************
 res@gsnLeftString               = ""
 res@gsnRightString              = ""
 u_phase = new ( (/2,nPhase,nlat,nlon/),float) 
 v_phase = new ( (/2,nPhase,nlat,nlon/),float)
 z_phase = new ( (/2,nPhase,nlat,nlon/),float)
 pr_phase = new ( (/2,nPhase,nlat,nlon/),float)
 olr_phase = new ( (/2,nPhase,nlat,nlon/),float)
 tas_phase = new ( (/2,nPhase,nlat,nlon/),float)
do nSeason=1,2
    if (nSeason.eq.1) then
        resP@gsnPanelMainString  = yrStrt+"-"+yrLast+": May to Oct"
    else
        resP@gsnPanelMainString  = yrStrt+"-"+yrLast+": Nov to Apr"
    end if
  
  do n=0,nPhase-1

    na = n+4             ; temporary adjustment for 0 <= ang < 45 represents MJO phase 5 not MJO phase 1 
    if(na.gt.7) then
      na = na - 8
    end if
;    print(na)

     if (nSeason.eq.1) then
         nt = ind(mjo_indx.gt.1.0    .and.                     \
                 (imon.ge.6          .and. imon.le.9).and.    \
                  ang.ge.angBnd(0,n) .and. ang.lt.angBnd(1,n))
     else
         nt = ind(mjo_indx.gt.1.0    .and.                     \
                 (imon.ge.12         .or.  imon.le. 3).and.    \
                  ang.ge.angBnd(0,n) .and. ang.lt.angBnd(1,n))
     end if
     if (.not.all(ismissing(nt))) then
         xAvg = dim_avg_Wrap( x(lat|:,lon|:,time|nt) )
         uAvg = dim_avg_Wrap( u(lat|:,lon|:,time|nt) )
         vAvg = dim_avg_Wrap( v(lat|:,lon|:,time|nt) )
         zAvg = dim_avg_Wrap( z(lat|:,lon|:,time|nt) )
         prAvg = dim_avg_Wrap( p(lat|:,lon|:,time|nt) )
         tasAvg = dim_avg_Wrap( t(lat|:,lon|:,time|nt) )
         olr_phase(nSeason-1,na,:,:) = xAvg
         u_phase(nSeason-1,na,:,:) = uAvg
         v_phase(nSeason-1,na,:,:) = vAvg
         z_phase(nSeason-1,na,:,:) = zAvg
         pr_phase(nSeason-1,na,:,:) = prAvg
         tas_phase(nSeason-1,na,:,:) = tasAvg
         nDays(na) = dimsizes(nt)
     
         res@tmXBLabelsOn         = False    ; do not draw lon labels
         res@tmXBOn               = False    ;             lon tickmarks
         if (n.eq.(nPhase-1)) then           ; 
             res@tmXBLabelsOn     = True     ; draw lon labels
             res@tmXBOn           = True     ;          tickmarks
         end if

         plot(na) = gsn_csm_vector_scalar_map(wks,uAvg,vAvg,xAvg,res)
     end if
     delete(nt)                  ; will change next iteration
  end do

  resP@gsnPanelFigureStrings= pLabel+nDays
  gsn_panel(wks,plot,(/nPhase,1/),resP)     ; now draw as one plot
 end do


;************************************************
; Write PC results to netCDF for further analysis
;************************************************

  olr           = olr_phase
  olr!0         = "season"
  olr!1         = "phase"
  olr!2         = "lat"
  olr!3         =  "lon"
  olr&season   = (/1,2/)
  olr&phase    = (/1,2,3,4,5,6,7,8/)
  olr&lat      = x&lat
  olr&lon      = x&lon
  olr@long_name = "phase-avergaed OLR"

  u850           = u_phase
  u850!0         = "season"
  u850!1         = "phase"
  u850!2         = "lat"
  u850!3         =  "lon"
  u850&season   = (/1,2/)
  u850&phase    = (/1,2,3,4,5,6,7,8/)
  u850&lat      = x&lat
  u850&lon      = x&lon
  u850@long_name = "phase-avergaed u850"

  w500           = v_phase
  w500!0         = "season"
  w500!1         = "phase"
  w500!2         = "lat"
  w500!3         = "lon"
  w500&season   = (/1,2/)
  w500&phase    = (/1,2,3,4,5,6,7,8/)
  w500&lat      = x&lat
  w500&lon      = x&lon
  w500@long_name = "phase-avergaed w500"

  zg250           = z_phase
  zg250!0         = "season"
  zg250!1         = "phase"
  zg250!2         = "lat"
  zg250!3         = "lon"
  zg250&season   = (/1,2/)
  zg250&phase    = (/1,2,3,4,5,6,7,8/)
  zg250&lat      = x&lat
  zg250&lon      = x&lon
  zg250@long_name = "phase-avergaed zg250"

  pr           = pr_phase
  pr!0         = "season"
  pr!1         = "phase"
  pr!2         = "lat"
  pr!3         = "lon"
  pr&season   = (/1,2/)
  pr&phase    = (/1,2,3,4,5,6,7,8/)
  pr&lat      = x&lat
  pr&lon      = x&lon
  pr@long_name = "phase-avergaed Pr"

  tas           = tas_phase
  tas!0         = "season"
  tas!1         = "phase"
  tas!2         = "lat"
  tas!3         = "lon"
  tas&season   = (/1,2/)
  tas&phase    = (/1,2,3,4,5,6,7,8/)
  tas&lat      = x&lat
  tas&lon      = x&lon
  tas@long_name = "phase-avergaed tas"

  diro = "./"+model+"/"
  filo = "MJO_phase_composite_"+expname+"_4mon_lag6_baseonREOF.nc"
  system("/bin/rm -f "+diro+filo)   ; remove any pre-existing file
  ncdf = addfile(diro+filo,"c")     ; open output netCDF file
                                    ; make time an UNLIMITED dimension
  ncdf->olr       = olr
  ncdf->u850      = u850
  ncdf->w500      = w500
  ncdf->pr        = pr
  ncdf->zg250     = zg250
  ncdf->tas       = tas
end
