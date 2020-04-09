;******************************************************
;
; mjoclivar_14.ncl
;
;***********************************************************
; Combined EOFs
; Latest Update: July, 2016: Eun-Pa Lim; Bureau of Meteorology, Australia
;***********************************************************
;;
;;      The following are automatically loaded from 6.2.0 onward
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
;;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 

undef("read_rename")
function read_rename(f[1]:file, varName[1]:string       \
                    ,iStrt[1]:integer, iLast[1]:integer \
                    ,latS[1]:numeric , latN[1]:numeric  )
; Utility to force specific named dimensions
; This is done for historical reasons (convenience) 
begin
   print(model)
   work    = f->$varName$(iStrt:iLast,{latS:latN},:)   ; (time,lat,lon)
   work!0  = "time"                                    ; CAM model names
   work!1  = "lat"
   work!2  = "lon"
   return(work)
end
; =========================>  MAIN  <==============================
begin
;   model = "GFDL-CM4"
;   expname = "rcp85"
;   ymdStrt = 20750101                         ; start yyyymmdd
;   ymdLast = 21001231                         ; last

   neof    =  2

   latS    = -15
   latN    =  15

   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000

   pltDir  = "./"+model+"/"                             ; plot directory
   pltType = "png" 
   pltName = "mjoclivar14_"+expname   

   diri    = "./"+model+"/anom/"                             ; input directory

   filolr  = "rlut.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filu200 = "u250.day.anomalies."+yrStrt+"-"+yrLast+".nc"
   filu850 = "u850.day.anomalies."+yrStrt+"-"+yrLast+".nc"

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
; Find the indices corresponding to the start/end times
;***********************************************************
   f       = addfile (diri+filolr , "r")                         
   TIME    = f->time                          ; days since ...

   YMD     = cd_calendar(TIME, -2)            ; entire (time,6)

   iStrt   = ind(YMD.eq.ymdStrt)              ; index start
   iLast   = ind(YMD.eq.ymdLast)              ; index last 
   delete([/ TIME, YMD /])

;***********************************************************
; Read anomalies
;***********************************************************

   work    = read_rename(f,"rlut_anom",iStrt,iLast,latS,latN) ; (time,lat,lon) 
   OLR     = dim_avg_n_Wrap(work, 1)                         ; (time,lon)

   f       = addfile (diri+filu850 , "r")                         
   work    = read_rename(f,"u850_anom",iStrt,iLast,latS,latN) ; (time,lat,lon) 
   U850    = dim_avg_n_Wrap(work, 1)          ; (time,lon)

   f       = addfile (diri+filu200 , "r")                         
   work    = read_rename(f,"u250_anom",iStrt,iLast,latS,latN) ; (time,lat,lon) 
   U200    = dim_avg_n_Wrap(work, 1)          ; (time,lon)

   dimw    = dimsizes( work )
   ntim    = dimw(0)
   nlat    = dimw(1)
   mlon    = dimw(2)
   delete(work)

   lon     = OLR&lon                                          
   time    = OLR&time                         
   date    = cd_calendar(time, -2)            ; yyyymmdd

;************************************************
; Apply the band pass filter to the original anomalies
;************************************************
   olr   = wgt_runave_n_Wrap ( OLR, wgt, 0, 0) ; (time,lon)
   u850  = wgt_runave_n_Wrap (U850, wgt, 0, 0)
   u200  = wgt_runave_n_Wrap (U200, wgt, 0, 0)

;************************************************
; remove temporal means of band pass series: *not* necessary 
;************************************************
   olr   = dim_rmvmean_n( olr, 0)              ; (time,lon)
   u850  = dim_rmvmean_n(u850, 0)
   u200  = dim_rmvmean_n(u200, 0)

;************************************************
; Compute the temporal variance at each lon
;************************************************
   var_olr  = dim_variance_n_Wrap( olr, 0)     ; (lon)
   var_u850 = dim_variance_n_Wrap(u850, 0)
   var_u200 = dim_variance_n_Wrap(u200, 0)

;************************************************
; Compute the zonal mean of the temporal variance
;************************************************
  zavg_var_olr  = dim_avg_n_Wrap( var_olr , 0)     
  zavg_var_u850 = dim_avg_n_Wrap( var_u850, 0)
  zavg_var_u200 = dim_avg_n_Wrap( var_u200, 0)

;************************************************
; Normalize by sqrt(avg_var*)
;************************************************
  olr   =  olr/sqrt(zavg_var_olr )          ; (time,lon)
  u850  = u850/sqrt(zavg_var_u850)
  u200  = u200/sqrt(zavg_var_u200)

;************************************************
; Combine the normalized data into one variable
;************************************************
  cdata     = new ( (/3*mlon,ntim/), typeof(olr), getFillValue(olr))
  do ml=0,mlon-1
     cdata(ml       ,:) = (/  olr(:,ml) /)
     cdata(ml+  mlon,:) = (/ u850(:,ml) /)
     cdata(ml+2*mlon,:) = (/ u200(:,ml) /)
  end do

;************************************************
; Compute **combined** EOF; Sign of EOF is arbitrary
;************************************************
  eof_cdata    = eofunc(cdata   , neof, False)      ; (neof,3*mlon)
;  print("==============")
;  printVarSummary(eof_cdata)
;  printMinMax(eof_cdata, True)

  eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False)   ; (neof,3*ntim)
;  print("==============")                                  
;  printVarSummary(eof_ts_cdata)
;  printMinMax(eof_ts_cdata, True)

;************************************************
; For clarity, explicitly extract each variable. Create time series 
;************************************************

  nvar = 3  ; "olr", "u850", "u200"
  ceof = new( (/nvar,neof,mlon/), typeof(cdata), getFillValue(cdata))

  do n=0,neof-1
     ceof(0,n,:) = eof_cdata(n,0:mlon-1)      ; olr
     ceof(1,n,:) = eof_cdata(n,mlon:2*mlon-1) ; u850
     ceof(2,n,:) = eof_cdata(n,2*mlon:)       ; u200
  end do

  ceof!0   = "var"
  ceof!1   = "eof"
  ceof!2   = "lon"   
  ceof&lon =  olr&lon

  ceof_ts        = new( (/nvar,neof,ntim/), typeof(cdata), getFillValue(cdata))
  ceof_ts(0,:,:) = eofunc_ts_Wrap( olr(lon|:,time|:),ceof(0,:,:),False)   ; (0,neof,ntim)
  ceof_ts(1,:,:) = eofunc_ts_Wrap(u850(lon|:,time|:),ceof(1,:,:),False)   ; (1,neof,ntim)
  ceof_ts(2,:,:) = eofunc_ts_Wrap(u200(lon|:,time|:),ceof(2,:,:),False)   ; (2,neof,ntim)

;**********************************************t*
; Add code contributed by Marcus N. Morgan, Florida Institute of Technology; Feb 2015
; Calculate % variance (pcv_ )accounted for by OLR, U850 and U200
;************************************************

    pcv_eof_olr  = new(neof,typeof(ceof))
    pcv_eof_u850 = new(neof,typeof(ceof))
    pcv_eof_u200 = new(neof,typeof(ceof))
      
    do n=0,neof-1
       pcv_eof_olr(n)  = avg((ceof(0,n,:)*sqrt(ceof@eval(n)))^2)*100
       pcv_eof_u850(n) = avg((ceof(1,n,:)*sqrt(ceof@eval(n)))^2)*100
       pcv_eof_u200(n) = avg((ceof(2,n,:)*sqrt(ceof@eval(n)))^2)*100
     ;;print("pcv: neof="+(n+1)+":  "+pcv_eof_olr(n)+"  "+pcv_eof_u850(n)+"  "+pcv_eof_u200(n))
    end do

;************************************************
; Change sign of EOFs for spatial structures of PC1 and PC2 
; to represent convection over the tropical Indian Ocean and the tropical western Pacific Ocean, respectively 
; (Ad hoc approach) 
;************************************************
  do_flip1 = 0
  do_flip2 = 0
  do_swip = 0
  imax_olr_eof1   = maxind(ceof(0,0,:))
  imax_olr_eof2   = maxind(ceof(0,1,:))
  lonmax_eof1 = ceof&lon(imax_olr_eof1)      ; longitude of max value (i.e. suppressed convection)
  lonmax_eof2 = ceof&lon(imax_olr_eof2)

  imin_olr_eof1   = minind(ceof(0,0,:))
  imin_olr_eof2   = minind(ceof(0,1,:))
  lonmin_eof1 = ceof&lon(imin_olr_eof1)      ; longitude of max value (i.e. suppressed convection)
  lonmin_eof2 = ceof&lon(imin_olr_eof2)

  if (lonmin_eof1.ge.60 .and. lonmin_eof1.le.120 .and. lonmin_eof2.ge.120 .and. lonmin_eof2.le.200) then
    print("no change") 
  else if (lonmin_eof1.ge.60 .and. lonmin_eof1.le.120 .and. lonmax_eof2.ge.120 .and. lonmax_eof2.le.200) then
     do_flip2 = 1
  else if (lonmin_eof1.ge.120 .and. lonmin_eof1.le.130 .and. lonmax_eof2.ge.140 .and. lonmax_eof2.le.200) then
     do_flip2 = 1
  else if  (lonmax_eof1.ge.60 .and. lonmax_eof1.le.120 .and. lonmin_eof2.ge.120 .and. lonmin_eof2.le.200) then
     do_flip1 = 1
  else if (lonmax_eof1.ge.60 .and. lonmax_eof1.le.120 .and. lonmax_eof2.ge.110 .and. lonmax_eof2.le.200) then
    do_flip1 = 1
    do_flip2 = 1
  else if (lonmin_eof1.ge.120 .and. lonmin_eof1.le.200 .and. lonmin_eof2.ge.60 .and. lonmin_eof2.le.120) then
    do_swip = 1
  else if (lonmin_eof1.ge.120 .and. lonmin_eof1.le.200 .and. lonmax_eof2.ge.60 .and. lonmax_eof2.le.120) then
    do_flip2 = 1
    do_swip =1
  else if (lonmax_eof1.ge.120 .and. lonmax_eof1.le.200 .and. lonmin_eof2.ge.60 .and. lonmin_eof2.le.120) then
    do_flip1 = 1
    do_swip = 1
  else if (lonmax_eof1.ge.120 .and. lonmax_eof1.le.200 .and. lonmax_eof2.ge.60 .and. lonmax_eof2.le.120) then
    do_flip1 = 1
    do_flip2 = 1
    do_swip = 1
  else
    print("exception meet")
  end if
  end if
  end if 
  end if
  end if
  end if
  end if
  end if
  end if
  if (do_flip1.eq.1) then  ; Change the sign of EOF1 
      ceof(:,0,:)       = -ceof(:,0,:)                  ; if OLR is positive
      ceof_ts(:,0,:)    = -ceof_ts(:,0,:)               ;  over the tropical Indian Ocean
      eof_cdata(0,:)    = -eof_cdata(0,:)
      eof_ts_cdata(0,:) = -eof_ts_cdata(0,:)
  end if

  if (do_flip2.eq.1) then
      ceof(:,1,:)       = -ceof(:,1,:)                   ; if OLR is positive 
      ceof_ts(:,1,:)    = -ceof_ts(:,1,:)                ; over the tropical western Pacific Ocean
      eof_cdata(1,:)    = -eof_cdata(1,:)
      eof_ts_cdata(1,:) = -eof_ts_cdata(1,:)
  end if


  if (do_swip.eq.1) then  ; swip eof1 and eof2
  
      tmp_8      = ceof@pcvar(0)
      tmp_9      = ceof@pcvar(1) 

      tmp_1       = ceof(:,0,:)                 
      tmp_2       = ceof_ts(:,0,:)               
      tmp_3       = eof_cdata(0,:)
      tmp_4       = eof_ts_cdata(0,:)
     
      tmp_5       = pcv_eof_u200(0) 
      tmp_6       = pcv_eof_u850(0)
      tmp_7       = pcv_eof_olr(0)
      
      ceof(:,0,:)       = ceof(:,1,::)                
      ceof_ts(:,0,:)    = ceof_ts(:,1,:)              
      eof_cdata(0,:)    = eof_cdata(1,:)
      eof_ts_cdata(0,:) = eof_ts_cdata(1,:)
      pcv_eof_u200(0) = pcv_eof_u200(1)
      pcv_eof_u850(0) = pcv_eof_u850(1)
      pcv_eof_olr(0)  = pcv_eof_olr(1)

      ceof(:,1,:)       = tmp_1
      ceof_ts(:,1,:)    = tmp_2
      eof_cdata(1,:)    = tmp_3
      eof_ts_cdata(1,:) = tmp_4
      pcv_eof_u200(1) = tmp_5
      pcv_eof_u850(1) = tmp_6
      pcv_eof_olr(1)  = tmp_7

      ceof@pcvar(0) = tmp_9
      ceof@pcvar(1) = tmp_8
  end if


;  print("==============")
;  printVarSummary(eof_cdata)
;  printMinMax(eof_cdata, True)

;************************************************
; Compute cross correlation of each variable's EOF time series at zero-lag
;************************************************
  r_olr_u850  = escorc(ceof_ts(0,:,:) , ceof_ts(1,:,:) )  ; (neof)
  r_olr_u200  = escorc(ceof_ts(0,:,:) , ceof_ts(2,:,:) )
  r_u850_u200 = escorc(ceof_ts(1,:,:) , ceof_ts(2,:,:) )

;  print("==============")
;  do n=0,neof-1
;     print("neof="+n \
;          +"  r_olr_u850=" +sprintf("%4.3f",r_olr_u850(n))  \
;          +"  r_olr_u200=" +sprintf("%4.3f",r_olr_u200(n))  \
;          +"  r_u850_u200="+sprintf("%4.3f",r_u850_u200(n)) )
;  end do
; print("==============")

;************************************************
; Compute cross correlation of the multivariate EOF; EOF 1 vs EOF 2
;************************************************

  mxlag     = 25
  rlag_01   = esccr(eof_ts_cdata(0,:),eof_ts_cdata(1,:), mxlag)   ; (N,mxlag+1)
  rlag_10   = esccr(eof_ts_cdata(1,:),eof_ts_cdata(0,:), mxlag)   ; (N,mxlag+1)
  ccr_12    = new ( (/2*mxlag+1/), float)    

  ccr_12(mxlag:)    = rlag_10(0:mxlag)   
  ccr_12(0:mxlag)   = rlag_01(::-1)       ; reverse order
;;print(ccr_12)


;************************************************
; Normalize the multivariate EOF 1&2 component time series
; Compute (PC1^2+PC2^2): values > 1 indicate "strong" periods
;************************************************
; stdev is already (by method) 1 no need to normalized (wenyu)
;  eof_ts_cdata(0,:) = eof_ts_cdata(0,:)/stddev(eof_ts_cdata(0,:))
;  eof_ts_cdata(1,:) = eof_ts_cdata(1,:)/stddev(eof_ts_cdata(1,:))

; print(stddev(eof_ts_cdata(0,:)))
; print(stddev(eof_ts_cdata(1,:)))
  mjo_ts_index      = eof_ts_cdata(0,:)^2 + eof_ts_cdata(1,:)^2 
  mjo_ts_index_smt  = runave(mjo_ts_index, 91, 0) ; 91-day running mean

  nGood   = num(.not.ismissing(mjo_ts_index))     ; # non-missing
  nStrong = num(mjo_ts_index .ge. 1.0)
;  print("nGood="+nGood+"   nStrong="+nStrong+"   nOther="+(nGood-nStrong))

;************************************************
; Write PC results to netCDF for use in another example.
;************************************************
  mjo_ts_index!0    = "time"
  mjo_ts_index&time = time 
  mjo_ts_index@long_name = "MJO PC INDEX" 
  mjo_ts_index@info      = "(PC1^2 + PC2^2)" 

  PC1           = eof_ts_cdata(0,:)
  PC1!0         = "time"
  PC1&time      =  time
  PC1@long_name = "PC1"
  PC1@info      = "PC1/stddev(PC1)"

  PC2           = eof_ts_cdata(1,:)
  PC2!0         = "time"
  PC2&time      =  time
  PC2@long_name = "PC2"
  PC2@info      = "PC2/stddev(PC2)"

  diro = "./"+model+"/"
  filo = "MJO_PC_"+expname+".nc"
  system("/bin/rm -f "+diro+filo)   ; remove any pre-existing file
  ncdf = addfile(diro+filo,"c")     ; open output netCDF file
                                    ; make time an UNLIMITED dimension 
  filedimdef(ncdf,"time",-1,True)   ; recommended  for most applications
                                    ; output variables directly
  ncdf->MJO_INDEX = mjo_ts_index    
  ncdf->PC1       = PC1     
  ncdf->PC2       = PC2     

;************************************************
; Write EOF results to netCDF for use in another example.
;************************************************

  EOF1_olr           = ceof(0,0,:)
  EOF1_olr!0         = "lon"
  EOF1_olr&lon      =  lon
  EOF1_olr@long_name = "1st EOF for OLR"

  EOF2_olr           = ceof(0,1,:)
  EOF2_olr!0         = "lon"
  EOF2_olr&lon      =  lon
  EOF2_olr@long_name = "2nd EOF for OLR"

  EOF1_u850           = ceof(1,0,:)
  EOF1_u850!0         = "lon"
  EOF1_u850&lon      =  lon
  EOF1_u850@long_name = "1st EOF for U850"

  EOF2_u850           = ceof(1,1,:)
  EOF2_u850!0         = "lon"
  EOF2_u850&lon      =  lon
  EOF2_u850@long_name = "2nd EOF for U850"

  EOF1_u250           = ceof(2,0,:)
  EOF1_u250!0         = "lon"
  EOF1_u250&lon      =  lon
  EOF1_u250@long_name = "1st EOF for U250"

  EOF2_u250           = ceof(2,1,:)
  EOF2_u250!0         = "lon"
  EOF2_u250&lon      =  lon
  EOF2_u250@long_name = "2nd EOF for U250"

  diro = "./"+model+"/"
  filo = "MJO_EOF_"+expname+".nc"
  system("/bin/rm -f "+diro+filo)   ; remove any pre-existing file
  ncdf = addfile(diro+filo,"c")     ; open output netCDF file
                                    ; make time an UNLIMITED dimension
  ncdf->EOF1_olr       = EOF1_olr
  ncdf->EOF2_olr       = EOF2_olr
  ncdf->EOF1_u850      = EOF1_u850
  ncdf->EOF2_u850      = EOF2_u850
  ncdf->EOF1_u250      = EOF1_u250
  ncdf->EOF2_u250      = EOF2_u250
 
;------------------------------------------------------------
; PLOTS
;------------------------------------------------------------

  yyyymmdd = cd_calendar(time, -2)
  yrfrac   = yyyymmdd_to_yyyyfrac(yyyymmdd, 0.0)
  delete([/ yrfrac@long_name, lon@long_name /])

  day      = ispan(-mxlag, mxlag, 1)
 ;day@long_name = "lag (day)"

  pltPath = pltDir+pltName

  wks = gsn_open_wks(pltType,pltPath)
  gsn_define_colormap(wks,"default")
  plot = new(3,graphic)                

;************************************************
; Multivariate EOF plots
;************************************************
  rts           = True
  rts@gsnDraw   = False       ; don't draw yet
  rts@gsnFrame  = False       ; don't advance frame yet
  rts@gsnScale  = True        ; force text scaling               

  rts@vpHeightF = 0.40        ; Changes the aspect ratio
  rts@vpWidthF  = 0.85
  rts@vpXF      = 0.10        ; change start locations
  rts@vpYF      = 0.75        ; the plot
  rts@xyLineThicknesses = (/2, 2, 2/)
  rts@xyLineColors      = (/"black","red","blue"/)
  rts@gsnYRefLine       = 0.                  ; reference line   
  rts@trXMaxF           = max(lon)
  rts@trXMinF           = min(lon)

  rts@pmLegendDisplayMode    = "Always"            ; turn on legend
  rts@pmLegendSide           = "Top"               ; Change location of 
  rts@pmLegendParallelPosF   = 1.16                ; move units right
  rts@pmLegendOrthogonalPosF = -0.50               ; move units down
  rts@pmLegendWidthF         = 0.15                ; Change width and
  rts@pmLegendHeightF        = 0.15                ; height of legend.
  rts@lgLabelFontHeightF     = 0.0175


  rtsP                       = True                ; modify the panel plot
;  rtsP@gsnMaximize           = True                ; large format
  rtsP@gsnPanelMainString     = "Multivariate EOF: 15S-15N: "+yrStrt+"-"+yrLast 
  
  do n=0,neof-1
    rts@xyExplicitLegendLabels = (/"OLR: "+sprintf("%4.1f", pcv_eof_u200(n)) +"%" \
                                  ,"U850: "+sprintf("%4.1f", pcv_eof_u850(n))+"%" \
                                  ,"U200: "+sprintf("%4.1f", pcv_eof_olr(n))+"%" /)
    rts@gsnLeftString  = "EOF "+(n+1)
    rts@gsnRightString = sprintf("%3.1f",ceof@pcvar(n))  +"%"
    plot(n) = gsn_csm_xy (wks,lon,ceof(:,n,:),rts)
  end do
  gsn_panel(wks,plot(0:1),(/2,1/),rtsP)     ; now draw as one plot

;-----------------------------------------
; The following doesn't work with some older versions of NCL
; With old versions, the user must delete each individually.
;----------------------------------------- 
  delete([/ rts@xyExplicitLegendLabels, rts@pmLegendDisplayMode \
          , rts@xyLineThicknesses     , rts@gsnLeftString       \
          , rts@gsnRightString        , rts@xyLineColors        \
          , rts@trXMaxF               , rts@trXMinF             /] )

  lag                        = ispan(-mxlag,mxlag,1)
  lag@long_name              = "lag (days)"

  plot(0)                    = gsn_csm_xy (wks, lag ,ccr_12,rts)
  rtsP@gsnPanelMainString    = "Cross Correlation: Multivariate EOF: 15S-15N: " \
                   +  yrStrt+"-"+yrLast 
  rtsP@gsnPaperOrientation   = "portrait"        ; force portrait
  gsn_panel(wks,plot(0),(/1,1/),rtsP)     ; now draw as one plot

;************************************************
; MJO "strong" index 
;************************************************
  rts@gsnYRefLine        = 1.0
  rts@gsnYRefLineColor   = "black"
  rts@xyMonoDashPattern  = True
  rts@xyLineColors       = (/"black", "blue"/)
  rts@xyLineThicknesses  = (/1, 2/)
  rts@pmLegendDisplayMode    = "Always"            ; turn on legend
  rts@pmLegendWidthF         = 0.12                ; Change width and
  rts@pmLegendHeightF        = 0.10                ; height of legend.
  rts@pmLegendParallelPosF   = 0.86                ; move units right
  rts@pmLegendOrthogonalPosF = -0.40               ; move units down
  rts@xyExplicitLegendLabels = (/"daily", "91-day runavg" /)

  mjo_ind_plt = new ( (/2,ntim/), typeof(mjo_ts_index))
  mjo_ind_plt(0,:) = mjo_ts_index
  mjo_ind_plt(1,:) = (/ mjo_ts_index_smt /)
  plot(0) = gsn_csm_xy(wks, yrfrac,mjo_ind_plt,rts)

  rtsP@gsnPanelMainString   = "MJO Index: (PC1^2+ PC2^2) : 15S-15N: "+yrStrt+"-"+yrLast 
  gsn_panel(wks,plot(0),(/1,1/),rtsP)     ; now draw as one plot

 end
