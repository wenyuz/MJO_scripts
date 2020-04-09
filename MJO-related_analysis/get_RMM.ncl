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
   pltName = "mjoclivar14_"+expname+"_baseonREOF"   

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
; Read in reanalysis EOFs
;************************************************

  eof_cdata = new( (/neof,3*mlon/), typeof(cdata), getFillValue(cdata))

   f       = addfile ("/Volumes/Zhou/daily/NCEP/MJO_EOF_historical.nc", "r")
   EOF1_olr    = f->EOF1_olr
   EOF2_olr    = f->EOF2_olr
   EOF1_u850   = f->EOF1_u850
   EOF2_u850   = f->EOF2_u850
   EOF1_u250   = f->EOF1_u250
   EOF2_u250   = f->EOF2_u250

   eof_cdata(0,0:mlon-1) = EOF1_olr
   eof_cdata(0,mlon:2*mlon-1) = EOF1_u850
   eof_cdata(0,2*mlon:) = EOF1_u250

   eof_cdata(1,0:mlon-1) = EOF2_olr
   eof_cdata(1,mlon:2*mlon-1) = EOF2_u850
   eof_cdata(1,2*mlon:) = EOF2_u250

;  print("==============")
;  printVarSummary(eof_cdata)
;  printMinMax(eof_cdata, True)

;************************************************
; compute PC1 and PC2
;************************************************

  eof_ts_cdata = eofunc_ts(cdata,eof_cdata,False)   ; (neof,3*ntim)
;  print("==============")
;  printVarSummary(eof_ts_cdata)
;  printMinMax(eof_ts_cdata, True)

;************************************************
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
  filo = "MJO_PC_"+expname+"_baseonREOF.nc"
  system("/bin/rm -f "+diro+filo)   ; remove any pre-existing file
  ncdf = addfile(diro+filo,"c")     ; open output netCDF file
                                    ; make time an UNLIMITED dimension 
  filedimdef(ncdf,"time",-1,True)   ; recommended  for most applications
                                    ; output variables directly
  ncdf->MJO_INDEX = mjo_ts_index    
  ncdf->PC1       = PC1     
  ncdf->PC2       = PC2     

 end
