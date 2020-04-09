;-------------------------------------------------------------
; mjoclivar_2.ncl
;-------------------------------------------------------------
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
;******************************************************

begin
;******************************************************
; User specifications
;******************************************************
   NC      = True                            ; create netCDF?       
   PLOT    = False                             ; sample plots?

  if (expname.eq."historical") then
     expdir = "hist"
  else
     expdir = "rcp85"
  end if
   yrStrt  = ymdStrt/10000
   yrLast  = ymdLast/10000
  
   nhar    = 4                                ; number of fourier comp      

   diri    = "./"+model+"/"+expdir+"/"

vNames = (/"rlut","u850","u250","w500","pr","zg250"/)
nvar = dimsizes(vNames)
do i =0, nvar-1
   vName = vNames(i)
   print(vName)
   if (vName.eq."pr") then 
      var     = "pr"                            ; name of file                   
      fili    = var+"_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."rlut") then
      var     = "rlut"                            ; name of file
      fili    = var+"_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."tas") then
      var     = "tas"                            ; name of file
      fili    = var+"_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."v850") then
      var     = "va"
      fili    = var+"850_day_"+model+"_"+expname+".regrid.nc"               ; input file       
   end if
   if (vName.eq."w500") then
      var     = "wap"
      fili    = var+"500_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."u850") then
      var     = "ua"
      fili    = var+"850_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."u250") then
      var     = "ua"
      fili    = var+"250_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if
   if (vName.eq."zg250") then
      var     = "zg"
      fili    = var+"250_day_"+model+"_"+expname+".regrid.nc"               ; input file
   end if

   if (NC) then
       diro= "./"+model+"/anom/"         ; output dir
       filo= vName+".day.anomalies."+yrStrt+"-"+yrLast+".nc"  ; output file  
   end if

   if (PLOT) then
       wksType = "png"                        ; send graphics to PNG file
       wksName = "mjoclivar"                  ; "mjo."+yrStrt+"_"+yrLast
   end if

;***********************************************************
; Read user specified time and create required yyyyddd                    
;***********************************************************
   f       = addfile (diri+fili , "r")                          

   time    = f->time                          ; all times on file
   ymd     = cd_calendar(time, -2)            ; yyyymmdd
   iStrt   = ind(ymd.eq.ymdStrt)              ; index start
   iLast   = ind(ymd.eq.ymdLast)              ; index last 
   delete(time)
   delete(ymd)

;***********************************************************
; Read user specified time and create required yyyyddd                    
;***********************************************************
   time    = f->time(iStrt:iLast)             ; time:units = "hours since"
   TIME    = cd_calendar(time, 0)             ; type float 
   year    = floattointeger( TIME(:,0) )
   month   = floattointeger( TIME(:,1) )
   day     = floattointeger( TIME(:,2) ) 
   ddd     = day_of_year(year, month, day) 
   yyyyddd = year*1000 + ddd                  ; needed for input

;***********************************************************
; Read data: short2flt                                     
;*********************************************************** 
   if ((var.eq."pr") .or. (var.eq."rlut")) then
     x       =  short2flt( f->$var$(iStrt:iLast,:,:) )    ; convert to float
   else
     x       =  short2flt( f->$var$(iStrt:iLast,0,:,:) )    ; convert to float 
   end if
;   printVarSummary( x ) 

;***********************************************************
; Compute daily climatology: raw and then 'smoothed'  
;***********************************************************
   xClmDay = clmDayTLL(x, yyyyddd)     ; daily climatology at each grid point
;   printVarSummary(xClmDay)            

;***********************************************************
; Compute smoothed daily climatology using 'nhar' harmonics
;***********************************************************
   xClmDay_sm = smthClmDayTLL(xClmDay, nhar)
;   printVarSummary(xClmDay_sm)

;***********************************************************
; Compute daily anomalies using raw and smoothed climatologies
;***********************************************************
    xAnom      = calcDayAnomTLL (x, yyyyddd, xClmDay)     
;    printVarSummary(xAnom)
;    printMinMax(xAnom, True)

    xAnom_sm   = calcDayAnomTLL (x, yyyyddd, xClmDay_sm)     
    xAnom_sm@long_name = "Anomalies from Smooth Daily Climatology"
;    printVarSummary(xAnom_sm)
;    printMinMax(xAnom_sm, True)

    delete( x )    ; no longer needed


;***********************************************************
; Create netCDF: convenience use 'simple' method
;***********************************************************

    dimx   = dimsizes(xAnom)
    ntim   = dimx(0)
    nlat   = dimx(1)
    mlon   = dimx(2)

    if (NC) then

        lat    = f->lat
        lon    = f->lon

        system("/bin/rm -f "+diro+filo)      ; rm any pre-exist file, if any
        fnc    = addfile (diro+filo, "c")

        filAtt = 0
        filAtt@title         = vName+": Daily Anomalies: "+yrStrt+"-"+yrLast  
	filAtt@source_file   = fili
        filAtt@creation_date = systemfunc("date")
        fileattdef( fnc, filAtt )         ; copy file attributes  

        setfileoption(fnc,"DefineMode",True)
        
        varNC      = vName+"_anom"
        varNC_sm   = vName+"_anom_sm"

        dimNames = (/"time", "lat", "lon"/)  
	dimSizes = (/ -1   ,  nlat,  mlon/) 
	dimUnlim = (/ True , False, False/)   
	filedimdef(fnc,dimNames,dimSizes,dimUnlim)

        filevardef(fnc, "time"  ,typeof(time),getvardims(time)) 
        filevardef(fnc, "lat"   ,typeof(lat) ,getvardims(lat)) 
        filevardef(fnc, "lon"   ,typeof(lon) ,getvardims(lon))
        filevardef(fnc, varNC   ,typeof(xAnom)   ,getvardims(xAnom))    
        filevardef(fnc, varNC_sm,typeof(xAnom)   ,getvardims(xAnom))    

        filevarattdef(fnc,"time"  ,time)                    ; copy time attributes
        filevarattdef(fnc,"lat"   ,lat)                     ; copy lat attributes
        filevarattdef(fnc,"lon"   ,lon)                     ; copy lon attributes
        filevarattdef(fnc,varNC   ,xAnom)                
        filevarattdef(fnc,varNC_sm,xAnom_sm)                

        fnc->time       = (/time/)     
        fnc->lat        = (/lat/)
        fnc->lon        = (/lon/)
        fnc->$varNC$    = (/xAnom   /)
        fnc->$varNC_sm$ = (/xAnom_sm/)
    end if
end do
;************************************************
; plotting parameters
;************************************************
    if (PLOT) then
        LAT   = (/ 60, 45,  5,  -5, -45, 60 /)
        LON   = (/270, 30, 90,  90, 180,  0 /)
        nPts  = dimsizes( LAT )

        plot  = new ( nPts, graphic)
        data  = new ( (/2,366/), typeof(xClmDay), getFillValue(xClmDay))

        wks   = gsn_open_wks (wksType,wksName)    

        res                   = True                      ; plot mods desired
        res@gsnDraw           = False
        res@gsnFrame          = False
        res@trXMinF           =   1
        res@trXMaxF           = 366
       ;res@tiMainString      = ""
       
        res@xyLineThicknesses = (/1.0, 2.0/)              ; make 2nd lines thicker
        res@xyLineColors      = (/"blue","red"/)          ; change line color
        res@xyMonoDashPattern = True                      ; all solid

        do np=0,nPts-1 
           data(0,:) = xClmDay(:,{LAT(np)},{LON(np)})
           data(1,:) = xClmDay_sm(:,{LAT(np)},{LON(np)})
           res@gsnCenterString = "lat="+LAT(np)+"  lon="+LON(np)
           plot(np)  = gsn_csm_y (wks,data,res) ; create plot
        end do
       
        resP                     = True                ; modify the panel plot
        resP@gsnPanelMainString  = vName+": Raw/Smooth Daily Climatology: "+yrStrt+"-"+yrLast
        resP@gsnMaximize         = True
        resP@gsnPaperOrientation = "portrait"
        gsn_panel(wks,plot,(/(nPts/2),2/),resP)               ; now draw as one plot

       ;==========
       ; Plot anomalies for an arbitrarily selected near equatorial location
       ; Time: Oct 1, 1996 to April 1,1997  [arbitrary selection]
       ;==========
        LATX    = 0
        LONX    = 90

        yyyymmdd  = cd_calendar(time, -2)
      ;;yrfrac    = yyyymmdd_to_yyyyfrac (yyyymmdd, 0)
      ;;delete(yrfrac@long_name)

        xAnom@long_name    = "Anomalies from Raw"   ; short labels for plot
        xAnom_sm@long_name = "Anomalies from Smooth"

        ntBegin   = ind(yyyymmdd.eq.19961001)
        ntEnd     = ind(yyyymmdd.eq.19970401)
       
        monthLabels    = (/1,4,7,10/)
        monNam = (/"Jan","Feb","Mar","Apr","May","Jun" \
                  ,"Jul","Aug","Sep","Oct","Nov","Dec" /)
        xVal   = new(ntim, typeof(xAnom&time) , "No_FillValue") ; bigger than
        xLab   = new(ntim, "string", "No_FillValue")        ; needed
        xValm  = new(ntim, typeof(xAnom&time) , "No_FillValue") ; bigger than

        ntm            = -1
        cr             = inttochar(10)                     ; carriage return
        do nt=ntBegin,ntEnd
         if (day(nt).eq.1) then
             ntm       = ntm + 1
             xVal(ntm) = xAnom&time(nt)
             xLab(ntm) = monNam(month(nt)-1)
             if (month(nt).eq.1) then
                 xLab(ntm) = xLab(ntm) + cr +sprinti("%0.4i", year(nt))
             end if
         end if
        end do

        rxy  = True      
        rxy@gsnDraw     = False
        rxy@gsnFrame    = False
        rxy@gsnYRefLine = 0.0                    ; create a reference line   
        rxy@gsnAboveYRefLineColor = "red"        ; above ref line fill red
        rxy@gsnBelowYRefLineColor = "blue"       ; below ref line fill blue
        rxy@xyLineThicknessF  = 2.0                                               
        rxy@vpHeightF  = 0.4                     ; resize            
        rxy@vpWidthF   = 0.8                  
        rxy@tmXBMode   = "Explicit"
        rxy@tmXBValues = xVal(0:ntm)
        rxy@tmXBLabels = xLab(0:ntm)

        plot(0)                  = gsn_csm_xy (wks,time(ntBegin:ntEnd) \
                              ,xAnom(ntBegin:ntEnd,{0},{90}),rxy) 
        plot(1)                  = gsn_csm_xy (wks,time(ntBegin:ntEnd) \
                              ,xAnom_sm(ntBegin:ntEnd,{0},{90}),rxy) 
        resP@gsnPanelMainString  = vName+": Anomalies: (0,90E)"
        gsn_panel(wks,plot(0:1),(/2,1/),resP)   

    end if
end
