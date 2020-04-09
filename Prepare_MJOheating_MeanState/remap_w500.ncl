begin

newlat = (/-87.8638,-85.0965,-82.3129, -79.5256, -76.7369, -73.9475, -71.1577, -68.3678, -65.5776, -62.7873, -59.9970, -57.2066, -54.4162, -51.6257, -48.8352,-46.0447, -43.2542, -40.4636, -37.6731, -34.8825, -32.0919, -29.3014, -26.510, -23.7202, -20.9296, -18.1390, -15.3484 ,-12.5578, -9.76715, -6.97653, -4.18592, -1.39531, 1.39531, 4.18592, 6.97653, 9.76715, 12.5578, 15.3484, 18.1390 ,20.9296, 23.7202, 26.5108, 29.3014, 32.0919, 34.8825, 37.6731, 40.4636, 43.2542, 46.0447, 48.8352, 51.6257, 54.4162, 57.2066, 59.9970, 62.7873, 65.5776, 68.3678, 71.1577, 73.9475, 76.7369, 79.5256, 82.3129, 85.0965, 87.8638/)
newlon = fspan(0,357.1875,128)
newlev = (/.99500,0.97999,0.94995,0.89988,0.82977,0.74468,0.64954,0.54946,0.45447,0.36948,0.29450,0.22953,0.17457,0.12440,0.0846830,0.0598005,0.0449337,0.0349146,0.0248800,0.00829901/)
newlat@units = "degrees_north"
newlon@units = "degrees_east"
time_org = (/3.5,9.5,15.5,21.5,27.5,33.5,39.5,45.5/)
newtime = fspan(1,48,48)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/w500/mme_hist.nc","r")
sst=f1->w500(1,:,:,:)
lon=f1->lon
lat=f1->lat
time=f1->phase
newsst1 =new((/8,64,128/),typeof(sst))
newsst2 =new((/48,64,128/),typeof(sst))
newsst =new((/48,20,64,128/),typeof(sst))
newsst =0.
newsst1 = linint2(lon,lat,sst,True,newlon,newlat,0)
printMinMax(newsst1,1)
newsst2 = linint1_n(time_org,newsst1,True,newtime,0,0)
printMinMax(newsst2,1)
i = 0
do while(i.le.19)
;   newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*exp(-20.*(newlev(i)-0.45 )^2)*864./20.*0.00002
;    newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*cos((newlev(i)-0.5)*get_pi("f"))*864*0.056/86400
    newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(1.-2.*abs(newlev(i)-0.5))*864*0.06/86400
;    if (newlev(i).gt.0.5) then
;      newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(1.-(newlev(i)-0.5)/0.5)*864*0.056/86400
;    else
;      newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(newlev(i)-0.1)/0.4*864*0.056/86400
;    end if
;   convert hpa/day to K/day * 0.56
; original unit pa/s, final unit K/s
   i=i+1
end do
printMinMax(newsst,1)
newsst3 = newsst(:,:,::-1,:)
newsst!0   ="time"
newsst!1   ="lev"
newsst!2   ="lat"
newsst!3   ="lon"
newsst&time = newtime
newsst&lev = newlev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/w500/mme_hist.winter.grd",newsst3)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/w500/mme_rcp.nc","r")
sst=f1->w500(1,:,:,:)
lon=f1->lon
lat=f1->lat
time=f1->phase
newsst1 =new((/8,64,128/),typeof(sst))
newsst2 =new((/48,64,128/),typeof(sst))
newsst =new((/48,20,64,128/),typeof(sst))
newsst =0.
newsst1 = linint2(lon,lat,sst,True,newlon,newlat,0)
printMinMax(newsst1,1)
newsst2 = linint1_n(time_org,newsst1,True,newtime,0,0)
printMinMax(newsst2,1)
i = 0
do while(i.le.19)
    newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(1.-2.*abs(newlev(i)-0.5))*864*0.06/86400*1.25
;    newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(1.-max(/1., 2.5*abs(newlev(i)-0.5)/))*864*0.056/86400*1.25
;    if (newlev(i).gt.0.5) then
;      newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(1.-(newlev(i)-0.48)/0.52)*864*0.056/86400*1.25
;    else
;      newsst(:,i,24:35,15:77)=-newsst2(:,24:35,15:77)*(newlev(i)-0.08)/0.4*864*0.056/86400*1.25
;    end if
   i=i+1
end do
printMinMax(newsst,1)
newsst3 = newsst(:,:,::-1,:)
newsst!0   ="time"
newsst!1   ="lev"
newsst!2   ="lat"
newsst!3   ="lon"
newsst&time = newtime
newsst&lev = newlev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/w500/mme_rcp.winter.grd",newsst3)

end
