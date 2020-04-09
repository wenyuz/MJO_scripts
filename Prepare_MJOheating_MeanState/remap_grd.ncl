begin

newlat = (/-87.8638,-85.0965,-82.3129, -79.5256, -76.7369, -73.9475, -71.1577, -68.3678, -65.5776, -62.7873, -59.9970, -57.2066, -54.4162, -51.6257, -48.8352,-46.0447, -43.2542, -40.4636, -37.6731, -34.8825, -32.0919, -29.3014, -26.510, -23.7202, -20.9296, -18.1390, -15.3484 ,-12.5578, -9.76715, -6.97653, -4.18592, -1.39531, 1.39531, 4.18592, 6.97653, 9.76715, 12.5578, 15.3484, 18.1390 ,20.9296, 23.7202, 26.5108, 29.3014, 32.0919, 34.8825, 37.6731, 40.4636, 43.2542, 46.0447, 48.8352, 51.6257, 54.4162, 57.2066, 59.9970, 62.7873, 65.5776, 68.3678, 71.1577, 73.9475, 76.7369, 79.5256, 82.3129, 85.0965, 87.8638/)
newlon = fspan(0,357.1875,128)
newlat@units = "degrees_north"
newlon@units = "degrees_east"

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/ps_mme.mon.regrid.nc","r")
sst=f1->ps/100.
lon=f1->lon
lat=f1->lat
time=f1->time
newsst =new((/12,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="lat"
newsst!2   ="lon"
newsst&time = time
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/ps.t42l20.grd",newsst)
delete(sst)
delete(newsst)


f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/hus_mme.mon.regrid.nc","r")
sst=f1->hus(:,0:7,:,:)
lon=f1->lon
lat=f1->lat
time=f1->time
plev=f1->plev(0:7)/100.
newsst =new((/12,8,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="plev"
newsst!2   ="lat"
newsst!3   = "lon"
newsst&time = time
newsst&plev = plev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/husnew.grd",newsst)
delete(sst)
delete(newsst)
delete(plev)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/ua_mme.mon.regrid.nc","r")
sst=f1->ua(:,0:16,:,:)
lon=f1->lon
lat=f1->lat
time=f1->time
plev=f1->plev(0:16)/100.
newsst =new((/12,17,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="plev"
newsst!2   ="lat"
newsst!3   = "lon"
newsst&time = time
newsst&plev = plev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/uanew.grd",newsst)
delete(sst)
delete(newsst)
delete(plev)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/va_mme.mon.regrid.nc","r")
sst=f1->va(:,0:16,:,:)
lon=f1->lon
lat=f1->lat
time=f1->time
plev=f1->plev(0:16)/100.
newsst =new((/12,17,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="plev"
newsst!2   ="lat"
newsst!3   = "lon"
newsst&time = time
newsst&plev = plev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/vanew.grd",newsst)
delete(sst)
delete(newsst)
delete(plev)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/ta_mme.mon.regrid.nc","r")
sst=f1->ta(:,0:16,:,:)
lon=f1->lon
lat=f1->lat
time=f1->time
plev=f1->plev(0:16)/100.
newsst =new((/12,17,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="plev"
newsst!2   ="lat"
newsst!3   = "lon"
newsst&time = time
newsst&plev = plev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/tanew.grd",newsst)
delete(sst)
delete(newsst)
delete(plev)

f1 = addfile("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/zg_mme.mon.regrid.nc","r")
sst=f1->zg(:,0:16,:,:)
lon=f1->lon
lat=f1->lat
time=f1->time
plev=f1->plev(0:16)/100.
newsst =new((/12,17,64,128/),typeof(sst))
newsst = linint2(lon,lat,sst,True,newlon,newlat,0)
newsst!0   ="time"
newsst!1   ="plev"
newsst!2   ="lat"
newsst!3   = "lon"
newsst&time = time
newsst&plev = plev
newsst&lat = newlat
newsst&lon = newlon
fbinwrite("/global/homes/w/wenyuz/ln_solver/bs/mme/hist/zgnew.grd",newsst)
delete(sst)
delete(newsst)
delete(plev)

end
