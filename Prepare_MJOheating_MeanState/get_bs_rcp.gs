'reinit'
'open /global/homes/w/wenyuz/ln_solver/bs/ncep/ncep.clim.y58-97.t42.ctl'
'open /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/zg.ctl'
'open /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/hus.ctl'
'open /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/ta.ctl'
'open /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/ua.ctl'
'open /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/va.ctl'
'set fwrite -sq -cl -be /global/homes/w/wenyuz/ln_solver/bs/mme/rcp/bs.t42l20.grd'
'set gxout fwrite'

'set x 1 128'
'set y 1 64'


ii=1
while(ii<=12)
'set t 'ii

k=1
while(k<=17)
'set z 'k
'd z.2'
*recl=recl+1
k=k+1
endwhile

k=1
while(k<=8)
'set z 'k
'd rh'
k=k+1
endwhile

k=1
while(k<=8)
'set z 'k
'd q.3'
k=k+1
endwhile

k=1
while(k<=17)
'set z 'k
'd t.4'
k=k+1
endwhile


k=1
while(k<=17)
'set z 'k
'd u.5'
k=k+1
endwhile

k=1
while(k<=17)
'set z 'k
'd v.6'
k=k+1
endwhile

k=1
while(k<=12)
'set z 'k
'd omg'
k=k+1
endwhile


ii=ii+1
endwhile


'disable fwrite'


