'reinit'
'open /global/homes/w/wenyuz/ln_solver/bs/w500/mme_hist.winter.ctl'
'set fwrite -sq -cl -be /global/homes/w/wenyuz/ln_solver/sample/frc.t42l20.mjo.hist.winter.grd'
'set gxout fwrite'

'set x 1 128'
'set y 1 64'


ii=1
'set z '1
'set t '1
'define aa = tmp*0.'
while(ii<=48)
'set t 'ii

k=1
while(k<=20)
'set z 'k
'd aa'
k=k+1
endwhile

k=1
while(k<=20)
'set z 'k
'd aa'
k=k+1
endwhile


k=1
while(k<=20)
'set z 'k
'd tmp'
k=k+1
endwhile

k=1
'set z 'k
'd aa'


ii=ii+1
endwhile


'disable fwrite'


