pro lnsp4_5_read_transit

t0=TEN(8,45,0)

path='C:\Users\magne\Desktop\LNSP4_5_2024\'
;file_name=DIALOG_PICKFILE(/READ,FILTER='*.fit',  GET_PATH=path)
file_name=path+'lnsp4_5ch_240201_083502_084944.fit';'lnsp4_5ch_240201_102047_103530.fit';'lnsp4_5ch_240201_090452_091951.fit';'lnsp4_5ch_240201_105028_110501.fit'
FXBOPEN, unit, file_name, 1, header1
;print, header1
FXBREAD, unit,ut_r,17;13;9;5;1
FXBREAD, unit, rcp,18;14;10;6;2
FXBREAD, unit,ut_l,19;15;11;7;3
FXBREAD, unit, lcp,20;16;12;8;4
FXBCLOSE,unit

stokes_i=450.

ut_=0.0006




lcp_=INTERPOL(lcp,ut_l, ut_r)

t1=t0-0.07
t2=t0+0.07
a1=WHERE( ut_r GT(t1-0.005)) AND (ut_r LT(t1+0.005))
a2=WHERE( ut_r GT(t2-0.005)) AND (ut_r LT(t2+0.005))
meanr1=MEAN(rcp(a1))
meanl1=MEAN(lcp_(a1))
meanr2=MEAN(rcp(a2))
meanl2=MEAN(lcp_(a2))
rcp0=INTERPOL([meanr1,meanr2],[t1,t2],ut_r)
lcp0=INTERPOL([meanl1,meanl2],[t1,t2],ut_r)

rcp=rcp-rcp0
lcp=lcp_-lcp0

x=-15*60.*(ut_r-t0)

;print, x

a=SORT(x)
x=x(a)

;print,x

rcp=rcp(a)
lcp=lcp(a)

window,13, xsize=600, ysize=800
device, decomposed=0
loadct,0

plot, x, (rcp+lcp)/stokes_i, xrange=[-40.,40.], xstyle=1, thick=3, yrange=[0., 2.5], ystyle=1,background=255, color=0;,
;oplot, ut_l, lcp



end
