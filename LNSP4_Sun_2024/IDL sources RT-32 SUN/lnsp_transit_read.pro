pro lnsp_transit_read

file_name=DIALOG_PICKFILE(/READ,FILTER='*.fit',  GET_PATH=path)


FXBOPEN, unit, file_name, 1, header1
print, header1
FXBREAD, unit,ut,1
FXBREAD, unit,a,2
FXBREAD, unit,pol,3
FXBCLOSE,unit

ut_r=0
ut_l=0
a_r=0.
a_l=0.
n=N_ELEMENTS(ut)
FOR i=0, n-1 DO BEGIN
	IF pol[i] EQ 0 THEN BEGIN
		ut_r=[ut_r, ut[i]]
		a_r=[a_r, a[i]]
	ENDIF ELSE BEGIN
		ut_l=[ut_r, ut[i]]
		a_l=[a_l, a[i]]
	ENDELSE
ENDFOR
ut_r=ut_r[3:N_ELEMENTS(ut_r)-10]
ut_l=ut_l[3:N_ELEMENTS(ut_l)-10]
rcp=a_r[3:N_ELEMENTS(a_r)-10]
lcp=a_l[3:N_ELEMENTS(a_l)-10]

rcp0=MEAN(rcp[20:40])
rcpmax=MAX(rcp-rcp0)
lcp0=MEAN(lcp[20:40])
lcpmax=MAX(lcp-lcp0)
window,13, xsize=800,ysize=600
device, decomposed=0
loadct,3
rcpmax=20

plot, ut_r, smooth(rcp-rcp0,2)+smooth(lcp-lcp0, 2), yrange=[-0.2*rcpmax, 2*rcpmax], COLOR=255
oplot, ut_l, smooth(lcp-lcp0, 2),COLOR=120


end