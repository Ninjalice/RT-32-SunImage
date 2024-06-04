pro data_to_sun_spiral5
lat_rt32=57.5535171694D0     ;RT32 geographic latitude 57d 33' 12''.661
lon_rt32=21.8545525000D0     ;RT32 geographic longitude 21d 51' 16''.389
;postroenie 2D kart dlya dannih iz SP3 ( stary priemnik) iz spiralnogo scanirovaniya
PRINT, " "
PRINT, "---------------------------------------------------------------------------"
PRINT, "---------                  SINTEZ 2D KART STOKES I V            -----------"
PRINT, " "
;chitatj fail programmi scanovirovaniya *.ptf
az_offset0=0.
el_offset0=0. ; offset nuley datchikov polozheniya antenni
;
fname_antenna=DIALOG_PICKFILE(/READ,FILTER='*.ptf', TITLE='Select scan program file sYYMMDDN.ptf', GET_PATH=path)
str='  '
OPENR, file_ant, fname_antenna,/GET_LUN
;chitatj stroki, poka ne poyavitsya [Table Dat a]
REPEAT BEGIN
	READF, file_ant, str
	IF STRCMP(str,'# Az offset',11) THEN BEGIN
		str1=STRSPLIT(STRCOMPRESS(STRTRIM(str,2)),' ',/EXTRACT)
		az_offset=FLOAT(str1[4])
		PRINT, "Az offset=", az_offset0, "  deg"
	ENDIF
	IF STRCMP(str,'# El offset',11) THEN BEGIN
		str1=STRSPLIT(STRCOMPRESS(STRTRIM(str,2)),' ',/EXTRACT)
		el_offset=FLOAT(str1[4])
		PRINT, "El offset=", el_offset0, "  deg"
	ENDIF
ENDREP UNTIL STRCMP(str,'[Table Data]',12)
READF, file_ant, str ; propustitj esche ognu stroku
	i=0
	WHILE NOT EOF(file_ant) DO BEGIN
		READF, file_ant, str
		str1=STRSPLIT(STRCOMPRESS(STRTRIM(str,2)),' ',/EXTRACT)
		str2=STRSPLIT(STRCOMPRESS(STRTRIM(str1[0],2)),'T',/EXTRACT)
		str3=STRSPLIT(STRCOMPRESS(STRTRIM(str2[0],2)),'-',/EXTRACT)
		str4=STRSPLIT(STRCOMPRESS(STRTRIM(str2[1],2)),':',/EXTRACT)
		IF i EQ 0 THEN BEGIN
			year=FLOAT(str3[0])
			month=FLOAT(str3[1])
			day=FLOAT(str3[2])
			print,'  '
			print, "Nachalo sessii "+str1[0]
		ENDIF
		utime=TEN(FLOAT(str4[0]),FLOAT(str4[1]), FLOAT(str4[2]))
		IF i EQ 0 THEN BEGIN
    		utime_ant=utime
    		az_ant=FLOAT(str1[1])
    		el_ant=FLOAT(str1[2])
    		i=i+1
  		ENDIF ELSE BEGIN
  			utime_ant=[utime_ant,utime]
    		az_ant=[az_ant,FLOAT(str1[1])]
    		el_ant=[el_ant,FLOAT(str1[2])]
    		i=i+1
		ENDELSE
	ENDWHILE
CLOSE, /ALL
print, "Konec sessii    "+str1[0]
PRINT, "Prochitano", N_ELEMENTS(utime_ant),' tochek polozheniy antenni'
utime_start=utime_ant[0]
utime_end=utime_ant[N_ELEMENTS(utime_ant)-1]
year=FLOAT(str3[0])
month=FLOAT(str3[1])
day=FLOAT(str3[2])
RT32_SUN_POSITION, year, month, day, utime_ant, az_sun, el_sun, /REFRA
q=RT32_SUN_PARA(year, month, day, utime_ant)

az_offset0=+0.02

az_ant=az_ant-az_offset0

;----------------------
el_offset0=+0.16
;--------------------

el_ant=el_ant-el_offset0
x0_ant=(az_ant-az_sun)*COS(!DTOR*el_ant)
y0_ant=(el_ant-el_sun)
;povernuty na q,esli +q -povertutj po chas strelke!
;koordinati antenni ot centa Solnca v kartinnoy ploskosti
x1_ant= 3600.*(x0_ant*COS(!DTOR*q) +y0_ant*SIN(!DTOR*q))
y1_ant= 3600.*(-x0_ant*SIN(!DTOR*q) +y0_ant*COS(!DTOR*q))
;
;
num_ch=3
n_sample=LNSP5CH_file_READ(num_ch, year, month, day, utime_start, utime_end, ut_r, ut_l, rcp, lcp, PATH=path)
print, n_sample
DEVICE, decomposed=0
WINDOW, 13
loadct, 3
plot, ut_r,rcp, YRANGE=[MIN(rcp), MAX(rcp)]
oplot, ut_l,lcp, COLOR=180
;

;otnositelnaya calibrovka
sky0_r=DBLARR(5)
sky0_l=DBLARR(5)
sun0_r=DBLARR(5)
sun0_l=DBLARR(5)
t_sun0=DBLARR(5)
t_sky0=DBLARR(5)
FOR i=0,4 DO BEGIN
	t_sun=utime_start +30./3600. + i*660./3600.
	t_sky=utime_start +(660.-30.-20.)/3600. +i*660./3600.
	t_sun0[i]=t_sun
	t_sky0[i]=t_sky
	a=WHERE ((ut_r GT t_sun-18./3600.) AND (ut_r LT t_sun+18./3600.))
	sun0_r[i]=MEAN(rcp(a))
	a=WHERE ((ut_l GT t_sun-18./3600.) AND (ut_l LT t_sun+18./3600.))
	sun0_l[i]=MEAN(lcp(a))
	b=WHERE ((ut_r GT t_sky-18./3600.) AND (ut_r LT t_sky+18./3600.))
	sky0_r[i]=MEAN(rcp(b))
	b=WHERE ((ut_l GT t_sky-18./3600.) AND (ut_l LT t_sky+18./3600.))
	sky0_l[i]=MEAN(lcp(b))
ENDFOR
rcp_=(rcp-INTERPOL(sky0_r,t_sky0,ut_r))/(INTERPOL(sun0_r,t_sun0,ut_r)-INTERPOL(sky0_r,t_sky0,ut_r))
lcp_=(lcp-INTERPOL(sky0_l,t_sky0,ut_l))/(INTERPOL(sun0_l,t_sun0,ut_l)-INTERPOL(sky0_l,t_sky0,ut_l))
;
;opredelim XY semplov
x_r=INTERPOL(x1_ant,utime_ant,ut_r)
x_l=INTERPOL(x1_ant,utime_ant,ut_l)
y_r=INTERPOL(y1_ant,utime_ant,ut_r)
y_l=INTERPOL(y1_ant,utime_ant,ut_l)

;RCP
;virezaem scan 1
ind=WHERE( (ut_r GT utime_start+60./3600.) AND (ut_r LT utime_start +(660.-20.-60.-20.)/3600.) )
xx_r=x_r(ind)
yy_r=y_r(ind)
rrcp_=rcp_(ind)
;virezaem scan 2
ind=WHERE( (ut_r GT (utime_start+660./3600.)+60./3600.) AND (ut_r LT (utime_start+660./3600.) +(660.-20.-60.-20.)/3600.) )
xx_r=[xx_r,x_r(ind)]
yy_r=[yy_r,y_r(ind)]
rrcp_=[rrcp_,rcp_(ind)]
ind=WHERE( (ut_r GT (utime_start+2.*660./3600.)+60./3600.) AND (ut_r LT (utime_start+2.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_r=[xx_r,x_r(ind)]
yy_r=[yy_r,y_r(ind)]
rrcp_=[rrcp_,rcp_(ind)]
ind=WHERE( (ut_r GT (utime_start+3.*660./3600.)+60./3600.) AND (ut_r LT (utime_start+3.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_r=[xx_r,x_r(ind)]
yy_r=[yy_r,y_r(ind)]
rrcp_=[rrcp_,rcp_(ind)]
ind=WHERE( (ut_r GT (utime_start+4.*660./3600.)+60./3600.) AND (ut_r LT (utime_start+4.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_r=[xx_r,x_r(ind)]
yy_r=[yy_r,y_r(ind)]
rrcp_=[rrcp_,rcp_(ind)]
dd=10. ; pixel 10x10 arcsec
rcp_map=TRI_SURF(rrcp_,xx_r/dd,yy_r/dd, GS=[1,1], BOUNDS=[-240.,-240.,240.,240.], /LINEAR);, MISSING=-1.)
;
;LCP
;RCP
;virezaem scan 1
ind=WHERE( (ut_l GT utime_start+60./3600.) AND (ut_l LT utime_start +(660.-20.-60.-20.)/3600.) )
xx_l=x_l(ind)
yy_l=y_l(ind)
llcp_=lcp_(ind)
;virezaem scan 2
ind=WHERE( (ut_l GT (utime_start+660./3600.)+60./3600.) AND (ut_l LT (utime_start+660./3600.) +(660.-20.-60.-20.)/3600.) )
xx_l=[xx_l,x_l(ind)]
yy_l=[yy_l,y_l(ind)]
llcp_=[llcp_,lcp_(ind)]
ind=WHERE( (ut_l GT (utime_start+2.*660./3600.)+60./3600.) AND (ut_l LT (utime_start+2.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_l=[xx_l,x_l(ind)]
yy_l=[yy_l,y_l(ind)]
llcp_=[llcp_,lcp_(ind)]
ind=WHERE( (ut_l GT (utime_start+3.*660./3600.)+60./3600.) AND (ut_l LT (utime_start+3.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_l=[xx_l,x_l(ind)]
yy_l=[yy_l,y_l(ind)]
llcp_=[llcp_,lcp_(ind)]
ind=WHERE( (ut_l GT (utime_start+4.*660./3600.)+60./3600.) AND (ut_l LT (utime_start+4.*660./3600.)+(660.-20.-60.-20.)/3600.) )
xx_l=[xx_l,x_l(ind)]
yy_l=[yy_l,y_l(ind)]
llcp_=[llcp_,lcp_(ind)]
dd=10. ; pixel 10x10 arcsec
lcp_map=TRI_SURF(llcp_,xx_l/dd,yy_l/dd, GS=[1,1], BOUNDS=[-240.,-240.,240.,240.], /LINEAR);, MISSING=-1.)

s=SIZE(rcp_map)
nx=s[1]
ny=s[2]
x_map=10.*(FINDGEN(nx)-FLOOR(nx/2.))
y_map=10.*(FINDGEN(ny)-FLOOR(ny/2.))
datetime=RT32_DATETIME_TO_STR(year, month, day, utime_start)
print, datetime
RT32_SUN_PBR, year, month, day, p, b0, rsun

st_i=(rcp_map+lcp_map)/2.
st_v=(rcp_map-lcp_map)/2.
stokes_i=make_map(st_i, XC=0., YC=dd, DX=dd, DY=dd, date=datetime, ROLL=-p)
stokes_v=make_map(st_v, XC=0., YC=dd, DX=dd, DY=dd, date=datetime, ROLL=-p)

r_map=make_map(rcp_map, XC=0., YC=dd, DX=dd, DY=dd, date=datetime, ROLL=-p)
l_map=make_map(lcp_map, XC=0., YC=dd, DX=dd, DY=dd, date=datetime, ROLL=-p)

WINDOW,16, XSIZE=800, YSIZE=800
loadct,3
plot_map, stokes_i, /square, grid=10., FOV=[40.,40.], TITLE=datetime, /COLORBAR, DMIN=0.10, DMAX=3.5, CHARSIZE=1.7, CHARTHICK=1.2
WINDOW,17, XSIZE=800, YSIZE=800
loadct,13
plot_map, stokes_v, /square, grid=10., FOV=[40.,40.], TITLE=datetime, /COLORBAR, DMIN=-0.10, DMAX=+0.1, CHARSIZE=1.7, CHARTHICK=1.2

end