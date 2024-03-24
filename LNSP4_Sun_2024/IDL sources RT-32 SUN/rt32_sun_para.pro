function RT32_SUN_PARA, year, month, day, utc
;vozvrascheat parallaksniy ugol Solnca
;
;INPUT:
;year, month, day - data
;utc - UT1 , hours  scalar or vector
;
;OUTPUT:
;paralaktichesky ugol v gradusah, scalar of vector of utc size
;
;ispozuet
;RT32_SUN_TO_HOR
;RT32_LST_TRUE
;
lat_rt32=57.5535171694D0     ;RT32 geographic latitude 57d 33' 12''.661
lon_rt32=21.8545525000D0     ;RT32 geographic longitude 21d 51' 16''.389
a=size(utc)                ;create array if utc is scalar
IF a[0] EQ 0 THEN BEGIN
    utime=DBLARR(1)
    utime[0]=utc
ENDIF ELSE utime=utc
RT32_SUN_POSITION, year, month, day, utime, az, el, ra_sun, dec_sun, /REFRA; efemeridi Solnca - gradusi
n=N_ELEMENTS(utime)
q=DBLARR(n)
FOR i=0, n-1 DO BEGIN
	RT32_LST_TRUE, year, month, day, utime[i], lst				; lokalnoe sidericheskoe vremya - chasi
	h=15.*lst-ra_sun[i]											; chasovoy ugol Solnca v gradusah
	q_tan=SIN(!DTOR*h)/(TAN(lat_rt32)*COS(!DTOR*dec_sun[i])-SIN(!DTOR*dec_sun[i])*COS(!DTOR*h))
	q[i]=-!RADEG*ATAN(q_tan)										; parallakticheskiy ugol, gradusi, "+" nalon vlevo protiv chasovoy strelki, "-" nalon vpravo po chasovoy strelke
ENDFOR
IF n EQ 1 THEN RETURN, q[0] ELSE RETURN, q
end