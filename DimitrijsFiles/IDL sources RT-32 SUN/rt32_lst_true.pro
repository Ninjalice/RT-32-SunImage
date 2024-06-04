
pro RT32_LST_TRUE, year_, month_, day_, utc, lst, DUT1=dut1
;
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; Calculation of Sun Local True Siderial time on date&time for Irbene RT-32 radio telescope
;
;INPUT:
; year, month, day - date
; utc - scalar or vector, Universal Time Coordinared , hours
; keyword DUT1 - correction of UTC to UT1,  UT1=UTC+DUT1,  seconds, from ftp://maia.usno.navy.mil/ser7/ser7.dat
;
;OUTPUT:
; lst - Local True Siderial time, hours
;
;
;-----------------------------------------------------------
year=year_
month=month_
day=day_
lat_rt32=57.5535171694D0     ;RT32 geographic latitude 57d 33' 12''.661
lon_rt32=21.8545525000D0     ;RT32 geographic longitude 21d 51' 16''.389
IF KEYWORD_SET(dut1) THEN d_ut1=dut1 ELSE d_ut1=0.
a=size(utc)                ;create array if utc is scalar
IF a[0] EQ 0 THEN BEGIN
    utime=DBLARR(1)
    utime[0]=utc+d_ut1
ENDIF ELSE utime=utc+d_ut1
;
IF month LE 2 THEN BEGIN    													;julian day on yyyy-mm-ddT00:00:00
    month=month+12
    year=year-1
ENDIF
jd2000=2451545.0d0                                        						; J2000 2000-1-1T12:00:00; julian day on J2000
c1=FIX(year/100)
c2=FIX(2-c1+FIX(c1/4))
jul_day0 =DOUBLE(LONG(365.25*(year+4716))+LONG(30.6001*(month+1))+day+c2-1524.5);julian day on yyyy-mm-ddT00:00:00
jul_day_now= DOUBLE(jul_day0)+DOUBLE(utime/24.d0)								;julian day on yyyy-mm-ddThh:mm:ss (array)
mjul_day0=jul_day0-jd2000                                                  		;julian days from 2000-01-01T12:00 to yyyy-mm-ddT00:00
;mjd_now=DOUBLE(mjul_day0)+DOUBLE(utime[0]/24.d0)                                ;julian days from 2000-01-01T12:00:00 to yyyy-mm-ddThh:mm:ss
t100=DOUBLE(mjul_day0/36525.0d0)												;centuries   from 2000-01-01T12:00:00 to yyyy-mm-ddT00:00:00
a1=DOUBLE(24110.54841d0)
a2=DOUBLE(8640184.812866d0*t100);
a3=DOUBLE(0.093104d0*t100*t100)
a4=DOUBLE(0.0000062d0*t100*t100*t100)
gmst0= DOUBLE(a1 +a2+ a3 - a4)
WHILE gmst0 GE 86400.0d0 DO gmst0=gmst0-86400.0d0 								;Greenwhich Mean Siderial Time in seconds on yyyy-mm-ddT00:00
lmst0=gmst0+3600.0d0*lon_rt32/15.0d0         									;local Mean Siderial Time in seconds on yyyy-mm-ddT00:00
;
lon_moon_mean=218.3165d0 + 481267.8813d0*t100    					           	;Moon mean longitude, degrees
lon_sun_mean =280.4665d0 +  36000.7698d0*t100              						;Sun mean longitude, degrees
omega=125.04452d0-1934.136261d0*t100+0.0020708d0*t100*t100+(t100*t100*t100)/450000.d0 ;longitude of ascending node of Moon mean orbit, degrees
obl=(23.d0+26.d0/60.d0+21.448d0/3600.d0)- (46.815D0/3600.d0)*t100-(0.00059d0/3600.d0)*t100*t100+ (0.001813d0/3600.d0)*t100*t100*t100; mean oliquity of equliptic, degrees
nut_lon=-17.20d0*SIN(!DTOR*omega) -1.32d0*SIN(2*!DTOR*lon_sun_mean) -0.23d0*SIN(2*!DTOR*lon_moon_mean) +0.21d0*SIN(2*!DTOR*omega); simplified
nut_obl=  9.20d0*COS(!DTOR*omega) +0.57d0*COS(2*!DTOR*lon_sun_mean) -0.10d0*COS(2*!DTOR*lon_moon_mean) +0.09d0*COS(2*!DTOR*omega); nutation on date&time, arc. sec.
;
ltst0=lmst0+nut_lon*COS(!DTOR*obl)/15.    								;local Apparent (True) Siderial Time in seconds on yyyy-mm-ddT00:00
ltst=DOUBLE(ltst0+1.00273790935d0*DOUBLE(utime*3600.0d0))			    ;local Apparent (True) Siderial Time in seconds on yyyy-mm-ddT utime
ltst=DOUBLE(ltst/3600.)													; L A S T in hours
FOR i=0, N_ELEMENTS(ltst)-1 DO WHILE ltst[i] GE 24. DO ltst[i]=ltst[i]-24.d0
;
IF N_ELEMENTS(utc) EQ 1 THEN BEGIN
    lst=ltst[0]
ENDIF ELSE BEGIN
    lst=ltst
ENDELSE
end