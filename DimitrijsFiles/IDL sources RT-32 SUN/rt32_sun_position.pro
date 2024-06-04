pro RT32_SUN_POSITION, year_, month_, day_, utc_, az, el, ra, dec, DUT1=dut1, REFRA=refra
;
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;
; Calculation of Sun local horizontal coordinates on date&time for Irbene RT-32 radio telescope
;
;INPUT:
; utc - scalar or vector, Universal Time Coordinared , hours
; keyword DUT1 - correction of UTC to UT1,  UT1=UTC+DUT1,  seconds, from ftp://maia.usno.navy.mil/ser7/ser7.dat
; keyword REFRA -correction of elevation for radio refrection (1.2 of optical refraction)
;
;OUTPUT:
; az, el - Sun horizontal coordinates of RT32 local horizon, degrees, scalar or vector, same dimension as utc
; ra, dec - true (apparent) Sun equatorial coordinates on current epoque, degrees, scalar or vector, same dimension as utc
;
;-----------------------------------------------------------
year=year_
month=month_
day=day_
utc=utc_
;
lat_rt32=57.5535171694D0     ;RT32 geographic latitude 57d 33' 12''.661
lon_rt32=21.8545525000D0     ;RT32 geographic longitude 21d 51' 16''.389
;-------------------correction of geographic latitude to geocentric latitude--------------------
lat_rt32_geocentric=!RADEG *ATAN(0.993305D0*TAN(!DTOR*lat_rt32))  ; b/a =0.996647 - elipcity of globe
;-------------------------------forming radio refraction = 1.2 mean optical---------------------------------------
r_refra_arg=1.+FINDGEN(90)            ; heights for radio refraction, degrees
r_refra=[1.0000,0.6461,0.3461,0.2597,0.2085,0.1747,$
0.1477,0.1307,0.1173,0.1060,0.0967,0.0890,$
0.0820,0.0763,0.0710,0.0667,0.0627,0.0590,$
0.0557,0.0527,0.0500,0.0477,0.0453,0.0433,$
0.0410,0.0397,0.0380,0.0363,0.0347,0.0333,$
0.0320,0.0310,0.0296,0.0287,0.0277,0.0267,$
0.0257,0.0247,0.0240,0.0233,0.0225,0.0217,$
0.0209,0.0201,0.0193,0.0187,0.0181,0.0175,$
0.0169,0.0163,0.0158,0.0157,0.0147,0.0142,$
0.0137,0.0132,0.0127,0.0127,0.0118,0.0113,$
0.0107,0.0104,0.0099,0.0094,0.0090,0.0086,$
0.0082,0.0078,0.0074,0.0070,0.0067,0.0063,$
0.0060,0.0056,0.0053,0.0049,0.0045,0.0041,$
0.0037,0.0033,0.0030,0.0027,0.0023,0.0020,$
0.0017,0.0013,0.0010,0.0007,0.0003,0.0000] ;radio refraction on 2-10 cm = 1.2 * mean optical refraction,, degrees
;
IF KEYWORD_SET(DUT1) THEN d_ut1=dut1/3600. ELSE d_ut1=0.
a=SIZE(utc)                				;create array if utc is scalar
IF a[0] EQ 0 THEN BEGIN
    utime=DBLARR(1)
    utime[0]=utc+d_ut1
ENDIF ELSE utime=utc+d_ut1
;
;forming julian day on yyyy-mm-ddT00:00:00
IF month LE 2 THEN BEGIN
    month=month+12
    year=year-1
ENDIF

c1=FIX(year/100)
c2=FIX(2-c1+FIX(c1/4))
jul_day0 =LONG(365.25D0*(year+4716)) + LONG(30.6001D0*(month+1)) +LONG(day) + LONG(c2) -1524.5D0;julian day on yyyy-mm-ddT00:00:00
;jul_day_now= DOUBLE(jul_day0)+DOUBLE(utime/24.D0)					;julian day on yyyy-mm-ddThh:mm:ss (array)
jul_cent_1900= ((jul_day0 + utime/24.D0) - 2415020.0D0)/36525.0D0					;julian centuries from 1900.0 1899-12-31 T 12:00 till yyyy-mm-ddThh:mm:ss
jul_cent_2000= (jul_day0 - 2451545.0D0)/36525.0D0					;julian centuries from 2000.0 2000-1-3 T 12:00
;
;---------------------nutation------------------------------------------------------------------------------
omega=125.04452D0-1934.136261D0*jul_cent_2000+0.0020708d0*jul_cent_2000*jul_cent_2000+(jul_cent_2000*jul_cent_2000*jul_cent_2000)/450000.D0 ;longitude of ascending node of Moon mean orbit, degrees
lon_moon_mean=218.3165D0 + 481267.8813D0*jul_cent_2000				        ;Moon mean longitude, degrees
lon_sun_mean =280.4665D0 +  36000.7698D0*jul_cent_2000     					;Sun mean longitude,  degrees
nut_lon=-17.20D0*SIN(!DTOR*omega) -1.32D0*SIN(2*!DTOR*lon_sun_mean) -0.23D0*SIN(2*!DTOR*lon_moon_mean) +0.21D0*SIN(2*!DTOR*omega); simplified
nut_obl=  9.20D0*COS(!DTOR*omega) +0.57D0*COS(2*!DTOR*lon_sun_mean) -0.10D0*COS(2*!DTOR*lon_moon_mean) +0.09D0*COS(2*!DTOR*omega); nutation on date&time, arc. sec.
;
;-----------mean and true obliquity---------------------------------------------------------------------
obl_mean=(23.d0+26.D0/60.D0+21.448D0/3600.D0)- (46.815D0/3600.d0)*jul_cent_2000-(0.00059D0/3600.D0)*jul_cent_2000*jul_cent_2000+ (0.001813D0/3600.D0)*jul_cent_2000*jul_cent_2000*jul_cent_2000; mean oliquity of equliptic, degrees
obl_true =23.452294D0-0.0130125D0*jul_cent_1900 + (9.2d0*COS(omega*!DTOR))/3600.0D0
; --------------------------local mean and apparent siderial time----------------------------------------
a1=24110.54841D0
a2=8640184.812866D0*jul_cent_2000;
a3=0.093104D0*jul_cent_2000*jul_cent_2000
a4=0.0000062D0*jul_cent_2000*jul_cent_2000*jul_cent_2000
gmst0= a1+a2+a3-a4
gmst0=gmst0 MOD 86400.0D0							;Greenwhich Mean Siderial Time on yyyy-mm-ddT00:00, seconds
gmst =gmst0 + 1.00273790935D0*utime*3600.0D0		;Greenwhich Mean Siderial Time on yyyy-mm-ddThh:mm, seconds (array)
gast =gmst  + nut_lon*COS(!DTOR*obl_true)/15.0D0	;Greenwhich Apparent Siderial Time on yyyy-mm-ddThh:mm, seconds (array)
ltst= gast  + 3600.0D0*lon_rt32/15.0D0  				;local Apparent Siderial Time on yyyy-mm-ddT hh:mm:ss,  seconds (attary)
;--------------------- sun position from astrolib, truncated Newcomb's algorithm ----------------------------
lon_sun_mean=DOUBLE((279.696678D0+((36000.768925D0*jul_cent_1900) MOD 360.0D0))*3600.0D0) 	;mean Sun longitude
;  allow for ellipticity of the orbit (equation of centre) using the Earth's mean anomaly ME
me = 358.475844D0 + ((35999.049750D0*jul_cent_1900) MOD 360.0D0)
ellcor  = (6910.1D0 - 17.2D0*jul_cent_1900)*SIN(!DTOR*me) + 72.3D0*SIN(2.0D0*!DTOR*me)
lon_sun_mean = lon_sun_mean + ellcor													;mean Sun longitude corrected for Earth orbit ellipcity
; allow for the Venus perturbations using the mean anomaly of Venus MV
mv = 212.603219D0 + ((58517.803875D0*jul_cent_1900) MOD 360.0D0)
vencorr = 4.8D0 * COS((299.1017D0 + mv - me)*!DTOR) + $
          5.5D0 * COS((148.3133D0 +  2.0D0 * mv  -  2.0D0 * me )*!DTOR) + 2.5D0 * COS((315.9433D0 +  2.0D0 * mv  -  3.0D0 * me )*!DTOR) + $
          1.6D0 * COS((345.2530D0 +  3.0D0 * mv  -  4.0D0 * me )*!DTOR) + 1.0D0 * COS((318.1500D0 +  3.0D0 * mv  -  5.0D0 * me )*!DTOR)
lon_sun_mean = lon_sun_mean + vencorr
;  Allow for the Mars perturbations using the mean anomaly of Mars MM
 mm = 319.529425D0  +  (( 19139.858500D0 * jul_cent_1900)  MOD  360.0D0 )
 marscorr = 2.0D0 * COS((343.8883D0 -  2.0D0 * mm  +  2.0D0 * me)*!DTOR ) + 1.8D0 * COS((200.4017D0 -  2.0D0 * mm  + me) * !DTOR)
 lon_sun_mean = lon_sun_mean + marscorr
; Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
mj = 225.328328D0  +  (( 3034.6920239D0 * jul_cent_1900) MOD 360.0D0 )
jupcorr = 7.2D0*COS((179.5317D0-mj+me)*!DTOR) + 2.6D0*COS((263.2167D0- mj)*!DTOR) + $
          2.7D0 * COS(( 87.1450D0  -  2.0D0 * mj  +  2.0D0 * me)*!DTOR) + 1.6D0 *COS((109.4933D0-2.0D0*mj+me)*!DTOR)
lon_sun_mean = lon_sun_mean + jupcorr
; Allow for the Moons perturbations using the mean elongation of the Moon from the Sun D
d = 350.7376814D0  + (( 445267.11422D0 * jul_cent_1900) MOD 360.0D0 )
mooncorr  = 6.5D0 * SIN(d*!DTOR)
lon_sun_mean = lon_sun_mean + mooncorr
; Allow for long period terms
longterm = 6.4D0 * SIN(( 231.19D0+20.20D0*jul_cent_1900)*!DTOR)
lon_sun_mean = lon_sun_mean + longterm
lon_sun_mean = (lon_sun_mean + 2592000.0D0) MOD 1296000.0D0 ; arc. sekundi
;longmed = lon_sun_mean/3600.0D0		;gradusi
lon_sun_mean = lon_sun_mean-20.5D0	;allow abberation, arc.sec.
; Allow for Nutation using the longitude of the Moons mean node OMEGA
omega = 259.183275D0 - ((1934.142008D0*jul_cent_1900)MOD 360.0D0)
lon_sun_mean = lon_sun_mean-17.2D0*SIN(omega*!DTOR)
; Form Right Ascension and Declination
lon_sun_mean=lon_sun_mean/3600.0D0		; gradusi
ra_sun=ATAN(SIN(lon_sun_mean*!DTOR)*COS(obl_true*!DTOR),COS(lon_sun_mean*!DTOR))
neg = WHERE(ra_sun LT 0.0D0, Nneg)
IF Nneg GT 0 THEN ra_sun[neg]=ra_sun[neg]+2.0d*!PI
dec_sun = ASIN(SIN(lon_sun_mean*!DTOR) * SIN(obl_true*!DTOR))
;----------- forming Sun horizontal position-----------------
hour_angle=DOUBLE(15.*DOUBLE(ltst/3600.0d0) - !RADEG*ra_sun)							;Sun hour angle, degrees
ell=DBLARR(N_ELEMENTS(hour_angle))
azz=DBLARR(N_ELEMENTS(hour_angle))
cos_dec=DOUBLE(COS(dec_sun))
sin_dec=DOUBLE(SIN(dec_sun))
cos_lat=DOUBLE(COS(!DTOR*lat_rt32))
;sin_lat=DOUBLE(SIN(!DTOR*lat_rt32_geocentric))
;cos_lat=DOUBLE(COS(!DTOR*lat_rt32_geocentric))
sin_lat=DOUBLE(SIN(!DTOR*lat_rt32))
FOR i= 0, N_ELEMENTS(hour_angle)-1 DO BEGIN
	IF hour_angle[i] LT 0. THEN hour_angle[i]=360.+hour_angle[i]
	cos_hu=DOUBLE(COS(!DTOR*hour_angle[i]))
	sin_hu=DOUBLE(SIN(!DTOR*hour_angle[i]))
	sin_el= DOUBLE(sin_dec[i]*sin_lat+cos_dec[i]*cos_lat*cos_hu)
	ele_rad=DOUBLE(ASIN(sin_el))
	ell[i]=!RADEG*ele_rad
	azi_rad=ATAN(sin_hu, cos_hu*sin_lat - (sin_dec[i]/cos_dec[i])*cos_lat)
	azz[i]=180.+!RADEG*azi_rad
	IF azz[i] LT 0. THEN azz[i]=360.+azz[i]
	IF azz[i] GE 360. THEN azz[i]=azz[i]-360.
	IF KEYWORD_SET(REFRA) THEN BEGIN
		d_refra=INTERPOL(r_refra, r_refra_arg,    ell[i])       ;correction to radio refraction, dergees
    	ell[i]=ell[i]+d_refra
    ENDIF
ENDFOR
;
;----------output array--------------------------
IF N_ELEMENTS(utc) EQ 1 THEN BEGIN
    az=azz[0]
    el=ell[0]
    ra =!RADEG* ra_sun[0]
    dec=!RADEG*dec_sun[0]
ENDIF ELSE BEGIN
    az=azz
    el=ell
    ra =!RADEG* ra_sun
    dec=!RADEG*dec_sun
ENDELSE
;
end

