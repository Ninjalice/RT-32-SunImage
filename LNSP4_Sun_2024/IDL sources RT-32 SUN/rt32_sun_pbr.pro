;--------------------------------------
pro rt32_sun_pbr, year, month, day, p, b0, r_sun
;
;vozvrascaet pozicionniy ugol, naklon osi Solnca, vidimiy radius Solnca
;na datu
;
; po J. Meeus Astronomical Algorithms
;
;ispolzuet
; JCNV, NUTATE

JDCNV, year,month,day,0., jd;      julianskij den na nachalo sutok

teta=DOUBLE(((jd-2398220.)*360./25.38)) MOD 360.;
i=7.25;                     naklon ekvatora Solnca k ploskosti ekliptiki, gradusi
k=73.6667+1.395833*(jd-2396758.)/36525; dolgota voshodyaschego uzla

t=DOUBLE((jd-2451545.)/36525.);           stoletij ot 2000

l0=(280.46645+DOUBLE(36000.76983*t)+DOUBLE(0.0003032*t*t)) MOD 360.;    geometricheskaia srednyaya dolgota Solnca otn. ravnodenstviya na tekuschuyu datu, gradusi
IF l0 LT 0. THEN l0=l0+360.

m0=(357.52910+35999.05030*t $
    -0.0001559*t*t-0.00000048*t*t*t) MOD 360.; srednyaya anomalija gradusi
IF m0 LT 0. THEN m0=m0+360.

ex=0.016708617-0.000042037*t;                         excentrisitet orbiti Zemli
c=(1.914600 -0.004817*t -0.0000014*t*t)*SIN(   !DTOR*m0) $
 +(0.019993 -0.000101*t               )*SIN(2.*!DTOR*m0) $
 + 0.000290                            *SIN(3.*!DTOR*m0);   uravnenie centra Solnca

lambda=(l0+c) MOD 360.;
IF lambda LT 0. THEN lambda=lambda+360. ;             istinnaya dolgota Solnca otn tochki ravnodentviya na tekuschuyu datu

anomaly=(m0+c) MOD 360.;
IF anomaly LT 0. THEN anomaly=anomaly+360.;           istinnaya anomaliya

rr=1.000001018*(1.-ex*ex)/(1.+ex*COS(!DTOR*anomaly));     radius vektor (a.e.)
r_sun=!RADEG*ATAN(0.004652, rr)*3600.;             vidimuy radius Solnca (ugl.sec)

abb=(-20.4898/rr)/3600.;                        popravka na abberaciyu s gradusah
lambda_app=(lambda+abb) MOD 360.;                vidimaya dolgota s uchetom abberacii, no bez nutacii

NUTATE, jd, nut_lon, nut_obliq;                    popravki na nutaciyu
lambda_app_=lambda_app+nut_lon/3600;                 vidimaya dolgota, no s uchetom nutacii po dolgote


e0=TEN(23,26,21.446)-DOUBLE(TEN(0,0, 46.815)*t) $
       -DOUBLE(TEN(0,0,0.00059)*t*t)+DOUBLE(TEN(0,0,0.001813)*t*t*t); sredniy naklon ekliptiki na tekuschuyu epohu v gradusah
e=e0+nut_obliq/3600.;                        istinny naklon s uchetom nutacii

tanx=-COS(!DTOR*lambda_app_   )*TAN(!DTOR*e)
tany=-COS(!DTOR*(lambda_app-k))*TAN(!DTOR*i)
x=!RADEG*ATAN(tanx)
y=!RADEG*ATAN(tany)
p=x+y;                             						   ;pozicionny ugol (gradusi)
b0=!RADEG*ASIN(SIN(!DTOR*(lambda_app-k))*SIN(!DTOR*i))   ; naklon osi Solnca (gradusi)

end