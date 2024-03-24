pro rt32_sun_spiral_scan_make
;
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++
;gotovit programmu spiralnogo skanirovaniya Solnca
;5 spiraley
;
;ispolzuet
;RT32_SUN_POSITION
;RT32_SUN_PARA
;RT32_DATETIME_TO_STR
;
;
;----------------------------------------------------


path='C:\Users\User\Desktop\LNSP4_5_2024\'
;ishodnie dannie dlya skanirovaniya - vvesti vruchnuyu
PRINT, "--------------------------------------------------------"
PRINT, "----------LNSP4 SUN SPIRAL SCANNING-----------------------"
;
;READ ,PROMPT='Start at yyyy mm dd , UTC hh, UTC mm', year, month, day, hour_start, minute_start
year=2024
month=2
day=15
hour_start=11
minute_start=25

az_offset0=+0.01
el_offset0=+0.18 	;offseti datchikov polozheniya antenni, gradusi
;
az_offset2=0. ;offset obluchatelya otnositelno ele osi antenni, gradusi
el_offset2=0.
;
num_scan=5
step1=6.; stupenki spirali arc.min.
step2=6.;
step3=6.;
step4=10.;
step5=12.
sky=16.*4	;chistoe nebo - 4 radiusa Solnca, arcmin
; intervali vremeni skanirovanija, sekundi
t_cal=60.; vremya kalibrovki
t1=40.
t2=100.;
t3=120;
t4=120.
t5=120.; vremya perehoda
t_slew=20.;
t_scan=t_cal+t1+t2+t3+t4+t5+t_slew+t_cal+t_slew;  dlitelnostj scana, sec
;
;sformitovatj imya faila
file_out=STRING(FIX(year-2000))
IF month LE 9 THEN file_out=file_out+'0'+STRING(FIX(month)) ELSE file_out=file_out+STRING(FIX(month))
IF day   LE 9 THEN file_out=file_out+'0'+STRING(FIX(day))   ELSE file_out=file_out+STRING(FIX(day))
file_out=STRCOMPRESS(file_out+"_",/REMOVE_ALL)
IF hour_start LE 9 THEN file_out=file_out+'0'+STRING(FIX(hour_start)) ELSE file_out=file_out+STRING(FIX(hour_start))
IF minute_start LE 9 THEN file_out=file_out+'0'+STRING(FIX(minute_start)) ELSE file_out=file_out+STRING(FIX(minute_start))
file_name1=STRCOMPRESS(('sun_scan_'+file_out+'.ptf'),/REMOVE_ALL)
OPENW, file1, path+file_name1, /GET_LUN
;
PRINTF, file1,"# Table of horizontal coordinates for Sun scanning"
PRINTF, file1,"# Celestial Object ID:  SUN"
PRINTF, file1,"# Spiral scanning, number of scans:      "+STRCOMPRESS(STRING(FIX(num_scan)),/REMOVE_ALL)
PRINTF, file1,"# Scan time schedule:"
PRINTF, file1,"# Sun center calibration 	60 sec"
PRINTF, file1,"# 1 spiral turn				40 sec  step 0->6 arcmin"
PRINTF, file1,"# 2 spiral					100 sec  step 6->12 arcmin"
PRINTF, file1,"# 3 spiral					120 sec step 12->18 arcmin"
PRINTF, file1,"# 4 spiral					120 sec step 18-28 arcmin"
PRINTF, file1,"# 5 spiral					120 sec step 28-40 arcmin"
PRINTF, file1,"# Slew to sky 4 Rsun			20 sec"
PRINTF, file1,"# Sky calibration			60 sec"
PRINTF, file1,"# Slew to Sun center			20 sec"
RT32_SUN_POSITION, year, month, day, TEN(hour_start,minute_start,0), az_start, el_start, /REFRA
PRINTF, file1,"# Az@Start Time      : "+STRING(az_start)+" degr."
PRINTF, file1,"# El@Start Time      : "+STRING(el_start)+" degr."
PRINTF, file1,"# Start Time (UTC+0) : "+RT32_DATETIME_TO_STR(year, month,day,TEN(hour_start,minute_start,0))
utime_start=TEN(hour_start, minute_start,0) ; nachalo sessii, chasi
time_session=num_scan*t_scan/3600. 				;dlitelnostj sessii, chasi
utime_mean= utime_start+time_session/2.		;privedennoe UT vremya sessii
utime_end = utime_start+time_session;		;konec UT sessii, chasi
PRINTF, file1,"# End Time, (UTC), hours:"+RT32_DATETIME_TO_STR(year, month,day,utime_end)
PRINTF, file1,"# Mean observation time (UTC):"+RT32_DATETIME_TO_STR(year, month,day,utime_mean)
PRINTF, file1,"# Az offset          : "+string(az_offset0)+" degr."
PRINTF, file1,"# El offset          : "+string(el_offset0)+" degr."
PRINTF, file1,"#!!! Offset of the antenna pointing system must be set 0 !!!
PRINTF, file1," "
PRINTF, file1,"[Interpolation]"
PRINTF, file1,"Newton"
PRINTF, file1," "
PRINTF, file1,"[Load Mode]"
PRINTF, file1,"New"
PRINTF, file1," "
PRINTF, file1,"[Start Time]"
PRINTF, file1,RT32_DATETIME_TO_STR(year, month,day,TEN(hour_start,minute_start,0))
PRINTF, file1," "
PRINTF, file1,"[Table Data]"
PRINTF, file1,"#  Time   Az source [degr.]     El source[degr.]"

; sdvig ot centra Solnca v techenii odnjgo scana
x=DBLARR(700)
y=DBLARR(700);
FOR i=0,FIX(t_cal)-1 DO BEGIN; calibration Sun center
	x[i]=0.
	y[i]=0.
ENDFOR
FOR i=0,FIX(t1)-1 DO BEGIN; 1 turn
	i0=FIX(t_cal)
	fi=i*360./t1
	r=i*step1/t1
	x[i+i0]=r*COS(!DTOR*fi)
	y[i+i0]=r*SIN(!DTOR*fi)
ENDFOR
FOR i=0,FIX(t2)-1 DO BEGIN; 2 turn
	i0=FIX(t_cal)+FIX(t1)
	fi=i*360./t2
	r=step1+i*(step2/t2)
	x[i+i0]=r*COS(!DTOR*fi)
	y[i+i0]=r*SIN(!DTOR*fi)
ENDFOR
FOR i=0,FIX(t3)-1 DO BEGIN; 3 turn
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)
	fi=i*360./t3
	r=step1+step2+i*(step3/t3)
	x[i+i0]=r*COS(!DTOR*fi)
	y[i+i0]=r*SIN(!DTOR*fi)
ENDFOR
FOR i=0,FIX(t4)-1 DO BEGIN; 4 turn
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)
	fi=i*360./t4
	r=step1+step3+step3+(i*step4/t4)
	x[i+i0]=r*COS(!DTOR*fi)
	y[i+i0]=r*SIN(!DTOR*fi)
ENDFOR
FOR i=0, FIX(t5)-1 DO BEGIN; 5 turn
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)+FIX(t4)
	fi=i*360./t5
	r=step1+step2+step3+step4+i*(step5/t5)
	x[i+i0]=r*COS(!DTOR*fi)
	y[i+i0]=r*SIN(!DTOR*fi)
ENDFOR
FOR i=0, FIX(t_slew)-1 DO BEGIN; slew to clibration, Sky
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)+FIX(t4)+FIX(t5)
	y[i+i0]=0.
	x0=step1+step2+step3+step4+step5
	x[i+i0]=x0+i*(sky-x0)/t_slew
ENDFOR
FOR i=0, FIX(t_cal)-1 DO BEGIN; calibration sky
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)+FIX(t4)+FIX(t5)+FIX(t_slew)
	x[i+i0]=sky
	y[i+i0]=0.
ENDFOR
FOR i=0, FIX(t_slew)-1 DO BEGIN; slew to Sun center
	i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)+FIX(t4)+FIX(t5)+FIX(t_slew)+FIX(t_cal)
	y[i+i0]=0.
	x[i+i0]=(FIX(t_slew)-i-1)*sky/t_slew
ENDFOR
i0=FIX(t_cal)+FIX(t1)+FIX(t2)+FIX(t3)+FIX(t4)+FIX(t5)+FIX(t_slew)+FIX(t_cal)+FIX(t_slew)

; sdelatj num_scan scanov po i0 tochek
xx=DBLARR(num_scan*i0)
yy=DBLARR(num_scan*i0)
FOR j=0,num_scan-1 DO BEGIN
	ff=j*(360./num_scan)
	FOR i=0,i0-1 DO BEGIN
		ii=j*i0+i
		xx[ii]= x[i]*COS(!DTOR*ff)- y[i]*SIN(!DTOR*ff)
		yy[ii]= x[i]*SIN(!DTOR*ff)+ y[i]*COS(!DTOR*ff)
	ENDFOR
ENDFOR

utc=utime_start+FINDGEN(num_scan*i0)/3600.	;massiv, utc cherez sekundu ot utime_start, chasi vremeni
q=rt32_sun_para(year, month, day, utc)		;massiv, paralaksniy ugol v kazhdoi tochke na utc
xx1= xx*COS(!DTOR*q) - yy*SIN(!DTOR*q)
yy1= xx*SIN(!DTOR*q) + yy*COS(!DTOR*q)		; povorot obeih koordinat smescheniya na paralaksniy ugol, + povorot vlevo protiv chasstrelki
RT32_SUN_POSITION, year, month, day, utc, az_sun, el_sun, /REFRA	;az-el koordinati centra solnca, s refrakciey!
az_anten=az_sun + xx1/COS(!DTOR*el_sun)/60.
el_anten=el_sun + yy1/60.



;zapisivatj tablicu v fail
FOR i=0, num_scan*i0-1 DO BEGIN
	PRINTF, file1, RT32_DATETIME_TO_STR(year,month,day, utc[i])+'.00'+STRING(az_anten[i])+STRING(el_anten[i])
ENDFOR



PRINTF, file1,"#---------------END---------------"
CLOSE,/ALL

print,'-------------------------------------------------------------'
print, 'Saved: ', file_name1, "  ", num_scan*i0, "  tochek"


window, 13, xsize=800,ysize=800
plot, az_anten, el_anten,xtitle='azimuth', ytitle='elevation',title='RT32 antenna sun scanning Start '+RT32_datetime_to_str(year,month,day,utime_start),/isotropic

end
