function RT32_DATETIME_tO_STR, year, month, day, utc
;vozvraschaet stroku vida yyy-mm-ddThh:mm:ss
;INPUT:
;year,month,day
;ut - UT vremya , chasi
aa=SIXTY(utc)
str=' '
date=STRTRIM(STRING(FIX(year)),2)+'-'
IF month LE 9 THEN date=date+'0'
date=date+STRTRIM(STRING(FIX(month)),2)+'-'
IF day LE 9 THEN date=date+'0'
date=date+STRTRIM(STRING(FIX(day)),2)
IF aa[0] LE 9 THEN time='0' ELSE time=''
time=time+STRTRIM(STRING(FIX(aa[0])),2)+':'
IF aa[1] LE 9 THEN time=time+'0'
time=time+STRTRIM(STRING(FIX(aa[1])),2)+':'
IF aa[2] LT 10. THEN time=time+'0'
;IF aa[2] LT 0.5 THEN a_sec=STRTRIM(STRING(FIX(FLOOR(aa[2]))),2) ELSE a_sec=STRTRIM(STRING(FIX(FLOOR(aa[2]+1.))),2)
time=time+STRTRIM(STRING(FIX(FLOOR(aa[2]+0.000501))),2)
str=date+'T'+time
RETURN, str
end