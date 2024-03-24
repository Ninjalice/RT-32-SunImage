function LNSP5CH_FILE_READ, num_ch, year, month, day, utime_begin, utime_end, utrcp, utlcp, rcp, lcp, PATH=path
;

;++++++++++++++++++++++++++++++++++++++++++
;chtenie FITS failov vida lnsp4_5ch_xxxxxx_xxxxxx_xxxxxx;
;
;vozvraschet chislo otschetov
;
; INPUT:
; num_ch 					-nomer kanala v rezhime 5 kanalov, 0-4
; year, month, date 		- data
; utime_begin, utime_end 	- UTC otrezok vremeni, chasi
;
;OUTPUT:
; utrcp,utlcp - vremja otschetov samolov obeih polarizacij, chasi
; rcp, lcp - otschet rcp i lcp, biti ACP
;
; KEYWORD:
;PATH - direktoriy, gde iskatj faili dannih, esli otsutstvuet - ischetsya v tekuschem
;;
;ispolzuet LNSP5CH_FILE_FIND, FXBOPEN, FXBREAD, FXBCLOSE iz astrolib
;
;----------------------------------------------------
; chislo failov s iskomimi dannimi
list=''
IF KEYWORD_SET(PATH) THEN BEGIN
	file_num=LNSP5CH_FILE_FIND( year, month, day, utime_begin, utime_end, list, PATH=path)
ENDIF ELSE BEGIN
	file_num=LNSP5CH_FILE_FIND( year, month, day, utime_begin, utime_end, list, PATH=path)
ENDELSE
;
IF file_num NE 0 THEN BEGIN
	FOR i=0, file_num-1 DO BEGIN
		FXBOPEN, unit, list[i], 1, header
		FXBREAD, unit,utrcp,4*num_ch+1
		FXBREAD, unit,rcp,4*num_ch+2
		FXBREAD, unit,utlcp,4*num_ch+3
		FXBREAD, unit,lcp, 4*num_ch+4
		FXBCLOSE, unit
	ENDFOR
	;otrezatj zadannuyu chastj po vremeni
	index=WHERE( (utrcp GE utime_begin) AND (utrcp LE utime_end))
	utrcp=utrcp[index]
	rcp=rcp[index]
	index=WHERE( (utlcp GE utime_begin) AND (utlcp LE utime_end))
	utlcp=utlcp[index]
	lcp=lcp[index]
	RETURN, N_ELEMENTS(utrcp)
ENDIF ELSE RETURN, 0
;
end