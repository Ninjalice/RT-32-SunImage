function LNSP5ch_file_find, year, month, day, time_begin, time_end,f_list, PATH=path
;
;+++++++++++++++++++++++++++++++++++++++++
;otiskivaet v papke faila tipa LNSP4_5ch_xxxxxx_xxxx_xxxx.fit na datu-vremya
;vozaraschaet chislo naidennih failov
;INPUT:
; year, month, day - data
; time_begin, time_end - vremya nachala i konca iskomogo perioda vremeni
;OUTPUT:
;f_list - massiv stro imen naidennyh failov
;chislo naidennih failov v tekuscheq ili ukazannoy papke
;keywords:
;FILE_NAME_LIST - spisok naidennih failov
;PATH - ukazannaya papka, esli ego net - ischetsya v tekuschei
;
;
; stroka, sootvetstvujuschaya date
str_date_=''
str_date_= STRCOMPRESS(STRING(FIX(year-2000)),/REMOVE_ALL)
IF month LE 9 THEN str_date_=str_date_+'0'
str_date_= str_date_+STRCOMPRESS(STRING(FIX(month)),/REMOVE_ALL)
IF day LE 9 THEN str_date_=str_date_+'0'
str_date_= str_date_+STRCOMPRESS(STRING(FIX(day)),/REMOVE_ALL)
;ischem vse faili v direktorii
file_name_date='lnsp4_5ch_'+str_date_+'*.fit
IF KEYWORD_SET(PATH) THEN file_name=path+file_name_date ELSE file_name=file_name_date
all_files=STRARR(1)
all_files[0]=''
all_files=FILE_SEARCH(file_name)
IF STRLEN(all_files[0]) NE 0 THEN BEGIN
;nenulevaya dlina stroki - naydeni faili
	;prosmatrivaem imena
	n_file=0
	FOR i=0, N_ELEMENTS (all_files)-1 DO BEGIN
		str_start=STRMID(all_files[i],16,6,/REVERSE_OFFSET)
		str_end  =STRMID(all_files[i], 9,6,/REVERSE_OFFSET)
		IF n_file EQ 0 THEN BEGIN
			files_=all_files[i]
			t_start_file=TEN(FIX(STRMID(str_start,0,2)), FIX(STRMID(str_start,2,2)), FIX(STRMID(str_start,4,2)))
			t_end_file  =TEN(FIX(STRMID(str_end,  0,2)), FIX(STRMID(str_end,  2,2)), FIX(STRMID(str_end,  4,2)))
			n_file=n_file+1
		ENDIF ELSE BEGIN
			files_=[files_,all_files[i]]
			t_start_file=[t_start_file,TEN(FIX(STRMID(str_start,0,2)), FIX(STRMID(str_start,2,2)), FIX(STRMID(str_start,4,2)))]
			t_end_file  =[t_end_file,  TEN(FIX(STRMID(str_end,  0,2)), FIX(STRMID(str_end,  2,2)), FIX(STRMID(str_end,  4,2)))]
			n_file=n_file+1
		ENDELSE
	ENDFOR
	a=SORT(t_start_file)
	files_=files_(a)
	t_start_file=t_start_file(a)
	t_end_file  =t_end_file(a)
	index=WHERE((t_end_file GE time_begin) AND (t_start_file LE time_end))
	aa=SIZE(index)
	IF aa[0] NE 0 THEN BEGIN
		f_list=files_[index]
		RETURN, N_ELEMENTS(index)
	ENDIF	ELSE RETURN, 0
ENDIF ELSE RETURN, 0
end

