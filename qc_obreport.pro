pro qc_obreport, filename=filename, obid=obid, $
 mjdrange=mjdrange, $
 survey=survey, $
 casu=casu, vircam=vircam, vhs=vhs, $
 verbose=verbose, debug=debug, tile=tile, $
 batch=batch, qcplots=qcplots, $
 qcstatus=qcstatus
;+
;$Id: qc_paw.pro,v 1.1 2009/12/12 22:08:04 rgm Exp rgm $
;
; NAME:
;   qc_obreport
;
; PURPOSE:
;   Generate report on executed OBs
;
;   Need to split OB report into wavebands and/or templates and give stats
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2009/12/12 22:08:04 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.1 $ 
;
;-
COMPILE_OPT IDL2

!x.style=1  & !y.style=1
!p.background=fsc_color('white')
!p.color=fsc_color('black')

;help, !version,/str
message,/inf, systime(/utc) 
message,/inf, getenv('USER')
message,/inf, getenv('PWD')
traceback=scope_traceback()
message,/inf, traceback[0]
message,/inf, traceback[1]

splog,'verbose: ',verbose

if not keyword_set(debug) then debug=0
if not keyword_set(verbose) then verbose=0
if keyword_set(debug) then verbose=2

splog,'verbose: ',verbose
MESSAGE,/INF,'Reading: '+ filename
t0=systime(1)
; read in the dqc data within mjd range and ra, dec range

dqc_read_data, data, filename=filename, $
 casu=casu, wsa=wsa, vsa=vsa, vircam=vircam, vhs=vhs, $
 mjdrange=mjdrange, survey=survey, $
 verbose=verbose

message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
print,'DQC file read in: ',n_elements(data),' records'

qcstatus=data.qcstatus
itest=UNIQ(qcstatus, SORT(qcstatus))
n_qcstatus=n_elements(itest)
message,/inf,'Number of unique QC status values: ' + string(n_qcstatus, '(i6)')
splog, qcstatus[itest]

itest=where(strpos(qcstatus, 'A') lt 0 and $ 
 strpos(qcstatus,'B') lt 0 , count)
data=data[itest]
splog, 'Numvber of QCSTATUS non-AB OB records: ', n_elements(itest)
pause


help, data
if keyword_set(debug) then help, data, /str
if keyword_set(verbose) then structure_stats, data, /sort

if not keyword_set(tile) then begin
  itest=where(data.is_tile eq 'F', count)
  data=data[itest]
endif

message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)

obid=data.obsid
ipoint=UNIQ(obid, SORT(obid))
n_obids=n_elements(ipoint)
message,/inf,'Number of unique OBids: ' + string(n_obids, '(i6)')
; generate a list of unique OBids
obid_unique=data[ipoint].obsid
obsname_unique=data[ipoint].obsname

message,/inf, 'Summary info for each unique obid '
print, 'i obid nfiles obname'
format='(i7, i8, i6, 3f10.4, 4i6, 8f8.2, 4x, a)'

; save the OB centres for each OB; check is this is also a keyword
ramean=make_array(n_obids,/double)
decmean=make_array(n_obids,/double)

if (n_obids gt 0) then begin
  for i=0, n_obids-1 do begin

    itest=where(data.obsid eq obid_unique[i],count)
    ntest=n_elements(itest)

    decmin=$
     min([data[itest].dec1,data[itest].dec2,data[itest].dec3,data[itest].dec4])
    decmax=$
     max([data[itest].dec1,data[itest].dec2,data[itest].dec3,data[itest].dec4])
    decmean[i]=total(data[itest].dec)/ntest
    decrange = decmax - decmin

    ramin=$
     min([data[itest].ra1,data[itest].ra2,data[itest].ra3,data[itest].ra4])
    ramax=$
     max([data[itest].ra1,data[itest].ra2,data[itest].ra3,data[itest].ra4])
    ramean[i]=total(data[itest].ra)/ntest
    rarange=(ramax-ramin)*cos(decmean[i]/!dradeg)
    ;rarange=(ramax-ramin)

    ;mjdobsmin=min(data[itest].mjdobs)
    ;mjdobsmax=max(data[itest].mjdobs)
    mjdobsmin=min(data[itest].mjd)
    mjdobsmax=max(data[itest].mjd)
    mjdobsrange=mjdobsmax-mjdobsmin
    tplnexpmin=min(data[itest].tplnexp)
    tplnexpmax=max(data[itest].tplnexp)
    tplexpnomin=min(data[itest].tplexpno)
    tplexpnomax=max(data[itest].tplexpno)

    datatest=data[itest]
    ; number of different template numbers

    print, format=format, $
     i+1, obid_unique[i], ntest, $
     mjdobsmin-55000, mjdobsmax-55000, mjdobsrange*24.0, $
     tplnexpmin, tplnexpmax, $
     tplexpnomin, tplexpnomax, $
     ramean[i], ramin, ramax, rarange, $
     decmean[i], decmin, decmax, decrange, $
     strtrim(data[ipoint[i]].obsname,2)

  endfor
endif

pause

; find OB that observe the same region of sky
radius=900.0

nearneigh_self, ramean, decmean, radius, xy=xy, $
 m1=m1, m2=m2, $
 dra=dra, ddec=ddec, dr_gcirc=dr_gcirc

nmatch=n_elements(m1)
message,/inf,'Number of overlaping OBs: ' + string(nmatch)
print,m1[0]
nmatch=0
if nmatch gt 0 then begin
for i=0, nmatch-1 do begin
  print, format='(i4,i8,i8,4f10.4,f8.4,2x,a,a,2x,a,4x,a)', $
   i+1, obid_unique[m1[i]], obid_unique[m2[i]], $
   ramean[m1[i]],decmean[m1[i]], $
   ramean[m2[i]],decmean[m2[i]], dr_gcirc[i], $
   obsname_unique[m1[i]], obsname_unique[m2[i]], $
   strtrim(adstring(ramean[m1[i]],decmean[m1[i]], precision=0),2), $
   strtrim(adstring(ramean[m2[i]],decmean[m2[i]], precision=0),2)
endfor
pause

nearneigh_plots, dra, ddec, dr_gcirc, $
 radius=radius, $
 /summary

radius=60.0
nearneigh_plots, dra, ddec, dr_gcirc, $
 radius=radius, $
 /summary

pause
endif

isort=sort(obid)
ndata=n_elements(obid)
obid_sorted=obid[isort]
message, /inf, 'Number of rows: ' + string(ndata)

; print info about each image in OBid order for debugging
if keyword_set(debug) then begin
  for i=0, ndata-1 do begin
    format='(i6,i8,f12.4,2x,a,4i10)' 
    print, format=format, i, $
     data[isort[i]].obsid, $
     data[isort[i]].mjd, $
     strtrim(mjd_isodate(data[isort[i]].mjd),2), $
     data[isort[i]].nightobs, $
     data[isort[i]].tplstart, $
     data[isort[i]].obsstart, $
     data[isort[i]].dateobs
 endfor
endif
message,/inf




; cycle through the list of unique OBs and derive some statistics
; obid_unique contains a list of unique obids

; print some OB info
debug=1
for i=0, n_obids-1 do begin

  itest=where(data.obsid eq obid_unique[i],count)
  ntest=n_elements(itest)
  splog, 'i, obid, count, ntest: ',i,obid_unique[i],count,ntest
  ;print, i+1, obid_unique[i], ntest, itest[0], obid[isort[itest[0]]]
 
    ;determine the unique nightobs for each OB and the LST range
    nightobs=data[itest].nightobs
    inightobs=UNIQ(nightobs, SORT(inightobs))
    print, i+1, n_elements(inightobs), nightobs[inightobs]
    print, i+1, $
     data[itest[inightobs]].nightobs, $
     data[itest[inightobs]].tplstart, $
     data[itest[inightobs]].obsstart, $
     data[itest[inightobs]].dateobs

    lst=data[itest].lst
    print, i+1, min(lst), max(lst), (max(lst)-min(lst))/60.0
    mjdobs=data[itest].mjd
    print, i+1, min(mjdobs), max(mjdobs), (max(mjdobs)-min(mjdobs))*24.0*60.0

    obid=obid_unique[i]
    if not keyword_set(batch) then batch=1
    verbose=0
    data_save=data
    qc_ob, filename=filename, /casu, /vircam, /vhs, obid=obid, $
     verbose=verbose, batch=batch, data=data
    data=data_save
 
endfor


end
