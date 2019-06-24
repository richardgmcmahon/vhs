pro qc_tilereport, filename=filename, obid=obid, $
 mjdrange=mjdrange, $
 survey=survey, $
 casu=casu, vircam=vircam, vhs=vhs, $
 verbose=verbose, debug=debug, $
 resultfile=resultfile, warningfile=warningfile
;+
;$Id: qc_tilereport.pro,v 1.2 2010/01/10 02:03:52 rgm Exp rgm $
;
; NAME:
;   qc_tilereport
;
; PURPOSE: 
;   Generate report on executed tiles. See also qc_obreport
;
;
; COMMENTS:
;   Since some tiles straddle 24hrs there is some 24hr wrap
;   handling. This is incomplete, so beware.
;
;   OBs that are paused and then restarted produce a warning
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2010/01/10 02:03:52 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.2 $ 
;
;
;  $Log: qc_tilereport.pro,v $
;  Revision 1.2  2010/01/10 02:03:52  rgm
;  *** empty log message ***
;
;
;
;-
COMPILE_OPT IDL2

;help, !version,/str
message,/inf, systime(/utc) 
message,/inf, getenv('USER')
message,/inf, getenv('PWD')
traceback=scope_traceback()
message,/inf, traceback[0]
message,/inf, traceback[1]
get_datestamp, datestamp

if not keyword_set(debug) then debug=0
if not keyword_set(verbose) then verbose=0
if keyword_set(debug) then verbose=2

resultfile='qc_tilereport_' + datestamp + '.results'
openw, ilun_resultfile, resultfile, /get_lun

warningfile='qc_tilereport_' + datestamp + '.warning'
openw, ilun_warningfile, warningfile, /get_lun

printf, ilun_resultfile, '# ', systime(/utc)
printf, ilun_resultfile, '# ' + getenv('USER') 
printf, ilun_resultfile, '# ' + getenv('PWD') 
printf, ilun_resultfile, '# ', SCOPE_TRACEBACK()
printf, ilun_resultfile, format='(a,a)','# ',filename

printf, ilun_warningfile, '# ', systime(/utc)
printf, ilun_warningfile, '# ' + getenv('USER') 
printf, ilun_warningfile, '# ' + getenv('PWD') 
printf, ilun_warningfile, '# ',SCOPE_TRACEBACK()
printf, ilun_warningfile, format='(a,a)','# ',filename

MESSAGE,/INF,'Reading: '+ filename
t0=systime(1)
; read in the dqc data within mjd range and ra, dec range

dqc_read_data, data, filename=filename, $
 casu=casu, wsa=wsa, vsa=vsa, vircam=vircam, vhs=vhs, $
 mjdrange=mjdrange, survey=survey, /verbose

outstring='# Elpased time(secs): ' + string(systime(1)-t0)
message,/inf, outstring
printf, ilun_resultfile, outstring

outstring='# DQC file read in: ' + string(n_elements(data)) + ' records'
message,/inf, outstring
printf, ilun_resultfile, outstring

help, data
if keyword_set(debug) then help, data, /str
if keyword_set(verbose) then structure_stats, data, /sort

mjdobs=data.mjdobs
chipno=data.chipno

;isort=sort(mjdobs,/info)
;isort=bsort(mjdobs,/info)
isort=multisort(mjdobs,chipno)

data=data[isort]
ndata=n_elements(data)

outstring='# Date range: ' + $
  mjd_isodate(min(data.mjdobs)) + '  ' + $
  mjd_isodate(max(data.mjdobs))
printf, ilun_resultfile, outstring
message, /inf, outstring

outstring='# Summary info for each tile in mjdobs order '
printf, ilun_resultfile, outstring
message,/inf, outstring

outstring='# i obid nfiles obname'
printf, ilun_resultfile, outstring
message,/inf, outstring


itile=0

outstring='# Number of processed images: ' + string(ndata)
printf, ilun_resultfile, outstring
message,/inf,outstring

;for i=1,1000 do begin
;  print, i, data[i].filtname
;endfor

;pause,/force
ra1_tile=0.0
ra2_tile=0.0
ra3_tile=0.0
ra4_tile=0.0

; loop through all the qc records already sorted by mjdobs
;  to locate and count complete tiles.
; analyze sequentially by OB and group by filter.
; loop tests for end of a tile via a change of filter or OB.
; This is not a sufficeint condition since sometimes a filter in the
; same OB is repeated. 
; Need to identify some other keywords that can be tested.
;  tplexpno seems to be the one. 

; set the OB/Tile warning counter to 0
iwarning=0
for i=0, ndata-1 do begin
  npaws=0
  itile=itile+1
  filter=data[i].filtname
  obid=data[i].obsid
  filter_same=1
  obid_same=1
  npaws=1
  istart=i
  while (filter_same eq 1) and (obid_same eq 1) and (i lt ndata-1) do begin
    ; this looks for the start of a new tile/filter
    if strpos(data[i+1].filtname, filter) lt 0 then filter_same=-1
    ; this looks for the start of a new OB
    if data[i+1].obsid ne obid then obid_same=-1
    ; tile but also finds OBs that are paused and then restarted
    if data[i+1].tplexpno lt data[i].tplexpno then obid_same=-1
    if filter_same ge 0 and obid_same ge 1 then npaws=npaws+1
    if filter_same ge 0 and obid_same ge 1 then i=i+1
  endwhile
  iend=i   
  irange=iend-istart

  npawprints=npaws/16.0
  format='(i6,2i8,2i7,2x,a,2i4,4i4,2f10.3,f8.4,f8.2,2x,a)'
  outstring= string(format=format, $
   itile, data[istart].obsid, data[iend].obsid, $
   istart, iend, filter, $
   npaws, npaws/16, $
   data[istart].tplnexp, data[iend].tplnexp, $
   data[istart].tplexpno, data[iend].tplexpno, $
   data[istart].mjdobs-50000, data[iend].mjdobs-50000, $
   data[iend].mjdobs-data[istart].mjdobs, $
   (data[iend].mjdobs-data[istart].mjdobs)*(24.0*60.0), $
   data[istart].dateobs)
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    iwarning=iwarning+1
    printf, ilun_warningfile, outstring
  endif


  format='(i6,i8,2x,a,4i10,2x,a)'
  outstring=string(format=format,$
   itile, data[istart].obsid, $
   strtrim(mjd_isodate(data[isort[i]].mjdobs),2), $
   data[istart].dateobs, data[iend].dateobs, $
   data[istart].nightobs, data[iend].nightobs, $
   strtrim(data[istart].obsname,2))
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
  endif

  ; some fits data filename info
  format='(i6,i8,i4,2x,a,2x,a)'
  outstring= string(format=format, $
   itile, data[istart].obsid, $
   data[istart].programme_id, $
   data[istart].filename, $
   data[istart].esoname)
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
  endif

  ; analyze the pawprint centres specified in header: RA, Dec
  ramin=min(data[istart:iend].ra)
  ramax=max(data[istart:iend].ra)
  rarange= rarange(ramin, ramax)
  ramean=ramean(ramin, ramax, rawrap=rawrap)

  decmin=min(data[istart:iend].dec)
  decmax=max(data[istart:iend].dec)
  decrange= (decmax-decmin)
  decmean= 0.5*(decmax+decmin)

  format='(i6,i8,2x,i1,2f9.3,2f9.3,f9.3,f7.3,f9.3,f8.3,2x,a)'
  outstring=string(format=format,itile, $
   data[istart].obsid, rawrap, $
   ramin, ramax, $
   decmin, decmax, $
   rarange, decrange, $
   ramean, decmean, $
   strtrim(adstring(ramean,decmean, precision=1),2))
  print, outstring
  printf, ilun_resultfile, outstring

  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
  endif

  ; analyze the chips centres derived from the data
  cenramin=min(data[istart:iend].cenra)
  cenramax=max(data[istart:iend].cenra)

  ; compute mean taking into account 24hr wrap
  cenrarange= rarange(cenramin, cenramax)
  cenramean=ramean(cenramin, cenramax, rawrap=rawrap)

  cendecmin=min(data[istart:iend].cendec)
  cendecmax=max(data[istart:iend].cendec)
  cendecrange= (cendecmax-cendecmin)
  cendecmean= 0.5*(cendecmax+cendecmin)

  format='(i6,i8,2x,i1,2f9.3,2f9.3,f9.3,f7.3,f9.3,f8.3,2x,a)'
  outstring=string(format=format, $
   itile, $
   data[istart].obsid, $
   rawrap, $
   cenramin, cenramax, $
   cendecmin, cendecmax, $
   cenrarange, cendecrange, $
   cenramean, cendecmean, $
   strtrim(adstring(cenramean,cendecmean, precision=1),2))
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
  endif

  ; analyze the chip corners that are computed from the processed data
  pawramin=min( $
   [data[istart:iend].ra1, data[istart:iend].ra2, $   
    data[istart:iend].ra3, data[istart:iend].ra4])    
  pawramax=max( $
   [data[istart:iend].ra1, data[istart:iend].ra2, $   
    data[istart:iend].ra3, data[istart:iend].ra4])    

  pawrarange= rarange(pawramin, pawramax)
  pawramean= ramean(pawramin, pawramax, rawrap=rawrap)

  pawdecmin=min( $
   [data[istart:iend].dec1,data[istart:iend].dec2, $   
    data[istart:iend].dec3,data[istart:iend].dec4])    
  pawdecmax=max( $
   [data[istart:iend].dec1,data[istart:iend].dec2, $   
    data[istart:iend].dec3,data[istart:iend].dec4])    
  pawdecrange= (pawdecmax-pawdecmin)
  pawdecmean= 0.5*(pawdecmax+pawdecmin)

  
  ; store the tile centres 
  if itile eq 1 then begin
    ratile_cen=pawramean
    dectile_cen=pawdecmean
  endif
  if itile gt 1 then begin
    ratile_cen=[ratile_cen, pawramean]
    dectile_cen=[dectile_cen, pawdecmean]
  endif


  format='(i6,i8,2x,i1,2f9.3,2f9.3,f9.3,f7.3,f9.3,f8.3,2x,a)'
  outstring=string(format=format,itile, $
   data[istart].obsid, $
   rawrap, $
   pawramin, pawramax, $
   pawdecmin, pawdecmax, $
   pawrarange, pawdecrange, $
   pawramean, pawdecmean, $
   strtrim(adstring(pawramean,pawdecmean, precision=1),2))
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
  endif


   ; determine pointing offsets
   delta_ra=ramean-pawramean
   ; deal with 24hr wrap
   if abs(delta_ra) gt 180.0 then delta_ra=(abs(delta_ra)-360.0)
   delta_dec=decmean-pawdecmean


  ; store the offsets between SADT tile centres and observed tile centres
  if itile eq 1 then begin
    tile_delta_ra=delta_ra
    tile_delta_dec=delta_dec
  endif
  if itile gt 1 then begin
    tile_delta_ra=[tile_delta_ra, delta_ra]
    tile_delta_dec=[tile_delta_dec, delta_dec]
  endif

  format='(i6,i8,2x,i1,2f9.3,2f10.2,f9.3,2f10.2, 2x,a,3x,a)'
  outstring=string(format=format,itile, $
   data[istart].obsid, $
   rawrap, $
   ramean, decmean, $
   delta_ra*3600, delta_dec*3600.0, $
   cos(decmean/!radeg), $
   (delta_ra*3600)*cos(decmean/!radeg), delta_dec*3600.0, $
   strtrim(adstring(delta_ra, precision=1),2), $
   strtrim(adstring(delta_dec, precision=1),2))
  print, outstring
  printf, ilun_resultfile, outstring
  if npaws ne 96 then begin
    printf, ilun_warningfile, outstring
    for idata=istart, iend do begin
      if data[idata].chipno eq 1 then begin
        format='(i6,i8,i8,i4,i8,i8,i8,i8,i8)'
        printf, format=format, ilun_warningfile, $
         itile, data[istart].obsid, $
         idata, $
         data[idata].chipno, $
         data[idata].expno, $
         data[idata].obsnum, $
         data[idata].image_id, $
         data[idata].tplexpno, data[idata].tplnexp
      endif
    endfor
  endif

endfor

ntiles=itile
nwarnings=iwarning

outstring='# Total number of chips analyzed   : ' + string(ndata)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

outstring='# Total number of tiles identified : ' + string(ntiles)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

outstring='# Nominal total number of pawprints: ' + string(ndata/16.0)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

outstring='# Nominal total number of tiles    : ' + string(ndata/96.0)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

outstring='# Elpased time (secs)              : ' + string(systime(1)-t0)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

outstring='# Number of warnings               : ' + string(nwarnings)
message,/inf, outstring
printf, ilun_resultfile, outstring
printf, ilun_warningfile, outstring

; close the files
close,ilun_resultfile
close,ilun_warningfile

pause,/force

window, xsize=600, ysize=600 
!p.multi = [0,1,1]
!x.style=1  
!y.style=1
!p.background=fsc_color('white')
!p.color=fsc_color('black')

xdata=ratile_cen
ydata=dectile_cen
psym=3

print,'Measured Tile centre RA range: ',min(xdata), max(xdata)
print,'Measured Tile centre Dec range: ',min(ydata), max(ydata)

title=trim_filename(filename)
xtitle='RA'
ytitle='Declination'
xrange=[min(xdata)-2.0, max(xdata)+2.0]
yrange=[min(ydata)-2.0, max(ydata)+2.0]
charsize=1.4
plot,xdata,ydata,psym=psym, charsize=charsize,$
 title=title, xtitle=xtitle, ytitle=ytitle, $
 xrange=xrange, yrange=yrange
pause

; plot the Delta RA and Delta Dec sepaartion of tiles 
; from SADT and pipeline

; convert from degrees to arc seconds
xdata=tile_delta_ra*3600.0
ydata=tile_delta_dec*3600.0
psym=1

title=trim_filename(filename)
xtitle='Delta RA'
ytitle='Delta Dec'
xrange=[-30.0,30.0]
yrange=[-30.0,30.0]
charsize=1.4
ndata=n_elements(xdata)
itest=where( $
 xdata ge xrange[0] and xdata le xrange[1] and $
 ydata ge yrange[0] and ydata le yrange[1], count)
ndata_plot=n_elements(itest)
xmedian=median(xdata)
ymedian=median(ydata)
plot,xdata,ydata,psym=psym, $
 title=title, xtitle=xtitle, ytitle=ytitle, $
 xrange=xrange, yrange=yrange, charsize=charsize
legend='n(all):  ' + string(ndata)
legend=[legend,'n(plot): ' + string(ndata_plot)]
legend=[legend, 'median : ' + string(xmedian,'(f6.2)') + ', ' + $
 string(xmedian,'(f6.2)')]
legend,legend
pause

plothist,xdata,xrange=xrange,xtitle=xtitle
pause

plothist,ydata,xrange=yrange,xtitle=ytitle
pause


xdata=xdata-10.0
ydata=ydata-10.0
psym=1

title=trim_filename(filename)
xtitle='Delta RA'
ytitle='Delta Dec'
xrange=[-30.0,30.0]
yrange=[-30.0,30.0]
charsize=1.4
ndata=n_elements(xdata)
itest=where( $
 xdata ge xrange[0] and xdata le xrange[1] and $
 ydata ge yrange[0] and ydata le yrange[1], count)
ndata_plot=n_elements(itest)
xmedian=median(xdata)
ymedian=median(ydata)
plot,xdata,ydata,psym=psym, $
 title=title, xtitle=xtitle, ytitle=ytitle, $
 xrange=xrange, yrange=yrange, charsize=charsize
legend='n(all):  ' + string(ndata)
legend=[legend,'n(plot): ' + string(ndata_plot)]
legend=[legend, 'median : ' + string(xmedian,'(f6.2)') + ', ' + $
 string(xmedian,'(f6.2)')]
legend,legend
pause

binsize=0.2
plothist,xdata,xrange=xrange,xtitle=xtitle, bin=binsize
pause

plothist,ydata,xrange=yrange,xtitle=ytitle, bin=binsize
pause

xdata=xdata[itest]
ydata=ydata[itest]

xmad=mad(xdata)
ymad=mad(ydata)
print, median(xdata), xmad
print, median(ydata), ymad

xmad=medabsdev(xdata)
ymad=medabsdev(ydata)
print, median(xdata), xmad
print, median(ydata), ymad

print,moment(xdata,mdev=mdev,sdev=sdev),sdev,mdev
print,moment(ydata,mdev=mdev,sdev=sdev),sdev,mdev

; determine the tile spacing from the measured tile centres

;qc_tileoffsets, data=data, ra=ra, dec=dec
;pro qc_tileoffsets, data=data, ra=ra, dec=dec

ra1=ratile_cen
dec1=dectile_cen
ra2=ratile_cen
dec2=dectile_cen

SearchRadiusDegrees=2

message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
allow=100
close_match_radec,ra1,dec1,ra2,dec2,m1,m2,SearchRadiusDegrees,$
 allow,miss1,/silent
print,'close_match_radec completed'
message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
nmatches=n_elements(m1)
nmiss1=n_elements(miss1)
print,'nmatches: ',nmatches,nmiss1,n_elements(ra1),n_elements(ra2)


index=where(m1 ne m2) ; eliminate the self matches
ramatch1=ra1[m1[index]]
decmatch1=dec1[m1[index]]
ramatch2=ra2[m2[index]]
decmatch2=dec2[m2[index]]

dra=ramatch1-ramatch2
ddec=decmatch1-decmatch2
dra=dra
ddec=ddec
xrange=[-2.0,2.0]
yrange=[-2.0,2.0]
xtitle='RA sep(degrees)'
ytitle='Dec sep(degrees)'

plot,dra,ddec,psym=1,$
 xrange=xrange,yrange=yrange,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle,/isotropic
plotfile='qc_tilesep_' + datestamp + '_a.png'
pause, plotfile=plotfile

xrange=[0.0,2.0]
yrange=[0.0,2.0]
xtitle='Abs[RA sep(degrees)]'
ytitle='Abs[Dec sep(degrees)]'

plot,abs(dra),abs(ddec),psym=1,$
 xrange=xrange,yrange=yrange,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle,/isotropic
plotfile='qc_tilesep_b.png'
pause, plotfile=plotfile


xrange=[0.0,2.0]
yrange=[1.4,1.5]
xtitle='Abs[RA sep(degrees)]'
ytitle='Abs[Dec sep(degrees)]'

plot,abs(dra),abs(ddec),psym=1,$
 xrange=xrange,yrange=yrange,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle
plotfile='qc_tilesep_c.png'
pause, plotfile=plotfile


xrange=[0.0,2.0]
yrange=[1.4*60.0,1.5*60.0]
xtitle='Abs[RA sep (degrees)]'
ytitle='Abs[Dec sep (arcmins)]'

plot,abs(dra),abs(ddec)*60.0,psym=1,$
 xrange=xrange,yrange=yrange,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle
plotfile='qc_tilesep_d.png'
pause, plotfile=plotfile


!p.multi = [0,1,1]
xtitle='Frequency'
ytitle='Delta RA [degrees]'
plothist,dra,bin=0.1,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle
plotfile='qc_tilesep_dRA_hist.png'
pause, plotfile=plotfile



xtitle='Frequency'
ytitle='Delta Dec [degrees]'
plothist,dra,bin=0.1,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle
plotfile='qc_tilesep_dDec_hist.png'
pause, plotfile=plotfile



units=0
gcirc,units,ramatch1/!radeg,decmatch1/!radeg,ramatch2/!radeg,decmatch2/!radeg,dr
dr=dr*!radeg ; convert from radians to degrees

dr=sphdist(ramatch1,decmatch1,ramatch2,decmatch2,/degrees)*60.0
print,'minmax(sphdist: ',min(dr),max(dr),n_elements(dr)

xtitle='Frequency'
ytitle='Delta R [degrees]'
plothist,dr,bin=0.1,charsize=2.0,$
  title=title,xtitle=xtitle,ytitle=ytitle
plotfile='qc_tilesep_dR_hist.png'
pause, plotfile=plotfile

end

