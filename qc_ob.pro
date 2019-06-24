pro qc_ob, filename=filename, obid=obid, $
 mjdrange=mjdrange, $
 casu=casu, vircam=vircam, vhs=vhs, $
 report=report, data=data, $
 verbose=verbose, debug=debug, batch=batch
;+
;$Id: qc_paw.pro,v 1.1 2009/12/12 22:08:04 rgm Exp rgm $
;
; NAME:
;  QC_OB
;
; PURPOSE:
;  Quality control analysis of a single OB
;
; TODO
;  Create a summery report file and also an exception report file
;  which identifies data that has QC problems.
;
;  also write:
;    qc_ob_list module
;    qc_ob_driver module
;
;
; CALL SEQUENCE
;
;
; INPUTS
;  filename: FITS format file with QC data from 
;
; OUTPUTS
; 
;
; EXAMPLE:
;
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2009/12/12 22:08:04 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.1 $ 
;
;-

splog, traceback(/verbose)
splog, 'batch: ', batch

if not keyword_set(debug) then debug=0
if not keyword_set(verbose) then verbose=0

MESSAGE,/INF,'Reading: '+ filename
; read in the dqc data within mjd range and ra, dec range
t0=systime(1)

dqc_read_data, data, filename=filename, obid=obid, $
 casu=casu, wsa=wsa, vsa=vsa, vircam=vircam, vhs=vhs, $
 mjdrange=mjdrange, verbose=verbose

message,/inf,'Elpased time(secs): ' + string(systime(1)-t0)
print,'DQC file read in: ',n_elements(data),' rqecords'
if keyword_set(verbose) then begin
  structure_info, data, /sort
  help, data
endif
if keyword_set(debug) then help, data, /str

itest=where(strpos(data.is_tile, 'T') lt 0, count)
data=data[itest]

isort=sort(data.filename)

data=data[isort]

nfiles=n_elements(data.(0))
runnum=long(strmid(data.filename,10,5))
tmpstr={runnum:0L}
newstr=replicate(tmpstr,nfiles)
newstr.runnum=runnum
data=struct_addtags(data, newstr)
print, n_elements(runnum)
print, runnum
; VCAM.2009-11-04T00:29:14.7693.fits
; v20091103_00093_st.fit
; nightobs 20091103

for i=0, nfiles-1 do begin
  format='(i6,i4,2x,a,2x,a,a,a,4x,f12.5,2x,a,2x,i6)'

  if tag_indx(data,'mjdobs') gt 0 then $
   print, format=format, i+1, data[i].chipno, $
   strtrim(data[i].filtname,2), $
   strtrim(data[i].filepath,2), '/', strtrim(data[i].filename,2), $
   data[i].mjdobs, mjd_isodate(data[i].mjdobs), runnum[i]

  if tag_indx(data,'mjd') gt 0 then $
   print, format=format, i+1, data[i].chipno, $
   strtrim(data[i].filtname,2), $
   strtrim(data[i].filepath,2), '/', strtrim(data[i].filename,2), $
   data[i].mjd, mjd_isodate(data[i].mjd), runnum[i]

endfor


j=0
for i=0, nfiles-1 do begin
  if data[i].chipno eq 1 then begin
    j=j+1
    format='(i6,i6,i4,2x,a,2x,i4)'
    print, format=format, $
     i+1,j, data[i].chipno, $
     strtrim(data[i].filtname,2), $
     data[i].tplexpno
  endif
endfor

; DES tplexpno 1, 3, 5, 7, 9


window, xsize=900, ysize=700 
!p.multi = [0,1,5]
!x.style=1  
!y.style=1
!p.background=fsc_color('white')
!p.color=fsc_color('black')

; using tagnames to define the columns

ytag='sec_seeing'
if tag_indx(data,ytag) lt 0 then begin
  ytag='seeing'
endif

qc_ob_plot_mjd_tag, data, ytag, obid
;plotid, /right

ytag='ellipticity'
qc_ob_plot_mjd_tag, data, ytag, obid


ytag='airmass'
if tag_indx(data,ytag) lt 0 then begin
  ;ytag=['amstart','amend']
  ytag=['amstart']
endif
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='magzpt'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='maglimit'
if tag_indx(data,ytag) lt 0 then ytag='maglim'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

plotfile ='dqc_ob_' + strtrim(string(obid),2) + '_a.png'
pause, plotfile=plotfile, batch=batch

ytag='az'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='alt'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='winddir'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='windsp'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

plotfile ='dqc_ob_' + strtrim(string(obid),2) + '_b.png'
pause, plotfile=plotfile, batch=batch

ytag='apcor3'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='tplexpno'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='tplnexp'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='obsnum'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2


plotfile ='dqc_ob_' + strtrim(string(obid),2) + '_c.png'
pause, plotfile=plotfile ,/force, batch=batch

ytag='amend'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='runnum'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

ytag='expno'
qc_ob_plot_mjd_tag, data, ytag, obid, psym=2

plotfile ='dqc_ob_' + strtrim(string(obid),2) + '_d.png'
pause, plotfile=plotfile ,/force, batch=batch

!p.multi = [0,2,2]

xdata=data.ra
ydata=data.dec

xmin=min(xdata)
xmax=max(xdata)
xrange=(xmax-xmin)
message, /inf, 'xrange: ' + string(xrange)
if xrange gt 180.0 then begin
  ndata=n_elements(xdata)
  message, /inf, 'Need to wrap RA '
  idata=where(xdata gt 180.0, count)
  xdata[idata]=xdata[idata]-360d0
endif
xmin=min(xdata)
xmax=max(xdata)
xrange=abs(xmax-xmin)
xmin = xmin - (0.1*xrange)
xmax = xmax + (0.1*xrange)
xrange = [xmin, xmax]

ymin=min(ydata)
ymax=max(ydata)
yrange=(ymax-ymin)
ymin = ymin - (0.1*yrange)
ymax = ymax + (0.1*yrange)
yrange = [ymin, ymax]

print, xrange, yrange
title='Pawprint centres (RA, Dec)'
xtitle='Right Ascension'
ytitle='Declination'
plot, xdata, ydata, psym=1, charsize=1.4, $
 xrange=xrange, yrange=yrange, /isotropic, $
 title=title, xtitle=xtitle, ytitle=ytitle
pause

xtag='cen_ra'
if tag_indx(data,xtag) ge 0 then xdata=data.cen_ra
if tag_indx(data,xtag) lt 0 then xdata=data.cenra

ytag='cen_dec'
if tag_indx(data,ytag) ge 0 then xdata=data.cen_dec
if tag_indx(data,ytag) lt 0 then xdata=data.cendec

xmin=min(xdata)
xmax=max(xdata)
xrange=(xmax-xmin)
message, /inf, 'xrange: ' + string(xrange)
if xrange gt 180.0 then begin
  ndata=n_elements(xdata)
  message, /inf, 'Need to wrap RA '
  idata=where(xdata gt 180.0, count)
  xdata[idata]=xdata[idata]-360d0
endif

xmin=min(xdata)
xmax=max(xdata)
xrange=abs(xmax-xmin)
xmin = xmin - (0.20*xrange)
xmax = xmax + (0.20*xrange)
xrange = [xmin, xmax]

ymin=min(ydata)
ymax=max(ydata)
yrange=(ymax-ymin)
ymin = ymin - (0.20*yrange)
ymax = ymax + (0.20*yrange)
yrange = [ymin, ymax]

print, xrange, yrange
title='Chip centres (cen_ra, cen_dec)'
plot, xdata, ydata, psym=1, charsize=1.4, $
 xrange=xrange, yrange=yrange, /isotropic, $
  title=title, xtitle=xtitle, ytitle=ytitle

xdata=data.ra1
xdata=[xdata, data.ra2]
xdata=[xdata, data.ra3]
xdata=[xdata, data.ra4]

xmin=min(xdata)
xmax=max(xdata)
range=(xmax-xmin)
message, /inf, 'range: ' + string(range)
if range gt 180.0 then begin
  ndata=n_elements(xdata)
  message, /inf, 'Need to wrap RA '
  idata=where(xdata gt 180.0, count)
  xdata[idata]=xdata[idata]-360d0
endif

ydata=data.dec1
ydata=[ydata, data.dec2]
ydata=[ydata, data.dec3]
ydata=[ydata, data.dec4]

title='Chip corners (ra[1-4], dec[1-4])'
plot, xdata, ydata, psym=1, charsize=1.4, $
 xrange=xrange, yrange=yrange, /isotropic, $
  title=title, xtitle=xtitle, ytitle=ytitle


xdata=data.crval1
ydata=data.crval2

xmin=min(xdata)
xmax=max(xdata)
xrange=abs(xmax-xmin)
xmin = xmin - (0.1*xrange)
xmax = xmax + (0.1*xrange)
xrange = [xmin, xmax]

ymin=min(ydata)
ymax=max(ydata)
yrange=(ymax-ymin)
ymin = ymin - (0.1*yrange)
ymax = ymax + (0.1*yrange)
yrange = [ymin, ymax]

print, xrange, yrange
title='CRVAL'
plot, xdata, ydata, psym=1, charsize=1.4, $
 xrange=xrange, yrange=yrange, /isotropic, $
  title=title, xtitle=xtitle, ytitle=ytitle

plotfile ='dqc_ob_' + strtrim(string(obid),2) + '_c.png'
pause, plotfile=plotfile, batch=batch


; determine the number of wavebands and which bands
filtname=data.filtname
ibands=uniq(filtname,sort(filtname))
nbands=n_elements(ibands)
message, /inf, 'Number of wavebands: ' + string(nbands,'(i6)')
wavebands=filtname[ibands[0]]
for i=1, nbands-1 do begin
 wavebands=[wavebands,filtname[ibands[i]]]
endfor
print, 'Wavebands: ',wavebands

; cycle through each filter and write out the ra and dec limits
wavebands=['J','H','K']
ramin=9999.9
ramax=-9999.9
decmin=9999.9
decmax=-9999.9
for iband=0, 2 do begin

  

endfor



;make_astr, ast, cd=cd


;qc_ob_wcs
; analyze the WCS info
;
; use make_astr to create an WCS structre that can be used with

; CD1_1, CD1_2
; CD2_1, CD2_2
; CRPIX1, CRPIX2
; CRVAL1, CRVAL2
; CTYPE1, CTYPE2
; PV2_1, PV2_2, PV2_3, PV2_4, PV2_5
; NAXIS1, NAXIS2

i=0
debug=0
ndata=n_elements(data.(0))
ndata=5
for i=0, ndata-1 do begin

print
message, /inf,'CENCOORDS: ' + strtrim(data[i].cencoords,2)
print, 'cenra, cendec: ', data[i].cenra, data[i].cendec, $
strtrim(adstring(data[i].cenra,data[i].cendec, precision=2),2)
print, 'coords: ', data[i].coords
print, 'ra, dec: ', data[i].ra, data[i].dec, $
 strtrim(adstring(data[i].ra,data[i].dec, precision=2),2)
print, 'RA1, Dec1:  ' ,data[i].ra1, data[i].dec1,'   ', $
 strtrim(adstring(data[i].ra1,data[i].dec1, precision=2),2)
print, 'RA2, Dec2:  ' ,data[i].ra2, data[i].dec2,'   ', $
 strtrim(adstring(data[i].ra2,data[i].dec2, precision=2),2)
print, 'RA3, Dec3:  ' ,data[i].ra3, data[i].dec3,'   ', $
 strtrim(adstring(data[i].ra3,data[i].dec3, precision=2),2)
print, 'RA4, Dec4:  ' ,data[i].ra4, data[i].dec4,'   ', $
 strtrim(adstring(data[i].ra4,data[i].dec4, precision=2),2)


cd1_1=data[i].cd1_1
cd1_2=data[i].cd1_2
cd2_1=data[i].cd2_1
cd2_2=data[i].cd2_2

crpix1=data[i].crpix1
crpix2=data[i].crpix2

crval1=data[i].crval1
crval2=data[i].crval2

ctype1=data[i].ctype1
ctype2=data[i].ctype2

pv2_1=data[i].pv2_1
pv2_2=data[i].pv2_2
pv2_3=data[i].pv2_3
pv2_4=data[i].pv2_4
pv2_5=data[i].pv2_5


cd = [[cd1_1, cd2_1], [cd1_2, cd2_2]]

cd = [[cd1_1, cd2_1], [cd1_2, cd2_2]]

;print, 'cd: ',cd

crpix=[crpix1, crpix2]
crval=[crval1, crval2]
ctype=[ctype1, ctype2]

pv2 = [0.0, pv2_1, pv2_2, pv2_3, pv2_4, pv2_5]

make_astr,astr, CD=cd, CRPIX = crpix, CRVAL = crval, $
                    CTYPE = ctype, $
                    PV2 = pv2
if keyword_set(debug) then help, astr, /str

naxis1=data[i].naxis1
naxis2=data[i].naxis2

x=0
y=0
xy2ad_rgm, x, y, astr, a, d, $
     xsi=xsi, eta=eta, rtheta=rtheta, phi=phi, xx=xx, yy=yy
PRINT, format='(2F10.3,2F10.5,1X,A)', x, y, a, d, $
 strtrim(adstring(a,d,precision=2),2)

x=naxis1-1
y=0
xy2ad_rgm, x, y, astr, a, d, $
     xsi=xsi, eta=eta, rtheta=rtheta, phi=phi, xx=xx, yy=yy
PRINT, format='(2F10.3,2F10.5,1X,A)', x, y, a, d, $
 strtrim(adstring(a,d,precision=2),2)


x=naxis1-1
y=naxis2-1
xy2ad_rgm, x, y, astr, a, d, $
     xsi=xsi, eta=eta, rtheta=rtheta, phi=phi, xx=xx, yy=yy
PRINT, format='(2F10.3,2F10.5,1X,A)', x, y, a, d, $
 strtrim(adstring(a,d,precision=2),2)

x=0
y=naxis2-1
xy2ad_rgm, x, y, astr, a, d, $
     xsi=xsi, eta=eta, rtheta=rtheta, phi=phi, xx=xx, yy=yy
PRINT, format='(2F10.3,2F10.5,1X,A)', x, y, a, d, $
 strtrim(adstring(a,d,precision=2),2)



endfor

end


