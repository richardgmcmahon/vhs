pro vhs_icrf_analysis, rfc=rfc, $
 ulas=ulas, waveband=waveband, $
 dr8=dr8, dr7qso=dr7qso, $
 viking=viking, vmc=vmc, $
 zone=zone, ukidssXiXn=ukidssXiXn, detection=detection, $
 publication=publication, psplot=psplot
;+
;
; added support for UKIDSS internal error analysis
;
;
;  publication=True: plot without annotation
;
;
;-

COMPILE_OPT idl2

window,0,xsize=600,ysize=600
!p.background=fsc_color('white')
!p.color=fsc_color('black')
!x.style=1
!y.style=1


;Open postscript device if necessary
if keyword_set(psplot) then begin
  psname='idl.ps'
  colorpsopen,psname,'device,xsize=11,ysize=5.17,'+$
  'yoffset=11.5,xoffset=2.0,/inches,/landscape'
endif

inpath='/data/vhs/icrf/'
filename='icrs_1_vhs.fits'

if keyword_set(rfc) then filename = 'rfc_2012b_vhs.fits'

if keyword_set(viking) then filename = 'rfc_2012b_viking.fits'
if keyword_set(vmc) then filename = 'rfc_2012b_vmcdr1.fits'

if keyword_set(dr7) then filename = 'rfc_2012b_sdssdr7.csv'
if keyword_set(dr8) then filename = 'rfc_2012b_sdssdr8.csv'

if keyword_set(dr7qso) then filename = 'dr7qsos_2mass_vizier.csv'

if keyword_set(ulas) then begin
  if not keyword_set(waveband) then filename='rfc_2012b_ukidssdr9_lasyjhk.fits'
  if keyword_set(waveband) then $
   filename='rfc_2012b_ukidssdr9_' + waveband + '.fits'
endif

if keyword_set(ukidssXiXn) then begin
  infile=ukidssXiXn
  filename=trim_filename(infile)
endif

if not keyword_set(ukidssXiXn) then infile= inpath + filename

test=file_test(infile)
if test eq 0 then begin
 message, infile + ' does no exist'
endif
message, /inf, 'Reading: ' + infile

if keyword_set(dr7) or keyword_set(dr8) then begin
; name, up_name, up_ra, up_dec, 
; objID, ra, dec, 
; run, rerun, camcol, field, type,
; modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z
  format='l, l, d, d, LL, d, d, l, l, l, l, a, f,f,f,f,f'
  readcol, infile, $
   name, up_name, up_ra, up_dec, $
   objID, ra, dec, run, rerun, camcol, field, type, $
   modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z, $
   format=format, skipline=1
  print, name[0], up_name[0]
  print, up_ra[0], up_dec[0]
  print, objid[0]
  print, ra[0], dec[0]
  print, run[0], rerun[0], camcol[0], field[0]
  print, type[0]

  ndata=n_elements(name)
  struct_def = { $
   id: 0L, $
   upload_ra: 0.0d, $
   upload_dec: 0.0d, $
   objid: 0LL, $
   ra: 0.0d, $
   dec: 0.0d, $
   distance: 0.0, $
   type: ''}
  data=replicate(struct_def,ndata)
  data.upload_ra=up_ra
  data.upload_dec=up_dec
  data.ra=ra
  data.dec=dec
  data.id=name
  data.objid=objid
  data.type=type

  ndata=n_elements(data)
  dra=make_array(ndata,/double)
  ddec=make_array(ndata,/double)
  dr=make_array(ndata,/double)
  dr_gcirc=make_array(ndata,/double)

  for i=0, ndata-1 do begin

    dra[i] = cos(data[i].dec/!dradeg) * 3600d0 * $
     (data[i].upload_ra - (data[i].ra))
    ddec[i]= 3600d0*(data[i].upload_dec - (data[i].dec))

    ra1_tmp = data[i].upload_ra
    dec1_tmp = data[i].upload_dec
    ra2_tmp = data[i].ra
    dec2_tmp = data[i].dec
    gcirc,2,ra1_tmp,dec1_tmp,ra2_tmp,dec2_tmp,dr_gcirc_tmp
    dr_gcirc[i]=dr_gcirc_tmp

  endfor

  print, minmax(dra)
  print, minmax(ddec)
  print, minmax(dr_gcirc)

  itest=where(abs(dra) lt 60 and abs(ddec) lt 60, count)
  print, 'count: ',count
  data=data[itest]
  dra=dra[itest]
  ddec=ddec[itest]
  dr_gcirc=dr_gcirc[itest]

  print, minmax(dra)
  print, minmax(ddec)
  print, minmax(dr_gcirc)

  data.distance=dr_gcirc

  structure_info, data

  pause

endif

if keyword_set(dr7qso) then begin
;d_arcsec,2MASS,RAJ2000,DEJ2000,
;errHalfMaj,errHalfMin,errPosAng,
;Jmag,Hmag,Kmag,e_Jmag,e_Hmag,e_Kmag,Qfl,Rfl,X,MeasureJD,
;_RAJ2000,_DEJ2000,SDSS,z,umag,gmag,rmag,imag,zmag,logX,DR7,RAJ2000,DEJ2000
  format='f, a, d, d, f, f, f, f, f, f, f,f,f, a, a, a, d,d,d'
  readcol, infile, $
   dr, name, ra, dec, $
   dum1, dum2, dum3, $
   dum4, dum5, dum6, $
   dum7, dum8, dum9, $
   dum10, dum11, dum12, $
   mjd_obs,ra_sdss, dec_sdss, $
   format=format, skipline=1
  print, name[0]
  print, ra[0], dec[0]
  print, ra_sdss[0], dec_sdss[0]

  ndata=n_elements(name)
  struct_def = { $
   id: 0L, $
   distance: 0.0, $
   ra: 0.0d, $
   dec: 0.0d, $
   upload_ra: 0.0d, $
   upload_dec: 0.0d}
  data=replicate(struct_def,ndata)
  data.distance=dr
  data.ra=ra
  data.dec=dec
  data.upload_ra=ra_sdss
  data.upload_dec=dec_sdss
  data.id=name
endif

if not keyword_set(dr7) and not keyword_set(dr8) $ 
 and not keyword_set(dr7qso) $
 and not keyword_set(ukidssXiXn) then begin

  data=mrdfits(infile,1,hdr)
  structure_info, data

  data.ra=data.ra*!dradeg
  data.dec=data.dec*!dradeg
  tag_test='framesetid'
  if keyword_set(waveband) then tag_test='multiframeid'
  itag=tag_indx(data, tag_test)
  message, /inf, 'tag, itag: ' + tag_test + ': ' + string(itag)

  itest=where(data.(itag) gt 0, count)
  data=data[itest]
  structure_info, data
  if keyword_set(zone) then begin
    if zone eq 1 then begin
      itest=where(data.ra ge 120 and data.ra le 240) 
      data=data[itest]
    endif
    if zone eq 2 then begin
      itest=where(data.ra le 120 or data.ra ge 240) 
      data=data[itest]
    endif
    if zone eq 3 then begin
      itest=where(data.dec le 45.0 and data.ra le 180)
      data=data[itest]
    endif
    if zone eq 4 then begin
      itest=where(data.dec le 45.0 and data.ra ge 180)
      data=data[itest]
    endif
  endif
endif


if keyword_set(ukidssXiXn) then begin
  data=mrdfits(infile,1,hdr)
  structure_info, data
endif

if not keyword_set(ukidssXiXn) then xdata=data.distance

if keyword_set(ukidssXiXn) then begin

  itest=where( $
   data.mergedclass eq -1 and $
   data.yapermag3 gt -1000.0 and data.j_1apermag3 gt -1000.0 and $
   data.hapermag3 gt -1000.0 and data.kapermag3 gt -1000.0, count)
  print, 'Number of objects detected in YJHK: ', count


  dra=data[itest].Kxi
  ddec=data[itest].Keta

  ydata=data[itest].kapermag3

  ; convert from pixels to arcsecs assuming no interleaving
  if keyword_set(detection) then begin
    dra=data[itest].Kerr_X*0.40
    ddec=data[itest].Kerr_Y*0.40
  endif

  ;xdata=data[itest].j_1xi
  ;dra=data[itest].J_1xi
  ;ddec=data[itest].J_1eta

  ;endif

  xdata=data[itest].Kxi
  print, minmax(xdata)
  itest=where(xdata gt -1000.0, count)
  xdata=xdata[itest]
  binsize=0.01
  plothist, xdata, bin=binsize, title='Kxi histogram'
  pause

  xdata=sqrt(dra^2 + ddec^2)

  xrange=[0, 1.0]
  yrange=[20.0, 10.0]
  print, minmax(xdata)

  if keyword_set(publication) then title=filename

  ytitle='Kapermag3'
  xtitle='sqrt(Kxi^2 + Kxn^2) arc secs (Source)'
  if keyword_set(detection) then xtitle='sqrt(dX^2 + dY^2) Pixels (Detection)'
  if keyword_set(detection) then xtitle='sqrt(dX^2 + dY^2) arc secs (Detection)'

  plot, xdata, ydata, psym=3, charsize=1.4, $
   xrange=xrange, yrange=yrange, $
   xtitle=xtitle, ytitle=ytitle, title=title
  ndata=n_elements(xdata)
  median=median(xdata)
  print, 'median = ',median
  legend='n = ' + string(ndata)
  legend=[legend, 'median = ' + string(median,'(F8.3)')]
  sigma=medabsdev(xdata,/sigma)
  legend=[legend, 'sigma(MAD) = ' + string(sigma,'(F8.3)')]

  al_legend, legend, /right, charsize=1.4
  plotid, /right

  pause

endif

if not keyword_set(dr7) and not keyword_set(dr8) and $
  not keyword_set(ukidssXiXn) then begin

tag_test='ksapermag3'
if keyword_set(ulas) then tag_test='kapermag3'

if keyword_set(waveband) then tag_test=waveband + 'apermag3'
if keyword_set(waveband) then tag_test='apermag3'

itag=tag_indx(data,tag_test)
message, /inf, 'tag, itag: ' + tag_test + ': ' + string(itag)
if itag ge 0 then begin
  itest=where(data.(itag) gt 10, count)
  binsize=0.1
  xdata=data[itest].(itag)
  plothist, xdata, bin=binsize
  pause

  xdata=data[itest].distance
  binsize=0.05
  plothist, xdata, bin=binsize
  pause
endif

ndata=n_elements(data)
dra=make_array(ndata,/double)
ddec=make_array(ndata,/double)
dr=make_array(ndata,/double)
dr_gcirc=make_array(ndata,/double)

for i=0, ndata-1 do begin

  dra[i] = cos(data[i].dec/!dradeg) * 3600d0 * $
   (data[i].upload_ra - data[i].ra)
  ddec[i]= 3600d0*(data[i].upload_dec - data[i].dec)

  ra1_tmp = data[i].upload_ra
  dec1_tmp = data[i].upload_dec
  ra2_tmp = data[i].ra
  dec2_tmp = data[i].dec
  gcirc,2,ra1_tmp,dec1_tmp,ra2_tmp,dec2_tmp,dr_gcirc_tmp
  dr_gcirc[i]=dr_gcirc_tmp

endfor

endif

binsize=0.05
if keyword_set(dr7) or keyword_set(dr8) then binsize=0.02

summary=1
radius=1.0
title=filename
charsize=1.2
psymsize=0.7
psym=8
plotsym, 0, psymsize, /fill
if keyword_set(ukidssXiXn) then psym=3
ndata=n_elements(dra)
if ndata gt 100000 then binsize=0.01
if publication then charsize=1.7
nearneigh_plots, dra, ddec, dr_gcirc, $
 radius=radius, binsize=binsize, label=label, $
 title=title, sidelabel=sidelabel, datestamp=datestamp, $
 nplots=nplots,summary=summary, batch=batch, $
 charsize=charsize, psym=psym, publication=publication, psplot=psplot





end
