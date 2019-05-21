;$Id: qc_compare_casu_vsa.pro,v 1.1 2012/07/21 15:59:31 rgm Exp rgm $
pro qc_compare_casu_vsa, casucat=casucat, vsacat=vsacat, $
                         verbose=verbose, tmass=tmass, xieta=xieta, $
                         waveband=waveband

;+
;
; NAME:
;       qc_compare_casu_vsa
;
; PURPOSE:
;  Compair a pair of VISTA CASU and VSA bandmerged catalogues
;
;  $Log: qc_compare_casu_vsa.pro,v $
;  Revision 1.1  2012/07/21 15:59:31  rgm
;  Initial revision
;
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2012/07/21 15:59:31 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.1 $ 
;
;-


print, 'WAVEBAND:', waveband
detection=1

compile_opt idl2

window, xsize=700, ysize=700 
!p.multi = [0,1,1]
!x.style=1  & !y.style=1
!p.background=fsc_color('white')
!p.color=fsc_color('black')


; open and read the vsa catalogue
message, /inf, 'Reading file: ' + vsacat
vsa = mrdfits(vsacat,1,hdr)
if keyword_set(verbose) then structure_info, vsa, /sort

if not keyword_set(tmass) then begin
  ; open and read the casu catalogue
  message,/inf,'Reading file: ' + casucat
  casu = mrdfits(casucat, 1, hdr)
  if keyword_set(verbose) then structure_info, casu, /sort
  print, minmax(casu.ra)
  print, minmax(casu.dec)
endif


vsa.ra=vsa.ra*!dradeg
vsa.dec=vsa.dec*!dradeg
print, minmax(vsa.ra)
print, minmax(vsa.dec)

radius=1.0
if not keyword_set(tmass) then begin

  itest1=where(casu.japercormag3 gt 0, count)
  itest2=where(vsa.japermag3 gt 0, count)

  ra1=casu[itest1].ra
  dec1=casu[itest1].dec

  ra2=vsa[itest2].ra
  dec2=vsa[itest2].dec

  allow=1

  max_radius_degrees = radius/3600.0
  close_match_radec, ra1, dec1, ra2, dec2, $
    m1, m2, max_radius_degrees, allow, miss1, silent=silent

endif

if keyword_set(tmass) then begin

  tag_test = waveband + '_X'
  itag = tag_indx(vsa, tag_test)
  print, tag_test, itag
  itest = where(abs(vsa.(itag)) lt 100000.0, count)

  if not keyword_set(detection) then begin
    splog, 'Source RA Dec: '
    ra1=vsa[itest].ra
    dec1=vsa[itest].dec
  endif

  if keyword_set(detection) then begin
    tag_test=waveband + '_RA'
    tag_test = 'RA'
    itag = tag_indx(vsa, tag_test)
    splog, 'Detection  RA: ',itag
    ra1 = vsa[itest].(itag)
    tag_test=waveband + '_dec'
    tag_test = 'DEC'
    ITAG = TAG_INDX(vsa, tag_test)
    splog, 'Detection Dec: ',itag
    dec1 = vsa[itest].(itag)
  endif

  ra2=vsa[itest].tmass_ra
  dec2=vsa[itest].tmass_dec

  ndata=n_elements(ra1)
  m1=lindgen(ndata)
  m2=m1

endif

if keyword_set(xieta) then qc_vsaxieta, data=vsa

if keyword_set(xieta) then qcplot_xieta, data=vsa

close_match_radec_offsets, $
   ra1,dec1,ra2,dec2, $
   m1,m2, $
   dra, ddec, dr, dr_gcirc, /vectorize

title=vsacat
label = 'qc_' + trim_filename(vsacat,/ext) + '_' + waveband
nearneigh_plots, dra, ddec, dr_gcirc, $
    title=title, label=label, $
    radius=radius, binsize=binsize, /summary, log=0

ra1 = ra1[m1]
dec1 = dec1[m1]
ra2 = ra2[m2]
dec2 = dec2[m2]


if keyword_set(tmass) then begin
  tag_test=waveband + '_X'
  itag = tag_indx(vsa, tag_test)
  ra1 = vsa[itest].(itag)

  tag_test=waveband + '_Y'
  itag = tag_indx(vsa, tag_test)
  dec1 = vsa[itest].(itag)
endif

ramin = min(ra1)
ramax = max(ra1)

decmin = min(dec1)
decmax = max(dec1)

rarange = ramax - ramin
decrange = decmax - decmin

nbins = 8
rabin = rarange/nbins
decbin = decrange/nbins

print, 'rabin:  ', rabin
print, 'decbin: ', decbin


xgrid=make_array(nbins, /double)
ygrid=make_array(nbins, /double)
ugrid=make_array(nbins, nbins, /double)
vgrid=make_array(nbins, nbins, /double)
print, 'Compute the systematics in a grid:', nbins, nbins
for i=0, nbins-1 do begin
  grid_ramin = ramin + (i*rabin)
  grid_ramax= grid_ramin + rabin
  xgrid[i]=(grid_ramin+grid_ramax)/2.0
  for j=0, nbins-1 do begin
    print, '2dbin: ', i, j
    grid_decmin = decmin + (j*decbin)
    grid_decmax= grid_decmin + decbin
    ygrid[j]=(grid_decmin+grid_decmax)/2.0
    itest=where(ra1 ge grid_ramin and ra1 le grid_ramax $
      and dec1 ge grid_decmin and dec1 le grid_decmax, count)
    ugrid[i,j] = median(dra[itest])
    vgrid[i,j] = median(ddec[itest])
    print, i, j, count, ugrid[i,j], vgrid[i,j], xgrid[i], ygrid[j]
  endfor
endfor

print, minmax(ugrid)
print, minmax(vgrid)


xtitle='RA (degrees)'
ytitle='Dec (degrees)'
if keyword_set(tmass) then xy=1

print, minmax(ra1)
print, minmax(dec1)

title=vsacat
plot_radec_systematics, ra1, dec1, ra2, dec2, dra, ddec, $
 title=title, waveband=waveband, xy=xy, label=label

end

