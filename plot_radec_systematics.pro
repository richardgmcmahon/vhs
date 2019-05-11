pro plot_radec_systematics, ra1, dec1, ra2, dec2, dra, ddec, $
 title=title, label=label, xy=xy, waveband=waveband


xtitle='RA (degrees)'
ytitle='Dec (degrees)'
if keyword_set(xy) then begin
  xtitle='X (pixel)'
  ytitle='Y (pixel)'
endif


nbins=8

charsize=1.4

print
print, minmax(ra1)
print, minmax(dec1)

print, minmax(ra2)
print, minmax(dec2)

print, minmax(dra)
print, minmax(ddec)

print

ramin=min(ra1)
ramax=max(ra1)

decmin=min(dec1)
decmax=max(dec1)

rarange = ramax - ramin
decrange = decmax - decmin

nbins=8
rabin=rarange/nbins
decbin=decrange/nbins

npts=nbins*nbins

xvect=make_array(npts, /double)
yvect=make_array(npts, /double)

xpos=make_array(npts,  /double)
ypos=make_array(npts,  /double)

k=-1
for i=0, nbins-1 do begin
  grid_ramin = ramin + (i*rabin)
  grid_ramax = grid_ramin + rabin

  for j=0, nbins-1 do begin
    k=k+1
    xpos[k]=(grid_ramin+grid_ramax)/2.0
    grid_decmin = decmin + (j*decbin)
    grid_decmax= grid_decmin + decbin
    ypos[k]=(grid_decmin+grid_decmax)/2.0
    itest=where(ra1 ge grid_ramin and ra1 le grid_ramax $
      and dec1 ge grid_decmin and dec1 le grid_decmax, count)
    xvect[k]=median(dra[itest])
    yvect[k]=median(ddec[itest])
    print, i, j, k, count, xpos[k], ypos[k], xvect[k], yvect[k]
  endfor

endfor

xrange=[ramin, ramax]
yrange=[decmin, decmax]

plot, xpos, ypos, /nodata, /isotropic, $
 xrange=xrange, yrange=yrange, charsize=charsize, $
 title=title, xtitle=xtitle, ytitle=ytitle
plotid, /right

;aspect=0.5
charsize=2.0
length=0.08
print, 'length = ',length
partvelvec, xvect, yvect, xpos, ypos, $
 charsize=charsize, /over, length=length, $
 aspect=aspect

rmedian=median( sqrt((xvect*xvect) + (yvect*yvect)) )
rmax=max( sqrt((xvect*xvect) + (yvect*yvect)) )

print, 'length = ',length
length=0.08*length

print, xvect[0], yvect[0], xpos[0], ypos[0]
print, xvect[0], yvect[0], xpos[0]+1000, ypos[0]+1000

;partvelvec, xvect[0], yvect[0], xpos[0]+1000, ypos[0]+1000, $
; charsize=charsize, /isotropic, /over, veccolor=fsc_color('red'), $
; aspect=aspect, length=length


;if keyword_set(tmass) then begin
  legend=""
  if keyword_set(waveband) then legend=waveband
  ndata=n_elements(ra1)  
  legend=[legend,'n = ' + string(ndata)]
  legend=[legend,'dr(median) [mas]= ' + string(rmedian*1000,'(F5.1)')]
  legend=[legend,'dr(max)    [mas]= ' + string(rmax*1000,'(F5.1)')]

  al_legend, legend, /clear, /top, /center, charsize=1.4
;endif

; plot a legend arrow
xvect= [0.10, 0.10]
yvect= [0.10, 0.10]

xvect=xvect*10000.0
yvect=yvect*10000.0

xpos=  [((ramax-ramin)/2.0)-1000.0, ((ramax-ramin)/2.0)+1000.0]
ypos=  [(decmax-(0.5*(decmax-decmin)))-1000.0, (decmax-(0.5*(decmax-decmin)))+1000.0]

print, (ramax-ramin)/2.0, decmax-(0.5*(decmax-decmin))

print
print, minmax(xpos)
print, minmax(ypos)
print, minmax(xvect)
print, minmax(yvect)

;partvelvec, xvect, yvect, xpos, ypos, /over, veccolor=fsc_color('red'), $
;  length=length

plotfile = label + '_radec_spatialsys.png'

pause, plotfile=plotfile

end
