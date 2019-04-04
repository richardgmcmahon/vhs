; $Id: template.pro,v 1.1 1997/04/01 19:41:45 doug Exp $
pro plot_vhs_footprints, pause=pause, $
 cfhtls=cfhtls, $
 vstatlas=vstatlas, $
 hermes=hermes, gaia=gaia, $
 dust=dust, des=des, $
 spt=spt
;+
;
; Plot VHS footprint and other related surveys
;
;
; MODIFICATION HISTORY:
;
; 
;  $Id$
;  
;  $Source$
;
;  $Log$
;
;
; Login name of author of last revision:   $Author$ 
; Date and time (UTC) of revision:         $Date$
; Login name of user locking the revision: $Locker$ 
; CVS revision number:                     $Revision$ 
;-
;-


; rgm added Round82 option in 20120904
if not keyword_set(des) then des='Round82'
;desyr1=2

thick=2

; plot some galactic latitude lines
oplot, [0.0,24.0], [0.0,0.0], linestyle=2
plot_galcoords, /oplot, b=30, color=fsc_color('red'), $
 linestyle=2, thick=thick
plot_galcoords, /oplot, b=5, color=fsc_color('red'), $
 linestyle=0, thick=thick

plot_galcoords, /oplot, b=0, color=fsc_color('black'), $
 linestyle=2, thick=1

plot_galcoords, /oplot, b=-5, color=fsc_color('red'), $
 linestyle=0, thick=thick

plot_galcoords, /oplot, b=-30, color=fsc_color('red'), $
 linestyle=2, thick=thick

; plot VVV bulge region
; l= 350 to 10, latitude -10 to +5
plot_galcoords, /oplot, b=-10, glrange=[350,360], $
 color=fsc_color('red'), linestyle=0
plot_galcoords, /oplot, b=-10, glrange=[0,10], $
 color=fsc_color('red'), linestyle=0

plot_galcoords, /oplot, b=5, glrange=[350,360], $
 color=fsc_color('red'), linestyle=0
plot_galcoords, /oplot, b=5, glrange=[0,10], $
 color=fsc_color('red'), linestyle=0


; longitude boundaries
plot_galcoords, /oplot, glong=10, gbrange=[0,5], $
 color=fsc_color('red'), linestyle=0
plot_galcoords, /oplot, glong=10, gbrange=[-10,0], $
 color=fsc_color('red'), linestyle=0

plot_galcoords, /oplot, glong=350, gbrange=[0,5], $
 color=fsc_color('red'), linestyle=0
plot_galcoords, /oplot, glong=350, gbrange=[-10,0], $
 color=fsc_color('red'), linestyle=0

thick=1

; plot the UKIDSS UDS and DXS
if keyword_set(cfhtls) then begin

 ; UKIDSS UDS
 ; 02 18 -05 00
 ; 02 16 -> 02 20; -05 30 -04 30 
  xdata=[2.267, 2.267, 2.333, 2.333, 2.267]
  ydata=[-5.5 ,-4.5  ,-4.5    , -5.5 , -5.5]
  oplot, xdata,ydata,color=fsc_color('red')

endif


if keyword_set(cfhtls) then begin

; cfhls wide fields
; w1
; 02h 18m 00s DEC:-07.d 00m 00s
; 8deg x 9deg
; 02 02 -> 02 38 
; -2.5 -> -11.5
  xdata=[02.033, 02.033, 02.633, 02.633, 02.033]
  ydata=[-11.5 ,-2.5  ,-2.5    , -11.5 , -11.5]
  oplot, xdata,ydata,color=fsc_color('blue')

; cfhls wide fields
; D1
; 02h 26m 00s DEC:-04 30m 00s
; 1deg x 1deg
; 02 24 -> 02 28
; -5.00 -> -40
  xdata=[2.4, 2.4, 2.467, 2.467, 2.4]
  ydata=[-5.0 ,-4.0  ,-4.0    , -5.0 , -5.0]
  oplot, xdata,ydata,color=fsc_color('blue')


; W2 - 7deg x 7deg
; 08:54:00 -04:15:00 2000
; 08 40 -> 09 08
; -0.75 -> -7.75
xdata=[08.666, 08.666, 09.133, 09.133,08.666]
ydata=[-7.75 ,-0.75  ,-0.75  ,-7.75  ,-7.75]
oplot, xdata,ydata,color=fsc_color('blue')

; W4 - 4deg x 4deg
; 22:13:18 +01:19:00 2000
; 21 57 - 22 25
; -1.5 - +5
xdata=[21.95, 21.95, 22.4166, 22.4166,21.95]
ydata=[-1.5 , 5.0  , 5.0  ,-1.5  ,-1.5]
oplot, xdata,ydata,color=fsc_color('blue')

endif


if keyword_set(spt) then begin
  ; SPT
  ; 22h00 > RA > 2h00, -30 < dec < -25
  xdata=[ 20.0, 20.00, 24.0, 24.0,20.0]
  ydata=[-65.0, -45.0, -45.0, -65.0,-65.0]
  oplot, xdata,ydata,color=fsc_color('Yellow'), $
    thick=2

  xdata=[ 0.0, 0.00, 7.0, 7.0,00.0]
  ydata=[-65.0, -45.0, -45.0, -65.0,-65.0]
  oplot, xdata,ydata,color=fsc_color('Yellow'), $
   thick=2
endif

; the DES footprint has been a moving target
; with the connecting region and VIKING region
; undergoing tweaks
;des_2012_old=1
if keyword_set(des) and des eq 'des_2012_old' then begin

  ; DES
  ; Stripe82
  color=fsc_color('Purple')
  xdata=[0.0,4.00,4.00,0.0]
  ;xdata=[0.0,3.67,3.67,0.0]
  ydata=[-1.0,-1.0,+1.0,+1.0]
  oplot, xdata,ydata,color=color, $
   thick=1

  xdata=[ 24.0, 20.67, 20.67,24.0]
  ydata=[-1.0,-1.0,+1.0,+1.0]
  oplot, xdata,ydata,color=color, $
   thick=1

  ; Connecting strip
  ; 02h00 < RA < 4h00, -30 < dec < -1 changed to 0340 on XX
  ; 02h00 < RA < 3h40, -30 < dec < -1 
  xdata=[ 2.0, 2.0]
  ydata=[-1.0, -25.0]
  oplot, xdata,ydata,color=color, $
   thick=2

  xdata=[3.67, 3.67]
  xdata=[4.00, 4.00]
  ydata=[-1.0, -30.0]
  oplot, xdata,ydata,color=color, $
   thick=2


  ; Cap
  ; 22h00 > RA > 2h00, -30 < dec < -25
  xdata=[ 22.0, 22.00, 24.0]
  ydata=[-30.0, -25.0, -25.0]
  oplot, xdata,ydata,color=color, $
   thick=2

  xdata=[ 00.0, 2.00]
  ydata=[-25.0, -25.0]
   oplot, xdata,ydata,color=color, $
   thick=2


  ; SPT
  ; 20h00 < RA < 07h00, -65< dec < -30 
  xdata=[3.67,7.00,7.00,0.0]
  xdata=[4.00,7.00,7.00,0.0]
  ydata=[-30.0,-30.0,-65.0,-65.0]
  oplot, xdata,ydata,color=color, $
   thick=2

  xdata=[ 22.0, 20.00, 20.00,24.0]
  ydata=[-30.0, -30.0, -65.0,-65.0]
  oplot, xdata,ydata,color=color, $
   thick=2

  xyouts, 7.0, -80.0, 'DES footprint: circa 2011)', $
   charsize=1.4, color=color



endif


if des eq 'Round82' then begin

  ; DES-Stripe82-2012
  ;xdata=[0.0,4.00,4.00,0.0]
  ;xdata=[0.0,3.67,3.67,0.0]
  ;ydata=[-1.0,-1.0,+1.0,+1.0]
  ;oplot, xdata,ydata,color=fsc_color('Blue'), $
  ; thick=2

  xdata=[ ((360.0-3.0)/15.0), ((360.0-43.0)/15.0), ((360.0-43.0)/15.0), ((360.0-3.0)/15.0)]
  ydata=[-1.0,-1.0,+1.0,+1.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1


  ; DES-Round82-2012
  xdata=[0.0, 45.00/15.0, 45.00/15.0]
  ydata=[3.0, 3.0, -25.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1

  xdata=[ 24.0, ((360.0-3.0)/15.0), ((360.0-3.0)/15.0)]
  ydata=[3.0, 3.0, -25.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1

  ; DES-VIKING-2012 
  linestyle=0
  xdata=[45.0/15.0, 60.00/15.0, 60.00/15]
  ydata=[-25.0, -25.0, -40.0 ]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1, linestyle=linestyle

  xdata=[ (360.0-3.0)/15.0, (360.0-30.00)/15.0, (360.0-30.00)/15.0]
  ydata=[-25.0            , -25.0             , -40.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1, linestyle=linestyle

  ; DES-SPT-2012 
  ; 20h00 < RA < 07h00, -65< dec < -40
  ; note change from -65< dec < -30 to -65< dec < -40
  linestyle=0
  xdata=[4.0,7.00,7.00,0.0]
  ydata=[-40.0,-40.0,-65.0,-65.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1, linestyle=linestyle

  xdata=[ 22.0, 20.00, 20.00,24.0]
  ydata=[-40.0, -40.0, -65.0,-65.0]
  oplot, xdata,ydata,color=fsc_color('Blue'), $
   thick=1, linestyle=linestyle

  xyouts, 7.0, -80.0, 'DES: Round82 (April 29, 2012)', $
   charsize=1.4, color=fsc_color('Blue')

endif


;desyr1=1
if keyword_set(desyr1) then begin
; DES Season 1
thick=1

; VVDS Deep
; 30 < RA < 45; -10 < dec < +5
xdata=[ 14.0/15.0, 40.0/15.0, 40.0/15.0, 14.0/15.0, 14.0/15.0]
ydata=[     -7.0,     -7.0,      +3.0,      +3.0,     -7.0] 
oplot, xdata,ydata,color=fsc_color('Blue'), $
 thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=-45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0


; MaxVis
; 75 < RA < 100; -60 < dec < -49 [v1]
; 60 < RA < 90; -60 < dec < -49  [accepted]
xdata=[ 60.0/15.0, 90.0/15.0, 90.0/15.0, 60.0/15.0, 60.0/15.0]
ydata=[     -60.0,     -60.0,     -49.0,     -49.0,     -60.0] 
oplot, xdata,ydata,color=fsc_color('Blue'), $
 thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=-45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0


; SPT Deep
; -17 < RA <  2, -60 < dec < -49
xdata=[ 24.0, (360.0-17.0)/15.0, (360.0-17.0)/15.0,  24.0]
ydata=[-60.0,             -60.0,             -49.0, -49.0] 
oplot, xdata,ydata,color=fsc_color('Blue'), $
 thick=thick, linestyle=2
POLYFILL, xdata, ydata, /line_fill, orientation=45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=-45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0


xdata= [   0.0, 2.0/15.0, 2.0/15.0, 0.0]
ydata= [ -60.0,    -60.0,      -49.0,     -49.0] 
oplot, xdata,ydata,color=fsc_color('Blue'), $
 thick=thick, linestyle=2
POLYFILL, xdata, ydata, /line_fill, orientation=45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0
POLYFILL, xdata, ydata, /line_fill, orientation=-45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0



; SPT Broad 20 < RA <  105, -65 < dec < -45
xdata=[ 20.0/15.0, 105.0/15.0, 105.0/15.0, 20.0/15.0, 20.0/15.0]
ydata=[     -65.0,      -65.0,      -45.0,     -45.0,  -65.0] 
oplot, xdata,ydata,color=fsc_color('Blue'), $
 thick=thick, linestyle=2
POLYFILL, xdata, ydata, /line_fill, orientation=-45.0, $
 color=fsc_color('Blue'), thick=thick, linestyle=0

endif




if keyword_set(pause) then pause

if keyword_set(dust) then begin
  message, /inf, 'Plot reddening'
  dustplot, /overplot, lrange=[0,360.0],brange=[-90.0,45], $
   rangestr = ['0.1', '1.0', 'A(g_s)'],csys='equ2000',/xterm
endif


if keyword_set(pause) then pause


if keyword_set(vstatlas) then begin
  color=fsc_color('Purple')
  color=fsc_color('Red')

  ; P88
  ; 10h -> RA -> 15h30, -2 -> Dec -> -11
  ; 21h30 -> RA -> ~02h, -40 -> Dec -> -35
  ; 21h30 -> RA -> ~0h, -14 -> Dec -> -9
  xdata=[ 10.0,  15.5,  15.5,  10.0,  10.0]
  ydata=[ -2.0 , -2.0, -11.0, -11.0,  -2.0]
  oplot, xdata,ydata,color=color, thick=2

  xdata=[ 0.0,   4.0,   4.0,  0.0,  0.0]
  ydata=[ -9.0,  -9.0, -14.0, -14.0,  -9.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2

  xdata=[  0.0,   4.0,  4.0,   0.0,   0.0]
  ydata=[-25.0 ,-25.0, -35.0, -35.0, -25.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2
  xdata=[ 21.5,  24.0,  24.0,  21.5,  21.5]
  ydata=[-25.0 ,-25.0, -35.0, -35.0, -25.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2

  ;xdata=[ 21.5,  24.0,  24.0,  21.5,  21.5]
  ;ydata=[ -9.0,  -9.0, -14.0, -14.0,  -9.0]
  ;oplot, xdata,ydata,color=fsc_color('purple')


  ; P89 
  ; 10h -> RA -> 15h30, -20 -> Dec -> -11
  ; 21h30 -> RA -> ~02h, -40 -> Dec -> -35
  ; 21h30 -> RA -> ~0h, -14 -> Dec -> -9
  xdata=[ 10.0,  15.5,  15.5,  10.0,  10.0]
  ydata=[-11.0 ,-11.0, -20.0, -20.0, -11.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2

  xdata=[  0.0,   2.0,  2.0,   0.0,   0.0]
  ydata=[-35.0 ,-35.0, -40.0, -40.0, -35.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2
  xdata=[ 21.5,  24.0,  24.0,  21.5,  21.5]
  ydata=[-35.0 ,-35.0, -40.0, -40.0, -35.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2
  xdata=[ 21.5,  24.0,  24.0,  21.5,  21.5]
  ydata=[ -9.0,  -9.0, -14.0, -14.0,  -9.0]
  oplot, xdata,ydata,color=fsc_color('purple'), thick=2

endif




if keyword_set(hermes) then begin
; Hermes ADFS
; 04:59:26  -52:55:44
; 04:53:07  -51:45:21
; 04:25:43  -54:21:12
; 04:32:55  -55:52:18
; 04:59:26  -52:55:44
xdata=[ 74.85/15.0, 73.28/15.0, 66.43/15.0, 68.22/15.0, 74.86/15.0]
ydata=[     -52.92,     -51.75,     -54.35,     -55.87,  -52.92] 
oplot, xdata,ydata,color=fsc_color('Red'), $
 thick=thick, linestyle=2
POLYFILL, xdata, ydata,  color=fsc_color('Red')
endif


if keyword_set(pause) then pause

; VIKING

; NGC
; 10h00 < RA < 15h30 , −5 < dec < +4 
xdata=[10.0,15.50,15.50,10.0,10.0]
ydata=[+4.0,+4.0,-5.0,-5.0,+4.0]
oplot, xdata,ydata,color=fsc_color('red')

; GAMA 09
xdata=[8.6 , 9.4, 9.4, 8.6, 8.6]
ydata=[+3.0,+3.0,-2.0,-2.0,+3.0]
oplot, xdata,ydata,color=fsc_color('red')

; SGC
; 22h00 < RA < 03h30 , −36 < dec < −26 
xdata=[0.0,3.50,3.50,0.0]
ydata=[-26.0,-26.0,-36.0,-36.0]
oplot, xdata,ydata,color=fsc_color('red')

xdata=[ 24.0, 22.00, 22.00,24.0]
ydata=[-26.0, -26.0, -36.0,-36.0]
oplot, xdata,ydata,color=fsc_color('red')


if keyword_set(gaia) then begin
; GAIA spec
plot_galcoords, /oplot, b=-20, glrange=[310,360], $
 color=fsc_color('red'), linestyle=1
plot_galcoords, /oplot, b=-20, glrange=[0,50], $
 color=fsc_color('red'), linestyle=1

plot_galcoords, /oplot, b=20, glrange=[310,360], $
 color=fsc_color('red'), linestyle=1
plot_galcoords, /oplot, b=20, glrange=[0,50], $
 color=fsc_color('red'), linestyle=1

plot_galcoords, /oplot, b=20, glrange=[220,270], $
 color=fsc_color('red'), linestyle=1
plot_galcoords, /oplot, b=20, glrange=[220,270], $
 color=fsc_color('red'), linestyle=1

plot_galcoords, /oplot, b=-20, glrange=[220,270], $
 color=fsc_color('red'), linestyle=1
plot_galcoords, /oplot, b=-20, glrange=[220,270], $
 color=fsc_color('red'), linestyle=1

endif

end
