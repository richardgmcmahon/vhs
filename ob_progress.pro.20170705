;$Id: ob_progress.pro,v 1.4 2011/11/12 12:46:58 rgm Exp rgm $
pro ob_progress, infile=infile, vsa=vsa, wfau=wfau, wsa=wsa, casu=casu, $
 verbose=verbose, debug=debug, pause=pause, $
 ps=ps, $
 incomplete=incomplete, survey=survey, runs=runs, $
 date=date, executed=executed, any=any, sdss=sdss, obstatus=obstatus, $
 batch=batch, reverse=reverse, wrfits=wrfits, $
 wise=wise, $
 ukidss=ukidss, overlap=overlap, $
 format=format, tiles=tiles, pawprints=pawprints, $
 rarange=rarange, decrange=decrange,$
 near=near, duplicity=duplicity, $
 exclude=exclude, polyfill=polyfill, title=title, $
 waveband=waveband, sgp=sgp, $
 clean_des=clean_des, $
 desfoot=desfoot, $
 desfile=desfile, $
 despolygon=despolygon,  overplot_despolygon=overplot_despolygon, $
 vstatlas=vstatlas, vstfile=vstfile, $
 plotpath=plotpath, label=label, $
 nofootprints=nofootprints, $
 publication=publication, $
 dust=dust, $
 hermes=hermes, $
 gaia=gaia, $
 spt=spt, $
 viking=viking, video=video, $
 pngsize=pngsize, $
 filtername=filtername, $
 dqc=dqc, $
 extra=extra

;+
;
;infile_vst=infile_vst
; clean_des: list the ATLAS footprint OBs with DES observations
;
; e.g. run=['A','B']
;
;ob_id status ob_name                          
;RA (hrs)      DEC (deg) 
;container_id seeing sky_transparency
;
;
; Status ’+’ (Accepted)
; Status ’-’ (Rejected)
; Status ’A’ (Aborted)
; Status ’C’ (Completed)
; Status ’D’ (Defined)
; Status 'K' (Cancelled)
; Status ’M’ (Must repeat)
; Status ’P’ (Partially defined)
; Status ’S’ (Started)
; Status 'X' 
;
; NAME:
;       ob_progress
;
; PURPOSE:
;       Plot VHS OBs
;
;
; MODIFICATION HISTORY:
;       July, 2010: Original version
;
; 
;  $Id: ob_progress.pro,v 1.4 2011/11/12 12:46:58 rgm Exp rgm $
;
;  $Log: ob_progress.pro,v $
;  Revision 1.4  2011/11/12 12:46:58  rgm
;  changed readcol format from short int to long int since ob_ids are long int
;
;  Revision 1.3  2011/04/26 21:45:19  rgm
;  added support for CASU dqc fits file
;
;  Revision 1.2  2011/03/23 16:09:12  rgm
;  first version
;
;  Revision 1.1  2010/07/29 16:00:03  rgm
;  Initial revision by Richard McMahon
;
;
;
; Login name of author of last revision:   $Author: rgm $ 
; Date and time (UTC) of revision:         $Date: 2011/11/12 12:46:58 $
; Login name of user locking the revision: $Locker: rgm $ 
; CVS revision number:                     $Revision: 1.4 $ 
;
;-

compile_opt idl2

common com_splog, loglun, fullprelog

if keyword_set(desfile) then begin
  print, 'desfile: ', desfile
  test=file_test(desfile)
  if test eq 0 then begin
    message,'desfile does not exist'
  endif
endif

despolygon=1
if keyword_set(despolygon) then begin
  splog, traceback()
  despolyfile = '/home/rgm/soft/des/round13-poly.txt'
  print, 'despolyfile: ', despolyfile
  test = file_test(despolyfile)
  print, 'test: ', test
  if test eq 0 then begin
    message,'despolygon file does not exist'
  endif
  format = '(d,d)'
  readcol, despolyfile, format=format, $
      ra_despolygon, dec_despolygon

  ra_despolygon=ra_despolygon/15.0
  print, 'RA range:  ', minmax(ra_despolygon)
  print, 'Dec range: ', minmax(dec_despolygon)
endif

pause

get_datestamp, datestamp

if not keyword_set(label) then label=""
if not keyword_set(dqc) then dqc=0
if not keyword_set(casu) then casu=0
if not keyword_set(vsa) then vsa=0
if not keyword_set(wfau) then wfau=0
if not keyword_set(filtername) then filtername=""


; read in some extra data
if keyword_set(extra) then begin
  format='(a,d,d)' 
  readcol, extra,format=format, $
   name_extra, ra_extra, dec_extra
  ra_extra=ra_extra/15.0
  print, 'RA extra range:  ', minmax(ra_extra)
  print, 'Dec extra range: ', minmax(dec_extra)
  pause, batch=batch
endif


if keyword_set(vstfile) then begin
  infile_vst=vstfile
  vsttiles=mrdfits(infile_vst,1,hr)
  structure_info, vsttiles
  pause, batch=batch
endif

if keyword_set(desfile) then begin
   ;desfile='/home/rgm/Projects/DES/Footprint/des_completed_20140207.fits'
   ;desfile='/home/rgm/Projects/DES/Footprint/des_footprint_20140207.fits'
endif

if keyword_set(desfile) then begin
    infile_des=desfile
    test=file_test(infile_des)
    if test eq 0 then infile_des='/data/des/SVA1/SVA1_COADD_GRIZY_tiles.fits'
    print, 'Reading: ', infile_des
    ipos=strpos(infile_des, 'txt')
    if ipos gt 0 then begin
      format='(d,d)' 
      readcol, infile_des, format=format, ra, dec
      struct = {ra:0d0,dec:0d0}
      destiles = replicate(struct, n_elements(ra))
      destiles.ra=ra
      destiles.dec=dec
      print, 'minmax(destiles.ra):  ', minmax(destiles.ra)
      print, 'minmax(destiles.dec): ', minmax(destiles.dec)
    endif


    ipos=strpos(infile_des, 'fits')
    if ipos gt 0 then begin
      destiles=mrdfits(infile_des,1,hr)
    endif

    itest=where(destiles.ra lt 0, count)
    if count gt 0 then destiles[itest].ra=destiles[itest].ra + 360.0
    structure_info, destiles

endif

; read in the VIKING DQC file
if keyword_set(viking) then begin

  vikpath = '/data/vhs/dqc/2012/'
  vikfile = 'vistaqc_20120630_pawprints_viking_thin.fits'

  vikpath = '/data/vhsardata/VIKING/'
  vikfile = 'VIKINGv20161202_thin_sources.fits'
  
  infile_viking = vikpath + vikfile
  viking=mrdfits(infile_viking,1,hr)
  print, minmax(viking.ra)
  print, minmax(viking.dec)
  pause, batch=batch

endif

if keyword_set(video) then begin

  vikpath='/data/vhs/dqc/2012/'
  vikfile='vistaqc_20120630_pawprints_video_thin.fits'
  infile_viking=vikpath + vikfile
  video=mrdfits(infile_viking,1,hr)
  print, minmax(video.ra)
  print, minmax(video.dec)
  pause, batch=batch

endif


if keyword_set(plotpath) then begin
  plotpath=plotpath + "/"
  test=file_test(plotpath,/directory)
  if test eq 0 then begin
      message, /inf, 'plotpath: ' + plotpath + ' does not exist'
      message, /inf, 'creating plotpath: ' + plotpath 
      file_mkdir, plotpath
  endif
endif
if not keyword_set(plotpath) then plotpath=""

des_veto=0
if keyword_set(clean_des) then des_veto=1

format='portal'
if format eq 'portal' then begin
  nruns=n_elements(runs)
  if n_elements(runs) eq 0 then runs=''
  if n_elements(survey) eq 0 then survey=''

  run_title=''
  if n_elements(runs)  gt 0 then begin
    for i=0, nruns-1 do begin
      run_title = run_title + runs[i]
    endfor
  endif
endif

splog, traceback()
splog, 'run_title: ', run_title
;pause
logfile='ob_progress' + $
 cgTimeStamp(11, RANDOM_DIGITS=2, /UTC) + '.logfile'
logger, logfile=logfile, /head, hdr=hdr
splog, filename=logfile, /append, systime(/utc)
flush, loglun

htmlfile= 'ob_progress_' + $
   '_' + datestamp + '.html'
openw, ilun_htmlfile, htmlfile, /get_lun

; ESO code
set_obstatus_codes, obstatus_codes
 
IF keyword_set(ps) THEN begin
  cmd='device,xsize=11.0,ysize=7.5,yoffset=11.0,xoffset=0.5,/inches,/landscape'
endif

IF not keyword_set(ps) THEN begin
  xsize=768
  ysize=512
  window, xsize=xsize, ysize=ysize
  if keyword_set(pngsize) then window, xsize=pngsize[0], ysize=pngsize[1]
  ;help, !x, /str
  ;help, !y, /str
  ;help, !p ,/str
  ;png_start, size=[768, 512] 
  ;help, !x, /str
  ;help, !y, /str
  ;help, !p, /str
  !p.background=fsc_color('white')
  !p.color=fsc_color('black')
  DEVICE, SET_FONT='Helvetica Bold Italic', /TT_FONT
endif
!x.style=1  & !y.style=1


if keyword_set(infile) then begin
  format='dqc'
  dqc=1
endif

if not keyword_set(infile) then begin
  infile='p85_vhs_2010-07-26.dat'
  infile='vhs_status_2010-07-26.dat'
endif

get_datestamp, datestamp
splog, 'datestamp: ', datestamp

if not keyword_set(date) then date=datestamp

;date='20110128'
;date='20110212'
;date='20110221'
;date='20110301'
path='/data/vhs/progress/'

t0=systime(1)
;read_struct, infile, struct_def, data, skipline=40

if format eq 'p2pp' then begin
  splog, traceback()
  format='(i,a,a,f,f,i,f,a)' 
  readcol,infile,format=format, $
   ob_id, ob_status, ob_name,ra,dec,container_id,seeing, sky_trans
endif

if format eq 'portal' then begin
  splog, traceback()
  inpath= path + date + '/'
  program='179A2010'

  infile= inpath + program + '.csv'
  ;if n_elements(runs) gt 0 then begin
  ;  if runs[0] eq 'A' then infile= inpath + program + 'A.csv'
  ;  if runs[0] eq 'B' then infile= inpath + program + 'B.csv'
  ;  if runs[0] eq 'C' then infile= inpath + program + 'C.csv'
  ;  if runs[0] eq 'D' then infile= inpath + program + 'D.csv'
  ;  if runs[0] eq 'BC' then infile= inpath + program + 'BC.csv'
  ;endif

  plotfile=plotpath + 'ob_progress_radec_' + date + survey + run_title + $
   '_' + label + filtername + datestamp 

  IF keyword_set(ps) THEN begin
    plotfile= plotfile + '.ps'
    colorpsopen,plotfile,cmd
  endif
  IF not keyword_set(ps) THEN plotfile = plotfile + '.png'
  
  splog, traceback()
  splog, 'Reading: ', infile

  ; pre-20110223
  ; run ID, OB ID, Execution time (s), 
  ; OB status,OB name,
  ; RA (hrs),DEC (deg),
  ; Container type, Container ID,Seeing,Sky transparency,Airmass,
  ; FLI, Moon distance (deg),Twilight (min)
  ; 179.A-2010(A), 431297, 2045 ,C ,Baseline_DES_Str1_RA060_1_1_1,
  ; 0.028893,-0.737300,G,431296,1.2,2CLR,1.7,1.0,30,0

  ; run ID,OB ID,Execution time (s),
  ; OB status,Status date,OB name,OD name,
  ; RA (hrs),DEC (deg),
  ; Container type, Container ID,Seeing,Sky transparency,Airmass,
  ; FLI, Moon distance (deg),Twilight (min)

  splog, strlen(date), n_elements(date)
  splog, traceback()
  if long(date) le long(20110222) then begin
    format='(a,i,f,a,a,f,f)'
    format='(a,l,f,a,a,f,f)'
    readcol,infile,format=format, $
     run_id, ob_id, exectime, ob_status, $
     ob_name, ra, dec, delimiter=','
    ; container_id, $
    ; seeing, sky_trans
  endif
  if long(date) ge long(20110223) then begin
    format='(a,i,f,a,a,a,a,f,f)' 
    format='(a,l,f,a,a,a,a,f,f)' 
    readcol,infile,format=format, $
     run_id, ob_id, exectime, ob_status, ob_status_date, $
     ob_name, od_name, ra, dec, delimiter=','
     ; container_id, $
     ; seeing, sky_trans
  endif
endif
if keyword_set(debug) then pause, /force

filehead=trim_filename(infile,/head)
outfile = 'ob_progress_' + filehead + '_' + datestamp + '.txt'
print, 'outfile: ', outfile
openw, ilun_sumfile, outfile, /get_lun
printf, ilun_sumfile,'# infile = ',infile
printf, ilun_sumfile, '# ',traceback()
printf, ilun_sumfile, format='(a)', '# i+1, i, m1sort, m2sort, ob_id1, ob_status1, ob_id2, obstatus2, dra, ddec, radec1, radec2, ob_name1, od_name1, ob_name2, od_name2'

print,'dqc: ', dqc
print,'vsa: ', vsa
splog, traceback()

if format eq 'dqc' then begin
  splog, traceback()
  plotfile= $
   plotpath + '/ob_progress_radec_' + filehead + label + filtername + '_' + datestamp 
  IF keyword_set(ps) THEN begin
    plotfile= plotfile + '.ps'
    colorpsopen,plotfile,cmd
  endif
  IF not keyword_set(ps) THEN plotfile = plotfile + '.png'

  wavebands_default=['Z','Y','J','H','K']
  survey_name=['ATLAS','DES','GPS']
  survey_code=['ATL','DES','GPS']

  ;survey=surveys[2]

  splog,'Reading: ',infile
  test=file_test(infile)
  if not test then begin
    splog, 'Input file does not exist: ',infile
    splog,traceback()
    message,'Exited'
  endif

  data=mrdfits(infile,1,hr)
  structure_info, data, /sort

  if keyword_set(vsa) then begin
    itag_ra=tag_indx(data,'ra')
    if itag_ra ge 0 then data.ra=data.ra*!radeg
    if itag_ra lt 0 then begin
      itag_rabase=tag_indx(data,'rabase')
      if itag_rabase ge 0 then begin
         addstruct = {ra:0d0}
         addstruct = replicate(addstruct, n_elements(data))
         newstruct = struct_addtags(data, addstruct)
         data=newstruct
         data.ra=data.rabase*15
      endif
    endif


    itag_dec=tag_indx(data,'dec')
    if itag_dec ge 0 then data.dec=data.dec*!radeg
    if itag_dec lt 0 then begin
      itag_decbase=tag_indx(data,'decbase')
      if itag_decbase ge 0 then begin
         addstruct = {dec:0d0}
         addstruct = replicate(addstruct, n_elements(data))
         newstruct = struct_addtags(data, addstruct)
         data=newstruct
         data.dec=data.decbase
      endif
    endif

  endif

  if keyword_set(pawprints) then begin
    label=label + 'pawprints_'
    itest=where(strpos(strupcase(data.is_tile), 'T') lt 0 , count)
    data=data[itest]
  endif

  if keyword_set(verbose) then structure_info, data, /sort
  ndata=n_elements(data)
  logstring='Number of data records: ' + string(ndata)
  printf, ilun_sumfile, logstring

  if not keyword_set(vsa) then begin
    itile=where(strpos(strupcase(data.is_tile), 'T') ge 0 , count)
    logstring='Number of Tiles: ' + string(count)
  endif

  if keyword_set(vsa) then count=n_elements(data)

  ntiles_all=count
  message, /inf, logstring
  printf, ilun_sumfile, logstring

  if not keyword_set(vsa) then begin
    ipawprint=where(strpos(strupcase(data.is_tile), 'T') lt 0 , count)
    logstring='Number of Pawprints: ' + string(count)
    printf, ilun_sumfile, logstring
    message, /inf,logstring
  endif

  if keyword_set(pawprints) then ntiles_all=count 


  ;itest=itile
  ;itest=ipawprint
  ;data=data[itest]

  ndata=n_elements(data)
  nruns=ndata

  if not keyword_set(vsa) then begin

  ; determine the ESO QCSTATUS codes used
  itest=UNIQ(data.qcstatus, SORT(data.qcstatus))
  QCStatusCodesUsed=data[itest].qcstatus
  nqcstatus_codes=-1
  nqcstatus_codes=n_elements(QCStatusCodesUsed)
  splog, 'Number of unique status codes used: ', nqcstatus_codes
  print, 'Status codes used: ', string(data[itest].qcstatus)
  legend_statuscode='ESO QC code'
  for i=0, nqcstatus_codes-1 do begin
    itest=where(data.qcstatus eq QCStatusCodesUsed[i], count)
    splog, 'Status codes: ', QCStatusCodesUsed[i], count
    legend_statuscode=[legend_statuscode, $
     QCStatusCodesUsed[i] + ': ' + string(count,'(i5)')]
  endfor

  endif

  if keyword_set(survey) then begin
    itest=where(strpos(data.obsname, survey) ge 0, count)
    splog,'Number in survey ',survey,': ',count
    print,'Number in survey ',survey,': ',count
    data=data[itest]
  endif

  splog, traceback()
  if not keyword_set(vsa) then begin
    ; unique filenames
    itest=UNIQ(data.filename, SORT(data.filename))
    n_filename_unique=n_elements(itest)
    splog, 'Total number of unique filenames:  ',n_filename_unique
   endif

  if keyword_set(vsa) then begin
    itest=UNIQ(data.yfilename, SORT(data.yfilename))
    n_filename_unique=n_elements(itest)
    splog, 'Total number of Y band unique filenames:  ',n_filename_unique
  endif

  ; unique OB IDS
  if not keyword_set(vsa) then begin
    itest=UNIQ(data.obsid, SORT(data.obsid))
    n_ob_unique=n_elements(itest)
    tmpdata=data[itest]
    logstring= 'Total number of unique OBs: ' + string(n_ob_unique)
    printf, ilun_sumfile, logstring
    splog, logstring

    itag_mjd=tag_indx(data,'mjd')
    if itag_mjd lt 0 then   itag_mjd=tag_indx(data,'mjdobs')
    mjd_min=min(data.(itag_mjd))  
    mjd_max=max(data.(itag_mjd))  

    message,/inf, $
     'MJD: ' + STRING(mjd_min) + STRING(mjd_max) 
    message,/inf, 'ISO date: ' + $
      MJD_ISODATE(mjd_min) + '  ' + MJD_ISODATE(mjd_max)
  endif 

  structure_info, data, /sort
  splog, traceback()
  splog, 'Number of records: ',n_elements(data)
  survey_test='DES'
  ob_select_survey, data=data, survey=survey_test, $
   ilun_sumfile=ilun_sumfile, filtername=filtername, $
   index=index, $
   n_ob=n_ob_des, ntiles_all=ntiles_des, $
   ntiles_y=ntiles_y_des, ntiles_j=ntiles_j_des, $
   ntiles_h=ntiles_h_des, ntiles_k=ntiles_k_des
  splog, traceback()
  splog, 'Number of records: ',n_elements(data)
  pause, batch=batch

  splog, traceback()
  survey_test='ATL'
  splog, 'Number of records: ',n_elements(data)
  ob_select_survey, data=data, survey=survey_test, $
   ilun_sumfile=ilun_sumfile, filtername=filtername, $
   index=index, $
   n_ob=n_ob_atlas, ntiles_all=ntiles_atlas, $
   ntiles_y=ntiles_y_atlas, ntiles_j=ntiles_j_atlas, $
   ntiles_h=ntiles_h_atlas, ntiles_k=ntiles_k_atlas
  pause, batch=batch

  splog, traceback()
  survey_test='GPS'
  splog, 'Number of records: ',n_elements(data)
  ob_select_survey, data=data, survey=survey_test, $
   ilun_sumfile=ilun_sumfile, filtername=filtername, $
   index=index, $
   n_ob=n_ob_gps, ntiles_all=ntiles_gps, $
   ntiles_y=ntiles_y_gps, ntiles_j=ntiles_j_gps, $
   ntiles_h=ntiles_h_gps, ntiles_k=ntiles_k_gps
  pause, batch=batch

  structure_info, data, /sort
  splog, traceback()
  splog, 'Total number of frames: ',ntiles_des + ntiles_atlas + ntiles_gps
  splog, 'Total number of OBs: ',n_ob_des + n_ob_atlas + n_ob_gps
  splog, 'Number of records: ',n_elements(data)
  pause, batch=batch

  ; unique filenames
  itest=UNIQ(data.filename, SORT(data.filename))
  n_filename_unique=n_elements(itest)

  splog, 'Number of unique filenames: ',n_filename_unique

  itag_filtername=tag_indx(data,'filtername')
  if itag_filtername lt 0 then itag_filtername=tag_indx(data,'filtname')


  itest=UNIQ(data.(itag_filtername), SORT(data.(itag_filtername)))
  wavebands=data[itest].(itag_filtername)
  splog, 'Wavebands: ',wavebands
  nbands=n_elements(wavebands)
  ntiles_y=0
  ntiles_j=0
  ntiles_h=0
  ntiles_k=0
  for i=0, nbands-1 do begin
    itest=where(data.(itag_filtername) eq wavebands[i], count)
    splog,i,': ',wavebands[i],': ',count
    if strpos(wavebands[i],'Y') ge 0 then ntiles_y=count
    if strpos(wavebands[i],'J') ge 0 then ntiles_j=count
    if strpos(wavebands[i],'H') ge 0 then ntiles_h=count
    if strpos(wavebands[i],'K') ge 0 then ntiles_k=count

    if strpos(wavebands[i],'u_SDSS') ge 0 then ndata_u=count
    if strpos(wavebands[i],'g_SDSS') ge 0 then ndata_g=count
    if strpos(wavebands[i],'r_SDSS') ge 0 then ndata_r=count
    if strpos(wavebands[i],'i_SDSS') ge 0 then ndata_i=count
    if strpos(wavebands[i],'z_SDSS') ge 0 then ndata_z=count

    tmp=data[itest]

    ; find the duplicate tiles based on repeat OB and waveband combos

  endfor


  ;endif


  if keyword_set(vst) then begin

  itest=UNIQ(data.prog, SORT(data.prog))
  programs=data[itest].prog
  splog, 'Programs: ',programs
  nprograms=n_elements(programs)
  for i=0, nprograms-1 do begin
    itest=where(data.prog eq programs[i], count)
    splog,i,': ',programs[i],': ',count
  endfor

  ; surveyname
  itest=UNIQ(data.surveyname, SORT(data.surveyname))
  surveynames=data[itest].surveyname
  splog, 'SurveyNames: ',surveynames
  nsurveynames=n_elements(surveynames)
  for i=0, nsurveynames-1 do begin
    itest=where(data.surveyname eq surveynames[i], count)
    splog,i,': ',surveynames[i],': ',count
  endfor

  itest=where(strpos(data.prog, '3011') ge 0, count)
  splog,'Numberof OB in survey ',survey,': ',count
  data=data[itest]

  endif

  ;pause

  splog, traceback()
endif

splog, traceback()

if format ne 'dqc' then begin
  splog, traceback()
  ; convert to a structure
  struct_def = { $
ra: 0.0d, $
dec: 0.0d, $
run_id: '', $
ob_id: 0L, $
exectime:0.0, $
ob_status:'', $
ob_status_date:'', $
ob_name:'', $
od_name:'', $
y:-1, $
j:1, $
h:-1, $
k:1 }

ndata=n_elements(ra)
data=replicate(struct_def,ndata)
data.ra=ra
data.dec=dec
data.run_id=run_id
data.ob_id=ob_id
data.ob_name=ob_name
data.exectime=exectime
data.ob_status=ob_status
  if n_elements(od_name) gt 0 then data.od_name=od_name
  if n_elements(ob_status_date) gt 0 then data.ob_status_date=ob_status_date

for i=0, ndata-1 do begin
 ; add Y and H band for ATLAS OBs
 if strpos(data[i].od_name, 'ATL') ge 0 then begin
   data[i].y=1
   data[i].h=1
 endif
 ; H band for DES OBs
 if strpos(data[i].od_name, 'DES') ge 0 then begin
   data[i].h=1
 endif
endfor

if keyword_set(sgp) then begin
  itest=where(data.ra gt 12.0, count)
  data[itest].ra = data[itest].ra - 24.0
endif

; check the DES 
ndata=n_elements(data)
splog, traceback()
splog,'DES OB check:',survey
;pause
for i=0, ndata-1 do begin
  if strpos(data[i].od_name, 'DES') ge 0 and $
   data[i].dec ge -10.0   and data[i].dec le -3.0 then begin 
   print, i, data[i].ra, data[i].dec, data[i].ob_name, data[i].od_name
  endif
endfor


if des_veto then begin
itest=where(strpos(data.od_name, 'DES') ge 0 and $
   data.dec ge -10.0   and data.dec le -3.0 and $
   data.ra le 2.0, count, complement=icomplement, ncomplement=ncomplement)
splog,'Number of DES OBS excluded: ', count, ncomplement, ndata
data=data[icomplement]
endif
;pause

endif

if not keyword_set(vsa) then begin
; exclude the cancelled OBs
if keyword_set(exclude) then begin
  ipos=strpos(data.ob_status, 'K')
  itest = where(ipos ge 0, count)
  message, /inf,'Number of Cancelled OBs: ' + string(count)
  itest = where(ipos lt 0, count)
  message, /inf,'Number of valid  OBs: ' + string(count)
  data=data[itest]
endif  
endif
splog, traceback()

if keyword_set(verbose) then structure_info, data, /sort

if keyword_set(debug) then pause, batch=batch

; determine the number of unique runs and keep the most recent
; Each OBs has one unique entry in the Portal files
if not keyword_set(dqc) then begin

if nruns gt 0 then begin
  for i=0,nruns-1 do begin
    run_test=strcompress('('+ runs[i] + ')')  
    itest=where(strpos(data.run_id, run_test) gt 0, count)
    splog,'Number of OBS in run ',runs[i], count
    if i eq 0 then itest_all=itest
    if i gt 0 then itest_all=[itest_all,itest]
  endfor
  if itest_all[0] ge 0 then begin
    splog,traceback()
    splog,'Number of OBS in runs ',n_elements(itest_all)
    data=data[itest_all]
  endif
endif
endif
splog, traceback()

if keyword_set(wrfits) then begin
  if not keyword_set(outfile) then $
   outfile='vhs_ob_status_' + datestamp + '.fits'
  splog,'Writing : ', outfile
  if keyword_set(verbose) then structure_info, data
  mwrfits, data, outfile, /create
  ;pause
endif

if keyword_set(survey) then begin
  if not keyword_set(dqc) then begin
    itest=where(strpos(data.od_name, survey) ge 0, count)
  endif
  if keyword_set(dqc) then begin
    itest=where(strpos(data.obsname, survey) ge 0, count)
  endif
  splog,'Number is survey ',survey,': ',count
  data=data[itest]
endif


if keyword_set(waveband) then begin
  itag=tag_indx(data,waveband)
  itest=where(data.(itag) gt 0, count)
  splog,'Number with band: ',waveband,': ',count
  data=data[itest]
endif

splog, traceback()

ncodes=0
print, 'dqc : ', dqc

if not keyword_set(dqc) then begin

t1=systime(1)
print, n_elements(data.ob_status)
message,/inf,'Elapsed time(secs): ' + string(t1-t0)

ndata=n_elements(data)
splog, 'Number of OB records: ', ndata
splog, filename=logfile, 'Number of OB records: ', ndata

itest=UNIQ(data.ob_id, SORT(data.ob_id))
n_unique=n_elements(itest)
splog, 'Number of unique  OBs: ', n_elements(itest)
splog, filename=logfile, 'Number of unique  OBs: ', n_elements(itest)

; unique OD_NAMES
itest=UNIQ(data.od_name, SORT(data.od_name))
ntest=n_elements(itest)
n_odname_unique=n_elements(itest)
splog, 'Total number of unique odnames: ',n_odname_unique
for i=0,ntest-1 do begin
  print, i, ': ',data[itest[i]].od_name
endfor


; unique OB_NAMES
itest=UNIQ(data.ob_name, SORT(data.ob_name))
n_obname_unique=n_elements(itest)
splog, 'Total number of unique obnames: ',n_obname_unique


; STATUS codes: + - A C M X
itest=UNIQ(data.ob_status, SORT(data.ob_status))
splog, 'Number of unique status codes used: ', n_elements(itest)
splog, 'Status codes used: ', data[itest].ob_status
StatusCodesUsed=data[itest].ob_status
ncodes=n_elements(StatusCodesUsed)

endif

if not keyword_set(rarange) then rarange=[0.0, 24.0]
xrange=rarange
if keyword_set(reverse) then xrange = [24.0, 0.0]
if keyword_set(sgp) then xrange = [-12.0, 12.0]
rarange=xrange

if keyword_set(rarange) then xrange=rarange

xtitle='RA (hours)'
ytitle='Declination (degrees)'
yrange=[-90.0, 45.0]

if keyword_set(decrange) then yrange=decrange

splog, traceback()

if not keyword_set(title) then $
  title = label + ' VHS Progress: ' + infile
charsize=1.4
plot, [0], [0], $
 charsize=charsize, $
 xrange=xrange, yrange=yrange, $
 title=title, xtitle=xtitle, ytitle=ytitle, $
 /nodata

if not keyword_set(dqc) then begin

for i=0, ncodes-1 do begin
  xdata=data.ra
  ydata=data.dec
  ipos=strpos(data.ob_status, StatusCodesUsed[i])
  itest = where(ipos gt -1, count)
  message, /inf,'Number of OBs: ' + StatusCodesUsed[i] + ' ' + string(count)  
  splog,'Number of OBs: ' + StatusCodesUsed[i] + ' ' + string(count)  
  if keyword_set(obstatus) and count gt 0 then begin
    xdata=data[itest].ra
    ydata=data[itest].dec
    psym=6
    if not keyword_set(title) then $
     title = label + ' VHS Progress: ' + infile
    if n_elements(runs)  gt 0 then title = title + '  ' + run_title
    if keyword_set(survey) then title = survey + '  ' + title
    splog, traceback()
    splog, 'run_title: ', run_title
    splog, 'title:     ', title
    ;pause, /force

    oplot, xdata, ydata, psym=psym, charsize=charsize, $
      xrange=xrange, yrange=yrange, $
      title=title, xtitle=xtitle, ytitle=ytitle
    IF not keyword_set(ps) THEN  plotid, /right
    ndata=n_elements(xdata)
    legend=StatusCodesUsed[i] + string(ndata,'(i6)')
    if not keyword_set(publication) then al_legend, legend
    ;pause, batch=batch
  endif

endfor

;pause, batch = batch

endif

splog, traceback()


if keyword_set(dqc) then begin

  splog, traceback()
  if keyword_set(debug) then structure_info, data, /sort
  xdata=data.ra/15.0
  ydata=data.dec
  psym=6
  if not keyword_set(title) then $
   title=label + ' VHS Progress: ' + infile
  if keyword_set(vst) then title='VST ATLAS Progress: ' + infile
  if keyword_set(survey) then title = survey + '  ' + title
  if n_elements(runs)  gt 0 then title = label + title + '  ' + run_title
  if keyword_set(publication) and not keyword_set(title) then title=''

  splog, traceback()
  splog, 'run_title: ', run_title
  splog, 'title:     ', title
  ;pause, /force

  xtitle='RA (hours)'
  ytitle='Declination (degrees)'
  yrange=[-90.0, 45.0]
  if keyword_set(decrange) then yrange=decrange
  charsize=1.4
  print, 'xdata range: ', min(xdata), max(xdata)
  print, 'ydata range: ', min(ydata), max(ydata)
  oplot, xdata, ydata, psym=psym, color=fsc_color('Green')
  ;charsize=charsize
  ;, $
  ;  xrange=xrange, yrange=yrange, $
  ;  title=title, xtitle=xtitle, ytitle=ytitle, /nodata

  if not keyword_set(vsa) then begin

  ipos=strpos(data.obsname, 'DES')
  itest = where(ipos gt -1, count)
  splog, 'Number of DES OBs: ',count

  xdata=data[itest].ra/15.0
  ydata=data[itest].dec
  oplot, xdata, ydata, psym=psym, color=fsc_color('blue')

  ipos=strpos(data.obsname, 'ATL')
  itest = where(ipos gt -1, count)
  tmpdata=data[itest]
  itest=where(abs(tmpdata.glat) ge 30.0, count)
  xdata=tmpdata[itest].ra/15.0
  ydata=tmpdata[itest].dec
  ;oplot, xdata, ydata, psym=psym, color=fsc_color('Dark Green')
  oplot, xdata, ydata, psym=psym, color=fsc_color('Green')

  ipos=strpos(data.obsname, 'GPS')
  itest = where(ipos gt -1, count)
  tmpdata=data[itest]
  itest=where(abs(tmpdata.glat) lt 30.0, count)
  xdata=tmpdata[itest].ra/15.0
  ydata=tmpdata[itest].dec
  oplot, xdata, ydata, psym=psym, color=fsc_color('red')

  endif

  if not keyword_set(ps) THEN  plotid, /right
  if not keyword_set(nofootprints) then $
   plot_vhs_footprints, des=des, pause=pause

  ndata=n_elements(xdata)
  if keyword_set(nqcstatus_codes) then begin
    if nqcstatus_codes gt 0 and not keyword_set(publication) then $
     al_legend, legend_statuscode, pos=[0.45,0.90], /norm, charsize=1.2 
  endif

  legend='Tile statistics ' 
  if keyword_set(pawprints) then  legend='Pawprint statistics ' 
  ;if not keyword_set(vst) and not keyword_set(dqc) then begin
  if not keyword_set(vst) and not keyword_set(vsa) then begin

    legend=[legend,'---ALL--DES-ATLAS--GPS']

    legend=[legend,'  Y : ' + $
     string(ntiles_y,'(i6)') + string(ntiles_y_des,'(i6)') + $
     string(ntiles_y_atlas,'(i6)') + string(ntiles_y_gps,'(i6)')]

    legend=[legend,'  J : ' + $
     string(ntiles_j,'(i6)') + string(ntiles_j_des,'(i6)') + $
     string(ntiles_j_atlas,'(i6)') + string(ntiles_j_gps,'(i6)')]

    legend=[legend,'  H : ' + $
     string(ntiles_h,'(i6)') + string(ntiles_h_des,'(i6)') + $
     string(ntiles_h_atlas,'(i6)') + string(ntiles_h_gps,'(i6)')]

    legend=[legend,'  Ks: ' + $
     string(ntiles_k,'(i6)') + string(ntiles_k_des,'(i6)') + $
     string(ntiles_k_atlas,'(i6)') + string(ntiles_k_gps,'(i6)')]

    legend=[legend,'  All: ' + $
     string(ntiles_all,'(i6)') + string(ntiles_des,'(i6)') + $
     string(ntiles_atlas,'(i6)') + string(ntiles_gps,'(i6)')]

  endif

  if keyword_set(vst) then begin
    legend=[legend,'  u : ' + string(ndata_u,'(i6)')]
    legend=[legend,'  g : ' + string(ndata_g,'(i6)')]
    legend=[legend,'  r : ' + string(ndata_r,'(i6)')]
    legend=[legend,'  i : ' + string(ndata_i,'(i6)')]
    legend=[legend,'  z : ' + string(ndata_z,'(i6)')]
  endif

  if not keyword_set(publication) then al_legend, legend, /clear, /left, charsize=1.2
  tile_legend = legend

  legend='OB statistics ' 
  legend=[legend,' Unique OBs : ' + string(n_ob_unique,'(i6)') ]
  legend=[legend,' DES    OBs : ' + string(n_ob_des,'(i6)') ]
  legend=[legend,' ATLAS  OBs : ' + string(n_ob_atlas,'(i6)') ]
  legend=[legend,' GPS    OBs : ' + string(n_ob_gps,'(i6)') ]
  legend=[legend,' Unique filenames: ' + string(n_filename_unique,'(i6)') ]
  if not keyword_set(publication) then al_legend, legend, /clear, /right, charsize=1.2
  ob_legend = legend

  date_legend = 'Date range: '+string(mjd_iso(mjd_min))+$
  ' : '+string(mjd_iso(mjd_max))
  if not keyword_set(publication) then $
    al_legend, date_legend, /bottom, /left, charsize=1.0

  if keyword_set(despolygon) then begin
    splog,traceback()
    plot_despolygon, ra_despolygon, dec_despolygon, $
     des_polyfill=des_polyfill
  endif

  pause, batch=batch, plotfile=plotfile

  goto, exit

endif

splog, traceback()



psym=6
psyms=[3,3,3,psym]
if not keyword_set(title) then $
 title='VHS progress: ' + infile


if n_elements(runs) gt 0 and not keyword_set(title) then begin
  if not keyword_set(title) then $
  title = title + ': '
  for i=0, nruns-1 do begin
    title = title + runs[i]
  endfor
endif

splog, traceback()
splog, 'run_title: ', run_title
splog, 'title:     ', title
;pause, /force


if keyword_set(survey) then title = survey + ' ' + title
if n_elements(runs)  gt 0 then title = title + '  ' + run_title
if keyword_set(publication) and not keyword_set(title) then title=''

xtitle='RA (hours)'
ytitle='Declination (degrees)'
yrange=[-90.0, 45.0]
if keyword_set(decrange) then yrange=decrange
charsize=1.4
xdata=xrange
ydata=yrange
plot, xdata, ydata, psym=psym, charsize=charsize, $
 xrange=xrange, yrange=yrange, $
 title=title, xtitle=xtitle, ytitle=ytitle, /nodata
IF not keyword_set(ps) THEN  plotid, /right


if keyword_set(desfile) then begin
  xdata=destiles.ra/15.0
  ydata=destiles.dec
  psym=6
  ;psym=3
  oplot, xdata, ydata, psym=psym, color=fsc_color('blue'), symsize=0.5
  ;oplot, xdata, ydata, psym=psym, color=fsc_color('red')
endif

if keyword_set(despolygon) then begin
    splog,traceback()
    plot_despolygon, ra_despolygon, dec_despolygon, $
     des_polyfill=des_polyfill
endif

if keyword_set(overplot_despolygon) then begin
  splog,traceback()
  plot_despolygon, ra_despolygon, dec_despolygon, $
   des_polyfill=des_polyfill
endif

if keyword_set(vstfile) then begin
  xdata=vsttiles.ra/15.0
  ydata=vsttiles.dec
  psym=6
  ;psym=3
  oplot, xdata, ydata, psym=psym, color=fsc_color('blue'), symsize=0.5
  ;oplot, xdata, ydata, psym=psym, color=fsc_color('red')
endif

splog, traceback()
splog, 'run_title: ', run_title
splog, 'title:     ', title
;pause, /force

xrange=rarange


if keyword_set(dust) then begin

  dustplot, /overplot, lrange=[0,360.0],brange=[-90.0,45], $
   rangestr = ['0.0', '0.1', 'A(B)'],csys='equ2000',/xterm

endif





if keyword_set(wise) then begin

  wisefile='/data/vhs/wise/wise_prelim_footprint.fits'
  splog,'Reading: ',wisefile
  wise=mrdfits(wisefile,1,hdr)
  ;structure_info, wise
  ;pause
  xdata=wise.ra/15.0
  ydata=wise.dec
  yrange=[-90.0, 45.0]
  ;title=sdssfile
  ;charsize=1.4
  ;plot, xdata, ydata, psym=psym, charsize=charsize, $
  ;  xrange=xrange, yrange=yrange, $
  ; title=title, xtitle=xtitle, ytitle=ytitle, /nodata
  ;pause
  ;oplot, [0.0,24.0], [0.0,0.0], linestyle=2
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Sea Green')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Yellow')
  psym=symcat(15)
  oplot, xdata, ydata, psym=psym, color=fsc_color('Pink')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Yellow')
  splog,traceback()
  ;pause,batch=batch

  ndata=n_elements(xdata)
  if keyword_set(polyfill) then begin
  fieldsize_width=1.56444 ; degrees
  for i=0, ndata-1 do begin
    xwidth =  fieldsize_width/(15.0*cos(ydata[i]/!radeg))
    ywidth =  fieldsize_width
    xpoly=[xdata[i]+(xwidth/2.0), xdata[i]+(xwidth/2.0), $
     xdata[i]-(xwidth/2.0),xdata[i]-(xwidth/2.0)]
    ypoly=[ydata[i]+(ywidth/2.0), ydata[i]-(ywidth/2.0), $
     ydata[i]-(ywidth/2.0),ydata[i]+(ywidth/2.0)]
    ; check for 24hr wrapping
    ;print, 'xpoly: ',i, xpoly
    ;print, 'ypoly: ',i, ypoly
    polyfill, xpoly, ypoly, col=fsc_color('pink'), noclip=0
    
    ; check for 24hr wrapping
    if xpoly[0] ge 24.0 or xpoly[2] le 0 then begin   
      splog, '24hr wrap detected xdata, ydata = ', xdata[i], ydata[i]
      splog, '24hr wrap detected xpoly = ', xpoly
      splog, '24hr wrap detected ypoly = ', ypoly
    endif

    wraptest=1
    if keyword_set(wraptest) then begin
    if xpoly[0] ge 24.0 then begin
      ;splog, '24hr wraptest'
      ;splog, '24hr wrap detected xdata, ydata = ', xdata[i], ydata[i]
      ;splog, '24hr wrap detected xpoly = ', xpoly
      ;splog, '24hr wrap detected ypoly = ', ypoly
      ;pause
      xsave=xpoly
      xpoly[0]=24.0 
      xpoly[1]=24.0 
      ;splog, '24hr wrap detected xpoly = ', xpoly
      polyfill, xpoly, ypoly, col=fsc_color('pink'), noclip=0
      ;pause
      xpoly[0]=xsave[0]-24.0 
      xpoly[1]=xsave[1]-24.0 
      xpoly[2]=0.0 
      xpoly[3]=0.0 
      ;splog, '24hr wrap detected xpoly = ', xpoly
      polyfill, xpoly, ypoly, col=fsc_color('pink'), noclip=0
      ;pause
    endif

    if xpoly[2] le 0.0 then begin
      ;splog, '24hr wrap detected xdata, ydata = ', xdata[i], ydata[i]
      ;splog, '24hr wrap detected xpoly = ', xpoly
      ;splog, '24hr wrap detected ypoly = ', ypoly
      ;pause
      xsave=xpoly
      xpoly[2]=0.0 
      xpoly[3]=0.0 
      ;splog, '24hr wrap detected xpoly = ', xpoly
      polyfill, xpoly, ypoly, col=fsc_color('pink'), noclip=0
      ;pause
      xpoly[0]=24.0 
      xpoly[1]=24.0
      xpoly[2]=xsave[2]+24 
      xpoly[3]=xsave[3]+24 
      ;splog, '24hr wrap detected xpoly = ', xpoly
      polyfill, xpoly, ypoly, col=fsc_color('pink'), noclip=0
      ;pause
    endif


    endif
    ;pause
  endfor


endif
endif

splog, traceback()




if keyword_set(overlap) then begin

  ;help,wise
  ;help,wise,/str
  nwise=n_elements(wise.ra)
  message, /inf,'Number of WISE frames: ' + string(nwise)

  ipos=strpos(data.ob_status, 'C')
  itest = where(ipos gt -1, count)
  message, /inf,'Number of completed OBs: ' + string(count)


  ra1=data[itest].ra*15d0
  dec1=data[itest].dec
  print, minmax(ra1)
  print, minmax(dec1)

  ra2=wise.ra
  dec2=wise.dec
  print, minmax(ra2)
  print, minmax(dec2)

  allow=1
  ep=1.4
  close_match_radec,ra1,dec1,ra2,dec2,m1,m2,ep,allow,miss1

  ;pause

  data=data[itest]
  data=data[m1]

endif





if keyword_set(sdss) then begin

  sdssfile='/data/rgm/boss/window_flist_radec.fits'
  splog,'Reading: ',sdssfile
  sdss=mrdfits(sdssfile,1,hdr)
  itest=where(sdss.rerun eq 301,count)  
  xdata=sdss[itest].ra/15.0
  ydata=sdss[itest].dec
  yrange=[-90.0, 45.0]
  ;title=sdssfile
  charsize=1.4
  ;plot, xdata, ydata, psym=psym, charsize=charsize, $
  ;  xrange=xrange, yrange=yrange, $
  ; title=title, xtitle=xtitle, ytitle=ytitle, /nodata
  ;oplot, [0.0,24.0], [0.0,0.0], linestyle=2
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Sea Green')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Yellow')
  oplot, xdata, ydata, psym=3, color=fsc_color('Powder Blue')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Yellow')
  splog,traceback()
  ;pause,batch=batch

endif


if keyword_set(ukidss) then begin

  ukidssfile='/data/rgm/ukidss/dr9/las/lasdr9_dqc_fs_footprint.fits'
  splog,'Reading: ',ukidssfile
  ukidss=mrdfits(ukidssfile,1,hdr)
  ;structure_info, ukidss
  xdata=ukidss.ra*!radeg/15.0
  ydata=ukidss.dec*!radeg
  yrange=[-90.0, 20.0]
  ;title=sdssfile
  charsize=1.4
  ;plot, xdata, ydata, psym=psym, charsize=charsize, $
  ;  xrange=xrange, yrange=yrange, $
  ; title=title, xtitle=xtitle, ytitle=ytitle, /nodata
  ;oplot, [0.0,24.0], [0.0,0.0], linestyle=2
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Sea Green')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Yellow')
  oplot, xdata, ydata, psym=3, color=fsc_color('Red')
  ;oplot, xdata, ydata, psym=3, color=fsc_color('Light Yellow')
  splog,traceback()
  ;pause,batch=batch
  ukidss_label='UKIDSS LAS DR9'
  ;legend, ukidss_label, /top, /center, text=fsc_color('red')
  al_legend, ukidss_label, pos=[0.5,0.93],/norm, text=fsc_color('red')
endif


if keyword_set(pause) then pause

; overplot various footprints
if not keyword_set(nofootprints) then plot_vhs_footprints, $
 vstatlas=vstatlas, hermes=hermes, /cfhtls, des=des, $
 pause=pause

if keyword_set(extra) then begin
  print, minmax(ra_extra)
  oplot, ra_extra, dec_extra, color=fsc_color('red'), psym=symcat(14)
endif

if keyword_set(pause) then pause

;if not keyword_set(dqc) then begin

ndata_all=n_elements(data.ob_status)

ipos=strpos(data.ob_status, 'K')
itest = where(ipos ge 0, count)
n_cancelled=count
message,/inf,traceback()
message, /inf,'Total number of OBs: ' + string(ndata_all)
message, /inf,'Number of Cancelled OBs: ' + string(n_cancelled)
n_valid=ndata_all  - n_cancelled

legend=        'All         OBs: ' + string(ndata_all,'(i6)')
legend=[legend,'Unique     OBs: ' + string(n_unique,'(i6)')]
legend=[legend,'Valid       OBs: ' + string(n_valid,'(i6)')]


ipos=strpos(data.ob_status, 'C')
itest = where(ipos gt -1, count)
message, /inf,'Number of completed OBs: ' + string(count)

if count gt 0 then begin
  xdata=data[itest].ra
  ydata=data[itest].dec
endif

if count eq 0 then begin
  xdata=[-999.9]
  ydata=[-999.9]
endif

if count gt 0 then ndata=n_elements(xdata)
if count eq 0 then ndata=0

legend=[legend, 'Completed OBs: ' + string(ndata,'(i6)')]

if not keyword_set(polyfill) then $
 ;oplot, xdata, ydata, psym=psym, color=fsc_color('Dark Green')
 oplot, xdata, ydata, psym=psym, color=fsc_color('Green')

if keyword_set(despolygon) then begin
    splog,traceback()
    plot_despolygon, ra_despolygon, dec_despolygon, $
     des_polyfill=des_polyfill
endif


if keyword_set(polyfill) then $
 plot_vistatile, xdata=xdata, ydata=ydata

if keyword_set(polyfill) and keyword_set(viking) then begin
 ;plot_vistatile, xdata=viking.ra, ydata=viking.dec, $
 ; color=fsc_color('red')
 oplot, viking.ra/15.0, viking.dec, $
  psym=symcat(15), color=fsc_color('red'), symsize=1.0
endif

if keyword_set(polyfill) and keyword_set(video) then begin
 ;plot_vistatile, xdata=viking.ra, ydata=viking.dec, $
 ; color=fsc_color('red')
 oplot, video.ra/15.0, video.dec, $
  psym=symcat(15), color=fsc_color('orange'), symsize=5.0
endif

pcolors=[fsc_color('black'), fsc_color('black'), fsc_color('black'), $
 fsc_color('Dark Green'),fsc_color('orange')]

if not keyword_set(executed) then begin
ipos=strpos(data.ob_status, 'C')
itest = where(ipos lt 0, count)
message, /inf,traceback()
message, /inf,'Number of incomplete OBs: ' + string(count)
n_incomplete=0
if count gt 0 then begin
  xdata=[data[itest].ra]
  ydata=[data[itest].dec]
  ;if count eq 1 then begin
  ;  xdata=[data[itest].ra]
  ;  ydata=[data[itest].dec]
  ;endif
  psym=6
  oplot, xdata, ydata, psym=psym, color=fsc_color('orange')
  n_incomplete=n_elements(xdata) - n_cancelled 
endif

legend=[legend, 'Pending    OBs: ' + string(n_incomplete,'(i6)')]
psyms=[psyms,psym]


ipos=strpos(data.ob_status, 'K')
itest = where(ipos ge 0, count)
message,/inf,traceback()
message, /inf,'Number of Cancelled OBs: ' + string(count)
print, itest[0]
n_cancelled=count
if count gt 0 then begin
  xdata=data[itest].ra
  ydata=data[itest].dec
  if count eq 1 then begin
    xdata=[data[itest].ra]
    ydata=[data[itest].dec]
  endif
  nplot=n_elements(xdata)
  message, /inf,'Number of points to be plotted: ' + string(nplot)
  psym=6
  oplot, xdata, ydata, psym=psym, color=fsc_color('red')
endif

legend=[legend, 'Cancelled  OBs: ' + string(n_cancelled,'(i6)')]
psyms=[psyms,psym]
pcolors=[pcolors, fsc_color('red')]

endif


if not keyword_set(nofootprints) then plot_vhs_footprints, des=des, $
 vstatlas=vstatlas, /cfhtls, pause=pause

;legend, legend, /clear, charsize=1.2

splog, 'n_elements(legend): ',n_elements(legend)
splog, 'n_elements(psym) : ',n_elements(psyms)
print, 'psyms: ',psyms
splog, 'n_elements(colors) : ',n_elements(pcolors)
print, 'pcolors: ',pcolors
;legend, legend, /clear, charsize=1.2
;legend, legend, psym=psyms, colors=pcolors, /clear, charsize=1.2
;legend, legend, psym=psyms, colors=pcolors, charsize=1.2
if not keyword_set(publication) then al_legend, legend, psym=psyms, colors=pcolors, charsize=1.2, /clear
;if not keyword_set(publications) then legend, legend, /clear, charsize=1.2


; Status ’+’ (Accepted)
ipos=strpos(data.ob_status, '+')
itest = where(ipos ge 0, count)
message, /inf,'Number of Accepted OBs: ' + string(count)
n_accepted=count
legend=['Accepted OBs: ' + string(n_accepted,'(i6)')]


; Status ’D’ (Defined)
ipos=strpos(data.ob_status, 'D')
itest = where(ipos ge 0, count)
message, /inf,'Number of Defined OBs: ' + string(count)
n_defined=count
legend=[legend,'Defined OBs: ' + string(n_defined,'(i6)')]


; Status ’M’ (Must repeat)
ipos=strpos(data.ob_status, 'M')
itest = where(ipos ge 0, count)
message, /inf,'Number of Must repeat OBs: ' + string(count)
n_must=count
legend=[legend,'Must repeat OBs: ' + string(n_must,'(i6)')]


ipos=strpos(data.ob_status, 'K')
itest = where(ipos ge 0, count)
message,/inf,traceback()
message, /inf,'Number of Cancelled OBs: ' + string(count)
n_cancelled=count
legend=[legend,'Kancelled OBs: ' + string(n_cancelled,'(i6)')]

ipos=strpos(data.ob_status, 'A')
itest = where(ipos ge 0, count)
message, /inf,'Number of Aborted OBs: ' + string(count)
n_cancelled=count
legend=[legend,'Aborted OBs: ' + string(n_cancelled,'(i6)')]

ipos=strpos(data.ob_status, 'X')
itest = where(ipos ge 0, count)
message, /inf,'Number of Status X OBs: ' + string(count)
n_status_x=count
legend=[legend,'Status X OBs: ' + string(n_status_x,'(i6)')]


if not keyword_set(publication) then al_legend, legend, /clear, charsize=1.2, /right_legend,/top_legend


overplot='True'
if keyword_set(overplot) then begin

  if keyword_set(desfile) then begin
    xdata=destiles.ra/15.0
    ydata=destiles.dec
    psym=6
    oplot, xdata, ydata, psym=psym, color=fsc_color('blue'), symsize=0.5
  endif

  if keyword_set(vstfile) then begin
    xdata=vsttiles.ra/15.0
    ydata=vsttiles.dec
    psym=6
    ;psym=3
    oplot, xdata, ydata, psym=psym, color=fsc_color('blue'), symsize=0.5
    ;oplot, xdata, ydata, psym=psym, color=fsc_color('red')
  endif

endif

if keyword_set(despolygon) and keyword_set(overplot_despolygon) then begin
  splog,traceback()
  polyfill, ra_despolygon, dec_despolygon, col=fsc_color('blue'), noclip=0
  ndata=n_elements(ra_despolygon)
  xdata=ra_despolygon*0.0
  for i=0, ndata-1 do begin
    xdata[i]=24.0+ra_despolygon[i]
  endfor
  polyfill, xdata, dec_despolygon, col=fsc_color('blue'), noclip=0
  splog,traceback()
endif


;png_write
splog,traceback()
pause, batch=batch, plotfile=plotfile

; check the spatial separation of OBs

ra=data.ra*15.0 ; convert for degrees
dec=data.dec
radius=1800.0
radius=60.0
radius=120.0
radius=900.0
radius=1800.0
;radius=2400.0
;radius=7200.0
message, /inf,traceback()
print, n_elements(ra)
nearneigh_self, ra, dec, radius, xy=xy, $
 m1=m1, m2=m2, nmatch=100, $
 dra=dra, ddec=ddec, dr_gcirc=dr_gcirc
print, n_elements(dra)
message, /inf,traceback()

if n_elements(dra) gt 0 then begin
  binsize=10.0

  label='ob_progress_nn_' + date + survey + run_title + $
    '_' + datestamp 

  splog, traceback()
  summary=1
  label='ob_progress_nn_' + date + survey + run_title 
  
  nearneigh_plots, dra, ddec, dr_gcirc, $
   radius=radius, binsize=binsize, label=label, $
   title=title, sidelabel=sidelabel, datestamp=datestamp, $
   nplots=nplots,summary=summary, batch=batch

endif

; reset plot region values
window, xsize=xsize, ysize=ysize

; look at the group statistics
splog,'n_elements(m1): ',n_elements(m1)
splog,'n_elements(m2): ',n_elements(m2)
nmatch=n_elements(m1)
if n_elements(dra) gt 0 then begin
  for i=0, nmatch-1 do begin
    print, i+1, i, m1[i], m2[i], dra[i], ddec[i]
  endfor
endif

close_match_group, m1, m2, x1, y1, x2, y2, xy=xy, $
 m1sort, m2sort, m1start, nmembers1, diag=diag, self=self, $
 batch=batch

printf, ilun_sumfile, '# ',traceback()
for i=0, nmatch-1 do begin
  dra_test=(data[m1sort[i]].ra*15.0 - data[m2sort[i]].ra*15.0) * $
    (cos(data[m1sort[i]].dec/!dradeg))*3600.0
  ddec_test=(data[m1sort[i]].dec - data[m2sort[i]].dec)*3600.0

  format='(i6,i6,i6,i6,i8,1x,a,i8,1x,a,1x,2f8.1,1x,a,1x,' $
   + 'a,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a)'
  string=string(format=format, $
   i+1, i, m1sort[i], m2sort[i], $
   data[m1sort[i]].ob_id, data[m1sort[i]].ob_status, $ 
   data[m2sort[i]].ob_id, data[m2sort[i]].ob_status, $
   dra_test, ddec_test, $
   adstring(data[m1sort[i]].ra*15.0, data[m1sort[i]].dec, 0), $
   adstring(data[m2sort[i]].ra*15.0, data[m2sort[i]].dec, 0), $
   data[m1sort[i]].ob_name, data[m1sort[i]].od_name, $ 
   data[m2sort[i]].ob_name, data[m2sort[i]].od_name)

   print, string

   printf, ilun_sumfile,string
   splog, /noname, string

endfor
printf, ilun_sumfile, '# ',traceback()

splog,'n_elements(m1): ',n_elements(m1)
splog,'n_elements(m2): ',n_elements(m2)
splog,'n_elements(m1sort): ',n_elements(m1sort)
splog,'n_elements(m2sort): ',n_elements(m2sort)
splog,'n_elements(m1start): ',n_elements(m1start)
splog,'n_elements(nmembers1): ',n_elements(nmembers1)
xdata=nmembers1

plotfile='ob_progress_groups_' + date + survey + run_title + $
   '_' + label + datestamp + filtername + '.png'

print, min(xdata), max(xdata) 
if max(xdata)-max(xdata) then begin
  plothist, xdata, charsize=1.4, $
   title=title, xtitle='Number of group members', ytitle='Frequency'
  plotid, /right
  pause, batch=batch, plotfile=plotfile
endif

ra_save=data.ra

itest=UNIQ(data.exectime, SORT(data.exectime))
n_exectimes=n_elements(itest)
print, '# Unique Execution times '
for i=0, n_exectimes - 1 do begin
    icount=where(data.exectime eq data[itest[i]].exectime, count)
    print, i-0, data[itest[i]].exectime, count
endfor

exectime_save=data.exectime

xdata=data.exectime
xtitle='OB Execution Time (seconds)'
ytitle='Frequency'
xmin=min(xdata)
xmax=max(xdata)
xspan=xmax-xmin
xrange=[xmin-(xspan/10.0), xmax+(xspan/10)]
print, 'xrange: ', xrange
if xspan gt 0.1 then begin
  plothist, xdata, charsize=1.4, $
    title=title, xtitle=xtitle, ytitle=ytitle, $
    xrange=xrange
  plotid, /right

  plotfile= plotpath + $
   'ob_progress_exectime_' + date + survey + run_title + '.png'
  pause, plotfile=plotfile, batch=batch
endif

ipos=strpos(data.ob_status, 'C')
itest = where(ipos gt -1, count)
message, /inf,'Number of completed OBs: ' + string(count)
if count gt 0 then begin
  xdata=data[itest].ra
  ydata=data[itest].exectime
  yrange=[0.8*min(ydata), 1.2*max(ydata)]
endif

if count eq 0 then begin
  xdata=[-999.9]
  ydata=[-999.9]
  yrange=[0.8*min(data.exectime), 1.2*max(data.exectime)]
endif

plotfile= plotpath + 'ob_progress_ra_exectime_' + date + survey + run_title 
IF keyword_set(ps) THEN begin
    plotfile=plotfile + '.ps'
    colorpsopen,plotfile,cmd
endif
IF not keyword_set(ps) THEN plotfile = plotfile + '.png'

xrange=rarange
xtitle='Right Ascension (hours)'
ytitle='Execution time (seconds)'
plot, xdata, ydata,psym=psym, charsize=charsize, $
 title=title, xtitle=xtitle, ytitle=ytitle, $
 xrange=xrange, yrange=yrange, /nodata
plotid, /right

; DES
if n_elements(data.od_name) gt 0 then begin
  splog, 'Total number of OBs: ',n_elements(data.od_name)
  ipos=strpos(data.od_name, 'DES')
  itest = where(ipos gt -1, count)
  n_des=count
  splog, 'Number of DES   OBs: ',n_des
  if count gt 0 then begin
    xdata=data[itest].ra
    ydata=data[itest].exectime
    oplot, xdata, ydata,psym=psym, color=fsc_color('Dark Green')
  endif

; GPS
  ipos=strpos(data.od_name, 'GPS')
  itest = where(ipos gt -1, count)
  n_gps=count
  splog, 'Number of GPS   OBs: ',n_gps
  if count gt 0 then begin
    xdata=data[itest].ra
    ydata=data[itest].exectime
    oplot, xdata, ydata,psym=psym, color=fsc_color('red')
  endif

; ATLAS
  ipos=strpos(data.od_name, 'AT')
  itest = where(ipos gt -1, count)
  n_atlas=count
  splog, 'Number of ATLAS OBs: ',n_atlas
  if count gt 0 then begin
    xdata=data[itest].ra
    ydata=data[itest].exectime
    oplot, xdata, ydata,psym=psym, color=fsc_color('orange')
  endif

  splog, 'Total number of OBs: ',n_des+n_atlas+n_gps

endif


;xdata=data[itest].ra
;ydata=data[itest].exectime
;oplot, xdata, ydata,psym=psym, color=fsc_color('orange')

splog,traceback()
plotfile= plotpath + 'ob_progress_ra_exectime_' + date + survey + run_title 
IF keyword_set(ps) THEN begin
  plotfile=plotfile + '.ps'
  colorpsopen,plotfile,cmd
endif
IF not keyword_set(ps) THEN plotfile = plotfile + '.png'
pause, batch=batch, plotfile=plotfile



; loop through and compute execution time by RA

ob_total=make_array(24, /long)
exectime_total=make_array(24, /float)
for i=0,23 do begin
   itest=where(fix(data.ra) eq i, count)
   print, i, count
   itest=where(data.ra ge i and data.ra lt i+1, count)
   print, i, count
   ob_total[i]=count
   if count gt 0 then exectime_total[i]=total(data[itest].exectime)
endfor

xdata=fix(data.ra)

plothist, xdata, xhist, yhist, /noplot, halfbin=0
yrange=[0.0,max(yhist)*1.1]
yrange=[0.0,250.0]
title='OB status: ' + title
xtitle='RA'
ytitle='Number of OBs per RA hour'
plothist, xdata, charsize=charsize, halfbin=0, $
 xrange=xrange, yrange=yrange, $
 title=title, xtitle=xtitle, ytitle=ytitle
plotid,/right
ndata=n_elements(xdata)
legend='Total OBs: ' + string(ndata,'(i5)')

ipos=strpos(data.ob_status, 'C')
itest = where(ipos gt -1, count)
message, /inf,'Total number of OBs     : ' + string(ndata)
message, /inf,'Number of completed  OBs: ' + string(count)
if count gt 0 then begin
  xdata=fix(data[itest].ra)
endif

if count eq 0 then begin
  xdata=[-999.9,-990.9]
endif


plothist, xdata, color=fsc_color('Dark Green'), halfbin=0, /overplot, thick=2
ndata=count
legend=[legend, 'Completed  OBs: ' + string(ndata,'(i5)')]


ipos=strpos(data.ob_status, 'C')
itest = where(ipos lt 0, count)
message, /inf,'Number of incomplete OBs: ' + string(count)
if count gt 1 then begin
  xdata=fix(data[itest].ra)
  plothist, xdata, color=fsc_color('red'), halfbin=0, /overplot
  ndata=n_elements(xdata)
endif
legend=[legend, 'Incomplete OBs: ' + string(count,'(i5)')]

pcolors=[fsc_color('black'),fsc_color('Dark Green'), fsc_color('red')]

splog, 'n_elements(legend): ',n_elements(legend)
splog, 'n_elements(psym) : ',n_elements(psyms)
splog, 'n_elements(colors) : ',n_elements(pcolors)
linestyles=[0,0,0]
al_legend, legend, linestyle=linestyles, colors=pcolors, /clear, charsize=1.2

splog,traceback()
plotfile= plotpath + 'ob_progress_ra_hist_' + date + survey + run_title 
IF keyword_set(ps) THEN begin
  plotfile=plotfile + '.ps'
  colorpsopen,plotfile,cmd
endif
IF not keyword_set(ps) THEN plotfile = plotfile + '.png'
pause, batch=batch, plotfile=plotfile

; analysis  of execution times
exectime_total=make_array(24, /float)

splog, 'Total allocated time (hours): ',total(data.exectime)/3600.0
splog, 'RA range: ',minmax(data.ra)
count_all=0
for i=0,23 do begin
   itest=where(fix(data.ra) eq i, count)
   print, i, count
   itest=where(data.ra ge i and data.ra lt i+1, count)
   print, i, count
   if count gt 0 then exectime_total[i]=total(data[itest].exectime)
   count_all=count_all + count
endfor
splog, traceback()
ndata=n_elements(data.ra)
splog, 'Number of OBs: ', ndata, count_all
splog, 'Total allocated time (hours): ',total(exectime_total)/3600.0


xdata=findgen(24)+0.5
print, xdata
ydata=exectime_total/3600.0

yrange=[0.0,max(ydata)*1.1]
;yrange=[0.0,100.0]
ytitle='Total execution time per hour of RA'
plot, xdata, ydata, charsize=charsize, $
 xrange=xrange, yrange=yrange, $
 title=title, xtitle=xtitle, ytitle=ytitle
plotid, /right
legend='Total allocated time (hrs):   ' + string(total(ydata),'(f6.1)')

splog, traceback()
ipos=strpos(data.ob_status, 'C')
idata = where(ipos gt -1, count)
message, /inf,'Number of completed OBs: ' + string(count)
if count gt 0 then message, /inf,'Total execution time: ' + $
 string(total(data[idata].exectime)/3600.0)
if count lt 0 then message, /inf,'Total execution time: ' + '0.0'

exectime_total=make_array(24, /float,value=0.0)
if count eq 0 then ydata=0
if count gt 0 then begin
  count_all=0
  for i=0,23 do begin
    itest=where(fix(data[idata].ra) eq i, count)
    print, i, count
    itest=where(data[idata].ra ge i and data[idata].ra lt i+1, count)
    print, i, count
    if count gt 0 then exectime_total[i]=total(data[idata[itest]].exectime)
    count_all=count_all+count
  endfor
  splog, traceback()
  splog,'Number of OB counted: ',count_all
  ydata=exectime_total/3600.0
  oplot, xdata, ydata, color=fsc_color('Dark Green'), thick=2
endif

legend=$
 [legend,'Total executed time (hrs):   ' + string(total(ydata),'(f6.1)')]
splog, 'Total execution time (hours): ',total(exectime_total)/3600.0

ipos=strpos(data.ob_status, 'C')
idata = where(ipos lt 0, count)
message, /inf,'Number of incomplete OBs: ' + string(count)
if count gt 0 then message, /inf,'Total execution time: ' $
 + string(total(data[idata].exectime)/3600.0)

exectime_total=make_array(24, /float,value=0.0)
ydata=exectime_total
count_all=0
if count gt 0 then begin
  for i=0,23 do begin
   itest=where(fix(data[idata].ra) eq i, count)
   print, i, count
    itest=where(data[idata].ra ge i and data[idata].ra lt i+1, count)
    print, i, count
    if count gt 0 then $
     exectime_total[i]=total(data[idata[itest]].exectime)
    count_all=count_all + count
  endfor
  splog, traceback()
  splog,'Number of OB counted: ',count_all
  ydata=exectime_total/3600.0
  oplot, xdata, ydata, color=fsc_color('red')

endif
legend=[legend,'Total unexecuted time (hrs): ' + string(total(ydata),'(f6.1)')]

al_legend, legend, $
 linestyle=linestyles, colors=pcolors, /clear, charsize=1.2

splog,traceback()
plotfile= plotpath + 'ob_progress_exectime_hist_' + date + survey + run_title 

IF keyword_set(ps) THEN begin
  plotfile=plotfile + '.ps'
  colorpsopen,plotfile,cmd
endif
IF not keyword_set(ps) THEN plotfile = plotfile + '.png'
pause, batch=batch, plotfile=plotfile

print, 'Total unexecuted time (hours): ',total(exectime_total)/3600.0


exit:

;if format eq 'dqc' then 

close,ilun_sumfile
close,ilun_htmlfile

end

