SELECT
/* 
   get DQC parameters for VHS framesets 

   TODO:

   could add ISNULL NULLIF to the maglimit calculation

   Notes
   Depth calculation uses abs() for some values that may be negative; 
   would be better to use ISNULL/NULLIF see above

   down stream checking will be needed

   see also separate query which gets the DQC parameters for
   each multiframe, multiframedetector. 
--
-- $Id: vhs_vsa_dqc_tiles_fs_metadata_sql92.sql,v 1.2 2014/04/26 17:09:49 rgm Exp rgm $	
-- $Log: vhs_vsa_dqc_tiles_fs_metadata_sql92.sql,v $
-- Revision 1.2  2014/04/26 17:09:49  rgm
-- added abs() to the depth calculation to deal with log of negative values
--
-- Revision 1.1  2014/04/26 16:58:36  rgm
-- Initial revision
--
-- Revision 1.1  2014/04/26 15:58:59  rgm
-- Initial revision
--
-- Revision 1.3  2012/09/21 15:46:00  rgm
-- Added join with EsoKeys to get OB id and execution time
--
-- Revision 1.2  2012/07/07 15:03:04  rgm
-- *** empty log message ***
--
-- Revision 1.1  2012/07/07 15:00:21  rgm
-- Initial revision
--
--
-- Login name of author of last revision:   $Author: rgm $ 
-- Date and time (UTC) of revision:         $Date: 2014/04/26 17:09:49 $
-- Login name of user locking the revision: $Locker: rgm $ 
-- CVS revision number:                     $Revision: 1.2 $ 
--

MultiFrame: a single image in any band with 4 chips
MultiFrameDetector: a chip in a MultiFrame
FrameSet: is a set of bandmerged MultiFrameDetectors

Note in a frameset the MFDs do not have EXACTLY the same
footprint since sometimes a different guidestar is used and 
the regions may have been retiled.

Extreme examples are shown here:

http://www.ast.cam.ac.uk/~rgm//ukidss/framesets/


The Astrometric limits for a MFD are in another table called
CuurentAstrometry.

Example queries based on UKIDSS LAS

select count(*) from lasMergeLog

Select count(*) from lasMergeLog WHERE ymfId gt 0

SELECT Count(*) FROM lasMergeLog 
WHERE 
  ymfId > 0 AND
  J_1mfId > 0 AND 
  HmfId > 0 AND
  KmfId > 0 


Notes:
   VSA/WSA transmutes some of the CASU keywords. In the case
   of WCS keywords eg CD1_1 to cd11, ctype1 to ctypeX. 
   Project and Object seem to uniquely define the multiframes within a Tile.

History  
   20160906 added ESO Object keyword
   20060826 added wcs information from table CurrentAstrometry
   20080701 DR4 version; removed some ZP columns since they do not exist

 $Source: /Users/rgm/Projects/vhs/vsa/footprint/RCS/vhs_vsa_dqc_tiles_fs_metadata_sql92.sql,v $
 $Id: vhs_vsa_dqc_tiles_fs_metadata_sql92.sql,v 1.2 2014/04/26 17:09:49 rgm Exp rgm $ 

*/

  /* Frameset data from MergeLog */ 
  ml.frameSetID,
  ml.cuEventID,
  ml.ra, ml.dec,

  /* OB ids */
  mfEsoKeysY.ObsID AS YObsID,   
  mfEsoKeysJ.ObsID AS JObsID,   
  mfEsoKeysH.ObsID AS HObsID,   
  mfEsoKeysKs.ObsID AS KsObsID,

  /* OB execution times (used to identify VHS subcomponent) */
  mfEsoKeysY.obsExectime AS YobsExectime,   
  mfEsoKeysJ.obsExectime AS JobsExectime,   
  mfEsoKeysH.obsExectime AS HobsExectime,   
  mfEsoKeysKs.obsExectime AS KsobsExectime,

  /* OB ids */
  mfEsoKeysY.ObsProgID AS YObsProgID,   
  mfEsoKeysJ.ObsProgID AS JObsProgID,   
  mfEsoKeysH.ObsProgID AS HObsProgID,   
  mfEsoKeysKs.ObsProgID AS KsObsProgID,

  /* ESO QC grades */
  mfy.esograde as yesograde,
  mfj.esograde as jesograde,
  mfh.esograde as hesograde,
  mfks.esograde as ksesograde,

  /* ESO Project Progid codes */
  mfy.project as yproject,
  mfj.project as jproject,
  mfh.project as hproject,
  mfks.project as ksproject,

  /* ESO Object keyword */
  mfy.object as yobject,
  mfj.object as jobject,
  mfh.object as hobject,
  mfks.object as ksobject,


  /*  casu version info */
  mfy.casuvers as ycasuversion,
  mfj.casuvers as jcasuversion,
  mfh.casuvers as hcasuversion,
  mfks.casuvers as kscasuversion,

  /* Stats per fs */
  nsources, 
  Ynsources, Jnsources, Hnsources, Ksnsources,
  JKsnsources, JHKsnsources, YJHKsnsources,
  onlyYnsources,   onlyJnsources,   onlyHnsources,   onlyKsnsources,  

  CONVERT(INT, ISNULL(onlyYnsources, 0)) 
    + CONVERT(INT, ISNULL(onlyJnsources, 0))
    + CONVERT(INT, ISNULL(onlyHnsources, 0)) 
    + CONVERT(INT, ISNULL(onlyKsnsources, 0)) 
    AS AllOrphans,

  /* avoid the division by zero by replacing 0 with NULL with NULLIF*/
  convert(float, onlyYnsources) / convert(float, nullif(Ynsources,0)) 
   AS YOrphansFraction,
  convert(float, onlyJnsources) / convert(float, nullif(Jnsources,0)) 
   AS JOrphansFraction,
  convert(float, onlyHnsources) / convert(float, nullif(Hnsources,0)) 
   AS HOrphansFraction,
  CONVERT(float, onlyKsnsources) / CONVERT(float, nullif(Ksnsources,0)) 
   AS KsOrphansFraction,
 
  /* avoid the division by zero by replacing 0 with NULL with NULLIF*/
  (CONVERT(float, ISNULL(onlyYnsources, 0)) 
    + CONVERT(float, ISNULL(onlyJnsources, 0))
    + CONVERT(float, ISNULL(onlyHnsources, 0)) 
    + CONVERT(float, ISNULL(onlyKsnsources, 0))) 
    / CONVERT(float, nullif(nsources, 0))
   AS AllOrphansFraction,

/* Y data */
/* Y Mulitiframe data*/
  mfy.multiframeID as yMultiframeID,
  mfy.Filtername as yFiltername,
  mfy.rabase*15 as yRaBase, 
  mfy.decbase as yDecBase,
  mfy.mjdobs as yMjdObs,
  mfy.exptime as yExpTime,
  rtrim(substring
   (mfy.filename,charindex("v2",mfy.filename,-1),32)) 
   as yfilename,
  rtrim(substring
   (mfy.catname,charindex("v2",mfy.catname,-1),32)) 
   as ycatname,
/* Y data from MultiFrameDetector */
  mfdy.extnum as yExtnum,
  mfdy.tablerows as ytablerows,
  mfdy.axis1length as ynaxis1,
  mfdy.axis2length as ynaxis2,
  mfdy.detrows as ydetrows,
  mfdy.detcols as ydetcols,
  /* mfdy.pixelScale as yPixelScale, */
  mfdy.skylevel as yskylevel,
  mfdy.skynoise as yskynoise,
  mfdy.seeing as yseeingPix,
  mfdy.seeing*cay.xPIxSize as yseeingArcSecs,
  mfdy.avStellarEll as yAvStellarEll,
  mfdy.apercor3 as yAperCor3,
  mfdy.PhotZPCat as yZP,
  mfdy.PhotZPCat as yphotzpcat,
/* Y CurrentAstrometry data */
  cay.centralRa as yCentralRA,
  cay.centralDec as yCentralDec,
  cay.stdCRms as ystdCRms,
  cay.numrms  as ynumrms,
  cay.xPixSize as yXPixSize,
  cay.yPixSize as yYPixSize,
  cay.minRA as yMinRA,
  cay.MaxRA as  YMaxRa,
  cay.minDec as YMinDec,
  cay.maxDec as YMaxDec,
  cay.cTypeX as yCtype1,
  cay.cTypeY as yCtype2,
  cay.crValX as yCRVAL1,
  cay.crValY as yCRVAL2,
  cay.crUnit1 as yCRUNIT1,
  cay.crUnit2 as yCRUNIT2,
  cay.crPixX as yCRPIX1,
  cay.crPixY as yCRPIX2,
  cay.cd11 as ycd1_1,
  cay.cd21 as ycd2_1,
  cay.cd12 as ycd1_2,
  cay.cd22 as ycd2_2,
  cay.pv21 as yPV2_1,
  cay.pv22 as yPV2_2,
  cay.pv23 as yPV2_3,

  /* 5sigma depth */
  /* commented out due to floating point problem */
  /* use abs to avoid -ve info; maybe should use NULLIF or ISNULL */
  /* factor of 1.2 is WFCAM effect of covariance or correlated noise between */ 
  /* pixels due to interpixel capacitance */
   mfdy.PhotZPCat - 2.5 * LOG10(abs(3.0 * mfdy.skyNoise * SQRT(1.2 * 3.141593) / (cay.xPixSize * mfy.expTime)))-mfdy.AperCor3 AS y_depth_dye2006,

  /* j data */
  /* J Mulitiframe data*/
  mfj.multiframeID as jMultiframeID,
  mfj.Filtername as jFiltername,
  mfj.rabase*15 as jRaBase, 
  mfj.decbase as jDecBase,
  mfj.mjdobs as jMjdObs,
  mfj.exptime as jExpTime,
  rtrim(substring
   (mfj.filename,charindex("v2",mfj.filename,-1),32)) 
   as jfilename,
  rtrim(substring
   (mfj.catname,charindex("v2",mfj.catname,-1),32)) 
   as jcatname,
/* j data from MultiFrameDetector */
  mfdj.extnum as jExtnum,
  mfdj.tablerows as jtablerows,
  mfdj.axis1length as jnaxis1,
  mfdj.axis2length as jnaxis2,
  mfdj.detrows as jdetrows,
  mfdj.detcols as jdetcols,
  /* mfdj.pixelScale as jPixelScale, */
  mfdj.skylevel as jskylevel,
  mfdj.skynoise as jskynoise,
  mfdj.seeing as jseeingPix,
  mfdj.seeing*caj.xPIxSize as jseeingArcSecs,
  mfdj.avStellarEll as jAvStellarEll,
  mfdj.apercor3 as jAperCor3,
  mfdj.PhotZPCat as jZP,
  mfdj.PhotZPCat as jphotzpcat,
/* j CurrentAstrometry data */
  caj.centralRa as jCentralRA,
  caj.centralDec as jCentralDec,
  caj.stdCRms as jstdCRms,
  caj.numrms  as jnumrms,
  caj.xPixSize as jXPixSize,
  caj.yPixSize as jYPixSize,
  caj.minRA as jMinRA,
  caj.MaxRA as  jMaxRa,
  caj.minDec as jMinDec,
  caj.maxDec as jMaxDec,
  caj.cTypeX as jCtype1,
  caj.cTypeY as jCtype2,
  caj.crValX as jCRVAL1,
  caj.crValY as jCRVAL2,
  caj.crUnit1 as jCRUNIT1,
  caj.crUnit2 as jCRUNIT2,
  caj.crPixX as jCRPIX1,
  caj.crPixY as jCRPIX2,
  caj.cd11 as jcd1_1,
  caj.cd21 as jcd2_1,
  caj.cd12 as jcd1_2,
  caj.cd22 as jcd2_2,
  caj.pv21 as jPV2_1,
  caj.pv22 as jPV2_2,
  caj.pv23 as jPV2_3,
  /* 5sigma depth */
  /* due to log of -of negative numbers using abs */
  mfdj.PhotZPCat - 2.5 * LOG10(3.0 * abs(mfdj.skyNoise) * SQRT(1.2 * 3.141593) / (abs(caj.xPixSize * mfj.expTime)))-mfdj.AperCor3 AS j_depth_dye2006, 
/* h data */
/* H Mulitiframe data*/
  mfh.multiframeID as hMultiframeID,
  mfh.Filtername as hFiltername,
  mfh.rabase*15 as hRaBase, 
  mfh.decbase as hDecBase,
  mfh.mjdobs as hMjdObs,
  mfh.exptime as hExpTime,
  rtrim(substring
   (mfh.filename,charindex("v2",mfh.filename,-1),32)) 
   as hfilename,
  rtrim(substring
   (mfh.catname,charindex("v2",mfh.catname,-1),32)) 
   as hcatname,
/* h data from MultiFrameDetector */
  mfdh.extnum as hExtnum,
  mfdh.tablerows as htablerows,
  mfdh.axis1length as hnaxis1,
  mfdh.axis2length as hnaxis2,
  mfdh.detrows as hdetrows,
  mfdh.detcols as hdetcols,
  /* mfdh.pixelScale as hPixelScale, */
  mfdh.skylevel as hskylevel,
  mfdh.skynoise as hskynoise,
  mfdh.seeing as hseeingPix,
  mfdh.seeing*cah.xPIxSize as hseeingArcSecs,
  mfdh.avStellarEll as hAvStellarEll,
  mfdh.apercor3 as hAperCor3,
  mfdh.PhotZPCat as hZP,
  mfdh.PhotZPCat as hphotzpcat,
/* h CurrentAstrometry data */
  cah.centralRa as hCentralRA,
  cah.centralDec as hCentralDec,
  cah.stdCRms as hstdCRms,
  cah.numrms  as hnumrms,
  cah.xPixSize as hXPixSize,
  cah.yPixSize as hYPixSize,
  cah.minRA as hMinRA,
  cah.MaxRA as  hMaxRa,
  cah.minDec as hMinDec,
  cah.maxDec as hMaxDec,
  cah.cTypeX as hCtype1,
  cah.cTypeY as hCtype2,
  cah.crValX as hCRVAL1,
  cah.crValY as hCRVAL2,
  cah.crUnit1 as hCRUNIT1,
  cah.crUnit2 as hCRUNIT2,
  cah.crPixX as hCRPIX1,
  cah.crPixY as hCRPIX2,
  cah.cd11 as hcd1_1,
  cah.cd21 as hcd2_1,
  cah.cd12 as hcd1_2,
  cah.cd22 as hcd2_2,
  cah.pv21 as hPV2_1,
  cah.pv22 as hPV2_2,
  cah.pv23 as hPV2_3,

  /* 5sigma depth */
  /* commented out due to floating point problem */
  mfdh.PhotZPCat - 2.5 * LOG10(3.0 * abs(mfdh.skyNoise) * SQRT(1.2 * 3.141593) / (abs(cah.xPixSize) * abs(mfh.expTime)))-mfdh.AperCor3 AS h_depth_dye2006, 

  /* ks data */

  /* Ks Mulitiframe data*/
  mfks.multiframeID as ksMultiframeID,
  mfks.Filtername as ksFiltername,
  mfks.rabase*15 as ksRaBase, 
  mfks.decbase as ksDecBase,
  mfks.mjdobs as ksMjdObs,
  mfks.exptime as ksExpTime,
  rtrim(substring
   (mfks.filename,charindex("v20",mfks.filename,-1),32)) 
   as ksfilename,
  rtrim(substring
   (mfks.catname,charindex("v20",mfks.catname,-1),32)) 
   as kscatname,

/* k data from MultiFrameDetector */
  mfdks.extnum as ksExtnum,
  mfdks.tablerows as kstablerows,
  mfdks.axis1length as ksnaxis1,
  mfdks.axis2length as ksnaxis2,
  mfdks.detrows as ksdetrows,
  mfdks.detcols as ksdetcols,
  /* mfdks.pixelScale as ksPixelScale, */
  mfdks.skylevel as ksskylevel,
  mfdks.skynoise as ksskynoise,
  mfdks.seeing as ksseeingPix,
  mfdks.seeing*caks.xPIxSize as ksseeingArcSecs,
  mfdks.avStellarEll as ksAvStellarEll,
  mfdks.apercor3 as ksAperCor3,
  mfdks.PhotZPCat as ksZP,
  mfdks.PhotZPCat as ksphotzpcat,

/* k CurrentAstrometry data */
  caks.centralRa as ksCentralRA,
  caks.centralDec as ksCentralDec,
  caks.stdCRms as ksstdCRms,
  caks.numrms  as ksnumrms,
  caks.xPixSize as ksXPixSize,
  caks.yPixSize as ksYPixSize,
  caks.minRA as ksMinRA,
  caks.MaxRA as  ksMaxRa,
  caks.minDec as ksMinDec,
  caks.maxDec as ksMaxDec,
  caks.cTypeX as ksCtype1,
  caks.cTypeY as ksCtype2,
  caks.crValX as ksCRVAL1,
  caks.crValY as ksCRVAL2,
  caks.crUnit1 as ksCRUNIT1,
  caks.crUnit2 as ksCRUNIT2,
  caks.crPixX as ksCRPIX1,
  caks.crPixY as ksCRPIX2,
  caks.cd11 as kscd1_1,
  caks.cd21 as kscd2_1,
  caks.cd12 as kscd1_2,
  caks.cd22 as kscd2_2,
  caks.pv21 as ksPV2_1,
  caks.pv22 as ksPV2_2,
  caks.pv23 as ksPV2_3,

  /* 5sigma depth */
  /* commented out due to floating point problem */
  mfdks.PhotZPCat - 2.5 * LOG10(3.0 * abs(mfdks.skyNoise) * SQRT(1.2 * 3.141593) / (abs(caks.xPixSize) * abs(mfks.expTime))) - mfdks.AperCor3 AS ks_depth_dye2006

FROM

  VHSmergelog as ml

  /* Join with Multiframe */
  LEFT OUTER JOIN MultiFrame as mfY
    ON ml.Ymfid = mfY.multiframeID

  LEFT OUTER JOIN MultiFrame as mfJ
    ON ml.Jmfid = mfJ.multiframeID 

  LEFT OUTER JOIN MultiFrame as mfh
    ON ml.hmfid = mfh.multiframeID

  LEFT OUTER JOIN MultiFrame as mfKs
    ON ml.Ksmfid = mfKs.multiframeID


  /* Join with MultiFrameDetector */
  LEFT OUTER JOIN MultiFrameDetector as mfdy
    ON ml.ymfid = mfdy.multiframeID 
    AND ml.yenum = mfdy.extnum 

  LEFT OUTER JOIN MultiFrameDetector as mfdj
    ON ml.jmfid = mfdj.multiframeID 
    AND ml.jenum = mfdj.extnum 

  LEFT OUTER JOIN MultiFrameDetector as mfdh
    ON ml.hmfid = mfdh.multiframeID
    AND ml.henum = mfdh.extnum 

  LEFT OUTER JOIN MultiFrameDetector as mfdks
    ON ml.ksmfid = mfdks.multiframeID
    AND ml.ksenum = mfdks.extnum 


  /* Join with CurrentAstrometry */
  LEFT OUTER JOIN CurrentAstrometry  as cay
    ON  ml.ymfID = cay.multiframeID 
    AND   ml.yenum = cay.extnum 

  LEFT OUTER JOIN CurrentAstrometry  as caj
    ON  ml.jmfID = caj.multiframeID 
    AND  ml.jenum = caj.extnum 

  LEFT OUTER JOIN CurrentAstrometry  as cah
    ON  ml.hmfID = cah.multiframeID 
    AND  ml.henum = cah.extnum 

  LEFT OUTER JOIN CurrentAstrometry  as caks
    ON  ml.ksmfID = caks.multiframeID 
    AND  ml.ksenum = caks.extnum 

  /*   MultiframeEsoKeys */
  LEFT OUTER JOIN MultiframeEsoKeys as mfesokeysY
    ON ml.ymfid  = mfesokeysy.multiframeID 

  LEFT OUTER JOIN MultiframeEsoKeys as mfesokeysJ
    ON ml.jmfid  = mfesokeysj.multiframeID 

  LEFT OUTER JOIN MultiframeEsoKeys as mfesokeysH
    ON ml.hmfid  = mfesokeysh.multiframeID 

  LEFT OUTER JOIN MultiframeEsoKeys as mfesokeysKs
    ON ml.ksmfid  = mfesokeysks.multiframeID 


  /*   MultiframeDetectorEsoKeys */
  LEFT OUTER JOIN MultiFrameDetectorEsoKeys as MFDEKy
    ON  ml.ymfID = mfdeky.multiframeID 
    AND ml.yenum = mfdeky.extnum 

  LEFT OUTER JOIN MultiFrameDetectorEsoKeys as MFDEKj
    ON  ml.jmfID = mfdekj.multiframeID 
    AND ml.jenum = mfdekj.extnum 

  LEFT OUTER JOIN MultiFrameDetectorEsoKeys as MFDEKh
    ON  ml.hmfID = mfdekh.multiframeID 
    AND ml.henum = mfdekh.extnum 

  LEFT OUTER JOIN MultiFrameDetectorEsoKeys as MFDEKks
    ON  ml.ksmfID = MFDEKks.multiframeID 
    AND ml.ksenum = MFDEKks.extnum 

  /* Start of the Counts */
  LEFT OUTER JOIN 
   (SELECT 
      framesetid, count(*) as nSources,

      SUM(CASE WHEN yAperMag3 > -1 THEN 1 ELSE 0 END) AS YnSources,
      SUM(CASE WHEN jAperMag3 > -1 THEN 1 ELSE 0 END) AS JnSources,
      SUM(CASE WHEN hAperMag3 > -1 THEN 1 ELSE 0 END) AS HnSources,
      SUM(CASE WHEN ksAperMag3 > -1 THEN 1 ELSE 0 END) AS KsnSources,

      SUM(CASE WHEN jAperMag3 > -1  AND
        ksAperMag3 > -1 THEN 1 ELSE 0 END) AS JKsnSources,
      SUM(CASE WHEN jAperMag3 > -1 and hapermag3 > -1 
        and ksAperMag3 > -1 THEN 1 ELSE 0 END) AS JHKsnSources,
      SUM(CASE WHEN yAperMag3 > -1 AND jAperMag3 > -1 
        and hapermag3 > -1 and ksAperMag3 > -1 
        THEN 1 ELSE 0 END) AS YJHKsnSources,

      SUM(CASE WHEN yAperMag3 > -1 AND jAperMag3 < -1 
        AND hAperMag3 < -1 AND ksAperMag3 < -1 
        THEN 1 ELSE 0 END) AS onlyYnSources,
      SUM(CASE WHEN yAperMag3 < -1 AND jAperMag3 > -1 
        AND hAperMag3 < -1 AND ksAperMag3 < -1 
        THEN 1 ELSE 0 END) AS onlyJnSources,
      SUM(CASE WHEN yAperMag3 < -1 AND jAperMag3 < -1 
        AND hAperMag3 > -1 AND ksAperMag3 < -1 
        THEN 1 ELSE 0 END) AS onlyHnSources,
      SUM(CASE WHEN yAperMag3 < -1 AND jAperMag3 < -1 
        AND hAperMag3 < -1 AND ksAperMag3 > -1 
        THEN 1 ELSE 0 END) AS onlyKsnSources

    FROM vhssource
    
    GROUP BY framesetid
    ) AS STATS ON ml.framesetid = STATS.framesetid

ORDER BY

  ml.framesetid ASC
