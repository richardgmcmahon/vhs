"""



"""
from __future__ import print_function, division

import glob
import os
import sys
import time
from time import strftime, gmtime

from os.path import basename

import tempfile, logging
import traceback

import matplotlib as mpl
print('matplotlib: ', mpl.__version__)
print('matplotlib.matplotlib_fname(): ',mpl.matplotlib_fname())
# Force matplotlib to not use any Xwindows backend:
try:
  mpl.use('Agg')
except:
  pass

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

import numpy as np

import astropy
print('astropy.__version__: ', astropy.__version__)
from astropy.table import Table, join

from astropy.coordinates import Angle
import astropy.units as u


# use sys.path.append of functions not in PYTHONPATH
sys.path.append('/home/rgm/soft/python/lib/')
from librgm.table_stats import table_stats
from librgm.plotid import plotid
from librgm.plot_radec import plot_radec
#from librgm import des

from rd_sadt import rd_sadt
import plot_radec_tile as plt_tile

import des_footprint

# functions for ESO progress files
import ob_progress as ob_progress

debug=False



def rd_obprogress(url=False, date=None):
  """
  Read in VISTA ESO progress file
  Only reads VHS for now

  specify url so can read over internet

  """

  url_default='http://www.ast.cam.ac.uk/~rgm/vhs/progress/data/current/'

  inpath='/data/vhs/progress/'
  progid='179A2010'
  # no changes should be needed below
  if date is None: datestamp = strftime("%Y%m%d", gmtime())
  if date is not None: datestamp = date

  datestamp = strftime("%Y%m%d", gmtime())

  infile = inpath + date + '/' + progid + '.fits'

  print('Current working directory: ', os.getcwd())
  print('Reading: ', infile)
  print('Path: ', os.path.dirname(infile))
  print('Absolute Path: ', os.path.abspath(infile))

  # read table
  ob_table=Table.read(infile)
  # add filename into the table metadata
  ob_table.meta['filename'] = infile

  # how to access the filename later
  # 1. filename=ob_table.meta['filename']
  # 2. print('plotting: ',ob_table.meta['filename'])
  # 3. print('plotting: ' os.path.basename(ob_table.meta['filename']))
  # where 3. shows how to trim of the path of the filename

  print(ob_table.colnames)
  print('table.meta: ', ob_table.meta)
  table_stats(infile)

  return ob_table


def plot_sadt_tile(sadtfile=None, sadtfiles=None, sadtfile_path=None,
                   xlimits=None, ylimits=None, raUnits=None,
                   color='orange'):
  """

  """

  from astropy.table import vstack

  print('plot_sadt_file:')
  if xlimits is None: xlimits=(-180.0,180.0)
  if ylimits is None: ylimits=(-90.0,2.0)
  print('xlimits: ', xlimits)
  print('ylimits: ', ylimits)


  print('sadtfile:', sadtfile)
  print('sadtfiles:', sadtfiles)
  print('sadtfile_path:', sadtfile_path)

  raw_input("Enter any key to continue: ")

  # read a single sadtfile
  if sadtfile_path is None and sadtfiles is None :
    sadt = rd_sadt(infile=sadtfile, csv=False)
    print('Path: ', os.path.dirname(sadtfile))
    sadt.meta['inpath']=os.path.dirname(sadtfile)
    sadt.meta['filename']=sadtfile

  if sadtfile_path is not None or sadtfiles is not None:

    if sadtfiles is None:
        filelist = os.listdir(sadtfile_path)
    if sadtfiles is not None:
        filelist = sadtfiles

    print('sadtfile_path: ', sadtfile_path)
    print('Number of files in filelist:', len(filelist))
    print('filelist:', filelist)
    ifile=0
    for file in filelist:
      ifile = ifile + 1
      print('reading file: ', file)
      #if sadtfile_path is not None:
      #    sadtfile = sadtfile_path + '/' + file
      #if sadtfile_path is None:
      #   sadtfile = file
      sadtfile = file


      if ifile == 1: sadt=rd_sadt(infile=sadtfile, csv=False)
      print('Number of SADT tiles: ', len(sadt))
      if ifile > 1:
          temp=rd_sadt(infile=sadtfile, csv=False)
          print('Number of SADT tiles: ', len(temp))
          sadt=vstack([sadt, temp])
          print('Number of SADT tiles: ', len(sadt))
          del temp

      print('Path: ', os.path.dirname(sadtfile))
      sadt.meta['inpath']=os.path.dirname(sadtfile)
      sadt.meta['filename'+str(ifile)]=sadtfile

  print('Number of SADT tiles: ', len(sadt))
  print(sadt.colnames)

  raw_input("Enter any key to continue: ")

  xdata=sadt['ra']
  ydata=sadt['dec']
  print('xrange: ', np.min(xdata), np.max(xdata))
  print('yrange: ', np.min(ydata), np.max(ydata))
  print('ndata: ', len(xdata))

  ntiles=len(xdata)
  plt.plot(xdata, ydata, '+k', ms=1.0, label='SADT Tiles: ' + str(ntiles))
  plt.title(desfile)
  plt.suptitle(sadtfile)
  plotid(progname=True)
  plt.grid()
  plt.legend(markerscale=2.0, fontsize='medium')

  # filled polygon
  ra=xdata
  dec=ydata
  for i in range(len(ra)):
      if i==0 or i == len(ra): print('i: ', i)

      ra_poly, dec_poly=plt_tile.mk_polytile(ra_cen=ra[i], dec_cen=dec[i],
       coverage='twice')

      xypolygon=np.column_stack((ra_poly, dec_poly))
      if i==0 or i == len(ra):
        print('xypolygon.shape: ',   xypolygon.shape)
        print('xypolygon: ',xypolygon)

      polygon = Polygon(xypolygon, True, color=color, alpha=0.2)
      plt.gca().add_patch(polygon)

      ra_poly, dec_poly=plt_tile.mk_polytile(ra_cen=ra[i], dec_cen=dec[i],
       coverage='single')

      xypolygon=np.column_stack((ra_poly, dec_poly))
      if i==0 or i == len(ra):
        print('xypolygon.shape: ',   xypolygon.shape)
        print('xypolygon: ',xypolygon)

      polygon = Polygon(xypolygon, True, color=color, alpha=0.1)
      plt.gca().add_patch(polygon)

  plt.xlim(xlimits)
  plt.ylim(ylimits)


  print('xlimits: ', xlimits)
  print('ylimits: ', ylimits)


  plotdir='./'
  figname='vhs_des_check'+ '_sadt_'+ datestamp + '.png'
  print('Saving: '+figname)
  plt.savefig(plotdir+figname)


def plot_viking_tiles(ra=None, dec=None, PA=0.0,
    wrap_24hr=True, overplot=True):
    """
    ra, dec in degrees

    """

    plot_polygon=True

    xdata=ra
    ydata=dec
    angle=Angle(xdata * u.deg)
    angle.wrap_at('180d', inplace=True)
    xdata=angle.degree

    print('plotting: ',vhs_obprogress.meta['filename'])
    filename=vhs_obprogress.meta['filename']

    if not plot_polygon:
      plt.plot(xdata, ydata, 's', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted)\n' + filename)

    if plot_polygon:
      plt.plot(xdata, ydata, '.', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted):' + filename)

      # filled polygon
      ra=xdata
      dec=ydata
      for i in range(len(ra)):
        if i==0 or i == len(ra): print('i: ', i)
        ra_poly, dec_poly=plt_tile.mk_polytile(ra_cen=ra[i], dec_cen=dec[i],
         coverage='twice', PA=PA)

        xypolygon=np.column_stack((ra_poly, dec_poly))
        if i==0 or i == len(ra):
          print('xypolygon.shape: ',   xypolygon.shape)
          print('xypolygon: ',xypolygon)

        polygon = Polygon(xypolygon, True, color='green', alpha=0.1)
        plt.gca().add_patch(polygon)

        #print(ob_table.meta)

      plt.suptitle(dqcfile_vhs)
      plt.legend(fontsize='small')

      figname='vhs_des_check_progress_vhs_obprogress_' + datestamp + '.png'
      print('Saving: '+figname)
      plt.savefig(plotdir+figname)
      plt.suptitle('')


def plot_vhs_obprogress(vhs_obprogress=None, wrap_ra24hr=True,
    overplot=True, savefig=True):
    """

    """

    plot_polygon=True

    PA=90.0

    xdata=vhs_obprogress['RA (hrs)']*15.0
    ydata=vhs_obprogress['DEC (deg)']
    angle=Angle(xdata * u.deg)
    angle.wrap_at('180d', inplace=True)
    xdata=angle.degree

    print('plotting: ',vhs_obprogress.meta['filename'])
    filename=vhs_obprogress.meta['filename']

    if not plot_polygon:
      plt.plot(xdata, ydata, 's', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted)\n' + filename)

    if plot_polygon:
      plt.plot(xdata, ydata, '.', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted):' + filename)

      # filled polygon
      ra=xdata
      dec=ydata
      for i in range(len(ra)):
        if i==0 or i == len(ra): print('i: ', i)
        ra_poly, dec_poly=plt_tile.mk_polytile(ra_cen=ra[i], dec_cen=dec[i],
         coverage='twice', PA=PA)

        xypolygon=np.column_stack((ra_poly, dec_poly))
        if i==0 or i == len(ra):
          print('xypolygon.shape: ',   xypolygon.shape)
          print('xypolygon: ',xypolygon)

        polygon = Polygon(xypolygon, True, color='green', alpha=0.1)
        plt.gca().add_patch(polygon)

        #print(ob_table.meta)

      plt.suptitle(dqcfile_vhs)
      plt.legend(fontsize='small')

      figname='vhs_des_check_progress_vhs_obprogress_' + datestamp + '.png'
      print('Saving: '+figname)
      plt.savefig(plotdir+figname)
      plt.suptitle('')


def plot_vhs_tiles(table=None,ra=None, dec=None,
    wrap_ra24hr=False, PA=90.0, overplot=True, savefig=True):
    """

    table is astropy table
    ra, dec in degrees

    """

    plot_polygon=True

    if wrap_ra24hr:
        angle=Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata=angle.degree

    # assumes filename is in table.meta
    print('plotting: ', table.meta['filename'])
    filename = table.meta['filename']

    if not plot_polygon:
      plt.plot(xdata, ydata, 's', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted)\n' + filename)

    if plot_polygon:
      plt.plot(xdata, ydata, '.', ms=7.0,
       color='yellow', markeredgecolor='yellow',
       alpha=0.1,
       label='OB Progress (submitted):' + filename)

      # filled polygon
      ra=xdata
      dec=ydata
      for i in range(len(ra)):
        if i==0 or i == len(ra): print('i: ', i)
        ra_poly, dec_poly=plt_tile.mk_polytile(ra_cen=ra[i], dec_cen=dec[i],
         coverage='twice', PA=PA)

        xypolygon=np.column_stack((ra_poly, dec_poly))
        if i==0 or i == len(ra):
          print('xypolygon.shape: ',   xypolygon.shape)
          print('xypolygon: ',xypolygon)

        polygon = Polygon(xypolygon, True, color='green', alpha=0.1)
        plt.gca().add_patch(polygon)

        #print(ob_table.meta)

      plt.suptitle(dqcfile_vhs)
      plt.legend(fontsize='small')

      figname='vhs_des_check_progress_vhs_obprogress_'+ datestamp + '.png'
      print('Saving: '+figname)
      plt.savefig(plotdir+figname)
      plt.suptitle('')


def plot_phase2_des_check():
    """


    """

    datestamp = strftime("%Y%m%d", gmtime())

    plotdir='./'
    figname='vhs_des_check'+ '_sadt_check_' + datestamp + '.png'
    print('Saving: '+figname)
    plt.savefig(plotdir+figname)


if __name__ == '__main__':

    import argparse

    from argparse import ArgumentParser

    global t0

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='check des footprint against SADT xml Tile layout')

    date_default = time.strftime("%Y%m%d", gmtime())
    parser.set_defaults(date=date_default)
    parser.add_argument(
        "--date", default=date_default,
        dest='date',
        help="date as string e.g. '20151101' [DEFAULT is today(gmt)]")

    parser.add_argument(
        "--rarange", default=[0.0, 360.0], type=float, nargs=2,
        help="RA range in degrees in form Deg Deg]")

    parser.add_argument(
        "--decrange", default=[-90.0, 90.0], type=float, nargs=2,
        help="Declination range in degrees in form Deg Deg")


    parser.add_argument(
        "--raunits", default="hour",
        help="RA units hour or degree; default = 'hour'")


    parser.set_defaults(path=None)
    parser.add_argument(
        "--path",
        help="Path to SADT xml files")


    print('Number of arguments:', len(sys.argv), 'arguments: ', sys.argv[0])
    args = parser.parse_args()
    date = args.date

    rarange = args.rarange
    print('RA range:', rarange)

    decrange = args.decrange
    print('Dec range:', decrange)

    raUnits = args.raunits
    print('raUnits:', raUnits)

    t0 = time.time()

    datestamp = strftime("%Y%m%d", gmtime())

    if debug: help(plt_tile)

    plotdir='./'


    t0=time.time()


    # extend North to give legend some space
    # show DES regions
    xlimits=(-60.0,120.0)
    xlimits=(120.0, -60.0)
    xlimits=(180.0, -180.0)
    ylimits=(-90.0,90.0)

    xlimits = rarange
    ylimits = decrange

    # DES related files
    desfile = '/home/rgm/soft/des/round13-poly.txt'
    descat = '/data/des/Y2Q1/mag_psf_i_leq_18p5_thin.fits'

    # VHS SADT files
    # maybe change to a directory
    sadtfile_path = '/home/rgm/Projects/VHS/Phase2/Period96/20150819'
    # 20160803
    sadtfile_path = '/data/mbanerji/Projects/VHS_Survey/SADT_outputs/Boundary_tests/'

    sadtfile_path = '/data/mbanerji/Projects/VHS_Survey/SADT_outputs/Boundary_tests/Dec2017/'

    sadtfiles = ['Str45_DES_RA_XXX_XXX.xml',
                 'Str46_DES_RA_XXX_XXX.xml',
                 'Str47_DES_RA_XXX_XXX.xml',
                 'Str48_DES_RA_XXX_XXX.xml']

    # Nov 2016 SADT xml files
    sadtfiles = ['Str08_DES_RA_XXX_XXX.xml',
                 'Str09_DES_RA_XXX_XXX.xml',
                 'Str10_DES_RA_XXX_XXX.xml',
                 'Str11_DES_RA_XXX_XXX.xml',
                 'Str12_DES_RA_XXX_XXX.xml']

    # sadtfile='Str02_DES_RA_XXX_XXX.xml'
    # sadtfile = sadtfiles[0]

    sadtfiles = sorted(glob.glob(sadtfile_path + '*.xml'))
    # filelist=os.listdir(sadtfile_path)

    # _GPS
    # sadtfiles = sorted(glob.glob(sadtfile_path + '*_GPS*.xml'))

    if args.path is not None: sadtfile_path = args.path
    print('sadtfile_path: ', sadtfile_path)
    sadtfiles = sorted(glob.glob(sadtfile_path + '*.xml'))
    for sadtfile in sadtfiles:
        print('SADT file: ', sadtfile)

    raw_input("Enter any key to continue: ")

    # VHS and other VISTA footprint files
    dqcfile_vhs='/data/vhs/dqc/casu/2015/vistaqc_20150531_tiles_vhs.fits'
    dqcfile_viking ='/data/vhs/dqc/casu/2015/vistaqc_20150531_tiles_viking.fits'

    vhs_dqc=False

    vhs_progress=False

    vhs_sadt = True

    viking_dqc = True

    des_y2q1 = False

    vhs_progress_dr3 = False
    if vhs_progress_dr3:
        eso_date='20130930'
        print('Reading vhs_progress_dr3: ' + eso_date)
        vhs_obprogress_dr3=rd_obprogress(date=eso_date)

        table=vhs_obprogress_dr3
        ExecutionTime = ob_progress.rd_ExecutionTimes(table=table, debug=False)

        DES = ob_progress.vhs_select_survey(ExecutionTime, survey='DES')
        vhs_obprogress_dr3 = vhs_obprogress_dr3[DES]

        table=vhs_obprogress_dr3
        OBstatus = ob_progress.rd_OBstatus(table=table, debug=False)

        completed=True
        if completed:
            itest = (OBstatus == 'C')
            print(len(itest), len(OBstatus[itest]))

            vhs_obprogress_dr3 = vhs_obprogress_dr3[itest]
            print('Number of completed VHS DES OBs: ', len(vhs_obprogress_dr3))

        raw_input("Enter any key to continue: ")

    vhs_progress_20150930 = False
    if vhs_progress_20150930:
        eso_date='20150930'
        print('Reading vhs_progress_dr3: ' + eso_date)
        vhs_obprogress_20150930=rd_obprogress(date=eso_date)

        table=vhs_obprogress_20150930
        ExecutionTime = ob_progress.rd_ExecutionTimes(table=table, debug=False)

        DES = ob_progress.vhs_select_survey(ExecutionTime, survey='DES')
        vhs_obprogress_20150930 = vhs_obprogress_20150930[DES]

        table=vhs_obprogress_20150930
        OBstatus = ob_progress.rd_OBstatus(table=table, debug=False)

        completed=True
        if completed:
            itest = (OBstatus == 'C')
            print(len(itest), len(OBstatus[itest]))

            vhs_obprogress_20150930 = vhs_obprogress_20150930[itest]

            print('Number of completed VHS DES OBs: ',
                len(vhs_obprogress_20150930))

        raw_input("Enter any key to continue: ")

    if vhs_progress:
        print('Read VHS progress from:', date)
        vhs_obprogress=rd_obprogress(date=date)
        raw_input("Enter any key to continue: ")

    if viking_dqc:
        print('Read VIKING DQC')
        viking=Table.read(dqcfile_viking)
        viking.meta['filename'] = dqcfile_viking
        print('Table column names: \n', viking.colnames)
        print('Table metadata: \n', viking.meta)
        viking.pprint


    if vhs_dqc:
        print('Read VHS DQC')
        vhsdqc=Table.read(dqcfile_vhs)
        vhsdqc.meta['filename'] = dqcfile_vhs
        print('Table column names: \n', vhsdqc.colnames)
        print('Table metadata: \n', vhsdqc.meta)
        vhsdqc.pprint

    despolygon = des_footprint.rd_despolygon(desfile)

    raw_input("Enter any key to continue: ")
    print('xlimits:', xlimits)
    print('ylimits:', ylimits)

    fig = plt.figure(figsize=(12,6))

    if vhs_sadt:
        print('sadtfile_path:', sadtfile_path)
        plot_sadt_tile(sadtfile=None, sadtfiles=sadtfiles,
                       sadtfile_path=sadtfile_path,
                       xlimits=xlimits, ylimits=ylimits)

    print('xlimits:', xlimits)
    print('ylimits:', ylimits)
    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.title(desfile)
    plt.xlabel('RA (degrees)')
    plt.ylabel('Declination (degrees)')

    xdata=despolygon['ra']
    ydata=despolygon['dec']

    plt.plot(xdata, ydata, label='DES Round13-poly')
    plotid(progname=True)

    figname='vhs_des_check_sadt_despolygon_' + datestamp + '.png'
    print('Saving: '+figname)
    plt.savefig(plotdir+figname)
    plt.suptitle('')


    if viking:
        xdata=viking['ra']
        ydata=viking['dec']
        angle=Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata=angle.degree

        plt.plot(xdata, ydata, 'sr',
            ms=7.0, markeredgecolor='r', alpha=0.2, label='VIKING')
        plt.suptitle(dqcfile_viking)

        figname='vhs_des_check_viking_1_' + datestamp + '.png'
        print('Saving: '+figname)
        plt.savefig(plotdir+figname)
        plt.suptitle('')

    if vhs_progress:
        print('Plot VHS obprogress data')
        plot_vhs_obprogress(vhs_obprogress=vhs_obprogress)
        raw_input("Enter any key to continue: ")

    fig = plt.figure(figsize=(12,6))

    plt.xlim(xlimits)
    plt.ylim(ylimits)
    plt.title(desfile)
    plt.xlabel('RA (degrees)')
    plt.ylabel('Declination (degrees)')

    xdata=despolygon['ra']
    ydata=despolygon['dec']

    plt.plot(xdata, ydata)
    plotid(progname=True)

    if viking:
        xdata=viking['ra']
        ydata=viking['dec']
        angle=Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata=angle.degree

        plt.plot(xdata, ydata, 'sr',ms=7.0, markeredgecolor='r', alpha=0.2)
        plt.suptitle(dqcfile_viking)

        figname='vhs_des_check_viking_1_' + datestamp + '.png'
        print('Saving: '+figname)
        plt.savefig(plotdir+figname)
        plt.suptitle('')




    raw_input("Enter any key to continue: ")

    if des_y2q1:
        print('Reading : ', descat)
        y2q1=Table.read(descat)
        y2q1.meta['filename'] = descat
        print(y2q1.colnames)

        xdata=y2q1['RA']
        print(np.min(xdata), np.max(xdata))
        angle=Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata=angle.degree
        ydata=y2q1['DEC']
        print(np.min(xdata), np.max(xdata))
        print(np.min(ydata), np.max(ydata))

        ndata=len(xdata)
        plt.plot(xdata, ydata, '.', ms=0.1, color='blue', alpha=0.1,
            markeredgecolor=None,
        label='Y2Q1 sources: ' +str(ndata))
        plt.title(desfile)
        plt.suptitle(sadtfile)

        xdata=despolygon['ra']
        ydata=despolygon['dec']

        plt.plot(xdata, ydata, label='DES Round13')
        plotid(progname=True)
        plt.grid()

        plt.legend(fontsize='small')
        plotid(progname=True)
        plotdir='./'
        figname='vhs_des_check'+ '_sadt_y2q1_' + datestamp + '.png'
        print('Saving: '+figname)
        plt.savefig(plotdir+figname)


    if vhs_dqc:
        xdata=vhsdqc['ra']
        ydata=vhsdqc['dec']
        angle=Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata=angle.degree

        plt.plot(xdata, ydata, 'sy',ms=7.0, markeredgecolor='y', alpha=0.2,
         label=os.path.basename(dqcfile_vhs))

        plt.suptitle(dqcfile_vhs)

        plt.legend(fontsize='small')

        figname='vhs_des_check_viking_2_' + datestamp + '.png'
        print('Saving: '+figname)
        plt.savefig(plotdir+figname)
        plt.suptitle('')


    if vhs_progress:
        ob_table=vhs_obprogress
        plot_vhs_obprogress(vhs_obprogress=ob_table)
