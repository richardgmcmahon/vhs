
"""

Read and Plot VISTA survey progress files based on IDL program ob_progress.pro

ISSUES: show we use figfile or plotfile consistently?

plt.savefig(figdir + figname)
plt.savefig(figpath + figname)
plt.savefig(plotpath + plotname)
plt.savefig(plotdir + plotname)

I propose that we use figpath and figname

path is more generic than dir or folder

https://matplotlib.org/faq/usage_faq.html#parts-of-a-figure


Input formats to be supported:

1: ESO Portal csv or FITS format files
2: CASU DQC FITS file
3: VSA based DQC FITS file
4. RA, Dec source lists [need to support various units]

Example(s):

python ob_progress.py --help

Current status including incomplete OBs
python ob_progress.py

Status of run A on 20151228
python ob_progress.py -d '20151228' -r 'A'

Completed OBs
python ob_progress.py --obstatus 'C'

Uses http://matplotlib.org/api/patches_api.html#matplotlib.patches.Polygon
ESO OB Status codes
   Status  + (Accepted)
   Status  - (Rejected)
   Status  A  (Aborted)
   Status  C (Completed)
   Status  D (Defined)
   Status  K  (Cancelled)
   Status  M (Must repeat)
   Status  P (Partially defined)
   Status  S  (Started)
   Status  X

"""

from __future__ import print_function, division
__version__ = "v0.0.1"

import os
import sys
import time
from time import strftime
from time import gmtime

import traceback

import numpy as np
print('numpy.__version__: ', np.__version__)

from matplotlib import pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon

import astropy
print('astropy.__version__: ', astropy.__version__)

from astropy.table import Table

from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import HeliocentricTrueEcliptic
from astropy.coordinates import BarycentricTrueEcliptic
from astropy.coordinates import GeocentricTrueEcliptic
from astropy.coordinates import Angle

import astropy.units as u

# use sys.path.append of functions not in PYTHONPATH
sys.path.append('/home/rgm/soft/python/lib/')
from librgm.table_stats import table_stats
from librgm.plotid import plotid
from librgm.plot_radec import plot_radec

import plot_radec_tile as plt_tile

import des_footprint


import ConfigParser
Config = ConfigParser.RawConfigParser()

# read config file
Config.read('ob_progress.cfg')
print(Config.sections())

dqcfile_viking = Config.get('obprogress', 'dqcfile_casu_viking')
print(Config.get('obprogress', 'dqcfile_casu_viking'))

sourcefile_viking = Config.get('obprogress', 'sourcefile_vsa_viking')
print(Config.get('obprogress', 'sourcefile_vsa_viking'))


# read fits progress file
# progresspath='/data/vhs/progress/'
progresspath = Config.get('obprogress', 'inpath') + '/'


def set_obstatus_codes():
    """

    """
    obstatus_codes = [['+', 'Accepted'],
                      ['-', 'Rejected'],
                      ['A', 'Aborted'],
                      ['C', 'Completed'],
                      ['D', 'Defined'],
                      ['K', 'Cancelled'],
                      ['M', 'Must repeat'],
                      ['P', 'Partially defined'],
                      ['S', 'Started'],
                      ['X', 'Status X']]

    return obstatus_codes


def plot_vista_tiles(table=None, ra=None, dec=None, PA=90.0,
                     raUnits='degrees', wrap_ra24hr=False, aitoff=False,
                     filename=None,
                     title=None, xlabel=None, ylabel=None, label=None,
                     verbose=False, debug=False,
                     color=None, alpha=None,
                     overplot=False,
                     savefig=True, plotdir=None, figfile=None,
                     suffix=None,
                     rarange=[0.0, 24.0], decrange=[-90.0, 10.0],
                     figsize=(12.0, 6.0)):
    """
    see also librgm.plot_radec
    table is an astropy table
    INPUT: ra, dec in degrees

    """

    trace = traceback.extract_stack()[-2]
    print(os.path.basename(trace[0]), ':', str(trace[1]))
    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ':', str(trace[1]))

    trace = traceback.extract_stack()
    print(trace)

    t0 = time.time()
    print('Elapsed time(secs): ', time.time() - t0)

    print('plot_vista_tiles: label;', label)
    plot_polygon = True

    # configure for aitoff
    if aitoff:
        raUnits = 'radians'
        wrap_ra24hr = False

    # matplotlib aitoff only support full sky
    if rarange is None and not aitoff:
        plt.xlim([0, 24.0])
        if raUnits == 'degrees':
            plt.xlim([0, 360.0])
    if rarange is not None and not aitoff:
        plt.xlim(rarange)

    if decrange is None and not aitoff:
        plt.ylim([-90, 30])
    if decrange is not None and not aitoff:
        plt.ylim(decrange)

    print('rarange:  ', rarange)
    print('decrange: ', decrange)

    if alpha is None:
        alpha = 1.0
    if color is None:
        color = 'green'

    if not overplot:
        # plt.clf()
        plt.figure(num=None, figsize=figsize)
        if aitoff:
            plt.subplot(111, projection="aitoff")

    if raUnits == 'degrees':
        plt.xlabel('RA (degrees)')
    if raUnits == 'hours':
        plt.xlabel('RA (hours)')
    if xlabel is not None:
        plt.xlabel(xlabel)
    plt.ylabel('Declination (degrees)')
    if ylabel is not None:
        plt.ylabel(ylabel)
    # if title != None: plt.title(title)

    xdata = ra
    ydata = dec
    ndata = len(ra)

    print('RA range: ', np.min(xdata), np.max(xdata))
    print('Dec range:', np.min(ydata), np.max(ydata))

    # convert to degrees for polygon and wrap calc
    if raUnits == 'hours':
        xdata = xdata * 15.0

    print('RA range: ', np.min(xdata), np.max(xdata))
    print('Dec range:', np.min(ydata), np.max(ydata))
    if wrap_ra24hr:
        angle = Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata = angle.degree

    print('RA range: ', np.min(xdata), np.max(xdata))
    print('Dec range:', np.min(ydata), np.max(ydata))

    print('Elapsed time(secs): ', time.time() - t0)

    # assumes filename is in table.meta
    if table is None and filename is None:
        filename = ''
    if table is not None:
        filename = table.meta['filename']
        print('plotting: ', filename)

    if not plot_polygon:
        if raUnits == 'hours':
            xdata = xdata / 15.0
        plt.plot(xdata, ydata, 's', ms=7.0,
                 color='yellow', markeredgecolor='yellow',
                 alpha=alpha,
                 label='OB Progress (submitted):' + str(ndata))
    print('Elapsed time(secs):', time.time() - t0)

    # phantom  plot to get the label onto legend
    if plot_polygon:
        if label is None:
            label = ''
        else:
            label = label + ': ' + str(ndata)

        plt.plot([-999], [-999], 's', ms=0.1,
                 color=color, markeredgecolor=color,
                 alpha=alpha, label=label)
        plt.grid()

    # do not plot the polygon when color = 'none' as used for adding info
    # to the legend
    if plot_polygon and color != 'none':
        if raUnits == 'hours':
            xdata = xdata / 15.0

        # filled polygon
        ra = xdata
        dec = ydata
        if raUnits == 'hours':
            ra = ra * 15.0
        print('Elapsed time(secs): ', time.time() - t0)
        for i in range(len(ra)):
            if i == 0 or i == len(ra):
                print('i: ', i)
            ra_poly, dec_poly = plt_tile.mk_polytile(
                ra_cen=ra[i], dec_cen=dec[i],
                coverage='twice', PA=PA)

            if raUnits == 'hours':
                ivertex = -1
                if debug:
                    print('ra_poly: ', ra_poly)
                for vertex in ra_poly:
                    ivertex = ivertex + 1
                    ra_poly[ivertex] = ra_poly[ivertex] / 15.0
                if debug:
                    print('ra_poly: ', ra_poly)

            if debug:
                print('Elapsed time(secs): ', time.time() - t0)
            if aitoff:
                # convert to wrapped radians
                ivertex = -1
                if debug:
                    print('ra_poly: ', ra_poly)
                for vertex in ra_poly:
                    ivertex = ivertex + 1
                    c = SkyCoord(ra=ra_poly[ivertex],
                                 dec=dec_poly[ivertex], unit=(u.deg, u.deg),
                                 frame='icrs')
                    ra_poly[ivertex] = c.ra.wrap_at(180 * u.deg).radian
                    dec_poly[ivertex] = c.dec.radian
                if debug:
                    print('ra_poly: ', ra_poly)
                    print('Elapsed time(secs): ', time.time() - t0)

            xypolygon = np.column_stack((ra_poly, dec_poly))
            if i == 0 or i == len(ra):
                print('xypolygon.shape:', xypolygon.shape)
                print('xypolygon:', xypolygon)

            # http://matplotlib.org/api/patches_api.html#matplotlib.patches.Polygon
            polygon = Polygon(xypolygon, closed=True,
                              edgecolor='none', facecolor=color, alpha=alpha)
            plt.gca().add_patch(polygon)
            # this does nothing
            plt.set_label = 'label test for polygon'

    print('Elapsed time(secs):', time.time() - t0)

    plt.title(filename)
    plt.legend(fontsize='small', markerscale=50,
               scatterpoints=1, numpoints=1, ncol=2, loc='upper center',
               prop={'family': 'monospace', 'size': 'small'})
    plt.grid()

    if savefig:
        datestamp = time.strftime("%Y%m%d", gmtime())
        if plotdir is None:
            plotdir = ''
        plotid(progname=True)
        prefix = os.path.splitext(__file__)[0]
        if suffix is None:
            suffix = ''
        if suffix is not None:
            suffix = suffix + '_'
        suffix = os.path.basename(filename) + '_' + suffix
        figname = prefix + '_' + suffix + \
            'radec_vista_tiles_' + datestamp + '.png'
        print('Saving:' + figname)
        plt.savefig(plotdir + figname)


def rd_ExecutionTimes(table=None, debug=False):
    """ reads Execution times from ESO portal file

    table is read in with astropy table

    """

    ExecutionTime = table['Execution time (s)']

    # unique_ExecutionTimes=set(ExecutionTime)
    unique_ExecutionTimes = np.unique(ExecutionTime)

    print('Number of unique Execution times: ', len(unique_ExecutionTimes))
    print(unique_ExecutionTimes)

    if debug:
        help(unique_ExecutionTimes)
        help(ExecutionTime)

    print(ExecutionTime.dtype)
    print(ExecutionTime.shape)

    # isort=np.argsort(unique_ExecutionTimes)
    # print(isort)

    for time in unique_ExecutionTimes:
        if debug:
            help(time)
        # print(time.shape)
        itest = (ExecutionTime == time)
        print(time, ':', len(ExecutionTime[itest]))

    return ExecutionTime


def rd_OBstatus(table=None, debug=True):
    """ reads OB status from ESO portal file

    table is read in with astropy table

    """

    OBstatus = table['OB status']

    # unique_ExecutionTimes=set(ExecutionTime)
    unique = np.unique(OBstatus)

    print('Number of unique OB status:', len(unique))
    print(unique)

    # isort=np.argsort(unique_ExecutionTimes)
    # print(isort)

    for status in unique:
        # print(time.shape)
        itest = (status == OBstatus)
        print(status + ':', len(OBstatus[itest]))

    return OBstatus


def rd_OBstatusDate(table=None, debug=True):
    """ reads OB status from ESO portal file

    table is read in with astropy table

    """

    OBstatusDate = table['Status date']

    print('OB status date range:',
          min(OBstatusDate), max(OBstatusDate))

    return OBstatusDate


def vhs_select_survey(ExecutionTime, survey='DES'):
    """Select VHS survey OBs via ExecutionTime

    """
    if survey == 'ATLAS':
        SELECTED = (
            (ExecutionTime == 1510) |
            (ExecutionTime == 1543) |
            (ExecutionTime == 1549) |
            (ExecutionTime == 1684) |
            (ExecutionTime == 1910) |
            (ExecutionTime == 1649)
        )
        print('Number of VHS ATLAS OBs:', len(ExecutionTime[SELECTED]))

    if survey == 'DES':
        SELECTED = (
            (ExecutionTime == 1809) |
            (ExecutionTime == 1811) |
            (ExecutionTime == 1829) |
            (ExecutionTime == 1845) |
            (ExecutionTime == 2045) |
            (ExecutionTime == 2051) |
            (ExecutionTime == 2129)
        )
        print('Number of VHS DES OBs:   ', len(ExecutionTime[SELECTED]))

    if survey == 'GPS':
        SELECTED = (
            (ExecutionTime == 829) |
            (ExecutionTime == 842) |
            (ExecutionTime == 887) |
            (ExecutionTime == 1005)
        )
        print('Number of VHS GPS OBs:   ', len(ExecutionTime[SELECTED]))

    # these are Dry Run test on the CDFS
    if survey == 'OTHER':
        SELECTED = \
            (ExecutionTime == 1717) | \
            (ExecutionTime == 1961) | \
            (ExecutionTime == 2111)
        print('Number of VHS OTHER OBs: ', len(ExecutionTime[SELECTED]))

    return SELECTED


def vhs_select_status(OBstatus, status='C'):
    """

    """
    STATUS = (OBstatus == status)

    return STATUS


def vhs_select_run(run='A'):
    """Select VHS survey Run via run

    """

    return RUN


def aitoff_test(ra, dec):
    """

    """
    # aitoff example
    ra = ra * 15.0
    c = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    ra_rad = c.ra.wrap_at(180.0 * u.degree).radian
    dec_rad = c.dec.radian

    plt.figure(figsize=(8.0, 6.0))
    plt.clf()
    plt.subplot(111, projection="aitoff")
    plt.title("Aitoff projection of our random data")
    plt.grid(True)
    plt.plot(ra_rad, dec_rad, 'o', markersize=2, alpha=0.5)
    plt.subplots_adjust(top=0.90, bottom=0.0)
    figname = 'ob_progress_aitoff.png'
    print('Saving: ' + figname)
    plt.savefig(figname)

    print('Elapsed time(secs):', time.time() - t0)

    aitoff = True
    PA = 90.0
    color = 'blue'
    alpha = 1.0
    debug = False
    plt.figure(figsize=(8.0, 6.0))
    plt.subplot(111, projection="aitoff")
    plt.title("Aitoff projection of our random data")
    plt.grid(True)
    for i in range(len(ra)):
        if i == 0 or i == len(ra):
            print('i:', i)
        (ra_poly, dec_poly) = plt_tile.mk_polytile(
            ra_cen=ra[i], dec_cen=dec[i],
            coverage='twice', PA=PA)
        if debug:
            print('Elapsed time(secs): ', time.time() - t0)

        if aitoff:
            # convert to wrapped radians
            c = SkyCoord(ra=ra_poly, dec=dec_poly,
                         unit=(u.deg, u.deg), frame='icrs')
            ra_poly = c.ra.wrap_at(180 * u.deg).radian
            dec_poly = c.dec.radian
            if debug:
                print('ra_poly:', i, len(ra_poly), ra_poly)
                print('Elapsed time(secs):', time.time() - t0)

        xypolygon = np.column_stack((ra_poly, dec_poly))
        if i == 0 or i == len(ra):
            print('xypolygon.shape:', xypolygon.shape)
            print('xypolygon:', xypolygon)

        # http://matplotlib.org/api/patches_api.html#matplotlib.patches.Polygon
        polygon = Polygon(
            xypolygon, closed=True,
            edgecolor='none', facecolor=color, alpha=alpha)
        plt.gca().add_patch(polygon)
        # this does nothing
        plt.set_label = 'label test for polygon'

    plt.subplots_adjust(top=0.90, bottom=0.0)
    figname = 'ob_progress_aitoff_polygons.png'
    print('Saving:' + figname)
    plt.savefig('./' + figname)

    print('Elapsed time(secs):', time.time() - t0)

    print('Plotting the shaded tiles')


def plot_status_notcompleted(table, raUnits='hours', wrap_ra24hr=False):
    """


    """
    # repeat for not completed OBs
    print('Number of rows in table:', len(table))
    ra = table['RA (hrs)']
    dec = table['DEC (deg)']

    ExecutionTime = rd_ExecutionTimes(table=table, debug=False)
    OBstatus = rd_OBstatus(table=table, debug=False)
    OBstatusDate = rd_OBstatusDate(table=table, debug=True)

    notcompleted = True
    if notcompleted:
        # ignore cancelled OBs
        itest = (OBstatus != 'C') & (OBstatus != 'K')
        print(len(itest), len(OBstatus[itest]))
        ra = ra[itest]
        dec = dec[itest]
        ExecutionTime = ExecutionTime[itest]
        OBstatus = OBstatus[itest]

    COMPLETED = vhs_select_status(OBstatus)
    NOTCOMPLETED = [not i for i in COMPLETED]

    ATLAS = vhs_select_survey(ExecutionTime, survey='ATLAS')
    DES = vhs_select_survey(ExecutionTime, survey='DES')
    GPS = vhs_select_survey(ExecutionTime, survey='GPS')
    OTHER = vhs_select_survey(ExecutionTime, survey='OTHER')

    total_atlas = len(ExecutionTime[ATLAS])
    total_des = len(ExecutionTime[DES])
    total_gps = len(ExecutionTime[GPS])
    total_obs = len(ExecutionTime[ATLAS]) + \
        len(ExecutionTime[DES]) + len(ExecutionTime[GPS])

    print('Number of ATLAS OBs:   ', total_atlas)
    print('Number of DES OBs:   ', total_des)
    print('Number of GPS OBs:   ', total_gps)
    print('Number of all (ATLAS, DES, GPS) OBs:   ', total_obs)

    color = 'orange'
    print('len(ra[DES]): ', len(ra[DES]))
    if len(ra[DES]) > 0:
        label = 'VHS-DES   (incomplete OBs)'
        plot_vista_tiles(
            table=table,
            ra=ra[DES], dec=dec[DES], aitoff=aitoff,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=1.0, color=color,
            figfile=figfile, title=title, raUnits=raUnits, overplot=True,
            rarange=rarange, decrange=decrange,
            suffix=date)

    if len(ra[ATLAS]) > 0:
        label = 'VHS-ATLAS (incomplete OBs)'
        plot_vista_tiles(
            table=table, ra=ra[ATLAS], dec=dec[ATLAS], aitoff=aitoff,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=0.5, color=color,
            figfile=figfile, title=title, raUnits=raUnits, overplot=True,
            rarange=rarange, decrange=decrange,
            suffix=date)

    if len(ra[GPS]) > 0:
        label = 'VHS-GPS   (incomplete OBs)'
        plot_vista_tiles(
            table=table, ra=ra[GPS], dec=dec[GPS], aitoff=aitoff,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=0.25, color=color,
            figfile=figfile, title=title, raUnits=raUnits, overplot=True,
            rarange=rarange, decrange=decrange,
            suffix=date)

def save_radec(table):
    """

    """
    ra = table['RA (hrs)'] * 15.0
    dec = table['DEC (deg)']

    # col1 = fits.Column(name='RA', format='E', array=ra)
    # col2 = fits.Column(name='DEC', format='E', array=dec)

    outTable = Table()
    outTable['RA'] = ra
    outTable['Dec'] = dec

    print('Savefile for DES non-DE paper')
    print('RA range: ', np.min(ra), np.max(ra))
    print('DEC range:', np.min(dec), np.max(dec))

    plt.figure(figsize=(8.0, 6.0))
    ndata = len(ra)
    plt.plot(ra, dec, '.', label=str(ndata))
    plt.legend()
    figname = appname + '_save_radec_' + datestamp + '.png'
    print('Saving: ' + figname)
    plt.savefig('./' + figname)

    savefile = 'tmp.fits'
    print('Saving: ', savefile)
    print('Number of rows:', ndata)
    outTable.write(savefile, overwrite=True)


def plot_despolygon(table=None,
                    raUnits='hours', wrap_ra24hr=False,
                    label='DES Round13-poly'):
    """

    """
    import traceback
    trace = traceback.extract_stack()[-1]
    print(os.path.basename(trace[0]), ':', str(trace[1]))

    print('plot_despolygon; label: ', label)

    datestamp = time.strftime("%Y%m%d", gmtime())

    # convert RA in degrees to hours
    if raUnits == 'hours':
        xdata = table['ra'] / 15.0
    if raUnits == 'degrees':
        xdata = table['ra']
    ydata = table['dec']

    if label is None:
        label = ''

    plt.plot(xdata, ydata, color='blue', label=label)

    if raUnits == 'hours':
        plt.plot(xdata + 24.0, ydata, color='blue', label=label)

    if wrap_ra24hr and raUnits == 'degrees':
        plt.plot(xdata - 180.0, ydata, color='blue', label=label)

    figname = appname + '_' + datestamp + '.png'
    print('Saving: ' + figname)
    plt.savefig('./' + figname)


def plot_skycoords(system='Galactic', axis='latitude',
                   units='degrees', wrap_ra24hr=False,
                   overplot=False, label=True):
    """
    could be generalised for ecliptic and lat/long

    based on IDL plot_galcoords.pro
    handles 24hr wrap

    """

    print('Plotting ' + system + ' Coordinates')
    print('plot_skycoords; label:', label)

    if not overplot:
        plt.figure(figsize=(8.0, 6.0))

    if system == 'Galactic':
        color = 'red'
        latitudes = [-30.0, -5.0, 0.0, 5.0, 30.0]
        linestyles = ['--', '--', '-', '--', '--']

        # latitudes= [-30.0, -20.0,-10.0, -5.0, 0.0, 5.0, 10.0, 20.0, 30.0]
        # linestyles = [ '--', '--', '--', '--', '-', '--','--', '--','--']

    if system == 'Ecliptic':
        color = 'blue'
        latitudes = [0.0]
        linestyles = ['-.']

        # latitudes= [-30.0, -20.0,-10.0, -5.0, 0.0, 5.0, 10.0, 20.0, 30.0]
        # linestyles = [ '--', '--', '--', '--', '-', '--','--', '--','--']

    npoints = 3600
    npoints = npoints + 1

    ilabel = 0
    ilat = -1
    # could be generalised for longitude and ecliptic coordinates
    for lat in latitudes:
        ilat = ilat + 1
        print('latitude:', lat, latitudes[ilat])

        latitude = np.full((npoints,), lat)
        longitude = np.linspace(0.0, 360.0, num=npoints)

        if system == 'Galactic':

            Galactic = SkyCoord(
                longitude * u.degree,
                latitude * u.degree, frame='galactic')
            Galactic_icrs = Galactic.icrs

            if wrap_ra24hr:
                ra = Galactic_icrs.ra.wrap_at(180 * u.deg).degree
            if not wrap_ra24hr:
                ra = Galactic_icrs.ra.degree
            dec = Galactic_icrs.dec.degree

        if system == 'Ecliptic':
            Ecliptic = SkyCoord(
                longitude * u.degree,
                latitude * u.degree, frame='geocentrictrueecliptic')
            Ecliptic_icrs = Ecliptic.icrs

            if wrap_ra24hr:
                ra = Ecliptic_icrs.ra.wrap_at(180 * u.deg).degree
            if not wrap_ra24hr:
                ra = Ecliptic_icrs.ra.degree

            dec = Ecliptic_icrs.dec.degree

        print('0,0: ', ra[0], dec[0])

        # determine if where there are 24hrs wraps
        if units == 'degrees':
            wrap_range = 240.0
        if units == 'hours':
            wrap_range = 18.0
            ra = ra / 15.0

        i = -1
        nwrap = 0
        iwrap = 0
        wrap = []
        ndata = len(ra)
        for i in range(0, ndata):
            if i <= ndata - 2:
                if abs(ra[i] - ra[i + 1]) >= wrap_range:
                    wrap.append(i)
                    iwrap = iwrap + 1
                    print('wrap:', iwrap, i, ra[i], ra[i + 1])
        nwrap = iwrap

        print('plot_skycoords; label:', label, ilabel)
        if label is not None:
            label = system + ' Coordinates'
        print('plot_skycoords; label:', label, ilabel)

        linestyle = linestyles[ilat]
        if nwrap == 0:
            print('No wrapping: ', latitudes[ilat])

            xdata = ra
            ydata = dec

            ilabel = plot_skycoords_plot(
                xdata, ydata, linestyle=linestyle,
                label=label, ilabel=ilabel, color=color,
                wrap_ra24hr=wrap_ra24hr, units=units)

        if nwrap == 1:
            iwrap = 0
            print()
            print('wrapping:', latitudes[ilat], nwrap)
            print('wrap at: ', wrap)
            print(iwrap, wrap[iwrap], ra[wrap[iwrap]], ra[wrap[iwrap] + 1])
            print(iwrap, wrap[iwrap], dec[wrap[iwrap]], dec[wrap[iwrap] + 1])

            ndata = len(ra)
            xdata = ra[0:wrap[iwrap]]
            ydata = dec[0:wrap[iwrap]]

            ilabel = plot_skycoords_plot(
                xdata, ydata, linestyle=linestyle,
                label=label, ilabel=ilabel, color=color,
                wrap_ra24hr=wrap_ra24hr, units=units)

            xdata = ra[wrap[iwrap] + 1:-1]
            ydata = dec[wrap[iwrap] + 1:-1]

            ilabel = plot_skycoords_plot(
                xdata, ydata, linestyle=linestyle,
                label=label, ilabel=ilabel, color=color,
                wrap_ra24hr=wrap_ra24hr, units=units)

        if nwrap >= 2:
            print()
            print('wrapping: ', latitudes[ilat], nwrap)
            print('wrap at: ', wrap)

            isegment = 0
            for iwrap in range(0, nwrap):
                print()
                print(iwrap, nwrap, wrap[iwrap],
                      ra[wrap[iwrap] - 1], ra[wrap[iwrap]],
                      ra[wrap[iwrap] + 1])
                print(iwrap, nwrap, wrap[iwrap],
                      dec[wrap[iwrap] - 1], dec[wrap[iwrap]],
                      dec[wrap[iwrap] + 1])

                if iwrap == 0:
                    isegment = isegment + 1
                    print()
                    print('segment: ', isegment)
                    print(iwrap, nwrap, 0, '-', wrap[iwrap])
                    print('segment (RA):  ', ra[0], '-', ra[wrap[iwrap]])
                    print('segment (Dec): ', dec[0], '-', dec[wrap[iwrap]])

                    xdata = ra[0:wrap[iwrap]]
                    ydata = dec[0:wrap[iwrap]]

                    ilabel = plot_skycoords_plot(
                        xdata, ydata,
                        linestyle=linestyle, color=color,
                        label=label, ilabel=ilabel,
                        wrap_ra24hr=wrap_ra24hr, units=units)

                if iwrap != 0:
                    isegment = isegment + 1
                    print()
                    print('segment: ', isegment)
                    print(iwrap, nwrap, wrap[iwrap - 1] + 1, '-', wrap[iwrap])
                    print('segment (RA):',
                          ra[wrap[iwrap - 1] + 1], '-', ra[wrap[iwrap]])
                    print('segment (Dec):',
                          dec[wrap[iwrap - 1] + 1], '-', dec[wrap[iwrap]])

                    xdata = ra[wrap[iwrap - 1] + 1:wrap[iwrap]]
                    ydata = dec[wrap[iwrap - 1] + 1:wrap[iwrap]]

                    ilabel = plot_skycoords_plot(
                        xdata, ydata,
                        linestyle=linestyle, color=color,
                        label=label, ilabel=ilabel,
                        wrap_ra24hr=wrap_ra24hr, units=units)

                if iwrap == nwrap - 1:
                    isegment = isegment + 1
                    print()
                    print('segment: ', isegment)
                    print(iwrap, nwrap, wrap[iwrap], npoints - 1)
                    print('segment (RA): ', ra[wrap[iwrap] + 1], '-', ra[-1])
                    print('segment (Dec):', dec[wrap[iwrap] + 1], '-', dec[-1])

                    xdata = ra[wrap[iwrap] + 1:-1]
                    ydata = dec[wrap[iwrap] + 1:-1]

                    ilabel = plot_skycoords_plot(
                        xdata, ydata,
                        linestyle=linestyle, color=color,
                        label=label, ilabel=ilabel,
                        wrap_ra24hr=wrap_ra24hr, units=units)

    plt.legend(fontsize='small')
    plt.grid()
    figname = appname + '_galactic_' + datestamp + '.png'
    print('Saving:' + figname)
    plt.savefig('./' + figname)


def plot_skycoords_plot(xdata, ydata,
                        label=None, ilabel=0,
                        color='red', linestyle=None,
                        wrap_ra24hr=False, units='degrees'):
    """


    """
    print('plot_skycoords_plot; label:', label, ilabel)

    # add plot legend label to first segment in series
    if ilabel == 0:
        plt.plot(xdata, ydata,
                 linestyle=linestyle,
                 color=color, label=label)
        ilabel = 1

    if ilabel != 0:
        plt.plot(xdata, ydata,
                 linestyle=linestyle,
                 color=color)

    if units == 'degrees':
        if wrap_ra24hr:
            plt.xlim([-180.0, 180.0])
        else:
            plt.xlim([0.0, 360.0])

    if not aitoff and units == 'hours':
        if wrap_ra24hr:
            plt.xlim([-12.0, 12.0])
        else:
            plt.xlim([0.0, 24.0])

    if not aitoff:
        plt.ylim([-90.0, +90.0])

    return ilabel


if __name__ == '__main__':

    import os
    import sys

    import doctest
    doctest.testmod()

    print('os.path.basename(sys.argv[0]): ', os.path.basename(sys.argv[0]))
    print('__file__: ', __file__)
    print(os.path.splitext(__file__)[0])
    # strip off file extension
    global appname
    appname = os.path.splitext(__file__)[0]
    print('Appname: ', appname)

    import ConfigParser

    from argparse import ArgumentParser

    t0 = time.time()

    parser = ArgumentParser(
        description='OB progress analysis (see also IDL version ob_progress)'
    )

    date_default = time.strftime("%Y%m%d", gmtime())
    parser.set_defaults(date=date_default)
    parser.add_argument(
        "--date", default=date_default,
        dest='date', help="date as string e.g. '20151101'")



    parser.add_argument(
        "--rarange", default=[0.0, 24.0], type=float, nargs=2,
        help="RA range in hours in form Deg Deg]")

    parser.add_argument(
        "--decrange", default=[-90.0, 30.0], type=float, nargs=2,
        help="Declination range in degrees in form Deg Deg")

    runs_default = 'ALL'
    parser.set_defaults(runs=runs_default)
    parser.add_argument(
        "--runs", default=runs_default,
        dest='runs', help="runs as string list e.g. 'ABC' ")

    survey_default = 'ALL'
    parser.set_defaults(survey=survey_default)
    parser.add_argument(
        "--survey", default=survey_default,
        dest='survey',
        help="survey as string list e.g. 'DES' or 'GPS' or 'ATLAS' ")

    obstatus_default = 'ALL'
    parser.set_defaults(obstatus=obstatus_default)
    parser.add_argument(
        "--obstatus",
        default=obstatus_default,
        help="obstatus as string list e.g 'AC'" +
        " (default is ALL excluding status K = cancelled OBs)")

    qcstatus_default = 'ALL'
    parser.set_defaults(default=qcstatus_default)
    parser.add_argument("--qcstatus",
                        help="qcstatus as string list e.g. 'D' ")

    viking_dqc_default = False
    parser.set_defaults(viking_dqc=viking_dqc_default)
    parser.add_argument("--viking_dqc",
                        action="store_true",
                        help="VIKING DQC file")

    viking_sources_default = False
    parser.set_defaults(viking_sources=viking_sources_default)
    parser.add_argument(
       "--viking_sources",
       action="store_true",
       help="""VIKING source list file (assume VSA format with RA,
               Dec in radians)""")

    # radian or degree; need to support parsing check
    viking_radec_format_default = "degree"
    parser.set_defaults(viking_radec_format=viking_radec_format_default)
    parser.add_argument(
       "--viking_radec_format",
       action="store_true",
       help="""VIKING dqc/source/tile list RA, Dec format (radian or degree)""")

    wrap_ra24hr_default = False
    parser.set_defaults(wrap_ra24hr=wrap_ra24hr_default)
    parser.add_argument("--wrap_ra24hr", action="store_true",
                        dest='wrap_ra24hr',
                        help="turn on 24hr RA wrap")

    aitoff_default = False
    parser.set_defaults(aitoff=aitoff_default)
    parser.add_argument("--aitoff", action="store_true",
                        dest='aitoff',
                        help="optional aitoff aspect")

    parser.add_argument("--verbose", action="store_true",
                        dest='verbose',
                        help="optional verbose mode")


    parser.add_argument("--debug", action="store_true",
                        dest='debug',
                        help="optional debug very verbose mode")

    # parser.add_argument('nums', nargs=2)
    parser.add_argument("--pause", action="store_true",
                        dest='pause', help="turn on pausing option")

    # parser.add_argument("-t", "--test", action="store_true", dest="test",
    #              default="", help="run tests")

    parser.add_argument("-v", "--version", action="store_true",
                        dest="version",
                        default="", help="print version number and  exit")

    print('Number of arguments:', len(sys.argv), 'arguments: ', sys.argv[0])
    args = parser.parse_args()

    if args.version:
        print('version:', __version__)
        sys.exit(0)

    pause = False
    if args.pause:
        pause = True

    iarg = 0
    for arg in vars(args):
        iarg = iarg + 1
        print('Argument:', str(iarg) + ':', arg, '=', getattr(args, arg))

    # only supports 1 or all for now until I add parsing of string like
    # in IDL version
    if args.runs == 'ALL':
        run = ''
    if args.runs != 'ALL':
        run = args.runs
    print('Runs: ', run)

    if args.obstatus == 'ALL':
        obstatus = ''
    if args.obstatus != 'ALL':
        obstatus = args.obstatus
    print('OBstatus: ', obstatus)

    if args.qcstatus == 'ALL':
        qcstatus = ''
    if args.qcstatus != 'ALL':
        qcstatus = args.qcstatus
    print('QCstatus: ', qcstatus)

    survey = args.survey
    print('survey: ', survey)

    # print('args.wrap_ra24hr:', args.wrap_ra24hr)
    wrap_ra24hr = args.wrap_ra24hr
    print('wrap_ra24hr: ', wrap_ra24hr)

    raUnits = 'hours'
    print('raUnits:', raUnits)

    print('args.date:', args.date)
    date = args.date

    rarange = args.rarange
    print('RA range:', rarange)

    decrange = args.decrange
    print('Dec range:', decrange)


    # no changes should be needed below

    viking_radec_format = args.viking_radec_format

    VIKING = False
    viking_dqc = False
    if args.viking_dqc:
        VIKING = True
        viking_dqc = True

    viking_sources = args.viking_sources
    if args.viking_sources:
        VIKING = True

    print('viking_radec_format:', viking_radec_format)

    zoom = True

    outpath = './'

    aitoff = args.aitoff
    debug = args.debug
    verbose = args.verbose

    if pause:
        raw_input("Enter any key to continue: ")

    datestamp = time.strftime("%Y%m%d", gmtime())

    progid = '179A2010'
    infile = progresspath + date + '/' + progid + run + '.fits'

    print('Input progress file:', infile)
    figpath = outpath

    if pause:
        raw_input("Enter any key to continue: ")

    print('Reading: ', infile)
    if not os.path.exists(infile):
        print(infile, 'does not exist; exiting')
        print('Check path and date of most recent progress file')
        sys.exit(0)

    table = Table.read(infile)
    table.meta['filename'] = infile
    print('table.colnames: ', table.colnames)
    print('Number of rows read in: ', len(table))

    desfile = '/home/rgm/soft/des/round13-poly.txt'
    print('Reading: ', desfile)
    despolygon = des_footprint.rd_despolygon(desfile)

    if pause:
        raw_input("Enter any key to continue: ")

    if VIKING:
        if viking_dqc:
            viking = Table.read(dqcfile_viking)
            viking.meta['filename'] = dqcfile_viking
            print('Table column names: \n', viking.colnames)
            print('Table metadata: \n', viking.meta)
            viking.pprint

        if viking_sources:
            viking = Table.read(sourcefile_viking)
            viking.meta['filename'] = dqcfile_viking
            print('Table column names: \n', viking.colnames)
            print('Table metadata: \n', viking.meta)
            viking.pprint


    ra = table['RA (hrs)']
    dec = table['DEC (deg)']

    ExecutionTime = rd_ExecutionTimes(table=table, debug=False)
    OBstatus = rd_OBstatus(table=table, debug=False)
    OBstatusDate = rd_OBstatusDate(table=table, debug=True)

    # not used yet
    COMPLETED = vhs_select_status(OBstatus, status='C')
    CANCELLED = vhs_select_status(OBstatus, status='K')
    # print(COMPLETED)
    print('Total number of OBs: ', len(table))
    print('Number of Completed OBs: ', len(table[COMPLETED]))
    NOTCOMPLETED = [not i for i in COMPLETED]
    NOTCOMPLETED = ~COMPLETED
    print('Number of NOT Completed OBs: ', len(table[NOTCOMPLETED]))
    print('Number of Cancelled(K)  OBs: ', len(table[CANCELLED]))
    NOTCOMPLETED = NOTCOMPLETED & ~CANCELLED
    print('Number of NOT Completed OBs: ', len(table[NOTCOMPLETED]))

    completed = False
    if completed:
        itest = (OBstatus == 'C')
        print('Completed OBs: ')
        print(len(itest), len(OBstatus[itest]))
        ra = ra[itest]
        dec = dec[itest]
        ExecutionTime = ExecutionTime[itest]
        OBstatus = OBstatus[itest]
        OBstatusDate = OBstatusDate[itest]

    ATLAS = vhs_select_survey(ExecutionTime, survey='ATLAS')
    DES = vhs_select_survey(ExecutionTime, survey='DES')
    GPS = vhs_select_survey(ExecutionTime, survey='GPS')
    OTHER = vhs_select_survey(ExecutionTime, survey='OTHER')

    total_obs = (len(ExecutionTime[ATLAS]) +
                 len(ExecutionTime[DES]) + len(ExecutionTime[GPS]))
    print('Number of all (ATLAS, DES, GPS) OBs:', total_obs)

    print('DES; All:', len(table[DES]))
    DES = DES & COMPLETED
    print('DES; Completed:', len(table[DES]))

    ATLAS = ATLAS & COMPLETED
    GPS = GPS & COMPLETED

    total_obs = (len(ExecutionTime[ATLAS]) +
                 len(ExecutionTime[DES]) + len(ExecutionTime[GPS]))
    print('Number of completed (ATLAS, DES, GPS) OBs:', total_obs)

    # save for DES non-DE paper
    savefile = False
    if savefile:
        save_radec(table[DES])
        raw_input("Enter any key to continue: ")

    figfile = figpath + '/' + 'ob_progress_radec_des_' + datestamp + '.png'
    title = 'VHS-DES Progress: ' + infile

    plot_radec(ra[DES], dec[DES], title=title,
               plotfile=figfile,
               rarange=rarange, decrange=decrange)

    print('Elapsed time(secs):', time.time() - t0)

    # aitoff_test(ra, dec)
    # raw_input("Enter any key to continue: ")

    # option to centre on 0 hrs with 24hr wrap

    if raUnits == 'hours':
        if wrap_ra24hr:
            rarange = [-12.0, 12.0]

    print()
    print('Plotting VHS tiles')
    print('RA range: ', np.min(ra), np.max(ra))
    print('Dec range:', np.min(dec), np.max(dec))

    overplot = False

    # plot_skycoords(overplot=False, wrap_ra24hr=True)
    # plot_skycoords(overplot=True, system='Ecliptic', wrap_ra24hr=True)
    # if pause: raw_input("Enter any key to continue: ")

    plot_skycoords(overplot=False,
                   wrap_ra24hr=True, units='hours', label=None)
    plot_skycoords(overplot=True, system='Ecliptic', wrap_ra24hr=True,
                   units='hours', label=None)
    if pause:
        raw_input("Enter any key to continue: ")

    # inefficent way to get the total into the legend via blank plot
    overplot = False
    label = 'VHS ALL   (completed OBs)'
    plot_vista_tiles(table=table,
                     ra=ra[COMPLETED], dec=dec[COMPLETED],
                     aitoff=aitoff, overplot=overplot,
                     verbose=verbose, debug=debug,
                     wrap_ra24hr=wrap_ra24hr, label=label,
                     alpha=1.0, color='none',
                     figfile=figfile, title=title, raUnits=raUnits,
                     rarange=rarange, decrange=decrange,
                     suffix=date)
    overplot = True

    raUnits = 'hours'
    plot_skycoords(overplot=overplot, wrap_ra24hr=wrap_ra24hr,
                   units=raUnits,
                   label=None)
    plot_skycoords(overplot=overplot, system='Ecliptic',
                   wrap_ra24hr=wrap_ra24hr,
                   units=raUnits, label=None)

    plot_despolygon(table=despolygon, label=None, raUnits=raUnits,
                    wrap_ra24hr=wrap_ra24hr)
    plt.grid()

    if pause:
        raw_input("Enter any key to continue: ")

    if survey == 'ALL' or survey == 'DES':
        label = 'VHS-DES   (completed OBs)'
        print('Plot completed VHS-DES OBs')
        plot_vista_tiles(
            table=table, ra=ra[DES], dec=dec[DES],
            aitoff=aitoff, overplot=overplot,
            verbose=verbose, debug=debug,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=1.0, color='green',
            figfile=figfile, title=title, raUnits=raUnits,
            rarange=rarange, decrange=decrange,
            suffix=date)
        overplot = True

    if pause:
        raw_input("Enter any key to continue: ")
    if survey == 'ALL' or survey == 'ATLAS':
        label = 'VHS-ATLAS (completed OBs)'
        print('Plot completed VHS-ATLAS OBs')
        plot_vista_tiles(
            table=table, ra=ra[ATLAS], dec=dec[ATLAS],
            aitoff=aitoff, overplot=overplot,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=0.5, color='green',
            figfile=figfile, title=title, raUnits=raUnits,
            rarange=rarange, decrange=decrange,
            suffix=date)
        overplot = True

    if pause:
        raw_input("Enter any key to continue: ")

    if survey == 'ALL' or survey == 'GPS':
        label = 'VHS-GPS   (completed OBs)'
        print('Plot completed VHS-ATLAS OBs')
        plot_vista_tiles(
            table=table, ra=ra[GPS], dec=dec[GPS],
            aitoff=aitoff, overplot=overplot,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=0.25, color='green',
            figfile=figfile, title=title, raUnits=raUnits,
            rarange=rarange, decrange=decrange,
            suffix=date)
        overplot = True

    if pause:
        raw_input("Enter any key to continue: ")

    # plot the incomplete OBs
    if obstatus != 'C':
        # inefficent? way to get the total into the legend via blank plot
        label = 'VHS ALL   (incomplete OBs)'
        plot_vista_tiles(
            table=table,
            ra=ra[NOTCOMPLETED], dec=dec[NOTCOMPLETED],
            aitoff=aitoff, overplot=overplot,
            verbose=verbose, debug=debug,
            wrap_ra24hr=wrap_ra24hr, label=label,
            alpha=1.0, color='none',
            figfile=figfile, title=title, raUnits=raUnits,
            rarange=rarange, decrange=decrange,
            suffix=date)

        overplot = True
        print('Number of tiles:', len(table))
        plot_status_notcompleted(table,
                                 raUnits=raUnits,
                                 wrap_ra24hr=wrap_ra24hr)

    if viking_dqc:
        xdata = viking['ra']
        ydata = viking['dec']
        angle = Angle(xdata * u.deg)
        angle.wrap_at('180d', inplace=True)
        xdata = angle.degree / 15.0

        # plt.plot(xdata, ydata, 'sr',
        #    ms=5.0, markeredgecolor='r', alpha=0.2, label='VIKING')

        print('Plotting VIKING tiles')
        plt.suptitle(dqcfile_viking)
        print('RA range:  ', np.min(viking['ra'] / 15.0),
              np.max(viking['ra'] / 15.0))
        print('Dec range:', np.min(viking['dec']), np.max(viking['dec']))
        label = 'VIKING   (completed OBs)'
        plot_vista_tiles(table=table,
                         ra=viking['ra'] / 15.0, dec=viking['dec'],
                         wrap_ra24hr=wrap_ra24hr, label=None,
                         PA=0.0, aitoff=aitoff,
                         alpha=0.1, color='red',
                         figfile=figfile, title=None,
                         raUnits=raUnits, overplot=True,
                         rarange=rarange, decrange=decrange,
                         suffix=date)

    if viking_sources:
        # VSA format
        viking.info('stats')
        xdata = viking['RA']
        ydata = viking['DEC']
        # angle = Angle(xdata * u.rad)
        # angle.wrap_at('180d', inplace=True)
        xdata = np.degrees(xdata) / 15.0
        ydata = np.degrees(ydata)

        # plt.plot(xdata, ydata, 'sr',
        #    ms=5.0, markeredgecolor='r', alpha=0.2, label='VIKING')

        print('Plotting VIKING sources')
        plt.suptitle(sourcefile_viking)
        print('RA range:  ', np.min(xdata), np.max(xdata))
        print('Dec range:', np.min(ydata), np.max(ydata))
        label = 'VIKING   (completed OBs)'
        plt.plot(xdata, ydata, 'r.')
        plt.grid()

        figfile = figpath + '/' + 'ob_progress_des_viking_' + datestamp + '.png'
        print('Saving: ' + figfile)
        plt.savefig(figfile)

    print('Completed plotting the shaded tiles')
    if pause:
        raw_input("Enter any key to continue: ")

    figfile = figpath + '/' + 'ob_progress_radec_atlas_' + datestamp + '.png'
    plot_radec(ra[ATLAS], dec[ATLAS], title='VHS-GPS Progress: ' + infile,
               plotfile=figfile,
               rarange=rarange, decrange=decrange)

    figfile = figpath + '/' + 'ob_progress_radec_gps_' + datestamp + '.png'
    plot_radec(ra[GPS], dec[GPS], title='VHS-GPS Progress: ' + infile,
               plotfile=figfile,
               rarange=rarange, decrange=decrange)

    figfile = outpath + '/' + 'ob_progress_ra_executiontime_' + \
        datestamp + '.png'

    print('wrap_ra24hr:', wrap_ra24hr)
    print('raUnits:', raUnits)

    # plot_ra_extime(ra, executionTime, title='VHS Progress: ' + fitsfile,
    # figfile=figfile,
    # rarange=[0.0, 24.0])
