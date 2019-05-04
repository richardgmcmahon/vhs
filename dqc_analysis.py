from __future__ import print_function, division

__version__="v0.0.1"

"""

TODO:

get filtername from VSA query and ESO program ID


"""

import os
import sys
import time
import traceback
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

import numpy as np

from scipy import stats
import astroML.stats as aml

import astropy
from astropy import coordinates as coord
from astropy import units as u
#from astropy.coordinates import ICRSCoordinates
from astropy.io import ascii

sys.path.append('/home/rgm/soft/python/lib/librgm')
from table_stats import *

import libdqc as libdqc

#import pyfits as pyfits
from astropy.io import fits as pyfits

now = time.localtime(time.time())
print('Current time: ',time.strftime("%Y-%m-%d %H:%M:%S %Z", now))
date=time.strftime("%Y%m%d", now)
print('day: ',date)
print('Current working directory: ',os.getcwd())
print( 'Executing: ',sys.argv[0])

now = time.localtime(time.time())
timestamp = time.strftime("%Y-%m-%dT%H:%M:%S",now)
datestamp = time.strftime("%Y%m%d",now)

print('timestamp: ', timestamp)
print('datestamp: ', datestamp)


if __name__=='__main__':

  print('Executing: ',sys.argv[0], __version__)

  trace = traceback.extract_stack()[-1]
  print(os.path.basename(trace[0]), ':', str(trace[1]))

  t0=time.time()

  debug=False
  verbose=False


  inpath= '/data/vhs/dqc/vsa/2014/VHSv20140409/'
  filename= 'vhs_vsa_dqc_tiles_fs_metadata.fits'

  dqcfile= inpath + filename

  print('Elpased time: %.2f seconds' %(time.time() - t0))

  print('Get table stats: ', dqcfile)

  dqc = libdqc.rd_dqc(infile=dqcfile)

  table_stats(dqcfile)

  print('Elpased time: %.2f seconds' %(time.time() - t0))


  #libdqc.dqc_stats(dqc=dqc)

  test=False
  if test:
    test = (data['nite'] >= '20120901') & (data['nite'] <= '20130401')
    data=data[test]


  interactive=True
  showplots=True

  t0=time.time()

  wavebands=['Y','J','H','Ks']

  colparam='NSOURCES'

  nsources_total=np.sum(dqc['nsources'])
  print('Total number of bandmerged sources: ', nsources_total)

  wavebands=['Y','J','H','Ks']
  for waveband in wavebands:
    nsources=np.sum(dqc[waveband+'nsources'])
    print('Total number of bandmerged sources(' + waveband + '): ', nsources)

  nsources = np.array(
   [dqc['Ynsources'], dqc['Jnsources'],
    dqc['Hnsources'],dqc['Ksnsources']]).flatten()

  nsources_max=max(nsources)
  xlimit_max=nsources_max

  libdqc.plot_byband(data=dqc, wavebands=wavebands, colparam=colparam,
   showplots=True, xlabel=colparam, filename=dqcfile,
   masklimit=0, xlimit_min=0, xlimit_max=xlimit_max, xscale='log')

  debug=True
  if debug: raw_input("Press ENTER to continue: ")

  ra=dqc['ra']
  dec=dqc['dec']
  libdqc.plot_radec(ra, dec, title=None, xlabel=None, ylabel=None,
   rarange=None, decrange=None, showplots=False, figfile=None)

  mask=dqc['KSSEEINGARCSECS'] > 0
  print('Unmasked: ', len(dqc), len(dqc[mask]))
  mask=dqc['KSSEEINGARCSECS'] < 3.0
  print('Unmasked: ', len(dqc), len(dqc[mask]))

  colparam='SEEINGARCSECS'
  wavebands=['Y','J','H','Ks']
  libdqc.plot_byband(data=dqc, wavebands=wavebands, colparam=colparam,
   masklimit=0.0,
   showplots=True, xlabel=colparam, filename=dqcfile)

  colparam='_DEPTH_DYE2006'
  libdqc.plot_byband(data=dqc, wavebands=wavebands, colparam=colparam,
    masklimit=0.0,
    showplots=True, xlabel=colparam, filename=dqcfile)


  #colnames=['KSSEEINGARCSECS']
  #libdqc.plot_byband(data=dqc[colnames], colnames=colnames,
  # showplots=True, filename=dqcfile)


  ksmaglimit=libdqc.get_maglimit(data=dqc, waveband='KS')

  libdqc.plot_cdf(data=ksmaglimit,
   showplots=True, filename=dqcfile, plotlabel="ksmaglimit")

  plt.figure(num=None, figsize=(8.0, 8.0))

  xdata=ksmaglimit
  ydata=xdata-dqc['KS_DEPTH_DYE2006']

  ms=1
  plt.plot(xdata, ydata, 'ob', markeredgecolor='b', ms=ms)

  plt.show()
