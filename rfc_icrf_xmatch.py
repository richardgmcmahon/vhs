
import sys

import numpy as np

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy import units as u

from astropy.coordinates import SkyCoord

from astroquery.vizier import Vizier
from astroquery.xmatch import XMatch

from astropy.stats import mad_std


def getconfig(configfile=None, debug=False, silent=False):
    """
    read config file

    Note the Python 2 ConfigParser module has been renamed to configparser
    in Python 3 so it better to use import configparser in Python 2 for
    future proofing

    see also getconfig.cfg

    TODO: catch exceptions

    Support for lists:

    see:

    https://stackoverflow.com/questions/335695/lists-in-configparser

    https://github.com/cacois/python-configparser-examples

    look in cwd, home and home/.config

    home/.config not implemented yet


    """
    import os
    import configparser

    # read the configuration file
    # config = configparser.RawConfigParser()
    config = configparser.SafeConfigParser()

    print('__file__', __file__)
    print('configfile:', configfile)
    configfile_default = os.path.splitext(__file__)[0] + '.cfg'
    print('configfile_default:', configfile_default)

    if configfile is None:
        configfile_default = os.path.splitext(__file__)[0] + '.cfg'
        if debug:
            print('__file__', __file__)
            print('configfile_default:', configfile_default)
        configfile = configfile_default

    print('Open configfile:', configfile)
    if debug:
        print('Open configfile:', configfile)

    try:
        if not silent:
            print('Reading config file', configfile)

        try:
            config.read(configfile)
        except IOError:
            print('config file', configfile, "does not exist")
            configfile = os.path.join(os.environ["HOME"], configfile)
            print('trying ', configfile)
            config.read(configfile)

    except Exception as e:
        print('Problem reading config file: ', configfile)
        print(e)

    if debug:
        print('configfile:', configfile)
        print('sections:', config.sections())
        for section_name in config.sections():
            print('Section:', section_name)
            print('Options:', config.options(section_name))
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
        print()

        for section_name in config.sections():
            print()
            print('Section:', section_name)
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
                print(section_name, ':',
                      name, config.get(section_name, name))
        print()

    return config


def getargs(verbose=False):
    """

    Template getargs function

    Usage

    python getargs.py --help


    def getargs():

    ....

    if __name__=='__main__':

        args = getargs()
        debug = args.debug()



    parse command line arguements

    not all args are active

    """
    import sys
    import pprint
    import argparse

    # there is probably a version function out there
    __version__ = '0.1'

    description = 'This is a template using getargs'
    epilog = """WARNING: Not all options may be supported
             """
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # the destination defaults to the option parameter
    # defaul=False might not be needed

    # default type is string
    parser.add_argument("--survey",
                        default = '2MASS',
                        help="Input survey")

    parser.add_argument("--float", type=float,
                        help="float input")

    parser.add_argument("--configfile",
                        default=None,
                        help="configuration file")

    parser.add_argument("--infile",
                        help="Input file name")

    parser.add_argument("--debug",
                        action='store_true',
                        help="debug option")

    parser.add_argument("--pause",
                        action='store_true',
                        help="Pausing option")

    parser.add_argument("--verbose", default=verbose,
                        action='store_true',
                        help="Verbose option")

    parser.add_argument("--version",
                        action='store_true',
                        help="verbose option")

    args = parser.parse_args()

    if args.debug or args.verbose:
        print()
        print('Number of arguments:', len(sys.argv),
              'arguments: ', sys.argv[0])

    if args.debug or args.verbose:
        print()
        for arg in vars(args):
            print(arg, getattr(args, arg))
        print()

    if args.debug or args.verbose:
        print()
        pprint.pprint(args)
        print()

    if args.version:
        print('version:', __version__)
        sys.exit(0)

    return args


if __name__ == '__main__':

    import time

    t0 = time.time()

    args = getargs(verbose=True)
    debug = args.debug
    survey_xmatch = args.survey

    # configfile defaults to
    configfile = 'rfc_icrf_xmatch.cfg'
    # config = getconfig(configfile=configfile, debug=debug)
    config = getconfig(debug=debug)
    infile = '/data/vhs/icrf/rfc_2019a_cat.radec'
    path_rfc = config.get('INPUTS', 'path_rfc')
    filename_rfc = config.get('INPUTS', 'filename_rfc')
    infile = path_rfc + filename_rfc


    # change the Vizier timeout in seconds
    print(Vizier.TIMEOUT)
    Vizier.TIMEOUT = 300
    print(Vizier.TIMEOUT)


    # help(XMatch)
    # help(XMatch.query)

    cat1 = Table.read(infile, format='ascii')
    cat1.info('stats')
    #cat1.rename_column('col1', 'ra')
    #cat1.rename_column('col2', 'dec')
    #cat1.rename_column('col3', 'id')

    cat1.info()

    title_prefix = 'VHSDR5'
    title_prefix = 'rfc_2019a'

    # xmatch_cat = '2MASS'
    # xmatch_cat = 'GAIA'
    # xmatch_cat = 'PS1'

    radec_colnames = ['ra', 'dec']

    # 2MASS
    if survey_xmatch == '2MASS':
        cat2 = 'vizier:' + 'II/246/out'
        radec_colnames = ['RAJ2000', 'DEJ2000']

    if survey_xmatch == 'UKIDSS':
        # UKIDSS DR8
        cat2 = 'vizier:' + 'II/314/las8'

    # VIKING
    if survey_xmatch == 'VIKING':
        cat2 = 'vizier:' + 'II/343/viking2'
        radec_colnames = ['RAdeg', 'DEdeg']

    # Gaia
    if survey_xmatch == 'GAIA':
        cat2 = 'vizier:I/345/gaia2'
        radec_colnames = ['ra_epoch2000', 'dec_epoch2000']

    # PS1
    if survey_xmatch == 'PS1':
        cat2 = 'vizier:' + 'II/349/ps1'
        radec_colnames = ['RAJ2000', 'DEJ2000']

    # VHS
    if survey_xmatch == 'VHS':
        infile2 = '/data/vhs/icrf/' + 'rfc_2019a_vhs_dr5eso.fits'
        radec_colnames = ['RA', 'DEC']

    #ra2 = result['ra']
    #dec2 = result['dec']


    #cat2 = 'vizier:' + 'II/294/sdss7'
    #cat2 = 'vizier:' + 'II/306/sdss8'
    #cat2 = 'vizier:' + 'V/139/sdss9'

    xmatch_catalogue = survey_xmatch
    title = title_prefix + ': ' + xmatch_catalogue
    print(survey_xmatch, title, radec_colnames)

    #table = XMatch.query(cat1=cat1,
    #                     cat2=cat2,
    #                     max_distance=5 * u.arcsec,
    #                     colRA1='col1',
    #                     colDec1='col2')
    #table.info()
    #sys.exit()


    cat1.rename_column('col1', 'ra_rfc')
    cat1.rename_column('col2', 'dec_rfc')
    cat1.rename_column('col3', 'id_rfc')
    cat1.info()


    if survey_xmatch == 'VHS':
        result = Table.read(infile2)
        itest = np.unique(result['UPLOAD_ID'])

    if survey_xmatch != 'VHS':
        print('cat1:', cat1)
        print('Vizier xmatch:', cat2)
        result = XMatch.query(cat1=cat1,
                              cat2=cat2,
                              max_distance=1.0 * u.arcsec,
                              colRA1='ra_rfc',
                              colDec1='dec_rfc')
        itest = np.unique(result['id_rfc'])

    print(len(itest), len(result))
    print('Elapsed time(secs): ',time.time() - t0)
    print(type(result), len(result))
    result.info()
    result.info('stats')

    #sys.exit()

    if survey_xmatch == 'VHS':
        ra_rfc = result['UPLOAD_RA']
        dec_rfc = result['UPLOAD_DEC']
        itest = result['SOURCEID'] > 0
        ra_rfc = ra_rfc[itest]
        dec_rfc = dec_rfc[itest]

    if survey_xmatch != 'VHS':
        ra_rfc = result['ra_rfc']
        dec_rfc = result['dec_rfc']

    ra2 = result[radec_colnames[0]]
    dec2 = result[radec_colnames[1]]

    radec1 = SkyCoord(ra=ra_rfc*u.degree, dec=dec_rfc*u.degree)

    print(ra2.unit, dec2.unit)
    if survey_xmatch == 'VHS':
        ra2 = ra2[itest]
        dec2 = dec2[itest]
        ra2.unit = 'radian'
        dec2.unit = 'radian'
        print(ra2.unit, dec2.unit)
        radec2 = SkyCoord(ra=ra2, dec=dec2)

    if survey_xmatch != 'VHS':
        radec2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)

    sep = radec1.separation(radec2).arcsec
    print(np.median(sep))

    dra, ddec = radec1.spherical_offsets_to(radec2)
    print(np.median(dra.arcsec), mad_std(dra.arcsec),
          np.min(dra.arcsec), np.max(dra.arcsec))
    print(np.median(ddec.arcsec), mad_std(ddec.arcsec),
          np.min(ddec.arcsec), np.max(ddec.arcsec))


    plt.figure(figsize=(6, 6))

    xdata = dra.arcsec
    ydata = ddec.arcsec
    ndata = len(xdata)
    plt.scatter(xdata, ydata, s=1, label=str(ndata))
    plt.legend(loc='upper left')

    ax = plt.gca()
    dra_median = np.median(dra.arcsec)
    dra_mad_std = mad_std(dra.arcsec)
    ddec_median = np.median(ddec.arcsec)
    ddec_mad_std = mad_std(ddec.arcsec)

    text = 'Median dRA: ' + '{:7.4f}'.format(dra_median) + '"'
    text = text + '\n' + 'Sigma(MAD) dRa: ' + \
           '{:7.4f}'.format(dra_mad_std) + '"'
    text = text + '\n' + 'Median dDec: ' + \
           '{:7.4f}'.format(ddec_median) + '"'
    text = text + '\n' + 'Sigma(MAD) dRa: ' + \
           '{:7.4f}'.format(ddec_mad_std) + '"'

    plt.text(0.95, 0.05,
             text,
             verticalalignment='bottom', horizontalalignment='right',
             transform=ax.transAxes,
             color='black', fontsize=10,
             bbox=dict(facecolor='white',edgecolor='green',
                       boxstyle='square'))


    plt.title(title)
    plt.xlabel('Delta RA (arcsec)')
    plt.ylabel('Delta Dec (arcsec)')
    plt.axes().set_aspect('equal')
    # plt.axes().set_aspect('equal', 'datalim')

    plt.xlim(-0.5, 0.5)
    plt.ylim(-0.5, 0.5)
    plt.grid()


    xmatch_catalogue = xmatch_catalogue.replace('/', '_')
    plotfile = title_prefix + '_' + xmatch_catalogue + '.png'
    # remove the protected chaacters
    plotfile = plotfile.replace('/', '_')
    plotfile = plotfile.replace(':', '_')
    print('Saving:', plotfile)
    plt.savefig(plotfile)

    plt.show()
