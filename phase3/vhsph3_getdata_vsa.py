"""

Download VHS bandmerged catalogue from VSA so that they can be delivered for
ESO Phase3

horus is the public website
djer is a test site

http://horus.roe.ac.uk:8080/vdfs/VSQL_form.jsp

Reads an input file and submits a query to the WSA archive.

Based on wsa_freeform_sql_rgm.py

Reads a vhsMergelog DQC file from VSA and is used to
generate generate the Phase3 files from the VSA


TODO

add logic that config file overrides the option defaults
and  that optional inputs overside the config file.

write a logfile of the inputs used.

support for fits, etc on command line
unzip the fits

do not use email address since then the query will background and
you loose state.

20120708; rgm original by rgm


TODO
option to specify the output filename
option to specify fits or csv on the command line


"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import sys
import re
import string
import time
import logging

import urllib
import urllib2
import cookielib

import argparse
import ConfigParser

import astropy
print('astropy.__version__: ', astropy.__version__)
from astropy.table import Table

# sys.path.append('/home/rgm/soft/libdqc/')
sys.path.append('/home/rgm/soft/python/lib/ssm/')
import AstroUtils

sys.path.append('/home/rgm/soft/python/lib/librgm/')
import table_stats

logging.basicConfig(format='%(asctime)s %(message)s',
                    datefmt='%m-%d-%Y %I:%M:%S %p:',
                    level=logging.INFO)


now = time.localtime(time.time())
print('Current time:', time.strftime("%Y-%m-%d %H:%M:%S %Z", now))
date = time.strftime("%Y%m%d", now)
print('day:', date)
print('Current working directory:', os.getcwd())
print('Executing:', sys.argv[0])

logging.info("Timestamp logging info")

# The optional fourth argument is a dict with members that will take
# precedence in interpolation.
# print config.get('VDFS', 'foo', 0, {'bar': 'Documentation',
#                                        'baz': 'evil'})
# Set the third, optional argument of get to 1 if you wish to use raw mode.

config = ConfigParser.ConfigParser()
config_vsa_login = ConfigParser.ConfigParser()

configfile_vsa_login = 'vhs_vsa_login.cfg'
config_vsa_login.read(configfile_vsa_login)

username = config_vsa_login.get('DEFAULT', 'username')
password = config_vsa_login.get('DEFAULT', 'password')
community = config_vsa_login.get('DEFAULT', 'community')


configfile = 'vhsph3_getdata_vsa.cfg'
config.read(configfile)

host = config.get('VSA-VHS', 'host')
archive = config.get('VSA-VHS', 'archive')


# outpath = '/data/vhs/vsa/VHSv20120417/tmp/'
# outpath = '/data/vhs/vsa/VHSv20120417/phase3_dr2_tmp1/'
# outpath = '/data/vhs/vsa/VHSv20120417/phase3_dr2_tmp2/'
# outpath = '/data/vhs/vsa/VHSv20120417/phase3_dr2_v1/'
# outpath = '/data/vhs/vsa/VHSv20120417/phase3_dr2_v2/'
# outpath = '/data/vhs/vsa/VHSv20140517/phase3_dr2_v1/'
# outpath = '/data/vhsardata/vsa/VHSv20160507/phase3_dr4_v1/'

outpath = config.get('DEFAULT', 'outpath')
dqcpath = config.get('DEFAULT', 'dqcpath')
dqcfilename = config.get('DEFAULT', 'dqcfilename')

description = '''Download VHS data tile by tile from VSA'''
epilog = '''Work in progress'''

parser = argparse.ArgumentParser(
    description=description, epilog=epilog,
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument("-t", "--test",
                    action="store_true", default=False,
                    help="run test")

parser.add_argument("--debug",
                    action="store_true", default=False,
                    help="run in debug mode")

parser.add_argument("-l", "--limit",
                    action="store_true", default=False,
                    help="limit SQL to 100 rows with TOP ")

parser.add_argument("-i", "--input", help="input metadata file")

parser.add_argument("-q", "--query",
                    help="specify query on the command line")

parser.add_argument("-f", "--format", default='.fits.gz',
                    help="Output format (csv/fits/fits.gz/html)")

parser.add_argument("-o", "--output", help="Output filename")

parser.add_argument("--config", help="Specify the config  file",
                    default='./')

parser.add_argument("--outpath", help="Output path", default='./')

parser.add_argument(
    "-e", "--email", default='',
    help="Email address for results of long queries; should not be used when batch processing is required")

# args.db = 'VHSv20120926'
# args.db = 'VHSv20130417'
# args.db = 'VHSDR2'
# args.db = 'VHSv20140517'
parser.set_defaults(db='VHSv20160507')
parser.add_argument("-d", "--db", dest="db",
                    help="VSA Database to use")

default_release = 'VHSDR2'
parser.add_argument("--dr", dest="release", default=default_release,
                    help="VSA Data release to use")

args = parser.parse_args()

debug = args.debug

outpath = args.outpath

if not args.test:
    dqcpath = '/data/vhs/vsa/VHSv20120417/'
    dqcfilename = 'vsa_fs_VHSv20120417.fits'

    # VHS_DR3
    dqcpath = '/data/vhs/vsa/dqc/VHSv20140517/'
    dqcfilename = 'vhs_vsa_dqc_tiles_fs_metadata.fits'

    # VHS_DR4
    dqcpath = '/data/vhs/dqc/vsa/2016/VHSv20160507/'
    dqcfilename = 'vhs_vsa_dqc_tiles_fs_metadata.fits'

    args.dqcfile = dqcpath + dqcfilename
    print('dqcfile:', args.dqcfile)

    dqc = Table.read(args.dqcfile)
    print('Number of rows in table:', len(dqc))
    dqc.info()
    dqc.info('stats')

    fsids = dqc['FRAMESETID']
    nsources = dqc['NSOURCES']

    print('outpath:', outpath)
    if not os.path.isdir(outpath):
        print('outpath:', outpath)
        os.mkdir(outpath)

    framesetid = '472446402561'
    framesetid = fsids[0]

t0 = time.time()

if args.test:
    print('Run test query and download')
    fsids = [472446402561]
    outpath = './'


def doquery(framesetid, debug=False):
    """
    Run SQL query on VSA and grab the results file.


    """

    t0 = time.time()

    outfile = outpath + 'fs' + str(framesetid) + '.fits'

    if os.path.exists(outfile):
        print('Skipping since output file already exists:', outfile)
        return

    lockfile = outpath + 'fs' + str(framesetid) + '.lock'

    if os.path.exists(lockfile):
        print('Skipping since lockfile already exists:', lockfile)
        return

    # create lockfile
    # create path if it does not exist
    # you might need to catch exception you cannot create the lockfile
    if not os.path.exists(outpath):
        try:
            os.makedirs(path)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
                return

    print('Create lockfile:', lockfile)
    open(lockfile, 'wt')

    query = """
        SELECT %s
            dbo.fIAUNameVHS(ra, dec) as SourceName, *
        FROM
            vhsSource
        WHERE
            framesetid = %s """

    if args.limit:
        query = query % ('TOP 100', str(framesetid))
    else:
        query = query % ('', str(framesetid))

    print('SQL query:', query)

    # Set the third, optional argument of get to 1 if you wish to use raw mode.

    URLDBLOGIN = "http://%s:8080/vdfs/DBLogin?archive=%s&user=%s&passwd=%s&community=%s&submit=Login"

    q = URLDBLOGIN % (host, archive, username, password, community)
    print('Connection:', q)
    logging.info(" Timestamp logging info")
    print('Connecting now....')

    response = urllib2.urlopen(q)
    request = urllib2.Request(q)

    # Extract the cookies from the response header
    # and use them for future connections
    cj = cookielib.CookieJar()
    cj.extract_cookies(response, request)
    opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))

    args.format = 'fits.gz'
    args.format = 'fits'
    # This is the base query.
    if args.format == 'fits' or args.format == 'fits.gz':
        dd = {'formaction': 'freeform', 'sqlstmt': query,
              'emailAddress': args.email,
              'database': args.db,
              'format': 'FITS', 'compress': 'GZIP', 'rows': 30}
    elif args.format == 'fits.gz':
        dd = {'formaction': 'freeform', 'sqlstmt': query,
              'emailAddress': args.email,
              'database': args.db,
              'format': 'FITS', 'compress': 'GZIP', 'rows': 30}
    elif args.format == 'csv':
        dd = {'formaction': 'freeform', 'sqlstmt': query,
              'emailAddress': args.email,
              'database': args.db,
              'format': 'CSV', 'compress': 'NONE', 'rows': 30}
    elif args.format == 'html':
        dd = {'formaction': 'freeform', 'sqlstmt': query,
              'emailAddress': args.email,
              'database': args.db,
              'format': 'HTML', 'compress': 'NONE', 'rows': 30}

    urlquery = "http://%s:8080/vdfs/WSASQL?archive=VSA&" % (host)

    basequery = urlquery + urllib.urlencode(dd)

    print('basequery:', basequery)

    t0 = time.time()
    res = opener.open(basequery).read()
    print('VSA query completed; Elapsed time(secs):', time.time() - t0)

    htmlresult = outpath + "wsa_freeform_result_fs" + str(framesetid) + ".html"

    f = open(htmlresult, 'w')
    f.write(res)
    f.close()
    print("%s file saved" % htmlresult)

    # Look for the csv/fits link in web page
    datalink = re.compile(
        "<a href=\"(\S+%s).+" % args.format).search(res).group(1)
    # Request the file
    print('datalink:', datalink)

    result = opener.open(datalink).read()
    print('type(result), len(result):', type(result), len(result))
    print('Elapsed time(secs):', time.time() - t0)

    # Save file to local disk
    if args.output:
        localcsv = args.output
    else:
        localcsv = "/tmp/" + datalink.split('/')[-1]

    f = open(outfile, 'w')
    f.write(result)
    f.close()

    print("%s file saved" % outfile)
    print('Remove lockfile:', lockfile)
    if os.path.exists(lockfile):
        os.remove(lockfile)

    print('Elapsed time(secs):', time.time() - t0)

    return

# main work loop
i = 0
for fsid in fsids:
    i = i + 1
    print('fsid:', fsid, i, 'out of:', len(fsids))
    print('Number of sources in frameset:', nsources[i - 1])

    result = doquery(fsid, debug=debug)

    print('Completed frameset:', fsid)
    print('Elapsed time(secs):', time.time() - t0)


print('Elapsed time(secs):', time.time() - t0)
