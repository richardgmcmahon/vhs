""" make Python 2.7 and 3.x compatible """
from __future__ import print_function, division

def set_mjdrange_vista(period=86,
                       verbose=True, debug=False, test=False):
  """ returns mjd range for ESO periods for VISTA public surveys

  based on set_mjdrange_vista.pro which is based on set_mjdrange_ukidss

  could be generalised to support ESO generically with a Dry Run as
  an exception

  To set the date from the start of the P84 to the end of P86, use
  two calls as follows.

  mjdrange0=set_mjdrange_vista(period='p84')
  mjdrange1=set_mjdrange_vista(period='p86')
  mjdrange=[mjdrange0[0],mjdrange1[1]]


  # data/run release date limits
  # year, month, day
  # Period 87:  Apr 2011 - Sep 2011

  # note that VISTA p85 Dry Run ended on Feb  28th


  """

  import os
  import sys

  from datetime import date
  from datetime import datetime
  import calendar
  from dateutil.relativedelta import relativedelta


  import numpy as np

  from astropy.time import Time

  #print('__name__: ',__name__)
  print('Period:', period)

  if period is None:
      period=86

  # use double precision to give mjd accuracy at the 1 second level
  mjdrange_all=np.array([50000.0,100000.0], dtype='d')

  print()
  print('ISO: ', \
     Time(mjdrange_all[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_all[1], scale='utc', format='mjd').isot)
  print()


  mjdrange_vhsdr1 = np.array([0,0], dtype='d')
  mjdrange_vhsdr2 = np.array([0,0], dtype='d')
  mjdrange_vhsdr3 = np.array([0,0], dtype='d')
  mjdrange_vhsdr4 = np.array([0,0], dtype='d')
  mjdrange_vhsdr5 = np.array([0,0], dtype='d')


  # VISTA periods
  p84_VISTA_date_range = [[2009, 10, 1], [2010, 2, 28]]
  dryrun_date_range = p84_VISTA_date_range
  # not that VISTA p85 started on March 1st 2010 and not April 1st 2010
  # so we call it p85v here
  p85vISTA_date_range = [[2010,3,1],[2010,9,30]]

  # Normal ESO periods
  p84_date_range=[[2009,10,1],[2010,3,31]]
  p85_date_range=[[2010,4,1],[2010,9,30]]
  p86_date_range=[[2010,10,1],[2011,3,31]]
  p87_date_range=[[2011,4,1],[2011,9,30]]
  p88_date_range=[[2011,10,1],[2012,3,31]]

  #
  year=1967
  month=9
  end_of_month=calendar.monthrange(year, month)[1]
  print('year, month, end_of_month:', year, month, end_of_month)

  date_plus_6months = date(year, month, end_of_month) + relativedelta(months=+6)
  end_of_month = calendar.monthrange(date_plus_6months.year,
      date_plus_6months.month)[1]
  print('End of Observing Period:',
        date_plus_6months.year, date_plus_6months.month, end_of_month)

  #datetime.date(2011, 1, 31)

  # Here we try to generate any period starting from ESO Period 1
  # ESO Observing Period 1
  year1 = 1968
  month1 = 10

  year0 = year1 - 1

  iperiod = 86
  period_year_start = year1 + int((iperiod - 1) / 2)
  period_month_start = (date(year1, month1, 1) +
      relativedelta(months =+ (iperiod * 6))).month

  print(iperiod, period_year_start, period_month_start, 1)

  if test:
    for iperiod in range(1, 101, 1):
      period_year_start = year0 + int(iperiod / 2)
      period_month_start = (date(year0 + 1, month1, 1) +
          relativedelta(months =+ (iperiod * 6))).month

      print()
      print('Period:', iperiod)
      print('Start of Observing Period:',
            iperiod, period_year_start, period_month_start, 1)

      year = period_year_start
      month = period_month_start - 1
      end_of_month=calendar.monthrange(year, month)[1]
      date_plus_6months = date(year, month, end_of_month) + relativedelta(months=+6)
      end_of_month = calendar.monthrange(date_plus_6months.year,
                                         date_plus_6months.month)[1]
      print('End of Observing Period:',
        date_plus_6months.year, date_plus_6months.month, end_of_month)

    raw_input("Enter any key to continue: ")

  # data/run release date limits
  # year, month, day
  #Period 87:  Apr 2011 - Sep 2011

  # note that VISTA p85 Dry Run ended on Feb  28th



  # VISTA periods
  p84v_date_range=[[2009,10,1],[2010,2,28]]
  dryrun_date_range=p84v_date_range
  # not that VISTA p85 started on March 1st
  p85v_date_range=[[2010,3,1],[2010,9,30]]

  # Normal ESO periods
  p84_date_range=[[2009,10,1],[2010,3,31]]
  p85_date_range=[[2010,4,1],[2010,9,30]]

  p86_date_range=[[2010,10,1],[2011,3,31]]

  p94_date_range=[[2014,10,1],[2015,3,31]]

  p96_date_range=[[2015,10,1],[2016,3,31]]

  p98_date_range=[[2016,10,1],[2017,3,31]]

  end_of_month=calendar.monthrange(year1, month1)[1]
  print(year1, month1, end_of_month)

  date_plus_6months=date(year1,month1,end_of_month)+relativedelta(months=+6)
  end_of_month=calendar.monthrange(date_plus_6months.year,
   date_plus_6months.month)[1]
  print(end_of_month)
  print(date_plus_6months.year,date_plus_6months.month,
   end_of_month)
  #datetime.date(2011, 1, 31)

  iperiod_start = 86
  year_start = year1 + int((iperiod_start-1)/2)
  month_start = (date(year1,month1,1) +
   relativedelta(months=+(iperiod*6))).month
  print(iperiod, year_start, month_start, 1)

  end_of_month=calendar.monthrange(year_start, month_start)[1]

  date_plus_5months = date(year_start,month_start,end_of_month) + \
      relativedelta(months=+5)
  end_of_month=calendar.monthrange(date_plus_5months.year,
                                   date_plus_5months.month)[1]
  print(date_plus_5months.year,date_plus_5months.month,
        end_of_month)

  print('p86_date_range: ',   p86_date_range)

  if debug:
    raw_input("Enter any key to continue: ")

  # ESO Period 1
  start_year = 1967
  start_month = 10
  end_of_month = calendar.monthrange(start_year, start_month)[1]

  even_period_start=(start_year, start_month, end_of_month)
  even_period_end=(start_year+1, 2, 28)

  p84_date_range=[[2009,10,1],[2010,3,31]]
 # date_range=

  print('Date range:')
  mjdrange = [-99.9, -99.9]
  for i in range(2):
      mjdrange[i] = Time(datetime( \
       p84_date_range[i][0], p84_date_range[i][1], p84_date_range[i][2]), \
       scale='utc').mjd

  print('MJD: ' + str(mjdrange[0]) + '  ' + str(mjdrange[1]))
  print('ISO: ', \
  Time(mjdrange[0], scale='utc', format='mjd').isot, \
  Time(mjdrange[1], scale='utc', format='mjd').isot)




  #period_date_range =   {p84_date_range: p84_date_range}
  # trawl through the data releases
  # note order needed for JULDAY
  # JULDAY(month,day,year,hour,min,sec)
  #mjdrange_p84=make_array(2, /double)
  mjdrange_p84=np.array([0,0], dtype='d')
  if period.lower() == 'p84' or period.lower() == 'dryrun' \
   or period.lower().find('dr1') >= 0 or period.lower().find('dr3'):


    print('P84 Dry Run  date range:')
    for i in range(2):
      mjdrange_p84[i]=Time(datetime( \
       p84_date_range[i][0], p84_date_range[i][1],p84_date_range[i][2]), \
       scale='utc').mjd

    print('MJD: ' + str(mjdrange_p84[0]) + '  ' + str(mjdrange_p84[1]))
    print('ISO: ', \
     Time(mjdrange_p84[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p84[1], scale='utc', format='mjd').isot)

    if period.lower().find('dr1'): mjdrange_vhsdr3[0]=mjdrange_p84[0]
    if period.lower().find('dr3'): mjdrange_vhsdr3[0]=mjdrange_p84[0]

  mjdrange_p85=np.array([0,0], dtype='d')
  if period.lower() == 'p85' or period.lower() == 'dr1':
    print('P85 date range:')
    for i in range(2):
      mjdrange_p85[i]=Time(p85_date_range[i][1], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p85[0]) + '  ' + str(mjdrange_p85[1]))
    print('ISO: ', \
     Time(mjdrange_p85[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p85[1], scale='utc', format='mjd').isot)


  mjdrange_p86=np.array([0,0], dtype='d')
  if period.lower() == 'p86' or period.lower() == 'dr2':
    print('P86 date range:')
    for i in range(2):
      mjdrange_p86[i]=Time(p86_date_range[1,i], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p86[0]) + '  ' + str(mjdrange_p86[1]))
    print('ISO: ', \
     Time(mjdrange_p86[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p86[1], scale='utc', format='mjd').isot)

  mjdrange_p87=np.array([0,0], dtype='d')
  if period.lower() == 'p87' or period.lower() == 'dr2':
    print('P87 date range:')
    for i in range(2):
      mjdrange_p87[i]=Time(p87_date_range[1,i], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p87[0]) + '  ' + str(mjdrange_p87[1]))
    print('ISO: ', \
     Time(mjdrange_p87[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p87[1], scale='utc', format='mjd').isot)

  mjdrange_p88=np.array([0,0], dtype='d')
  if period.lower() == 'p88':
    print('P88 date range:')
    for i in range(2):
      mjdrange_p88[i]=Time(p88_date_range[1,i], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p88[0]) + '  ' + str(mjdrange_p88[1]))
    print('ISO: ', \
     Time(mjdrange_p88[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p88[1], scale='utc', format='mjd').isot)


  mjdrange_p89=np.array([0,0], dtype='d')
  if period.lower() == 'p89':
    print('P89 date range:')
    for i in range(2):
      mjdrange_p89[i]=Time(p89_date_range[1,i], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p89[0]) + '  ' + str(mjdrange_p89[1]))
    print('ISO: ', \
     Time(mjdrange_p89[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p89[1], scale='utc', format='mjd').isot)

  mjdrange_p90=np.array([0,0], dtype='d')
  if period.lower() == 'p90':
    print('P90 date range:')
    for i in range(2):
      mjdrange_p90[i]=Time(p90_date_range[1,i], scale='utc', format='mjd').mjd

    print('MJD: ' + str(mjdrange_p90[0]) + '  ' + str(mjdrange_p90[1]))
    print('ISO: ', \
     Time(mjdrange_p90[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p90[1], scale='utc', format='mjd').isot)

  mjdrange_p91=np.array([0,0], dtype='d')
  if period.lower() == 'p91':
    print('P91 date range:')
    for i in range(2):
      mjdrange_p91[i]=Time(datetime( \
       p91_date_range[i][0], p91_date_range[i][1],p91_date_range[i][2]), \
       scale='utc').mjd

    print('MJD: ' + str(mjdrange_p91[0]) + '  ' + str(mjdrange_p91[1]))
    print('ISO: ', \
     Time(mjdrange_p91[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_p91[1], scale='utc', format='mjd').isot)


  mjdrange_dr3=np.array([0,0], dtype='d')
  if period.lower() == 'dr3':
    print('dr3 date range:')

    mjdrange_dr3[0]=\
     Time(datetime( \
     p84_date_range[1][0], p84_date_range[1][1],p84_date_range[1][2]), \
     scale='utc').mjd

    mjdrange_dr3[1]=\
     Time(datetime( \
     p91_date_range[1][0], p91_date_range[1][1],p91_date_range[1][2]), \
     scale='utc').mjd

    print('MJD: ' + str(mjdrange_dr3[0]) + '  ' + str(mjdrange_dr3[1]))
    print('ISO: ', \
     Time(mjdrange_dr3[0], scale='utc', format='mjd').isot, \
     Time(mjdrange_dr3[1], scale='utc', format='mjd').isot)

  if period != None:

    if period.lower() == 'p84': mjdrange=mjdrange_p84
    if period.lower() == 'p85': mjdrange=mjdrange_p85
    if period.lower() == 'p86': mjdrange=mjdrange_p86
    if period.lower() == 'p87': mjdrange=mjdrange_p87
    if period.lower() == 'p88': mjdrange=mjdrange_p88
    if period.lower() == 'p89': mjdrange=mjdrange_p89
    if period.lower() == 'p90': mjdrange=mjdrange_p90
    if period.lower() == 'p91': mjdrange=mjdrange_p91
    if period.lower() == 'p92': mjdrange=mjdrange_p92
    if period.lower() == 'p93': mjdrange=mjdrange_p93
    if period.lower() == 'p94': mjdrange=mjdrange_p94
    if period.lower() == 'p95': mjdrange=mjdrange_p95
    if period.lower() == 'p96': mjdrange=mjdrange_p96
    if period.lower() == 'p97': mjdrange=mjdrange_p97

    if period.lower() == 'dr1': mjdrange=mjdrange_dr1
    if period.lower() == 'dr2': mjdrange=mjdrange_dr2
    if period.lower() == 'dr3': mjdrange=mjdrange_dr3
    if period.lower() == 'dr4': mjdrange=mjdrange_dr4

    if period.lower()  == 'vhsdr1_dr2' or period.lower() == 'vhsdr1dr2':
      mjdrange=[mjdrange_vhsdr1[0],mjdrange_vhsdr2[1]]

  return mjdrange


if __name__ == '__main__':

  """

  run unit tests

  """


  i =- 1
  for iperiod in range(84, 101, 1):

    i += 1
    print('i, Period:', i, iperiod)
    period = 'p'+str(iperiod)
    print('Period:', period)

    mjdrange = set_mjdrange_vista(period=period, verbose=True, debug=True)


  mjdrange=set_mjdrange_vista(period='p84',verbose=True, debug=True)

  mjdrange=set_mjdrange_vista(period='p91',verbose=True, debug=True)

  mjdrange=set_mjdrange_vista(period='dr3',verbose=True, debug=True)
