c23456789012345678901234567890123456789012345678901234567890123456789012
c
c this include files contains general declrations and assigns values
c used in the ECA&D fortran modules
c
c GvdS Jul 2010 the choices for start years of the periods have been changed from:
c               data periods / 1901,1946,1951,1961,1976,1979/
c               to:
c               data periods / 1901,1950,1976/
c GvdS Aug 2010 year 1961 is added as startyear
c GvdS May 2011 data periods / 1901,1950,1961,1976/ is changed to:
c                          /1851, 1901, 1951, 1979
c GvdS Aug 2011 percentage and percnormal changed to 0.70 (used to be 0.80 and 0.90
c               respectively
c               values pyear=362,p6month=181,pseason=86,pmonth=25 changed to:
c               pyear=350,p6month=175,pseason=85,pmonth=25

      integer      yrbeg,yrend,calyrbeg,calyrend,maxstat
      parameter    (yrbeg=1700,yrend=2020,calyrbeg=1961,calyrend=1990)
      parameter    (maxstat=5000)
      real*8       pi
      parameter    (pi = 3.1415926535897932384626433832795d0)
      integer      absent
      parameter    (absent=-9999)
      integer      nstat
      real*8       statsN(maxstat,5)
      character*40 statsC(maxstat,2)
c
c this part sets parameters etc.
c
c
c number of days which must be present before a year/6month/season/month index can be calculated
c note that the pmonth value is KNMI-based: it relates to 80% of the data
      integer      pyear,p6month,pseason,pmonth
      parameter    (pyear=350,p6month=175,pseason=85,pmonth=25)
c
c percentage of data which must be present for trend and homogeneity calculation
      real*8       percentage
      parameter    (percentage=0.70D0)
c percentage of data which must be present for climatology calculation
      real*8       percnormal
      parameter    (percnormal=0.70D0)
c
c number of "seasons" used: 0=annual, 1=ONDJFM, 2=AMJJAS, 3=DJF, 4=MAM, 5=JJA, 6=SON
c 7=Jan, 8=Feb, 9=Mar, 10=Apr, 11=May, 12=Jun, 13=Jul, 14=Aug, 15=Sep, 16=Oct, 17=Noc, 18=Dec.
      integer      nseason
      parameter    (nseason=19)
c
c minrunmean: number of datapoints which must be present before a running mean is calculated (relates to the whole
c length of the array: used for the special time filter?
c maxstep: maximum number of datapoints in a (simple) weighted running mean
c npoint: length of window in a weighted running mean
      integer      minrunmean, maxstep,npoint
      parameter    (minrunmean=25,maxstep=101,npoint=21)
c
c nwl: in the calculation of the calender mean a nwl length window is used centered on each calender day
c note: nwl needs to be odd
      integer      nwl
      parameter    (nwl=5)
c
c fancy: .true. : distributions are determined by fits to a three- or two-parameter Gamma distribution
c        .false. : distributions are determined by simple sorting of the data
      logical      fancy
      parameter    (fancy=.false.)
c
c when calculating percentiles of precipitation of wet days (R > 1mm): one might see dates with no rain (or < 1 mm)
c in the calibration period. Similar problems occur for temperature
c if # days < pdayspresent then percentiles are set to "absent"
      integer      pdayspresent
      parameter    (pdayspresent=12)
c 
c the significance testing of the trends is done with a simple t-test (null-hypothesis: trend is indistinguishable
c from zero). In this test, one needs to estimate the number of degrees of freedom. In uncorrelated data, this 
c number is equal to the number of samples. When the data are correlated, this number is lower. See \S 6.6.8 of
c Von Storch and Zwiers for theory. 
c if vonstorch = .false., the data are asumed to be independent, 
c if vonstorch = .true., the autocorrelation between data values is calculated and an estimate for the number 
c                        of degrees of freedom is made.
      logical      vonstorch
      parameter    (vonstorch=.true.)
c
c calculating exceedences from a percentile-based threshold may introduce inhomogeneities in the index.
c See Zhang et al. (2005) and rep070827 for theory.
c if zhang = .false., the `ordinary' way of calculating exceedences is invoked.
c if zhang = .true., the Zhang et al. approach is used
      logical      zhang
      parameter    (zhang=.true.)
c
c define an array which links ind_abbr to ind_id
      integer      mxindices
      parameter    (mxindices=83)
      character*10 indexnames(mxindices)
      integer      indexids(mxindices)
      data indexnames / 'TG','TN','TX','DTR','ETR','GD4','GSL',
     +                  'vDTR','CFD','FD','HD17',
     +                  'ID','CSDI','TG10p','TN10p','TX10p',
     +                  'SU','TR','WSDI',
     +                  'TG90p','TN90p','TX90p','RR','RR1','SDII',
     +                  'CDD','CWD','R10mm',
     +                  'R20mm','RX1day','RX5day','R75p',
     +                  'R75pTOT','R95p','R95pTOT','R99p',
     +                  'R99pTOT','PP','SPI3','SPI6','SD',
     +                  'SS','TXx','TNx','TXn','TNn','SSp','PET','SD1',
     +                  'SD5cm','SD50cm',
     +                  'CD','CW','WD','WW','WGF',
     +                  'FXx','FG6Bft','FGcalm',
     +                  'DDnorth','DDsouth','DDwest','DDeast',
     +                  'FG','vPP','CSU','RH','CC','CC2','CC6',
     +                  'PRCPTOT','TCI','SSp0','SSp1','SSp2',
     +                  'SSp6','SSp7','SSp8','UTCI','HI','BEDD','TCI60',
     +                  'TCI80'/
      data indexids / 1,2,3,4,5,6,7,
     +                8,9,10,11,
     +                12,13,15,16,17,
     +                18,19,20,
     +                22,23,24,25,26,27,
     +                28,29,30,
     +                31,32,33,34,
     +                35,36,37,38,
     +                39,40,49,50,52,
     +                53,54,55,56,57,58,59,60,
     +                61,62,
     +                64,65,66,67,68,
     +                70,73,75,
     +                76,77,78,79,
     +                80,81,82,83,84,85,86,87,88,89,
     +                90,91,92,93,94,95,96,97,98,99/

c
c define an array which sets the intervals over which the trends and the homogeneity are calculated.
      integer      mxperiods
      parameter    (mxperiods=4)
      integer      periods(mxperiods)
      data periods / 1851,1901,1951,1979/
