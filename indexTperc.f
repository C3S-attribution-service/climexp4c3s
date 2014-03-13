c23456789012345678901234567890123456789012345678901234567890123456789012
      program indexTpercentiles
c-----------------------------------------------------------------------
c     This module calculates indices based on the temperatures
c     TX, TG and TG (not combined) which make use of thresholds 
c     related to the percentiles of the frequency distribution
c
c     which set of indices is calculated depends on the input cint
c
c     cint   : name of temperature parameter (TX, TG, TN)
c-----------------------------------------------------------------------
c
c     2011/07/04 GvdS Noted a bug in subroutine fillValue. calls
c                     to this subroutine commented out
c     2012/06/05 GvdS Bugfix in calcFreqDistr: if the minimum data availability
c                     is not exceeded , then the freq. distribution 
c                     values are set to absent
c
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       i,j,k,n,numdat,stagrp,fyear,lyear,case
      character*128 ainput,infile
      character     cint*2
      real*8        p10(calyrbeg:calyrend+1,12,31)
      real*8        p90(calyrbeg:calyrend+1,12,31)
      logical       ex,verbatim

      verbatim=.false.

      if(verbatim)
     + write(6,*) '----------> indexT percentiles <---------------'

      call inputt(5,numdat)

      if (numdat.ne.1) call error('wrong amount of data on input line')

      cint = ainput()

      if(verbatim) write(6,11) cint
11    format('Index name             : ',a8)

c change to lower case
      if(cint.eq.'TX') cint='tx'
      if(cint.eq.'TN') cint='tn'
      if(cint.eq.'TG') cint='tg'

      if((cint.eq.'tx').or.(cint.eq.'TX')) then
        case = 1
        if(verbatim) write(6,*) '--- calculated are the 
     +     indices TX10p, TX90p and WSDI'
      elseif((cint.eq.'tn').or.(cint.eq.'TN')) then
        case = 2
        if(verbatim) write(6,*) '--- calculated are the 
     +     indices TN10p, TN90p and CSDI'
      elseif((cint.eq.'tg').or.(cint.eq.'TG')) then
        case = 3
        if(verbatim) write(6,*) 
     +   '--- calculated are the indices TG10p, TG90p'
      else
        call error('input not recognized')
      endif
       

      infile = cint(1:2)//'_infile_indices.txt'
      inquire(file=infile,exist=ex)

      if (.not.ex) then
        write(6,*) 'Cannot find file ',infile
      else

        call readData(infile,fyear,lyear,stagrp,a,qc)

c set the data to the unit degrees C
        call xscale(fyear,lyear,0.1d0,a)

c calculate the frequency distribution for each calender date for data
c in the period calyrbeg - calyrend 

        call calcFreqDistr(a,qc,p10,p90)

c determine which set of indices are to be calculated
        if (case.eq.1) then

          call calcTX(fyear,lyear,p10,p90,stagrp,a,qc)

        elseif (case.eq.2) then

          call calcTN(fyear,lyear,p10,p90,stagrp,a,qc)

        else

          call calcTG(fyear,lyear,p10,p90,stagrp,a,qc)

        endif

      endif
   
      end
c
c---------------------------------------------------------------------------
c
      subroutine calcTX(fyear,lyear,p10,p90,stagrp,a,qc)
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,nseason)
      real*8        p10(calyrbeg:calyrend+1,12,31)
      real*8        p90(calyrbeg:calyrend+1,12,31)
      integer       fyear,i,j,stagrp,lyear,findindexid,indid
      character     locid*6
      character     newdata*128,cindex*10

c calculate TX10p: number of days < 10th percentile
      call calcTp10(fyear,lyear,a,qc,p10,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TX10p.txt'
      cindex='TX10p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

c calculate TX90p: number of days > 90th percentile
      call calcTp90(fyear,lyear,a,qc,p90,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TX90p.txt'
      cindex='TX90p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

c calculate WSDI: warm spell duration index 
      call calcWSDI(fyear,lyear,a,qc,p90,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_WSDI.txt'
      cindex='WSDI'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine calcTN(fyear,lyear,p10,p90,stagrp,a,qc)
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,nseason)
      real*8        p10(calyrbeg:calyrend+1,12,31)
      real*8        p90(calyrbeg:calyrend+1,12,31)
      integer       fyear,i,j,stagrp,lyear,findindexid,indid
      character     newdata*128,cindex*10

c calculate TN10p: number of days < 10th percentile
      call calcTp10(fyear,lyear,a,qc,p10,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TN10p.txt'
      cindex='TN10p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

c calculate TN90p: number of days > 90th percentile
      call calcTp90(fyear,lyear,a,qc,p90,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TN90p.txt'
      cindex='TN90p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

c calculate CSDI: cold spell duration index
      call calcCSDI(fyear,lyear,a,qc,p10,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_CSDI.txt'
      cindex='CSDI'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine calcTG(fyear,lyear,p10,p90,stagrp,a,qc)
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      real*8        c(yrbeg:yrend,nseason)
      real*8        p10(calyrbeg:calyrend+1,12,31)
      real*8        p90(calyrbeg:calyrend+1,12,31)
      integer       fyear,lyear,i,j,stagrp,findindexid,indid
      character     newdata*128,cindex*10

c calculate TG10p: number of days < 10th percentile
      call calcTp10(fyear,lyear,a,qc,p10,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TG10p.txt'
      cindex='TG10p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

c calculate TG90p: number of days > 90th percentile
      call calcTp90(fyear,lyear,a,qc,p90,b)

      call lowessfilter(fyear,lyear,b,c)
 
c write to file
      newdata = 'index_nonprec_TG90p.txt'
      cindex='TG90p'
      indid = findindexid(cindex)
      call write2file(newdata,indid,stagrp,fyear,lyear,b,c)

      return
      end
c
c---------------------------------------------------------------------------
c
