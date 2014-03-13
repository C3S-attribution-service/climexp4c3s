      subroutine calcR10mm(fyear,lyear,a,qc,b)
c this routine calculates the number of days with heavy precip (RR => 10 mm)
      implicit none
      include 'comgeneral.h'

      real*8        a(yrbeg:yrend,12,31)
      integer       qc(yrbeg:yrend,12,31)
      real*8        b(yrbeg:yrend,nseason)
      integer       c(yrbeg:yrend,12,2)
      integer       fyear,lyear,i,j,k,length,absen,nsum,npre
      real*8        absenr,absentr,thresh
      parameter     (thresh=10.0d0)

      call calcRnnmm(fyear,lyear,a,qc,thresh,b)

      return
      end
c
c---------------------------------------------------------------------------
c
