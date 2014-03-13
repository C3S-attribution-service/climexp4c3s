subroutine climexp2extreme(data,indx,npermax,fyr,lyr,yr1,yr2,nperyear,npernew, &
&	minindx,var,units,newunits,climdex,lvar,lwrite)
!
!	glue routine to convert my time series format to the extreme format,
!	call the correct subroutine and convert the results back
!
!	Geert Jan van Oldenborgh, July 2013
!
	implicit none
	include 'extreme.h'

	integer fyr,lyr,yr1,yr2,npermax,nperyear,npernew
	real data(npermax,fyr:lyr),indx(npernew,fyr:lyr)
	real minindx
	character var*(*),lvar*(*),units*(*),newunits*(*),climdex*(*)
	logical lwrite
	real*8 a(yrbeg:yrend,12,31),b(yrbeg:yrend,nseason)
	real*8 p10(calyrbeg:calyrend+1,12,31),p90(calyrbeg:calyrend+1,12,31)
	real*8 p75(calyrbeg:calyrend+1,nseason),p95(calyrbeg:calyrend+1,nseason), & 
	&      p99(calyrbeg:calyrend+1,nseason)
	real*8 rsum(yrbeg:yrend,nseason),p(calyrbeg:calyrend+1,nseason)
	real*8 thresh
	integer qc(yrbeg:yrend,12,31)
	integer mo,yr,dy,dd,i,j,k,l,fyear,lyear
	!
	!	init
	!
	if ( yr1.lt.yrbeg ) yr1 = yrbeg
	if ( yr2.gt.yrend ) yr2 = yrend
	!
	!	actually, Gerard's routines do not use yr1,yr2 to select a time 
	!	period, but to indicate the first, last year with data. Find these.
	!
	do yr=yr1,yr2
		do mo=1,nperyear
			if ( data(mo,yr).lt.1e33 ) then
				goto 101
			end if
		end do
	end do
	write(0,'(a)') 'climate2extreme: no valid data found'
	write(*,'(a)') '# climate2extreme: no valid data found'
	stop
101	continue
	yr1 = yr
	do yr=yr2,yr1,-1
		do mo=1,nperyear
			if ( data(mo,yr).lt.1e33 ) then
				goto 102
			end if
		end do
	end do
102 continue
	yr2 = yr
	!
	a = -9999
	qc = 0
	indx = 3e33
	if ( lwrite ) then
		print *,'climexp2extreme: input'
		print *,'  climdex = ',trim(climdex)
		print *,'  minindx = ',minindx
		print *,'  nperyear,npernew = ',nperyear,npernew
		print *,'  var,units = ',trim(var),' [',trim(units),']'
	end if
	if ( npernew.eq.-1 ) then
		write(0,*) 'climexp2extreme: cannot handle Jul-Jun years yet'
	end if
	if ( nperyear.eq.360 ) then
		call setcalendar('360_day')
	else if ( nperyear.eq.365 ) then
		call setcalendar('noleap')
	else if ( nperyear.eq.366 ) then
		call setcalendar('proleptic_gregorian')
	else
		write(0,*) 'climexp2extreme: error: expecting daily data, not ',nperyear
		call abort
	end if
	!
	!	copy input
	!
	fyear = yr1
	lyear = yr2
	if ( units.eq.'K' ) then
		call makestandardseries(data,npermax,fyr,lyr,nperyear,var,units,lwrite)
	end if
	a = -9999
	do yr=yr1,yr2
		dd = 0
		do mo=1,12
			call lengthofmonth(2000,mo,l) ! our storage has all leap years for nperyear 366
			do dy=1,l
				dd = dd + 1
				if ( data(dd,yr).lt.1e33 ) then
					a(yr,mo,dy) = data(dd,yr)
				else
					qc(yr,mo,dy) = 9
				end if
			end do
		end do
	end do
	!
	!	call appropriate routine, see list at http://www.climdex.org/indices.html
	!
	call tolower(var)
	thresh = minindx ! convert from real*4 to real*8
	if ( var.eq.'tn' .or. index(var,'min').ne.0 ) then ! minimum temperature
		if ( climdex.eq.'FD' ) then
			lvar = 'Frost days (TN < 0C)'
			call calcFD(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'CFD' ) then
			lvar = 'Maximum number of consecutive frost days (TN < 0C)'
			call calcCFD(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'TR' ) then
			lvar = 'Tropical nights (TN > 20C)'
			call calcTR(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'TNx' ) then
			lvar = 'Maximum value of daily minimum temperature'
			call calcMAX(fyear,lyear,a,qc,b)
			newunits = 'Celsius'
		else if ( climdex.eq.'TNn' ) then
			lvar = 'Minimum value of daily minimum temperature'
			call calcMIN(fyear,lyear,a,qc,b)
			newunits = 'Celsius'
		else if ( climdex.eq.'TN10p' ) then
			lvar = 'Days with TN < 10th percentile of daily minimum temperature (cold nights)'
			call calcFreqDistr(a,qc,p10,p90)
			call calcTp10(fyear,lyear,a,qc,p10,b)
			newunits = 'dy'
		else if ( climdex.eq.'TN90p' ) then
			lvar = 'Days with TN > 90th percentile of daily minimum temperature (warm nights)'
			call calcFreqDistr(a,qc,p10,p90)
			call calcTp90(fyear,lyear,a,qc,p90,b)
			newunits = 'dy'
		else if ( climdex.eq.'CSDI' ) then
			lvar = 'Cold-spell duration index'
			call calcFreqDistr(a,qc,p10,p90)
			call calcCSDI(fyear,lyear,a,qc,p10,b)
			newunits = 'dy'
		else
			write(0,*) 'unknown Tmin index ',trim(climdex)
			call abort
		end if

	else if ( var.eq.'tx' .or. index(var,'max').ne.0 ) then ! maximum temperature
		if ( climdex.eq.'SU' ) then
			lvar = 'Summer days (TX > 25C)'
			call calcSU(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'CSU' ) then
			lvar = 'Maximum number of consecutive summer days (TX > 25C)'
			call calcCSU(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'ID' ) then
			lvar = 'Ice days (TX < 0C)'
			call calcID(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'TXx' ) then
			lvar = 'Maximum value of daily maximum temperature'
			call calcMAX(fyear,lyear,a,qc,b)
			newunits = 'Celsius'
		else if ( climdex.eq.'TXn' ) then
			lvar = 'Minimum value of daily maximum temperature'
			call calcMIN(fyear,lyear,a,qc,b)
			newunits = 'Celsius'
		else if ( climdex.eq.'TX10p' ) then
			lvar = 'Days with TX < 10th percentile of daily maximum temperature (cold day-times)'
			call calcFreqDistr(a,qc,p10,p90)
			call calcTp10(fyear,lyear,a,qc,p10,b)
			newunits = 'dy'
		else if ( climdex.eq.'TX90p' ) then
			lvar = 'Days with TX > 90th percentile of daily maximum temperature (warm day-times)'
			call calcFreqDistr(a,qc,p10,p90)
			call calcTp90(fyear,lyear,a,qc,p90,b)
			newunits = 'dy'
		else if ( climdex.eq.'WSDI' ) then
			lvar = 'Warm-spell duration index'
			call calcFreqDistr(a,qc,p10,p90)
			call calcCSDI(fyear,lyear,a,qc,p10,b)
			newunits = 'dy'
		else
			write(0,*) 'unknown Tmax index ',trim(climdex)
			call abort
		end if

	else if ( var(1:1).eq.'t' ) then ! let's hope it is mean temperature
		if ( climdex.eq.'GSL' ) then
			lvar = 'Growing season length'
			call calcGSL(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'GD4' ) then
			lvar = 'Growing degree days (sum of TG > 4C)'
			call calcGD4(fyear,lyear,a,qc,b)
			newunits = 'K dy'
		else if ( climdex.eq.'HD17' ) then
			lvar = 'Heating degree days (sum of 17C - TG)'
			call calcHD17(fyear,lyear,a,qc,b)
			newunits = 'K dy'
		else
			write(0,*) 'unknown Tmean index ',trim(climdex)
			call abort
		end if

	else if ( var(1:2).eq.'rr' .or. var(1:2).eq.'rh' .or. &
	&         var(1:1).eq.'p' .and. var(1:4).ne.'pres' ) then ! let's hope it is precipitation
		if ( units.ne.'mm/dy' .and. units.ne.'mm dy-1' ) then
			call makestandardseries(data,npermax,fyr,lyr,nperyear,var,units,lwrite)
		end if
		if ( climdex.eq.'RX1day' ) then
			lvar = 'Highest 1-day precipitation amount'
			call calcRXday(fyear,lyear,a,qc,b)
			newunits = 'mm'
		else if ( climdex.eq.'RX5day' ) then
			lvar = 'Highest 5-day precipitation amount'
			call calcRX5day(fyear,lyear,a,qc,b)
			newunits = 'mm'
		else if ( climdex.eq.'SDII' ) then
			lvar = 'Simple daily intensity index'
			call calcSDII(fyear,lyear,a,qc,b)
			newunits = 'mm/dy'
		else if ( climdex.eq.'SDIInn' ) then
			lvar = 'Simple daily intensity index'
			call calcSDIInn(fyear,lyear,a,qc,thresh,b)
			newunits = 'mm/dy'
		else if ( climdex.eq.'RR1' ) then ! wet days
			lvar = 'Wet days (RR >= 1mm)'
			call calcRR1(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'R10mm' ) then
			lvar = 'Heavy precipitation days (precipitation >= 10mm)'
			call calcR10mm(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'R20mm' ) then
			lvar = 'Very heavy precipitation days (precipitation >= 20mm)'
			call calcR20mm(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'Rnnmm' ) then
			lvar = 'Heavy precipitation days (precipitation >= XXmm)'
			call calcRnnmm(fyear,lyear,a,qc,thresh,b)
			newunits = 'dy'
		else if ( climdex.eq.'CDD' ) then
			lvar = 'Maximum number of consecutive dry days (RR < 1mm)'
			call calcCDD(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'CWD' ) then
			lvar = 'Maximum number of consecutive wet days (RR >= 1mm)'
			call calcCWD(fyear,lyear,a,qc,b)
			newunits = 'dy'
		else if ( climdex.eq.'R95pTOT' ) then
			lvar = 'Precipitation fraction due to very wet days (> 95th percentile)'
			call calcFreqDistrRR(a,qc,p75,p95,p99)
			call calcRXXp(fyear,lyear,a,qc,p95,rsum,b)
			call calcRXXptot(fyear,lyear,a,qc,rsum,b)
			newunits = '%'
		else if ( climdex.eq.'R99pTOT' ) then
			lvar = 'Precipitation fraction due to extremely wet days (> 99th percentile)'
			call calcFreqDistrRR(a,qc,p75,p95,p99)
			call calcRXXp(fyear,lyear,a,qc,p99,rsum,b)
			call calcRXXptot(fyear,lyear,a,qc,rsum,b)
			newunits = '%'
		else if ( climdex.eq.'RnnTOT' ) then
			lvar = 'Precipitation fraction due to days above threshold'
			p99 = thresh ! cheat to use a standard routine
			call calcRXXp(fyear,lyear,a,qc,p99,rsum,b)
			call calcRXXptot(fyear,lyear,a,qc,rsum,b)
			newunits = '%'
		else if ( climdex.eq.'SPI3' ) then
			lvar = '3-Month Standardized Precipitation Index'
			call calcSPI3(fyear,lyear,a,qc,b)
			newunits = '1'
		else if ( climdex.eq.'SPI6' ) then
			lvar = '6-Month Standardized Precipitation Index'
			call calcSPI6(fyear,lyear,a,qc,b)
			newunits = '1'
		else if ( climdex.eq.'PRCPTOT' ) then
			lvar = 'Total precipitation on wet days'
			call calcPRCPTOT(fyear,lyear,a,qc,b)
			if ( npernew.eq.1 ) then
				newunits = 'mm/yr'
			else if ( npernew.eq.2 ) then
				newunits = 'mm/6mo'
			else if ( npernew.eq.4 ) then
				newunits = 'mm/3mo'
			else if ( npernew.eq.12 ) then
				newunits = 'mm/mo'
			else
				newunits = 'mm'
			end if
		else
			write(0,*) 'unknown precip index ',trim(climdex)
			call abort
		end if		

	else
		write(0,*) 'climexp2extreme: error: cannot recognise variable ', &
		&	trim(var),' [',trim(units),']'
		call abort
	end if
	!
	!	copy output back, note that the numbers are one higher than in the documentation :-)
	!
	indx = 3e33
	if ( npernew.eq.1 ) then
		do yr=yr1,yr2
			if ( b(yr,1).gt.-9998 ) indx(1,yr) = b(yr,1)
		end do
	else if ( npernew.eq.2 ) then
		do yr=yr1,yr2
			if ( b(yr,2).gt.-9998 ) indx(1,yr) = b(yr,2)
			if ( b(yr,3).gt.-9998 ) indx(2,yr) = b(yr,3)
		end do
	else if ( npernew.eq.4 ) then
		do yr=yr1,yr2
			do i=1,4
				if ( b(yr,3+i).gt.-9998 ) indx(i,yr) = b(yr,3+i)
			end do
		end do
	else if ( npernew.eq.12 ) then
		do yr=yr1,yr2
			do i=1,12
				if ( b(yr,7+i).gt.-9998 ) indx(i,yr) = b(yr,7+i)
			end do
		end do
	else
		write(0,*) 'climexp2extreme: error: unknown value for npernew ',npernew
		call abort
	end if

end subroutine climexp2extreme	
 