program maskseries

!   mask out part of a series based on values of another series
!   e.g. all months in the NAO when NINO3>0

    implicit none
    include 'param.inc'
    include 'getopts.inc'
    integer :: i,j,nperyear,my1,my2,n(npermax),nmetadata1
    real ::  data(npermax,yrbeg:yrend),mask(npermax,yrbeg:yrend)
    character :: file1*1024,var1*40,units1*120,lvar1*120,svar1*120,history1*50000, &
        metadata1(2,100)*2000
    character :: file2*1024,var2*40,units2*120,lvar2*120,svar2*120,history2*50000, &
        metadata2(2,100)*2000        
    character :: line*255
    integer,external :: rindex

    if ( command_argument_count() < 2 ) then
        print *,'usage: maskseries series1 series2 options'
        call exit(-1)
    endif
            
    lwrite = .false. 
    call get_command_argument(1,file1)
    call readseriesmeta(file1,data,npermax,yrbeg,yrend,nperyear,var1,units1,lvar1,svar1, &
        history1,metadata1,lstandardunits,lwrite)
    call get_command_argument(2,file2)
    call readseriesmeta(file2,mask,npermax,yrbeg,yrend,i,var2,units2,lvar2,svar2, &
        history2,metadata2,lstandardunits,lwrite)
    if ( i /= nperyear ) then
        write(0,*) 'maskseries: error: cannot inperpolate in time (yet)'
        write(*,*) 'maskseries: error: cannot inperpolate in time (yet)'
        call exit(-1)
    endif
    call getopts(3,command_argument_count(),nperyear,yrbeg,yrend,.false.,0,0)
    i = rindex(file2,'/')
    j = index(file2,'.dat')
    if ( j == 0 ) then
        j = index(file2,'.nc')
        if ( j == 0 ) then
            j = len_trim(file2) + 1
        endif
    endif
    if ( maxindx < 1e33 ) then
        if ( minindx > -1e33 ) then
            print '(a,g10.2,3a,g10.3)','# only when ',minindx &
                ,' &llt; ',trim(file2(i+1:j-1)),' &lt; ',maxindx
        else
            print '(3a,g10.3)','# for ',trim(file2(i+1:j-1)),' &lt; ' &
                ,maxindx
        endif
    else
        if ( minindx > -1e33 ) then
            print '(3a,g10.3)','# for ',trim(file2(i+1:j-1)),' &gt; ',minindx
        else
            write(*,*) 'maskseries: error: expected ''lt val'' or ''gt val'''
            write(0,*) 'maskseries: error: expected ''lt val'' or ''gt val'''
            call exit(-1)
        endif
    endif
    if ( mdiff > 0 ) then
        if ( lwrite ) print '(a)','# Taking monthly anomalies'
        call mdiffit(mask,npermax,nperyear,yrbeg,yrend,mdiff)
    endif
    if ( lsum > 1 ) then
        if ( lwrite ) print '(a,i3)','# Summing ',lsum
        call sumit(mask,npermax,nperyear,yrbeg,yrend,lsum,oper)
    endif
    if ( ldetrend ) then
        if ( lwrite ) print *,'Detrending'
        call detrend(mask,npermax,nperyear,yrbeg,yrend,yr1,yr2,m1,m2,lsel)
    endif
    if ( ndiff /= 0 ) then
        if ( lwrite ) print *,'Taking differences - data'
        call diffit(mask,npermax,nperyear,yrbeg,yrend,ndiff)
    endif
    if ( anom .or. (lsel > 1 .or. nfittime > 0 ) .and. ndiff <= 0 ) then
        if ( lwrite ) print *,'Taking anomalies'
        call anomal(mask,npermax,nperyear,yrbeg,yrend,yr1,yr2)
    endif

    do i=yr1,yr2
        do j=1,nperyear
            if ( mask(j,i) < maxindx .neqv. mask(j,i) > minindx ) then
                data(j,i) = 3e33
            endif
        enddo
    enddo
    call printvar(6,var1,units1,lvar1)
    call add_varnames_metadata(var2,lvar2,svar2,metadata2,'mask')
    call merge_metadata(metadata1,nmetadata1,metadata2,' ',history2,'mask_')
    call printmetadata(6,file1,' ',' ',history1,metadata1)
    call printdatfile(6,data,npermax,nperyear,yrbeg,yrend)
end program maskseries
