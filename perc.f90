subroutine val2frac(pdf,npdf,val,frac)
    implicit none
    integer,intent(in) :: npdf
    real,intent(in) :: pdf(npdf),val
    real,intent(out) :: frac
    integer :: i

!   computes the fraction value of val in the PDF described by the
!   sorted array pdf(1:npdf).  For values outside the array the end
!   points are returned.

    call checkpdf(pdf,npdf)
    if ( val > 1e33 ) then
        frac = 3e33
    elseif ( val <= pdf(1) ) then
        frac = 1/real(npdf+1)
    elseif ( val >= pdf(npdf) ) then
        frac = npdf/real(npdf+1)
    else
        do i=1,npdf-1
            if ( val < pdf(i+1) ) then
                frac = ((val-pdf(i))*(i+1) + (pdf(i+1)-val)*i) &
                    /((pdf(i+1)-pdf(i))*(npdf+1))
                exit
            endif
        enddo
    endif
end subroutine val2frac

subroutine frac2val(frac,pdf,npdf,val)
    implicit none
    integer :: npdf
    real :: pdf(npdf),frac,val
    integer :: i
    real :: s

!       computes the value val corresponding to the fraction frac in
!       the PDF described by the sorted array pdf(1:npdf).  For values
!       outside the interval 1/(npdf+1),npdf/(npdf+1) the end points are
!       returned.

    if ( frac > 1e33 ) then
        val = 3e33
    elseif ( frac < 0 .or. frac > 1 ) then
        write(0,*) 'frac2val: error: invalid fraction: ',frac
        call exit(-1)
    elseif ( frac <= 1/real(npdf+1) ) then
        val = pdf(1)
    elseif ( frac >= npdf/real(npdf+1) ) then
        val = pdf(npdf)
    else
        s = (npdf+1)*frac
        i = int(s)
        s = s - i
        val = pdf(i+1)*s + pdf(i)*(1-s)
    endif
end subroutine frac2val

subroutine checkpdf(pdf,npdf)
    implicit none
    integer :: npdf
    real :: pdf(npdf)
    integer :: i
    if ( npdf < 1 ) then
        write(0,*) 'checkpdf: no PDF present ',npdf
        call exit(-1)
    endif
    if ( pdf(1) > 1e33 .or. pdf(npdf) > 1e33 ) then
        write(0,*) 'checkpdf: invalid PDF: '
        write(0,'(i4,g16.4)') (i,pdf(i),i=1,npdf)
        call exit(-1)
    endif
end subroutine checkpdf
