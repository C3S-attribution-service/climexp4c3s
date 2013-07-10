        subroutine val2frac(pdf,npdf,val,frac)
        implicit none
        integer npdf
        real pdf(npdf),val,frac
        integer i
*
*       computes the fraction value of val in the PDF described by the
*       sorted array pdf(1:npdf).  For values outside the array the end
*       points are returned.
*
        call checkpdf(pdf,npdf)
        if ( val.gt.1e33 ) THEN
            frac = 3e33
        elseif ( val.le.pdf(1) ) then
            frac = 1/real(npdf+1)
        elseif ( val.ge.pdf(npdf) ) then
            frac = npdf/real(npdf+1)
        else
            do i=1,npdf-1
                if ( val.lt.pdf(i+1) ) then
                    frac = ((val-pdf(i))*(i+1) + (pdf(i+1)-val)*i)
     +                   /((pdf(i+1)-pdf(i))*(npdf+1))
                    exit
                endif
            enddo
        endif
        end

        subroutine frac2val(frac,pdf,npdf,val)
        implicit none
        integer npdf
        real pdf(npdf),frac,val
        integer i
        real s
*
*       computes the value val corresponding to the fraction frac in
*       the PDF described by the sorted array pdf(1:npdf).  For values
*       outside the interval 1/(npdf+1),npdf/(npdf+1) the end points are
*       returned.
*
        if ( frac.gt.1e33 ) then
            val = 3e33
        elseif ( frac.lt.0 .or.frac.gt.1 ) then
            write(0,*) 'frac2val: invalid fraction: ',frac
            call abort
        elseif ( frac.le.1/real(npdf+1) ) then
            val = pdf(1)
        elseif ( frac.ge.npdf/real(npdf+1) ) then
            val = pdf(npdf)
        else
            s = (npdf+1)*frac
            i = int(s)
            s = s - i
            val = pdf(i+1)*s + pdf(i)*(1-s)
        endif
        end

        subroutine checkpdf(pdf,npdf)
        implicit none
        integer npdf
        real pdf(npdf)
        integer i
        if ( npdf.lt.1 ) then
            write(0,*) 'checkpdf: no PDF present ',npdf
            call abort
        endif
        if ( pdf(1).gt.1e33 .or. pdf(npdf).gt.1e33 ) then
            write(0,*) 'checkpdf: invalid PDF: '
            write(0,'(i4,g16.4)') (i,pdf(i),i=1,npdf)
            call abort
        endif
        end
