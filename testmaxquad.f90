program testmaxquad
    real :: x(10),y(10)
  1 do i=1,10
        x(i) = i
        call random_number(y(i))
    enddo
    call maxquad(xmax,ymax,x,y,10)
    read *,i
    goto 1

end program testmaxquad
