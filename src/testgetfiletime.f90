program testgetfiletime
    character *255 file
    integer time
    integer getfiletime
    print *,'file?'
    read '(a)',file
    time = getfiletime(file)
    print *,'time = ',time
end

