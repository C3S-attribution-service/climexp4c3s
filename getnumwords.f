        integer function getnumwords(string)
*
*       count the number of space-delimited words in string
*
        implicit none
        character string*(*)
        integer i,n
*
*       hey, this is Fortran!
*
        i = 0
        n = 0
  100   continue
        if ( i.ge.len(string) ) goto 800
        i = i + 1
		! make sure tabs and windows-style newlines (\r\n) are not accidentally counted as a separate word
        if ( string(i:i).eq.' ' .or. string(i:i).eq.char(9)
     +		.or. string(i:i).eq.char(13) ) goto 100
*       found the beginning of a word
        n = n + 1
  200   continue
        if ( i.ge.len(string) ) goto 800
        i = i + 1
        if ( string(i:i).ne.' ' .and. string(i:i).ne.char(9)
     +		.and. string(i:i).ne.char(13) ) goto 200
        goto 100
  800   continue
        getnumwords = n
        end
