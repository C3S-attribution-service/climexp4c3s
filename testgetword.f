        program testgetword
        character string*80,word*80
        integer n

        print *,'give string'
        read '(a)',string
 100    continue
        print *,'N?'
        read *,n
        call getword(n,string,',',word)
        print *,trim(word)
        goto 100
        end
