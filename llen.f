        integer function llen(a)
        character*(*) a
        do 10 i=len(a),1,-1
            if ( a(i:i).ne.'?' .and. a(i:i).ne.' ' .and. 
     +           a(i:i).ne.char(0) ) goto 20
   10   continue
   20   continue
        llen = max(i,1)
        end
        
