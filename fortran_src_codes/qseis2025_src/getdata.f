      subroutine getdata(unit,line)
      implicit none
      integer*4 unit
      character line*180,char*1
c
      integer*4 i
c
c     this subroutine reads over all comment lines starting with "#".
c
      char='#'
100   continue
      if(char.eq.'#')then
        read(unit,'(a)')line
        i=1
        char=line(1:1)
200     continue
        if(char.eq.' ')then
          i=i+1
          char=line(i:i)
          goto 200
        endif
        goto 100
      endif
c
      return
      end
