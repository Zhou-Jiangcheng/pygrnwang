      subroutine qsbsj(dk,nk)
      implicit none
      integer*4 nk
      real*8 dk
c
      include 'qsglobal.h'
c
      integer*4 i,ir,ik
      real*8 k,x
      real*8 bessj0,bessj1,bessj
c
      do ik=1,nk
        k=dble(ik)*dk
        do ir=1,nr
          x=k*r(ir)
          bsj(ik,0,ir)=bessj0(x)
          bsj(ik,1,ir)=bessj1(x)
          if(x.gt.2.d0)then
            bsj(ik,2,ir)=bsj(ik,1,ir)*2.d0/x-bsj(ik,0,ir)
          else
            bsj(ik,2,ir)=bessj(2,x)
          endif
          if(x.gt.3.d0)then
            bsj(ik,3,ir)=bsj(ik,2,ir)*4.d0/x-bsj(ik,1,ir)
          else
            bsj(ik,3,ir)=bessj(3,x)
          endif
          if(x.gt.4.d0)then
            bsj(ik,4,ir)=bsj(ik,3,ir)*6.d0/x-bsj(ik,2,ir)
          else
            bsj(ik,4,ir)=bessj(4,x)
          endif
          bsj(ik,-1,ir)=-bsj(ik,1,ir)
          bsj(ik,-2,ir)=bsj(ik,2,ir)
          do i=-2,4
            bsj(ik,i,ir)=bsj(ik,i,ir)*geospr(ir)
          enddo
        enddo
      enddo
c
      return
      end