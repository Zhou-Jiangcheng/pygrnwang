      subroutine qpsprop0(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of response to psv source (l=0)
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,istp,ly,key
      complex*16 y0(2),c(2),yup(2),ylw(2)
      complex*16 wave(2),coef(2,2),b(2,2)
c
c===============================================================================
c
c     propagation from surface to source
c
      yup(1)=(1.d0,0.d0)
      yup(2)=(0.d0,0.d0)
c
      if(lyr.eq.lyup)call cmemcpy(yup,y0,2)
c
      do ly=lyup,lys-1
        call caxcb(mas2x2inv(1,1,ly),yup,2,2,1,c)
        wave(1)=cdexp( cps(1,ly))
        wave(2)=cdexp(-cps(2,ly))
        if(ly.ge.lyr)then
          y0(1)=y0(1)*wave(2)
          y0(2)=y0(2)*wave(2)
        endif
c
	  c(1)=c(1)*wave(1)*wave(2)
c       c(2)=c(2)
c
        call caxcb(mas2x2lw(1,1,ly),c,2,2,1,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,2)
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c
c===============================================================================
c
c     propagation from bottom to source
c
      ylw(1)=mas2x2up(1,1,lylw)
      ylw(2)=mas2x2up(2,1,lylw)
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,2)
c
      do ly=lylw-1,lys,-1
        call caxcb(mas2x2inv(1,1,ly),ylw,2,2,1,c)
        wave(1)=cdexp(-cps(1,ly))
        wave(2)=cdexp( cps(2,ly))
        if(ly.lt.lyr)then
          y0(1)=y0(1)*wave(1)
          y0(2)=y0(2)*wave(1)
        endif
c
c       c(1)=c(1)
        c(2)=c(2)*wave(1)*wave(2)
c
        call caxcb(mas2x2up(1,1,ly),c,2,2,1,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,2)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c
      y0(1)=y0(1)/crrup(lyr)
      y0(2)=y0(2)/crrup(lyr)**2
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qpsprop0: anormal exit from cdgemp!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,2
            ypsv(i,istp)=b(2,istp)*y0(i)
          enddo
        enddo
      endif
      return
      end
