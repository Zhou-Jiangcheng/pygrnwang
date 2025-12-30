      subroutine qppsvkerng(f,ldeg,ypsv)
      use qpalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ypsv(6,4): psv solution vector (complex) with gravity effect
c
      integer*4 ldeg
      real*8 f
      complex*16 ypsv(6,4)
c
      integer*4 i,istp,lyup,lylw
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      lyup=lyupatm(ldeg)
      lylw=lylwpsv(ldeg)
      if(lylw.lt.max0(lys,lyr))return
      comi=dcmplx(PI2*f,PI2*fi)
      comi2=comi*comi
c
      if(ldeg.eq.0)then
        call qpspropg0(ypsv,lyup,lylw)
      else if(ldeg.eq.1)then
        call qpspropg1(ypsv,lyup,lylw)
      else
        call qpspropg(ypsv,ldeg,lyup,lylw)
      endif
      return
      end