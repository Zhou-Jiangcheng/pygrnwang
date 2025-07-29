      subroutine brunewvlet(tau,df,mm2,wvf)
      implicit none
      integer*4 mm2
      real*8 tau,df,fi
      complex*16 wvf(mm2)
c
      integer*4 l
      real*8 f
c
      real*8 pi2
      complex*16 c1
      data pi2/6.28318530717959d0/
      data c1/(1.d0,0.d0)/
c
c     for wavelet: Brune's source time function
c
      do l=1,mm2
        f=df*dble(l-1)
        wvf(l)=c1/dcmplx(1.d0,pi2*f*tau)**2
      enddo
      return
      end