      subroutine wavelet(tau,df,mm2,wvf)
      implicit none
      integer mm2
      double precision tau,df
      double complex wvf(mm2)
c
      integer l
      double precision f,x
c
      double precision pi2,eps
      data pi2,eps/6.28318530717959d0,1.0d-04/
c
      do l=1,mm2
        f=df*dble(l-1)
        x=f*tau
        if(x.eq.0.d0)then
          wvf(l)=(1.d0,0.d0)
        else if(x.ge.1.d0-eps.and.x.le.1.d0+eps)then
          wvf(l)=dcmplx(-1.d0/x/(1+x),0.d0)
        else
          wvf(l)=dcmplx(0.d0,1.d0/(pi2*x*(1.d0+x)*(1.d0-x)))
     &             *(cdexp(dcmplx(0.d0,-pi2*x))-(1.d0,0.d0))
        endif
      enddo
      return
      end
