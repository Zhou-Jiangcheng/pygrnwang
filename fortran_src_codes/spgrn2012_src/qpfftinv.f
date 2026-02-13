      subroutine qpfftinv(sgrn,tgrn,tred)
      implicit none
c
      include 'qpglobal.h'
c
      double precision tred
      double complex sgrn(nfmax)
      real*4 tgrn(2*nfmax)
c
      integer lf,mf,it
      double precision t,a,b
      double complex cswap(2*nfmax)
      double precision dswap(4*nfmax)
      equivalence (cswap,dswap)
c
      do lf=1,nf
        cswap(lf)=sgrn(lf)*cdexp(dcmplx(-fi,dble(lf-1)*df)
     &                         *dcmplx(PI2*tred,0.d0))
      enddo
c
      mf=1
      do lf=2*nf,nf+2,-1
        mf=mf+1
        cswap(lf)=dconjg(cswap(mf))
      enddo
      cswap(nf+1)=(0.d0,0.d0)
c
c     convention for Fourier transform:
c     f(t)=\int F(f) exp(i2\pi f t) df
c
      call four1(dswap,2*nf,+1)
      do lf=1,2*nf
        t=dble(lf-1)*dt
        cswap(lf)=cswap(lf)*dcmplx(dexp(-PI2*fi*t)*df,0.d0)
      enddo
c
      do it=1,ntout
        t=dble(it-1)*dtout
        lf=1+int(t/dt)
        if(lf.lt.2*nf)then
          b=t/dt-dble(lf-1)
          a=1.d0-b
          tgrn(it)=sngl(a*dreal(cswap(lf))+b*dreal(cswap(lf+1)))
        else
          tgrn(it)=sngl(dreal(cswap(lf)))
        endif
      enddo
      return
      end