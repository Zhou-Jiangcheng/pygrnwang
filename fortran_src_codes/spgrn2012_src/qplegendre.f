      subroutine qplegendre(ldegcut)
      use qpalloc
      implicit none
      include 'qpglobal.h'
      integer ldegcut
c
c     calculate associated Plm(l,m,x)/(1-x^2)^(m/2)
c     where Plm are Legendre polynomials
c
      integer l,m,ir
      double precision x
c
      do ir=1,nr
        x=dcos(disrad(ir))
        plm(0,0,ir)=1.d0
        plm(0,1,ir)=0.d0
        plm(1,1,ir)=1.d0
        plm(0,2,ir)=0.d0
        plm(1,2,ir)=0.d0
        plm(2,2,ir)=3.d0
        do m=0,2
          plm(m+1,m,ir)=dble(2*m+1)*x*plm(m,m,ir)
          do l=m+2,ldegcut
            plm(l,m,ir)=(dble(2*l-1)*x*plm(l-1,m,ir)
     &               -dble(l+m-1)*plm(l-2,m,ir))/dble(l-m)
          enddo
        enddo
      enddo
      return
      end
