      subroutine legendre(ldegcut,raddis,dp)
      implicit none
      integer ldegcut
      double precision raddis,dp(0:ldegcut)
c
c     x = cos(raddis)
c     define: P(l,x) = Legendre polynomials
c     This subroutine calculate dp(l,x) = p(l,x) - 1
c
      integer l
      double precision x,y
c
      x=dcos(raddis)
      y=-2.d0*(dsin(0.5d0*raddis))**2
c
      dp(0)=0.d0
      dp(1)=y
      do l=2,ldegcut
        dp(l)=(dble(2*l-1)*y+dble(2*l-1)*x*dp(l-1)
     &        -dble(l-1)*dp(l-2))/dble(l)
      enddo
c
      return
      end
