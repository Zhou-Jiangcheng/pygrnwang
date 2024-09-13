      subroutine qpdifmat2(ly,ldeg,rr,mat)
      use qpalloc
      implicit none
      integer*4 ly,ldeg
      real*8 rr
      complex*16 mat(6,6)
c
      real*8 mass,ro1,dro
      complex*16 cup,clw,crr,drr,crorr,cxirr,clarr,cmurr,cgrrr
c
      complex*16 c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      crr=dcmplx(rr,0.d0)
      cup=(crr-rrlw(ly))/(crrup(ly)-crrlw(ly))
      clw=c1-cup
      crorr=cup*croup(ly)+clw*crolw(ly)
      clarr=cup*claup(ly)+clw*clalw(ly)
      cmurr=cup*cmuup(ly)+clw*cmulw(ly)
      cxirr=clarr+c2*cmurr
c
      dro=(roup(ly)-rolw(ly))/(rrup(ly)-rrlw(ly))
      ro1=rolw(ly)-dro*rrlw(ly) 
      mass=PI*(rr-rrlw(ly))
     &    *((4.d0/3.d0)*ro1*(rr**2+rr*rrlw(ly)+rrlw(ly)**2)
     &    +dro*(rr**3+rr**2*rrlw(ly)+rr*rrlw(ly)**2+rrlw(ly)**3))
      cgrrr=cgrlw(ly)*(crrlw(ly)/crr)**2
     &     +dcmplx(BIGG*mass/rr**2,0.d0)
c
      mat(1,1)=(c1-c2*clarr/cxirr)/crr
      mat(1,2)=c1/cxirr/crr
c
      mat(2,1)=c4*(cmurr*(c1+c2*clarr/cxirr)/crr-crorr*cgrrr)
     &        -crorr*comi2*crr
      mat(2,2)=c2*clarr/cxirr/crr
c
      return
      end