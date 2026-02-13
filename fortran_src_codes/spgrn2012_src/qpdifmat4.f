      subroutine qpdifmat4(ly,ldeg,rr,mat)
      implicit none
      integer ly,ldeg
      double precision rr
      double complex mat(6,6)
c
      include 'qpglobal.h'
c
      double precision alfa
      double complex cldeg,clp1,cll1
      double complex crr,crorr,clarr,cgrrr,cgarr
c
      double complex c1,c2,c4
      data c1,c2,c4/(1.d0,0.d0),(2.d0,0.d0),(4.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clp1=cldeg+c1
      cll1=cldeg*clp1
c
      crr=dcmplx(rr,0.d0)
c
      alfa=(rr-rrlw(ly))/(rrup(ly)-rrlw(ly))
      crorr=dcmplx(rolw(ly)
     &     *dexp((dlog(roup(ly))-dlog(rolw(ly)))*alfa),0.d0)
      clarr=clalw(ly)+(claup(ly)-clalw(ly))*dcmplx(alfa,0.d0)
      cgarr=cgalw(ly)+(cgaup(ly)-cgalw(ly))*dcmplx(alfa,0.d0)
      cgrrr=cgrlw(ly)+(cgrup(ly)-cgrlw(ly))*dcmplx(alfa,0.d0)
c
      mat(1,1)=crorr*cgrrr/clarr-c1/crr
      mat(1,2)=cll1/crr-crorr*comi2*crr/clarr
      mat(1,3)=-crorr*crr/clarr
      mat(1,4)=(0.d0,0.d0)
c
      mat(2,1)=c1/crr
      mat(2,2)=(0.d0,0.d0)
      mat(2,3)=(0.d0,0.d0)
      mat(2,4)=(0.d0,0.d0)
c
      mat(3,1)=cgarr/crr
      mat(3,2)=(0.d0,0.d0)
      mat(3,3)=-clp1/crr
      mat(3,4)=c1/crr
c
      mat(4,1)=cgarr*clp1/crr
      mat(4,2)=-cgarr*cll1/crr
      mat(4,3)=(0.d0,0.d0)
      mat(4,4)=cldeg/crr
      return
      end