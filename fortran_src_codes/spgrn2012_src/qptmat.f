      subroutine qptmat(ldeg,ly,lwup)
      use qpalloc
      implicit none
c
      include 'qpglobal.h'
c
      integer ldeg,ly,lwup
c
      double complex cldeg,clm1,cdet,cx2,c2lp1,c2lp3
      double complex ph0jup,ph1jup,ph0yup,ph1yup
      double complex ph0jlw,ph1jlw,ph0ylw,ph1ylw
      double complex spb(2,2)
      double complex spbphj,spbphy
c
      integer i,j
      double complex ci,c1,c2,c3
      data ci,c1,c2,c3/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      clm1=dcmplx(dble(ldeg-1),0.d0)
      c2lp1=c2*cldeg+c1
      c2lp3=c2*cldeg+c3
c
      if(ksmallsh(ly))then
        ph0jup=spbphj(ldeg,kt(ly)*crrup(ly))
        ph1jup=spbphj(ldeg+1,kt(ly)*crrup(ly))
c
        ph0yup=spbphy(ldeg,kt(ly)*crrup(ly))
        ph1yup=spbphy(ldeg+1,kt(ly)*crrup(ly))
c
        cx2=(kt(ly)*crrup(ly))**2
        mat2x2up(1,1,ly)=crrup(ly)*ph0jup
        mat2x2up(2,1,ly)=crrup(ly)*cmu(ly)*(clm1*ph0jup
     &                  -cx2*ph1jup/c2lp3)
        mat2x2up(1,2,ly)=crrup(ly)*ph0yup
        mat2x2up(2,2,ly)=crrup(ly)*cmu(ly)*(clm1*ph0yup
     &                  -ph1yup*c2lp1)
      else
        spb(1,1)=c1-zjup(ldeg,ly,3)
        spb(2,1)=-ci*zjup(ldeg,ly,3)
        spb(1,2)=c1
        spb(2,2)=zhup(ldeg,ly,3)
c
        do j=1,2
          mat2x2up(1,j,ly)=crrup(ly)*spb(1,j)
          mat2x2up(2,j,ly)=crrup(ly)*cmu(ly)*(clm1*spb(1,j)-spb(2,j))
        enddo
      endif
c
      if(ly.eq.lylwsh(ldeg))return
c
      if(ksmallsh(ly))then
        ph0jlw=spbphj(ldeg,kt(ly)*crrlw(ly))
        ph1jlw=spbphj(ldeg+1,kt(ly)*crrlw(ly))
c
        ph0ylw=spbphy(ldeg,kt(ly)*crrlw(ly))
        ph1ylw=spbphy(ldeg+1,kt(ly)*crrlw(ly))
c
        cx2=(kt(ly)*crrlw(ly))**2
        mat2x2lw(1,1,ly)=crrlw(ly)*ph0jlw
        mat2x2lw(2,1,ly)=crrlw(ly)*cmu(ly)*(clm1*ph0jlw
     &                  -cx2*ph1jlw/c2lp3)
        mat2x2lw(1,2,ly)=crrlw(ly)*ph0ylw
        mat2x2lw(2,2,ly)=crrlw(ly)*cmu(ly)*(clm1*ph0ylw
     &                  -ph1ylw*c2lp1)
      else
        spb(1,1)=c1-zjlw(ldeg,ly,3)
        spb(2,1)=-ci*zjlw(ldeg,ly,3)
        spb(1,2)=c1
        spb(2,2)=zhlw(ldeg,ly,3)
c
        do j=1,2
          mat2x2lw(1,j,ly)=crrlw(ly)*spb(1,j)
          mat2x2lw(2,j,ly)=crrlw(ly)*cmu(ly)*(clm1*spb(1,j)-spb(2,j))
        enddo
      endif
c
      if(lwup.eq.0)then
c
c       calculate inverse matrix at upper radius
c
        cdet=mat2x2up(1,1,ly)*mat2x2up(2,2,ly)
     &      -mat2x2up(1,2,ly)*mat2x2up(2,1,ly)
        mat2x2inv(1,1,ly)=mat2x2up(2,2,ly)/cdet
        mat2x2inv(1,2,ly)=-mat2x2up(1,2,ly)/cdet
        mat2x2inv(2,1,ly)=-mat2x2up(2,1,ly)/cdet
        mat2x2inv(2,2,ly)=mat2x2up(1,1,ly)/cdet
      else
c
c       calculate inverse matrix at lower radius
c
        cdet=mat2x2lw(1,1,ly)*mat2x2lw(2,2,ly)
     &      -mat2x2lw(1,2,ly)*mat2x2lw(2,1,ly)
        mat2x2inv(1,1,ly)=mat2x2lw(2,2,ly)/cdet
        mat2x2inv(1,2,ly)=-mat2x2lw(1,2,ly)/cdet
        mat2x2inv(2,1,ly)=-mat2x2lw(2,1,ly)/cdet
        mat2x2inv(2,2,ly)=mat2x2lw(1,1,ly)/cdet
      endif
c
      return
      end