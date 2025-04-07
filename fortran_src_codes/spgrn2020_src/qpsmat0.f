      subroutine qpsmat0(ly,lwup)
      use qpalloc
      implicit none
c
c     calculate 2x2 spheroidal layer matrix for a solid shell in case of degree 2
c
      integer*4 ly,lwup
c
      integer*4 i,j
      complex*16 cmudr,c2mu,cxi,cx2,cr2,cdet
      complex*16 ph0jup,ph1jup,ph0yup,ph1yup
      complex*16 ph0jlw,ph1jlw,ph0ylw,ph1ylw
      complex*16 spb(2,2)
      complex*16 spbphj,spbphy
c
      complex*16 ci,c1,c2,c3
      data ci,c1,c2,c3/(0.d0,1.d0),(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
      c2mu=c2*cmu(ly)
      cxi=cla(ly)+c2mu
c
      cx2=(kp(ly)*crrup(ly))**2
      if(ksmallpsv(ly))then
        ph0jup=spbphj(0,kp(ly)*crrup(ly))
        ph1jup=spbphj(1,kp(ly)*crrup(ly))
c
        ph0yup=spbphy(0,kp(ly)*crrup(ly))
        ph1yup=spbphy(1,kp(ly)*crrup(ly))
c
        cr2=crrup(ly)**2
        mas2x2up(1,1,ly)=cr2*ph1jup/c3
        mas2x2up(2,1,ly)=cxi*cr2*ph0jup-c2*c2mu*cr2*ph1jup/c3
        mas2x2up(1,2,ly)=ph1yup
        mas2x2up(2,2,ly)=cxi*cx2*ph0yup-c2*c2mu*ph1yup
      else
        spb(1,1)=c1-zjup(0,ly,1)
        spb(2,1)=-ci*zjup(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhup(0,ly,1)
c
        do j=1,2
          mas2x2up(1,j,ly)=spb(2,j)
          mas2x2up(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
        enddo
      endif
c
      if(ly.eq.lylwpsv(0))return
c
      cx2=(kp(ly)*crrlw(ly))**2
      if(ksmallpsv(ly))then
        ph0jlw=spbphj(0,kp(ly)*crrlw(ly))
        ph1jlw=spbphj(1,kp(ly)*crrlw(ly))
c
        ph0ylw=spbphy(0,kp(ly)*crrlw(ly))
        ph1ylw=spbphy(1,kp(ly)*crrlw(ly))
c
        cr2=crrlw(ly)**2
        mas2x2lw(1,1,ly)=cr2*ph1jlw/c3
        mas2x2lw(2,1,ly)=cxi*cr2*ph0jlw-c2*c2mu*cr2*ph1jlw/c3
        mas2x2lw(1,2,ly)=ph1ylw
        mas2x2lw(2,2,ly)=cxi*cx2*ph0ylw-c2*c2mu*ph1ylw
      else
        spb(1,1)=c1-zjlw(0,ly,1)
        spb(2,1)=-ci*zjlw(0,ly,1)
        spb(1,2)=c1
        spb(2,2)=zhlw(0,ly,1)
c
        do j=1,2
          mas2x2lw(1,j,ly)=spb(2,j)
          mas2x2lw(2,j,ly)=cxi*cx2*spb(1,j)-c2*c2mu*spb(2,j)
        enddo
      endif
c
      if(lwup.eq.0)then
c
c       calculate inverse masrix at upper radius
c
        cdet=mas2x2up(1,1,ly)*mas2x2up(2,2,ly)
     &      -mas2x2up(1,2,ly)*mas2x2up(2,1,ly)
        mas2x2inv(1,1,ly)=mas2x2up(2,2,ly)/cdet
        mas2x2inv(1,2,ly)=-mas2x2up(1,2,ly)/cdet
        mas2x2inv(2,1,ly)=-mas2x2up(2,1,ly)/cdet
        mas2x2inv(2,2,ly)=mas2x2up(1,1,ly)/cdet
      else
c
c       calculate inverse masrix at lower radius
c
        cdet=mas2x2lw(1,1,ly)*mas2x2lw(2,2,ly)
     &      -mas2x2lw(1,2,ly)*mas2x2lw(2,1,ly)
        mas2x2inv(1,1,ly)=mas2x2lw(2,2,ly)/cdet
        mas2x2inv(1,2,ly)=-mas2x2lw(1,2,ly)/cdet
        mas2x2inv(2,1,ly)=-mas2x2lw(2,1,ly)/cdet
        mas2x2inv(2,2,ly)=mas2x2lw(1,1,ly)/cdet
      endif
c
      return
      end