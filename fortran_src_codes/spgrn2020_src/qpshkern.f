      subroutine qpshkern(f,ldeg,ysh)
      use qpalloc
      implicit none
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ysh(2,2): sh solution vector (complex)
c
      integer*4 ldeg
      real*8 f
      complex*16 ysh(2,2)
c
      integer*4 i,j,ly,istp,lylw,lwup
      real*8 xmin
      complex*16 cldeg,cll1,cruplw,delta
      complex*16 spbphj,spbphy
c
      complex*16 c1,c2
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      lylw=lylwsh(ldeg)
      if(lyr.lt.lyob.or.lylw.lt.max0(lys,lyr))return
c
      if(ldeg.eq.0.or.ldeg.le.1.and.f.le.0.d0)return
      comi=dcmplx(PI2*f,PI2*fi)
      comi2=comi*comi
      xmin=dsqrt(2.d0*dble(2*ldeg+1))
      cldeg=dcmplx(dble(ldeg),0.d0)
      cll1=cldeg*(cldeg+c1)
c
      do ly=lyob,min0(lylwsh(ldeg),lycm-1)
        if(ldeg.le.ldegsh(ly))then
          kt(ly)=comi/cvs(ly)
          if(cdabs(kt(ly))*rrup(ly).le.xmin)then
            ksmallsh(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,3)=cldeg*cruplw
            wh(ldeg,ly,3)=-(cldeg+c1)*cruplw
          else if(ksmallsh(ly))then
            ksmallsh(ly)=.false.
            call spbjh(ldegsh(ly),kt(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,3),zjlw(0,ly,3),wj(0,ly,3),
     &                 zhup(0,ly,3),zhlw(0,ly,3),wh(0,ly,3))
          endif
        endif
      enddo
c
      do ly=lyob,min0(lylwsh(ldeg),lycm-1)
        if(ly.lt.lys)then
          lwup=0
          cpt(1,ly)=-wj(ldeg,ly,3)
          cpt(2,ly)=-wh(ldeg,ly,3)
        else
          lwup=1
          cpt(1,ly)=wj(ldeg,ly,3)
          cpt(2,ly)=wh(ldeg,ly,3)
        endif
        call qptmat(ldeg,ly,lwup)
      enddo
      if(ldeg.eq.1)then
        call qptprop1(ysh,lylw)
      else
        call qptprop(ysh,lylw)
      endif
      return
      end
