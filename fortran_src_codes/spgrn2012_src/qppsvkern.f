      subroutine qppsvkern(f,ldeg,ypsv)
      use qpalloc
      implicit none
c
      include 'qpglobal.h'
c
c     calculation of response function in frequency-wavelength domain
c     ldeg: harmonic degree
c     ypsv(6,4): psv solution vector (complex)
c
      integer ldeg
      double precision f
      double complex ypsv(6,4)
c
      integer i,ly,istp,lyup,lylw,lwup
      double precision xmin
      double complex cldeg,cruplw,delta
c
      double complex c1,c2
      data c1,c2/(1.d0,0.d0),(2.d0,0.d0)/
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      lyup=lyupatm(ldeg)
      lylw=lylwpsv(ldeg)
      if(lylw.lt.max0(lys,lyr).or.f.le.0.d0.and.ldeg.eq.1)return
c
      if(f.le.0.d0)then
        lylw=min0(lylw,lycm)
        lyup=lyob
      endif
c
      cldeg=dcmplx(dble(ldeg),0.d0)
      comi=dcmplx(PI2*f,PI2*fi)
      comi2=comi*comi
      xmin=dsqrt(2.d0*dble(2*ldeg+1))
c
      do ly=lyup,min0(lyob-1,lylw)
        if(ldeg.le.ldegpsv(ly))then
          kp(ly)=comi/cvp(ly)
          if(cdabs(kp(ly))*rrup(ly).le.xmin)then
            ksmallpsv(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,1)=cldeg*cruplw
            wh(ldeg,ly,1)=-(cldeg+c1)*cruplw
          else if(ksmallpsv(ly))then
            ksmallpsv(ly)=.false.
            call spbjh(ldegpsv(ly),kp(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,1),zjlw(0,ly,1),wj(0,ly,1),
     &                 zhup(0,ly,1),zhlw(0,ly,1),wh(0,ly,1))
          endif
        endif
      enddo
      do ly=lyob,min0(lycm-1,lylw)
        if(ldeg.le.ldegpsv(ly))then
          kp(ly)=comi/cvp(ly)
          ks(ly)=comi/cvs(ly)
          if(cdabs(ks(ly))*rrup(ly).le.xmin)then
            ksmallpsv(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,1)=cldeg*cruplw
            wh(ldeg,ly,1)=-(cldeg+c1)*cruplw
            wj(ldeg,ly,2)=(cldeg+c2)*cruplw
            wh(ldeg,ly,2)=-(cldeg-c1)*cruplw
          else if(ksmallpsv(ly))then
            ksmallpsv(ly)=.false.
            call spbjh(ldegpsv(ly),kp(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,1),zjlw(0,ly,1),wj(0,ly,1),
     &                 zhup(0,ly,1),zhlw(0,ly,1),wh(0,ly,1))
            call spbjh(ldegpsv(ly),ks(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,2),zjlw(0,ly,2),wj(0,ly,2),
     &                 zhup(0,ly,2),zhlw(0,ly,2),wh(0,ly,2))
          endif
        endif
      enddo
      do ly=lycm,min0(lycc-1,lylw)
        if(ldeg.le.ldegpsv(ly))then
          kp(ly)=comi/cvp(ly)
          if(cdabs(kp(ly))*rrup(ly).le.xmin)then
            ksmallpsv(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,1)=cldeg*cruplw
            wh(ldeg,ly,1)=-(cldeg+c1)*cruplw
          else if(ksmallpsv(ly))then
            ksmallpsv(ly)=.false.
            call spbjh(ldegpsv(ly),kp(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,1),zjlw(0,ly,1),wj(0,ly,1),
     &                 zhup(0,ly,1),zhlw(0,ly,1),wh(0,ly,1))
          endif
        endif
      enddo
      do ly=lycc,min0(ly0-1,lylw)
        if(ldeg.le.ldegpsv(ly))then
          kp(ly)=comi/cvp(ly)
          ks(ly)=comi/cvs(ly)
          if(cdabs(ks(ly))*rrup(ly).le.xmin)then
            ksmallpsv(ly)=.true.
            if(rrlw(ly).gt.0.d0)then
              cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
            else
              cruplw=(0.d0,0.d0)
            endif
            wj(ldeg,ly,1)=cldeg*cruplw
            wh(ldeg,ly,1)=-(cldeg+c1)*cruplw
            wj(ldeg,ly,2)=(cldeg+c2)*cruplw
            wh(ldeg,ly,2)=-(cldeg-c1)*cruplw
          else if(ksmallpsv(ly))then
            ksmallpsv(ly)=.false.
            call spbjh(ldegpsv(ly),kp(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,1),zjlw(0,ly,1),wj(0,ly,1),
     &                 zhup(0,ly,1),zhlw(0,ly,1),wh(0,ly,1))
            call spbjh(ldegpsv(ly),ks(ly),rrup(ly),rrlw(ly),
     &                 zjup(0,ly,2),zjlw(0,ly,2),wj(0,ly,2),
     &                 zhup(0,ly,2),zhlw(0,ly,2),wh(0,ly,2))
          endif
        endif
      enddo
c
      if(ldeg.eq.0)then
        do ly=lyup,lylw
          if(ly.lt.lys)then
            lwup=0
            cps(1,ly)=-wj(ldeg,ly,1)
            cps(2,ly)=-wh(ldeg,ly,1)
          else
            lwup=1
            cps(1,ly)=wj(ldeg,ly,1)
            cps(2,ly)=wh(ldeg,ly,1)
          endif
          call qpsmat0(ly,lwup)
        enddo
        call qpsprop0(ypsv,lyup,lylw)
      else
        do ly=lyup,min0(lyob-1,lylw)
          if(rrlw(ly).gt.0.d0)then
            cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
          else
            cruplw=(0.d0,0.d0)
          endif
          lwup=0
          cps(1,ly)=-wj(ldeg,ly,1)
          cps(2,ly)=-wh(ldeg,ly,1)
          cps(3,ly)=-cldeg*cruplw
          cps(4,ly)=(cldeg+c1)*cruplw
          call qpsmatc(ldeg,ly,lylw,lwup)
        enddo
        do ly=lyob,min0(lycm-1,lylw)
          if(rrlw(ly).gt.0.d0)then
            cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
          else
            cruplw=(0.d0,0.d0)
          endif
          if(ly.lt.lys)then
            lwup=0
            cps(1,ly)=-wj(ldeg,ly,1)
            cps(2,ly)=-wh(ldeg,ly,1)
            cps(3,ly)=-wj(ldeg,ly,2)
            cps(4,ly)=-wh(ldeg,ly,2)
            cps(5,ly)=-cldeg*cruplw
            cps(6,ly)=(cldeg+c1)*cruplw
          else
            lwup=1
            cps(1,ly)=wj(ldeg,ly,1)
            cps(2,ly)=wh(ldeg,ly,1)
            cps(3,ly)=wj(ldeg,ly,2)
            cps(4,ly)=wh(ldeg,ly,2)
            cps(5,ly)=cldeg*cruplw
            cps(6,ly)=-(cldeg+c1)*cruplw
          endif
          call qpsmat(ldeg,ly,lylw,lwup)
        enddo
        do ly=lycm,min0(lycc-1,lylw)
          if(rrlw(ly).gt.0.d0)then
            cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
          else
            cruplw=(0.d0,0.d0)
          endif
          lwup=1
          cps(1,ly)=wj(ldeg,ly,1)
          cps(2,ly)=wh(ldeg,ly,1)
          cps(3,ly)=cldeg*cruplw
          cps(4,ly)=-(cldeg+c1)*cruplw
          call qpsmatc(ldeg,ly,lylw,lwup)
        enddo
        do ly=lycc,min0(ly0,lylw)
          if(rrlw(ly).gt.0.d0)then
            cruplw=dcmplx(dlog(rrup(ly)/rrlw(ly)),0.d0)
          else
            cruplw=(0.d0,0.d0)
          endif
          lwup=1
          cps(1,ly)=wj(ldeg,ly,1)
          cps(2,ly)=wh(ldeg,ly,1)
          cps(3,ly)=wj(ldeg,ly,2)
          cps(4,ly)=wh(ldeg,ly,2)
          cps(5,ly)=cldeg*cruplw
          cps(6,ly)=-(cldeg+c1)*cruplw
          call qpsmat(ldeg,ly,lylw,lwup)
        enddo
        if(ldeg.eq.1)then
          call qpsprop1(ypsv,lyup,lylw)
        else
          call qpsprop(ypsv,ldeg,lyup,lylw)
        endif
      endif
      return
      end
