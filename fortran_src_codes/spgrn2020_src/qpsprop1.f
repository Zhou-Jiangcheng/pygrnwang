      subroutine qpsprop1(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of response to psv source (l=1)
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,j,j0,istp,ly,key
      real*8 y4max
      complex*16 cdet,cyswap
      complex*16 c0(6,3),c1(6,3),cc0(4,2),cc1(4,2)
      complex*16 y0(6,3),y1(6,3),yc(6,3),uc(6)
      complex*16 yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      complex*16 wave(6),orth(3,3),orthc(2,2)
      complex*16 coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
c
      complex*16 c3
      data c3/(3.d0,0.d0)/
c
c===============================================================================
c
c     propagation from surface to atmosphere/ocean bottom
c
      if(lyr.lt.lyob.and.dreal(comi).le.0.d0)return
      if(lyob.gt.1.and.dreal(comi).gt.0.d0)then
        do j=1,2
          do i=1,4
            yupc(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yupc(1,1)=(1.d0,0.d0)
        yupc(3,2)=(1.d0,0.d0)
        do j=1,3
          do i=1,6
            yc(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yc(1,1)=yupc(1,1)
        yc(5,2)=yupc(3,2)
        if(lyr.eq.lyup)call cmemcpy(yc,y0,18)
c
        do ly=lyup,min0(lys-1,lyob-1)
          do j=1,3,2
            wave(j)=cdexp(cps(j,ly))
          enddo
          do j=2,4,2
            wave(j)=cdexp(-cps(j,ly))
          enddo
c
          call caxcb(mas4x4inv(1,1,ly),yupc,4,4,2,cc0)
c
c         orthonormalization of the p-sv modes
c
          cdet=cc0(2,1)*cc0(4,2)-cc0(4,1)*cc0(2,2)
          orthc(1,1)=cc0(4,2)/cdet
          orthc(1,2)=-cc0(2,2)/cdet
          orthc(2,1)=-cc0(4,1)/cdet
          orthc(2,2)=cc0(2,1)/cdet
c
          call caxcb(cc0,orthc,4,2,2,cc1)
c
          call caxcb(yc,orthc,6,2,2,y1)
          call cmemcpy(y1,yc,12)
          do j=1,2
            do i=1,6
              yc(i,j)=yc(i,j)*wave(2*j)
            enddo
          enddo
c
          if(ly.ge.lyr)then
c
c           orthonormalization of the receiver vectors
c
            call caxcb(y0,orthc,6,2,2,y1)
            call cmemcpy(y1,y0,12)
            do j=1,2
              do i=1,6
                y0(i,j)=y0(i,j)*wave(2*j)
              enddo
            enddo
          endif
c
	    cc1(1,1)=cc1(1,1)*wave(1)*wave(2)
 	    cc1(2,1)=(1.d0,0.d0)
	    cc1(3,1)=cc1(3,1)*wave(3)*wave(2)
          cc1(4,1)=(0.d0,0.d0)
c
          cc1(1,2)=cc1(1,2)*wave(1)*wave(4)
          cc1(2,2)=(0.d0,0.d0)
          cc1(3,2)=cc1(3,2)*wave(3)*wave(4)
          cc1(4,2)=(1.d0,0.d0)
c
          call caxcb(mas4x4lw(1,1,ly),cc1,4,4,2,yupc)
c
          if(ly.eq.lyr-1)then
            do j=1,2
              y0(1,j)=yupc(1,j)
              y0(2,j)=yupc(2,j)
              y0(3,j)=-yupc(2,j)/(crolw(ly)*comi2*crrlw(ly)**2)
              y0(4,j)=(0.d0,0.d0)
              y0(5,j)=yupc(3,j)
              y0(6,j)=yupc(4,j)
            enddo
            do i=1,6
              y0(i,3)=(0.d0,0.d0)
            enddo
          endif
        enddo
c
        ly=min0(lys-1,lyob-1)
        do j=1,2
          yup(1,j)=yupc(1,j)
          yup(2,j)=yupc(2,j)
          yup(3,j)=-yupc(2,j)/(crolw(ly)*comi2*crrlw(ly)**2)
          yup(4,j)=(0.d0,0.d0)
          yup(5,j)=yupc(3,j)
          yup(6,j)=yupc(4,j)
        enddo
        if(lys.ge.lyob)then
          yup(1,3)=(0.d0,0.d0)
          yup(2,3)=(0.d0,0.d0)
          yup(3,3)=(1.d0,0.d0)
          yup(4,3)=(0.d0,0.d0)
          yup(5,3)=(0.d0,0.d0)
          yup(6,3)=(0.d0,0.d0)
        endif
        if(lyr.eq.lyob.and.ly.eq.lyob-1)call cmemcpy(yup,y0,18)
      else
        do j=1,3
          do i=1,6
            yup(i,j)=(0.d0,0.d0)
          enddo
        enddo
        yup(1,1)=(1.d0,0.d0)
        yup(3,2)=(1.d0,0.d0)
        yup(5,3)=(1.d0,0.d0)
        call cmemcpy(yup,yc,18)
      endif
c
      if(lyr.eq.lyob)call cmemcpy(yup,y0,18)
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source
c
      do ly=lyob,lys-1
        do j=1,5,2
          wave(j)=cdexp(cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),yup,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(2,1)*c0(4,2)*c0(6,3)
     &      +c0(4,1)*c0(6,2)*c0(2,3)
     &      +c0(6,1)*c0(2,2)*c0(4,3)
     &      -c0(6,1)*c0(4,2)*c0(2,3)
     &      -c0(4,1)*c0(2,2)*c0(6,3)
     &      -c0(2,1)*c0(6,2)*c0(4,3)
        orth(1,1)=(c0(4,2)*c0(6,3)-c0(4,3)*c0(6,2))/cdet
        orth(2,1)=(c0(4,3)*c0(6,1)-c0(4,1)*c0(6,3))/cdet
        orth(3,1)=(c0(4,1)*c0(6,2)-c0(4,2)*c0(6,1))/cdet
        orth(1,2)=(c0(2,3)*c0(6,2)-c0(2,2)*c0(6,3))/cdet
        orth(2,2)=(c0(2,1)*c0(6,3)-c0(2,3)*c0(6,1))/cdet
        orth(3,2)=(c0(2,2)*c0(6,1)-c0(2,1)*c0(6,2))/cdet
        orth(1,3)=(c0(2,2)*c0(4,3)-c0(2,3)*c0(4,2))/cdet
        orth(2,3)=(c0(2,3)*c0(4,1)-c0(2,1)*c0(4,3))/cdet
        orth(3,3)=(c0(2,1)*c0(4,2)-c0(2,2)*c0(4,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
c
        call caxcb(yc,orth,6,3,3,y1)
        call cmemcpy(y1,yc,18)
        do j=1,3
          do i=1,6
            yc(i,j)=yc(i,j)*wave(2*j)
          enddo
        enddo
c
        if(ly.ge.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j)
            enddo
          enddo
        endif
c
	  c1(1,1)=c1(1,1)*wave(1)*wave(2)
 	  c1(2,1)=(1.d0,0.d0)
	  c1(3,1)=c1(3,1)*wave(3)*wave(2)
        c1(4,1)=(0.d0,0.d0)
	  c1(5,1)=c1(5,1)*wave(5)*wave(2)
        c1(6,1)=(0.d0,0.d0)
c
        c1(1,2)=c1(1,2)*wave(1)*wave(4)
        c1(2,2)=(0.d0,0.d0)
        c1(3,2)=c1(3,2)*wave(3)*wave(4)
        c1(4,2)=(1.d0,0.d0)
        c1(5,2)=c1(5,2)*wave(5)*wave(4)
        c1(6,2)=(0.d0,0.d0)
c
        c1(1,3)=c1(1,3)*wave(1)*wave(6)
        c1(2,3)=(0.d0,0.d0)
        c1(3,3)=c1(3,3)*wave(3)*wave(6)
        c1(4,3)=(0.d0,0.d0)
        c1(5,3)=c1(5,3)*wave(5)*wave(6)
        c1(6,3)=(1.d0,0.d0)
c
        call caxcb(mas6x6lw(1,1,ly),c1,6,6,3,yup)
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,18)
      enddo
c
      do j=1,3
        yup(1,j)=yup(1,j)/crrup(lys)
        yup(2,j)=yup(2,j)/crrup(lys)**2
        yup(3,j)=yup(3,j)/crrup(lys)
        yup(4,j)=yup(4,j)/crrup(lys)**2
c        yup(5,j)=yup(5,j)
        yup(6,j)=yup(6,j)/crrup(lys)
c
        yc(1,j)=yc(1,j)/crrup(lyup)
        yc(2,j)=yc(2,j)/crrup(lyup)**2
        yc(3,j)=yc(3,j)/crrup(lyup)
        yc(4,j)=yc(4,j)/crrup(lyup)**2
c        yc(5,j)=yc(5,j)
        yc(6,j)=yc(6,j)/crrup(lyup)
      enddo
c
c===============================================================================
c
c     propagation within inner core
c
      if(lylw.ge.lycc)then
c
c       lowest layer is within inner core
c
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6up(i,2*j-1,lylw)
          enddo
        enddo
      endif
c
      if(lylw.eq.lyr.and.lylw.gt.lys)call cmemcpy(ylw,y0,18)
c
      do ly=lylw-1,lycc,-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c===============================================================================
c
c     propagation within outer core
c
      if(lylw.ge.lycc)then
c
c       interface conditions: solid to liquid
c
        y4max=cdabs(ylw(4,3))
        j0=3
        do j=1,2
          if(y4max.lt.cdabs(ylw(4,j)))then
            y4max=cdabs(ylw(4,j))
            j0=j
          endif
        enddo
        do i=1,6
          cyswap=ylw(i,j0)
          ylw(i,j0)=ylw(i,3)
          ylw(i,3)=cyswap
        enddo
        do j=1,2
          ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
          ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
          ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
          ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
        enddo
        if(lycc.le.lyr.and.lycc.gt.lys)then
          do i=1,6
            cyswap=y0(i,j0)
            y0(i,j0)=y0(i,3)
            y0(i,3)=cyswap
          enddo
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
            enddo
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      else if(lylw.ge.lycm)then
c
c       lowest layer is within the liquid core
c
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4up(i,2*j-1,lylw)
          enddo
        enddo
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            if(cdabs(comi2).gt.0.d0)then
              y0(3,j)=-ylwc(2,j)/(croup(lylw)*comi2*crrup(lylw)**2)
            else
              y0(3,j)=(0.d0,0.d0)
            endif
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      endif
c
      do ly=min0(lylw-1,lycc-1),lycm,-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            if(cdabs(comi2).gt.0.d0)then
              y0(3,j)=-ylwc(2,j)/(croup(ly)*comi2*crrup(ly)**2)
            else
              y0(3,j)=(0.d0,0.d0)
            endif
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      enddo
c
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lylw.ge.lycm)then
c
c       interface conditions: liquid to solid
c
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=ylwc(2,j)
          if(cdabs(comi2).gt.0.d0)then
            ylw(3,j)=-ylwc(2,j)/(croup(lycm)*comi2*crrup(lycm)**2)
          else
            ylw(3,j)=(0.d0,0.d0)
          endif
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        ylw(1,3)=(0.d0,0.d0)
        ylw(2,3)=(0.d0,0.d0)
        ylw(3,3)=(1.d0,0.d0)
        ylw(4,3)=(0.d0,0.d0)
        ylw(5,3)=(0.d0,0.d0)
        ylw(6,3)=(0.d0,0.d0)
c
        if(lycm.eq.lyr.and.lycm.gt.lys)then
          call cmemcpy(ylw,y0,18)
        endif
      else if(lylw.ge.lyob)then
c
c       lowest layer is within the mantle
c
        do j=1,3
          do i=1,6
            ylw(i,j)=mas6x6up(i,2*j-1,lylw)
          enddo
        enddo
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          call cmemcpy(ylw,y0,18)
        endif
      endif
c
      do ly=min0(lylw-1,lycm-1),max0(lys,lyob),-1
        do j=1,5,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,6,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas6x6inv(1,1,ly),ylw,6,6,3,c0)
c
c       orthonormalization of the p-sv modes
c
        cdet=c0(1,1)*c0(3,2)*c0(5,3)
     &      +c0(3,1)*c0(5,2)*c0(1,3)
     &      +c0(5,1)*c0(1,2)*c0(3,3)
     &      -c0(5,1)*c0(3,2)*c0(1,3)
     &      -c0(3,1)*c0(1,2)*c0(5,3)
     &      -c0(1,1)*c0(5,2)*c0(3,3)
        orth(1,1)=(c0(3,2)*c0(5,3)-c0(3,3)*c0(5,2))/cdet
        orth(2,1)=(c0(3,3)*c0(5,1)-c0(3,1)*c0(5,3))/cdet
        orth(3,1)=(c0(3,1)*c0(5,2)-c0(3,2)*c0(5,1))/cdet
        orth(1,2)=(c0(1,3)*c0(5,2)-c0(1,2)*c0(5,3))/cdet
        orth(2,2)=(c0(1,1)*c0(5,3)-c0(1,3)*c0(5,1))/cdet
        orth(3,2)=(c0(1,2)*c0(5,1)-c0(1,1)*c0(5,2))/cdet
        orth(1,3)=(c0(1,2)*c0(3,3)-c0(1,3)*c0(3,2))/cdet
        orth(2,3)=(c0(1,3)*c0(3,1)-c0(1,1)*c0(3,3))/cdet
        orth(3,3)=(c0(1,1)*c0(3,2)-c0(1,2)*c0(3,1))/cdet
c
        call caxcb(c0,orth,6,3,3,c1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orth,6,3,3,y1)
          call cmemcpy(y1,y0,18)
          do j=1,3
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
        c1(1,1)=(1.d0,0.d0)
        c1(2,1)=c1(2,1)*wave(2)*wave(1)
        c1(3,1)=(0.d0,0.d0)
        c1(4,1)=c1(4,1)*wave(4)*wave(1)
        c1(5,1)=(0.d0,0.d0)
        c1(6,1)=c1(6,1)*wave(6)*wave(1)
c
        c1(1,2)=(0.d0,0.d0)
        c1(2,2)=c1(2,2)*wave(2)*wave(3)
        c1(3,2)=(1.d0,0.d0)
        c1(4,2)=c1(4,2)*wave(4)*wave(3)
        c1(5,2)=(0.d0,0.d0)
        c1(6,2)=c1(6,2)*wave(6)*wave(3)
c
        c1(1,3)=(0.d0,0.d0)
        c1(2,3)=c1(2,3)*wave(2)*wave(5)
        c1(3,3)=(0.d0,0.d0)
        c1(4,3)=c1(4,3)*wave(4)*wave(5)
        c1(5,3)=(1.d0,0.d0)
        c1(6,3)=c1(6,3)*wave(6)*wave(5)
c
        call caxcb(mas6x6up(1,1,ly),c1,6,6,3,ylw)
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,18)
      enddo
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lylw.ge.lyob.and.lys.lt.lyob)then
c
c       interface conditions: solid to liquid
c
        y4max=cdabs(ylw(4,3))
        j0=3
        do j=1,2
          if(y4max.lt.cdabs(ylw(4,j)))then
            y4max=cdabs(ylw(4,j))
            j0=j
          endif
        enddo
        do i=1,6
          cyswap=ylw(i,j0)
          ylw(i,j0)=ylw(i,3)
          ylw(i,3)=cyswap
        enddo
        do j=1,2
          ylwc(1,j)=ylw(1,j)-ylw(4,j)*ylw(1,3)/ylw(4,3)
          ylwc(2,j)=ylw(2,j)-ylw(4,j)*ylw(2,3)/ylw(4,3)
          ylwc(3,j)=ylw(5,j)-ylw(4,j)*ylw(5,3)/ylw(4,3)
          ylwc(4,j)=ylw(6,j)-ylw(4,j)*ylw(6,3)/ylw(4,3)
        enddo
        if(lyob.le.lyr.and.lyob.gt.lys)then
          do i=1,6
            cyswap=y0(i,j0)
            y0(i,j0)=y0(i,3)
            y0(i,3)=cyswap
          enddo
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)-ylw(4,j)*y0(i,3)/ylw(4,3)
            enddo
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      else if(lys.lt.lyob)then
c
c       lowest layer is within the atmosphere
c
        do j=1,2
          do i=1,4
            ylwc(i,j)=mas4x4up(i,2*j-1,lylw)
          enddo
        enddo
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            if(cdabs(comi2).gt.0.d0)then
              y0(3,j)=-ylwc(2,j)/(croup(lylw)*comi2*crrup(lylw)**2)
            else
              y0(3,j)=(0.d0,0.d0)
            endif
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      endif
c
      do ly=min0(lylw-1,lyob-1),lys,-1
        do j=1,3,2
          wave(j)=cdexp(-cps(j,ly))
        enddo
        do j=2,4,2
          wave(j)=cdexp(cps(j,ly))
        enddo
c
        call caxcb(mas4x4inv(1,1,ly),ylwc,4,4,2,cc0)
c
c       orthonormalization of the p-sv modes
c
        cdet=cc0(1,1)*cc0(3,2)-cc0(3,1)*cc0(1,2)
        orthc(1,1)=cc0(3,2)/cdet
        orthc(1,2)=-cc0(1,2)/cdet
        orthc(2,1)=-cc0(3,1)/cdet
        orthc(2,2)=cc0(1,1)/cdet
c
        call caxcb(cc0,orthc,4,2,2,cc1)
        if(ly.lt.lyr)then
c
c         orthonormalization of the receiver vectors
c
          call caxcb(y0,orthc,6,2,2,y1)
          call cmemcpy(y1,y0,12)
          do j=1,2
            do i=1,6
              y0(i,j)=y0(i,j)*wave(2*j-1)
            enddo
          enddo
        endif
c
        cc1(1,1)=(1.d0,0.d0)
        cc1(2,1)=cc1(2,1)*wave(2)*wave(1)
        cc1(3,1)=(0.d0,0.d0)
        cc1(4,1)=cc1(4,1)*wave(4)*wave(1)
c
        cc1(1,2)=(0.d0,0.d0)
        cc1(2,2)=cc1(2,2)*wave(2)*wave(3)
        cc1(3,2)=(1.d0,0.d0)
        cc1(4,2)=cc1(4,2)*wave(4)*wave(3)
c
        call caxcb(mas4x4up(1,1,ly),cc1,4,4,2,ylwc)
        if(ly.eq.lyr.and.ly.gt.lys)then
          do j=1,2
            y0(1,j)=ylwc(1,j)
            y0(2,j)=ylwc(2,j)
            if(cdabs(comi2).gt.0.d0)then
              y0(3,j)=-ylwc(2,j)/(croup(ly)*comi2*crrup(ly)**2)
            else
              y0(3,j)=(0.d0,0.d0)
            endif
            y0(4,j)=(0.d0,0.d0)
            y0(5,j)=ylwc(3,j)
            y0(6,j)=ylwc(4,j)
          enddo
          do i=1,6
            y0(i,3)=(0.d0,0.d0)
          enddo
        endif
      enddo
      if(lys.lt.lyob)then
        do j=1,2
          ylw(1,j)=ylwc(1,j)
          ylw(2,j)=ylwc(2,j)
          if(cdabs(comi2).gt.0.d0)then
            ylw(3,j)=-ylwc(2,j)/(croup(lys)*comi2*crrup(lys)**2)
          else
            ylw(3,j)=(0.d0,0.d0)
          endif
          ylw(4,j)=(0.d0,0.d0)
          ylw(5,j)=ylwc(3,j)
          ylw(6,j)=ylwc(4,j)
        enddo
        do i=1,6
          ylw(i,3)=(0.d0,0.d0)
        enddo
      endif
c
      do j=1,3
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
c
      do j=1,3
        y0(1,j)=y0(1,j)/crrup(lyr)
        y0(2,j)=y0(2,j)/crrup(lyr)**2
        y0(3,j)=y0(3,j)/crrup(lyr)
        y0(4,j)=y0(4,j)/crrup(lyr)**2
c        y0(5,j)=y0(5,j)
        y0(6,j)=y0(6,j)/crrup(lyr)
      enddo
c
c===============================================================================
c     source function
c===============================================================================
c
      if(lys.lt.lyob)then
        do istp=1,2
          do i=1,4
            b4(i,istp)=(0.d0,0.d0)
          enddo
          b4(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,2
          do i=1,2
            coef4(i,j)=yup(i,j)
            coef4(i,j+2)=-ylw(i,j)
          enddo
          do i=3,4
            coef4(i,j)=yup(i+2,j)
            coef4(i,j+2)=-ylw(i+2,j)
          enddo
        enddo
        key=0
        call cdsvd(coef4,b4,4,2,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop1: anormal exit from cdgemp!'
          return
        endif
c
c       uc = translational motion of earth's centre of mass
c
        do istp=1,2
          uc(istp)=(0.d0,0.d0)
          do j=1,2
            uc(istp)=uc(istp)+b4(j,istp)*yc(5,j)
          enddo
          uc(istp)=uc(istp)/cgrup(lyup)
        enddo
        if(lyr.le.lys)then
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,2
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,2
                ypsv(i,istp)=ypsv(i,istp)+b4(j+2,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
        do istp=3,4
          uc(istp)=(0.d0,0.d0)
          do i=1,6
            ypsv(i,istp)=(0.d0,0.d0)
          enddo
        enddo
      else
        do istp=1,4
          do i=1,6
            b6(i,istp)=(0.d0,0.d0)
          enddo
          b6(istp,istp)=(1.d0,0.d0)
        enddo
        do j=1,3
          do i=1,6
            coef6(i,j)=yup(i,j)
            coef6(i,j+3)=-ylw(i,j)
          enddo
        enddo
        key=0
        call cdsvd(coef6,b6,6,4,0.d0,key)
        if(key.eq.0)then
          print *,' Warning in qpsprop1: anormal exit from cdgemp!'
          return
        endif
c
c       uc = translational motion of earth's centre of mass
c
        do istp=1,4
          uc(istp)=(0.d0,0.d0)
          do j=1,3
            uc(istp)=uc(istp)+b6(j,istp)*yc(5,j)
          enddo
          uc(istp)=uc(istp)/cgrup(lyup)
        enddo
        if(lyr.le.lys)then
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j,istp)*y0(i,j)
              enddo
            enddo
          enddo
        else
          do istp=1,4
            do i=1,6
              ypsv(i,istp)=(0.d0,0.d0)
              do j=1,3
                ypsv(i,istp)=ypsv(i,istp)+b6(j+3,istp)*y0(i,j)
              enddo
            enddo
          enddo
        endif
      endif
c
c     transform to motion relative to the centre of mass
c
      do istp=1,4
        ypsv(1,istp)=ypsv(1,istp)-uc(istp)
        ypsv(3,istp)=ypsv(3,istp)-uc(istp)
        ypsv(5,istp)=ypsv(5,istp)-uc(istp)*cgrup(lyr)
      enddo
      return
      end