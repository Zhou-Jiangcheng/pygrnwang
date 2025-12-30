      subroutine qpspropg1(ypsv,lyup,lylw)
      use qpalloc
      implicit none
c
c     calculation of response (with gravity) to psv source (l=1)
c     ypsv(6,4): solution vector (complex)
c
      integer*4 lyup,lylw
      complex*16 ypsv(6,4)
c
c     work space
c
      integer*4 i,j,i0,j0,istp,ly,key
      real*8 y4max
      complex*16 cdet,alf,bet,cyabs,cyswap,ca,cb
      complex*16 y0(6,3),yc(6,3),uc(6)
      complex*16 yup(6,3),ylw(6,3),yupc(4,2),ylwc(4,2)
      complex*16 coef6(6,6),b6(6,4),coef4(4,4),b4(4,2)
      external qpdifmat4,qpdifmat6
c
      complex*16 c3
      data c3/(3.d0,0.d0)/
c
c===============================================================================
c
c     propagation from surface to atmosphere/ocean bottom
c
      if(lyob.gt.lyup)then
        do j=1,2
          do i=1,4
            yupc(i,j)=(0.d0,0.d0)
          enddo
        enddo
c
        yupc(1,1)=comi2*crrup(lyup)
        yupc(2,1)=cgrup(lyup)
c
        yupc(2,2)=(1.d0,0.d0)
        yupc(3,2)=-comi2
c
        if(lyr.eq.lyup)then
          do j=1,3
            do i=1,6
              y0(i,j)=(0.d0,0.d0)
            enddo
          enddo
c
          y0(1,1)=yupc(1,1)
          y0(3,1)=yupc(2,1)
c
          y0(3,2)=yupc(2,2)
          y0(5,2)=yupc(3,2)
c
        endif
c
        do j=1,3
          do i=1,6
            yc(i,j)=(0.d0,0.d0)
          enddo
        enddo
c
        yc(1,1)=yupc(1,1)
        yc(3,1)=yupc(2,1)
c
        yc(3,2)=yupc(2,2)
        yc(5,2)=yupc(3,2)
c
        do ly=lyup,min0(lys-1,lyob-1)
          if(ly.gt.lyup.and.ly.eq.lyos)then
c
c           interface ocean-atmosphere
c
            ca=yupc(1,1)*cgrup(ly)/crrup(ly)-yupc(3,1)
            cb=yupc(1,2)*cgrup(ly)/crrup(ly)-yupc(3,2)
            if(cdabs(ca).gt.cdabs(cb))then
              do i=1,4
                yupc(i,2)=yupc(i,1)*cb-yupc(i,2)*ca
              enddo
              yupc(2,2)=yupc(2,2)*crolw(ly-1)/croup(ly)
c
              yupc(1,1)=yupc(1,1)*comi2
              yupc(2,1)=(yupc(2,1)*crolw(ly-1)*comi2
     &                 +ca*(croup(ly)-crolw(ly-1)))/croup(ly)
              yupc(3,1)=yupc(3,1)*comi2
              yupc(4,1)=yupc(4,1)*comi2
              do i=1,6
                yc(i,2)=yc(i,1)*cb-yc(i,2)*ca
                yc(i,1)=yc(i,1)*comi2
              enddo
              if(ly.ge.lyr)then
                do i=1,6
                  y0(i,2)=y0(i,1)*cb-y0(i,2)*ca
                  y0(i,1)=y0(i,1)*comi2
                enddo
              endif
            else
              do i=1,4
                yupc(i,1)=yupc(i,1)*cb-yupc(i,2)*ca
              enddo
              yupc(2,1)=yupc(2,1)*crolw(ly-1)/croup(ly)
c
              yupc(1,2)=yupc(1,2)*comi2
              yupc(2,2)=(yupc(2,2)*crolw(ly-1)*comi2
     &                 +cb*(croup(ly)-crolw(ly-1)))/croup(ly)
              yupc(3,2)=yupc(3,2)*comi2
              yupc(4,2)=yupc(4,2)*comi2
              do i=1,6
                yc(i,1)=yc(i,1)*cb-yc(i,2)*ca
                yc(i,2)=yc(i,2)*comi2
              enddo
              if(ly.ge.lyr)then
                do i=1,6
                  y0(i,1)=y0(i,1)*cb-y0(i,2)*ca
                  y0(i,2)=y0(i,2)*comi2
                enddo
              endif
            endif
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,1)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,1)=yupc(i,1)*cyabs
          enddo
          do i=1,6
            yc(i,1)=yc(i,1)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*cyabs
            enddo
          endif
c
          alf=(0.d0,0.d0)
          do i=1,4
            alf=alf+yupc(i,2)*dconjg(yupc(i,1))/cypnorm(i,ly)**2
          enddo
          do i=1,4
            yupc(i,2)=yupc(i,2)-alf*yupc(i,1)
          enddo
          do i=1,6
            yc(i,2)=yc(i,2)-alf*yc(i,1)
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)-alf*y0(i,1)
            enddo
          endif
c
          cyabs=(0.d0,0.d0)
          do i=1,4
            cyabs=cyabs+yupc(i,2)*dconjg(yupc(i,2))/cypnorm(i,ly)**2
          enddo
          cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
          do i=1,4
            yupc(i,2)=yupc(i,2)*cyabs
          enddo
          do i=1,6
            yc(i,2)=yc(i,2)*cyabs
          enddo
          if(ly.ge.lyr)then
            do i=1,6
              y0(i,2)=y0(i,2)*cyabs
            enddo
          endif
c
          call ruku(yupc,4,2,ly,1,qpdifmat4,rrup(ly),rrlw(ly))
          if(ly.eq.lyr-1)then
            do j=1,3
              do i=1,6
                y0(i,j)=(0.d0,0.d0)
              enddo
            enddo
            do j=1,2
              y0(1,j)=yupc(1,j)
              y0(2,j)=crolw(ly)*crrlw(ly)**2
     &                *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &                -comi2*yupc(2,j)-yupc(3,j))
              y0(3,j)=yupc(2,j)
              y0(5,j)=yupc(3,j)
              y0(6,j)=yupc(4,j)
            enddo
          endif
        enddo
c
        ly=min0(lyob-1,lys-1)
        do j=1,2
          yup(1,j)=yupc(1,j)
          yup(2,j)=crolw(ly)*crrlw(ly)**2
     &            *(yupc(1,j)*cgrlw(ly)/crrlw(ly)
     &            -comi2*yupc(2,j)-yupc(3,j))
          yup(3,j)=(0.d0,0.d0)
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
        if(lyr.eq.lyup)call cmemcpy(yup,y0,18)
      endif
c
c===============================================================================
c
c     propagation from atmosphere/ocean bottom to source
c
      do ly=lyob,lys-1
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+yup(i,1)*dconjg(yup(i,1))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          yup(i,1)=yup(i,1)*cyabs
        enddo
        do i=1,6
          yc(i,1)=yc(i,1)*cyabs
        enddo
        if(ly.ge.lyr)then
          do i=1,6
            y0(i,1)=y0(i,1)*cyabs
          enddo
        endif
c
        alf=(0.d0,0.d0)
        do i=1,6
          alf=alf+yup(i,2)*dconjg(yup(i,1))/cypnorm(i,ly)**2
        enddo
        do i=1,6
          yup(i,2)=yup(i,2)-alf*yup(i,1)
        enddo
        do i=1,6
          yc(i,2)=yc(i,2)-alf*yc(i,1)
        enddo
        if(ly.ge.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)-alf*y0(i,1)
          enddo
        endif
c
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+yup(i,2)*dconjg(yup(i,2))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          yup(i,2)=yup(i,2)*cyabs
        enddo
        do i=1,6
          yc(i,2)=yc(i,2)*cyabs
        enddo
        if(ly.ge.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)*cyabs
          enddo
        endif
c
        alf=(0.d0,0.d0)
        bet=(0.d0,0.d0)
        do i=1,6
          alf=alf+yup(i,3)*dconjg(yup(i,1))/cypnorm(i,ly)**2
          bet=bet+yup(i,3)*dconjg(yup(i,2))/cypnorm(i,ly)**2
        enddo
        do i=1,6
          yup(i,3)=yup(i,3)-alf*yup(i,1)-bet*yup(i,2)
        enddo
        do i=1,6
          yc(i,3)=yc(i,3)-alf*yc(i,1)-bet*yc(i,2)
        enddo
        if(ly.ge.lyr)then
          do i=1,6
            y0(i,3)=y0(i,3)-alf*y0(i,1)-bet*y0(i,2)
          enddo
        endif
c
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+yup(i,3)*dconjg(yup(i,3))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          yup(i,3)=yup(i,3)*cyabs
        enddo
        do i=1,6
          yc(i,3)=yc(i,3)*cyabs
        enddo
        if(ly.ge.lyr)then
          do i=1,6
            y0(i,3)=y0(i,3)*cyabs
          enddo
        endif
c
        call ruku(yup,6,3,ly,1,qpdifmat6,rrup(ly),rrlw(ly))
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
        call qpstart6(1,lylw,ylw)
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          call cmemcpy(ylw,y0,12)
        endif
      endif
c
      do ly=lylw-1,lycc,-1
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          ylw(i,1)=ylw(i,1)*cyabs
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,1)=y0(i,1)*cyabs
          enddo
        endif
c
        alf=(0.d0,0.d0)
        do i=1,6
          alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
        enddo
        do i=1,6
          ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)-alf*y0(i,1)
          enddo
        endif
c
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          ylw(i,2)=ylw(i,2)*cyabs
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)*cyabs
          enddo
        endif
c
        call ruku(ylw,6,2,ly,1,qpdifmat6,rrlw(ly),rrup(ly))
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,12)
      enddo
c===============================================================================
c
c     propagation within outer core
c
      if(lylw.ge.lycc)then
c
c       interface conditions: solid to liquid
c
        ylwc(1,1)=ylw(4,2)*ylw(1,1)-ylw(4,1)*ylw(1,2)
        ylwc(2,1)=ylw(4,2)*ylw(2,1)-ylw(4,1)*ylw(2,2)
        ylwc(3,1)=ylw(4,2)*ylw(5,1)-ylw(4,1)*ylw(5,2)
        ylwc(4,1)=ylw(4,2)*ylw(6,1)-ylw(4,1)*ylw(6,2)
c
c       y2 = Ur
c
        ylwc(2,1)=-ylwc(2,1)/crolw(lycc-1)/crrlw(lycc-1)**2
     &         +ylwc(1,1)*cgrlw(lycc-1)/crrlw(lycc-1)-ylwc(3,1)
        ylwc(1,1)=ylwc(1,1)*comi2
        ylwc(3,1)=ylwc(3,1)*comi2
        ylwc(4,1)=ylwc(4,1)*comi2
        if(lycc.le.lyr.and.lycc.gt.lys)then
          do i=1,6
            y0(i,1)=(ylw(4,2)*y0(i,1)-ylw(4,1)*y0(i,2))*comi2
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      else if(lylw.ge.lycm)then
c
c       lowest layer is within the liquid core
c
        if(cdabs(comi).gt.0.d0)then
c
c         compressible liquid for dynamic response
c
          call qpstart4(1,lylw,ylwc)
          do i=1,4
            ylwc(i,2)=(0.d0,0.d0)
          enddo
        else
c
c         incompressible liquid used for static response
c
          do j=1,2
            do i=1,4
              ylwc(i,j)=(0.d0,0.d0)
            enddo
          enddo
c
          ylwc(3,1)=(1.d0,0.d0)
          ylwc(4,1)=(3.d0,0.d0)
        endif
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          y0(1,1)=ylwc(1,1)
          y0(2,1)=croup(lylw)*crrup(lylw)**2
     &            *(ylwc(1,1)*cgrup(lylw)/crrup(lylw)
     &            -comi2*ylwc(2,1)-ylwc(3,1))
          y0(3,1)=ylwc(2,1)
          y0(4,1)=(0.d0,0.d0)
          y0(5,1)=ylwc(3,1)
          y0(6,1)=ylwc(4,1)
          do i=1,6
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      endif
c
      do ly=min0(lylw-1,lycc-1),lycm,-1
        call ruku(ylwc,4,1,ly,1,qpdifmat4,rrlw(ly),rrup(ly))
        if(ly.eq.lyr.and.ly.gt.lys)then
          y0(1,1)=ylwc(1,1)
          y0(2,1)=croup(ly)*crrup(ly)**2
     &           *(ylwc(1,1)*cgrup(ly)/crrup(ly)
     &           -comi2*ylwc(2,1)-ylwc(3,1))
          y0(3,1)=ylwc(2,1)
          y0(4,1)=(0.d0,0.d0)
          y0(5,1)=ylwc(3,1)
          y0(6,1)=ylwc(4,1)
          do i=1,6
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      enddo
c===============================================================================
c
c     propagation from core-mantle boundary to source or ocean bottom
c
      if(lylw.ge.lycm)then
c
c       interface conditions: liquid to solid
c
        ylw(1,1)=ylwc(1,1)
        ylw(2,1)=croup(lycm)*crrup(lycm)**2
     &          *(ylwc(1,1)*cgrup(lycm)/crrup(lycm)
     &          -comi2*ylwc(2,1)-ylwc(3,1))
        ylw(3,1)=(0.d0,0.d0)
        ylw(4,1)=(0.d0,0.d0)
        ylw(5,1)=ylwc(3,1)
        ylw(6,1)=ylwc(4,1)
c
        ylw(1,2)=(0.d0,0.d0)
        ylw(2,2)=(0.d0,0.d0)
        ylw(3,2)=(1.d0,0.d0)
        ylw(4,2)=(0.d0,0.d0)
        ylw(5,2)=(0.d0,0.d0)
        ylw(6,2)=(0.d0,0.d0)
c
        if(lycm.eq.lyr.and.lycm.gt.lys)then
          call cmemcpy(ylw,y0,12)
        endif
      else if(lylw.ge.lyob)then
c
c       lowest layer is within the mantle
c
        call qpstart6(1,lylw,ylw)
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          call cmemcpy(ylw,y0,12)
        endif
      endif
c
      do ly=min0(lylw-1,lycm-1),max0(lys,lyob),-1
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+ylw(i,1)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          ylw(i,1)=ylw(i,1)*cyabs
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,1)=y0(i,1)*cyabs
          enddo
        endif
c
        alf=(0.d0,0.d0)
        do i=1,6
          alf=alf+ylw(i,2)*dconjg(ylw(i,1))/cypnorm(i,ly)**2
        enddo
        do i=1,6
          ylw(i,2)=ylw(i,2)-alf*ylw(i,1)
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)-alf*y0(i,1)
          enddo
        endif
c
        cyabs=(0.d0,0.d0)
        do i=1,6
          cyabs=cyabs+ylw(i,2)*dconjg(ylw(i,2))/cypnorm(i,ly)**2
        enddo
        cyabs=(1.d0,0.d0)/cdsqrt(cyabs)
        do i=1,6
          ylw(i,2)=ylw(i,2)*cyabs
        enddo
        if(ly.lt.lyr)then
          do i=1,6
            y0(i,2)=y0(i,2)*cyabs
          enddo
        endif
c
        call ruku(ylw,6,2,ly,1,qpdifmat6,rrlw(ly),rrup(ly))
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,12)
      enddo
c
c===============================================================================
c
c     propagation from ocean bottom to source in atmosphere
c
      if(lylw.ge.lyob.and.lys.lt.lyob)then
c
c       interface conditions: solid to liquid
c
        ylwc(1,1)=ylw(4,2)*ylw(1,1)-ylw(4,1)*ylw(1,2)
        ylwc(2,1)=ylw(4,2)*ylw(2,1)-ylw(4,1)*ylw(2,2)
        ylwc(3,1)=ylw(4,2)*ylw(5,1)-ylw(4,1)*ylw(5,2)
        ylwc(4,1)=ylw(4,2)*ylw(6,1)-ylw(4,1)*ylw(6,2)
c
c       y2 = Ur
c
        ylwc(2,1)=-ylwc(2,1)/crolw(lycc-1)/crrlw(lycc-1)**2
     &         +ylwc(1,1)*cgrlw(lycc-1)/crrlw(lycc-1)-ylwc(3,1)
        ylwc(1,1)=ylwc(1,1)*comi2
        ylwc(3,1)=ylwc(3,1)*comi2
        ylwc(4,1)=ylwc(4,1)*comi2
        if(lyob.le.lyr.and.lyob.gt.lys)then
          do i=1,6
            y0(i,1)=ylw(4,2)*y0(i,1)-ylw(4,1)*y0(i,2)
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      else if(lys.lt.lyob)then
c
c       lowest layer is within atmosphere
c
        if(cdabs(comi).gt.0.d0)then
c
c         compressible liquid for dynamic response
c
          call qpstart4(1,lylw,ylwc)
          do i=1,4
            ylwc(i,2)=(0.d0,0.d0)
          enddo
        else
c
c         incompressible liquid used for static response
c
          do j=1,2
            do i=1,4
              ylwc(i,j)=(0.d0,0.d0)
            enddo
          enddo
c
          ylwc(3,1)=(1.d0,0.d0)
          ylwc(4,1)=(3.d0,0.d0)
        endif
c
        if(lylw.eq.lyr.and.lylw.gt.lys)then
          y0(1,1)=ylwc(1,1)
          y0(2,1)=croup(lylw)*crrup(lylw)**2
     &            *(ylwc(1,1)*cgrup(lylw)/crrup(lylw)
     &            -comi2*ylwc(2,1)-ylwc(3,1))
          y0(3,1)=ylwc(2,1)
          y0(4,1)=(0.d0,0.d0)
          y0(5,1)=ylwc(3,1)
          y0(6,1)=ylwc(4,1)
          do i=1,6
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      endif
c
      do ly=min0(lylw-1,lyob-1),lys,-1
        if(ly.eq.lyos-1)then
c
c         interface ocean-atmosphere
c
          ylwc(2,1)=(ylwc(2,1)*croup(ly+1)*comi2
     &          +(cgrlw(ly)*ylwc(1,1)/crrlw(ly)-ylwc(3,1))
     &          *(crolw(ly)-croup(ly+1)))/crolw(ly)
          ylwc(1,1)=ylwc(1,1)*comi2
          ylwc(3,1)=ylwc(3,1)*comi2
          ylwc(4,1)=ylwc(4,1)*comi2
          if(ly.lt.lyr)then
            do i=1,6
              y0(i,1)=y0(i,1)*comi2
            enddo
          endif
        endif
c
        call ruku(ylwc,4,1,ly,1,qpdifmat4,rrlw(ly),rrup(ly))
        if(ly.eq.lyr.and.ly.gt.lys)then
          y0(1,1)=ylwc(1,1)
          y0(2,1)=croup(ly)*crrup(ly)**2
     &            *(ylwc(1,1)*cgrup(ly)/crrup(ly)
     &            -comi2*ylwc(2,1)-ylwc(3,1))
          y0(3,1)=ylwc(2,1)
          y0(4,1)=(0.d0,0.d0)
          y0(5,1)=ylwc(3,1)
          y0(6,1)=ylwc(4,1)
          do i=1,6
            y0(i,2)=(0.d0,0.d0)
          enddo
        endif
      enddo
      if(lys.lt.lyob)then
        ylw(1,1)=ylwc(1,1)
        ylw(2,1)=croup(lys)*crrup(lys)**2
     &          *(ylwc(1,1)*cgrup(lys)/crrup(lys)
     &          -comi2*ylwc(2,1)-ylwc(3,1))
        ylw(3,1)=ylwc(2,1)
        ylw(4,1)=(0.d0,0.d0)
        ylw(5,1)=ylwc(3,1)
        ylw(6,1)=ylwc(4,1)
        do i=1,6
          ylw(i,2)=(0.d0,0.d0)
        enddo
      endif
c
      if(lys.lt.lyob)then
        j0=2
      else
        j0=3
      endif
      do j=1,j0-1
        ylw(1,j)=ylw(1,j)/crrup(lys)
        ylw(2,j)=ylw(2,j)/crrup(lys)**2
        ylw(3,j)=ylw(3,j)/crrup(lys)
        ylw(4,j)=ylw(4,j)/crrup(lys)**2
c        ylw(5,j)=ylw(5,j)
        ylw(6,j)=ylw(6,j)/crrup(lys)
      enddo
      ylw(1,j0)=(1.d0,0.d0)
      ylw(2,j0)=(0.d0,0.d0)
      ylw(3,j0)=(1.d0,0.d0)
      ylw(4,j0)=(0.d0,0.d0)
      ylw(5,j0)=cgrup(lys)-comi2*crrup(lys)
      ylw(6,j0)=-c3*comi2
c
      if(lyr.le.lys)then
        do j=1,j0
          y0(1,j)=y0(1,j)/crrup(lyr)
          y0(2,j)=y0(2,j)/crrup(lyr)**2
          y0(3,j)=y0(3,j)/crrup(lyr)
          y0(4,j)=y0(4,j)/crrup(lyr)**2
c          y0(5,j)=y0(5,j)
          y0(6,j)=y0(6,j)/crrup(lyr)
        enddo
      else
        do j=1,j0-1
          y0(1,j)=y0(1,j)/crrup(lyr)
          y0(2,j)=y0(2,j)/crrup(lyr)**2
          y0(3,j)=y0(3,j)/crrup(lyr)
          y0(4,j)=y0(4,j)/crrup(lyr)**2
c          y0(5,j)=y0(5,j)
          y0(6,j)=y0(6,j)/crrup(lyr)
        enddo
        y0(1,j0)=(1.d0,0.d0)
        y0(2,j0)=(0.d0,0.d0)
        y0(3,j0)=(1.d0,0.d0)
        y0(4,j0)=(0.d0,0.d0)
        y0(5,j0)=cgrup(lyr)-comi2*crrup(lyr)
        y0(6,j0)=-c3*comi2
      endif
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
          print *,' Warning in qpspropg1: anormal exit from cdgemp!'
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
          print *,' Warning in qpspropg1: anormal exit from cdgemp!'
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