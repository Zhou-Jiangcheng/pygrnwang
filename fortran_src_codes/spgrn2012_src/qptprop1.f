      subroutine qptprop1(ysh,lylw)
      use qpalloc
      implicit none
c
      include 'qpglobal.h'
c
c     calculation of response to sh source
c     ysh(2,2): solution vector (complex)
c
      integer lylw
      double complex ysh(2,2)
c
c     work space
c
      integer i,istp,ly,key
      double complex amom,amomup,amomlw
      double complex y0(2),yup(2),ylw(2),c(2)
      double complex coef(2,2),b(2,2)
c
c===============================================================================
c
c     propagation from surface to source
c
      if(lys.lt.lyob)then
        do istp=1,2
          do i=1,2
            ysh(i,istp)=(0.d0,0.d0)
          enddo
        enddo
        return
      endif
c
      yup(1)=(1.d0,0.d0)
      yup(2)=(0.d0,0.d0)
c
      if(lyr.eq.lyob)call cmemcpy(yup,y0,2)
c
      amomup=(0.d0,0.d0)
      do ly=lyob,lys-1
        amomup=amomup+yup(1)*dcmplx(0.5d0*mshell(ly),0.d0)/crrup(ly)**2
        call caxcb(mat2x2inv(1,1,ly),yup,2,2,1,c)
c
	  c(1)=c(1)*cdexp(cpt(1,ly))
        c(2)=c(2)*cdexp(cpt(2,ly))
c
        call caxcb(mat2x2lw(1,1,ly),c,2,2,1,yup)
        amomup=amomup+yup(1)*dcmplx(0.5d0*mshell(ly),0.d0)/crrlw(ly)**2
        if(ly.eq.lyr-1)call cmemcpy(yup,y0,2)
      enddo
      yup(1)=yup(1)/crrup(lys)
      yup(2)=yup(2)/crrup(lys)**2
c
c===============================================================================
c
c     propagation from bottom to source
c
      if(lylw.eq.lycm)then
        ylw(1)=(1.d0,0.d0)
        ylw(2)=(0.d0,0.d0)
      else
        ylw(1)=mat2x2up(1,1,lylw)
        ylw(2)=mat2x2up(2,1,lylw)
      endif
      if(lylw.eq.lyr.and.lylw.gt.lys)then
        call cmemcpy(ylw,y0,2)
      endif
c
      amomlw=(0.d0,0.d0)
      do ly=lylw-1,lys,-1
        amomlw=amomlw+ylw(1)*dcmplx(0.5d0*mshell(ly),0.d0)/crrlw(ly)**2
        call caxcb(mat2x2inv(1,1,ly),ylw,2,2,1,c)
c
        c(1)=c(1)*cdexp(cpt(1,ly))
        c(2)=c(2)*cdexp(cpt(2,ly))
c
        call caxcb(mat2x2up(1,1,ly),c,2,2,1,ylw)
        amomlw=amomlw+ylw(1)*dcmplx(0.5d0*mshell(ly),0.d0)/crrup(ly)**2
        if(ly.eq.lyr.and.ly.gt.lys)call cmemcpy(ylw,y0,2)
      enddo
      ylw(1)=ylw(1)/crrup(lys)
      ylw(2)=ylw(2)/crrup(lys)**2
c
      y0(1)=y0(1)/crrup(lyr)
      y0(2)=y0(2)/crrup(lyr)**2
c
c===============================================================================
c     source function
c===============================================================================
c
      b(1,1)=(1.d0,0.d0)
      b(2,1)=(0.d0,0.d0)
      b(1,2)=(0.d0,0.d0)
      b(2,2)=(1.d0,0.d0)
      do i=1,2
        coef(i,1)=yup(i)
        coef(i,2)=-ylw(i)
      enddo
      key=0
      call cdsvd500(coef,b,2,2,0.d0,key)
      if(key.eq.0)then
        print *,' Warning in qptprop1: anormal exit from cdgemp!'
        return
      endif
      if(lyr.le.lys)then
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(1,istp)*y0(i)
          enddo
        enddo
      else
        do istp=1,2
          do i=1,2
            ysh(i,istp)=b(2,istp)*y0(i)
          enddo
        enddo
      endif
c
c     conservation of angular momentum
c
      do istp=1,2
        amom=amomup*b(1,istp)+amomlw*b(2,istp)
        ysh(1,istp)=ysh(1,istp)-amom*crrup(lyr)/dcmplx(mmantle,0.d0)
      enddo
      return
      end