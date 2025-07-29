      subroutine qpqmodel(f)
      use qpalloc
      implicit none
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      real*8 f
c
      integer*4 ly
      real*8 fac,pup,plw,sup,slw
      complex*16 cqpup,cqsup,cqplw,cqslw
c
      if(.not.dispersion.or.f.ge.FSBREF)then
        fac=0.d0
      else if(f.le.FSBLW)then
        fac=dlog(FSBLW/FSBREF)/PI
      else
        fac=dlog(f/FSBREF)/PI
      endif
c
      if(1.d0+fac/qsmin.le.0.d0)then
        stop 'Error in qpqmodel: too small Qs value!'
      endif
      do ly=1,ly0
        if(f.gt.0.d0)then
          cqpup=dcmplx(1.d0,0.5d0/qpup(ly))
          cqplw=dcmplx(1.d0,0.5d0/qplw(ly))
        else
          cqpup=(1.d0,0.d0)
          cqplw=(1.d0,0.d0)
        endif
        pup=1.d0+fac/qpup(ly)
        plw=1.d0+fac/qplw(ly)
c
        cvpup(ly)=dcmplx(vpup(ly)*pup/cdabs(cqpup),0.d0)*cqpup
        cvplw(ly)=dcmplx(vplw(ly)*plw/cdabs(cqplw),0.d0)*cqplw
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
c
        if(ly.ge.lyob.and.ly.lt.lycm.or.ly.ge.lycc)then
          sup=1.d0+fac/qsup(ly)
          slw=1.d0+fac/qslw(ly)
          if(f.gt.0.d0)then
            cqsup=dcmplx(1.d0,0.5d0/qsup(ly))
            cqslw=dcmplx(1.d0,0.5d0/qslw(ly))
          else
            cqsup=(1.d0,0.d0)
            cqslw=(1.d0,0.d0)
          endif
          cvsup(ly)=dcmplx(vsup(ly)*sup/cdabs(cqsup),0.d0)*cqsup
          cvslw(ly)=dcmplx(vslw(ly)*slw/cdabs(cqslw),0.d0)*cqslw
          cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        else
          cvsup(ly)=(0.d0,0.d0)
          cvslw(ly)=(0.d0,0.d0)
          cvs(ly)=(0.d0,0.d0)
        endif
c
        cmuup(ly)=croup(ly)*cvsup(ly)**2
        claup(ly)=croup(ly)*cvpup(ly)**2-(2.d0,0.d0)*cmuup(ly)
        cmulw(ly)=crolw(ly)*cvslw(ly)**2
        clalw(ly)=crolw(ly)*cvplw(ly)**2-(2.d0,0.d0)*cmulw(ly)
        cmu(ly)=(0.5d0,0.d0)*(cmuup(ly)+cmulw(ly))
        cla(ly)=(0.5d0,0.d0)*(claup(ly)+clalw(ly))
      enddo
      return
      end
