      subroutine qpqmodel(f)
      implicit none
c
      include 'qpglobal.h'
c
c     calculate q based on the constant q model
c
c     f = frequency
c
      double precision f
c
      integer ly
      double precision fac,pup,plw,sup,slw,eup,elw
c
      if(f.le.0.d0)return
      if(.not.dispersion)then
        fac=0.d0
      else if(f.le.FSBLW)then
        fac=dlog(FSBLW/FSBREF)/PI
      else if(f.ge.FSBUP)then
        fac=dlog(FSBUP/FSBREF)/PI
      else
        fac=dlog(f/FSBREF)/PI
      endif
c
      if(1.d0+fac/qsmin.le.0.d0)then
        stop ' Error in qpqmodel: too small Qs value!'
      endif
      do ly=1,ly0
        if(ly.ge.lyob.and.ly.lt.lycm.or.ly.ge.lycc)then
          sup=1.d0+fac/qsup(ly)
          slw=1.d0+fac/qslw(ly)
          eup=(vsup(ly)/vpup(ly))**2*4.d0/3.d0
          elw=(vslw(ly)/vplw(ly))**2*4.d0/3.d0
          pup=1.d0+fac*eup/qsup(ly)
          plw=1.d0+fac*elw/qslw(ly)
          cvsup(ly)=dcmplx(vsup(ly)*sup,0.d0)
     &             *dcmplx(1.d0,0.5d0/qsup(ly))
          cvslw(ly)=dcmplx(vslw(ly)*slw,0.d0)
     &             *dcmplx(1.d0,0.5d0/qslw(ly))
          cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        else
          pup=1.d0
          plw=1.d0
          cvsup(ly)=(0.d0,0.d0)
          cvslw(ly)=(0.d0,0.d0)
          cvs(ly)=(0.d0,0.d0)
        endif
        cvpup(ly)=dcmplx(vpup(ly)*pup,0.d0)
     &           *dcmplx(1.d0,0.5d0/qpup(ly))
        cvplw(ly)=dcmplx(vplw(ly)*plw,0.d0)
     &           *dcmplx(1.d0,0.5d0/qplw(ly))
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
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
