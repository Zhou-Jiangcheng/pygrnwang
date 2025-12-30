      subroutine qpmaxdeg(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
c     work space
c
      integer*4 ig,ly
      real*8 ray,raymin,raymax
c
      if(fullwave)then
        ldegmin=ldeg0
        raymax=0.d0
        do ig=1,ngrn
          lys=lygrn(ig)
          raymin=0.d0
          do ly=min0(lyr,lys),max0(lyr,lys)
            if(vsup(ly).gt.0.d0)then
              ray=rrup(ly)/vsup(ly)
            else
              ray=rrup(ly)/vpup(ly)
            endif
            if(raymin.le.0.d0.or.raymin.gt.ray)raymin=ray
            if(ly.gt.1)then
              if(vslw(ly-1).gt.0.d0)then
                ray=rrup(ly)/vslw(ly-1)
              else
                ray=rrup(ly)/vplw(ly-1)
              endif
              if(raymin.le.0.d0.or.raymin.gt.ray)raymin=ray
            endif
          enddo
          slwmax=dmax1(slwmax,raymin/rrup(lys))
        enddo
        slwmax=slwfac*slwmax
      else
        ldegmin=ldeg1+ndmax
      endif
c
      ldegmax=1+ldegmin+idnint(REARTH*PI2*fcut*slwmax)
c
      allocate(yr(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpmaxdeg: yr not allocated!'
      allocate(yt(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpmaxdeg: yp not allocated!'
      allocate(yp(0:ldegmax,4,0:ndmax),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpmaxdeg: yt not allocated!'
c
      return
      end
