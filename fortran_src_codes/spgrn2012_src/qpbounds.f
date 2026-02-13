      subroutine qpbounds(ldegcut)
      implicit none
      integer ldegcut
c
      include 'qpglobal.h'
c
c     work space
c
      integer ig,il,ldeg,ldeg0,ldeg1,ldegup,ly
      double precision dxs,dys2,omi,dll1,ksp2,expo,rr0min,depsmin
c
      rr0min=rr0(1)
      depsmin=grndep(1)
      lys=1
      do ig=2,ngrn
        if(rr0min.gt.rr0(ig))then
          rr0min=rr0(ig)
        endif
        if(depsmin.gt.grndep(ig))then
          lys=lygrn(ig)
          depsmin=grndep(ig)
        endif
      enddo
c
c     tappered disk source:
c     2*[cos(theta)-cos(alpha)]/[1-cos(alpha)]^2
c
      if(rr0min.le.0.d0)then
        do ldeg=0,ldegmax
          fdisk(ldeg)=1.d0
        enddo
      else
        dxs=dcos(rr0min/rrup(lys))
        dys2=(dsin(0.5d0*rr0min/rrup(lys)))**2
        call legendre(ldegmax,rr0min/rrup(lys),dplm1)
        fdisk(0)=1.d0
        fdisk(1)=(2.d0+dxs)/3.d0
        do ldeg=2,ldegmax
          fdisk(ldeg)=0.5d0/dys2**2
     &             /dble(ldeg)/dble(ldeg+1)/dble(ldeg+2)
     &             *(2.d0*dys2*(2.d0*dble(ldeg)*dys2+1.d0)
     &              +dble(ldeg)*dxs**2*dplm1(ldeg)
     &              -dble(2*ldeg+1)*dxs*dplm1(ldeg-1)
     &              +dble(ldeg+1)*dplm1(ldeg-2))
        enddo
c
        do ldeg=0,ldegmax
          fdisk(ldeg)=dabs(fdisk(ldeg))
          do il=ldeg+1,ldegmax
            fdisk(ldeg)=dmax1(fdisk(ldeg),dabs(fdisk(il)))
          enddo
        enddo
      endif
c
      ldegcut=0
      do ig=1,ngrn
        lys=lygrn(ig)
c
        omi=PI2*fcut
        ldeg0=ldegmin+ndmax
        if(slwmax.le.0.d0)then
          ldeg1=ldegmax
        else
          ldeg1=min0(ldegmax,idint(REARTH*omi*slwmax))
        endif
        ldegup=ldeg0
        do ldegup=ldeg0,ldeg1
          dll1=dble(ldegup)*dble(ldegup+1)
          expo=0.d0
          do ly=min0(lys,lyr),max0(lys,lyr)-1
            if(vsup(ly).gt.0.d0)then
              ksp2=(omi*rrup(ly)/vsup(ly))**2
            else
              ksp2=(omi*rrup(ly)/vpup(ly))**2
            endif
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
          enddo
          if(expo.gt.expos+dlog(fdisk(ldegup)))goto 20
        enddo
20      continue
        ldegcut=max0(ldegcut,ldegup)
      enddo
      ldegcut=ldegcut+1+ndmax
      if(ldegcut.ge.ldegmax)then
        write(*,'(2(a,i6))')' Warning in qpbounds: '
     &     //'ldegmax reguired/defined:',ldegcut,'/',ldegmax
        ldegcut=ldegmax
      endif
      return
      end
