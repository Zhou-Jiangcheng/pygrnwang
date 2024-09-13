      subroutine qpgrnspec(ig)
      use qpalloc
      implicit none
      integer*4 ig
c
      integer*4 i,il,istp,ly,lf,ldeg,ldeglw,ldegup,ldegf,ierr
      real*8 f,ksp2,dll1,omi,expo,slwgrn,ray,raymin
      real*8 fl,depst1,depst2,dys2,dxs
      complex*16 ca,cb,cag,cs1,cs2,cs3,cs4,ct1,ct2,cll1
      complex*16 ys1(4,0:2),ys2(4,0:2),yt1(4,0:2)
      complex*16 ypsv(6,4),ypsvg(6,4),ysh(2,2)
c
      complex*16 c2,c3,c4
      data c2,c3,c4/(2.d0,0.d0),(3.d0,0.d0),(4.d0,0.d0)/
c
      if(ig.eq.igfirst)then
        allocate(ksmallpsv(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ksmallpsv not allocated!'
        allocate(ksmallsh(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ksmallsh not allocated!'
c
        allocate(ldegpsv(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ldegpsv not allocated!'
        allocate(ldegsh(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ldegsh not allocated!'
c
        allocate(kp(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: kp not allocated!'
        allocate(ks(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: ks not allocated!'
        allocate(kt(ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: kt not allocated!'
        allocate(cps(6,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: cps not allocated!'
        allocate(cpt(2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: cpt not allocated!'
c
        allocate(mat2x2up(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2up not allocated!'
        allocate(mat2x2lw(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2lw not allocated!'
        allocate(mat2x2inv(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mat2x2inv not allocated!'
        allocate(mas2x2up(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas2x2up not allocated!'
        allocate(mas2x2lw(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas2x2lw not allocated!'
        allocate(mas2x2inv(2,2,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas2x2inv not allocated!'
        allocate(mas4x4up(4,4,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4up not allocated!'
        allocate(mas4x4lw(4,4,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4lw not allocated!'
        allocate(mas4x4inv(4,4,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas4x4inv not allocated!'
        allocate(mas6x6up(6,6,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6up not allocated!'
        allocate(mas6x6lw(6,6,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6lw not allocated!'
        allocate(mas6x6inv(6,6,ly0),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: mas6x6inv not allocated!'
c
        allocate(disk(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: disk not allocated!'
        allocate(fdisk(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: fdisk not allocated!'
c
        allocate(lylwsh(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: lylwsh not allocated!'
        allocate(lylwpsv(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: lylwpsv not allocated!'
        allocate(lyupatm(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: lylwatm not allocated!'
c
        allocate(zjup(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjup not allocated!'
        allocate(zjlw(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zjlw not allocated!'
        allocate(zhup(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhup not allocated!'
        allocate(zhlw(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: zhlw not allocated!'
        allocate(wj(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: wj not allocated!'
        allocate(wh(0:ldegmax,ly0,3),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: wh not allocated!'
c
        if(.not.nogravity)then
          allocate(zjupg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: zjupg not allocated!'
          allocate(zjlwg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: zjlwg not allocated!'
          allocate(zhupg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: zhupg not allocated!'
          allocate(zhlwg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: zhlwg not allocated!'
          allocate(wjg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: wjg not allocated!'
          allocate(whg(0:ldegmax),stat=ierr)
          if(ierr.ne.0)stop 'Error in qpgrnspec: whg not allocated!'
        endif
c
        allocate(dplm1(0:ldegmax),stat=ierr)
        if(ierr.ne.0)stop 'Error in qpgrnspec: dplm1 not allocated!'
      endif
c
c     tappered disk source:
c     2*[cos(theta)-cos(alpha)]/[1-cos(alpha)]^2
c
      if(rr0(ig).le.0.d0)then
        do ldeg=0,ldegmax
          fdisk(ldeg)=1.d0
        enddo
      else
        dxs=dcos(rr0(ig)/rrup(lys))
        dys2=(dsin(0.5d0*rr0(ig)/rrup(lys)))**2
        call legendre(ldegmax,rr0(ig)/rrup(lys),dplm1)
        fdisk(0)=1.d0
        fdisk(1)=(2.d0+dxs)/3.d0
        do ldeg=2,ldegmax
          fdisk(ldeg)=0.5d0/dys2**2
     &               /dble(ldeg)/dble(ldeg+1)/dble(ldeg+2)
     &               *(2.d0*dys2*(2.d0*dble(ldeg)*dys2+1.d0)
     &                +dble(ldeg)*dxs**2*dplm1(ldeg)
     &                -dble(2*ldeg+1)*dxs*dplm1(ldeg-1)
     &                +dble(ldeg+1)*dplm1(ldeg-2))
        enddo
      endif
      do ldeg=0,ldegmax
        disk(ldeg)=fdisk(ldeg)*dble(2*ldeg+1)/(4.d0*PI*rrup(lys)**2)
      enddo
c
      do ldeg=0,ldegmax
        fdisk(ldeg)=dabs(fdisk(ldeg))
        do il=ldeg+1,ldegmax
          fdisk(ldeg)=dmax1(fdisk(ldeg),dabs(fdisk(il)))
        enddo
      enddo
c
c     Initiation
c
      do istp=1,4
        do i=1,6
          ypsv(i,istp)=(0.d0,0.d0)
        enddo
      enddo
      do istp=1,2
        do i=1,2
          ysh(i,istp)=(0.d0,0.d0)
        enddo
      enddo
c
      open(21,file=rgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(22,file=tgrnfile(ig),
     &     form='unformatted',status='unknown')
      open(23,file=pgrnfile(ig),
     &     form='unformatted',status='unknown')
c
      if(fullwave)then
        do ldeg=1+ndmax,ldegmin-1
          dll1=dble(ldeg)*dble(ldeg+1)
          expo=0.d0
          do ly=min0(lys,lyr),max0(lys,lyr)-1
            expo=expo+dsqrt(dll1)*dlog(rrup(ly)/rrlw(ly))
          enddo
          if(expo.gt.expol+dlog(fdisk(ldeg)))goto 10
        enddo
10      ldeglw=ldeg
c
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
        slwgrn=dmin1(slwmax,slwfac*raymin/rrup(lys))
      else
        ldeglw=ldegmin
        slwgrn=slwmax
      endif
      ldegup=min0(ldegmax-1,
     &       ldeglw+idint(rrup(min0(lyr,lys))*PI2*fcut*slwgrn))
c
      write(*,*)' '
      write(*,'(a,i3,a,f7.2,a)')' ... calculate Green functions for ',
     &        ig,'. source at depth ',(grndep(ig)-depatmos)/KM2M,' km'
	write(*,'(a,i5)')'   max. harmonic degree: L_max = ',ldegup
c
      write(21)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(22)nt,ntcut,dt,nf,nfcut,df,ldegup
      write(23)nt,ntcut,dt,nf,nfcut,df,ldegup
c
      do istp=1,4
        do ldeg=0,1
          yr(ldeg,istp,0)=(0.d0,0.d0)
          yt(ldeg,istp,0)=(0.d0,0.d0)
          yp(ldeg,istp,0)=(0.d0,0.d0)
        enddo
      enddo
c
      do lf=1,nfcut
        f=dble(lf-1)*df
        omi=PI2*f
        ldegf=min0(ldegmax-1,ldeglw+idnint(dble(ldegup-ldeglw)*f/fcut))
c
        call qpqmodel(f)
        do ldeg=0,ldegf+1
          dll1=dble(ldeg)*dble(ldeg+1)
c
c         determine degree dependent starting layer number
c         of sh solution
c
          lylwsh(ldeg)=min0(lycm,ly0)
          expo=0.d0
          do ly=max0(lys,lyr,lyob)+1,min0(lycm,ly0)-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwsh(ldeg)=ly
              goto 100
            endif
          enddo
100       continue
c
c         determine degree dependent starting layer number
c         of psv solution
c
          lylwpsv(ldeg)=ly0
          expo=0.d0
          do ly=max0(lys,lyr,lyob)+1,min0(lycm,ly0)-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
          do ly=lycm,min0(lycc,ly0)-1
            ksp2=(omi*rrup(ly)/vpup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
          do ly=lycc,ly0-1
            ksp2=(omi*rrup(ly)/vsup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expos)then
              lylwpsv(ldeg)=ly
              goto 200
            endif
          enddo
200       continue
c
          lyupatm(ldeg)=1
          expo=0.d0
          do ly=lyob-1,1,-1
            ksp2=(omi*rrup(ly)/vpup(ly))**2
            if(dll1.gt.ksp2)then
              expo=expo+dsqrt(dll1-ksp2)*dlog(rrup(ly)/rrlw(ly))
            endif
            if(expo.gt.expoa)then
              lyupatm(ldeg)=ly
              goto 300
            endif
          enddo
300       continue
        enddo
c
        if(selpsv.or.lys.lt.lyob)then
          depst1=(rratmos-depatmos-rrup(lylwpsv(0)))/KM2M
          depst2=(rratmos-depatmos-rrup(lylwpsv(ldegf+1)))/KM2M
        else
          depst1=(rratmos-depatmos-rrup(lylwsh(1)))/KM2M
          depst2=(rratmos-depatmos-rrup(lylwsh(ldegf+1)))/KM2M
        endif
c
c       determine layer dependent max. harmonic degree
c       of sh solution
c
        do ly=lyob,min0(lycm-1,ly0)
          ldegsh(ly)=1
          do ldeg=1,ldegf+1
            if(lylwsh(ldeg).ge.ly)then
              ldegsh(ly)=ldeg
            endif
          enddo
        enddo
c
c       determine layer dependent max. harmonic degree
c       of psv solution
c
        do ly=1,lyob-1
          ldegpsv(ly)=0
          do ldeg=0,ldegf+1
            if(lyupatm(ldeg).le.ly)then
              ldegpsv(ly)=ldeg
            endif
          enddo
        enddo
c
        do ly=lyob,ly0
          ldegpsv(ly)=0
          do ldeg=0,ldegf+1
            if(lylwpsv(ldeg).ge.ly)then
              ldegpsv(ly)=ldeg
            endif
          enddo
        enddo
c
        do ly=1,ly0
          ksmallpsv(ly)=.true.
          ksmallsh(ly)=.true.
        enddo
c
        do il=0,2
          do istp=1,4
            ys1(istp,il)=(0.d0,0.d0)
            ys2(istp,il)=(0.d0,0.d0)
            yt1(istp,il)=(0.d0,0.d0)
          enddo
        enddo
c
        do ldeg=0,ldegf+1
          do istp=1,4
            ys1(istp,0)=ys1(istp,1)
            ys2(istp,0)=ys2(istp,1)
c
            ys1(istp,1)=ys1(istp,2)
            ys2(istp,1)=ys2(istp,2)
c
            yt1(istp,0)=yt1(istp,1)
c
            yt1(istp,1)=yt1(istp,2)
          enddo
c
          if(.not.selpsv)then
            do istp=1,4
	        do i=1,6
                ypsv(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else if(nogravity)then
            call qppsvkern(f,ldeg,ypsv)
          else
            fl=fgr*dble(ldeg)/dble(ldeggr)
            if(f.le.fgr-fl)then
              call qppsvkerng(f,ldeg,ypsv)
            else if(f.ge.fgr*(1.d0+FLTAPER)-fl)then
              call qppsvkern(f,ldeg,ypsv)
            else
              call qppsvkerng(f,ldeg,ypsvg)
              call qppsvkern(f,ldeg,ypsv)
              ca=dcmplx(dsin(((f+fl)/fgr-1.d0)/FLTAPER)**2,0.d0)
              cb=(1.d0,0.d0)-ca
              do istp=1,4
                do i=1,6
                  ypsv(i,istp)=ca*ypsv(i,istp)+cb*ypsvg(i,istp)
                enddo
              enddo
            endif
          endif
c
          if(.not.selsh.or.lys.lt.lyob)then
            do istp=1,2
	        do i=1,2
                ysh(i,istp)=(0.d0,0.d0)
              enddo
            enddo
          else
            call qpshkern(f,ldeg,ysh)
          endif
c
c         1. Explosion (M11=M22=M33=1)
c
          cs1=dcmplx( disk(ldeg)/(roup(lys)*vpup(lys)**2),0.d0)
          cs2=dcmplx(-disk(ldeg)*4.d0*(vsup(lys)/vpup(lys))**2
     &                /rrup(lys),0.d0)
          cs4=dcmplx( disk(ldeg)*2.d0*(vsup(lys)/vpup(lys))**2
     &                /rrup(lys),0.d0)
          ys1(1,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
          ys2(1,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          yt1(1,2)=(0.d0,0.d0)
c
c         2. Strike-slip (M12=M21=1)
c
          if(ldeg.lt.2.or.lys.lt.lyob)then
            ys1(2,2)=(0.d0,0.d0)
            ys2(2,2)=(0.d0,0.d0)
            yt1(2,2)=(0.d0,0.d0)
          else
            ct2=dcmplx(disk(ldeg)
     &                /(dble(ldeg)*dble(ldeg+1)*rrup(lys)),0.d0)
            cs4=-ct2
            ys1(2,2)=cs4*ypsv(1,4)
            ys2(2,2)=cs4*ypsv(3,4)
            yt1(2,2)=ct2*ysh(1,2)
          endif
c
c         3. Dip-slip (M13=M31=1)
c
          if(ldeg.lt.1.or.lys.lt.lyob)then
            ys1(3,2)=(0.d0,0.d0)
            ys2(3,2)=(0.d0,0.d0)
            yt1(3,2)=(0.d0,0.d0)
          else
            ct1=dcmplx(disk(ldeg)/(dble(ldeg)*dble(ldeg+1)
     &                *roup(lys)*vsup(lys)**2),0.d0)
            cs3=ct1
            ys1(3,2)=cs3*ypsv(1,3)
            ys2(3,2)=cs3*ypsv(3,3)
            yt1(3,2)=ct1*ysh(1,1)
          endif
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c
          if(lys.lt.lyob)then
            ys1(4,2)=(0.d0,0.d0)
            ys2(4,2)=(0.d0,0.d0)
          else
            cs1=dcmplx( disk(ldeg)/(roup(lys)*vpup(lys)**2),0.d0)
            cs2=dcmplx( disk(ldeg)*(3.d0-4.d0*(vsup(lys)/vpup(lys))**2)
     &                  /rrup(lys),0.d0)
            cs4=-(0.5d0,0.d0)*cs2
            ys1(4,2)=cs1*ypsv(1,1)+cs2*ypsv(1,2)+cs4*ypsv(1,4)
            ys2(4,2)=cs1*ypsv(3,1)+cs2*ypsv(3,2)+cs4*ypsv(3,4)
          endif
          yt1(4,2)=(0.d0,0.d0)
c
c         convert to (r,t,p,g)-system
c         ===========================
c
c         1. Explosion (M11=M22=M33=1)
c            yr normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,1,cos(t)) = -dPlm/dt
c
          yr(ldeg,1,0)=ys1(1,2)
          yt(ldeg,1,0)=-ys2(1,2)
c
c         2. Strike-slip (M12=M21=1)
c            yr normalized by Plm(l,2,cos(t))*sin(2p)
c            yt normalized by Plm(l,2,cos(t))*sin(2p)/sin(t)
c            yp normalized by Plm(l,2,cos(t))*cos(2p)/sin(t)
c
          yr(ldeg,2,0)=ys1(2,2)
          if(ldeg.ge.3)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-3)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+2)/dble(2*ldeg+1),0.d0)
            yt(ldeg-1,2,0)= ca*ys2(2,0)-cb*ys2(2,2)-(2.d0,0.d0)*yt1(2,1)
            yp(ldeg-1,2,0)=-ca*yt1(2,0)+cb*yt1(2,2)+(2.d0,0.d0)*ys2(2,1)
          endif
c
c         3. Dip-slip (M13=M31=1)
c            yr normalized by Plm(l,1,cos(t))*cos(p)
c            yt normalized by Plm(l,1,cos(t))*cos(p)/sin(t)
c            yp normalized by Plm(l,1,cos(t))*sin(p)/sin(t)
c
          yr(ldeg,3,0)=ys1(3,2)
          if(ldeg.ge.2)then
            ca=dcmplx(dble(ldeg-2)*dble(ldeg-2)/dble(2*ldeg-3),0.d0)
            cb=dcmplx(dble(ldeg+1)*dble(ldeg+1)/dble(2*ldeg+1),0.d0)
            yt(ldeg-1,3,0)= ca*ys2(3,0)-cb*ys2(3,2)+yt1(3,1)
            yp(ldeg-1,3,0)=-ca*yt1(3,0)+cb*yt1(3,2)-ys2(3,1)
          endif
c
c         4. CLVD (M33=1,M11=M22=-0.5)
c            yr normalized by Plm(l,0,cos(t))
c            yt normalized by Plm(l,0,cos(t)) = -dPlm/dt
c
          yr(ldeg,4,0)=ys1(4,2)
          yt(ldeg,4,0)=-ys2(4,2)
        enddo
c
        write(*,'(i6,a,f12.4,a,i5,a,2(f7.2,a))')lf,'.',1.0d+03*f,
     &       ' mHz: cut-off degree = ',ldegf,
     &       ', start depth = ',depst1,' - ',depst2,' km'
c
        write(21)ldegf
        write(22)ldegf
        write(23)ldegf
        write(21)((yr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(22)((yt(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
        write(23)((yp(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
      enddo
      close(21)
      close(22)
      close(23)
c
      if(ig.eq.iglast)then
        deallocate(ksmallpsv,ksmallsh,kp,ks,cps,cpt,
     &    mat2x2up,mat2x2lw,mat2x2inv,mas2x2up,mas2x2lw,mas2x2inv,
     &    mas4x4up,mas4x4lw,mas4x4inv,mas6x6up,mas6x6lw,mas6x6inv,
     &    disk,fdisk,lylwsh,lylwpsv,lyupatm,
     &    zjup,zjlw,zhup,zhlw,wj,wh,dplm1)
        if(.not.nogravity)then
          deallocate(zjupg,zjlwg,zhupg,zhlwg,wjg,whg)
        endif
      endif
c
      return
      end
