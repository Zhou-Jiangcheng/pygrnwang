      subroutine tpssub(nz,h,vps,zs,rearth,head,diff,
     &                  nt,distance,traveltime,takeoff,slowness)
      implicit none
c-------------------------------------------------------------------------------
c     Calculte first arrival time of seismic waves based on ray theory
c
c     Input:
c     nz = number of layers of crust-mantle structure
c     h(nz) = layer thickness in m
c     vps(nz) = p or s wave velocity in m
c     zs = source depth in m
c     rearth = earth radius in m
c     head = .true. for including head waves
c     diff = .true. for including diffracted waves
c     nt = number of travel times to be returned
c     distance(nt) = distance in m
c
c     Output:
c     traveltime(nt) = travel time in second
c     takeoff(nt) = takeoff angle in degree
c     slowness(nt) = slowness in s/m
c-------------------------------------------------------------------------------
      integer*4 nz,nt
      real*8 zs,rearth
      logical*2 head, diff
      real*8 h(nz),vps(nz)
      real*8 distance(nt),traveltime(nt),takeoff(nt),slowness(nt)
c
      integer*4 npmax,npmin
      parameter(npmax=20001,npmin=100)
      integer*4 i,k,kc,m,iz,ir,jr,irs,ip,np,nr,ierr,nokay
      real*8 pi,piplus,deg2rad
      real*8 rs,a,b,c,ray,raymax,raymin,ins0,ins1,ins2,insup,inslw,dins
      real*8 dh,swap,dswap,tswap,slw,sgrd,fgrd,tup,tlw,dup,dlw
      logical*2 ninety,test
c
      logical*2, allocatable:: okay(:)
      integer*4, allocatable:: i1(:),i2(:)
      real*8, allocatable:: r(:),vel(:),tkf(:),ddg(:),tps(:),disrad(:)
c
      if(nt.le.1)then
        print *,'error in tpssub: nt <= 1!'
        stop
      endif
c
      test=.false.
c
      pi=4.d0*datan(1.d0)
      piplus=pi*(1.d0+1.0d-12)
      deg2rad=pi/180.d0
c
      ierr=0
      allocate(okay(nt),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: okay not allocated!'
      allocate(disrad(nt),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: disrad not allocated!'
c
      nr=nz+2
      allocate(r(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: r not allocated!'
      allocate(vel(nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: v not allocated!'
c
      allocate(i1(2*nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: i1 not allocated!'
      allocate(i2(2*nr),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: i2 not allocated!'
c
      np=npmax+(npmin+2)*nr
      allocate(tkf(np),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: tkf not allocated!'
      allocate(ddg(np),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: ddg not allocated!'
      allocate(tps(np),stat=ierr)
      if(ierr.ne.0)stop ' Error in tpssub: tps not allocated!'
c
      nr=nz
      r(1)=rearth
      vel(1)=vps(1)
      do iz=2,nz
        r(iz)=r(iz-1)-h(iz-1)
        vel(iz)=vps(iz)
        if(r(iz).lt.0.d0)then
          nr=iz
          r(nr)=0.d0
          goto 50
        endif
      enddo
50    continue
      if(vel(nr).gt.0.d0)then
        nr=nr+1
        r(nr)=0.d0
        vel(nr)=vel(nr-1)
      endif
c
c     delete pseudo interfaces
c
      m=nr
      do ir=2,m-1
        if(vel(ir).le.0.d0)then
          nr=ir
          goto 60
        else if(vel(ir).eq.vel(ir-1))then
          do jr=ir,nr-1
            r(jr)=r(jr+1)
            vel(jr)=vel(jr+1)
          enddo
          nr=nr-1
        endif
      enddo
60    continue
c
c     add a pseudo interface at source
c
      if(zs.le.0.d0)then
        rs=rearth
        irs=1
      else
        rs=rearth-zs
c
        if(rs.le.r(nr))then
          print *,'error in tpssub: too deep source!'
          stop
        endif
c
        do ir=2,nr
          if(rs.ge.r(ir))then
            irs=ir
            if(rs.gt.r(ir))then
              do jr=nr,ir,-1
                r(jr+1)=r(jr)
                vel(jr+1)=vel(jr)
              enddo
              r(ir)=rs
              vel(ir)=vel(ir-1)
              nr=nr+1
            endif
            goto 100
          endif
        enddo
      endif
100   continue
      swap=0.d0
      do ir=1,irs-1
        swap=swap+(r(ir)-r(ir+1))/vel(ir)
      enddo
c
      do i=1,np
        tkf(i)=0.d0
        ddg(i)=0.d0
        tps(i)=0.d0
      enddo
c
      kc=1
      i1(kc)=1
      k=0
      dh=rearth/dble(npmax-1)
c
c     wave phase with takeoff >= 90 deg
c
      if(irs.eq.1)then
        k=k+1
        tkf(k)=0.5d0*pi
        ddg(k)=0.d0
        tps(k)=0.d0
c
        k=k+1
        tkf(k)=0.5d0*pi
        ddg(k)=(r(1)-r(2))/rearth
        tps(k)=ddg(k)*rearth/vel(1)
c
        i2(kc)=k
        raymin=r(irs)/vel(irs)
      else
        k=k+1
        tkf(k)=pi
        ddg(k)=0.d0
        tps(k)=0.d0
        do ir=1,irs-1
          tps(k)=tps(k)+(r(ir)-r(ir+1))/vel(ir)
        enddo
c
        np=max0(npmin,1+idint((r(1)-r(irs))/dh))
c
        raymin=r(irs)/vel(irs-1)
        ninety=.true.
        do ir=2,irs-1
          if(raymin.gt.r(ir)/vel(ir-1))then
            raymin=r(ir)/vel(ir-1)
            ninety=.false.
          endif
        enddo
c
        dins=dasin(dmin1(1.d0,raymin*vel(irs-1)/r(irs)))/dble(np)
c
        do ip=1,np
          ins0=dble(ip)*dins
          if(ip.lt.np)then
            ray=dsin(ins0)*r(irs)/vel(irs-1)
          else
            ray=raymin
          endif
c
          dswap=0.d0
          tswap=0.d0
          do ir=irs-1,1,-1
            swap=dmin1(1.d0,ray*vel(ir)/r(ir+1))
            ins1=dasin(swap)
            ins2=dasin(swap*r(ir+1)/r(ir))
            dswap=dswap+ins1-ins2
            tswap=tswap+(r(ir)*dcos(ins2)-r(ir+1)*dcos(ins1))/vel(ir)
          enddo
          k=k+1
          tkf(k)=pi-ins0
          ddg(k)=dswap
          tps(k)=tswap
        enddo
        i2(kc)=k
        if(diff.and.ninety.and.vel(irs).lt.vel(irs-1))then
c
c         add diffracted wave
c
          slw=r(irs)/vel(irs-1)
          k=k+1
          ddg(k)=piplus
          tps(k)=tps(k-1)+(ddg(k)-ddg(k-1))*slw
          tkf(k)=tkf(k-1)
          i2(kc)=k
        else if(ninety.and.vel(irs).gt.vel(irs-1))then
c
c         add pn wave
c
          ins0=dasin(dmin1(1.d0,vel(irs-1)/vel(irs)))
          ray=dsin(ins0)*r(irs)/vel(irs-1)
c
          if(head)then
            dswap=0.d0
            tswap=0.d0
            do ir=irs-1,1,-1
              swap=dmin1(1.d0,ray*vel(ir)/r(ir+1))
              ins1=dasin(swap)
              ins2=dasin(swap*r(ir+1)/r(ir))
              dswap=dswap+ins1-ins2
              tswap=tswap+(r(ir)*dcos(ins2)-r(ir+1)*dcos(ins1))/vel(ir)
            enddo
            kc=kc+1
            slw=r(irs)/vel(irs)
            i1(kc)=i2(kc-1)+1
            k=i1(kc)
            ddg(k)=dswap
            tps(k)=tswap
            tkf(k)=0.5d0*pi
            k=k+1
            ddg(k)=piplus
            tps(k)=tps(k-1)+(ddg(k)-ddg(k-1))*slw
            tkf(k)=tkf(k-1)
            i2(kc)=k
          endif
          raymin=dmin1(raymin,r(irs)/vel(irs))
        endif
      endif
c
c     wave phase with takeoff < 90 deg
c
      raymax=raymin
c
      do jr=irs+1,nr
        swap=r(jr-1)/vel(jr-1)
        if(raymax.ge.swap)then
          insup=dasin(dmin1(1.d0,swap*vel(irs)/r(irs)))
        else
          insup=dasin(dmin1(1.d0,raymax*vel(irs)/r(irs)))
        endif
        swap=r(jr)/vel(jr-1)
        if(raymax.gt.swap)then
          inslw=dasin(dmin1(1.d0,swap*vel(irs)/r(irs)))
        else
          goto 500
        endif
        raymax=swap
c
        kc=kc+1
        i1(kc)=i2(kc-1)+1
        k=i2(kc-1)
c
        np=max0(npmin,1+idint((r(jr-1)-r(jr))/dh))
        dins=(inslw-insup)/dble(np)
c
        do ip=0,np
          ins0=insup+dble(ip)*dins
          ray=dsin(ins0)*r(irs)/vel(irs)
c
c         upgoing part of the current wave phase
c
          dup=0.d0
          tup=0.d0
          do ir=irs-1,1,-1
            swap=dmin1(1.d0,ray*vel(ir)/r(ir+1))
            ins1=dasin(swap)
            ins2=dasin(swap*r(ir+1)/r(ir))
            dup=dup+ins1-ins2
            tup=tup+(r(ir)*dcos(ins2)-r(ir+1)*dcos(ins1))/vel(ir)
          enddo
c
c         downgoing part of the current wave phase
c
          dlw=0.d0
          tlw=0.d0
          do ir=irs+1,jr
            swap=dmin1(1.d0,ray*vel(ir-1)/r(ir-1))
            ins1=dasin(swap)
            if(ir.lt.jr)then
              ins2=dasin(dmin1(1.d0,swap*r(ir-1)/r(ir)))
              dlw=dlw+ins2-ins1
              tlw=tlw+(r(ir-1)*dcos(ins1)-r(ir)*dcos(ins2))/vel(ir-1)
            else
              dlw=dlw+0.5d0*pi-ins1
              tlw=tlw+r(ir-1)*dcos(ins1)/vel(ir-1)
            endif
          enddo
c
          k=k+1
          tkf(k)=ins0
          ddg(k)=2.d0*dlw+dup
          tps(k)=2.d0*tlw+tup
        enddo
c
        i2(kc)=k
c
        if(diff.and.vel(jr).lt.vel(jr-1))then
c
c         add diffracted wave
c
          slw=r(jr)/vel(jr-1)
          k=k+1
          ddg(k)=piplus
          tps(k)=tps(k-1)+(ddg(k)-ddg(k-1))*slw
          tkf(k)=tkf(k-1)
          i2(kc)=k
        else if(jr.lt.nr.and.vel(jr).gt.vel(jr-1))then
c
c         add pn wave
c
          ray=r(jr)/vel(jr)
          ins0=dasin(dmin1(1.d0,ray*vel(irs)/r(irs)))
          if(head)then
            dup=0.d0
            tup=0.d0
            do ir=irs-1,1,-1
              swap=dmin1(1.d0,ray*vel(ir)/r(ir+1))
              ins1=dasin(swap)
              ins2=dasin(swap*r(ir+1)/r(ir))
              dup=dup+ins1-ins2
              tup=tup+(r(ir)*dcos(ins2)-r(ir+1)*dcos(ins1))/vel(ir)
            enddo
            dlw=0.d0
            tlw=0.d0
            do ir=irs+1,jr
              ins1=dasin(dmin1(1.d0,ray*vel(ir-1)/r(ir-1)))
              ins2=dasin(dmin1(1.d0,ray*vel(ir-1)/r(ir)))
              dlw=dlw+ins2-ins1
              tlw=tlw+(r(ir-1)*dcos(ins1)-r(ir)*dcos(ins2))/vel(ir-1)
            enddo
c
            dswap=2.d0*dlw+dup
            tswap=2.d0*tlw+tup
            slw=r(jr)/vel(jr)
c
            kc=kc+1
            i1(kc)=i2(kc-1)+1
            k=i1(kc)
            ddg(k)=dswap
            tps(k)=tswap
            tkf(k)=ins0
            k=k+1
            ddg(k)=piplus
            tps(k)=tps(k-1)+(ddg(k)-ddg(k-1))*slw
            tkf(k)=tkf(k-1)
            i2(kc)=k
          endif
          raymax=ray
        endif
500     continue
      enddo
c
      if(test)then
        open(31,file='testing.dat',status='unknown')
        write(31,'(a)')' distance[deg] traveltime[s]'
     &               //' takeoff[deg]'
        do m=1,kc
          do k=i1(m),i2(m)
            write(31,'(2f14.6,f13.6)')ddg(k)/deg2rad,tps(k),
     &          tkf(k)/deg2rad
          enddo
          write(31,'(a)')'            '
        enddo
      endif
c
c     interpolate and output
c
      do i=1,nt
        disrad(i)=dmod(distance(i)/rearth,2.d0*pi)
        if(disrad(i).gt.pi)disrad(i)=2.d0*pi-disrad(i)
      enddo
c
      nokay=0
c
      do i=1,nt
        traveltime(i)=0.d0
        okay(i)=.false.
        do m=1,kc
          do k=i1(m),i2(m)-1
            if(disrad(i).ge.dmin1(ddg(k),ddg(k+1)).and.
     &         disrad(i).le.dmax1(ddg(k),ddg(k+1)))then
              if(ddg(k).eq.ddg(k+1))then
                a=1.d0
                b=0.d0
              else
                a=(ddg(k+1)-disrad(i))/(ddg(k+1)-ddg(k))
                b=1.d0-a
              endif
              tswap=a*tps(k)+b*tps(k+1)
              if(.not.okay(i).or.traveltime(i).gt.tswap)then
                traveltime(i)=tswap
                takeoff(i)=a*tkf(k)+b*tkf(k+1)
                okay(i)=.true.
              endif
            endif
          enddo
        enddo
        if(okay(i))nokay=nokay+1
      enddo
c
      if(nokay.ge.nt)goto 1000
c
c     add reflected waves
c
      if(irs.eq.1)then
        m=npmin+nr
        dins=0.5d0*pi/dble(m)
        do k=1,m
          tkf(k)=dble(k-1)*dins
          ddg(k)=0.d0
          tps(k)=0.d0
        enddo
      else
        m=i2(1)
        do k=1,m
          tkf(k)=pi-tkf(k)
        enddo
      endif
c
      do ir=irs+1,nr-1
        do k=1,m
          ray=dsin(tkf(k))*r(irs)/vel(irs)
          if(ray.ge.r(ir)/vel(ir-1))goto 600
          swap=ray*vel(ir-1)/r(ir-1)
          ins1=dasin(swap)
          ins2=dasin(swap*r(ir-1)/r(ir))
          dlw=ins2-ins1
          tlw=(r(ir-1)*dcos(ins1)-r(ir)*dcos(ins2))/vel(ir-1)
          ddg(k)=ddg(k)+2.d0*dlw
          tps(k)=tps(k)+2.d0*tlw
        enddo
600     m=k-1
        if(test)then
          do k=1,m
            write(31,'(2f14.6,f13.6)')ddg(k)/deg2rad,tps(k),
     &          tkf(k)/deg2rad
          enddo
          write(31,'(a)')'            '
        endif
        nokay=0
        do i=1,nt
          do k=1,m-1
            if(disrad(i).ge.dmin1(ddg(k),ddg(k+1)).and.
     &         disrad(i).le.dmax1(ddg(k),ddg(k+1)))then
              if(ddg(k).eq.ddg(k+1))then
                a=1.d0
                b=0.d0
              else
                a=(ddg(k+1)-disrad(i))/(ddg(k+1)-ddg(k))
                b=1.d0-a
              endif
              tswap=a*tps(k)+b*tps(k+1)
              if(.not.okay(i).or.traveltime(i).gt.tswap)then
                traveltime(i)=tswap
                takeoff(i)=a*tkf(k)+b*tkf(k+1)
                okay(i)=.true.
              endif
            endif
          enddo
          if(okay(i))nokay=nokay+1
        enddo
      enddo
c
1000   continue
c
c     convert unit
c
      do i=1,nt
        if(irs.gt.1.and.takeoff(i).ge.0.5d0*pi)then
          slowness(i)=dsin(takeoff(i))/vel(irs-1)
        else
          slowness(i)=dsin(takeoff(i))/vel(irs)
        endif
        takeoff(i)=takeoff(i)/deg2rad
      enddo
c
      deallocate(okay,disrad,i1,i2,r,vel,tkf,ddg,tps)
c
      if(test)then
        close(31)
        pause
      endif
c
      if(nokay.lt.nt)then
        print *,zs,nokay,nt
        stop 'error in tpssub: a problem occurred!'
      endif
c
      return
      end