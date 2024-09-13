      subroutine qpgetinp(inpunit)
      use qpalloc
      implicit none
      integer*4 inpunit
c
c     work space
c
      integer*4 i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo,ierr
      real*8 suppress,r1,r2,dr,dr1,dr2,swap
      character*80 fswap
c
c     uniform receiver depth
c     ======================
c
      call skipdoc(inpunit)
      read(inpunit,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call skipdoc(inpunit)
      read(inpunit,*)twindow,dt
      ntcut=1+idnint(twindow/dt)
      nt=2
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      nf=nt/2
      df=1.d0/(dble(nt)*dt)
c
      call skipdoc(inpunit)
      read(inpunit,*)fcut
      nfcut=min0(nf,1+idnint(fcut/df))
      fcut=dble(nfcut-1)*df
      call skipdoc(inpunit)
      read(inpunit,*)slwmax
      slwmax=slwmax/KM2M
      fullwave=slwmax.le.0.d0
c
      call skipdoc(inpunit)
      read(inpunit,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        fi=0.d0
      else
        fi=dlog(suppress)*df/PI2
      endif
c
c     cutoffs of spectra
c     ==================
c
      call skipdoc(inpunit)
      read(inpunit,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop 'Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      call skipdoc(inpunit)
      read(inpunit,*)i,j
      selpsv=i.eq.1
      selsh=j.eq.1
      if(.not.(selpsv.or.selsh))then
        stop 'Error: none of PSV and SH is selected!'
      endif
c
c     Green's function files
c     ======================
c
      call skipdoc(inpunit)
      read(inpunit,*)spcgrndir
      call skipdoc(inpunit)
      read(inpunit,*)ngrn
      if(ngrn.le.0)then
        stop 'bad number of source depths!'
      endif
c
      ierr=0
      allocate(rr0(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: rr0 not allocated!'
      allocate(lygrn(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: lygrn not allocated!'
      allocate(grnsel(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: grnsel not allocated!'
      allocate(grndep(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: grndep not allocated!'
      allocate(grnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: grnfile not allocated!'
      allocate(rgrnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: rgrnfile not allocated!'
      allocate(tgrnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: tgrnfile not allocated!'
      allocate(pgrnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: pgrnfile not allocated!'
      allocate(stdgrnfile(ngrn),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: stdgrnfile not allocated!'
c
      do ig=1,ngrn
        call skipdoc(inpunit)
        read(inpunit,*)grndep(ig),rr0(ig),grnfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop 'bad Green function selection!'
        endif
        grndep(ig)=grndep(ig)*KM2M
        rr0(ig)=rr0(ig)*KM2M
      enddo
c
c     sort green function files by source depth
c
      do i=1,ngrn
        do j=i+1,ngrn
          if(grndep(j).lt.grndep(i))then
            swap=grndep(i)
            fswap=grnfile(i)
            iswap=grnsel(i)
c
            grndep(i)=grndep(j)
            grnfile(i)=grnfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=swap
            grnfile(j)=fswap
            grnsel(j)=iswap
          endif
        enddo
      enddo
c
      do flen=80,1,-1
        if(spcgrndir(flen:flen).ne.' ')goto 200
      enddo
200   continue
      do ig=1,ngrn
        rgrnfile(ig)=spcgrndir(1:flen)//'R_'//grnfile(ig)
        tgrnfile(ig)=spcgrndir(1:flen)//'T_'//grnfile(ig)
        pgrnfile(ig)=spcgrndir(1:flen)//'P_'//grnfile(ig)
      enddo
c
c     receiver parameters
c     ===================
c
      call skipdoc(inpunit)
      read(inpunit,*)stdgrndir
c
      do flen=80,1,-1
        if(stdgrndir(flen:flen).ne.' ')goto 300
      enddo
300   continue
      do ig=1,ngrn
        stdgrnfile(ig)=stdgrndir(1:flen)//grnfile(ig)
      enddo
c
      call skipdoc(inpunit)
      read(inpunit,*)infofile,tptable,tstable
      infofile=stdgrndir(1:flen)//infofile
      tptable=stdgrndir(1:flen)//tptable
      tstable=stdgrndir(1:flen)//tstable
c
      call skipdoc(inpunit)
      read(inpunit,*)twinout,dtout
      if(twinout.gt.twindow)then
        stop 'Error: too large time window for space-time domain!'
      endif
      if(dtout.le.0.d0)then
        stop 'Error: wrong time sampling for space-time domain!'
      else if(dtout.lt.dt)then
        stop 'Error: too small time interval for space-time domain!'
      endif
      ntout=1+idnint(twinout/dtout)
c
      call skipdoc(inpunit)
      read(inpunit,*)tstart
c
      call skipdoc(inpunit)
      read(inpunit,*)tau
      call skipdoc(inpunit)
      read(inpunit,*)nbpf,f1corner,f2corner
      call skipdoc(inpunit)
      read(inpunit,*)r1,r2,dr1,dr2
      if(r1.lt.0.d0.or.r1.gt.r2.or.dr1.le.0.d0.or.dr1.gt.dr2)then
        stop 'Error: wrong epidistance sampling for Green functions!'
      endif
      r1=r1*KM2M
      r2=r2*KM2M
      dr1=dr1*KM2M
      dr2=dr2*KM2M
      nr=1+idnint(2.d0*(r2-r1)/(dr1+dr2))
c
      allocate(dism(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: dism not allocated!'
      allocate(disrad(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: disrad not allocated!'
c
      dism(1)=r1
      dism(nr)=r2
      if(nr.gt.2)then
        dr=(r2-r1)/dble(nr-1)
        dr1=dr1*2.d0*dr/(dr1+dr2)
        dr2=dr2*2.d0*dr/(dr1+dr2)
        do i=2,nr-1
          dr=dr1+(dr2-dr1)*dble(i-2)/dble(nr-2)
          dism(i)=dism(i-1)+dr
        enddo
      endif
      do i=1,nr
c
c       rounded to meter
c
        dism(i)=dble(idnint(dism(i)))
        disrad(i)=dism(i)/REARTH
      enddo
c
c     multilayered model parameters
c     =============================
c
      call skipdoc(inpunit)
      read(inpunit,*)l,i
      if(l.le.0)then
        stop 'Error: bad number of data lines of earth model!'
      endif
      dispersion=i.eq.1
c
      l0=l+1
c
      allocate(dp0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: dp0 not allocated!'
      allocate(vp0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vp0 not allocated!'
      allocate(vs0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vs0 not allocated!'
      allocate(ro0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: ro0 not allocated!'
      allocate(qp0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qp0 not allocated!'
      allocate(qs0(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qs0 not allocated!'
c
      allocate(dp0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: dp0up not allocated!'
      allocate(vp0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vp0up not allocated!'
      allocate(vs0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vs0up not allocated!'
      allocate(ro0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: ro0up not allocated!'
      allocate(qp0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qp0up not allocated!'
      allocate(qs0up(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qs0up not allocated!'
c
      allocate(dp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: dp0lw not allocated!'
      allocate(vp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vp0lw not allocated!'
      allocate(vs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: vs0lw not allocated!'
      allocate(ro0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: ro0lw not allocated!'
      allocate(qp0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qp0lw not allocated!'
      allocate(qs0lw(l0),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpgetinp: qs0lw not allocated!'
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=REARTH
      do i=1,l
        call skipdoc(inpunit)
        read(inpunit,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input inpunits:    -,km,  km/s, km/s, g/cm^3,-,-
c
        if(i.eq.1)then
          depatmos=-KM2M*dp0(1)
          rratmos=REARTH+depatmos
        endif
        dp0(i)=KM2M*dp0(i)+depatmos
        vp0(i)=KM2M*vp0(i)
        vs0(i)=KM2M*vs0(i)
        ro0(i)=KM2M*ro0(i)
        if(vs0(i).gt.0.d0)qsmin=dmin1(qsmin,qs0(i))
        if(i.gt.1)then
          if(dp0(i).lt.dp0(i-1))then
            stop 'Error: bad layering of earth model!'
          endif
        endif
      enddo
c
      dpr=dpr+depatmos
      do i=1,ngrn
        grndep(i)=grndep(i)+depatmos
      enddo
c
      if(dp0(1).ne.0.d0)then
        stop 'Error: bad start depth!'
      else if(dp0(l).gt.rratmos)then
        stop 'Error: bad definition of earth radius!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        dp0(l)=rratmos
        vp0(l)=vp0(l-1)
        vs0(l)=vs0(l-1)
        ro0(l)=ro0(l-1)
        qp0(l)=qp0(l-1)
        qs0(l)=qs0(l-1)
      endif
c
      l0=0
      do i=2,l
        if(dp0(i).gt.dp0(i-1))then
          l0=l0+1
          dp0up(l0)=dp0(i-1)
          vp0up(l0)=vp0(i-1)
          vs0up(l0)=vs0(i-1)
          ro0up(l0)=ro0(i-1)
          qp0up(l0)=qp0(i-1)
          qs0up(l0)=qs0(i-1)
c
          dp0lw(l0)=dp0(i)
          vp0lw(l0)=vp0(i)
          vs0lw(l0)=vs0(i)
          ro0lw(l0)=ro0(i)
          qp0lw(l0)=qp0(i)
          qs0lw(l0)=qs0(i)
        endif
      enddo
c
      deallocate(dp0,vp0,vs0,ro0,qp0,qs0)
c
c     end of inputs
c     =============
c
      return
      end
