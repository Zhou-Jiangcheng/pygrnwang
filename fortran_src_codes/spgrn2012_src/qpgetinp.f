      subroutine qpgetinp(unit)
      implicit none
      integer unit
c
      include 'qpglobal.h'
c
c     work space
c
      integer i,j,l,ir,ig,isg,is,is1,flen,iswap,nhypo
      double precision twindow,twinout,suppress,r1,r2,dr,dr1,dr2
      double precision dswap(11)
      character*128 fswap
      character*180 comments
c
c     uniform receiver depth
c     ======================
c
      call getdata(unit,comments)
      read(comments,*)dpr
      dpr=KM2M*dpr
c
c     time (frequency) sampling
c     =========================
c
      call getdata(unit,comments)
      read(comments,*)twindow,dt
      ntcut=1+idnint(twindow/dt)
      nt=2
100   nt=2*nt
      if(nt.lt.ntcut)goto 100
      nf=nt/2
      if(nf.gt.nfmax)then
        stop 'Error: nfmax defined in qpglobal.h too small!'
      endif
      df=1.d0/(dble(nt)*dt)
c
      call getdata(unit,comments)
      read(comments,*)fcut
      nfcut=min0(nf,1+idnint(fcut/df))
      fcut=dble(nfcut-1)*df
      call getdata(unit,comments)
      read(comments,*)slwmax
      slwmax=slwmax/KM2M
c
      call getdata(unit,comments)
      read(comments,*)suppress
      if(suppress.le.0.d0.or.suppress.ge.1.d0)then
        fi=0.d0
      else
        fi=dlog(suppress)*df/PI2
      endif
c
c     cutoffs of spectra
c     ==================
c
      call getdata(unit,comments)
      read(comments,*)fgr,ldeggr
      if(fgr.lt.0.d0)fgr=0.d0
      if(ldeggr.lt.0)ldeggr=0
      if(fgr.gt.0.d0.and.ldeggr.le.0.or.
     &   fgr.le.0.d0.and.ldeggr.gt.0)then
        stop ' Bad fgr and ldeggr combination!'
      endif
      nogravity=fgr*dble(ldeggr).le.0.d0
c
      call getdata(unit,comments)
      read(comments,*)i,j
      selpsv=i.eq.1
      selsh=j.eq.1
      if(.not.(selpsv.or.selsh))then
        stop ' Error: none of PSV and SH is selected!'
      endif
c
c     Green's function files
c     ======================
c
      call getdata(unit,comments)
      read(comments,*)spcgrndir
      call getdata(unit,comments)
      read(comments,*)ngrn
      if(ngrn.le.0)then
        stop ' bad number of source depths!'
      else if(ngrn.gt.ngrnmax)then
        stop ' number of source depths exceeds the maximum!'
      endif
c
      do ig=1,ngrn
        call getdata(unit,comments)
        read(comments,*)grndep(ig),rr0(ig),grnfile(ig),grnsel(ig)
        if(grnsel(ig).lt.0.or.grnsel(ig).gt.1)then
          stop ' bad Green function selection!'
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
            dswap(1)=grndep(i)
            fswap=grnfile(i)
            iswap=grnsel(i)
c
            grndep(i)=grndep(j)
            grnfile(i)=grnfile(j)
            grnsel(i)=grnsel(j)
c
            grndep(j)=dswap(1)
            grnfile(j)=fswap
            grnsel(j)=iswap
          endif
        enddo
      enddo
c
      do flen=128,1,-1
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
      call getdata(unit,comments)
      read(comments,*)stdgrndir
c
      do flen=128,1,-1
        if(stdgrndir(flen:flen).ne.' ')goto 300
      enddo
300   continue
      do ig=1,ngrn
        stdgrnfile(ig)=stdgrndir(1:flen)//grnfile(ig)
      enddo
c
      call getdata(unit,comments)
      read(comments,*)infofile
      infofile=stdgrndir(1:flen)//infofile
c
      call getdata(unit,comments)
      read(comments,*)twinout,dtout
      if(twinout.gt.twindow)then
        stop ' Error: too large time window for space-time domain!'
      endif
      if(dtout.le.0.d0)then
        stop ' Error: wrong time sampling for space-time domain!'
      else if(dtout.lt.dt)then
        stop ' Error: too small time interval for space-time domain!'
      endif
      ntout=1+idnint(twinout/dtout)
c
      call getdata(unit,comments)
      read(comments,*)tstart,vreduc
      if(vreduc.le.0.d0)then
        stop ' Error in qpgetinp: wrong time reduction parameter!'
      endif
      vreduc=vreduc*KM2M
c
      call getdata(unit,comments)
      read(comments,*)tau
      call getdata(unit,comments)
      read(comments,*)nlpf,fcorner
      call getdata(unit,comments)
      read(comments,*)r1,r2,dr1,dr2
      if(r1.lt.0.d0.or.r1.gt.r2.or.dr1.le.0.d0.or.dr1.gt.dr2)then
        stop ' Error: wrong epidistance sampling for Green functions!'
      else
        r1=r1*KM2M
        r2=r2*KM2M
        dr1=dr1*KM2M
        dr2=dr2*KM2M
        nr=1+idnint(2.d0*(r2-r1)/(dr1+dr2))
        if(nr.gt.nrmax)then
          stop ' Error in qpgetinp: nrmax defined too small!'
        endif
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
c         rounded to meter
c
          dism(i)=dble(idnint(dism(i)))
          disrad(i)=dism(i)/REARTH
        enddo
      endif
c
c     multilayered model parameters
c     =============================
c
      call getdata(unit,comments)
      read(comments,*)l,i
      if(l.ge.lymax-2)then
        stop ' Error: lymax defined too small!'
      endif
      if(i.eq.1)then
        dispersion=.true.
      else
        dispersion=.false.
      endif
c
      qsmin=10000.d0
      depatmos=0.d0
      rratmos=REARTH
      do i=1,l
        call getdata(unit,comments)
        read(comments,*)j,dp0(i),vp0(i),vs0(i),ro0(i),qp0(i),qs0(i)
c
c       input units:    -,km,  km/s, km/s, g/cm^3,-,-
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
            stop ' Error: bad layering of earth model!'
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
        stop ' Error: bad start depth!'
      else if(dp0(l).gt.rratmos)then
        stop ' Error: bad definition of earth radius!'
      else if(dp0(l).lt.rratmos)then
        l=l+1
        if(l.ge.lymax-2)stop ' Error: lymax defined too small!'
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
c     end of inputs
c     =============
c
      return
      end
