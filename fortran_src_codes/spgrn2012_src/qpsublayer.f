	subroutine qpsublayer(ierr)
	implicit none
	integer ierr
c
      include 'qpglobal.h'
c
c	work space
c
	integer i,j,i0,l,ly,ig,di0
      double precision h,dh,z,zz,wvlen,up,lw,uplw4
	double precision rrs,rrr,dvp,dvs,ro1,dro,dqp,dqs,mass
c
	ly0=0
c
      zz=0.d0
	do l=1,l0
        if(vs0up(l).gt.0.d0)then
          wvlen=vs0up(l)/fcut
        else
          wvlen=vp0up(l)/fcut
        endif
	  h=dp0lw(l)-dp0up(l)
	  dvp=2.d0*dabs(vp0lw(l)-vp0up(l))/(vp0lw(l)+vp0up(l))
        if(vs0lw(l)+vs0up(l).gt.0.d0)then
	    dvs=2.d0*dabs(vs0lw(l)-vs0up(l))/(vs0lw(l)+vs0up(l))
        else
          dvs=0.d0
        endif
        dro=2.d0*dabs(dlog(ro0lw(l))-dlog(ro0up(l)))
     &     /(dlog(ro0lw(l))+dlog(ro0up(l)))
	  i0=1+idint(dmax1(dvp/RESOLUT,dvs/RESOLUT,dro/RESOLUT))
        i0=min0(i0,1+idint(4.d0*h/wvlen))
        if(fgr.gt.0.d0)then
          i0=max0(i0,1+idint(4.d0*h/(wvlen*fcut/fgr)))
        endif
        dvp=(vp0lw(l)-vp0up(l))/h
	  dvs=(vs0lw(l)-vs0up(l))/h
	  dro=(dlog(ro0lw(l))-dlog(ro0up(l)))/h
	  dqp=(qp0lw(l)-qp0up(l))/h
	  dqs=(qs0lw(l)-qs0up(l))/h
	  dh=h/dble(i0)
	  do i=1,i0
	    ly0=ly0+1
	    if(ly0.ge.lymax-2-ngrn)then
	      stop ' Error in qpsublayer: lymax too small!'
	    endif
	    z=dble(i-1)*dh
	    rrup(ly0)=rratmos-(zz+z)
	    vpup(ly0)=vp0up(l)+dvp*z
	    vsup(ly0)=vs0up(l)+dvs*z
	    roup(ly0)=ro0up(l)*dexp(dro*z)
	    qpup(ly0)=qp0up(l)+dqp*z
	    qsup(ly0)=qs0up(l)+dqs*z
          z=z+dh
	    rrlw(ly0)=rratmos-(zz+z)
	    vplw(ly0)=vp0up(l)+dvp*z
	    vslw(ly0)=vs0up(l)+dvs*z
	    rolw(ly0)=ro0up(l)*dexp(dro*z)
	    qplw(ly0)=qp0up(l)+dqp*z
	    qslw(ly0)=qs0up(l)+dqs*z
	  enddo
        zz=zz+h
	enddo
c
c     add source layers
c
      do ig=1,ngrn
        rrs=rratmos-grndep(ig)
        if(rrs.lt.0.d0)then
          stop ' Wrong source depth!'
        endif
        do ly=1,ly0
          if(rrs.ge.rrlw(ly))then
            lys=ly
            goto 100
          endif
        enddo
100     continue
        if(rrs.lt.rrup(lys).and.rrs.gt.rrlw(lys))then
          do ly=ly0,lys,-1
            rrup(ly+1)=rrup(ly)
	      vpup(ly+1)=vpup(ly)
	      vsup(ly+1)=vsup(ly)
	      roup(ly+1)=roup(ly)
	      qpup(ly+1)=qpup(ly)
	      qsup(ly+1)=qsup(ly)
            rrlw(ly+1)=rrlw(ly)
	      vplw(ly+1)=vplw(ly)
	      vslw(ly+1)=vslw(ly)
	      rolw(ly+1)=rolw(ly)
	      qplw(ly+1)=qplw(ly)
	      qslw(ly+1)=qslw(ly)
          enddo
          lys=lys+1
          up=(rrs-rrlw(lys))/(rrup(lys-1)-rrlw(lys))
          lw=1.d0-up
          rrlw(lys-1)=rrs
	    vplw(lys-1)=up*vpup(lys-1)+lw*vplw(lys)
	    vslw(lys-1)=up*vsup(lys-1)+lw*vslw(lys)
	    rolw(lys-1)=up*roup(lys-1)+lw*rolw(lys)
	    qplw(lys-1)=up*qpup(lys-1)+lw*qplw(lys)
	    qslw(lys-1)=up*qsup(lys-1)+lw*qslw(lys)
          rrup(lys)=rrs
	    vpup(lys)=vplw(lys-1)
	    vsup(lys)=vslw(lys-1)
	    roup(lys)=rolw(lys-1)
	    qpup(lys)=qplw(lys-1)
	    qsup(lys)=qslw(lys-1)
          ly0=ly0+1       
        endif
      enddo
c
c     add receiver layer
c
      rrr=rratmos-dpr
      if(rrr.lt.0.d0)then
        stop ' Wrong receiver depth!'
      endif
      do ly=1,ly0
        if(rrr.ge.rrlw(ly))then
          lyr=ly
          goto 200
        endif
      enddo
200   continue
      if(rrr.lt.rrup(lyr).and.rrr.gt.rrlw(lyr))then
        do ly=ly0,lyr,-1
          rrup(ly+1)=rrup(ly)
	    vpup(ly+1)=vpup(ly)
	    vsup(ly+1)=vsup(ly)
	    roup(ly+1)=roup(ly)
	    qpup(ly+1)=qpup(ly)
	    qsup(ly+1)=qsup(ly)
          rrlw(ly+1)=rrlw(ly)
	    vplw(ly+1)=vplw(ly)
	    vslw(ly+1)=vslw(ly)
	    rolw(ly+1)=rolw(ly)
	    qplw(ly+1)=qplw(ly)
	    qslw(ly+1)=qslw(ly)
        enddo
        lyr=lyr+1
        up=(rrr-rrlw(lyr))/(rrup(lyr-1)-rrlw(lyr))
        lw=1.d0-up
        rrlw(lyr-1)=rrr
	  vplw(lyr-1)=up*vpup(lyr-1)+lw*vplw(lyr)
	  vslw(lyr-1)=up*vsup(lyr-1)+lw*vslw(lyr)
	  rolw(lyr-1)=up*roup(lyr-1)+lw*rolw(lyr)
	  qplw(lyr-1)=up*qpup(lyr-1)+lw*qplw(lyr)
	  qslw(lyr-1)=up*qsup(lyr-1)+lw*qslw(lyr)
        rrup(lyr)=rrr
	  vpup(lyr)=vplw(lyr-1)
	  vsup(lyr)=vslw(lyr-1)
	  roup(lyr)=rolw(lyr-1)
	  qpup(lyr)=qplw(lyr-1)
	  qsup(lyr)=qslw(lyr-1)
        ly0=ly0+1      
      endif
c
c     determine indices of main interfaces
c
      lyos=1
      lyob=1
      if(vsup(1).le.0.d0)then
        do ly=2,ly0
          if(vsup(ly).gt.0.d0.and.vslw(ly-1).le.0.d0)then
            lyob=ly
            goto 301
          endif
        enddo
301     continue
        do ly=2,lyob-1
          if(roup(ly)-rolw(ly-1).gt.rolw(ly-1))then
            lyos=ly
            goto 302
          endif
        enddo
302     continue
      endif
      lycm=ly0+1
      do ly=max0(2,lyob+1),ly0
        if(vsup(ly).le.0.d0.and.vslw(ly-1).gt.0.d0)then
          lycm=ly
          goto 400
        endif
      enddo
400   lycc=ly0+1
      do ly=max0(2,lycm+1),ly0
        if(vsup(ly).gt.0.d0.and.vslw(ly-1).le.0.d0)then
          lycc=ly
          goto 500
        endif
      enddo
500   continue
c
c     determine indices of receiver layer
c
      rrr=rratmos-dpr
      do ly=1,ly0
        if(rrr.ge.rrup(ly))then
          lyr=ly
          goto 600
        endif
      enddo
600   continue
c
c     determine indices of source layers
c
      do ig=1,ngrn
        rrs=rratmos-grndep(ig)
        do ly=1,ly0
          if(rrs.ge.rrup(ly))then
            lygrn(ig)=ly
            goto 700
          endif
        enddo
700     continue
c        if(lygrn(ig).le.lyob.or.lygrn(ig).ge.lycm)then
c          stop ' Source in a liquid layer!'
c        endif
      enddo
c
      mass=0.d0
      do ly=ly0,1,-1
        dro=(roup(ly)-rolw(ly))/(rrup(ly)-rrlw(ly))
        ro1=rolw(ly)-dro*rrlw(ly) 
        mass=mass+PI*(rrup(ly)-rrlw(ly))*((4.d0/3.d0)*ro1
     &      *(rrup(ly)**2+rrup(ly)*rrlw(ly)+rrlw(ly)**2)
     &      +dro*(rrup(ly)**3+rrup(ly)**2*rrlw(ly)
     &      +rrup(ly)*rrlw(ly)**2+rrlw(ly)**3))
        cgrup(ly)=dcmplx(BIGG*mass/rrup(ly)**2,0.d0)
      enddo
      freeairgrd=-dreal(cgrup(lyr))*2.d0/rrup(lyr)
      if(lyr.gt.1)then
        freeairgrd=freeairgrd+4.d0*PI*BIGG*rolw(lyr-1)
      endif
c
      do ly=1,ly0-1
        cgrlw(ly)=cgrup(ly+1)
      enddo
      cgrlw(ly0)=(0.d0,0.d0)
c
	write(*,'(9a)')'  No','      R(km)','   Vp(km/s)',
     &    '   Vs(km/s)',' Ro(g/cm^3)','      Qp',
     &    '      Qs',' g(m/s^2)'
	do ly=1,ly0
	  write(*,1001)ly,rrup(ly)/1.d3,vpup(ly)/1.d3,vsup(ly)/1.d3,
     &               roup(ly)/1.d3,qpup(ly),qsup(ly),dreal(cgrup(ly))
        j=0
        do ig=1,ngrn
          if(lygrn(ig).eq.ly)then
            j=j+1
            write(*,'(a3,$)')' S '
          endif
        enddo
        if(j.eq.0)write(*,'(a3,$)')'   '
        if(lyr.eq.ly)then
          write(*,'(a3)')' R '
        else
          write(*,'(a3)')'   '
        endif
	  write(*,1002)rrlw(ly)/1.d3,vplw(ly)/1.d3,vslw(ly)/1.d3,
     &             rolw(ly)/1.d3,qplw(ly),qslw(ly),dreal(cgrlw(ly))
	enddo
c
      do ly=1,ly0
        crrup(ly)=dcmplx(rrup(ly),0.d0)
        croup(ly)=dcmplx(roup(ly),0.d0)
        claup(ly)=dcmplx(roup(ly)*(vpup(ly)**2-2.d0*vsup(ly)**2),0.d0)
        cmuup(ly)=dcmplx(roup(ly)*vsup(ly)**2,0.d0)
        cvpup(ly)=dcmplx(vpup(ly),0.d0)
        cvsup(ly)=dcmplx(vsup(ly),0.d0)
        cgaup(ly)=dcmplx(2.d0*PI2*BIGG*roup(ly),0.d0)
c
        crrlw(ly)=dcmplx(rrlw(ly),0.d0)
        crolw(ly)=dcmplx(rolw(ly),0.d0)
        clalw(ly)=dcmplx(rolw(ly)*(vplw(ly)**2-2.d0*vslw(ly)**2),0.d0)
        cvplw(ly)=dcmplx(vplw(ly),0.d0)
        cvslw(ly)=dcmplx(vslw(ly),0.d0)
        cmulw(ly)=dcmplx(rolw(ly)*vslw(ly)**2,0.d0)
        cgalw(ly)=dcmplx(2.d0*PI2*BIGG*rolw(ly),0.d0)
c
        cro(ly)=dcmplx(dsqrt(roup(ly)*rolw(ly)),0.d0)
        cla(ly)=(0.5d0,0.d0)*(claup(ly)+clalw(ly))
        cmu(ly)=(0.5d0,0.d0)*(cmuup(ly)+cmulw(ly))
        cvp(ly)=(0.5d0,0.d0)*(cvpup(ly)+cvplw(ly))
        cvs(ly)=(0.5d0,0.d0)*(cvsup(ly)+cvslw(ly))
        cgr(ly)=(0.5d0,0.d0)*(cgrup(ly)+cgrlw(ly))
        cga(ly)=(0.5d0,0.d0)*(cgaup(ly)+cgalw(ly))
      enddo
c
      mmantle=0.d0
      do ly=lyob,lycm-1
        ro1=0.5d0*(roup(ly)+rolw(ly))
        uplw4=rrup(ly)**4+rrup(ly)**3*rrlw(ly)
     &       +rrup(ly)**2*rrlw(ly)**2+rrup(ly)*rrlw(ly)**3
     &       +rrlw(ly)**4
        mshell(ly)=ro1*(rrup(ly)-rrlw(ly))*uplw4*8.d0*PI/15.d0
        mmantle=mmantle+mshell(ly)
      enddo
c
      do ly=1,lyob-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)
        cypnorm(3,ly)=cga(ly)
        cypnorm(4,ly)=cga(ly)
      enddo
      do ly=lyob,lycm-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(5,ly)=cga(ly)
        cypnorm(6,ly)=cga(ly)
      enddo
      do ly=lycm,lycc-1
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)
        cypnorm(3,ly)=cga(ly)
        cypnorm(4,ly)=cga(ly)
      enddo
      do ly=lycc,ly0
        cypnorm(1,ly)=(1.d0,0.d0)
        cypnorm(2,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(3,ly)=(1.d0,0.d0)
        cypnorm(4,ly)=cla(ly)+(2.d0,0.d0)*cmu(ly)
        cypnorm(5,ly)=cga(ly)
        cypnorm(6,ly)=cga(ly)
      enddo
c
1001	format(i4,f11.4,3f11.4,2f8.1,f9.4,$)
1002	format(f15.4,3f11.4,2f8.1,f9.4)
      return
	end
