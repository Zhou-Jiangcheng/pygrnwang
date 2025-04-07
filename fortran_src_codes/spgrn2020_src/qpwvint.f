      subroutine qpwvint(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,m,id,ir,it,ig,nd,nt0,nf0,ntcut0,nfcut0
      integer*4 lf,istp,ldeg,ldegf,ldegup,ldegcri
      real*4 tredr4
      real*8 depsarc,f,dt0,df0,rn,re,azi,t,tp2s
      real*8 a,b,slwcut,discri,twinred
      complex*16 cp0,cp1,cp2
      complex*16 cfac,ca,cb
c
      integer*4, allocatable:: idr(:),ldegtap(:,:)
      real*8, allocatable:: tap(:,:)
      complex*16, allocatable:: ssd(:),ssf(:),wvf(:),bpf(:)
c
      integer*4 fseek,ftell
c
      twinred=0.65d0*(dble(nt-1)*dt-tstart)
c
      ierr=0
      allocate(idr(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: idr not allocated!'
      allocate(ldegtap(4,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: ldegtap not allocated!'
      allocate(ssd(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: ssd not allocated!'
      allocate(ssf(nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: ssf not allocated!'
      allocate(wvf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: wvf not allocated!'
      allocate(bpf(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: ssd not allocated!'
      allocate(tap(0:ldegmax,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: tap not allocated!'
      allocate(plm(0:ldegmax,0:2,nr),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: plm not allocated!'
      allocate(sgrn(nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: sgrn not allocated!'
      allocate(tgrn(2*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: tgrn not allocated!'
      allocate(vgrnr(nf,nr,4),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: vgrnr not allocated!'
      allocate(vgrnt(nf,nr,4),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: vgrnt not allocated!'
      allocate(vgrnp(nf,nr,4),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: vgrnp not allocated!'
      allocate(cswap(2*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: cswap not allocated!'
      allocate(dswap(4*nf),stat=ierr)
      if(ierr.ne.0)stop 'Error in qpwvint: dswap not allocated!'
c
      call swavelet(tau,df,fi,nf,wvf)
c
c     band-pass filter
c
      if(nbpf.gt.0.and.f2corner.gt.f1corner.and.f1corner.gt.0.d0)then
        if(f1corner.le.0.d0)then
          call butterworth(nbpf,f2corner,df,nf,cswap)
        else
          call bandpass(nbpf,f1corner,f2corner,df,nf,cswap)
        endif
c
        m=1
        do lf=2*nf,nf+2,-1
          m=m+1
          cswap(lf)=dconjg(cswap(m))
        enddo
        cswap(nf+1)=(0.d0,0.d0)
        call four1w(cswap,dswap,2*nf,+1)
        do lf=1,2*nf
          t=dble(lf-1)*dt
          cswap(lf)=dcmplx(dswap(2*lf-1)*dexp(PI2*fi*t)*df,0.d0)
        enddo
        call four1w(cswap,dswap,2*nf,-1)
        do lf=1,nf
          bpf(lf)=cswap(lf)*dcmplx(dt,0.d0)
        enddo
        do lf=1,nf
          wvf(lf)=wvf(lf)*bpf(lf)
        enddo
      endif
c
      do ir=1,nr
        ssd(ir)=dcmplx(dsin(disrad(ir)),0.d0)
        ssf(ir)=dcmplx(2.d0*dsin(0.5d0*disrad(ir))**2,0.d0)
      enddo
c
      call qplegendre(ierr)
c
      do ig=1,ngrn
        lys=lygrn(ig)
        depsarc=grndep(ig)/REARTH
        do ir=1,nr
          if(slwmax.le.0.d0)then
            idr(ir)=0
          else
            idr(ir)=min0(ndmax,max0(0,
     &            idnint(dlog(disrad(ir)/depsarc)/dlog(50.d0))))
          endif
        enddo
c
        write(*,'(a)')' '
        write(*,'(a)')' open Green function data base: '
     &              //grnfile(ig)(1:40)
        write(*,'(a)')' ... please wait ...'
c
        open(21,file=rgrnfile(ig),
     &       form='unformatted',status='old')
        open(22,file=tgrnfile(ig),
     &       form='unformatted',status='old')
        open(23,file=pgrnfile(ig),
     &       form='unformatted',status='old')
c
        open(30,file=stdgrnfile(ig),
     &       form='unformatted',status='unknown')
c
        read(21)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(22)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
        read(23)nt0,ntcut0,dt0,nf0,nfcut0,df0,ldegup
c
        if(ntcut0.ne.ntcut.or.dt0.ne.dt.or.
     &     nfcut0.ne.nfcut.or.df0.ne.df)then
          print *,' Error in qpwvint: t/f sampling'
     &          //' inconsistent with Green functions!'
          write(*,'(a)')'               '//'  ntcut             dt'
     &                                   //'  nfcut             df'
          write(*,'(a,2(i7,f16.8))')' Current input:',
     &                               ntcut,dt,nfcut,df
          write(*,'(a,2(i7,f16.8))')'     Data base:',
     &                               ntcut0,dt0,nfcut0,df0
          stop
        endif
c
        nd=0
        do ir=1,nr
          nd=max0(nd,idr(ir))
        enddo
c
        do istp=1,4
          do ir=1,nr
            do lf=1,nf
              vgrnr(lf,ir,istp)=(0.d0,0.d0)
              vgrnt(lf,ir,istp)=(0.d0,0.d0)
            enddo
          enddo
        enddo
        do istp=2,3
          do ir=1,nr
            do lf=1,nf
              vgrnp(lf,ir,istp)=(0.d0,0.d0)
            enddo
          enddo
        enddo
c
        do lf=1,nfcut
          f=dble(lf-1)*df
          read(21)ldegf
          read(22)ldegf
          read(23)ldegf
c
          read(21)((yr(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(22)((yt(ldeg,istp,0),ldeg=0,ldegf),istp=1,4)
          read(23)((yp(ldeg,istp,0),ldeg=0,ldegf),istp=2,3)
c
          do ir=1,nr
            ldegtap(1,ir)=0
            ldegtap(2,ir)=0
          enddo
c
          m=0
          do ir=1,nr
            tp2s=ts(ir,ig)-tp(ir,ig)
            if(tp2s.lt.twinred)then
              slwcut=0.d0
              ldegtap(4,ir)=ldegf-nd
              m=ir
            else
              slwcut=slwtp(ir,ig)
     &              +(slwts(ir,ig)-slwtp(ir,ig))*twinred/tp2s
              ldegtap(4,ir)=max0(10*nd,
     &                       nd+idint(rrup(min0(lyr,lys))*PI2*f*slwcut))
            endif
          enddo
c
          if(m.lt.nr)then
            do ir=1,m
              b=dism(ir)/dism(m+1)
              a=1.d0-b
              ldegtap(4,ir)=idnint(a*dble(ldegf-nd)
     &                            +b*dble(ldegtap(4,m+1)))
            enddo
          endif
c
          do ir=1,nr
            ldegtap(3,ir)=ldegtap(4,ir)*4/5
            call taper(ldegtap(1,ir),ldegf,tap(0,ir))
          enddo
c
c         use differential filter to suppress spatial aliasing
c
          do id=0,nd
            yt(0,1,id)=(0.d0,0.d0)
            yt(0,4,id)=(0.d0,0.d0)
c
            yr(0,3,id)=(0.d0,0.d0)
            yt(0,3,id)=(0.d0,0.d0)
            yp(0,3,id)=(0.d0,0.d0)
c
            yr(0,2,id)=(0.d0,0.d0)
            yt(0,2,id)=(0.d0,0.d0)
            yp(0,2,id)=(0.d0,0.d0)
c
            yr(1,2,id)=(0.d0,0.d0)
            yt(1,2,id)=(0.d0,0.d0)
            yp(1,2,id)=(0.d0,0.d0)
          enddo
c
          do id=1,nd
c
c           m = 0
c
            ca=dcmplx(1.d0/3.d0,0.d0)
            yr(0,1,id)=yr(0,1,id-1)-ca*yr(1,1,id-1)
            yr(0,4,id)=yr(0,4,id-1)-ca*yr(1,4,id-1)
            do ldeg=1,ldegf-id
              ca=dcmplx(dble(ldeg+1)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg)/dble(2*ldeg-1),0.d0)
              yr(ldeg,1,id)=yr(ldeg,1,id-1)-ca*yr(ldeg+1,1,id-1)
     &                     -cb*yr(ldeg-1,1,id-1)
              yr(ldeg,4,id)=yr(ldeg,4,id-1)-ca*yr(ldeg+1,4,id-1)
     &                     -cb*yr(ldeg-1,4,id-1)
            enddo
c
c           m = 1
c
            do ldeg=1,ldegf-id
              ca=dcmplx(dble(ldeg+2)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-1)/dble(2*ldeg-1),0.d0)
              yt(ldeg,1,id)=yt(ldeg,1,id-1)-ca*yt(ldeg+1,1,id-1)
     &                     -cb*yt(ldeg-1,1,id-1)
              yt(ldeg,4,id)=yt(ldeg,4,id-1)-ca*yt(ldeg+1,4,id-1)
     &                     -cb*yt(ldeg-1,4,id-1)
              yr(ldeg,3,id)=yr(ldeg,3,id-1)-ca*yr(ldeg+1,3,id-1)
     &                     -cb*yr(ldeg-1,3,id-1)
              yt(ldeg,3,id)=yt(ldeg,3,id-1)-ca*yt(ldeg+1,3,id-1)
     &                     -cb*yt(ldeg-1,3,id-1)
              yp(ldeg,3,id)=yp(ldeg,3,id-1)-ca*yp(ldeg+1,3,id-1)
     &                     -cb*yp(ldeg-1,3,id-1)
            enddo
c
c           m = 2
c
            do ldeg=2,ldegf-id
              ca=dcmplx(dble(ldeg+3)/dble(2*ldeg+3),0.d0)
              cb=dcmplx(dble(ldeg-2)/dble(2*ldeg-1),0.d0)
              yr(ldeg,2,id)=yr(ldeg,2,id-1)-ca*yr(ldeg+1,2,id-1)
     &                     -cb*yr(ldeg-1,2,id-1)
              yt(ldeg,2,id)=yt(ldeg,2,id-1)-ca*yt(ldeg+1,2,id-1)
     &                     -cb*yt(ldeg-1,2,id-1)
              yp(ldeg,2,id)=yp(ldeg,2,id-1)-ca*yp(ldeg+1,2,id-1)
     &                     -cb*yp(ldeg-1,2,id-1)
            enddo
          enddo
c
          do ir=1,nr
            id=idr(ir)
            do ldeg=ldegtap(1,ir),ldegtap(4,ir)
              cp0=dcmplx(tap(ldeg,ir)*plm(ldeg,0,ir),0.d0)
              cp1=dcmplx(tap(ldeg,ir)*plm(ldeg,1,ir),0.d0)
              cp2=dcmplx(tap(ldeg,ir)*plm(ldeg,2,ir),0.d0)
              vgrnr(lf,ir,1)=vgrnr(lf,ir,1)+yr(ldeg,1,id)*cp0
              vgrnr(lf,ir,2)=vgrnr(lf,ir,2)+yr(ldeg,2,id)*cp2*ssd(ir)**2
              vgrnr(lf,ir,3)=vgrnr(lf,ir,3)+yr(ldeg,3,id)*cp1*ssd(ir)
              vgrnr(lf,ir,4)=vgrnr(lf,ir,4)+yr(ldeg,4,id)*cp0
c
              vgrnt(lf,ir,1)=vgrnt(lf,ir,1)+yt(ldeg,1,id)*cp1*ssd(ir)
              vgrnt(lf,ir,2)=vgrnt(lf,ir,2)+yt(ldeg,2,id)*cp2*ssd(ir)
              vgrnt(lf,ir,3)=vgrnt(lf,ir,3)+yt(ldeg,3,id)*cp1
              vgrnt(lf,ir,4)=vgrnt(lf,ir,4)+yt(ldeg,4,id)*cp1*ssd(ir)
c
              vgrnp(lf,ir,2)=vgrnp(lf,ir,2)+yp(ldeg,2,id)*cp2*ssd(ir)
              vgrnp(lf,ir,3)=vgrnp(lf,ir,3)+yp(ldeg,3,id)*cp1
            enddo
            if(id.gt.0)then
              cfac=(1.d0,0.d0)/ssf(ir)**id
              do istp=1,4
                vgrnr(lf,ir,istp)=vgrnr(lf,ir,istp)*cfac
                vgrnt(lf,ir,istp)=vgrnt(lf,ir,istp)*cfac
              enddo
              do istp=2,3
                vgrnp(lf,ir,istp)=vgrnp(lf,ir,istp)*cfac
              enddo
            endif
          enddo
          write(*,'(i6,a,f10.4,a,i5)')lf,'.',1.0d+03*f,
     &                        ' mHz: spectra read: ',ldegf
        enddo
        close(21)
        close(22)
        close(23)
        write(*,'(i6,a)')lf-1,' spectra read from '
     &                      //grnfile(ig)(1:40)
c
c       FFT
c
        do lf=nfcut+1,nf
          sgrn(lf)=(0.d0,0.d0)
        enddo
c
        do ir=1,nr
          tred=dble(idnint(tstart+tp(ir,ig)))
          tredr4=sngl(tred)
          write(30)tredr4
c
c         time domain Green's functions
c
c---------explosion source
c         vertival component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnr(lf,ir,1)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         radial component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnt(lf,ir,1)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c---------strike-slip source
c         vertival component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnr(lf,ir,2)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         radial component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnt(lf,ir,2)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         tangential component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnp(lf,ir,2)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c---------dip-slip source
c         vertival component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnr(lf,ir,3)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         radial component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnt(lf,ir,3)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         tangential component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnp(lf,ir,3)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c---------clvd source
c         vertival component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnr(lf,ir,4)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
c
c         radial component
c
          do lf=1,nfcut
            sgrn(lf)=wvf(lf)*vgrnt(lf,ir,4)
          enddo
          call qpfftinv(ierr)
          write(30)(tgrn(it),it=1,ntout)
        enddo
c
        close(30)
c
      enddo
c
      deallocate(idr,ssd,ssf,wvf,bpf,tap,plm,
     &           sgrn,tgrn,vgrnr,vgrnt,vgrnp)
c
      return
      end