      subroutine qpdocu(ierr)
      use qpalloc
      implicit none
      integer*4 ierr
c
      integer*4 i,j,k,ly,ig,flen,slen
      real*8 f
c      real, dimension (:,:,:), allocatable :: wdarray
      logical*2 head,diff
      character*24 fdate
c
      open(20,file=infofile,status='unknown')
      write(20,'(a)')'#  Space-time domain Green function database'
      write(20,'(a)')'#  calculated on '//fdate()
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  data_type size: REAL*4'
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  time_window[s] sampling[s] number_of_samples'
      write(20,'(2f10.2,i6)')dble(ntout-1)*dtout,dtout,ntout
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  number of distances, list of distances[km]'
      write(20,'(i8)')nr
      j=0
      do i=1,nr/5
        write(20,'(5f10.3)')(dism(k)/KM2M,k=j+1,j+5)
        j=j+5
      enddo
      if(j.lt.nr)then
        write(20,'(5f10.3)')(dism(k)/KM2M,k=j+1,nr)
      endif
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  each binary file includes all fundamental'
      write(20,'(a)')'#  velocity Green functions for a given source'
      write(20,'(a)')'#  depth data are stored in the following order:'
      write(20,'(a)')'#  {distance'
      write(20,'(a)')'#   {time_reduction,'
      write(20,'(a)')'#    mechanism(expl.,strike-slip,dip-slip,clvd)'
      write(20,'(a)')'#     {component(vertical,radial,tangential*)'
      write(20,'(a)')'#      {time_series}}}}'
      write(20,'(a)')'#  note: time reduction [s] needs (4)+4+(4)'
      write(20,'(a)')'#        bytes, each time series (seimogram)'
      write(20,'(a)')'#        needs (4)+number_of_samples*4+(4)'
      write(20,'(a)')'#        bytes, where (4) bytes are for the'
      write(20,'(a)')'#        reserved places before and after each'
      write(20,'(a)')'#        record'
      write(20,'(a)')'#        *)no tangential component for explosion'
      write(20,'(a)')'#          and clvd sources'
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  number of source_depths / files,'
      write(20,'(a)')'#  list of source_depths / filenames'
      write(20,'(i10)')ngrn
      do ig=1,ngrn
        do flen=80,1,-1
          if(grnfile(ig)(flen:flen).ne.' ')goto 100
        enddo
100     continue
        write(20,'(f10.1,a)')grndep(ig)/KM2M,
     &                     '  '//''''//grnfile(ig)(1:flen)//''''
      enddo
      write(20,'(a)')'#  layered earth model:'
      write(20,'(a)')'#  number of layers'
      write(20,'(i5)')ly0
      write(20,'(a)')'#  no thickness[km]'
     &             //'      vp[km/s]      vs[km/s]'
     &             //'   rho[g/cm^3]        qp        qs'
      f=FSBREF
      call qpqmodel(f)
      do ly=1,ly0
        write(20,'(i5,4f14.4,2f10.2)')ly,(rrup(ly)-rrlw(ly))/km2m,
     &    dreal(cvp(ly))/km2m,dreal(cvs(ly))/km2m,dreal(cro(ly))/km2m,
     &    0.5d0*(qpup(ly)+qplw(ly)),0.5d0*(qsup(ly)+qslw(ly))
      enddo
c
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#  binary files of arrival time [s], takeoff'
      write(20,'(a)')'#  angle [deg] and slowness [s/m] for P and S'
      write(20,'(a)')'#  waves. data are stored in the following order:'
      write(20,'(a)')'#  {green_function_depth'
      write(20,'(a)')'#   {green_function_distance,'
      write(20,'(a)')'#      {onset,takeoff, slowness}}}'
      write(20,'(a)')'#  note: each table needs (4)+number_of_distances'
      write(20,'(a)')'#        number_of_depths*12+(4)bytes, where (4)'
      write(20,'(a)')'#        bytes are for the reserved places before'
      write(20,'(a)')'#         and after each record'
      write(20,'(a)')'#================================================'
      do slen=80,1,-1
        if(stdgrndir(slen:slen).ne.' ')goto 200
      enddo
200   continue
      do flen=80,1,-1
        if(tptable(flen:flen).ne.' ')goto 300
      enddo
300   continue
      write(20,'(a)')'  '//''''//tptable(slen+1:flen)//''''
      do flen=80,1,-1
        if(tstable(flen:flen).ne.' ')goto 400
      enddo
400   continue
      write(20,'(a)')'  '//''''//tstable(slen+1:flen)//''''
      close(20)
c
      head=.true.
      diff=.true.
c
      allocate (hpmod(ly0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: hpmod not allocated!'
      allocate (vpmod(ly0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: vpmod not allocated!'
      allocate (vsmod(ly0),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: vsmod not allocated!'
      do ly=1,ly0
        hpmod(ly)=rrup(ly)-rrlw(ly)
        vpmod(ly)=dreal(cvp(ly))
        vsmod(ly)=dreal(cvs(ly))
      enddo
      
      allocate (tp(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: tp not allocated!'
      allocate(tkftp(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: tkftp not allocated!'
      allocate(slwtp(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: slwtp not allocated!'
      allocate (ts(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: ts not allocated!'
      allocate(tkfts(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: tkfts not allocated!'
      allocate(slwts(nr,ngrn),stat=ierr)
      if(ierr.ne.0)stop ' Error in qpdocu: slwts not allocated!'
c
      do ig=1,ngrn
        call tpssub(ly0,hpmod,vpmod,grndep(ig),REARTH,head,diff,
     &              nr,dism,tp(1,ig),tkftp(1,ig),slwtp(1,ig))
        call tpssub(ly0,hpmod,vsmod,grndep(ig),REARTH,head,diff,
     &              nr,dism,ts(1,ig),tkfts(1,ig),slwts(1,ig))
      enddo
      
C      allocate(wdarray(3,nr,ngrn))
C     
C      do i=1,nr
C          do ig=1,ngrn
C              wdarray(1,i,ig)=sngl(tp(i,ig))
C              wdarray(2,i,ig)=sngl(tkftp(i,ig))
C              wdarray(3,i,ig)=sngl(slwtp(i,ig))
C          enddo
C      enddo
C     
C      open(21,file=tptable,form='unformatted',status='unknown')
C      write(21)wdarray
C      close(21)
C
C      do i=1,nr
C          do ig=1,ngrn
C              wdarray(1,i,ig)=sngl(ts(i,ig))
C              wdarray(2,i,ig)=sngl(tkfts(i,ig))
C              wdarray(3,i,ig)=sngl(slwts(i,ig))
C          enddo
C      enddo
C     
C      open(22,file=tstable,form='unformatted',status='unknown')
C      write(22)wdarray
C      close(22)
C     
C      deallocate(wdarray)
      
      open(21,file=tptable,form='unformatted',status='unknown')
      write(21)((sngl(tp(i,ig)),sngl(tkftp(i,ig)),
     &            sngl(slwtp(i,ig)),i=1,nr),ig=1,ngrn)
      close(21)
      
      
      open(22,file=tstable,form='unformatted',status='unknown')
      write(22)((sngl(ts(i,ig)),sngl(tkfts(i,ig)),
     &            sngl(slwts(i,ig)),i=1,nr),ig=1,ngrn)
      close(22)
      
      return
      end
