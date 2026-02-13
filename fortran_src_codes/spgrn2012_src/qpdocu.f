      subroutine qpdocu(ierr)
      implicit none
      integer ierr
c
      include 'qpglobal.h'
c
      integer i,j,k,ig,flen
      character*24 fdate
c
      open(20,file=infofile,status='unknown')
      write(20,'(a)')'#   Space-time domain Green function database'
      write(20,'(a)')'#   calculated on '//fdate()
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   data_type size: REAL*4'
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   time_window[s] sampling[s] number_of_samples'
      write(20,'(2f10.2,i6)')dble(ntout-1)*dtout,dtout,ntout
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   number of distances, list of distances'
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
      write(20,'(a)')'#   each file includes all fundamental velocity'
      write(20,'(a)')'#   Green functions for a given source depth'
      write(20,'(a)')'#   data are stored in the following order:'
      write(20,'(a)')'#   distance'
      write(20,'(a)')'#    {time_reduction,'
      write(20,'(a)')'#     mechanism(expl.,strike-slip,dip-slip,clvd)'
      write(20,'(a)')'#      {component(vertical,radial,tangential)'
      write(20,'(a)')'#       {time_series}}}'
      write(20,'(a)')'#   note: time reduction [s] needs (4)+4+(4)'
      write(20,'(a)')'#         bytes, each time series (seimogram)'
      write(20,'(a)')'#         needs (4)+number_of_samples*4+(4)'
      write(20,'(a)')'#         bytes, where (4) bytes are for the'
      write(20,'(a)')'#         reserved places before and after each'
      write(20,'(a)')'#         record'
      write(20,'(a)')'#================================================'
      write(20,'(a)')'#   number of source_depths / files,'
      write(20,'(a)')'#   list of source_depths / filenames'
      write(20,'(i10)')ngrn
      do ig=1,ngrn
        do flen=128,1,-1
          if(stdgrnfile(ig)(flen:flen).ne.' ')goto 100
        enddo
100     continue
        write(20,'(f10.1,a)')grndep(ig)/KM2M,
     &                     '  '//''''//stdgrnfile(ig)(1:flen)//''''
      enddo
      close(20)
      return
      end