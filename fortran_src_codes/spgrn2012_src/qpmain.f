      program qpmain
      use qpalloc
      implicit none
c
      include 'qpglobal.h'
c
c     work space
c
      integer i,j,ig,ldegcut,ierr,runtime
      integer time
      character*128 arg,inputindex,inputfile,leftpart,rightpart
c     
c
c
      write(*,'(a)') 'Please enter the input file name:'
      read(*,'(a)') inputfile
      open(10,file=inputfile,status='old')
      call qpgetinp(10)
      close(10)
c
      call qpsublayer(ierr)
c
      call qpdocu(ierr)
c
      call qpbounds(ldegcut)
      allocate (plm(0:ldegcut,0:2,nr),stat=ierr)
      allocate (lylwsh(0:ldegcut),stat=ierr)
      allocate (lylwpsv(0:ldegcut),stat=ierr)
      allocate (lyupatm(0:ldegcut),stat=ierr)
      allocate (zjup(0:ldegcut,ly0,3),stat=ierr)
      allocate (zjlw(0:ldegcut,ly0,3),stat=ierr)
      allocate (zhup(0:ldegcut,ly0,3),stat=ierr)
      allocate (zhlw(0:ldegcut,ly0,3),stat=ierr)
      allocate (wj(0:ldegcut,ly0,3),stat=ierr)
      allocate (wh(0:ldegcut,ly0,3),stat=ierr)
      allocate (zjupg(0:ldegcut),stat=ierr)
      allocate (zjlwg(0:ldegcut),stat=ierr)
      allocate (zhupg(0:ldegcut),stat=ierr)
      allocate (zhlwg(0:ldegcut),stat=ierr)
      allocate (wjg(0:ldegcut),stat=ierr)
      allocate (whg(0:ldegcut),stat=ierr)
      allocate (yr(0:ldegcut,4,0:ndmax),stat=ierr)
      allocate (yt(0:ldegcut,4,0:ndmax),stat=ierr)
      allocate (yp(0:ldegcut,4,0:ndmax),stat=ierr)
      allocate (vgrnr(nf,nr,4),stat=ierr)
      allocate (vgrnt(nf,nr,4),stat=ierr)
      allocate (vgrnp(nf,nr,4),stat=ierr)
c
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          call qpgrnspec(ig,ldegcut)
        endif
      enddo
      call qpwvint(ldegcut)
c
      deallocate (plm,stat=ierr)
      deallocate (lylwsh,stat=ierr)
      deallocate (lylwpsv,stat=ierr)
      deallocate (lyupatm,stat=ierr)
      deallocate (zjup,stat=ierr)
      deallocate (zjlw,stat=ierr)
      deallocate (zhup,stat=ierr)
      deallocate (zhlw,stat=ierr)
      deallocate (wj,stat=ierr)
      deallocate (wh,stat=ierr)
      deallocate (zjupg,stat=ierr)
      deallocate (zjlwg,stat=ierr)
      deallocate (zhupg,stat=ierr)
      deallocate (zhlwg,stat=ierr)
      deallocate (wjg,stat=ierr)
      deallocate (whg,stat=ierr)
      deallocate (yr,stat=ierr)
      deallocate (yt,stat=ierr)
      deallocate (yp,stat=ierr)
      deallocate (vgrnr,stat=ierr)
      deallocate (vgrnt,stat=ierr)
      deallocate (vgrnp,stat=ierr)
c
      stop
      end
