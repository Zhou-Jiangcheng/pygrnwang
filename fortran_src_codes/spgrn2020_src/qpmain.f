      program qpmain
      use qpalloc
      implicit none
c
c     work space
c
      integer*4 i,j,ig,ierr,runtime
      integer*4 time
      character*128 inputfile
c
c     read input file file
c
      print *,'######################################################'
      print *,'#                                                    #'
      print *,'#               Welcome to the program               #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#      SSSS    PPPP      GGG     RRRR     N   N      #'
      print *,'#     S        P   P    G        R   R    NN  N      #'
      print *,'#      SSS     PPPP     G  GG    RRRR     N N N      #'
      print *,'#         S    P        G   G    R R      N   N      #'
      print *,'#     SSSS     P         GGG     R  R     N   N      #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                        for                         #'
      print *,'#              synthetic Green functions             #'
      print *,'#                      based on                      #'
      print *,'#      a spherical self-gravitating earth model      #'
      print *,'#                                                    #'
      print *,'#                  (Version 2020)                    #'
      print *,'#                                                    #'
      print *,'#                                                    #'
      print *,'#                      by                            #'
      print *,'#                 Rongjiang Wang                     #'
      print *,'#              (wang@gfz-potsdam.de)                 #'
      print *,'#                                                    #'
      print *,'#     GFZ German Research Centre for Geosciences     #'
      print *,'#            last modified: Jan 2020                 #'
      print *,'######################################################'
      print *,'                          '
      write(*,'(a,$)')' the input data file is '
      read(*,'(a)')inputfile
      runtime=time()
c
      open(10,file=inputfile,status='old')
      call qpgetinp(10)
      close(10)
c
      call qpsublayer(ierr)
c
      call qpdocu(ierr)
c
      call qpmaxdeg(ierr)
c
      igfirst=0
      do ig=ngrn,1,-1
        if(grnsel(ig).eq.1)igfirst=ig
      enddo
      iglast=ngrn+1
      do ig=1,ngrn
        if(grnsel(ig).eq.1)iglast=ig
      enddo
c
      do ig=1,ngrn
        if(grnsel(ig).eq.1)then
          lys=lygrn(ig)
          call qpgrnspec(ig)
        endif
      enddo
c
      call qpwvint(ierr)
c
      runtime=time()-runtime
      write(*,'(a)')' #############################################'
      write(*,'(a)')' #                                           #'
      write(*,'(a)')' #      End of computations with spgrn       #'
      write(*,'(a)')' #                                           #'
      write(*,'(a,i10,a)')' #       Run time: ',runtime,
     +                                           ' sec            #'
      write(*,'(a)')' #############################################'
      stop
      end
