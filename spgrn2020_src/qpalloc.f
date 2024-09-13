      module qpalloc
c
c     CONSTANTS
c     =========
      real*8 PI,PI2
      parameter(PI=3.14159265358979d0,PI2=6.28318530717959d0)
      real*8 DEG2RAD,KM2M
      parameter(DEG2RAD=1.745329251994328d-02,KM2M=1.0d+03)
      real*8 REARTH,BIGG
      parameter(REARTH=6.371d+06,BIGG=6.6732d-11)
      real*8 RESOLUT
      parameter(RESOLUT=0.01d0)
      real*8 FLTAPER
      parameter(FLTAPER=0.2d0)
      real*8 FSBLW,FSBREF
      parameter(FSBLW=0.25d-03,FSBREF=1.0d+00)
c
      real*8 expol,expos,expoa
      parameter(expol=6.0d+00,expos=16.0d+00,expoa=16.0d+00)
      real*8 slwfac
      parameter(slwfac=1.25d+00)
c
      integer*4 ldeg0,ldeg1,ndmax
      parameter(ldeg0=2500,ldeg1=50,ndmax=4)
c
      logical*2 nogravity,selpsv,selsh,dispersion,fullwave
c
      integer*4 ldegmin,ldegmax,igfirst,iglast
      integer*4 nr,ntout,l0,lys,lyr,lyos,lyob,lycm,lycc,ly0
      integer*4 ngrn,nt,ntcut,ntcutout,nf,nfcut,nbpf,ldeggr
c
      real*8 dt,df,fi,fcut,fgr,rratmos,depatmos
      real*8 slwmax,f1corner,f2corner,qsmin,mmantle
      real*8 dpr,freeairgrd,dtout,tau,tstart,tred,twindow,twinout
c
      complex*16 comi,comi2,cua(8,6)
c
      character*80 spcgrndir,stdgrndir,infofile,tptable,tstable
c
      logical*2, allocatable:: ksmallpsv(:),ksmallsh(:)
c
      integer*4, allocatable:: lygrn(:),grnsel(:),i0sub(:),
     &       ldegpsv(:),ldegsh(:),lylwsh(:),lylwpsv(:),lyupatm(:)
c
      real*4, allocatable:: tgrn(:)
c
      real*8, allocatable:: mshell(:),dswap(:),hpmod(:),
     &            vpmod(:),tp(:,:),tkftp(:,:),slwtp(:,:),
     &            vsmod(:),ts(:,:),tkfts(:,:),slwts(:,:)
      real*8, allocatable:: grndep(:),rr0(:),dism(:),disrad(:),
     &       dp0(:),dp0up(:),dp0lw(:),vp0(:),vp0up(:),vp0lw(:),
     &       vs0(:),vs0up(:),vs0lw(:),ro0(:),ro0up(:),ro0lw(:),
     &       qp0(:),qp0up(:),qp0lw(:),qs0(:),qs0up(:),qs0lw(:),
     &       rrup(:),rrlw(:),vpup(:),vplw(:),vsup(:),vslw(:),
     &       roup(:),rolw(:),qpup(:),qplw(:),qsup(:),qslw(:)
      real*8, allocatable:: disk(:),fdisk(:),dplm1(:),plm(:,:,:)
c
      complex*16, allocatable:: yr(:,:,:),yt(:,:,:),yp(:,:,:),sgrn(:),
     &                          cswap(:)
c
      complex*16, allocatable:: vgrnr(:,:,:),vgrnt(:,:,:),vgrnp(:,:,:)
      complex*16, allocatable:: crrup(:),crrlw(:),
     &       cro(:),croup(:),crolw(:),cla(:),claup(:),clalw(:),
     &       cmu(:),cmuup(:),cmulw(:),cvp(:),cvpup(:),cvplw(:),
     &       cvs(:),cvsup(:),cvslw(:),cgr(:),cgrup(:),cgrlw(:),
     &       cga(:),cgaup(:),cgalw(:)
      complex*16, allocatable:: kp(:),ks(:),kt(:),cps(:,:),cpt(:,:),
     &       mat2x2up(:,:,:),mat2x2lw(:,:,:),mat2x2inv(:,:,:),
     &       mas2x2up(:,:,:),mas2x2lw(:,:,:),mas2x2inv(:,:,:),
     &       mas4x4up(:,:,:),mas4x4lw(:,:,:),mas4x4inv(:,:,:),
     &       mas6x6up(:,:,:),mas6x6lw(:,:,:),mas6x6inv(:,:,:)
      complex*16, allocatable:: cypnorm(:,:)
      complex*16, allocatable:: zjup(:,:,:),zjlw(:,:,:),zhup(:,:,:),
     &       zhlw(:,:,:),wj(:,:,:),wh(:,:,:),zjupg(:),
     &       zjlwg(:),zhupg(:),zhlwg(:),wjg(:),whg(:)
c
      character*80, allocatable:: grnfile(:),rgrnfile(:),
     &             tgrnfile(:),pgrnfile(:),stdgrnfile(:)

c
      end module
