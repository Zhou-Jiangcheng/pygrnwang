      module qpalloc
c===================================================================
c     allocatable variables
c===================================================================
c
c     LEGENDRE POLYNOMIALS TABLES
c     ===========================
c
c     plm = associated Legendre polynomials divided by sin(x)**m
c
      double precision, allocatable:: plm(:,:,:)
c
c     SPHERICAL BESSEL FUNCTIONS
c     ==========================
c
c     lylwpsv = starting layer number for psv solution
c     lylwsh = starting layer number for sh solution
c
c     integer lylwsh(0:ldegmax),lylwpsv(0:ldegmax)
c     integer lyupatm(0:ldegmax)
c     zj(x)=x*j_(n+1)(x)/[j_n(x)+ix*j_(n+1)(x)]
c     wj(xa,xb)=ln{[j_n(xa)+i*xa*j_(n+1)(xa)]/[j_n(xb)+i*xb*j_(n+1)(xb)]}
c     zh(x)=x*h_(n+1)(x)/h_n(x), where h_n(x) = j_n(x) + sign[im(x)]iy_n(x)
c     wh(xa,xb)=ln[h_n(xa)/h_n(xb)]
c
c      double complex zjup(0:ldegmax,lymax,3),zjlw(0:ldegmax,lymax,3)
c      double complex zhup(0:ldegmax,lymax,3),zhlw(0:ldegmax,lymax,3)
c      double complex wj(0:ldegmax,lymax,3),wh(0:ldegmax,lymax,3)
      integer, allocatable:: lylwsh(:)
      integer, allocatable:: lylwpsv(:)
      integer, allocatable:: lyupatm(:)
      double complex, allocatable:: zjup(:,:,:)
      double complex, allocatable:: zjlw(:,:,:)
      double complex, allocatable:: zhup(:,:,:)
      double complex, allocatable:: zhlw(:,:,:)
      double complex, allocatable:: wj(:,:,:)
      double complex, allocatable:: wh(:,:,:)
c
c     same as above but with gravity effect
c
c      double complex zjupg(0:ldegmax),zjlwg(0:ldegmax)
c      double complex zhupg(0:ldegmax),zhlwg(0:ldegmax)
c      double complex wjg(0:ldegmax),whg(0:ldegmax)
      double complex, allocatable:: zjupg(:)
      double complex, allocatable:: zjlwg(:)
      double complex, allocatable:: zhupg(:)
      double complex, allocatable:: zhlwg(:)
      double complex, allocatable:: wjg(:)
      double complex, allocatable:: whg(:)
c
c     GREEN FUNCTION PARAMETERS
c     =========================
c
c     yr,yt,yp = velocity Green's function spectra
c     (r,t,p = spherical coordinate system: radius, theta, phi)
c
c      complex*8 yr(0:ldegmax,4,0:ndmax),yt(0:ldegmax,4,0:ndmax)
c      complex*8 yp(0:ldegmax,4,0:ndmax)
c      double complex vgrnr(nfmax,nrmax,4),vgrnt(nfmax,nrmax,4)
c      double complex vgrnp(nfmax,nrmax,4)
      complex*8, allocatable:: yr(:,:,:)
      complex*8, allocatable:: yt(:,:,:)
      complex*8, allocatable:: yp(:,:,:)
      double complex, allocatable:: vgrnr(:,:,:)
      double complex, allocatable:: vgrnt(:,:,:)
      double complex, allocatable:: vgrnp(:,:,:)
c
      end module
