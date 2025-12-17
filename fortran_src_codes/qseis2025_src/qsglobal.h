c     GLOBAL INDEX PARAMETERS FOR DEFINING ARRAYS
c     ===========================================
c     nzmax: max. interface index;
c     lmax: max. number of total homogeneous layers (lmax <= nzmax-2);
c     nrmax: max. number of traces;
c     nfmax: max. number of frequency samples
c
      integer*4 nzmax,lmax,nrmax,nfmax,ndtransmax
      parameter(lmax=500)
      parameter(nzmax=lmax+2)
      parameter(nrmax=101,nfmax=4096)
      parameter(ndtransmax=4)
c
c     INDEX PARAMETERS FOR SEISMOMETER CHARACTERISTICS
c     ================================================
c     (max. number of roots and poles)
c
      integer*4 nrootmax,npolemax
      parameter(nrootmax=10,npolemax=10)
c
c     EARTH RADIUS IN METER
c     =====================
c
      real*8 rr0,km2m
      parameter(rr0=6.371d+06,km2m=1.0d+03)
c
c     ATMOSPHERIC PARAMETERS
c     ======================
c     real*8 roair,vpair
c     parameter(roair=0.1300d+01,vpair=0.3318d+03)
      real*8 roair,vpair
      parameter(roair=0.d0,vpair=0.d0)
c
c     THE MININUM VS/VP RATIO: VSPMIN
c     ===============================
c     (if vs/vp < vspmin, then fluid medium is assumed)
c
      real*8 vspmin
      parameter(vspmin=0.05d0)
c
c     FOR FLAT-EARTH-TRANSFORMATION
c     =============================
c
      integer*4 ndens
      parameter(ndens=1)
c
c     PARAMETERS FOR ESTIMATE WAVENUMBER
c     ===========================================
c
      integer*4 nk0max,nbsjmax
      parameter(nk0max=10000,nbsjmax=nk0max*8)
      real*8 epswv,rd2r
c      parameter(epswv=1.0d-06,rd2r=5.0d-02)
      common /wavenumber/ epswv,rd2r
c
      complex*16 accair,cvpair,kpair,comega
      common /airpara/ accair,cvpair,kpair,comega
c
c     zr: receiver depth
c     lzr: sublayer no of receiver
c
      integer*4 lzr,lzrrs
      real*8 zr,zrrs
      common /receiver/ zr,zrrs,lzr,lzrrs
c
      integer*4 nr
      real*8 r(nrmax)
      common /distance/ r,nr
c
      integer*4 lp,nno(nzmax)
      real*8 hp(nzmax)
      common /sublayer/ hp,lp,nno
c
      integer*4 lprs,nnors(nzmax)
      real*8 hprs(nzmax)
      common /rssublayer/ hprs,lprs,nnors
c
c     original model parameters
c
      integer*4 l0
      real*8 z1(lmax),z2(lmax),ro1(lmax),ro2(lmax)
      real*8 vp1(lmax),vp2(lmax),vs1(lmax),vs2(lmax)
      real*8 qp1(lmax),qp2(lmax),qs1(lmax),qs2(lmax)
      common /model0/z1,z2,ro1,ro2,vp1,vp2,vs1,vs2,qp1,qp2,qs1,qs2,l0
c
c     original model parameters (receiver site)
c
      integer*4 l0rs
      real*8 z1rs(lmax),z2rs(lmax),ro1rs(lmax),ro2rs(lmax)
      real*8 vp1rs(lmax),vp2rs(lmax),vs1rs(lmax),vs2rs(lmax)
      real*8 qp1rs(lmax),qp2rs(lmax),qs1rs(lmax),qs2rs(lmax)
      common /rsmodel0/z1rs,z2rs,ro1rs,ro2rs,vp1rs,vp2rs,vs1rs,vs2rs,
     +                 qp1rs,qp2rs,qs1rs,qs2rs,l0rs
c       
c     layered model parameter:
c     n0: number of homogeneous layers
c
      integer*4 n0
      real*8 h(lmax),ro(lmax),vp(lmax),vs(lmax)
      real*8 qp(lmax),qs(lmax)
      common /model/ h,ro,vp,vs,qp,qs,n0
c
      complex*16 acc(lmax),kp(lmax),ks(lmax),cla(lmax),cmu(lmax)
      complex*16 cvp(lmax),cvs(lmax),wa(lmax),wb(lmax)
      common /cpara/ acc,kp,ks,cla,cmu,cvp,cvs,wa,wb
c       
c     layered model parameter:
c     n0rs: number of homogeneous layers (receiver site)
c
      integer*4 n0rs
      real*8 hrs(lmax),rors(lmax),vprs(lmax),vsrs(lmax)
      real*8 qprs(lmax),qsrs(lmax)
      common /rsmodel/ hrs,rors,vprs,vsrs,qprs,qsrs,n0rs
c
      complex*16 accrs(lmax),kprs(lmax),ksrs(lmax),
     +               clars(lmax),cmurs(lmax)
      complex*16 cvprs(lmax),cvsrs(lmax),wars(lmax),wbrs(lmax)
      common /rscpara/ accrs,kprs,ksrs,clars,cmurs,cvprs,cvsrs,wars,wbrs
c
c     for partial solution
c
      integer*4 ipartial,npar,ipsv(nzmax)
      real*8 zup(nzmax),zlow(nzmax)
      common /partial/ zup,zlow,ipartial,npar,ipsv
      logical*2 pup(nzmax),pdw(nzmax),svup(nzmax),svdw(nzmax),sh(nzmax)
      common /psvfilter/ pup,pdw,svup,svdw,sh
c
c     slowness cut-offs
c
      real*8 slw(4)
      logical*2 fullwave,calsh
      common /slwcutoffs/ slw,fullwave,calsh
c
c     source parameters
c
      integer*4 ls,ms(6),ics(6)
      real*8 zs
      real*8 sfct0(6,6),sfct1(6,6)
      common /source/ zs,sfct0,sfct1,ls,ms,ics
c
c     path filtering
c
      integer*4 iflat,ipath,lpath,isurf,ndtrans
      real*8 pathdepth
      common /pathfilter/ pathdepth,iflat,ipath,lpath,isurf,ndtrans
c
      integer*4 nt,nf
      real*8 dt,df,fi
      common /sampling/ dt,df,fi,nt,nf
c
      real*8 tstart,twindow,tau,v0
      common /tparas/ tstart,twindow,tau,v0
c
      real*8 mtensor(6),azimuth(nrmax)
      common /eqparas/ mtensor,azimuth
c
      integer*4 nnmax,nn0,iexist,wdeg
      parameter(nnmax=1024)
      real*8 wv0(nnmax)
      common /wavelets/ wv0,nn0,iexist,wdeg
c
c     seismometer filtering
c
      integer*4 nroot,npole
      complex*16 asm
      complex*16 root(nrootmax),pole(npolemax) 
      common /seismometer/ root,pole,asm,nroot,npole
c
c     table of J_n(x), n = -1, 0, 1, 2, 3
c
      real*8 bsj(nbsjmax,-1:3,nrmax),geospr(nrmax)
      common /bessels/ bsj,geospr
c
c     green's functions
c
      complex*16 grns(nfmax,19,nrmax,7)
      common /grnfcts/ grns
c
c     title text
c
      character*1 comptxt(4),varbtxt
      character*4 rcvtxt(nrmax)
      common /title/ comptxt,varbtxt,rcvtxt
c
c     input and output data files
c
      character*110 inputfile
      common /inputdata/ inputfile
      integer*4 ssel(7),fsel(19,7),flen(19,7),outsel(5)
      character*113 outfile(19,7)
      common /outsel/ ssel,fsel,flen,outsel
      common /outdata/ outfile
