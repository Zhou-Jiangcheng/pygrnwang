      subroutine qswvint(srate)
      implicit none
c
      real*8 srate
c
      include 'qsglobal.h'
c
      integer*4 istp,n,nrec,l,lf,lf1,i,ir,nrr,nrs
      integer*4 ik,ik1,ik2,nk1,nk2,nbsj,idtrans
      real*8 f,fcut,k,kc,dk,slwn
      real*8 pi,pi2,rr,rs,delta,cmax,ymax,yabs
      real*8 fac,zdis,rsdis,thickness,wvlen
      real*8 kcut(4),kcut1(4),kcut2(4)
      real*8 rdisk(nrmax)
      complex*16 carec,cbrec,clar,cmur,cr
      complex*16 ck,ck2,cdk,c2dk,cdk2,cfac,czdis2,cm
      complex*16 swap(nbsjmax+ndtransmax)
      complex*16 y0(9,6,nbsjmax+ndtransmax)
      complex*16 y(10,6,nbsjmax+ndtransmax)
      complex*16 cy(6,6),yb(11),cics(6),cm2(9,6)
      complex*16 jmm1,jm0,jmp1 ! bsj
      real*8 taper
c
      complex*16 c1,c2,c3
      data c1,c2,c3/(1.d0,0.d0),(2.d0,0.d0),(3.d0,0.d0)/
c
c     ics = 1  when the azmuth-factor is cos(ms*theta) for poloidal mode
c             (psv) and sin(ms*theta) for the toroidal mode (sh);
c     ics = -1 otherwise.
c
      pi=4.d0*datan(1.d0)
      pi2=2.d0*pi
c
      do istp=1,6
        cics(istp)=dcmplx(dble(ics(istp)),0.d0)
        cm2(1,istp)=dcmplx(dble(ms(istp)**2),0.d0)
        cm2(2,istp)=dcmplx(dble((ms(istp)-1)**2),0.d0)
        cm2(3,istp)=dcmplx(dble((ms(istp)+1)**2),0.d0)
        cm2(4,istp)=cm2(1,istp)
        cm2(5,istp)=dcmplx(dble(ms(istp)**2),0.d0)
        cm2(6,istp)=dcmplx(dble((ms(istp)-1)**2),0.d0)
        cm2(7,istp)=dcmplx(dble((ms(istp)+1)**2),0.d0)
        cm2(8,istp)=cm2(1,istp)
        cm2(9,istp)=cm2(1,istp)
      enddo
c
      nrec=nno(lzr)
      carec=dcmplx(1.d0/(ro(nrec)*vp(nrec)**2),0.d0)
      cbrec=dcmplx(2.d0*(vs(nrec)/vp(nrec))**2,0.d0)
      cmur=dcmplx(ro(nrec)*vs(nrec)**2,0.d0)
      clar=dcmplx(ro(nrec)*vp(nrec)**2,0.d0)-2*cmur
c
      l=min0(ls,lzr)
      n=nno(l)
      if(vs(n).le.vspmin*vp(n).or.
     &   .not.(svup(l).or.svdw(l)))then
        cmax=vp(n)
      else
        cmax=vs(n)
      endif
      do l=min0(ls,lzr)+1,max0(ls,lzr,lpath)-1
        n=nno(l)
        if(vs(n).le.vspmin*vp(n).or.
     &     .not.(svup(l).or.svdw(l)))then
          cmax=dmax1(cmax,vp(n))
        else
          cmax=dmax1(cmax,vs(n))
        endif
      enddo
c
      rsdis=0.d0
      do l=1,max0(ls,lzr,lpath)-1
        rsdis=rsdis+hp(l)
      enddo
c
      fcut=dble(nf)*df
c
      zdis=dabs(zr-zs)
      if(zdis.le.0.d0.and.r(1).le.0.d0)then
        stop 'Error in qswvint: identical '
     &     //'location of receiver and source!'
      endif
      czdis2=dcmplx(zdis*zdis,0.d0)
c
      thickness=0.d0
      do l=1,lp-1
        thickness=thickness+hp(l)
      enddo
c
      dk=pi/dmax1(srate*r(nr),3.d0*thickness)
      cdk=dcmplx(dk,0.d0)
      c2dk=dcmplx(2.d0*dk,0.d0)
      cdk2=dcmplx(dk*dk,0.d0)
c
      if(fullwave)then
        f=0.d0
        rdisk(1)=rd2r*dsqrt(zdis*zdis+r(1)*r(1))
        call qsqmodel(f)
        kc=0.d0
        ik1=1
        ymax=0.d0
101     k=dble(ik1)*dk
        call qskern(cy,f,k)
        yabs=0.d0
        do istp=1,4
          do i=1,3
            yabs=yabs+cdabs(cy(2*i-1,istp))**2
          enddo
        enddo
        yabs=dsqrt(yabs)*k*dexp(-0.5d0*(k*rdisk(1))**2)
        ymax=dmax1(ymax,yabs)
        if(yabs.gt.epswv*ymax.and.yabs.gt.0.d0.and.
     &     ik1.lt.nk0max)then
          ik1=ik1+1
          goto 101
        endif
        kcut1(1)=0.d0
        kcut1(2)=0.d0
        kcut1(3)=0.8d0*k
        kcut1(4)=k
        ! kcut1(3)=k
        ! kcut1(4)=1.25*k
c
        f=fcut
        n=nno(ls)
        wvlen=vp(n)/(f+df)
        rdisk(1)=rd2r*dmin1(dsqrt(zdis*zdis+r(1)*r(1)),wvlen)
        ik2=ik1
        ymax=0.d0
102     k=dble(ik2)*dk
        call qskern(cy,f,k)
        yabs=0.d0
        do istp=1,4
          do i=1,3
            yabs=yabs+cdabs(cy(2*i-1,istp))**2
          enddo
        enddo
        yabs=dsqrt(yabs)*k*dexp(-0.5d0*(k*rdisk(1))**2)
        ymax=dmax1(ymax,yabs)
        if(yabs.gt.epswv*ymax.and.yabs.gt.0.d0.and.
     &     ik2.lt.(nbsjmax+nk0max)/2)then
          ik2=ik2+1
          goto 102
        endif
c
        kcut2(1)=0.d0
        kcut2(2)=0.d0
        kcut2(3)=0.8d0*k
        kcut2(4)=k
        ! kcut2(3)=k
        ! kcut2(4)=1.25*k
        lf1=1
      else
        kcut1(1)=0.d0
        kcut1(2)=0.d0
        kcut1(3)=0.d0
        kcut1(4)=0.d0
        do i=1,4
          kcut2(i)=pi2*fcut*slw(i)
        enddo
        lf1=2
      endif
c
      nbsj=2+idint(dmax1(kcut1(4),kcut2(4))/dk)+ndtrans
      print *, kcut1(4)/dk, kcut2(4)/dk, ndtrans
      print *, nbsj
c
      if(nbsj.gt.nbsjmax)then
        stop ' parameter nbsjmax defined too small'
      else
        print *,' Calculate Bessel functions for x up to ',
     &       dble(nbsj)*dk*r(nr)
        print *,'nbsj',nbsj
        do ir=1,nr
          geospr(ir)=1.d0/(zdis*zdis+r(ir)*r(ir))**ndtrans
        enddo
        call qsbsj(dk,nbsj)
      endif
c
      do lf=1,nf
        do istp=1,6
          do i=1,19
            do ir=1,nr
              grns(lf,i,ir,istp)=(0.d0,0.d0)
             enddo
          enddo
        enddo
      enddo
c
      write(*,'(a,2(f10.7,a))')' Min./max. slowness at f_cut: ',
     &     1000.d0*kcut2(1)/(pi2*fcut),' / ',
     &     1000.d0*kcut2(4)/(pi2*fcut),' s/km'
c
      do lf=lf1,nf
        f=dble(lf-1)*df
        n=nno(ls)
        wvlen=vp(n)/(f+df)
        do ir=1,nr
          rdisk(ir)=rd2r*dmin1(dsqrt(zdis*zdis+r(ir)*r(ir)),wvlen)
        enddo
        call qsqmodel(f)
c
        do i=1,4
          kcut(i)=kcut1(i)
     &           +(kcut2(i)-kcut1(i))*dsqrt(f**2+(pi*fi)**2)/fcut
        enddo
c
        nk2=min0(2+idint(kcut(4)/dk),nbsj-ndtrans)
        nk1=min0(1+idint(kcut(1)/dk),nk2)
c
        ik1=max0(1,nk1-ndtrans)
        ik2=nk2+ndtrans
        do ik=ik1,ik2
          k=dble(ik)*dk
          ck=dcmplx(k,0.d0)
          call qskern(cy,f,k)
          do istp=1,6
            y0(1,istp,ik)=cy(1,istp)
            y0(2,istp,ik)=( cy(3,istp)+cics(istp)*cy(5,istp))/c2
            y0(3,istp,ik)=( cy(3,istp)-cics(istp)*cy(5,istp))/c2
            y0(4,istp,ik)=carec*cy(2,istp)-cbrec*ck*cy(3,istp)
            y0(5,istp,ik)=cy(2,istp)
            y0(6,istp,ik)=( cy(4,istp)+cics(istp)*cy(6,istp))/c2
            y0(7,istp,ik)=( cy(4,istp)-cics(istp)*cy(6,istp))/c2
            y0(8,istp,ik)=ck*cy(3,istp)
            y0(9,istp,ik)=ck*cy(5,istp)
          enddo
        enddo
c
        do ir=1,nr
          cr=dcmplx(r(ir),0.d0)
          ik1=max0(1,nk1-ndtrans)
          ik2=nk2+ndtrans
          do ik=ik1,ik2
            k=dble(ik)*dk
            do istp=1,6
              do i=1,9
                y(i,istp,ik)=y0(i,istp,ik)
     &                   *dcmplx(dexp(-0.5d0*(k*rdisk(ir))**2),0.d0)
              enddo
            enddo
          enddo
          do idtrans=1,ndtrans
            ik1=max0(1,nk1-ndtrans)+idtrans
            ik2=nk2+ndtrans-idtrans
            do istp=1,6
              do i=1,9
                do ik=ik1-1,ik2+1
                  swap(ik)=y(i,istp,ik)
                enddo
                do ik=ik1,ik2
                  ck=dcmplx(dble(ik)*dk,0.d0)
                  ck2=ck*ck
                  y(i,istp,ik)=swap(ik)*(czdis2+cm2(i,istp)/ck2)
     &              -(swap(ik+1)-swap(ik-1))/c2dk/ck
     &              -(swap(ik+1)-c2*swap(ik)+swap(ik-1))/cdk2
                enddo
              enddo
            enddo
          enddo
c
          ik1=max0(1,nk1+ndtrans)
          ik2=nk2
          do ik=ik1,ik2
            k=dble(ik)*dk
            ck=dcmplx(k,0.d0)
            cfac=dcmplx(k*dk
     &                 *taper(k,kcut(1),kcut(2),kcut(3),kcut(4)),0.d0)
c
            do istp=1,6
              do i=1,9
                y(i,istp,ik)=y(i,istp,ik)*cfac
              enddo
              jmm1 = dcmplx(bsj(ik,ms(istp)-1,ir),0.d0)
              jm0  = dcmplx(bsj(ik,ms(istp),ir),0.d0)
              jmp1 = dcmplx(bsj(ik,ms(istp)+1,ir),0.d0)
              cm=dcmplx(ms(istp), 0.d0)
c
              yb(1)=y(1,istp,ik)*jm0
              yb(2)=y(2,istp,ik)*jmm1
              yb(3)=y(3,istp,ik)*jmp1
              yb(4)=y(4,istp,ik)*jm0
              yb(5)=y(5,istp,ik)*jm0
              yb(6)=y(6,istp,ik)*jmm1
              yb(7)=y(7,istp,ik)*jmp1
              yb(8)=y(8,istp,ik)*jm0
              yb(9)=y(9,istp,ik)*jm0
              yb(10)=y(1,istp,ik)*jmm1
              yb(11)=y(1,istp,ik)*jmp1
c
c tz,tr,tt disp
              if(outsel(1).eq.1)then
                grns(lf,1,ir,istp)=grns(lf,1,ir,istp)+yb(1)
                grns(lf,2,ir,istp)=grns(lf,2,ir,istp)+yb(2)-yb(3)
                grns(lf,3,ir,istp)=grns(lf,3,ir,istp)
     &           -cics(istp)*(yb(2)+yb(3))
              endif
c tv volume
              if(outsel(2).eq.1)then
                grns(lf,4,ir,istp)=grns(lf,4,ir,istp)+yb(4)
              endif
c ezz,ezr,ezt strain
              if(outsel(3).eq.1)then
                grns(lf,5,ir,istp)=grns(lf,5,ir,istp)
     &           +(yb(5)-clar*yb(4))/(c2*cmur)
                grns(lf,6,ir,istp)=grns(lf,6,ir,istp)
     &           +(yb(6)-yb(7))/(c2*cmur)
                grns(lf,7,ir,istp)=grns(lf,7,ir,istp)
     &           -cics(istp)*(yb(6)+yb(7))/(c2*cmur)
c ett, err, ert strain
                if (r(ir).gt.0.d0)then
                  grns(lf,8,ir,istp)=grns(lf,8,ir,istp)
     &             -((cm-c1)*yb(2)+(cm+c1)*yb(3))/cr
                  grns(lf,9,ir,istp)=grns(lf,9,ir,istp)
     &             -yb(8)+((cm-c1)*yb(2)+(cm+c1)*yb(3))/cr
                  grns(lf,10,ir,istp)=grns(lf,10,ir,istp)+yb(9)/c2
     &                 -(cics(istp)*cm*(yb(2)-yb(3))
     &                 +(-cics(istp))*(yb(2)+yb(3)))/cr
                else if(ms(istp).eq.0.d0)then
                  grns(lf,8,ir,istp)=grns(lf,8,ir,istp)
     &             -ck/c2*(y(2,istp,ik)+y(3,istp,ik))
                  grns(lf,9,ir,istp)=grns(lf,9,ir,istp)
     &             -y(8,istp,ik)+ck/c2*(y(2,istp,ik)+y(3,istp,ik))
                  grns(lf,10,ir,istp)=grns(lf,10,ir,istp)
     &             +y(9,istp,ik)/c2
     &             +cics(istp)*ck/c2*(y(2,istp,ik)+y(3,istp,ik))
                else if(ms(istp).eq.1.d0)then
                  grns(lf,8,ir,istp)=grns(lf,8,ir,istp)+(0.d0,0.d0)
                  grns(lf,9,ir,istp)=grns(lf,9,ir,istp)+(0.d0,0.d0)
                  grns(lf,10,ir,istp)=grns(lf,10,ir,istp)+(0.d0,0.d0)
                else if(ms(istp).eq.2.d0)then
                  grns(lf,8,ir,istp)=grns(lf,8,ir,istp)
     &             -ck/c2*y(2,istp,ik)
                  grns(lf,9,ir,istp)=grns(lf,9,ir,istp)
     &             +ck/c2*y(2,istp,ik)
                  grns(lf,10,ir,istp)=grns(lf,10,ir,istp)
     &             -cics(istp)*ck/c2*y(2,istp,ik)
                endif
              endif
c szz,szr,szt stress
              if(outsel(4).eq.1)then
                grns(lf,11,ir,istp)=grns(lf,11,ir,istp)+yb(5)
                grns(lf,12,ir,istp)=grns(lf,12,ir,istp)+yb(6)-yb(7)
                grns(lf,13,ir,istp)=grns(lf,13,ir,istp)
     &           -cics(istp)*(yb(6)+yb(7))
c stt, srr, srt stress
                if (r(ir).gt.0.d0)then
                  grns(lf,14,ir,istp)=grns(lf,14,ir,istp)+clar*yb(4)
     &             +2*cmur*(-((cm-c1)*yb(2)+(cm+c1)*yb(3))/cr)
                  grns(lf,15,ir,istp)=grns(lf,15,ir,istp)+clar*yb(4)
     &             +2*cmur*(-yb(8)+((cm-c1)*yb(2)+(cm+c1)*yb(3))/cr)
                  grns(lf,16,ir,istp)=grns(lf,16,ir,istp)+2*cmur*(yb(9)/c2
     &              -(cics(istp)*cm*(yb(2)-yb(3))
     &              +(-cics(istp)*(yb(2)+yb(3))))/cr)
                else if(ms(istp).eq.0.d0)then
                  grns(lf,14,ir,istp)=grns(lf,14,ir,istp)+clar*yb(4)
     &              +2*cmur*(-ck/c2*(y(2,istp,ik)+y(3,istp,ik)))
                  grns(lf,15,ir,istp)=grns(lf,15,ir,istp)+clar*yb(4)
     &          +2*cmur*(-y(8,istp,ik)+ck/c2*(y(2,istp,ik)+y(3,istp,ik)))
                  grns(lf,16,ir,istp)=grns(lf,16,ir,istp)+2*cmur*
     &              (y(9,istp,ik)/c2+cics(istp)*ck/c2*(y(2,istp,ik)
     &              +y(3,istp,ik)))
                else if(ms(istp).eq.1.d0)then
                  grns(lf,14,ir,istp)=grns(lf,14,ir,istp)+dcmplx(0.d0,0.d0)
                  grns(lf,15,ir,istp)=grns(lf,15,ir,istp)+dcmplx(0.d0,0.d0)
                  grns(lf,16,ir,istp)=grns(lf,16,ir,istp)+dcmplx(0.d0,0.d0)
                else if(ms(istp).eq.2.d0)then
                  grns(lf,14,ir,istp)=grns(lf,14,ir,istp)+2*cmur*
     &               (-ck/c2*y(2,istp,ik))
                  grns(lf,15,ir,istp)=grns(lf,15,ir,istp)+2*cmur*
     &               (+ck/c2*y(2,istp,ik))
                  grns(lf,16,ir,istp)=grns(lf,16,ir,istp)+2*cmur*
     &              (-cics(istp)*ck/c2*y(2,istp,ik))
                endif
              endif
c oz,or,ot rotation
            if(outsel(5).eq.1)then
              grns(lf,17,ir,istp)=grns(lf,17,ir,istp)+yb(9)/c2
              grns(lf,18,ir,istp)=grns(lf,18,ir,istp)
     &         -c1/c2*cics(istp)*ck*(yb(10)+yb(11))
     &         +cics(istp)*(yb(6)+yb(7))/(2*cmur)
              grns(lf,19,ir,istp)=grns(lf,19,ir,istp)
     &         -c1/c2*ck*(yb(10)-yb(11))
     &         +(yb(6)-yb(7))/(2*cmur)
            endif
            enddo
          enddo
        enddo
c
        write(*,'(i6,a,E13.6,a,i7)')lf,'.',f,
     &      'Hz: slowness samples = ',1+nk2-nk1-ndtrans
      enddo
c
      if(iflat.eq.1)then
c
c       amplitude correction when using the flat-earth transform
c       see Mueller (1977) for n = -2
c
        rs=rr0*dexp(-zs/rr0)
        rr=rr0*dexp(-zrrs/rr0)
        nrs=5-ndens
        nrr=3-ndens
        do ir=1,nr
          if(r(ir).gt.0.d0)then
            delta=r(ir)/rr0
            fac=delta/dsin(delta)
          else
            fac=1.d0
          endif
          cfac=dcmplx(dsqrt((rr0/rr)**nrr*(rr0/rs)**nrs*fac),0.d0)
          do istp=1,6
            do i=1,19
              do lf=lf1,nf
                grns(lf,i,ir,istp)=grns(lf,i,ir,istp)*cfac
              enddo
            enddo
          enddo
        enddo
      endif
c
      return
      end