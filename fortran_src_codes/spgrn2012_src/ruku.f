      subroutine ruku(y,i0,j0,ly,ldeg,difmat,rr1,rr2)
      implicit none
      integer i0,j0,ly,ldeg
      double precision rr1,rr2
      double complex y(i0,j0)
      external difmat
c
      integer nrrmax
      parameter(nrrmax=4096)
c
      integer i,j,k,irr,nrr
      double precision rr,drr
      double complex cdrr,ra
      double complex y0(6,3),y1(6,3),y2(6,3),ya(6,3),yb(6,3)
      double complex yk(6,4),mat(6,6,3)
      logical again
c
      double precision eps
      double complex c1,c2,c6
      data eps/1.0d-03/
      data c1,c2,c6/(1.d0,0.d0),(2.d0,0.d0),(6.d0,0.d0)/
c
      do j=1,j0
        do i=1,i0
          y1(i,j)=(0.d0,0.d0)
          y2(i,j)=y(i,j)
          yb(i,j)=y(i,j)
        enddo
      enddo
      nrr=1
100   continue
      drr=(rr2-rr1)/dble(nrr)
      cdrr=dcmplx(drr,0.d0)
      do j=1,j0
        do i=1,i0
          y0(i,j)=y1(i,j)
          y1(i,j)=y2(i,j)
          y2(i,j)=y(i,j)
          ya(i,j)=yb(i,j)
        enddo
      enddo
      do irr=1,nrr
        rr=rr1+dble(irr-1)*drr
        call difmat(ly,ldeg,rr,mat(1,1,1))
        call difmat(ly,ldeg,rr+0.5d0*drr,mat(1,1,2))
        call difmat(ly,ldeg,rr+drr,mat(1,1,3))
        do j=1,j0
          do i=1,i0
            yk(i,1)=(0.d0,0.d0)
            do k=1,i0
              yk(i,1)=yk(i,1)+cdrr*mat(i,k,1)*y2(k,j)
            enddo
          enddo
c
          do i=1,i0
            yk(i,2)=(0.d0,0.d0)
            do k=1,i0
              yk(i,2)=yk(i,2)+cdrr*mat(i,k,2)*(y2(k,j)+yk(k,1)/c2)
            enddo
          enddo
c
          do i=1,i0
            yk(i,3)=(0.d0,0.d0)
            do k=1,i0
              yk(i,3)=yk(i,3)+cdrr*mat(i,k,2)*(y2(k,j)+yk(k,2)/c2)
            enddo
          enddo
c
          do i=1,i0
            yk(i,4)=(0.d0,0.d0)
            do k=1,i0
              yk(i,4)=yk(i,4)+cdrr*mat(i,k,3)*(y2(k,j)+yk(k,3))
            enddo
          enddo
c
          do i=1,i0
            y2(i,j)=y2(i,j)+(yk(i,1)+c2*(yk(i,2)+yk(i,3))+yk(i,4))/c6
          enddo
        enddo
      enddo
c
      again=.false.
      if(nrr.le.64)then
        do j=1,j0
          do i=1,i0
            yb(i,j)=y2(i,j)
            if(cdabs(yb(i,j)-ya(i,j)).gt.eps*cdabs(yb(i,j)))again=.true.
          enddo
        enddo
      else if(nrr.le.nrrmax)then
        do j=1,j0
          do i=1,i0
            ra=c2*y1(i,j)-y0(i,j)-y2(i,j)
            if(cdabs(ra).gt.eps*cdabs(y2(i,j)))then
              yb(i,j)=y1(i,j)+(y2(i,j)-y1(i,j))*(y1(i,j)-y0(i,j))/ra
            else
              yb(i,j)=y2(i,j)
            endif
            if(cdabs(yb(i,j)-ya(i,j)).gt.eps*cdabs(yb(i,j)))again=.true.
          enddo
        enddo
      else
        print *, ' Warning in ruku: Convergence problem!'
      endif
      if(again)then
        nrr=2*nrr
        goto 100
      endif
c
      do j=1,j0
        do i=1,i0
          y(i,j)=yb(i,j)
        enddo
      enddo
c
      return
      end