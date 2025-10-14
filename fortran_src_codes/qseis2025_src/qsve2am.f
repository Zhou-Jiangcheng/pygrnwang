      subroutine qsve2am(n,ck,z,dynamic,y,c,nj,rsite)
      implicit none
c
c     converting displacement-stress vectors to wave amplitudes
c     by Langer block-diagonal decomposition technique for com-
c     putational efficiency
c
      integer*4 n,nj
      real*8 z
      complex*16 ck
      complex*16 y(4,nj),c(4,nj)
      logical*4 dynamic,rsite
c
      include 'qsglobal.h'
c
      integer*4 i,j,m,key
      complex*16 ca,cadp,cads
      complex*16 swap(4),a(4,4),b(4,4)
c
      if(rsite)then
        ca=(0.5d0,0.d0)/accrs(n)
        cadp=ca/kprs(n)
        cads=ca/ksrs(n)
        do j=1,nj
          swap(1)=(-wars(n)*y(1,j)+ck*y(4,j))*cadp
          swap(2)=(-y(2,j)+wbrs(n)*y(3,j))*ca
          swap(3)=(-ck*y(2,j)+wars(n)*y(3,j))*cads
          swap(4)=(wbrs(n)*y(1,j)-y(4,j))*ca
c
          c(1,j)= swap(1)+swap(2)
          c(2,j)=-swap(1)+swap(2)
          c(3,j)=-swap(3)+swap(4)
          c(4,j)= swap(3)+swap(4)
        enddo
      else if(dynamic)then
        ca=(0.5d0,0.d0)/acc(n)
        cadp=ca/kp(n)
        cads=ca/ks(n)
        do j=1,nj
          swap(1)=(-wa(n)*y(1,j)+ck*y(4,j))*cadp
          swap(2)=(-y(2,j)+wb(n)*y(3,j))*ca
          swap(3)=(-ck*y(2,j)+wa(n)*y(3,j))*cads
          swap(4)=(wb(n)*y(1,j)-y(4,j))*ca
c
          c(1,j)= swap(1)+swap(2)
          c(2,j)=-swap(1)+swap(2)
          c(3,j)=-swap(3)+swap(4)
          c(4,j)= swap(3)+swap(4)
        enddo
      else
        call qsmatrix0(b,ck,z,n)
        do i=1,4
          do j=1,4
            a(i,j)=(0.d0,0.d0)
          enddo
          a(i,i)=(1.d0,0.d0)
        enddo
        key=0
        call cdgemp(b,a,4,4,1.d-99,key)
        if(key.eq.0)then
          print *,'warning in qsve2am: anormal exit from cdgemp!'
          return
        endif
        do j=1,nj
          do i=1,4
            c(i,j)=(0.d0,0.d0)
            do m=1,4
              c(i,j)=c(i,j)+a(i,m)*y(m,j)
            enddo
          enddo
        enddo
      endif
c
      return
      end
