c----------------------------------------------------------------------c
c                                                                      c
c                                                                      c
c                  PROGRAM       J O Y C E                             c
c                                                                      c
c                                                                      c
c                                                                      c
c The energy and its first and second geometrical derivatives obtained c
c by DFT calculations (performed e.g. with the Gaussian package) for a c
c number of conformations of a single molecule are used to parametrize c
c intramolecular force fields, suitable for computer simulations.      c
c                                                                      c
c                                                                      c
c                                                                      c 
c   Copyright 2007-2025 G. Prampolini, J. Cerezo, S. Giannini,         c
c   N. De Mitri and I. Cacelli
c                                                                      c 
c   Joyce is free software: you can redistribute it and/or modify      c
c   it under the terms of the GNU General Public License as            c
c   published by  the Free Software Foundation, either version 3 of    c
c   the License, or (at your option) any later version.                c
c                                                                      c
c   Joyce is distributed in the hope that it will be useful,           c
c   but WITHOUT ANY WARRANTY; without even the implied warranty of     c 
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      c
c   GNU General Public License for more details.                       c
c                                                                      c
c   You should have received a copy of the GNU General Public License  c
c   along with Joyce. If not, see <http://www.gnu.org/licenses/>.      c
c                                                                      c
c   For information about new develepments and debugs you can contact: c
c                                                                      c
c   javier.cerezo@uam.es                                               c
c   giacomo.prampolini@pi.iccom.cnr.it                                 c
c   samuele.giannin@cnr.it                                             c
c                                                                      c
c----------------------------------------------------------------------c

c----------------------------------------------------------------------c
      Subroutine Eigen (a,r,n)
c----------------------------------------------------------------------c
 
 
c     diagonalization of symmetric matrices with Jacobi method
c     rearrange the eigenvalues (and e.vec) in increasing order
 
      implicit real*8 (a-h,o-z)
      dimension a(*),r(*)
      range=1.0d-12
      nn=n*n
      do 12 i=1,nn
   12 r(i)=0.d0
      do 20 i=1,n
   20 r(i*n-n+i)=1.d0
      anorm=0.d0
      do 35 i=1,n
      do 34 j=1,i
      ia=j+(i*i-i)/2
   34 anorm=anorm+a(ia)*a(ia)
   35 anorm=anorm-a(ia)*a(ia)
      if (anorm) 165,165,40
   40 anorm=2.d0*dsqrt(anorm)
      anrmx=anorm*range/dble(n)
      ind=0
      thr=anorm
   45 thr=thr/dble(n)
   50 l=1
   55 m=l+1
   60 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
   62 if (dabs(a(lm))-thr) 130,65,65
   65 ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if (x.lt.0.d0) y=-y
      sinx=y/dsqrt(2.0d0*(1.0d0+dsqrt(1.0d0-y*y)))
      sinx2=sinx*sinx
      cosx2=1.0d0-sinx2
      cosx=dsqrt(cosx2)
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if (i.eq.l) goto 120
      if (i-m) 85,120,90
   85 im=i+mq
      goto 95
   90 im=m+iq
   95 if (i-l) 100,105,105
  100 il=i+lq
      goto 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  120 ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
  125 r(ilr)=x
      x=2.0d0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
  130 if (m.eq.n) goto 140
      m=m+1
      goto 60
  140 if (l-(n-1)) 145,150,145
  145 l=l+1
      goto 55
  150 if (ind.ne.1) goto 160
      ind=0
      goto 50
  160 if (thr.gt.anrmx) goto 45
  165 do 166 i=1,n
  166 a(i)=a((i*i+i)/2)

c     rearranges the eigenvalues

      iq=-n
      do 185 i=1,n
      iq=iq+n
      jq=n*(i-2)
      do 185 j=i,n
      jq=jq+n
      if(a(j)-a(i)) 170,185,185
  170 x=a(i)
      a(i)=a(j)
      a(j)=x
      do 180 k=1,n
      ilr=iq+k
      imr=jq+k
      x=r(ilr)
      r(ilr)=r(imr)
  180 r(imr)=x
  185 continue
      return
      end



c----------------------------------------------------------------------c
      Subroutine Eigen0 (a,r,n)
c----------------------------------------------------------------------c
 
 
c     diagonalization of symmetric matrices with Jacobi method
c     with no rearrangement
 
      implicit real*8 (a-h,o-z)
      dimension a(*),r(*)
      range=1.0d-12
      nn=n*n
      do 12 i=1,nn
   12 r(i)=0.d0
      do 20 i=1,n
   20 r(i*n-n+i)=1.d0
      anorm=0.d0
      do 35 i=1,n
      do 34 j=1,i
      ia=j+(i*i-i)/2
   34 anorm=anorm+a(ia)*a(ia)
   35 anorm=anorm-a(ia)*a(ia)
      if (anorm) 165,165,40
   40 anorm=2.d0*dsqrt(anorm)
      anrmx=anorm*range/dble(n)
      ind=0
      thr=anorm
   45 thr=thr/dble(n)
   50 l=1
   55 m=l+1
   60 mq=(m*m-m)/2
      lq=(l*l-l)/2
      lm=l+mq
   62 if (dabs(a(lm))-thr) 130,65,65
   65 ind=1
      ll=l+lq
      mm=m+mq
      x=0.5d0*(a(ll)-a(mm))
      y=-a(lm)/dsqrt(a(lm)*a(lm)+x*x)
      if (x.lt.0.d0) y=-y
      sinx=y/dsqrt(2.0d0*(1.0d0+dsqrt(1.0d0-y*y)))
      sinx2=sinx*sinx
      cosx2=1.0d0-sinx2
      cosx=dsqrt(cosx2)
      sincs=sinx*cosx
      ilq=n*(l-1)
      imq=n*(m-1)
      do 125 i=1,n
      iq=(i*i-i)/2
      if (i.eq.l) goto 120
      if (i-m) 85,120,90
   85 im=i+mq
      goto 95
   90 im=m+iq
   95 if (i-l) 100,105,105
  100 il=i+lq
      goto 110
  105 il=l+iq
  110 x=a(il)*cosx-a(im)*sinx
      a(im)=a(il)*sinx+a(im)*cosx
      a(il)=x
  120 ilr=ilq+i
      imr=imq+i
      x=r(ilr)*cosx-r(imr)*sinx
      r(imr)=r(ilr)*sinx+r(imr)*cosx
  125 r(ilr)=x
      x=2.0d0*a(lm)*sincs
      y=a(ll)*cosx2+a(mm)*sinx2-x
      x=a(ll)*sinx2+a(mm)*cosx2+x
      a(lm)=(a(ll)-a(mm))*sincs+a(lm)*(cosx2-sinx2)
      a(ll)=y
      a(mm)=x
  130 if (m.eq.n) goto 140
      m=m+1
      goto 60
  140 if (l-(n-1)) 145,150,145
  145 l=l+1
      goto 55
  150 if (ind.ne.1) goto 160
      ind=0
      goto 50
  160 if (thr.gt.anrmx) goto 45
  165 do 166 i=1,n
  166 a(i)=a((i*i+i)/2)
 
      return
      end



c----------------------------------------------------------------------c
      Subroutine Eigen2 (a,r,n,ndim)
c----------------------------------------------------------------------c
 
 
c  diagonalization of symmetric matrices with Jacobi method
c  rearrange the eigenvalues (and e.vec) in increasing order
c  version for 2-dim arrays
 
      implicit real*8 (a-h,o-z)
      dimension a(ndim,*),r(ndim,*)
      range=1.0d-12
 
c     symmetrize the input matrix (no check !)
 
      tol=range*100.d0
      do j=1,n
        do i=1,j
          diff=A(i,j)-A(j,i)
          if(dabs(diff).gt.tol) goto 200
          A(i,j)=0.5d0*(A(i,j)+A(j,i))
          A(j,i)=A(i,j)
        enddo
      enddo
 
      do 14 i=1,n
      do 12 j=1,n
   12 r(i,j)=0.d0
   14 r(i,i)=1.d0
      anorm=0.d0
      do 34 i=2,n
      do 34 j=1,i-1
   34 anorm=anorm+a(i,j)*a(i,j)
 
      if (anorm) 165,165,40
   40 anorm=2.d0*dsqrt(anorm)
      anrmx=anorm*range/dble(n)
      ind=0
      thr=anorm
   45 thr=thr/dble(n)
   50 l=1
   55 m=l+1
   60 continue
   62 if (dabs(a(l,m))-thr) 130,65,65
   65 ind=1
      x=0.5d0*(a(l,l)-a(m,m))
      y=-a(l,m)/dsqrt(a(l,m)*a(l,m)+x*x)
      if (x.lt.0.d0) y=-y
      sinx=y/dsqrt(2.0d0*(1.0d0+dsqrt(1.0d0-y*y)))
      sinx2=sinx*sinx
      cosx2=1.0d0-sinx2
      cosx=dsqrt(cosx2)
      sincs=sinx*cosx

c     rotation for the l,m elements

      do 125 i=1,n
      if (i.eq.l) goto 120
      if (i-m) 85,120,85
   85 continue
      x=a(i,l)*cosx-a(i,m)*sinx
      a(i,m)=a(i,l)*sinx+a(i,m)*cosx
      a(i,l)=x
      a(l,i)=x
      a(m,i)=a(i,m)
  120 x=r(i,l)*cosx-r(i,m)*sinx
      r(i,m)=r(i,l)*sinx+r(i,m)*cosx
  125 r(i,l)=x

c     special for l,m

      x=2.0d0*a(l,m)*sincs
      y=a(l,l)*cosx2+a(m,m)*sinx2-x
      x=a(l,l)*sinx2+a(m,m)*cosx2+x
      a(l,m)=(a(l,l)-a(m,m))*sincs+a(l,m)*(cosx2-sinx2)
      a(m,l)=a(l,m)
      a(l,l)=y
      a(m,m)=x

  130 if (m.eq.n) goto 140
      m=m+1
      goto 60
  140 if (l-(n-1)) 145,150,145
  145 l=l+1
      goto 55
  150 if (ind.ne.1) goto 160
      ind=0
      goto 50
  160 if (thr.gt.anrmx) goto 45
  165 continue
 
c     rearrange the eigenvalues
 
      do 185 i=1,n
      do 185 j=i,n
      if(a(j,j)-a(i,i)) 170,185,185
  170 x=a(i,i)
      a(i,i)=a(j,j)
      a(j,j)=x
      do 180 k=1,n
      x=r(k,i)
      r(k,i)=r(k,j)
  180 r(k,j)=x
  185 continue
      return
 
  200 write(6,202) i,j,A(i,j),A(j,i),diff,tol
      stop
  202 format(//5x,'ERROR in EIGEN2: input matrix is not symmetric'
     &     ,/10x,'i j .................',2i5
     &     ,/10x,'A(i,j) ..............',f20.10
     &     ,/10x,'A(j,i) ..............',f20.10
     &     ,/10x,'A(i,j)-A(j,i) .......',f20.10
     &     ,/10x,'tolerance ...........',f20.10,///)
      end



c----------------------------------------------------------------------c
      Subroutine Rearr2 (a,r,n,ndim,s)
c----------------------------------------------------------------------c
 
c     rearrange eigenvalues and eigenvectors 
c     in increasing/decrasing order
 
      implicit real*8 (a-h,o-z)
      dimension a(ndim,*),r(ndim,*)
      character*1 s
 
      if(s.eq.'d' .or. s.eq.'D' .or. s.eq.'-') goto 15
 
c     increasing order
 
      do 10 i=1,n-1
      do 10 j=i+1,n
      if(a(j,j).gt.a(i,i)) goto 10
      x=a(i,i)
      a(i,i)=a(j,j)
      a(j,j)=x
      do k=1,n
        x=r(k,i)
        r(k,i)=r(k,j)
        r(k,j)=x
      enddo
   10 continue
      return
 
c     decreasing order
 
  15  do 20 i=1,n-1
      do 20 j=i+1,n
      if(a(j,j).lt.a(i,i)) goto 20
      x=a(i,i)
      a(i,i)=a(j,j)
      a(j,j)=x
      do k=1,n
        x=r(k,i)
        r(k,i)=r(k,j)
        r(k,j)=x
      enddo
   20 continue
      return
 
      end


c----------------------------------------------------------------------c
      Subroutine Lowdin (S,avt,thr0,nfb,iop,zprt)
c----------------------------------------------------------------------c
 
 
c                         computes
c      the S(1/2) symmetric matrix such that S=S(1/2)*S(1/2)
c   or the S(-1/2) symmetric matrix such that S(-1)=S(-1/2)*S(-1/2)
c   which are useful in the Lowdin orthogonalization
 
c  S     inp   metric matrix (triangular) - destroyed
c        out   S(1/2) or S(-1/2) depending on iop
c  avt   out   metric eigenvector matrix by column (square)
c  thr0  inp   threshold for selection of the metric eigenvalues
c  nfb   inp   problem dimension
c  iop   inp   >=0 compute S(1/2)   <0 compute S(-1/2)  
c  zprt  inp   (logical)  t=print some results  f=noprint
 
 
      implicit real*8(a-h,o-z)
      parameter (maxdim=3600)
      dimension avt(*),S(*),y(maxdim)
      logical zprt,zpr
 
      if(nfb.gt.maxdim) then
        write(6,2) nfb,maxdim
        stop
      endif
    2 format(/6x,'Error in subr. LOWDIN: too large dimens.',i6
     &      ,/6x,'                      allowed max =',i5,/)
      if(thr0.lt.0.d0) then
        write(6,'(/6x,"Error in LOWDIN: thr0<0",1pd12.3,/)') thr0
        stop
      endif
      if(nfb.le.1) then
        write(6,'(/6x,"Error in LOWDIN: nfb < 1",/)')
        stop
      endif
 
c     diagonalization of the metric matrix
 
      ndim=nfb
      Call Eigen (S,avt,nfb)
      smax=S(nfb)
      thre=thr0*Dmin1(smax,1.0d0)
      zero=-thre
 
      n1=0
      n0=0
      do j=1,nfb
        y(j)=S(j)
        if(y(j).lt.thr0) n1=n1+1
        if(y(j).lt.zero) n0=n0+1
      enddo
      nuve=nfb-n1
      if(n0.gt.0) zpr=.true.
 
      if(zpr) then
        nwr=n1+5
        if(nwr.gt.nfb) nwr=nfb
        write(6,600)  nfb,n1,nuve,thr0
        if(n1.gt.0) write(6,605)  (y(k),k=1,n1)
        write(6,606)  n1+1,(y(k),k=n1+1,nwr)
        iwr=n1   
        jwr=n1+1
        if(n1.eq.0) iwr=jwr
        do j=iwr,jwr
          jj2=j*nfb
          jj1=jj2-nfb+1
          write(6,610)  j,y(j),(avt(k),k=jj1,jj2)
        enddo
      endif
 
      if(n0.gt.0) then
        do j=1,n0
          jj2=j*nfb
          jj1=jj2-nfb+1
          if (levprt .gt. 0) write(6,610)  j,S(j),(avt(k),k=jj1,jj2)
        enddo
        write(6,615)
        stop
      endif
 
c    computes the eigen**1/2 or eigen**-1/2
 
      do i=1,nfb
        if(i.gt.n1) then
          if(iop.ge.0) y(i)=dsqrt(y(i))
          if(iop.lt.0) y(i)=1.0d0/dsqrt(y(i))
        endif
      enddo
 
c   computes the matrix:  s=avt*Y*(avt+)
c   with Y = metric eigenvalues **1/2 or **(-1/2)
 
      mk=0
      do k=1,nfb
        do m=1,k
          ki=k
          mi=m
          aaa=0.0d0
          do i=1,nfb
            if(i.gt.n1) then
              aaa=aaa+avt(ki)*avt(mi)*y(i)
            endif
            ki=ki+nfb
            mi=mi+nfb
          enddo
          mk=mk+1
          S(mk)=aaa
        enddo
      enddo
 
 
      return
  600 format(/,/7x,'----------------------------------------------'
     &        ,/7x,'-    Symmetric (Lowdin) Orthogonalization    -'
     &        ,/7x,'-   trhough diagonalization of the Metric    -'
     &        ,/7x,'----------------------------------------------'
     &   ,/7x,'starting dimension of the eigenv. problem',i5
     &   ,/7x,'n. of vanishing metric eigenvalues .. ...',i5
     &   ,/7x,'n. of very small metric eigenvalues .....',i5
     &   ,/7x,'threshold for the metric eigenvalues ....',1pd11.2)
  605 format(' discarded metric eigenvalues'/(3x,1p5d11.3))
  606 format(' first kept metric eigenvalues from',i6/(3x,1p5d11.3))
  610 format(/' metric eigensolution n. ',i3,g15.7/(1x,10f7.2))
  615 format(/,/9x,'Irreparable Error in subr. LOWDIN:'
     &        ,/9x,'The metric matrix is not positive definite'
     &        ,/9x,'I"m not able to continue, STOP',/)



      end



c----------------------------------------------------------------------c
      Subroutine Inv_SpecDec (S,avt,thr0,nfb,nzeri,zprt)
c----------------------------------------------------------------------c


c       invert the S matrix, by spectral decomposition
 
c  S     inp   metric matrix (triangular) - destroyed
c        out   S**(-1)
c  avt   out   metric eigenvector matrix by column (square)
c  thr0  inp   threshold for selection of the metric eigenvalues
c              Eigenvectors with eigenvalues > thr0*Max(eigenval) are used for inversion
c  nfb   inp   problem dimension
c  nzeri out   number of eigenvectors excluded in the inversion. They correspond
c              to the first nzeri eigenvectors in the avt array.
c  zprt  inp   (logical)  t=print some results  f=noprint
 
 
      implicit real*8(a-h,o-z)
      parameter (maxdim=3600)
      dimension avt(*),S(*),y(maxdim)
      logical zprt,zpr
      zpr=zprt
 
      if(nfb.gt.maxdim) then
        write(6,2) nfb,maxdim
        stop
      endif
    2 format(/6x,'Error in subr. Inv_Eig: too large dimens.',i6
     &      ,/6x,'                      allowed max =',i5,/)
      if(thr0.lt.0.d0) then
        write(6,'(/6x,"Error in Inv_Eig: thr0<0",1pd12.3,/)') thr0
        stop
      endif
      if(nfb.le.1) then
        write(6,'(/6x,"Error in Inv_Eig: nfb < 1",/)')
        stop
      endif
 
c     diagonalization of the metric matrix
 
      ndim=nfb
      Call Eigen (S,avt,nfb)
      smax=S(nfb)
      thre=thr0*Dabs(smax)           
      zero=-thre
 
      n1=0
      n0=0
      do j=1,nfb
        y(j)=S(j)
        if(y(j).lt.thre) n1=n1+1
        if(y(j).lt.zero) n0=n0+1
      enddo
      nuve=nfb-n1
      if(n0.gt.0) zpr=.true.
 
      if(zpr) then
        nwr=n1+5
        if(nwr.gt.nfb) nwr=nfb
        if (levprt .gt. 0 ) write(6,600)  nfb,n1,nuve,thr0,smax,thre
        if(n1.gt.0 .and. levprt .gt. 0) write(6,605)  (y(k),k=1,n1)
        write(6,606)  n1+1,(y(k),k=n1+1,nwr)
        iwr=n1   
        jwr=n1+1
        if(n1.eq.0) iwr=jwr
        do j=iwr,jwr
          jj2=j*nfb
          jj1=jj2-nfb+1
          write(6,610)  j,y(j),(avt(k),k=jj1,jj2)
        enddo
      endif
 
      if(n0.gt.0) then
        do j=1,n0
          jj2=j*nfb
          jj1=jj2-nfb+1
          write(6,610)  j,S(j),(avt(k),k=jj1,jj2)
        enddo
        write(6,615)
        stop
      endif
 
c     computes the eigen**(-1)
 
      do i=1,nfb
        if(i.gt.n1) then
          y(i)=1.0d0/y(i)
        endif
      enddo
 
c     computes the matrix:  S=avt*Y*(avt+)
c     with Y = inverse metric eigenvalues
 
      mk=0
      do k=1,nfb
        do m=1,k
          ki=k
          mi=m
          aaa=0.0d0
          do i=1,nfb
            if(i.gt.n1) then
              aaa=aaa+avt(ki)*avt(mi)*y(i)
            endif
            ki=ki+nfb
            mi=mi+nfb
          enddo
          mk=mk+1
          S(mk)=aaa
        enddo
      enddo
      nzeri=n1
 
 
      return
  600 format(/,/10x,'----------------------------------------------'
     &        ,/10x,'-     Inversion of a symmetric matrix        -'
     &        ,/10x,'-     trhough spectral decomposition         -'
     &        ,/10x,'----------------------------------------------'
     &   ,/10x,'starting dimension of the eigenv. problem',i5
     &   ,/10x,'n. of vanishing metric eigenvalues .. ...',i5
     &   ,/10x,'n. of taken eigenvalues .................',i5
     &   ,/10x,'factor to mult max(eig) to get thr ......',1pd11.2
     &   ,/10x,'max eigenvalue ..........................',1pd11.2
     &   ,/10x,'threshold used for eigenvalues ..........',1pd11.2)
  605 format(10x,' discarded metric eigenvalues'/(3x,1p5d11.3))
  606 format(10x,'first kept metric eigenvalues from',i6/(10x,1p5d11.3))
  610 format(/10x,'metric eigensolution n. ',i3,g15.7/(10x,10f7.2))
  615 format(/,/10x,'Irreparable Error in subr. Inv_Eig:'
     &        ,/10x,'The metric matrix is not positive definite'
     &        ,/10x,'I"m not able to continue, STOP',/)



      end
