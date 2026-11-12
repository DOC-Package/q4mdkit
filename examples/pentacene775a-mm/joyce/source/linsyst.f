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
      Subroutine Lin_Eq_SVD (A,B,X,n,neq,thr0,Sens,ndim,ierr)
c----------------------------------------------------------------------c

c     solve the linear equation system with simmetric matrix
c     using the Singular Value Decomposition method adapted
c     for symmetric matrices i.e. Spectral Decomposition
c
c                     A(n,n) * X(n,m) = B(n,m)
c
c     for m values of B (columns)

c     A       i  characteristic matrix of the system
c             o  covariance matrix, = 1/A
c     X       o  solutions (unknown quantities)
c     B       i  known vectors
c     n       i  problem dimension
c     neq     i  n. of equations to be solved
c     thr0    i  threshold for zero eigenvalues (in the inversion)
c     Sens    o  sensitivity index for each unknown quantity
c                Sens(i)=1 the unknown X(i) is non redundant, its value
c                          is well determined, unique and active
c                Sens(i)=0 the unknown X(i) is redundant or not sensitive
c                          for the current equation system
c     ndim    i  dimensions of A,B,X in the caller
c     ierr    o  error parameter (0=OK)

      implicit real*8(a-h,o-z)
      parameter (nmax=2000, ntri=nmax*(nmax+1)/2, nqua=nmax*nmax)
      dimension A(ndim,*),B(ndim,*),X(ndim,*),Sens(*)
      dimension F(ntri),C(nqua),ndsym(4),bnew(nmax)
      character title*60
      logical zprt
c
      if(n.gt.nmax) then
        write(6,10) n
        stop
      endif
   10 format(//' ERROR in Lin_Eq_SVD: increase nmax to',i7,/)
   12 format(//' ERROR in Lin_Eq_SVD: A matr not symmetric'
     &       ,/' element',2i5,f20.10
     &       ,/' element',2i5,f20.10,/)
c
      amax=0.d0
      do i=1,n
        do j=1,n
          if(dabs(A(i,j)).gt.amax) amax=dabs(A(i,j))
        enddo
      enddo
c
      tolsym=1.d-07*amax
      tolzer=tolsym*1.d-03
      ij=0
      do j=1,n
        do i=1,j
          errsym=0.5d0*(A(i,j)-A(j,i))
          summat=0.5d0*(A(i,j)+A(j,i))
          sumabs=0.5d0*(dabs(A(i,j))+dabs(A(j,i)))
          if(sumabs.gt.tolzer) then
            if(dabs(errsym/summat).gt.tolsym) then
              write(6,12) i,j,A(i,j),j,i,A(j,i)
              stop
            endif
          endif
          ij=ij+1
          F(ij)=summat
          A(i,j)=summat
          A(j,i)=summat
        enddo
      enddo

c     invert the A matrix by the spectral decomposition

      zprt=.true.
      Call Inv_SpecDec (F,C,thr0,n,nzeri,zprt)
      ierr=1
      if(n.eq.nzeri) return
      ierr=0
      if (levprt .gt. 0) 
     &     write(6,*) ' n. of low eigenvalues in the inversion',nzeri

c     test on the inversion

      errijm=0.d0
      do j=1,n
        do i=1,j
          xij=0
          do k=1,n
            xij=xij+A(i,k)*F(ijtri(j,k))
          enddo
          ref=0.0d0
          if(i.eq.j) ref=1.0d0
          errij=dabs(xij-ref)
          if(errij.gt.errijm) then
            imi=i
            jmi=j
            errijm=errij
          endif
        enddo
      enddo

c     compute the X

      do im=1,neq
        do i=1,n
          xi=0.0d0
          do j=1,n
            xi=xi+F(ijtri(i,j))*B(j,im)
          enddo
          X(i,im)=xi
        enddo
      enddo

c     test on the accuracy of the solutions

      errmax=0.0d0
      do im=1,neq
        psnb=0.0d0
        psbb=0.0d0
        ermx=0.0d0
        do i=1,n
          bi=0.0d0
          do j=1,n
            bi=bi+A(i,j)*X(j,im)
          enddo
          bnew(i)=bi
          diff=bnew(i)-B(i,im)
          psnb=psnb+diff**2
          psbb=psbb+B(i,im)**2
          if(dabs(diff).gt.ermx) then
            ii=i
            ermx=dabs(diff)
          endif
        enddo
        errmed=dsqrt(psnb/psbb)
        if (zlevprt .gt. 0 ) write(6,72) ermx,ii,errmed,dsqrt(psbb)
        if(errmed.gt.errmax) errmax=errmed
      enddo
      write(6,74) errijm,imi,jmi
      if(errmax.gt.0.01d0) write(6,76)
   72 format(/10x,'ACCURACY OF THE LINEAR EQUATION SOLUTIONS'
     &      ,/10x,'Maximum Error in |A*x-B|    =  ',f13.6,'  on elem',i4
     &      ,/10x,'  |A*x-B|/|B|               =  ',f13.6
     &      ,/10x,'  |B|                       =  ',f13.6)
   74 format( 10x,'Max err in the matrix inversion',g13.3
     &      ,/10x,'for the matrix element           ',2i5)
   76 format(8(/16x,'W A R N I N G :  ERROR IN LIN EQ SYST  VERY HIGH'))

c     computes the sensitivity index

      kbeg=nzeri*n
      do i=1,n
        ki=kbeg+i
        aaa=0.0d0
        do k=nzeri+1,n
          aaa=aaa+C(ki)*C(ki)
          ki=ki+n
        enddo
        Sens(i)=aaa
      enddo

c     copy the inverse matrix on the A array

      ij=0
      do j=1,n
        do i=1,j
          ij=ij+1
          A(i,j)=F(ij)
          A(j,i)=F(ij)
        enddo
      enddo
c
      return
      end

c----------------------------------------------------------------------c
      Function ijtri (i,j)
c----------------------------------------------------------------------c

      if(i.gt.j) then
        ijtri=i*(i-1)/2+j
      else
        ijtri=j*(j-1)/2+i
      endif
      return
      end

c----------------------------------------------------------------------c
      Subroutine Lin_Eq (A,B,X,n,m,ndim,ierr)
c----------------------------------------------------------------------c

c     solves linear eq. system  A(.,.)*X(.,j)=B(.,j)
c                               j=1,2,..,m .

c     method gauss&jordan with pivot search
c     written by Roberto Moccia

      implicit real*8(a-h,o-z)
      parameter (nmax=2000)
      dimension A(ndim,*),B(ndim,*),x(ndim,*),is(nmax),js(nmax)
      if(n.gt.nmax) stop 'ERROR in Lin_Eq: increase nmax'
      if(n.le.1)    stop 'ERROR in Lin_Eq: n<2'
      if(m.lt.1)    stop 'ERROR in Lin_Eq: m<1'
      sog=1.d-12
c
      aso=1.d0
      ierr=0
      ls=1
      do 50 i=1,n-1
      aa=dabs(A(i,i))
      ind=i
      do 10 j=i+1,n
      bb=dabs(A(i,j))
      if(bb.lt.aa) goto 10
      ind=j
      aa=bb
   10 continue
      if(ind.eq.i) goto 25
      do 20 j=1,n
      aa=A(j,i)
      A(j,i)=A(j,ind)
   20 A(j,ind)=aa
      is(ls)=i
      js(ls)=ind
      ls=ls+1
   25 aso=aso*A(i,i)
      if(dabs(aso).lt.sog) goto 90
      do 30 j=i+1,n
      aa=A(j,i)/A(i,i)
      do lm=1,m
      B(j,lm)=B(j,lm)-B(i,lm)*aa
      enddo
      do 30 l=i+1,n
   30 A(j,l)=A(j,l)-A(i,l)*aa
   50 continue
c
      ls=ls-1
      do lm=1,m
      x(n,lm)=B(n,lm)/A(n,n)
      do 70 l=1,n-1
      aa=0.d0
      do 65 j=1,l
   65 aa=aa+A(n-l,n-l+j)*x(n-l+j,lm)
   70 x(n-l,lm)=(B(n-l,lm)-aa)/A(n-l,n-l)
      enddo
c
      if(ls.eq.0) goto 80
      do lm=1,m
      do 75 l=1,ls
      i=is(ls+1-l)
      j=js(ls+1-l)
      aa=x(i,lm)
      x(i,lm)=x(j,lm)
   75 x(j,lm)=aa
      enddo
   80 return
   90 ierr=1
      write(6,92) dabs(aso),sog
   92 format(//8x,'ERROR in Lin_Eq: determ too small'
     &       ,/8x,'determ , toler ',1p2d13.3,//)
      return
c     stop
      end
