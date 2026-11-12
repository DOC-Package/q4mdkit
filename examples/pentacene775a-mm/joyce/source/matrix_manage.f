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
      Subroutine Trsml (A,C,X,N,M,zover)
c----------------------------------------------------------------------c
 
c   similitude transformation of the real symmetric a matrix :
c        A(out,dim=N) = (C+) * A(inp,dim=M) * C     if zover=T
c        X(out,dim=N) = (C+) * A(inp,dim=M) * C     if zover=F
c   the transformation requires a triangular work array
c   rather than a rectangular array as in common routines.
 
c  A = triangular symmetric matrix stored in dictionary order :
c      i.e the i,j element is found in A(ij) where
c      ij=i+j*(j-1)/2 if i<=j
c      at input has dim=M
c      if(zover=T) at output it is the transformed matrix of dim=N
c      if(zover=F) at output it is unchanged
c  C = N vectors of dimension M stored by column 
c  X = work matrix of dimension N*(N+1)/2
c      at output it is the transformed matrix of dim=N
c  N = dimension of the transformed output matrix
c  M = dimension of the input matrix
c  zover =  T  at return A=X = transformed matrix (A input is overwritten)
c           F  at return X = transformed matrix , A = unchanged
 
c     the subroutine employes the vector Y (dim>=M) whose
c     dimension is set by the parameter max which must be
c     greater than the input M value.
 
c               written by ivo cacelli  -  Pisa 
c               first version        30.01.1985
c               last update          14.10.1996
 
      implicit real*8(a-h,o-y),logical(z)
      parameter (max=5200)
      dimension A(*),C(*),X(*),Y(max) 
 
      if(N.le.0)             goto 100
      if(M.le.0.or.M.gt.max) goto 100 
      ij=0
      jj=-M
      do 70 j=1,N
      jj=jj+M
      ll=0
      do 50 l=1,M
      ll=ll+l-1
      ff=0.d0
      do 40 k=1,M
      if(k-l) 20,20,25
   20 kl=ll+k
      goto 30
   25 kl=kl+k-1
   30 ff=ff+C(jj+k)*A(kl)
   40 continue
   50 Y(l)=ff

      ii=-M
      do 70 i=1,j
      ii=ii+M
      ff=0.d0
      do 60 l=1,M
   60 ff=ff+Y(l)*C(ii+l)
      ij=ij+1
   70 X(ij)=ff

      if(zover) then
        do j=1,(N*N+N)/2
          A(j)=X(j)
        enddo
      endif

      return
  100 write(6,200) N,M,max
  200 format(/,/5x,'subr. Trsml : dimension',2i5,'  out-of-range'
     &      ,'  0-',i3/,/)
      stop 'Error in Trsml: out-of-range dimension'
      end

c----------------------------------------------------------------------c
      Subroutine Trsml2 (A,C,X,N,M,ndim,zover)
c----------------------------------------------------------------------c
 
c   similitude transformation of the real A matrix :
c        X(out,dim=N) = (C+) * A(inp,dim=M) * C
 
c  A    = input 2-dim array matrix  (M,M)
c  C    = input vectors             (M,N)
c  X    = output matrix             (N,N)
c  N    = dimension of the transformed output matrix
c  M    = dimension of the input matrix
c  ndim = first array dimension in the caller
c  zover =  T  at return A=X = transformed matrix (A input is overwritten)
c           F  at return X = transformed matrix , A = unchanged
 
      implicit real*8(a-h,o-z)
      parameter (max=5200)
      dimension A(ndim,*),C(ndim,*),X(ndim,*),Y(max) 
      logical zover
 
      if(N.le.0)             goto 100
      if(M.le.0.or.M.gt.max) goto 100 

      do 70 j=1,N

      do 50 l=1,M
      ff=0.d0
      do 40 k=1,M
      ff=ff+C(k,j)*A(k,l)
   40 continue
   50 Y(l)=ff

      do 70 i=1,N
      ff=0.d0
      do 60 l=1,M
   60 ff=ff+Y(l)*C(l,i)
   70 X(i,j)=ff

      if(zover) then
        do i=1,N
          do j=1,N
            A(i,j)=X(i,j)
          enddo
        enddo
      endif

      return
  100 write(6,200) N,M,max
  200 format(/,/5x,'subr. Trsml2 : dimension',2i5,'  out-of-range'
     &      ,'  0 -',i5/,/)
      stop
      end

c----------------------------------------------------------------------c
      Subroutine MatMt2 (A,C,X,N,M,ndim,zover)
c----------------------------------------------------------------------c
 
c   Perform the following transformation on the A matrix:
c        X(out,dim=N,M) = (C+) * A(inp,dim=N,M)
 
c  A    = input 2-dim array matrix  (N,M)
c  C    = input vectors             (N,N)   square !!!
c  X    = output matrix             (N,M)
c  N    = dimension of the transforming matrix C (# row of A,X)
c  M    = # columns (vectors) to br trasnf. of A,X
c  ndim = first array dimension in the caller
c  zover =  T  at return A=X = transformed matrix (A input is overwritten)
c           F  at return X = transformed matrix , A = unchanged
 
      implicit real*8(a-h,o-z)
      dimension A(ndim,*),C(ndim,*),X(ndim,*)
      logical zover
 
      if(N.le.0 .or. M.le.0) goto 100
      do 50 j=1,N
      do 50 l=1,M
      ff=0.d0
      do 40 k=1,N
      ff=ff+C(k,j)*A(k,l)
   40 continue
   50 X(j,l)=ff

      if(zover) then
        do i=1,N
          do j=1,M
            A(i,j)=X(i,j)
          enddo
        enddo
      endif

      return
  100 write(6,200) N,M
  200 format(/,/5x,'subr. Matvc2 : dimension',2i5,'  wrong !!',//)
      stop
      end

c----------------------------------------------------------------------c
      Subroutine MatM2 (A,C,X,N,M,ndim,zover)
c----------------------------------------------------------------------c
 
c   Perform the following transformation on the A matrix:
c        X(out,dim=N,M) = A(inp,dim=N,M) * C(inp,dim=M,M)
 
c  A    = input 2-dim array matrix  (N,M)
c  C    = input vectors             (M,M)   square !!!
c  X    = output matrix             (N,M)
c  ndim = first array dimension in the caller
c  zover =  T  at return A=X = transformed matrix (A input is overwritten)
c           F  at return X = transformed matrix , A = unchanged
 
      implicit real*8(a-h,o-z)
      dimension A(ndim,*),C(ndim,*),X(ndim,*)
      logical zover
 
      if(N.le.0 .or. M.le.0) goto 100
      do 50 j=1,N
      do 50 l=1,M
      ff=0.d0
      do 40 k=1,M
      ff=ff+A(j,k)*C(k,l)
   40 continue
   50 X(j,l)=ff

      if(zover) then
        do i=1,N
          do j=1,M
            A(i,j)=X(i,j)
          enddo
        enddo
      endif

      return
  100 write(6,200) N,M
  200 format(/,/5x,'subr. CvMat2 : dimension',2i5,'  wrong !!',//)
      stop
      end
