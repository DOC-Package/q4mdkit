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
      Subroutine Invert (c,n,ndim,ier)
c----------------------------------------------------------------------c

c     matrix inversion (written by R. Ambrosetti)

      implicit real*8 (a-h,o-z)
      parameter (max=250, errmax=1.d-11)
      dimension c(ndim,*),al(max,max),bl(max),cl(max),lz(max)
c
      if(n.gt.max) goto 200
      do 10 i=1,n
      lz(i)=i
      do 10 j=1,n
      al(i,j)=c(i,j)
   10 continue
c
      do 20 i=1,n
      k=i
      yl=al(i,i)
      l=i-1
      lp=i+1
      if(n.lt.lp)  goto 14
      do 12 j=lp,n
      w=al(i,j)
      if(dabs(w).le.dabs(yl)) goto 12
      k=j
      yl=w
   12 continue
   14 do 16 j=1,n
      cl(j)=al(j,k)
      al(j,k)=al(j,i)
      al(j,i)=-cl(j)/yl
      al(i,j)=al(i,j)/yl
   16 bl(j)=al(i,j)
      al(i,i)=1.0/yl
      j=lz(i)
      lz(i)=lz(k)
      lz(k)=j
      do 19 k=1,n
      if (i.eq.k)  goto 19
      do 18 j=1,n
      if (i.eq.j)  goto 18
      al(k,j)=al(k,j)-bl(j)*cl(k)
   18 continue
   19 continue
   20 continue
c
      do 30 i=1,n
      if (i.eq.lz(i))  goto 30
      k=i+1
      do 28 j=k,n
      if (i.ne.lz(j))  goto 28
      m=lz(i)
      lz(i)=lz(j)
      lz(j)=m
      do  26  l=1,n
      cl(l)=al(i,l)
      al(i,l)=al(j,l)
   26 al(j,l)=cl(l)
   28 continue
   30 continue

c     checks the inversion error

      err=0.d0
      do 50 j=1,n
      do 48 k=1,n
      yl=0.d0
      do 46 l=1,n
   46 yl=yl+c(j,l)*al(l,k)
      if(j.eq.k) yl=yl-1.d0
   48 if(dabs(yl).gt.err) err=dabs(yl)
   50 continue
      if(dabs(err).gt.errmax) goto 220
c
      do 60 i=1,n
      do 60 j=1,n
   60 c(i,j)=al(i,j)
      ier=0
      return
c
  200 write(6,202) n,max
  202 format(/,/4x,'error: insuff. dimens. in subr. Invert'
     &        ,/4x,'dim. of the input matrix ....',i5
     &        ,/4x,'max allowed dimension .......',i5,/,/)
      stop ' ERROR in subr. Invert'
  220 write(6,222) err,errmax
  222 format(//4x,'error: loss of precision in subr. Invert'
     &       ,/4x,'evaluated inversion error ...',d12.3
     &       ,/4x,'max allowed error ...........',d12.3,/,/)
      ier=1
      return
      end

c----------------------------------------------------------------------c
      Subroutine Minv (A,N,D)
c----------------------------------------------------------------------c

c     invert a real non-symmetric matrix

      IMPLICIT DOUBLE PRECISION (a-h,o-z)
      parameter (max=2800)
      DIMENSION a(*),l(max),m(max)
      if(n.gt.max) stop 'ERROR: increase max in MINV'
      if(n.le.0  ) stop 'ERROR: N<1          in MINV'
c
      if(n.eq.1) go to 1999
c
c     search for largest element
c
      d=1.0
      nk=-n
      do 80 k=1,n
      nk=nk+n
      l(k)=k
      m(k)=k
      kk=nk+k
      biga=a(kk)
      do 20 j=k,n
      iz=n*(j-1)
      do 20 i=k,n
      ij=iz+i
   10 if(abs(biga)-abs(a(ij))) 15,20,20
   15 biga=a(ij)
      l(k)=i
      m(k)=j
   20 continue

c     interchange rows

      j=l(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      ki=ki+n
      hold=-a(ki)
      ji=ki-k+j
      a(ki)=a(ji)
   30 a(ji) =hold

c     interchange columns

   35 i=m(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      jk=nk+j
      ji=jp+j
      hold=-a(jk)
      a(jk)=a(ji)
   40 a(ji) =hold

c     divide column by minus pivot (value of pivot element is
c     contained in biga)

   45 if(biga) 48,46,48
   46 d=0.0
      write(6,'(//5x,"ERROR in DMINV: singular matrix",/)')
      return
   48 do 55 i=1,n
      IF(i-k) 50,55,50
   50 ik=nk+i
      a(ik)=a(ik)/(-biga)
   55 continue

c     reduce matrix

      do 65 i=1,n
      ik=nk+i
      hold=a(ik)
      ij=i-n
      do 65 j=1,n
      ij=ij+n
      IF(i-k) 60,65,60
   60 IF(j-k) 62,65,62
   62 kj=ij-i+k
      a(ij)=hold*a(kj)+a(ij)
   65 continue

c     divide row by pivot

      kj=k-n
      do 75 j=1,n
      kj=kj+n
      IF(j-k) 70,75,70
   70 a(kj)=a(kj)/biga
   75 continue

c     product of pivots

      d=d*biga

c     replace pivot by reciprocal

      a(kk)=1.d0/biga
   80 continue

c     final row and column interchange

      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=l(k)
      if(i-k) 120,120,108
  108 jq=n*(k-1)
      jr=n*(i-1)
      do 110 j=1,n
      jk=jq+j
      hold=a(jk)
      ji=jr+j
      a(jk)=-a(ji)
  110 a(ji) =hold
  120 j=m(k)
      if(j-k) 100,100,125
  125 ki=k-n
      do 130 i=1,n
      ki=ki+n
      hold=a(ki)
      ji=ki-k+j
      a(ki)=-a(ji)
  130 a(ji) =hold
      go to 100
  150 return
 1999 a(1)=1.d0/a(1)
      return
      end

c----------------------------------------------------------------------c
      Subroutine Minv2 (A,N,D,ndim)
c----------------------------------------------------------------------c

c     invert a real non-symmetric (two-index) matrix

      implicit double precision (a-h,o-z)
      parameter (max=2800)
      dimension A(ndim,*),indl(max),indm(max)
      if(n.gt.max) stop 'ERROR: increase max in MINV2'
      if(n.le.0  ) stop 'ERROR: N<1          in MINV2'
c
      if(n.eq.1) then
        a(1,1)=1.d0/a(1,1)
        return
      endif

c     search for largest element

      d=1.0
      do 80 k=1,n
      indl(k)=k
      indm(k)=k
      biga=a(k,k)
      do 20 j=k,n
      do 20 i=k,n
      if(abs(biga)-abs(a(i,j))) 15,20,20
   15 biga=a(i,j)
      indl(k)=i
      indm(k)=j
   20 continue

c     interchange rows

      j=indl(k)
      if(j-k) 35,35,25
   25 ki=k-n
      do 30 i=1,n
      hold=-a(k,i)
      a(k,i)=a(j,i)
   30 a(j,i)=hold

c     interchange columns

   35 i=indm(k)
      if(i-k) 45,45,38
   38 jp=n*(i-1)
      do 40 j=1,n
      hold=-a(j,k)
      a(j,k)=a(j,i)
   40 a(j,i)=hold

c      divide column by minus pivot (value of pivot element is
c      contained in biga)

   45 if(biga) 48,46,48
   46 d=0.0
      write(6,'(//5x,"ERROR in DMINV2: singular matrix",/)')
      return
   48 do 55 i=1,n
      if(i.ne.k) a(i,k)=a(i,k)/(-biga)
   55 continue

c     reduce matrix

      do 65 i=1,n
      hold=a(i,k)
      do 65 j=1,n
      if(i-k) 60,65,60
   60 if(j-k) 62,65,62
   62 continue
      a(i,j)=hold*a(k,j)+a(i,j)
   65 continue

c     divide row by pivot

      do 75 j=1,n
      if(j.ne.k) a(k,j)=a(k,j)/biga
   75 continue

c     product of pivots

      d=d*biga

c     replace pivot by reciprocal

      a(k,k)=1.d0/biga
   80 continue

c     final row and column interchange

      k=n
  100 k=(k-1)
      if(k) 150,150,105
  105 i=indl(k)
      if(i-k) 120,120,108
  108 do 110 j=1,n
      hold=a(j,k)
      a(j,k)=-a(j,i)
  110 a(j,i) =hold
  120 j=indm(k)
      if(j-k) 100,100,125
  125 do 130 i=1,n
      hold=a(k,i)
      a(k,i)=-a(j,i)
  130 a(j,i) =hold
      go to 100
  150 return
      end

c----------------------------------------------------------------------c
      Subroutine Invertc (c,n,ndim,ier,err)
c----------------------------------------------------------------------c

c     complex matrix inversion 
c     original code written by R. Ambrosetti, for real matrices
c     extended to complex matrices by Ivo Cacelli  march 2006

      implicit complex*16 (a-h,o-z)
      real*8 err,errmax,ar,ai
      parameter (max=550)
      dimension c(ndim,*),al(max,max),bl(max),cl(max),lz(max)
      errmax=1.d-10
c
      if(n.gt.max) goto 200
      do 10 i=1,n
      lz(i)=i
      do 10 j=1,n
      al(i,j)=c(i,j)
   10 continue
c
      do 20 i=1,n
      k=i
      yl=al(i,i)
      l=i-1
      lp=i+1
      if(n.lt.lp)  goto 14
      do 12 j=lp,n
      w=al(i,j)
      if(cdabs(w).le.cdabs(yl)) goto 12
      k=j
      yl=w
   12 continue
   14 do 16 j=1,n
      cl(j)=al(j,k)
      al(j,k)=al(j,i)
      al(j,i)=-cl(j)/yl
      al(i,j)=al(i,j)/yl
   16 bl(j)=al(i,j)
      al(i,i)=1.0/yl
      j=lz(i)
      lz(i)=lz(k)
      lz(k)=j
      do 19 k=1,n
      if (i.eq.k)  goto 19
      do 18 j=1,n
      if (i.eq.j)  goto 18
      al(k,j)=al(k,j)-bl(j)*cl(k)
   18 continue
   19 continue
   20 continue
c
      do 30 i=1,n
      if (i.eq.lz(i))  goto 30
      k=i+1
      do 28 j=k,n
      if (i.ne.lz(j))  goto 28
      m=lz(i)
      lz(i)=lz(j)
      lz(j)=m
      do  26  l=1,n
      cl(l)=al(i,l)
      al(i,l)=al(j,l)
   26 al(j,l)=cl(l)
   28 continue
   30 continue

c     checks the inversion error

      err=0.d0
      do 50 j=1,n
      do 48 k=1,n
      yl=0.d0
      do 46 l=1,n
   46 yl=yl+c(j,l)*al(l,k)
      ar=Dreal(yl)
      ai=Dimag(yl)
      if(j.eq.k) ar=ar-1.d0
      if(dabs(ar).gt.err) err=dabs(ar)
   48 if(dabs(ai).gt.err) err=dabs(ai)
   50 continue
c     write(6,*)   ' invert error = ',err,c(1,1),al(1,1),yl
      if(dabs(err).gt.errmax) goto 220
c
      do 60 i=1,n
      do 60 j=1,n
   60 c(i,j)=al(i,j)
      ier=0
      return
c
  200 write(6,202) n,max
  202 format(/,/4x,'error: insuff. dimens. in subr. Invertc'
     &        ,/4x,'dim. of the input matrix ....',i5
     &        ,/4x,'max allowed dimension .......',i5,/,/)
      stop ' ERROR in subr. Invert'
  220 write(6,222) err,errmax
  222 format(//4x,'error: loss of precision in subr. Invertc'
     &       ,/4x,'evaluated inversion error ...',d12.3
     &       ,/4x,'max allowed error ...........',d12.3,/,/)
      ier=1
      return
      end
