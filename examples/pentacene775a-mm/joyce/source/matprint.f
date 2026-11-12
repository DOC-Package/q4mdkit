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
c                          matprint.f
c            package of subroutines to print matrices
c----------------------------------------------------------------------c

c----------------------------------------------------------------------c
      Subroutine Pr_Tri (a,nsym,ndsym,title,ncol)
c----------------------------------------------------------------------c

c     prints a matrix stored as triangular matrix in dictionary order
c     only the upper triangular part of the matrix is showed.

c     a        matrix to be printed
c     nsym     number of diagonal blocks
c     ndsym    dimension of diagonal blocks
c     title    (char) title to be printed
c     ncol     number of columns to be printed at the
c              same time for each row

      implicit real*8(a-h,o-z)
      parameter (nfix=15)
      character*8 strin(nfix),blank,stri(2)
      character*(*) title
      dimension a(*),ndsym(*),nkk(nfix),www(nfix)
      data blank /'   . .  '/
c
      nl=len(title)
      k1=1
      dowhile ( title(k1:k1).eq.' ' .and. k1.lt.nl )
        k1=k1+1
      enddo
      k2=nl
      dowhile ( title(k2:k2).eq.' ' .and. k2.gt.0 )
        k2=k2-1
      enddo
      if(k1.le.k2) write(6,2) title(k1:k2)
    2 format(/,/3x,'-- ',a,' --')
      if(nsym.le.0) then
      write(6,4) nsym
    4 format(/,/6x,'wrong numb. of diag. blocks in Subr. Pr_tri'
     &       ,i6,/)
      return
      endif
c
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      nskipt=0
      nskipl=0
      do 30 isym=1,nsym
c
      nt=ndsym(isym)
      if(nt.le.0) then
      write(6,6) (ndsym(k),k=1,nsym)
    6 format(/,/6x,'Wrong diag. block dimens. in Subr. Pr_tri'
     &        ,/(4X,10I5))
      return
      endif
      num=nt/nco
      nres=nt-num*nco
      if(nres.gt.0) num=num+1
      if(nres.eq.0) nres=nco
      isu=0
      do 20 ks=1,num
      inf=isu+1
      isu=isu+nco
      if(isu.gt.nt) isu=nt
      ndim=isu-inf+1
      do 14 i=1,ndim
   14 nkk(i)=inf+i-1+nskipl
      if(nsym.gt.1) write(6,38) isym
      write(6,40)  (nkk(i),i=1,ndim)
c
      ijs=(inf-2)*(inf-1)/2+nskipt
      do 18 i=1,isu
      ij=ijs+i
      ja=0
      do 16 j=inf,isu
      ja=ja+1
      stri(1)=blank
      ij=ij+j-1
      if(i.gt.j) goto 16
      www(1)=a(ij)
      Call Fmtsel (stri,www,.false.,ix,1,1,1)
   16 strin(ja)=stri(1)
      ia=i+nskipl
      write(6,42) ia,(strin(j),j=1,ndim)
   18 continue
   20 continue
c
      nskipl=nskipl+nt
   30 nskipt=nskipt+nt*(nt+1)/2
      return
   38 format(4x,'diagonal block number',i3)
   40 format(4x,15i8)
   42 format(1x,i4,1x,15a8)
      end

c----------------------------------------------------------------------c
      Subroutine Pr_Rec (x,nv,nh,title,ncol)
c----------------------------------------------------------------------c

c     prints a rectangular matrix stored by columns
c            in an one-dimensional array
c
c     X       matrix to be printed
c     nh      number of columns   (horizontal dimensions)
c     nv      number of rows      (vertical dimensions)
c     ncol    number of columns to be printed in a same row
c     title   (char) title to be printed

      implicit real*8(a-h,o-z)
      parameter (nfix=15)
      character strin(nfix)*8,stri*8,title*(*)
      dimension X(*)
c
      nl=len(title)
      k1=1
      dowhile ( title(k1:k1).eq.' ' .and. k1.lt.nl )
        k1=k1+1
      enddo
      k2=nl
      dowhile ( title(k2:k2).eq.' ' .and. k2.gt.0 )
        k2=k2-1
      enddo
      if(k1.le.k2) write(6,2) title(k1:k2)
    2 format(/,/3x,'-- ',a,' --')
c
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      lsup=0
      dowhile (lsup.lt.nh)
        linf=lsup+1
        lsup=lsup+nco
        if(lsup.gt.nh) lsup=nh
        nw=lsup-linf+1
        write(6,40)  (i,i=linf,lsup)
        kinf=(linf-1)*nv+1
        ksup=(lsup-1)*nv+1
        do l=1,nv
          Call Fmtsel (strin,X,.false.,ja,kinf,ksup,nv)
          write(6,42) l,(strin(j),j=1,ja)
          kinf=kinf+1
          ksup=ksup+1
        enddo
      enddo
c
      return
   40 format(4x,15i8)
   42 format(1x,i4,1x,15a8)
      end

c----------------------------------------------------------------------c
      Subroutine Pr_Rec2 (x,nv,nh,title,ncol,ndim)
c----------------------------------------------------------------------c

c     prints a two-dimensional array 
c
c     X       matrix to be printed
c     nh      n. of columns to be printed (horizontal)
c     nv      n. of rows to be printed (vertical)
c     ncol    n. of columns to be printed in a same row
c     title   (char) title to be printed
c     ndim    first dimension of the X array in the caller

      implicit real*8(a-h,o-z)
      parameter (nfix=15)
      character strin(nfix)*8,stri*8,title*(*)
      dimension x(ndim,*),xx(nfix)
      logical zok

c     set an automatic threshold

      nn=0
      aa=0.0d0
      do i=1,nv
        do j=1,nh
          if(dabs(x(i,j)).gt.1.d-12) then
            nn=nn+1
            aa=aa+dabs(x(i,j))
          endif
        enddo
      enddo
      am=aa/dble(nn)
      thr=am*1.d-10
c
      nl=len(title)
      k1=1
      dowhile ( title(k1:k1).eq.' ' .and. k1.lt.nl )
        k1=k1+1
      enddo
      k2=nl
      dowhile ( title(k2:k2).eq.' ' .and. k2.gt.0 )
        k2=k2-1
      enddo
      if(k1.le.k2) write(6,2) title(k1:k2)
    2 format(/,/3x,'-- ',a,' --')

c     test for null matrix

      sum=0.d0
      do i=1,nv
        do j=1,nh
          sum=sum+dabs(X(i,j))
        enddo
      enddo
      if(sum/dble(nv*nh).lt.1.d-09) then
        write(6,'(8x," ----- this is a null matrix ------")')
        return
      endif
c
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      lsup=0
      dowhile (lsup.lt.nh)
        linf=lsup+1
        lsup=lsup+nco
        if(lsup.gt.nh) lsup=nh
        nw=lsup-linf+1
        write(6,40)  (i,i=linf,lsup)
        do k=1,nv
          zok=.false.
          do l=linf,lsup
            xx(l-linf+1)=X(k,l)
            if(dabs(X(k,l)).gt.thr) zok=.true.
          enddo
          if(zok) Call Fmtsel (strin,xx,.false.,ja,1,nw,1)
          if(zok) write(6,42) k,(strin(j),j=1,ja)
        enddo
      enddo
c
      return
   40 format(4x,15i8)
   42 format(1x,i4,1x,15a8)
      end

c----------------------------------------------------------------------c
      Subroutine Fmtsel (strin,xx,zexp,ind,min,max,incr)
c----------------------------------------------------------------------c

c     write the value of (xx(i),i=min,max,incr) on the
c     array character (string(j),j=1,ind)
c       zexp = t  if only zero's are contained in *string* then
c                 try to write in exponential form.

      double precision xx(*),x
      character*(*) strin(*)
      character strfmt*13, stk*1
      logical zexp
c
      nlt=Len(strin(1))
      if(nlt.lt.4) goto 85
      nl=nlt
      if(nl.gt.20) nl=20
      n=nl-1
c
      ind=0
      do 80 ix=min,max,incr
      x=xx(ix)
      ind=ind+1
c
      do 12 k=1,nlt
   12 strin(ind)(k:k)=' '

c     tries to print with floating point format

      ndec=n-3
   20 write(strfmt,90) n,ndec
      write(strin(ind),fmt=strfmt) x
      if(strin(ind)(2:2).ne.'*') goto 50
      ndec=ndec-1
      if(ndec.gt.0) goto 20

c     exponential format is employed

   30 ndec=n-7
      if(ndec.lt.0) goto 40
      write(strfmt,92) n,ndec
      write(strin(ind),fmt=strfmt) x
      if(dabs(x).lt.1.d-12) write(strin(ind),'("    0   ")')
      goto 80

c     null string is returned

   40 do 45 k=1,nl
   45 strin(ind)(k:k)=' '
      goto 80

c     check on the significant digits .ne. 0

   50 if(.not.zexp) goto 80
      do 55 k=1,nl
      stk=strin(ind)(k:k)
      if(stk.eq.'.'.or.stk.eq.'-') goto 55
      if(stk.eq.'0'.or.stk.eq.' ') goto 55
      goto 80
   55 continue
      ndec=n-7
      write(strfmt,92) n,ndec
      write(strin(ind),fmt=strfmt) x
      if(dabs(x).lt.1.d-13) write(strin(ind),'("    0   ")')
   80 continue
   85 return
   90 format('(1x,f',i2,'.',i2,')  ')
   92 format('(1x,1pd',i2,'.',i2,')')
      end

c----------------------------------------------------------------------c
      Subroutine Pr_Orb (v,e,n,nf,label,thrs,ncol)
c----------------------------------------------------------------------c

c     print vectors and energies
c
c     v              vectors (by column)
c     e              energies 
c     n              number of vectors to be printed
c     nf             vector dimension
c     label  (char)  label of each vector component
 
      implicit real*8 (a-h,o-z) 
      parameter (nfix=14)
      dimension  v(*),e(*),ntx(nfix)
      character*(*) label(*)
      character*40  string
c 
      nchar=len(label(1))
      do i=1,nchar
        string(i:i)=' '
      enddo
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      nv=nf
      nh=n 
      nwri=nh/nco
      if(nco*nwri.lt.nh) nwri=nwri+1 
c
      lsu=0 
      kin=0 
      do 40 iw=1,nwri 
      lin=lsu+1 
      lsu=lsu+nco 
      if(lsu.gt.nh) lsu=nh
      write(6,100) string(1:nchar),(e(k),k=lin,lsu) 
      write(6,101) string(1:nchar),(k,k=lin,lsu) 
      ninc=(lsu-lin)*nv 
      do 38 l=1,nf 
      kin=kin+1 
      ksu=kin+ninc
      iopr=0
      do 36 k=kin,ksu,nv
   36 if(dabs(v(k)).gt.thrs) iopr=1
      if(iopr.eq.0) goto 38
      write(6,102)  l,label(l),(v(k),k=kin,ksu,nv) 
   38 continue
      kin=ksu 
   40 continue
      return
  100 format(/1x,5x,a,14f8.1) 
  101 format( 1x,5x,a,i6,13i8) 
  102 format( 1x,i4,1x,a,1x,14f8.3) 
      end 

c----------------------------------------------------------------------c
      Subroutine Pr_Orb2 (v,e,n,nf,label,thrs,ncol,ndim)
c----------------------------------------------------------------------c

c     print vectors and energies
c
c     v              vectors (by column)
c     e              energies 
c     n              number of vectors to be printed
c     nf             vector dimension
c     label  (char)  label of each vector component
 
      implicit real*8 (a-h,o-z) 
      parameter (nfix=14)
      dimension  v(ndim,*),e(*),ntx(nfix)
      character*(*) label(*)
      character*40  string
c 
      nchar=len(label(1))
      do i=1,nchar
        string(i:i)=' '
      enddo
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      nv=nf
      nh=n 
      nwri=nh/nco
      if(nco*nwri.lt.nh) nwri=nwri+1 
c
      ksu=0 
      do 40 iw=1,nwri 
      kin=ksu+1 
      ksu=ksu+nco 
      if(ksu.gt.nh) ksu=nh
      write(6,100) string(1:nchar),(e(k),k=kin,ksu) 
      write(6,101) string(1:nchar),(k,k=kin,ksu) 
      do 38 i=1,nf 
      iopr=0
      do 36 k=kin,ksu
   36 if(dabs(v(i,k)).gt.thrs) iopr=1
      if(iopr.eq.0) goto 38
      write(6,102)  i,label(i),(v(i,k),k=kin,ksu) 
   38 continue
   40 continue
      return
  100 format(/1x,5x,a,14f8.1) 
  101 format( 1x,5x,a,i6,13i8) 
  102 format( 1x,i4,1x,a,1x,14f8.3) 
      end 

c----------------------------------------------------------------------c
      Subroutine Pr1_Orb (v,e,n,nf,thrs,ncol)
c----------------------------------------------------------------------c

c     print vectors and energies
c
c     v              vectors (by column)
c     e              energies 
c     n              number of vectors to be printed
c     nf             vector dimension
 
      implicit real*8 (a-h,o-z) 
      parameter (nfix=14)
      dimension v(*),e(*),ntx(nfix)
c 
      nco=ncol
      if(nco.gt.nfix.or.nco.le.0) nco=nfix
      nv=nf
      nh=n 
      nwri=nh/nco
      if(nco*nwri.lt.nh) nwri=nwri+1 
c
      lsu=0 
      kin=0 
      do 40 iw=1,nwri 
      lin=lsu+1 
      lsu=lsu+nco 
      if(lsu.gt.nh) lsu=nh
      write(6,100) (e(k),k=lin,lsu) 
      write(6,101) (k,k=lin,lsu) 
      ninc=(lsu-lin)*nv 
      do 38 l=1,nf 
      kin=kin+1 
      ksu=kin+ninc
      iopr=0
      do 36 k=kin,ksu,nv
   36 if(dabs(v(k)).gt.thrs) iopr=1
      if(iopr.eq.0) goto 38
      write(6,102)  l,(v(k),k=kin,ksu,nv) 
   38 continue
      kin=ksu 
   40 continue
      return
  100 format(/1x,4x,14f8.4) 
  101 format( 1x,4x,i6,13i8) 
  102 format( 1x,i3,1x,14f8.3) 
      end 

c----------------------------------------------------------------------c
      Subroutine Matric (a,nt,title,nfix)
c----------------------------------------------------------------------c

c     prints a lower triangular matrix (a) (dimension nt) in nfix columns
c        1
c        2  3
c        4  5  6
c       .. .. .. ..

      real*8 a(*)
      integer nkb(10)
      character*(*) title
c
      nl=len(title)
      k1=1
      dowhile ( title(k1:k1).eq.' ' .and. k1.lt.nl )
        k1=k1+1
      enddo
      k2=nl
      dowhile ( title(k2:k2).eq.' ' .and. k2.gt.0 )
        k2=k2-1
      enddo
      if(k1.le.k2) write(6,28) title(k1:k2)
c
      num=nt/nfix
      nres=nt-num*nfix
      if(nres)  2,2,1
    1 num=num+1
    2 m=0
      do 9 ks=1,num
      n=m+1
      m=m+nfix
      if(ks-num) 4,3,4
    3 m=nt
    4 do 5 i=1,nfix
      ik=n+i-1
    5 nkb(i)=ik
      nwr=m-n+1
      write(6,30) (nkb(i),i=1,nwr)
      jj=(n*n-n)/2
      do 8 j=n,nt
      min=jj+n
      mm=j
      if(j-m) 7,7,6
    6 mm=m
    7 max=jj+mm
      write(6,32) j,(a(k),k=min,max)
    8 jj=jj+j
    9 continue
c
      return
   28 format(/,/3x,'-- ',a,' --')
   30 format(/3x,10i8)
   32 format(1x,i3,2x,10f8.4)
      end
