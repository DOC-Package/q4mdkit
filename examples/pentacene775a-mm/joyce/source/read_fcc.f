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
      Subroutine Read_FCC (nat,geom,aname,etot,FC,Grad
     &                     ,iunr,maxat,max3n)
c----------------------------------------------------------------------c
 
c     reads from the fcc state file from FCclasses3:
c        1. geometry
c        2. energy gradients
c        3. force constant matrix
 
 
c  WARNING: common printout !! 
 
 
      implicit real*8(a-h,o-y) , logical(z)
      parameter (maxnuc=300, maxcart=maxnuc*3, maxatlab=54 )
      dimension geom(maxat,*),FC(max3n,*),Grad(*), iord(maxnuc) 
     &         ,nuc(maxnuc),vec(maxcart*maxcart),ivec(maxnuc*3)
      character title*72, jtype*10, method*30, basis*30
     &         ,aname(maxnuc)*5, keyword(4)*35
     &         ,labat(maxatlab)*2, str5*5, cardf*90
     &         ,cvec(maxnuc)*5
      common /printout/ levprt
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      data labat /'H ','He'
     &           ,'Li','Be',  'B ','C ','N ','O ','F ','Ne'
     &           ,'Na','Mg',  'Al','Si','P ','S ','Cl','Ar'
     &           ,'K ','Ca',  'Sc','Ti','V ','Cr','Mn','Fe','Co'
     &      ,'Ni','Cu','Zn',  'Ga','Ge','As','Se','Br','Kr'
     &           ,'Rb','Sr',  'Y ','Zr','Nb','Mo','Tc','Ru','Rh'
     &      ,'Pd','Ag','Cd',  'In','Sn','Sb','Te','I ','Xe'/
      save iord
 
      keyword(1)='ENER '
      keyword(2)='GRAD '
      keyword(3)='HESS '
 
      do i=1,maxat
        do k=1,4
          geom(i,k)=0.0d0
        enddo
      enddo
 
c     READING
 
      Call FccRd_Info(iunr,title)
     
      zreal = .true.
      n = 1
      Call FccRd (iunr,keyword(1),zreal,n,vec,ivec)
      etot=vec(1)

      Call FccRd_Geom (iunr,nat,cvec,vec)

      write(6,20) title(1:lenw(title)),nat,etot
   20 format(/3x,'============ reading FCC data ==========='
     &      ,/3x,'title: ',a 
     &      ,/3x,'n.atoms ...........................',i6
     &      ,/3x,'E(tot) ............................',f15.8)

c    reads data to rearrange the atom ordering
 
      Call RearrangeFCC (iord,nat)
      iww=0
      do i=1,nat
        if(iord(i).ne.i) iww=iww+1
      enddo
      if(iww.gt.0) write(6,522) (iord(k),k=1,nat)
  522 format(/8x,'the atomic sequence is rearranged according to:'
     &     ,/(8x,5i4,4x,5i4))


c     arrange atomic labels

      do i=1,nat
        ik=iord(i)
        str5='X   '
        str5(1:2)=cvec(ik)(1:2)
        write(str5(3:5),'(i3)') i
        Call Compact (str5,lw)
        aname(i)=str5
      enddo

c     get nuclear charges from names
 
      do i=1,nat
        ik=iord(i)
        call AtNum_from_AtName((cvec(ik)),nuc(i))
        geom(i,4)=dble(nuc(i))
      enddo

 
c     arrange nuclear coordinates
 
      jk=0
      do j=1,nat
        jj=iord(j)*3-3
        do k=1,3
          geom(j,k)=vec(jj+k)
        enddo
      enddo
      if(levprt.gt.0) write(6,48) (i,aname(i),(geom(i,k),k=1,3)
     &                           ,nuc(i),geom(i,4),   i=1,nat)
   48 format(/12x,'---------- GEOMETRY DATA (angstrom) ----------'
     &      ,/18x,'x',11x,'y',11x,'z',5x,'atom.n  Nucl.ch'
     &      ,/(1x,i3,2x,a5,3f12.7,1x,i5,f8.1))
 
c    reads cartesian force constants (hessian)
 
      do j=1,(nat*3)**2
        vec(j)=0.0d0
      enddo
      n = 3*nat*(3*nat+1)/2
      zreal = .true.
      rewind(iunr)
      Call FccRd (iunr,keyword(3),zreal,n,vec,ivec)
      if(n.le.0) then
        write(6,'(/10x,"Hessian not found on file fcc")')
      endif
      n=nat*3
      ij=0
      do j=1,n
        Call Get_at_xyz (j,ja,jx)
        jj=iord(ja)
        jk=(jj-1)*3+jx
        do i=1,j
          ij=ij+1
          Call Get_at_xyz (i,ia,ix)
          ii=iord(ia)
          ik=(ii-1)*3+ix
          ijk=ik*(ik-1)/2+jk
          if(jk.gt.ik) ijk=jk*(jk-1)/2+ik
          FC(i,j)=vec(ijk)
          FC(j,i)=vec(ijk)
        enddo
      enddo
      ncol=10
      ncart=n
      cardf='Force Constant Matrix (cartes) from FCC har/bohr^2'
      if(levprt.ge.2) Call Pr_Rec2 (FC,n,ncart,cardf,ncol,max3n)
 
c     reads cartesian gradients
 
      do j=1,nat*3
        vec(j)= -99999.d0
      enddo
      n = 3*nat
      zreal = .true.
      Call FccRd (iunr,keyword(2),zreal,n,vec,ivec)
      if(n.le.0) then
        write(6,'(/8x,"Gradient not found on file fcc")')
      endif
      do i=1,nat*3
        Call Get_at_xyz (i,ia,ix)
        ii=iord(ia)
        ik=(ii-1)*3+ix
        Grad(i)=vec(ik)
      enddo
      if(levprt.gt.0)then
        write(6,52) n
        do i=1,nat
          l=3*i-2
          m=3*i
          write(6,54) i,aname(i),(1.d3*Grad(k),k=l,m),(FC(k,k),k=l,m)
        enddo
      endif
   52 format(/9x,'Cartesian Energy Gradients  (mH/bohr)',i5
     &          ,'    diagonal Hessian(H/bohr)')
   54 format(i5,2x,a,3f11.5,3x,3f11.5)
 
      return
      end



c ---------------------------------------------------------
      Subroutine FccRd (IFUnit,Label,IsReal,N,Value,IValue)
c ---------------------------------------------------------
 
      Implicit Real*8(A-H,O-Z)
 
c     Find the data specified by Label on the FCC
c     file on unit IFUnit and return it:
c
c     IFUnit ... Unit for formatted data.
c     Label  ... Label for data.
c     IsReal ... Returned with whether the data is real (input).
c     N      ... Number of data items (input):
c                -1 ... none found.
c                >0 ... a vector of length N.
c     Value  ... Holds returned values if real.
c     IValue ... Holds returned values if integer.
c
      Character Label*(*), Type*1, Flag*2, ValStr*22, VFlag*(*)
      Character LabelI*40, cnull*1
      Parameter (VFlag='N=')
      Logical IsReal, IsVect
      Dimension Value(*), IValue(*)
 1000 Format(A40,3X,A1,3X,A2,A22)
 1010 Format(A1)
 1020 Format((5(1PD16.8)))
 1030 Format(6I12)
 
c     First, find the data.  Start from wherever we are, then
c     if necessary rewind the file.

      If (label .eq. 'GEOM ') then
          write(6,*) 'Error in FccRd:'
          write(6,*) ' This routine cannot be used to read GEOM'
          stop
      endif
 
      Do 100 IPass = 1, 2
   10   Read(IFUnit,1000,End=20) LabelI
        if(index(labeli(1:lenw(labeli)),
     &           label(1:lenw(label))).gt.0) goto 200
        Goto 10
   20   Rewind IFUnit
  100   Continue
 
c     We fell through and didn't find the data.
 
      N = -1
      Return
 
c     Label was found.  Determine N and (if scalar) the value.
 
  200 If(IsReal) then
        Read(IFUnit,*) (Value(I),I=1,N)
      else
        Read(IFUnit,*) (IValue(I),I=1,N)
      endIf
      Return
      End

c ---------------------------------------------------------
      Subroutine FccRd_Geom (IFUnit,Nat,AtName,Geo)
c ---------------------------------------------------------
 
      Implicit Real*8(A-H,O-Z)
 
c     Find the data specified by Label on the FCC
c     file on unit IFUnit and return it:
c
c     IFUnit ... Unit for formatted data.
c     Nat    ... Number of atoms (output)
c     AtName ... Holds returned atomnames (char)
c     Geo    ... Holds returned coordinates
c
      Character LabelI*40, cnull*1
      character AtName(*)*(*)
      Dimension Geo(*)
 1000 Format(A40,3X,A1,3X,A2,A22)
 1010 Format(A1)
 1020 Format((5(1PD16.8)))
 1030 Format(6I12)
 
c     First, find the data.  Start from wherever we are, then
c     if necessary rewind the file.

      Do 100 IPass = 1, 2
   10   Read(IFUnit,'(A)',End=20) LabelI
        if(index(labeli(1:lenw(labeli)),'GEOM').gt.0) goto 200
        Goto 10
   20   Rewind IFUnit
  100   Continue
 
c     We fell through and didn't find the data.
 
      Nat = -1
      Return
 
c     Label was found.  Determine N and (if scalar) the value.
 
  200 Read(IFUnit,*) Nat
      Read(IFUnit,*) cnull
      do i = 1, Nat
        j = 3*(i-1)
        Read(IFUnit,*) AtName(i), Geo(j+1), Geo(j+2), Geo(j+3)
      enddo
      Return
      End

c ---------------------------------------------------------
      Subroutine FccRd_Info (IFUnit,Title)
c ---------------------------------------------------------
 
      Implicit Real*8(A-H,O-Z)
 
c     Find the data specified by Label on the FCC
c     file on unit IFUnit and return it:
c
c     IFUnit ... Unit for formatted data.
c     Title  ... Title (output)
c
      Character LabelI*40
      character Title*(*)
 1000 Format(A40,3X,A1,3X,A2,A22)
 1010 Format(A1)
 1020 Format((5(1PD16.8)))
 1030 Format(6I12)
 
c     First, find the data.  Start from wherever we are, then
c     if necessary rewind the file.

      Do 100 IPass = 1, 2
   10   Read(IFUnit,'(A)',End=20) LabelI
        if(index(labeli(1:lenw(labeli)),'INFO').gt.0) goto 200
        Goto 10
   20   Rewind IFUnit
  100   Continue
 
c     We fell through and didn't find the data.
 
      N = -1
      Return
 
c     Label was found.  Determine N and (if scalar) the value.
 
  200 Read(IFUnit,'(A)') Title
      Return
      End

c ---------------------------------------------------------
      Subroutine AtNum_from_AtName(atname,atnum)
c ---------------------------------------------------------
 
      Implicit Real*8(A-H,O-Z)
 
c     Get the atom number for a given atom name
c
c     atname ... Atom name   (input)
c     atnum  ... Atom number (output)

      Character atname*(*)
      integer   atnum

      character  atname_from_atnum(103)*10
      parameter (atname_from_atnum = 
     &(/'H ','He',
     &  'Li','Be','B ','C ','N ','O ','F ','Ne',
     &  'Na','Mg','Al','Si','P ','S ','Cl','Ar',
     &  'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga'
     & ,'Ge','As','Se','Br','Kr',
     &  'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In'
     & ,'Sn','Sb','Te','I ','Xe',
     &  'Cs','Ba','La', 
c       Lantanides:  
c       ---------------------------------------------------
     &  'Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     &  'Tb','Dy','Ho','Er','Tm','Yb','Lu',
c       ---------------------------------------------------
     &  'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po'
     & ,'At','Rn',
     &  'Fr','Ra','Ac',
c       Actinides:
c       ---------------------------------------------------
     &  'Th','Pa','U ','Np','Pu','Am','Cm',
     &  'Bk','Cf','Es','Fm','Md','No','Lr'
c       ---------------------------------------------------
     & /))

      n = size(atname_from_atnum)
      atnum = 0
      do i=1,n
          if (adjustl(atname)==atname_from_atnum(i)) then
              atnum = i
              exit
          endif
      enddo

      Return
      End

c----------------------------------------------------------------c
      Subroutine RearrangeFCC (iord,nat)
c----------------------------------------------------------------c
 
c      rearrange the atom ordering
c      es: rearrange 6 18 19 20 21 22 23 
c      atom 6 in first position, atom 18 in second pos. .....
 
 
      dimension iord(*)
      character card*500
      data ireadFCC/0/
      save ireadFCC
      ireadFCC=ireadFCC+1
      if(ireadFCC.gt.1) return

      do i=1,nat
        iord(i)=i
      enddo
      iline=0
      rewind 5
      Call Cerca_L ('$rearr',5,card,ier,.true.,iline)
      if(ier.ne.0) return
 
      Call Count_W (card,nword)
      read(card,*) (iord(k),k=1,nword)
      do i=1,nword
        if(iord(i).gt.nat .or. iord(i).lt.1) then
          write(6,524) iord(i),(iord(k),k=1,nat)
          stop
        endif
      enddo
  524 format(/9x,'REARRANGE ERROR: wrong atom',i4,/(15x,10i4))
 
      do j=2,nat
        do i=1,j-1
          if(iord(j).eq.iord(i)) then
            write(6,526) iord(j),(iord(k),k=1,nat)
            stop
          endif
        enddo
      enddo
  526 format(/9x,'REARR ERROR: atom',i4,' repeated',/(15x,10i4))
 
      return
      end
