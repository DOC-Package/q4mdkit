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
      Subroutine Read_fchk (nat,geom,aname,etot,FC,Grad
     &                     ,iunr,maxat,max3n)
c----------------------------------------------------------------------c
 
c     reads from the formatted checkpoint file of Gaussian :
c        1. gaussian general information
c        2. geometry
c        3. energy gradients
c        4. force constant matrix
 
 
c  WARNING: common printout !! 
 
c  WARNING: if Hessian (= Force Constants) is not found,
c           then the written gradients are zero
c           i.e. they are written on file fchk but they are not
c                computed in the g09 run
 
 
      implicit real*8(a-h,o-y) , logical(z)
      parameter (maxnuc=300, maxcart=maxnuc*3, maxatlab=54 )
      dimension geom(maxat,*),FC(max3n,*),Grad(*), iord(maxnuc) 
     &         ,nuc(maxnuc),vec(maxcart*maxcart),ivec(maxnuc*3)
      character title*72, jtype*10, method*30, basis*30
     &         ,aname(maxnuc)*5, keyword(18)*35
     &         ,labat(maxatlab)*2, str5*5, cardf*90
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
 
      keyword(01)='Number of atoms                    '
      keyword(02)='Charge                             '
      keyword(03)='Multiplicity                       '
      keyword(04)='Number of electrons                '
      keyword(05)='Number of alpha electrons          '
      keyword(06)='Number of beta electrons           '
      keyword(07)='Number of basis functions          '
      keyword(08)='Number of contracted shells        '
      keyword(09)='Highest angular momentum           '
      keyword(10)='Largest degree of contraction      '
      keyword(11)='Number of primitive shells         '
      keyword(12)='SCF Energy                         '
      keyword(13)='Total Energy                       '
      keyword(14)='Atomic numbers                     '
      keyword(15)='Nuclear charges                    '
      keyword(16)='Current cartesian coordinates      '
      keyword(17)='Cartesian Gradient                 '
      keyword(18)='Cartesian Force Constants          '
 
      do i=1,maxat
        do k=1,4
          geom(i,k)=0.0d0
        enddo
      enddo
 
c     READING
 
      rewind iunr
      read (iunr,'(a)') title
      read (iunr,'(a10,a30,30x,a30)') jtype,method,basis
 
      Call Fchkrd (iunr,keyword(1),zreal,n,vec,ivec)
      nat=ivec(1)
      Call ChkDim (nat,maxnuc,'READ_FCHK: increase maxnuc')
      Call ChkDim (nat,maxat, 'READ_FCHK: increase maxat in main')
 
      Call Fchkrd (iunr,keyword(2),zreal,n,vec,ivec)
      molchr=ivec(1)
 
      Call Fchkrd (iunr,keyword(3),zreal,n,vec,ivec)
      mul=ivec(1)
 
      Call Fchkrd (iunr,keyword(4),zreal,n,vec,ivec)
      ne=ivec(1)
 
      Call Fchkrd (iunr,keyword(5),zreal,n,vec,ivec)
      na=ivec(1)
 
      Call Fchkrd (iunr,keyword(6),zreal,n,vec,ivec)
      nb=ivec(1)
 
      Call Fchkrd (iunr,keyword(7),zreal,n,vec,ivec)
      ndim=ivec(1)
 
      Call Fchkrd (iunr,keyword(8),zreal,n,vec,ivec)
      nshell=ivec(1)
 
      Call Fchkrd (iunr,keyword(10),zreal,n,vec,ivec)
      maxcn=ivec(1)
 
      Call Fchkrd (iunr,keyword(11),zreal,n,vec,ivec)
      nprim=ivec(1)
 
      Call Fchkrd (iunr,keyword(12),zreal,n,vec,ivec)
      escf=vec(1)
 
      Call Fchkrd (iunr,keyword(13),zreal,n,vec,ivec)
      etot=vec(1)
 
 
      write(6,20) title(1:lenw(title)),jtype(1:lenw(jtype))
     &           ,method(1:lenw(method)),basis(1:lenw(basis))
     &           ,nat,molchr,mul,nshell,ndim,nprim,maxcn
     &           ,na,nb,ne,escf,etot
   20 format(/10x,'============ reading Gaussian data ==========='
     &      ,/10x,'title: ',a 
     &      ,/10x,'Job.type=',a,6x,'Method=',a,6x,'Basis=',a
     &      ,/10x,'n.atoms mol.charge multipl ........',3i6
     &      ,/10x,'n.shell n.dim n.prim. max.contr ...',4i6
     &      ,/10x,'n. of electrons (alpha, beta, tot) ',3i6
     &      ,/10x,'E(scf)  E(tot) ....................',2f14.8)
 
c    reads data to rearrange the atom ordering
 
      Call Rearrange (iord,nat)
      iww=0
      do i=1,nat
        if(iord(i).ne.i) iww=iww+1
      enddo
      if(iww.gt.0) write(6,522) (iord(k),k=1,nat)
  522 format(/8x,'the atomic sequence is rearranged according to:'
     &     ,/(8x,5i4,4x,5i4))
 
c     reads atomic numbers

      Call Fchkrd (iunr,keyword(14),zreal,n,vec,ivec)
      do i=1,nat
        ik=iord(i)
        nuc(i)=ivec(ik)
      enddo
      do i=1,nat
        str5='X   '
        if(nuc(i).le.maxatlab) str5(1:2)=labat(nuc(i))
        write(str5(3:5),'(i3)') i
        Call Compact (str5,lw)
        aname(i)=str5
      enddo
 
c     reads nuclear charges
 
      Call Fchkrd (iunr,keyword(15),zreal,n,vec,ivec)
      do j=1,nat
c       geom(j,4)=vec(iord(j))  !! GP 5/5/2016
        geom(j,4)=dble(nuc(j))
      enddo
 
c     reads nuclear coordinates
 
      Call Fchkrd (iunr,keyword(16),zreal,n,vec,ivec)
      jk=0
      do j=1,nat
        jj=iord(j)*3-3
        do k=1,3
          geom(j,k)=vec(jj+k)*boh2ang
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
      Call Fchkrd (iunr,keyword(18),zreal,n,vec,ivec)
      if(n.le.0) then
        write(6,'(/10x,"Hessian not found on file fchk")')
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
      cardf='Force Constant Matrix (cartes) from FCHK har/bohr^2'
      if(levprt.ge.2) Call Pr_Rec2 (FC,n,ncart,cardf,ncol,max3n)
 
c     WARNING: if Hessian (= Force Constants) is not found,
c           then the written gradients are zero
c           i.e. they are written on file fchk but they are not
c           computed in the g09 run
 
 
c     reads cartesian gradients
 
      do j=1,nat*3
        vec(j)= -99999.d0
      enddo
      Call Fchkrd (iunr,keyword(17),zreal,n,vec,ivec)
      if(n.le.0) then
        write(6,'(/8x,"Gradient not found on file fchk")')
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
      Subroutine FChkRd (IFUnit,Label,IsReal,N,Value,IValue)
c ---------------------------------------------------------
 
      Implicit Real*8(A-H,O-Z)
 
c     Find the data specified by Label on the formatted checkpoint
c     file on unit IFUnit and return it:
c
c     IFUnit ... Unit for formatted data.
c     Label  ... Label for data.
c     IsReal ... Returned with whether the data is real.
c     N      ... Number of data items:
c                -1 ... none found.
c                 0 ... a scalar.
c                >0 ... a vector of length N.
c     Value  ... Holds returned values if real.
c     IValue ... Holds returned values if integer.
c
      Character Label*(*), Type*1, Flag*2, ValStr*22, VFlag*(*),
     &  LabelI*40, cnull*1
      Parameter (VFlag='N=')
      Logical IsReal, IsVect
      Dimension Value(*), IValue(*)
 1000 Format(A40,3X,A1,3X,A2,A22)
 1010 Format(A1)
 1020 Format((5(1PD16.8)))
 1030 Format(6I12)
 
c     First, find the data.  Start from wherever we are, then
c     if necessary rewind the file.
 
      Do 100 IPass = 1, 2
   10   Read(IFUnit,1000,End=20) LabelI, Type, Flag, ValStr
        Call FCGetN(Type,Flag,ValStr,IVal,Val)
        IsVect = Flag.eq.VFlag
        if(index(label,labeli(1:lenw(labeli))).gt.0) goto 200
        If(IsVect) then
          If(Type.eq.'I' .or. Type.eq.'C') then
            Call SkipFF(IFUnit,(IVal+5)/6)
          else if(Type.eq.'R') then
            Call SkipFF(IFUnit,(IVal+4)/5)
          endIf
        endIf
        If(Type.eq.'C') then
c         First check that we are not at EOF (i.e. if the last
c         section is of Type=C)
          Read(IFUnit,'(a)',End=20) cnull
          backspace IFUnit
          Call SkipFF(IFUnit,0)
        endif
        Goto 10
   20   Rewind IFUnit
        Call SkipFF(IFUnit,2)
  100   Continue
 
c     We fell through and didn't find the data.
 
      N = -1
      Return
 
c     Label was found.  Determine N and (if scalar) the value.
 
  200 IsReal = Type.eq.'R'
      If(IsVect) then
        N = IVal
        If(IsReal) then
          Read(IFUnit,1020) (Value(I),I=1,N)
        else
          Read(IFUnit,1030) (IValue(I),I=1,N)
          endIf
      else if(IsReal) then
        N = 0
        Value(1) = Val
      else
        N = 0
        IValue(1) = IVal
        endIf
      Return
      End


c ---------------------------------------------------------
      Subroutine FCGetN (Type,Flag,ValStr,IVal,Val)
c ---------------------------------------------------------

      Implicit Real*8(A-H,O-Z)
 
c     Given the type ('R' or 'I'), the flag ('N=' or blank) and the
c     value string from a label line in a formatted checkpoint file,
c     return the value.
 
      Character*(*) Type, Flag, ValStr
      Save Zero
      Data Zero/0.0d0/
 1000 Format(I12)
 1010 Format(D22.15)
 
      If(Type.eq.'I'.or.Flag.eq.'N=') then
        Read(ValStr,1000) IVal
        Val = Zero
      else
        Read(ValStr,1010) Val
        IVal = 0
        endIf
      Return
      End

c ---------------------------------------------------------
      Subroutine SkipFF(IFUnit,Num)
c ---------------------------------------------------------


      Implicit Integer(A-Z)
 
c     Skip N records on a fortran formatted unit.
      
 
      Character*1 CJunk
 1000 Format(A1)
 
      if(Num.eq.0) goto 200
 
      Do 100 I = 1, Num
        Read(IFUnit,1000) CJunk
  100 Continue
      Return
 
  200 Read(IFUnit,'(a)') CJunk
      if(CJunk.eq.' ') goto 200
      backspace IFUnit
      Return
      End
c----------------------------------------------------------------c
      Subroutine ChkDim (need,nreal,string)
c----------------------------------------------------------------c
 
c     checks dimension, printouts and stops calling routine
 
      character*(*) string
      if(need.gt.nreal) then
        ll=Lenw(string)
        write(6,8) string,need,nreal
        stop
      endif
      return
    8 format(//8x,'DIMENSION ERROR in ',a
     &       ,/8x,'requested dimension .............',i6
     &       ,/8x,'effective dimension .............',i6,//)
      end




c----------------------------------------------------------------c
      Subroutine Rearrange (iord,nat)
c----------------------------------------------------------------c
 
c      rearrange the atom ordering
c      es: rearrange 6 18 19 20 21 22 23 
c      atom 6 in first position, atom 18 in second pos. .....
 
 
      dimension iord(*)
      character card*500
      data iread/0/
      save iread
      iread=iread+1   ! re-inseerted by GP & SG 23/2/2024
      if(iread.gt.1) return

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
 
c----------------------------------------------------------------c
      Subroutine Get_at_xyz (numb,iato,ixyz)
c----------------------------------------------------------------c
 
      if(numb.le.0) goto 40
      iato=(numb+2)/3
      ixyz=numb-3*(iato-1)
      return
   40 write(6,'(//" ERROR: numb <= 0 in GET_AT_XYZ",//)')
      stop
      end
