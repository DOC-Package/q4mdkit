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
      Subroutine Get_Bonds (nat,geom,aname,dist,kbond,ndim,iout)
c----------------------------------------------------------------------c

c     computes the internuclear distances and the bond order

      implicit real*8(a-h,o-y), logical(z)
      parameter (max=300)
      dimension geom(ndim,*),dist(ndim,*),kbond(ndim,*)
     &         ,ibon(max),rbon(max),izat(max),radi(max)
     &         ,obond(max,max)
      character aname(*)*(*),abon(max)*1,At_Symb*5
      logical zdetail
      common /printout/ levprt
      zdetail=.false.

      data dmax /5.0d0/, tol/0.0001d0/
      if(nat.lt.1)   stop 'Error: Get_Bonds: nat<1'
      if(nat.gt.max) stop 'Error: Get_Bonds: max<nat'
      zdeb=.true.
      iou=6
      if(iout.gt.0 .and. iout.le.9) iou=iout
      if(nat.eq.1) then
        write(iou,'(   4x,"--- Cartes. coordinates (Angstrom) ---")') 
        write(iou,'(4x,"Atom      X",8x,"Y",8x,"Z")') 
        write(iou,42) 1,aname(1),(geom(1,k),k=1,3)
      endif

c     special case: atom

      if(nat.eq.1) then
        dist(1,1)=0.0d0
        kbond(1,1)=0
        return
      endif

c     set the interatomic distances (stored in /dist/)

      dsmin=1.d6
      do i=1,nat
        iz=geom(i,4)+tol
        if(iz.lt.0) iz=0
        izat(i)=iz
        do j=1,i
          rij=0.d0
          do k=1,3
            rij=rij+(geom(i,k)-geom(j,k))**2
          enddo
          dist(i,j)=dsqrt(rij)
          dist(j,i)=dist(i,j)
          if(i.ne.j .and. dist(j,i).lt.dsmin) then
            iii=i
            jjj=j
            dsmin=dist(j,i)
          endif
        enddo
      enddo

c     set the bond order for the atom-atom pairs

      do i=1,nat
        iz=izat(i)
        do j=1,nat
          kbond(i,j)=0
          obond(i,j)=0.d0
          if(i.ne.j .and. dist(i,j).lt.dmax) then
            jz=izat(j)
            rij=dist(i,j)
            Call Bond_Ord (iz,jz,rij,bondor,ipbo)
            kbond(i,j)=ipbo
            obond(i,j)=bondor
          endif
        enddo
      enddo

c     print geometry and distance between bonded atoms
      if(levprt .gt. 0) then 
        write(iou,'(/,/10x,"---     BOND ORDER CALCULATION     ---")') 
        write(iou,'(   10x,"--- Cartes. coordinates (Angstrom) ---")') 
        write(iou,'(10x,"Atom      X",8x,"Y",8x,"Z",9x,"bonds")') 
      endif
      do i=1,nat
        nb=0
        do k=1,nat
          if(kbond(i,k).gt.0) then
            nb=nb+1
            ibon(nb)=k
            rbon(nb)=dist(i,k)
            abon(nb)='/'
            if(kbond(i,k).gt.1) abon(nb)='#'
          endif
        enddo
        if(nb.gt.0 .and. levprt .gt. 0) then
          write(iou,42) i,aname(i),(geom(i,k),k=1,3)
     &              ,(ibon(k),abon(k),rbon(k),k=1,nb)
        elseif (levprt .gt. 0) then 
          write(iou,42) i,aname(i),(geom(i,k),k=1,3)
        endif
      enddo
      if (levprt .gt. 0)
     & write(iou,'(14x," ------- min distance  ",2i4,f9.3)')
     &       iii,jjj,dsmin
   42 format(14x,i3,2x,a,3f9.5,1x,4(i3,a,f5.3):
     &                  ,/(38x,4(i3,a,f5.3)))

c     print covalent atomic radii

      jz=0
      nt=0
      do i=1,nat
        iz=izat(i)
        iii=0
        if(nt.gt.0) then
          do k=1,nt
            if(ibon(k).eq.iz) iii=iii+1
          enddo
        endif
        if(iii.eq.0) then
          Call Bond_Ord (iz,jz,ri,BondOr,ipbo)
          nt=nt+1
          ibon(nt)=iz
          radi(nt)=ri
        endif
      enddo
      if (levprt .gt. 0) write(iou,54)
      do i=1,nt
        if (levprt .gt. 0)
     &  write(iou,56) ibon(i),At_Symb(ibon(i)),radi(i)
      enddo
   54 format(/'  atomic covalent radii used for B.O. calc')
   56 format(i6,2x,a5,f10.3)

c     print the bond orders in some details

      if(.not.zdetail) goto 70
      write(iou,60)
      do i=1,nat
        it=0
        bt=0.0d0
        zn=.true.
        do j=1,nat
          bij=obond(i,j)
          if(bij.gt.0.5d0) then
            kij=kbond(i,j)
            if(kbond(i,j).gt.0) bt=bt+obond(i,j)
            it=it+kbond(i,j)
            if(zn) write(iou,62) i,aname(i),j,aname(j),bij,kij,bt,it
            if(.not.zn) write(iou,63)       j,aname(j),bij,kij,bt,it
            zn=.false.
          endif
        enddo
      enddo
   60 format(/'   Details of the Bond Order calculation')
   62 format(i4,1x,a4,t12,i4,1x,a4,f9.4,'  -->',i3,f10.3,i3)
   63 format(t12,i4,1x,a4,f9.4,'  -->',i3,f10.3,i3)

c     search for aromatic bonds

   70 tolord=0.15d0
      tolbon=10.d0
      do 75 i=1,nat
        if(izat(i).lt.6 .or. izat(i).gt.8) goto 75

c       look for two bonds with BO>1.3 and angle = [50-70]

        m1=0
        m2=0
        do j=1,nat
          if(dabs(obond(i,j)-1.5d0).lt.tolord) then
            if(m1.gt.0 .and. m2.eq.0) m2=j
            if(m1.eq.0) m1=j
          endif
        enddo
        if(m1.eq.0 .or. m2.eq.0) goto 75
        m3=0
        Call Value_Ang (geom,ndim,m1,i,m2,m3,val,'d')
        if(dabs(val-120.d0).gt.tolbon) goto 75
        kbond(i,m1)=6
        kbond(i,m2)=6
   75 continue

c     print table of bonds

      if (levprt .gt. 0) write(iou,80)
   80   format(/4x,'-----------  Bonds between atoms -------------'
     &        ,/7x,' n.bonds |0=nobond  1,2,3=sing,doub,tripl  6=arom')
      do i=1,nat
        at=0.d0
        bt=0.d0
        do j=1,nat
          it=kbond(i,j)
          if(it.lt.4) at=at+kbond(i,j)
          if(it.eq.6) at=at+1.5d0
          if(it.gt.0) bt=bt+obond(i,j)
        enddo
        it=at+0.49d0
        ibon(i)=it
        if (levprt .gt. 0)
     &  write(iou,82) i,aname(i),it,bt,(kbond(i,k),k=1,nat)
      enddo
   82 format(i3,'-',a3,i3,f5.2,t18,4(10i1,1x),/(17x,4(10i1,1x)))

c     Check the number of bonds for some atom types

      if (levprt .gt. 0) then 
      write(6,*) '  '
      do i=1,nat
        nbi=ibon(i)
        if(nbi.eq.0)                    write(iou,100) i,aname(i)
        if(izat(i).eq.1 .and. nbi.ne.1) write(iou,101) i,aname(i)
        if(izat(i).eq.6 .and. nbi.ne.4) write(iou,106) i,aname(i)
        if(izat(i).eq.7 .and. nbi.ne.3) write(iou,107) i,aname(i)
        if(izat(i).eq.8 .and. nbi.ne.2) write(iou,108) i,aname(i)
        if(izat(i).eq.9 .and. nbi.ne.1) write(iou,109) i,aname(i)
      enddo
      endif
  100 format('  WARNING: atom ',i2,1x,a,' has no bond',/)
  101 format('  WARNING: atom ',i2,1x,a,' has not 1 bond',/)
  106 format('  WARNING: atom ',i2,1x,a,' has not 4 bonds',/)
  107 format('  WARNING: atom ',i2,1x,a,' has not 3 bonds',/)
  108 format('  WARNING: atom ',i2,1x,a,' has not 2 bonds',/)
  109 format('  WARNING: atom ',i2,1x,a,' has not 1 bond',/)
c
      return
      end

c----------------------------------------------------------------------c
      Subroutine Bond_Ord (IA,IB,RAB,BondOr,ipbo)
c----------------------------------------------------------------------c

c     Compute the bond order between atoms of atomic numbers IA
c     and IB, separated by RAB. Generate connectivity based on 
c     bond distances alone.  The criteria is whether the distances
c     is no more than 30% longer than the sum of covalent radii of
c     involved atoms. All bond types are determined using Pauling
c     bond orders.

      implicit Real*8 (A-H,O-Z)
      Save Pt2,Pt3,Pt5
      Data Pt2,Pt3,Pt5/0.2d0,0.3d0,0.5d0/
c
      if(ib.eq.0) then
        RAB = RCov97(IA,IB)
        return
      endif
      IPBO   = 0
      R1AB   = RAB
      R0AB   = RCov97(IA,IB)
      BondOr = Exp((R0AB-R1AB)/Pt3)
      if(BondOr.ge.Pt2) Then
        IBondO = Int(BondOr+Pt5)
        IBondO = Max(IBondO,1)
        IBondO = Min(IBondO,3)
        IPBO   = IBondO
      endif
c
      return
      end

c----------------------------------------------------------------------c
      Function RCov97 (IA,IB)
c----------------------------------------------------------------------c

c     This function returns an estimated covalent bond distance (in Ang)
c     between atoms of atomic numbers IA and IB. Setting IB to 0 returns
c     the covalent radius of IA. Parameters for atoms heavier than At are
c     taken from the UFF force field (A.K.Rappe',C.J.Casewit,K.S.Colwell,
c     W.A.Goddard III and W.M.Skiff, J.Am.Chem.Soc. 114,10024 (1992)).

      Implicit Real*8(A-H,O-Z)
      Parameter (MaxI=104)
      Dimension Rii(0:MaxI)
      Save Rii
      Data Rii/0.0d0,
     &  0.354D0,0.849D0,
     &  1.336D0,1.010D0,0.838D0,0.757D0,0.700D0,0.658D0,0.668D0,0.920D0,
     &  1.539D0,1.421D0,1.244D0,1.117D0,1.101D0,1.064D0,1.044D0,1.032D0,
     &  1.953D0,1.761D0,
     &  1.513D0,1.412D0,1.402D0,1.345D0,1.382D0,
     &  1.300D0,1.241D0,1.164D0,1.302D0,1.193D0,
     &                  1.260D0,1.197D0,1.211D0,1.190D0,1.192D0,1.147D0,
     &  2.260D0,2.052D0,
     &  1.698D0,1.564D0,1.473D0,1.467D0,1.322D0,
     &  1.478D0,1.332D0,1.338D0,1.386D0,1.403D0,
     &                  1.459D0,1.398D0,1.407D0,1.386D0,1.382D0,1.267D0,
     &  2.570D0,2.277D0,
     &  1.943D0,1.841D0,1.823D0,1.816D0,1.801D0,1.780D0,1.771D0,
     &  1.735D0,1.732D0,1.710D0,1.696D0,1.673D0,1.660D0,1.637D0,
     &  1.671D0,1.611D0,1.511D0,1.392D0,1.372D0,
     &  1.372D0,1.371D0,1.364D0,1.262D0,1.340D0,
     &                  1.518D0,1.459D0,1.512D0,1.500D0,1.545D0,1.420D0,
     &  2.880D0,2.512D0,1.983D0,1.721D0,1.711D0,1.684D0,1.666D0,
     &  1.657D0,1.660D0,1.801D0,1.761D0,1.750D0,1.724D0,1.712D0,
     &  1.689D0,1.679D0,1.698D0,
     &  1.850D0/
c
      RCov97 = Rii(Min(Max(IA,0),MaxI)) + Rii(Min(Max(IB,0),MaxI))
c
      return
      end

c----------------------------------------------------------------------c
      Function At_Symb (iz)
c----------------------------------------------------------------------c

c     given the atomic number, returns the atomic label

      implicit real*8(a-h,o-z)
      parameter (ntype=86)
      character symb(0:ntype)*2,At_Symb*5
      data symb/   'X '
     &  ,'H ','He'
     &  ,'Li','Be','B ','C ','N ','O ','F ','Ne'
     &  ,'Na','Mg','Al','Si','P ','S ','Cl','Ar'
     &  ,'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn'
     &                      ,'Ga','Ge','As','Se','Br','Kr'
     &  ,'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd'
     &                      ,'In','Sn','Sb','Te','I ','Xe'
     &  ,'Cs','Ba','La'
     &       ,'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy'
     &                           ,'Ho','Er','Tm','Yb','Lu'  !! Lu
     &                 ,'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg'
     &                      ,'Tl','Pb','Bi','Po','At','Rn'/
c
      if(iz.lt.0 .or. iz.gt.ntype) then
        write(6,'(//" ERROR in At_Symb: wrong Z",f10.3,//)') iz
        stop ' ERROR in At_Symb: wrong Z'
      endif
c
      At_Symb=' '
      At_Symb(1:2)=Symb(iz)
      return
      end

c----------------------------------------------------------------------c
      Function At_numb (name)
c----------------------------------------------------------------------c

c     given the atomic label, returns the atomic number

      implicit real*8(a-h,o-z)
      parameter (ntype=86)
      character symb(0:ntype)*2,name*(*),aa2*2,bb2*2
      data symb/   'X ','H ','He'
     &  ,'Li','Be','B ','C ','N ','O ','F ','Ne'
     &  ,'Na','Mg','Al','Si','P ','S ','Cl','Ar'
     &  ,'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn'
     &                      ,'Ga','Ge','As','Se','Br','Kr'
     &  ,'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd'
     &                      ,'In','Sn','Sb','Te','I ','Xe'
     &  ,'Cs','Ba','La'
     &       ,'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy'
     &                           ,'Ho','Er','Tm','Yb','Lu'  !! Lu
     &                 ,'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg'
     &                      ,'Tl','Pb','Bi','Po','At','Rn'/
c
      aa2=name(1:2)
      Call Via_Numer (aa2)
      Call LowCas (aa2)
      do i=1,ntype
        bb2=symb(i)
        Call LowCas (bb2)
        if(aa2.eq.bb2) goto 20
      enddo
      write(6,'(//" ERROR in At_numb: wrong name ",a,//)') name
      stop ' ERROR in At_numb: wrong atom name'
   20 At_numb=i
      return
      end

c----------------------------------------------------------------------c
      Subroutine Value_1IC (geom,maxat,i,j,k,l,val,st1)
c----------------------------------------------------------------------c

c     computes value of the IC defined by atoms i,j,k,l
c
c     k=0     distance (int units of c)
c     k.ne.0  angle in radians, if st1='d' angle in degrees

      implicit double precision (a-h,o-z)
      character*1 st1
      dimension geom(maxat,*)
      
      if(k.le.0) then
        val=dsqrt( (geom(i,1)-geom(j,1))**2
     &            +(geom(i,2)-geom(j,2))**2
     &            +(geom(i,3)-geom(j,3))**2 )
      else
        Call Value_Ang (geom,maxat,i,j,k,l,val,st1)
      endif
      return
      end

c----------------------------------------------------------------------c
      Subroutine Value_Ang (c,maxat,i,j,k,l,val,st1)
c----------------------------------------------------------------------c

c     computes the angle between three atoms (in degrees)
c           or the dihedral between four (in degrees)           

c           val is in radians, if st1='d' val is in degrees

      implicit double precision (a-h,o-z)
      data rad2deg /57.2957795130823229d0/
      character*1 st1
      dimension c(maxat,*)
c
      if((i*j*k).le.0) goto 80
      if(l.le.0) Call Angolo (c,maxat,i,j,k,val)
      if(l.gt.0) Call Diedro (c,maxat,i,j,k,l,val)
      if(st1.eq.'d' .or. st1.eq.'D') val=val*rad2deg
      return
   80 write(6,82) i,j,k,l
      stop
   82 format(//'  ERROR in Value_Ang: i j k l',4i4,//)
      end

c----------------------------------------------------------------------c
      Subroutine Angolo (c,maxat,i,j,k,val)
c----------------------------------------------------------------------c

c     computes the angle between three atoms (in degrees)

      implicit double precision (a-h,o-z)
      dimension c(maxat,*)
c
      x1=c(i,1)-c(j,1)
      y1=c(i,2)-c(j,2)
      z1=c(i,3)-c(j,3)
      x2=c(k,1)-c(j,1)
      y2=c(k,2)-c(j,2)
      z2=c(k,3)-c(j,3)
      b =dsqrt(x2*x2 + y2*y2 + z2*z2)
      a =dsqrt(x1*x1 + y1*y1 + z1*z1)
      ab=      x1*x2 + y1*y2 + z1*z2
      aabb=a*b
      if(dabs(aabb).lt.1.d-09) goto 100
      coseno=ab/aabb
      if(coseno.gt.+1.0d0) coseno=+1.d0
      if(coseno.lt.-1.0d0) coseno=-1.d0
c
      v1=y1*z2-z1*y2
      v2=z1*x2-x1*z2
      v3=x1*y2-y1*x2
      seno=dsqrt(v1*v1+v2*v2+v3*v3)/aabb
      if(seno.gt.+1.0d0) seno=+1.d0
      if(seno.lt.-1.0d0) seno=-1.d0
      val=datan2(seno,coseno)
c
      return
  100 write(6,*) ' ERROR in angolo: coincident atoms !!!'
      stop
      end

c----------------------------------------------------------------------c
      Subroutine Diedro (c,maxat,i,j,k,l,val)
c----------------------------------------------------------------------c

c     computes the dihedral betewwn four (in degrees)           

      implicit double precision (a-h,o-z)
      dimension c(maxat,*)
      tol=1.d-06
c
      xjk=c(j,1)-c(k,1)
      yjk=c(j,2)-c(k,2)
      zjk=c(j,3)-c(k,3)
      a=dsqrt(xjk**2 + yjk**2 + zjk**2)
      if(a.lt.tol) then
        write(6,'(//" ERROR in DIEDRO: atoms 2 3 too much close")')
        val = 0.d0
        return
      endif
      xjk=xjk/a
      yjk=yjk/a
      zjk=zjk/a
c
      xji=c(j,1)-c(i,1)
      yji=c(j,2)-c(i,2)
      zji=c(j,3)-c(i,3)
      scal=xji*xjk+yji*yjk+zji*zjk
      xji=xji-xjk*scal
      yji=yji-yjk*scal
      zji=zji-zjk*scal
      a=dsqrt(xji**2 + yji**2 + zji**2)
      if(a.lt.tol) then
        write(6,'(//" ERROR in DIEDRO: atoms 1 2 too much close")')
        val = 0.d0
        return
      endif
      xji=xji/a
      yji=yji/a
      zji=zji/a
c
      xkl=c(k,1)-c(l,1)
      ykl=c(k,2)-c(l,2)
      zkl=c(k,3)-c(l,3)
      a=dsqrt(xkl**2 + ykl**2 + zkl**2)
      if(a.lt.tol) then
        write(6,'(//" ERROR in DIEDRO: atoms 3 4 too much close")')
        val = 0.d0
        return
      endif
c
      scal=xkl*xjk+ykl*yjk+zkl*zjk
      xkl=xkl-xjk*scal
      ykl=ykl-yjk*scal
      zkl=zkl-zjk*scal
      a=dsqrt(xkl**2 + ykl**2 + zkl**2)
      if(a.lt.tol) then
        write(6,'(//" ERROR in DIEDRO: three atoms are colinear")')
        val = 0.d0
        return
      endif
      xkl=xkl/a
      ykl=ykl/a
      zkl=zkl/a
c
      coseno=xji*xkl+yji*ykl+zji*zkl
      seno  =xjk*(yji*zkl-ykl*zji)+yjk*(xkl*zji-xji*zkl)+
     &      zjk*(xji*ykl-yji*xkl)
      val=-atan2(seno,coseno)
c
      return
  140 write(6,144) '1',i,(c(i,iii),iii=1,3),'2',j,(c(j,iii),iii=1,3)
     &            ,'3',k,(c(k,iii),iii=1,3),'4',l,(c(l,iii),iii=1,3)
  144 format(/'  atoms for the dihedral',4(/3x,a,i5,3f14.7),//)
      stop
      end

c----------------------------------------------------------------------c
      Subroutine Curl (Va,Vb,Vab)
c----------------------------------------------------------------------c

c     computes vector product  Vab = Va x Vb

      implicit real*8(a-h,o-z)
      dimension Va(*),Vb(*),Vab(*)
c
      Vab(1)=Va(2)*Vb(3)-Va(3)*Vb(2)
      Vab(2)=Va(3)*Vb(1)-Va(1)*Vb(3)
      Vab(3)=Va(1)*Vb(2)-Va(2)*Vb(1)
      return
      end

c----------------------------------------------------------------------c
      Subroutine Atom_Count (nat,geom,aname,zord,maxat)
c----------------------------------------------------------------------c

c     given the atomic numbers, creates the atomic labels
c     with consecutive numbers after the symbol.

c     zord   .false.   numerazione interna per ciascun tipo di atomo
c            .true.    numerazione progressiva, ignora il tipo di atomo

c     Warning: if a '*' is already present in /aname/ it is preserved

      implicit real*8(a-h,o-z)
      parameter (max=130)
      dimension geom(maxat,*),inum(0:max)
      character str5*5,aname(*)*(*),symbi*2,At_Symb*5
      logical zord,znumb
      znumb=.true.
c
      tol=1.d-08
      la=Lenw(aname(1))
      do i=0,max
        inum(i)=0
      enddo
c
      do i=1,nat
        iz=geom(i,4)+tol
        if(iz.lt.0) then
          geom(i,4)=0.0d0
          iz=geom(i,4)+tol
        endif
        inum(iz)=inum(iz)+1
        ii=inum(iz)
        if(zord) ii=i
        symbi=At_Symb(iz)
        str5=symbi(1:2)
        if(znumb) write(str5(3:5),'(i3)') ii
        Call Compact (str5,lw)
        istar=index(aname(i),'*')
        ll=Lenw(str5)
        if(ll.lt.5 .and. istar.gt.0) str5(ll+1:ll+1)='*'
        aname(i)=str5
      enddo
   22 format(i4,4x,a,f7.1,3f14.7)
c
      return
      end

c----------------------------------------------------------------------c
      Function Atom_Mass (iat)
c----------------------------------------------------------------------c

c     given the atomic numbers, returns the atomic mass

      implicit real*8(a-h,o-z)
      dimension amas(107)
      data  amas /  1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,
     & 10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,
     & 20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,
     & 30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,
     & 40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,
     & 54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,
     & 65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,
     & 79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,
     & 91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,
     & 102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,
     & 118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,
     & 132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,
     & 183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,
     & 196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0,
     & 24*0.000D0 /
c
      if(iat.lt.1 .or. iat.gt.107) then
        Atom_Mass=0.0d0
      else
        Atom_Mass=amas(iat)
      endif
c
      return
      end

c----------------------------------------------------------------------c
      Subroutine Inertia (c,nat,baric,amom,vmom,maxat)
c----------------------------------------------------------------------c

c     compute the inertia moment and principal axes

c     c,nat   I  molecular geometry (not changed)
c     baric   O  center of mass
c     amom    O  inertia moments amd their vectors (princ axes)


c     WARNING: from 16/9/06  the masses are to be stored in c(i,5) !!
c              This allows for a larger flexibility (United Atom ...)

      implicit real*8(a-h,o-z)
      parameter (zero=0.0d+00)
      dimension c(maxat,*),baric(*),amom(*),vmom(*)
      dimension tmom(6),trot(3,3)

c     compute the center of mass

      totwt=zero
      baric(1)=zero
      baric(2)=zero
      baric(3)=zero
c
      do i=1,nat
        wi=c(i,5)
        if(wi.lt.0.01d0) goto 40
        totwt=totwt+wi
        do j=1,3
          baric(j)=baric(j)+wi*c(i,j)
        enddo
      enddo
      do i=1,3
        baric(i)=baric(i)/totwt
      enddo

c     compute the moment of inertia tensor

      xx=zero
      yy=zero
      zz=zero
      xy=zero
      xz=zero
      yz=zero
      do i=1,nat
        wi=c(i,5)
        xc=c(i,1)-baric(1)
        yc=c(i,2)-baric(2)
        zc=c(i,3)-baric(3)
        xx=xx+wi*(yc*yc+zc*zc)
        yy=yy+wi*(xc*xc+zc*zc)
        zz=zz+wi*(xc*xc+yc*yc)
        xy=xy-wi*xc*yc
        xz=xz-wi*xc*zc
        yz=yz-wi*yc*zc
      enddo
c
      tmom(1) = xx
      tmom(2) = xy
      tmom(3) = yy
      tmom(4) = xz
      tmom(5) = yz
      tmom(6) = zz
      ierr = 0
      n=3
      Call Eigen (tmom,vmom,n)
      amom(1)=tmom(1)
      amom(2)=tmom(2)
      amom(3)=tmom(3)
c
      return
   40 write(6,42)
   42 format(//' EEROR in INERTIA : one or more mass is too low',///)
      stop
      end
