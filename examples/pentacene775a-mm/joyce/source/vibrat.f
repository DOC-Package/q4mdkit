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

c =================================================================
c
      Subroutine Comp_Vibr (FCC,U,F,E,Freq,C,geom,scelta,zpic,zgrah
     #                     ,ncint,nat,nvib,labint
     #                     ,max3n,maxat)
c
c  compute the normal vibrational modes from the hessian matrix
c  given in cartesian coordinates (bohr). Internal coordinates are also
c  given and used, warning: they must be not redundant IC !!!
c
c  FCC     Force constant matrix in the cartesian coordinate basis
c          Units: (Hartree/bohr^2), (Hartree/bohr*rad), (Hartree/rad^2)
c  U       transformation matrix from internal to cartesian coordinates.
c          First index: internal, second index: cartesian (bohr)
c  F       (out) force const. matrix in internal coordinates (if zpic=t)
c  E       Hessian in the NM basis (only diagonal terms)
c  Freq    Computed frequencies
c  C       MN vectors in the mass weighted CC (MWCC)
c  geom    current molecular geometry in CC
c  scelta  label for printouts
c  zpic    T = try to project NMs onto ICs
c          F = do not project NMs onto ICs
c  zgrah   T = compute the vibr. modes by working in a space orthogonal to
c              the vibr+rot modes (use Graham-Schmidt orthogonalization)
c              (recommended choice)
c          F = compute the vibr. modes by analyzing the null eigenvectors
c              (not recommended)
c  ncint   n. of internal coordinates
c  nat     n. of atoms
c  nvib    n. of vibrational normal modes (usually 3*nat-6)
c  labint  label of the internal coordinates
c  max3n   dimension of two-dim arrays in the caller
c  maxat   dimension of two-dim array /geom/ in the caller
c
      implicit real*8(a-h,o-z)
      parameter (maxf=2000, maxq=maxf*maxf)
      dimension FCC(max3n,*),U(max3n,*),Freq(*),geom(maxat,*)
     &         ,F(*),E(*),C(max3n,*),Utrot(6,maxf)
      dimension V(maxq),D(maxq),UU(maxf,maxf),TT(maxf,maxf)
     &         ,wmass(maxf)
     &         ,ndsym(5),DD(maxf,maxf),VV(maxf,maxf),iadd(maxf)
      character labint(*)*(*),labc(maxf)*4,scelta*(*),title*50
      logical zpic,zfina,zgrah
      equivalence (V(1),VV(1,1)),(D(1),DD(1,1))
      common /printout/ levprt
      data har2cm1 /219474.6156d0/    !! 1 hartree = 221474 1/cm
      data amu/1822.8880d0/  !!   1 amu = 1822.888 electron mass
c
      n=3*nat
      ncart=n
      Call Check_Dim (n,maxf,'COMP_NORM_MOD: increase max')
      ncol=10
      thrs=0.01d0
c
      m=0
      do i=1,nat
        aa=geom(i,5)
        wmass(m+1)=dsqrt(aa*amu)
        wmass(m+2)=wmass(m+1)
        wmass(m+3)=wmass(m+1)
        m=m+3
      enddo
c
c  compute cartesian labels for print purposes
c
      j=0
      do i=1,nat
        j=j+1
        write(labc(j),'(i3,"x")') i
        j=j+1
        write(labc(j),'(i3,"y")') i
        j=j+1
        write(labc(j),'(i3,"z")') i
      enddo
c
      lw=Lenw(scelta)
      nb=(50-lw)/2
      title=' '
      title(nb+1:)=scelta
      nsym=1
      write(6,15) title
   15 format(
     &   /10x,'====================================================='
     &  ,/10x,'   N O R M A L   V I B R A T I O N A L   M O D E S'
     &  ,/10x,a
     &  ,/10x,'=====================================================',/)
c
c  print input data
c
      if(levprt.ge.2) then
        ij=0
        do i=1,n
          do j=1,i
            ij=ij+1
            D(ij)=FCC(i,j)
          enddo
        enddo
        title='FC matrix (H) in cart displac.s'
        ndsym(1)=n
        Call Pr_Tri (D,nsym,ndsym,title,ncol)
      endif
c
c  transforms F.C. matrix form cart. to mass weighted cart. coord.
c  and save it in triangular form
c
      write(6,'( 10x,"Compute the Mass weighted Hessian")')
      do i=1,n
        do j=1,n
          VV(i,j)=FCC(i,j)/wmass(i)/wmass(j)
        enddo
      enddo
      title='HESSIAN matrix (mH) in mass weig cart coord.s'
      ndsym(1)=n
      if(levprt.ge.3) then
        ij=0
        do j=1,n
          do i=1,j
            ij=ij+1
            D(ij)=VV(i,j)*1000.d0
          enddo
        enddo
        Call Pr_Tri (D,nsym,ndsym,title,ncol)
      endif
c
c  compute the translational and rotational modes in MWCC
c  they form an orthonormal 6-dim vector space
c
      Call U_for_Trot (geom,nat,Utrot,ncart,maxat)
c
c ---------------------------------------------------------
c --- use graham-schmidt orthogonalization to separate  ---
c --- the transl+rot modes from the vibrational modes   ---
c ---------------------------------------------------------
c
      if(.not.zgrah) goto 100
c
c   transform the Hessian in a vector space given by
c   3N-6 vectors, orthogonal to the transl+rot modes
c
      do i=1,6   !! copy transl+rot
        do k=1,ncart
          UU(k,i)=Utrot(i,k)
        enddo
      enddo
      do i=7,ncart+6
        do k=1,ncart
          UU(k,i)=0.0d0
        enddo
          UU(i-6,i)=1.0d0
      enddo
      Call Graham2 (UU,ncart+6,ncart,maxf)
      nvib=ncart-6
      do i=1,nvib
        do k=1,ncart
          UU(k,i)=UU(k,i+6)
        enddo
      enddo
      zfina=.false.
      Call Trsml2 (VV,UU,DD,nvib,ncart,maxf,zfina)
c
c   save the projected Hessian in MWCC in a triangular array
c
      ij=0
      do i=1,nvib
        do j=1,i
          ij=ij+1
          F(ij)=DD(i,j)
        enddo
      enddo
c
c  diagonalize the mass weighted Force Constant matrix (cart coord)
c
      write(6,'(10x,"Diagonalize the Mass weighted Hessian")')
c
      Call Eigen (F,D,nvib)
c
c  retrieve the original vectors in the MWCC space
c
      do i=1,nvib
        ii=(i-1)*nvib
        amax=0.d0
        do j=1,n
          aaa=0.0d0
          do k=1,nvib
            aaa=aaa+UU(j,k)*D(ii+k)
          enddo
          C(j,i)=aaa
          if(dabs(aaa).gt.amax) then
            amax=dabs(aaa)
            jmax=j
          endif
        enddo
        if(C(jmax,i).lt.0.0d0) then
          do j=1,n
            C(j,i)=-C(j,i)
          enddo
        endif
      enddo
c
c  compute frequencies and check
c
      write(6,56) scelta,(i,1.d3*F(i),i=1,nvib)
   56 format(/10x,'Eigenvalues of the (M-1/2)*F*M(-1/2) matrix (mH)',
     &        a,/10x,(5(i6,g11.4)))
c
      if(levprt.ge.3) write(6,58) scelta
   58 format(/'  All the NORMAL MODES in mass-weighted cartes ',a)
      do i=1,nvib
        if(F(i).ge.0.0d0) frau=dsqrt(F(i))
        if(F(i).lt.0.0d0) frau=-dsqrt(-F(i))
        Freq(i)=frau*har2cm1
        E(i)=F(i)
      enddo
      write(6,64) scelta,(i,Freq(i),i=1,nvib)
   64 format(/10x,'    Frequencies in 1/cm ',2x,a,/(5(i6,f10.3)))
      goto 140
c
c ---------------------------------------------------------
c --- use Hessian matrix projected outside the space    ---
c --- spanned by the transl+rot modes. The vibrational  ---
c --- modes are recognized by their projection onto     ---
c --- the transl+rot modes                              ---
c ---------------------------------------------------------
c
c
c  project the FC out of the transl+rot modes
c
  100 continue
      write(6,'(/8x,"Project FC matrix outside the transl,rot modes")')
      do j=1,n
        do i=1,n
          DD(i,j)=VV(i,j)
          UU(i,j)=0.0d0
          if(i.eq.j) UU(i,j)=1.0d0
          do itr=1,6
            UU(i,j)=UU(i,j)-Utrot(itr,i)*Utrot(itr,j)
          enddo
        enddo
      enddo
      zfina=.false.
      Call Trsml2 (DD,UU,VV,ncart,ncart,maxf,zfina)
c
c   save the projected Hessian in MWCC in a triangular array
c
      ij=0
      do i=1,n
        do j=1,i
          ij=ij+1
          F(ij)=DD(i,j)
        enddo
      enddo
      if(levprt.ge.2) then
        do ij=1,n*(n+1)/2
          D(ij)=F(ij)*1000.d0
        enddo
        title='(1-P)FC(1-P) matrix (mH) in mass weig cart displac.s'
        Call Pr_Tri (D,nsym,ndsym,title,ncol)
      endif
c
c  diagonalize the mass weighted Force Constant matrix (cart coord)
c
      write(6,'(/8x,"Diagonalize the Mass weighted Hessian")')
c
      Call Eigen (F,D,n)
c
c  compute frequencies and check
c
      write(6,116) scelta,(i,1.d3*F(i),i=1,n)
  116 format(/'    Eigenvalues of the (M-1/2)*F*M(-1/2) matrix (mH) ',a
     #      ,/(5(i6,g11.4)))
c
      if(levprt.ge.3) write(6,118) scelta
  118 format(/'  All the NORMAL MODES in mass-weighted cartes ',a)
      do i=1,n
        if(F(i).ge.0.0d0) frau=dsqrt(F(i))
        if(F(i).lt.0.0d0) frau=-dsqrt(-F(i))
        Freq(i)=frau*har2cm1
      enddo
      nn=n
      if(levprt.ge.3) Call Pr_Orb (D,Freq,nn,n,labc,thrs,ncol)
c
      Call Sort_Vibr (F,D,wmass,Utrot,Freq,V,iadd,n,nvib)
c
      write(6,120) scelta,(i,Freq(i),i=1,n)
  120 format(/'    Frequencies in 1/cm ',2x,a,/(5(i6,f10.3)))
      nf=n
      write(6,124) scelta
  124 format(/'  VIBRATIONAL NORMAL MODES in mass-weighted cartes ',a)
      if(levprt.ge.2) Call Pr_Orb (D,Freq,nvib,n,labc,thrs,ncol)
c
c  save the eigenvectors on the C 2-dim array,
c  the zero energy modes are on the higher positions
c
      inm=0
      do nm=1,n
        E(nm)=F(nm)
        do i=1,n
          inm=inm+1
          C(i,nm)=D(inm)
        enddo
      enddo
c
c ---------------------------------------------------------
c ---------------------------------------------------------
c ---      analysis of the vibrational modes            ---
c ---------------------------------------------------------
c ---------------------------------------------------------
c
c  compute vectors proportional to cartesian displacements
c           X(t) =  M(-1/2) * C * Q(t)
c
  140 continue
      j=0
      do nm=1,nvib
        aa=0.0d0
        j0=j
        do i=1,n
          j=j+1
          V(j)=C(i,nm)/wmass(i)
          aa=aa+V(j)**2
        enddo
        aa=dsqrt(aa)
        j=j0
        do i=1,n
          j=j+1
          V(j)=V(j)/aa
        enddo
      enddo
      nf=n
      if (levprt .gt. 0 ) write(6,132) scelta
  132 format(/'  Normalized VIBR NORMAL MODES: cart. displac.s ',a)
      if(levprt.ge.3) Call Pr_Orb (V,Freq,nvib,nf,labc,thrs,ncol)
c
c  compute vectors proportional to displacements of internal coord.s
c
      inm=0
      do nm=1,nvib
        anorm=0.0d0
        do int=1,ncint
          aa=0.d0
          knm=n*(nm-1)
          do k=1,n
            knm=knm+1
            aa=aa+U(int,k)*V(knm)
          enddo
          inm=inm+1
          D(inm)=aa
          anorm=anorm+aa*aa
        enddo
        if(anorm.gt.1.d-03) then
          anorm=dsqrt(anorm)
          do k=inm-ncint+1,inm
            D(k)=D(k)/anorm
          enddo
        endif
      enddo
      write(6,135) scelta
  135 format(/10x,'  VIBRATIONAL NORMAL MODES: Int.Coord. displac.s ',a)
      Call Pr_Orb (D,Freq,nvib,ncint,labint,thrs,ncol)
c
c  project normal modes onto the internal coordinates
c
      inm=0
      do nm=1,nvib
        anorm=0.0d0
        do int=1,ncint
          aa=0.d0
          do k=1,n
            aa=aa+U(int,k)*C(k,nm)*wmass(k)
          enddo
          inm=inm+1
          D(inm)=aa
        enddo
      enddo
c
c  skip projection if zpic=false
c
      if(.not.zpic) goto 200
c
c  compute the inverse of the U*(U+) matrix
c
      do jnt=1,ncint
        do int=1,jnt
          fij=0.0d0
          do k=1,n
            fij=fij+U(jnt,k)*U(int,k)
          enddo
          TT(int,jnt)=fij
          TT(jnt,int)=fij
        enddo
      enddo
      Call Invert (TT,ncint,max3n,ier)
c
c  in the case of redundant intern. coord.s : skip projection
c
      if(ier.ne.0) then
        write(6,'(//" WARNING: Redund ICs found: no NM projection",/)')
        goto 200
      endif
c
c  if the intern. coord.s are not redundant:
c  project the normal modes onto internal coordinates
c
      do i=1,nvib*ncint
        V(i)=D(i)
      enddo
      inm=0
      do nm=1,nvib
        do int=1,ncint
          aa=0.d0
          knm=ncint*(nm-1)
          do jnt=1,ncint
            knm=knm+1
            aa=aa+TT(int,jnt)*V(knm)
          enddo
          inm=inm+1
          D(inm)=aa
        enddo
      enddo
c
c  if the intern. coord.s are not redundant:
c  transforms the force constant matrix form cartesian norm. modes
c  to internal coordinates
c
      inm=0
      do mu=1,ncint
        do nu=1,mu
          aa=0.d0
          muk=mu
          nuk=nu
          do k=1,nvib
            aa=aa+D(muk)*E(k)*D(nuk)
            muk=muk+ncint
            nuk=nuk+ncint
          enddo
          inm=inm+1
          F(inm)=aa
        enddo
      enddo
      ndsym(1)=ncint
      title='Force Const matr in Int.C (recomputed) '//scelta
      Call Pr_Tri (F,nsym,ndsym,title,ncol)
c
  200 write(6,'(/10x,"End Normal Modes Calculation")')
      return
      end

      Subroutine Comp_Eckart (Freq,C,geom,scelta
     #                       ,nat,nvib,max3n,maxat)
c
c  compute vibrational space basis in MWCC as the orthogonal space
c  wrt translations and rotations
c
c  Freq    Frequencies set to zero (to avoid using the random values in memory
c  C       MN vectors in the mass weighted CC (MWCC)
c  geom    current molecular geometry in CC
c  scelta  label for printouts
c  nat     n. of atoms
c  nvib    n. of vibrational normal modes (usually 3*nat-6)
c  max3n   dimension of two-dim arrays in the caller
c  maxat   dimension of two-dim array /geom/ in the caller
c
      implicit real*8(a-h,o-z)
      parameter (maxf=2000, maxq=maxf*maxf)
      dimension Freq(*),geom(maxat,*)
     #         ,C(max3n,*),Utrot(6,maxf)
      dimension V(maxq),D(maxq),UU(maxf,maxf),TT(maxf,maxf),wmass(maxf)
     #         ,ndsym(5),DD(maxf,maxf),VV(maxf,maxf),iadd(maxf)
      character labc(maxf)*4,scelta*(*),title*50
      logical zpic,zfina,zgrah
      equivalence (V(1),VV(1,1)),(D(1),DD(1,1))
      common /printout/ levprt
      data har2cm1 /219474.6156d0/    !! 1 hartree = 221474 1/cm
      data amu/1822.8880d0/  !!   1 amu = 1822.888 electron mass
c
      n=3*nat
      ncart=n
      Call Check_Dim (n,maxf,'COMP_ARB_VIBR: increase max')
      ncol=10
      thrs=0.01d0
c
      m=0
      do i=1,nat
        aa=geom(i,5)
        wmass(m+1)=dsqrt(aa*amu)
        wmass(m+2)=wmass(m+1)
        wmass(m+3)=wmass(m+1)
        m=m+3
      enddo
c
c  compute cartesian labels for print purposes
c
      j=0
      do i=1,nat
        j=j+1
        write(labc(j),'(i3,"x")') i
        j=j+1
        write(labc(j),'(i3,"y")') i
        j=j+1
        write(labc(j),'(i3,"z")') i
      enddo
c
      lw=Lenw(scelta)
      nb=(50-lw)/2
      title=' '
      title(nb+1:)=scelta
      nsym=1
      write(6,15) title
   15 format(
     #  //12x,'====================================================='
     #  ,/12x,'====================================================='
     #  ,/12x,'     A R B I T R A R Y   V I B R A T I O N A L       '
     #  ,/12x,'                    B A S I S                        '
     #  ,/12x,a
     #  ,/12x,'====================================================='
     #  ,/12x,'=====================================================',/)
c
c
c  compute the translational and rotational modes in MWCC
c  they form an orthonormal 6-dim vector space
c
      Call U_for_Trot (geom,nat,Utrot,ncart,maxat)
c
c ---------------------------------------------------------
c --- use graham-schmidt orthogonalization to separate  ---
c --- the transl+rot modes from the vibrational modes   ---
c ---------------------------------------------------------
c
c   transform the Hessian in a vector space given by
c   3N-6 vectors, orthogonal to the transl+rot modes
c
      do i=1,6   !! copy transl+rot
        do k=1,ncart
          UU(k,i)=Utrot(i,k)
        enddo
      enddo
      do i=7,ncart+6
        do k=1,ncart
          UU(k,i)=0.0d0
        enddo
          UU(i-6,i)=1.0d0
      enddo
      Call Graham2 (UU,ncart+6,ncart,maxf)
      nvib=ncart-6
      do i=1,nvib
        do k=1,ncart
          UU(k,i)=UU(k,i+6)
        enddo
        Freq(i) = 0.d0
      enddo

c
c  retrieve the original vectors in the MWCC space
c
      do i=1,nvib
        ii=(i-1)*nvib
        amax=0.d0
        do j=1,n
          aaa=UU(j,i)
          C(j,i)=aaa
          if(dabs(aaa).gt.amax) then
            amax=dabs(aaa)
            jmax=j
          endif
        enddo
        if(C(jmax,i).lt.0.0d0) then
          do j=1,n
            C(j,i)=-C(j,i)
          enddo
        endif
      enddo

c
c ---------------------------------------------------------
c ---------------------------------------------------------
c ---  analysis as if they were vibrational modes  ---
c ---------------------------------------------------------
c ---------------------------------------------------------
c
c  compute vectors proportional to cartesian displacements
c           X(t) =  M(-1/2) * C * Q(t)
c
      j=0
      do nm=1,nvib
        aa=0.0d0
        j0=j
        do i=1,n
          j=j+1
          V(j)=C(i,nm)/wmass(i)
          aa=aa+V(j)**2
        enddo
        aa=dsqrt(aa)
        j=j0
        do i=1,n
          j=j+1
          V(j)=V(j)/aa
        enddo
      enddo
      write(6,232) scelta
  232 format(/'  Normalized VIBR COORDINATES: cart. displac.s ',a)
      if(levprt.ge.3) Call Pr_Orb (V,Freq,nvib,nf,labc,thrs,ncol)


      write(6,'(/8x,"End Generation Of Arb Vibrational Space",//)')
      return
      end
c
c =================================================================
c
      Subroutine U_for_Trot (geom,nat,Utrot,ncart,maxat)
c
c  Compute the projections of the translational and rotational
c  modes onto the mass weighted cartesian coordinates (MWCC)
c
c        The transl and rotat modes are orthogonal !!!
c
c
      implicit real*8(a-h,o-z)
      dimension geom(maxat,*),Utrot(6,*)
     #         ,Va(3),Vb(3),Vab(3),cmass(3),W(6),UU(6,6)
     #         ,amom(9),vmom(9),baric(3)
      logical zdeb
      character*64 tit
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      common /printout/ levprt
      zdeb=.false.
c
c               transl+rot are expressed in CC
c     imode=1   rotat modes in the new (Wilson Decius Cross) manner
c               transl+rot are expressed in MassWeighted-CC
c
      imode=1   !!! WARNING use imode=1  This gives 6 orthonormal vectors !!
c
      ndim=6
      if (levptr .gt. 0 ) write(6,12)
   12 format(/10x,'Compute the Transl+Rot modes in the MW-CC.s')
c
c  compute the center of mass of the molecule
c
      totwt=0.0d0
      cmass(1)=0.0d0
      cmass(2)=0.0d0
      cmass(3)=0.0d0
      do i=1,nat
        wi=geom(i,5)
        totwt=totwt+wi
        do j=1,3
          cmass(j)=cmass(j)+wi*geom(i,j)
        enddo
      enddo
      cmass(1)=cmass(1)/totwt
      cmass(2)=cmass(2)/totwt
      cmass(3)=cmass(3)/totwt
      if (levptr .gt. 0) write(6,16)  (boh2ang*cmass(k),k=1,3)
   16 format(/8x,'The Center of mass (in Angs) is',3f10.5)
c
      if(imode.ne.0) goto 30
c
c--------------------------------------------------------
c--    translations and rotations (in CC, old manner)  --
c--------------------------------------------------------
c
      do itr=1,3
        iccw=mod(itr,3)
        do k=1,ncart
          Utrot(itr,k)=0.0d0
          if(mod(k,3).eq.iccw) Utrot(itr,k)=1.0d0
        enddo
      enddo
c
      do iro=1,3
        i=iro+3
        do k=1,3
          Vb(k)=0.0d0
        enddo
        Vb(iro)=1.0d0
        do iat=1,nat
          va(1)=(geom(iat,1)-cmass(1))
          va(2)=(geom(iat,2)-cmass(2))
          va(3)=(geom(iat,3)-cmass(3))
          Call Curl (Va,Vb,Vab)
          iat3=(iat-1)*3
          Utrot(i,iat3+1)=Vab(1)
          Utrot(i,iat3+2)=Vab(2)
          Utrot(i,iat3+3)=Vab(3)
        enddo
      enddo
      goto 40
c
c--------------------------------------------------------
c--       translations and rotations (in MW-CC)        --
c--       according to Wilson, Decius, Cross           --
c
c   note that the rotations are only orthogonal to each other
c   in the frame defined by the principal inertia axes
c
c--------------------------------------------------------
c
   30 continue
      Call Inertia (geom,nat,baric,amom,vmom,maxat)
      if (levprt .gt. 0) write(6,32) (amom(i),i=1,3)
      do i=1,3
         if (levprt .gt. 0) write(6,34) (vmom(j),j=i,9,3)
      enddo
   32 format(/10x,'principal inertia axes (au)',/(6x,3f12.5))
   34 format(5x,3f12.5)
c
      do itr=1,6
        do k=1,ncart
          Utrot(itr,k)=0.0d0
        enddo
      enddo
      do i=1,nat
        ii=i*3-3
        amw=dsqrt(geom(i,5))
        Utrot(1,ii+1)=amw
        Utrot(2,ii+2)=amw
        Utrot(3,ii+3)=amw
        do k=1,3
          Vb(k)=Geom(i,k)-cmass(k)
        enddo
        m=0
        do ir=4,6
          do k=1,3
            m=m+1
            Va(k)=vmom(m)
          enddo
          Call Curl (Va,Vb,Vab)
          do k=1,3
            Utrot(ir,ii+k)=Vab(k)*amw
          enddo
        enddo
      enddo
c
   40 continue
c
c--      normalization
c
      do itr=1,6
        aaa=0.0d0
        do k=1,ncart
          aaa=aaa+Utrot(itr,k)**2
        enddo
        aaa=dsqrt(aaa)
        do k=1,ncart
          Utrot(itr,k)=Utrot(itr,k)/aaa
        enddo
      enddo
c
c--      printout
c
      i=ncart
      n=6
      ncol=10
      tit='Matrix of transl+rot modes in CCs (bohr,rad) '
      if(imode.ne.0) 
     #tit='Matrix of transl+rot modes in MWCCs (bohr,rad) '
      if(levprt.ge.2 .or. zdeb) Call Pr_Rec2 (Utrot,n,i,tit,ncol,ndim)
c
c--     ortonormality check
c
      ermax=0.d0
      do j=1,n
        do i=1,j
          uij=0.0d0
          do k=1,ncart
            uij=uij+Utrot(i,k)*Utrot(j,k)
          enddo
          UU(i,j)=uij
          UU(j,i)=uij
          erij=dabs(uij)
          if(i.eq.j) erij=dabs(uij-1.0d0)
          if(erij.gt.ermax) then
            ii=i
            jj=j
            ermax=erij
          endif
        enddo
      enddo
      if (levprt .gt. 0 ) write(6,55) ii,jj,ermax
      if(zdeb) then
        n=6
        i=n
        ndim=6
        tit='Check Orthonormality of the transl+rot modes'
        Call Pr_Rec2 (UU,n,i,tit,ncol,ndim)
      endif
   55 format(/10x,'Max err in the U(transl+rot)',2i4,1pd12.3)
c
      return
      end
c
c =================================================================
c
      Subroutine Sort_Vibr (F,D,wmass,Utrot,Freq,V,iadd,n,nvib)
c
c  analyze the normal modes and sort the transl+rot modes, which
c  will be put in the first position of the F,D arrays.
c  This routine works even in the case of negative frequencies.
c
c
c  Utrot   transformation matrix from internal to cartesian coordinates.
c          First index: internal, second index: cartesian (bohr)
c  Freq    vibrational frequencies (out)
c  V       work array (n*n)
c  F,D     eigensolutions of the mass weighted hessian matrix M(-1/2)*H*M(-1/2)
c          At return these are rearranged: vibrat + (transl+rot)
c  wmass   triples of sqrt(atomic mass) for transformation in mass
c          weighted cartesian coordinates
c  Utrot   (6,CC) projection of the rot+transl in the CC
c  n       3*nat = n. of cartesian coordinates
c  nvib    n. of vibrational normal modes (usually 3*nat-6=n-6)
c
      implicit real*8(a-h,o-z)
      dimension F(*),D(*),wmass(*),Utrot(6,*),Freq(*),V(*),iadd(*)
      character isok(2)*5
      data har2cm1 /219474.6156d0/    !! 1 hartree = 221474 1/cm
      nrt=6
c
c   compute the projection of the normal modes onto the transl+rot modes
c
      j=0
      do nm=1,n
        ptr=0.0d0
        do itr=1,nrt
          pp=0.0d0
          j=(nm-1)*n
          do i=1,n
            j=j+1
            pp=pp+D(j)*Utrot(itr,i)
          enddo
          ptr=ptr+pp*pp
        enddo
        Freq(nm)=ptr
        iadd(nm)=nm
      enddo
      write(6,20) (i,Freq(i),i=1,n)
   20 format(/'    Projections of the NMs onto the Transl+Rot modes'
     #      ,/(5(i6,f10.2)))
c
c   sort the transl + rot modes on the basis of the computed projection
c   One more value is selected for test purposes
c
      do k=1,nrt+1
        aa=0.0d0
        do nm=1,n
          if(Freq(nm).gt.aa) then !! FREQ here is the projectiona on tr+rot
            aa=Freq(nm)
            inm=nm
          endif
        enddo
        iadd(k)=inm
        Freq(inm)=-Freq(inm)
      enddo
      ivib=inm
      evib=F(ivib)
      pvib=-Freq(ivib)
      if(pvib.lt.1.d-12) pvib=1.d-12
c
      do i=1,nrt-1
        do j=i+1,nrt
          if(iadd(i).gt.iadd(j)) then
            ii=iadd(i)
            iadd(i)=iadd(j)
            iadd(j)=ii
          endif
        enddo
      enddo
      write(6,'(8x,"transl+rot Modes are",6i4)') (iadd(k),k=1,nrt)
      write(6,'(8x,"max proj of vibr NM onto transl+rot",f8.2)') pvib
c
c   check the transl + rot modes
c
      erm=0.d0
      do i=1,nrt
        irt=iadd(i)
        Freq(irt)=-Freq(irt)
        prj=Freq(irt)
        eig=F(irt)
        isok(1)='WRONG'
        isok(2)='WRONG'
        itest=0
        if(prj.gt.0.8d0) then
          itest=itest+1
          isok(1)=' OK  '
        endif
        if(eig.lt.1.d-06) then
          itest=itest+1
          isok(2)=' OK  '
        endif
        if(itest.ne.2) write(6,24) irt,isok
        if(dabs(eig).gt.erm) erm=dabs(eig)
      enddo
      write(6,'(8x,"max error on null EigenValues",1pd12.4)') erm
   24 format(//10x,'WARNING: the transl + rot mode',i4,' is critical'
     #       ,/19x,'proj. onto transl+rot modes   > 0.8  ',4x,a
     #       ,/19x,'eigenvalue                    < 1.d-6',4x,a,/)
c
c   rearrange the eigensolutions
c
      do i=1,n*n
        V(i)=D(i)
      enddo
      do i=1,n
        Freq(i)=F(i)
      enddo
c
      nv=0
      iv=0
      nr=n-nrt
      ir=(nr-1)*n
      iupp=0
      do nm=1,n
        ilow=iupp+1
        iupp=iupp+n
        ivib=1
        do irt=1,nrt
          if(iadd(irt).eq.nm) ivib=0
        enddo
        if(ivib.gt.0) then
          nv=nv+1
          F(nv)=Freq(nm)
          do i=ilow,iupp
            iv=iv+1
            D(iv)=V(i)
          enddo
        else
          nr=nr+1
          F(nr)=Freq(nm)
          do i=ilow,iupp
            ir=ir+1
            D(ir)=V(i)
          enddo
        endif
      enddo
      nvib=nv
      nzeri=n-nvib
c
c  compute frequencies and check
c
      write(6,'(8x,"vibrations are moved to lower positions")')
      fconv=har2cm1
      do i=1,n
        if(F(i).ge.0.0d0) frau=dsqrt(F(i))
        if(F(i).lt.0.0d0) frau=-dsqrt(-F(i))
        Freq(i)=frau*fconv
      enddo
c
      return
      end
c
c =================================================================
c
      Subroutine Graham (U,n,ndim)
c
c  A set of vectors U(*,i) is orthonormalized by Graham-Schmidt method
c  for metric = indentity.
c  A special treatment is deserved for the last vector in the case
c  its norma (after orthogonalization) is very small.
c
      implicit real*8(a-h,o-z)
      dimension U(ndim,*)
      tol=1.d-08
c
      i=1
      Call Normaliz_C (U,n,i,anr,ndim)
      if(anr.lt.tol) goto 120
c
      do i=2,n
        Call Normaliz_C (U,n,i,anr,ndim)
        if(anr.lt.tol) goto 120
        do j=1,i-1
          aij=0.d0
          do k=1,n
            aij=aij+U(k,i)*U(k,j)
          enddo
          do k=1,n
            U(k,i)=U(k,i)-aij*U(k,j)
          enddo
        enddo
        Call Normaliz_C (U,n,i,anr,ndim)
        if(anr.lt.tol) then
          if(i.eq.n) goto 80
          goto 120
        endif
      enddo
c
c  special treatment for the last (critical) vector
c
   80 toln=1.d-20
      n1=n-1
      uknmax=0.0d0
      do k=1,n
        akk=0.0d0
        do j=1,n1
          akk=akk+U(k,j)**2
        enddo
        ukn=0.0d0
        if(akk.gt.1.0d0) akk=1.0d0
        ukn=dsqrt(1.d0-akk)
        U(k,n)=ukn
        if(ukn.gt.uknmax) then
          uknmax=ukn
          kk=k
        endif
      enddo
      Call Normaliz_C (U,n,i,anr,ndim)
      i=n
      if(anr.lt.tol) goto 110
c---  lasct vector sign choice
      do k=1,n
        aa=0.0d0
        do j=1,n
          aa=aa+U(k,j)*U(kk,j)
        enddo
        apiu=U(k,n)*U(kk,n)/2.0d0
        if(k.eq.kk) aa=aa-1.d0
        if(aa.gt.+apiu) U(k,n)=-U(k,n)
        if(aa.lt.-apiu) U(k,n)=-U(k,n)
      enddo
c--- Graham-Schmidt for the last vector
        do j=1,n-1
          aij=0.d0
          do k=1,n
            aij=aij+U(k,n)*U(k,j)
          enddo
          do k=1,n
            U(k,n)=U(k,n)-aij*U(k,j)
          enddo
        enddo
        i=n
        Call Normaliz_C (U,n,i,anr,ndim)
        if(anr.lt.tol) goto 120
c---
c---   orthonormalization test on the final vector space
c---
      errm=0.0d0
      do j=1,n
        do i=1,j
          aaa=0.0d0
          do k=1,n
          aaa=aaa+U(k,i)*U(k,j)
          enddo
          if(i.eq.j) aaa=1.d0-aaa
          if(dabs(aaa).gt.errm) then
            iii=i
            jjj=j
            errm=dabs(aaa)
          endif
        enddo
      enddo
      if(errm.gt.1.d-10) goto 130
c
      return
c
  110 write(6,112) i,anr,tol
  112 format(//' ERROR in Graham: last vector wrong',i4,2d12.3//)
      stop
c
  120 write(6,122) i,anr,tol
  122 format(//' ERROR in Graham: too small norma for row',
     &   i4,d12.3,d13.3//)
      stop
c
  130 write(6,132) iii,jjj,errm
  132 format(//' ERROR in Graham: inaccurate vectors  ',2i6,d12.3//)
      stop
      end
c
c =================================================================
c
      Subroutine Graham2 (U,n,nc,ndim)
c
c  A set of vectors U(*,i) is orthonormalized by Graham-Schmidt method
c  for metric = identity.
c  The number of vectors may be larger than the space dimension and the
c  redundancies are removed along the process. The output matrix is
c  square othonormal.
c  On input : n  is the number of vectors
c             nc is the number of Cartesian coords (=space dim)
c  On output: n  is the number of non-redundant vectors (=space dim)
c  NO special treatment for the last vector is applied
c
      implicit real*8(a-h,o-z)
      dimension U(ndim,*)
      tol=1.d-08
c
      i=1
      Call Normaliz_C (U,n,i,anr,ndim)
      if(anr.lt.tol) goto 120
c
c     Track the non redundant vector with index ii
      ii  = 1
      nnr = n
      do i=2,n
        ii = ii + 1
        Call Normaliz_C (U,n,ii,anr,ndim)
        if(anr.lt.tol) goto 120
        do j=1,ii-1
          aij=0.d0
          do k=1,nc
            aij=aij+U(k,ii)*U(k,j)
          enddo
          do k=1,nc
            U(k,ii)=U(k,ii)-aij*U(k,j)
          enddo
        enddo
        Call Normaliz_C (U,n,ii,anr,ndim)
        if(anr.lt.tol) then
          !Shift vectors in U
          do j=ii,nnr-1
            do k=1,nc
              U(k,j) = U(k,j+1)
            enddo
          enddo
          ii  = ii - 1
          nnr = nnr- 1
        endif
      enddo
      if (ii.ne.nc) goto 140
      n = ii
c---
c---   orthonormalization test on the final vector space
c---
      errm=0.0d0
      do j=1,n
        do i=1,j
          aaa=0.0d0
          do k=1,n
          aaa=aaa+U(k,i)*U(k,j)
          enddo
          if(i.eq.j) aaa=1.d0-aaa
          if(dabs(aaa).gt.errm) then
            iii=i
            jjj=j
            errm=dabs(aaa)
          endif
        enddo
      enddo
      if(errm.gt.1.d-10) goto 130
c
      return
c
  110 write(6,112) i,anr,tol
  112 format(//' ERROR in Graham2: last vector wrong',i4,2d12.3//)
      stop
c
  120 write(6,122) i,anr,tol
  122 format(//' ERROR in Graham2: too small norma for row',
     &   i4,d12.3,d13.3//)
      stop
c
  130 write(6,132) iii,jjj,errm
  132 format(//' ERROR in Graham2: inaccurate vectors  ',2i6,d12.3//)
      stop
c
  140 write(6,142) ii, nc, n
  142 format(//' ERROR in Graham2: wrong dimension  ',3i6//)
      stop
      end
c
c =================================================================
c
      Subroutine Normaliz_C (U,n,i,aaa,ndim)
c
c  Normalize the i-th column of U
c
      implicit real*8(a-h,o-z)
      dimension U(ndim,*)
      tol=1.d-09
c
      aaa=0.0d0
      do k=1,n
        aaa=aaa+U(k,i)**2
      enddo
      aaa=dsqrt(aaa)
      do k=1,n
        U(k,i)=U(k,i)/aaa
      enddo
      return
      end
