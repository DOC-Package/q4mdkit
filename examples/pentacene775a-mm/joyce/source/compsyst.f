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
      Subroutine Comp_Syst (iflag,parfun,parfun2,escf,G,H,npar,string
     &                     ,alpha,beta,weig,we,wg,whd,who,labpt,npt
     &                     ,wfreq,igeom,valic0,escf0,zgrah
     &                     ,geom,nat,labpar
     &                     ,iflagr,labric,nric
     &                     ,maxpt,maxat,max3n,maxpar)
c----------------------------------------------------------------------c
 
 
c----------------------------------------------------------------------c
c   computes the part of the alpha and beta matrices pertinent to a
c   given geometry
c   WARNING: distances must be in bohr !!!    
 
c   The second derivatives (d2F / dQa dQb) 
c   i.e. derivative of the functions
c   with respect to the normal modes,
c   are computed by a numerical method.   
 
 
c  iflag   driving array for the definition of the internal coordinate
c          to which each function is referred.
c          it is assumed that the coupling terms 
c          come after the diag. terms

c  parfun  equilibrium values for the chosen IC's 
c          (default: set from QM)
c  parfun2 equilibrium values for the second IC involved in a coupling (if any)
c          (default: set from QM)

c  escf    energy of the current conformation
c  G       cartesian Gradients (3*nat) of the current conformation
c  H       cartesian Hessian matrix for the current conformation
 
c  npar    n. of parameters (linear) to be optimized

c  string  a label for the current conformation

c  alpha   alpha(mu,a) = a-th function in the mu point. 
c          A mu point may be energy, or first deriv. 
c          with respect to a normal mode, or
c          a second deriv. with respect to a pair of normal modes.

c  beta    exact value of the energy, gradient or 
c          hessian in the mu-point

c  weig    weigths of each point

c  we,wg,whd,who   weights for: energy, grad, 
c                  diag hessian, off-diag hessian, respectively

c  labpt   a label for each point

c  npt     n. of points

c  wfreq   constant to rescale the weigths of the frequencies, 
c          in order to emphasize the lower frequencies

c  igeom   index of the current geometry: 
c          0= equil 1,2,3 ...= other geometries

c  valic0  value of the internal coordinates 
c          at the equilibrium geometry

c  escf0   value of the energy at the equilibrium geometry

c  zgrah    a parameter for the Comp_Vibr routine

c  geom    current atomic coordinates (in bohr) 
c          with nucl charge(4) and mass(5))

c  nat     n. of atoms

c  labpar  label of each parameter to be optimized

c  iflagr  iflag array for the internal coordinates

c  labric  label for each redundant internal coordinate

c  nric    n. of internal redundant coordinates

c  maxpt   max n. of points (ener+grad+hess) (in the caller)
c  maxat   max n. of atoms (in the caller)
c  max3n   maxat*3 (in the caller)
c  maxpar  max n. of parameters to be optimized (in the caller)
c----------------------------------------------------------------------c
 
 
      implicit real*8(a-h,o-z)
      parameter (maxf=2000, maxtri=maxf*(maxf+1)/2)
      dimension iflag(maxpar,0:8),parfun(*),parfun2(*),G(*),H(max3n,*)
     &         ,alpha(maxpt,*),beta(*),weig(*)
     &         ,valic0(*),geom(maxat,*),iflagr(maxpar,0:8)
      dimension iflagx(maxf,0:8),wmass(maxf)
     &         ,fun0(maxf),fun1(maxf,2),fun2(maxf,3),valic(maxf)
      dimension U(maxf,maxf),C(maxf,maxf),T(maxf,maxf)
     &         ,E(maxf),F(maxtri),Freq(maxf),Gnm(maxf),Cm(maxf,maxf)
     &         ,Z(maxf,maxf),ncc(maxf),indcc(maxf,maxf),mflag(0:8)
      character string*(*),labpar(*)*(*),labpt(*)*(*),labric(*)*(*)
     &         ,scelta*10,card*80
      logical zene,zgrd,zhes,zpic,zgrah,z_LJassign
      common /printout/ levprt
      common /ljassign/z_LJassign
      data amu/1822.8880d0/  !!   1 amu = 1822.888 electron mass


      if(levprt.ge.1) write(6,7) string(1:Lenw(string))
      if(levprt.lt.1) write(6,8) string(1:Lenw(string))
      write(6,9) we,wg,whd,who,wfreq
    7 format(/6x,'=============================================='
     &      ,/19x,a
     &      ,/16x,'     increment the Alpha , Beta matrices'
     &      ,/16x,'==============================================',/)
    8 format(/10x,'INCREMENT Alpha, Beta matrices for ',a)
    9 format( 10x,'Energy weight ......................',f10.4
     #      ,/10x,'Gradient weights ...................',f10.4
     #      ,/10x,'diag Hessian weights ...............',f10.4
     #      ,/10x,'off-diag Hessian weights ...........',f10.4
     #      ,/10x,'further freq dep. weight for Hessian',f10.4)
      if(levprt.ge.2) then
        write(6,4) 'IFLAG ',(i,(iflag (i,k),k=0,8),i=1,npar)
        write(6,4) 'IFLAGR',(i,(iflagr(i,k),k=0,8),i=1,nric)
      endif
    4 format(/9x,'array ',a,/(2i5,2x,8i5))
 
      ncol=10
      ncart=nat*3
      need=Max0(ncart,npar)
      Call Check_Dim (need,maxf,'COMP_SYST: increase maxf')
      thrdo=1.0d-05
      zene=(we.gt.thrdo)
      zgrd=(wg.gt.thrdo)
      zhes=(whd.gt.thrdo)
 
c     computes array: sqrt(mass)
 
      m=0
      do i=1,nat
        aa=geom(i,5)
        wmass(m+1)=dsqrt(aa*amu)
        wmass(m+2)=wmass(m+1)
        wmass(m+3)=wmass(m+1)
        m=m+3
      enddo
 
c     copy the iflag array (for convenience)
c     and compute:   npar0 = n. of param.s without coupling terms
 
      npar0=0
      do i=1,npar
        do k=0,8
          iflagx(i,k)=iflag(i,k)
        enddo
        if(iflagx(i,0).lt.8) npar0=i
      enddo
 
c     compute the value, gradient and 2nd deriv of the FF functions
c     with respect to the internal coordinates (RICOs)


      Call Comp_Func (iflagx,fun0,fun1,fun2,parfun,parfun2,npar
     &               ,valic,nat,geom,maxat,maxf)
      if(levprt.gt.0) then
        write(6,20)
        icoup=0
        do i=1,npar
          it=iflagx(i,0)
          if(it.lt.8) then
            write(6,22) i,labpar(i),fun0(i),fun1(i,1),fun2(i,1)
     &                 ,valic(i),valic0(i),valic(i)-valic0(i)
          else
            icoup=icoup+1
            if(icoup.eq.1) write(6,23)
            ia=iflagx(i,1)
            ib=iflagx(i,2)
            write(6,24) i,labpar(i),fun0(i),ia,ib,fun1(i,1),fun1(i,2)
     &                 ,fun2(i,1),fun2(i,2),fun2(i,3)
          endif
        enddo
      endif
   20 format(/8x,'Value of functions and deriv. at the current geometry'
     &      ,/'   #',28x,'fun       df/dIC   d2f/dIC2      '
     &                  ,'RIC       RIC(0)    change'
     &      ,/'    ',28x,'Hart      Hart/b    Hart/bb      '
     &                  ,'b/rad     b/rad)')
   22 format(i4,2x,a20,6f11.5)
   23 format( '  #',28x,'fun    a   b    df/da    df/db   d2f/daa',
     &        '  d2f/dab  d2f/dbb')
   24 format(i4,2x,a20,f9.3,2i4,5f9.3)
 
c     prepare unit = 48 (freq, normal modes, for successive use)
 
      if(igeom.eq.0) then
        open (unit=48,file='provvi.48',form='unformatted')
        rewind 48
        kn=-1
        write(48) kn
      endif
 
 
c--------------------------------------------------------------
c       P R E L I M I N A R Y     C A L C U L A T I O N S
c--------------------------------------------------------------
 
 
      if(.not.zgrd .and. (.not.zhes)) goto 100
 
c     computes the transformation matrix (RIC = U * CC)
 
      Call Trasf2Int (geom,nat,U,iflagr,nric,maxat,maxf,maxpar)
      ncart=nat*3
      card='Transformation U Matrix: RIC = U * CC (bohr,rad)'
      if(levprt.ge.2) Call Pr_Rec2 (U,nric,ncart,card,ncol,maxf)
 
c     compute the normal modes
 
      do j=1,ncart
        do i=1,ncart
          T(i,j)=H(i,j)
        enddo
      enddo
      write(scelta,'("GEOMETRY",i2)') igeom
      zpic=.false.
      if ( dabs(who+whd) > 1.d-12 ) then
      Call Comp_Vibr (T,U,F,E,Freq,C,geom,scelta,zpic,zgrah
     &               ,nric,nat,nvib,labric,maxf,maxat)
      else
      Call Comp_Eckart (Freq,C,geom,scelta
     &                 ,nat,nvib,maxf,maxat)
      endif
 
c     write on unit 48 for successive use
 
      if(igeom.eq.0) then
        rewind 48
        write(48) nvib,ncart,(E(i),i=1,nvib)
        write(48) ((C(i,nm),i=1,ncart),nm=1,nvib)
      endif
 
c  compute the Cm matrix: Cm(i,K)=C(i,K)/sqrt(mass.i)
 
      do K=1,nvib
        do i=1,ncart
          Cm(i,K)=C(i,K)/wmass(i)
        enddo
      enddo
 
c     compute the T matrix: T(mu,K)=dq(mu)/dQK  QK=normal mode
c     WARNING: in this case some RICO's may be repeated (for torsions)
 
      maxfl=maxf
      Call Trasf2Int (geom,nat,U,iflagx,npar0,maxat,maxf,maxfl)
      do nm=1,nvib
        do int=1,npar0
          aa=0.d0
          do k=1,ncart
            aa=aa+U(int,k)*C(k,nm)/wmass(k)
          enddo
          T(int,nm)=aa
        enddo
      enddo
      card='Matrix T = U*M(-1/2)*C = dRICO/dNM (RICO may be repeated)'
      if(levprt.ge.3) Call Pr_Rec2 (T,npar0,nvib,card,ncol,maxf)
  100 continue
      nptbeg=npt
 
 
c     E N E R G Y
 
 
      if(.not.zene) goto 200
      npt=npt+1
      Call Check_Dim (npt,maxpt,'COMP_SYST: increase maxpt')
      do ip=1,npar
        alpha(npt,ip)=fun0(ip)
      enddo
      beta(npt)=(escf-escf0)              !!! in Hartree
      weig(npt)=we
      write(labpt(npt),'("Energy-",i2)') igeom
      write(6,'(10x,"Computed Alpha,Beta for ENERGY  point",i5)') npt
      if(levprt.ge.3) write(6,120) igeom,(i,alpha(npt,i),i=1,npar)
  120 format(/5x,'value of the functions at geom',i4
     #      ,/(5(i5,f9.4)))
 
 
c     G R A D I E N T
 
 
  200 if(.not.zgrd) goto 300
 
c     g09 gradients in the normal mode basis
 
      angrad=nvib
      do k=1,nvib
        gg=0.0d0
        do i=1,ncart
          gg=gg+G(i)*C(i,k)/wmass(i)
        enddo
        Gnm(k)=gg
      enddo
      write(6,'(10x,"Computed Alpha,Beta for GRADIENT points",2i5)')
     &      npt+1,npt+nvib-1
      write(6,212) (k,Freq(k),1000.d0*Gnm(k),k=1,nvib)
  212 format(/8x,'En Gradients with resp. to NM (mH/bohr)'
     &      ,/(i5,f9.1,f11.6))
 
c     gradient of each function with resp. to NM
 
      do K=1,nvib
        npt=npt+1
        Call Check_Dim (npt,maxpt,'COMP_SYST: increase maxpt')
        beta(npt)=Gnm(K)
        weig(npt)=wg/angrad
        write(labpt(npt),'("grd",i2,i4)') igeom,K
        do ip=1,npar0
          alpha(npt,ip)=fun1(ip,1)*T(ip,K)
        enddo
        if(npar.gt.npar0) then
          do ip=npar0+1,npar
            ia=iflag(ip,1)
            ib=iflag(ip,2)
            alpha(npt,ip)=fun1(ip,1)*T(ia,K)+fun1(ip,2)*T(ib,K)
          enddo
        endif
      enddo
 
 
c     H E S S I A N
 
 
  300 if(.not.zhes) goto 400
 
c     computes the scaling factor for the frequencies,
c     if requested through the wfreq number.
c     Notice that the sum of the hessian weigths does not change 
 
      anhess=dble(nvib*(nvib+1)/2)
      if(wfreq.gt.0.0d0) then
        worig=dble(nvib)*whd+dble(nvib*(nvib-3)/2)*who
        wtot=0.0d0
        do K=1,nvib
          do L=1,K
            www=whd
            if(K.ne.L) www=who
            Wkl=www/dsqrt((wfreq+Freq(K))*(wfreq+Freq(L)))
            wtot=wtot+Wkl
          enddo
        enddo
        wplus=worig/wtot
      endif
 
c     computes second derivatives with respect to NM numerically

      npt0=npt
      do 320 K=1,nvib
        do 310 L=1,K   !! the K,L pair identifies a point
          www=whd
          if(K.ne.L) www=who
          if(www.lt.1.d-05) goto 310
          if(wfreq.gt.0.0d0) then
            www=wplus*www/dsqrt((wfreq+Freq(K))*(wfreq+Freq(L)))
          endif
          npt=npt+1
          Call Check_Dim (npt,maxpt,'COMP_SYST: increase maxpt')
          weig(npt)=www/anhess
          beta(npt)=0.0d0
          if(K.eq.L) beta(npt)=E(K)
          write(labpt(npt),'("hes",i2,2i3)') igeom,K,L
          do ip=1,npar  !! ip identifies a function
            Call D2f_dQ2 (ip,K,L,ader2,iflagx,geom,Cm,parfun,parfun2,
     &                    nat,npar,maxf,maxat)
            alpha(npt,ip)=ader2
          enddo
  310   continue
  320 continue
      write(6,'(10x,"Computed Alpha,Beta for HESSIAN points",2i5)')
     &      npt0+1,npt
 
 
  400 continue
      nnew=npt-nptbeg
      write(6,410) nnew,string(1:Lenw(string))
  410 format(10x,'Computed',i5,'  new points for geom: ',a)
      return
      end



c----------------------------------------------------------------------c
      Subroutine Get_Uder (der,mflag,ii,jj,nat,geom,maxat)
c----------------------------------------------------------------------c
 
c  computes the derivative of the transformation matrix U,
c  with respect to a cartesian coordinate
c     d2 q(mu) / dxi dxj   =   d/dxi U(mu,j)   =   d/dxj U(mu,i)
c   where the details of the mu-th RICO are specified in /mflag/ 
 
c  Angle derivatives are in radians/bohr  (case: q=angle or dihedral)
c  Coordinate derivatives are c-numbers   (case: q=distance)
 
c  WARNING: the atomic coordinates in /geom/ are assumed to be in Angstrom
 
      implicit real*8(a-h,o-z)
      dimension geom(maxat,*),mflag(0:8)
 
      dr=2.d-04
      dnorm=12.0d0*dr
      iat=(jj-1)/3+1
      ico=jj-3*(iat-1)
      gsave=geom(iat,ico)
 
      geom(iat,ico)=gsave-dr-dr
      Call Get_1_U (geom,nat,mflag,ii,v1,dr,maxat)
      der=0.0d0
      if(dabs(v1).lt.1.d-12) return
      geom(iat,ico)=gsave-dr
      Call Get_1_U (geom,nat,mflag,ii,v2,dr,maxat)
      geom(iat,ico)=gsave+dr
      Call Get_1_U (geom,nat,mflag,ii,v4,dr,maxat)
      geom(iat,ico)=gsave+dr+dr
      Call Get_1_U (geom,nat,mflag,ii,v5,dr,maxat)
      geom(iat,ico)=gsave
      der=(v1-8.d0*v2+8.d0*v4-v5)/dnorm
      return
      end



c----------------------------------------------------------------------c
      Subroutine Get_1_U (geom,nat,mflag,ii,Uqx,dr,maxat)
c----------------------------------------------------------------------c
 
c  compute one element of the transformation matrix dq/dxi
c       q   = RICO determined by the mflag array
c       x   = CC identified by the ii index
 
c  Angle derivatives are in radians/bohr  (case: q=angle or dihedral)
c  Coordinate derivatives are c-numbers   (case: q=distance)
 
c  WARNING: the atomic coordinates in /geom/ are assumed to be in Angstrom
 
      implicit real*8(a-h,o-z)
      parameter (max=99)
      dimension geom(maxat,*),mflag(0:8)
     &         ,faux(max,5),gaux(max,5),ddr(4),valu(4)
      data rad2deg /57.2957795130823229d0/
 
      dnorm=12.0d0*dr
      ddr(1)=-dr-dr
      ddr(2)=-dr
      ddr(3)=+dr
      ddr(4)=+dr+dr
      do i=1,nat
        do k=1,4
          gaux(i,k)=geom(i,k)
          faux(i,k)=geom(i,k)
        enddo
      enddo
 
c     identifies atom and x,y,z component of the CartCoord=ii
 
      iato=(ii-1)/3+1
      ixyz=ii-(iato-1)*3
 
      itype=mflag(0)
      if(itype.gt.3) itype=3   !! 1=bond  2=angle  3=dihedral  4=dihedral
      Uqx=0.0d0
 
c      DISTANCES 
 
      if(itype.ne.1) goto 20
      ia=mflag(1)
      ib=mflag(2)
      if(iato.ne.ia .and. iato.ne.ib) return
      do ispo=1,4
        gaux(iato,ixyz)=faux(iato,ixyz)+ddr(ispo)
        val=dsqrt((gaux(ia,1)-gaux(ib,1))**2
     &           +(gaux(ia,2)-gaux(ib,2))**2
     &           +(gaux(ia,3)-gaux(ib,3))**2)
        valu(ispo)=val
      enddo
      der=(valu(1)-8.d0*valu(2)+8.d0*valu(3)-valu(4))/dnorm
      Uqx=der
      return
 
c      ANGLES AND DIHEDRALS   
 
   20 continue
      ia=mflag(1)
      ib=mflag(2)
      ic=mflag(3)
      id=mflag(4)
      nt=Min0(4,itype+1)
      if(iato.eq.ia) goto 30
      if(iato.eq.ib) goto 30
      if(iato.eq.ic) goto 30
      if(nt.eq.4 .and. iato.eq.id) goto 30
      return
   30 do ispo=1,4
        gaux(iato,ixyz)=faux(iato,ixyz)+ddr(ispo)
        Call Value_Ang (gaux,max,ia,ib,ic,id,val,'d')
        valu(ispo)=val
        if(ispo.gt.1) then   !!! test for +/-180 case
          if(valu(ispo)-valu(ispo-1).gt.20.d0) then
            valu(ispo)=valu(ispo)-360.d0
          elseif(valu(ispo-1)-valu(ispo).gt.20.d0) then
            valu(ispo)=valu(ispo)+360.d0
          endif
        endif
      enddo
      der=(valu(1)-8.d0*valu(2)+8.d0*valu(3)-valu(4))/dnorm
      Uqx=der/rad2deg
 
      return
      end


c----------------------------------------------------------------------c
      Subroutine D2f_dQ2 (fu,q1,q2,ader2,iflag,geom,Cm,parfun,parfun2,
     &                    nat,npar,maxf,maxat)
c----------------------------------------------------------------------c
 
c  computes second derivates w.r.o normal coordinates q1 and q2
 
      implicit real*8(a-h,o-z)
      parameter (max=2000)
      dimension iflag(maxf,0:8),parfun(*),parfun2(*),ifx(max,0:8)
     &         ,geom(maxat,*),geo(max,3),Cm(maxf,*)
     &         ,fun0(max),fun1(max,2),fun2(max,3)
     &         ,valu(-2:2),val(-2:2,-2:2),ppp(max)
      integer fu,q1,q2
      data icall/0/
      save icall
      logical zcheck
 
      zcheck=.false.
      icall=icall+1
      if(icall.eq.1 .and. zcheck) write(6,4)
    4 format(20(/5x,"WARNING: zcheck in D2f_dQ2: activated"))
 
      dq=0.0010d0   !!  increment for finite differences
      dfac=1.d0/(12.d0*dq)
      nf=nat*3
      Call Check_Dim (nf,max,'D2F_DQ2: increase max')
      mmx=max
 
c     copy iflag (for array problem)
 
      do i=1,npar
        do k=0,8
          ifx(i,k)=iflag(i,k)
        enddo
      enddo
 
c      computes:   d2.fu / d.Q(q1) d.Q(q2)
 
      do i=-2,2
        dqi=dq*dble(i)
        do j=-2,2
          dqj=dq*dble(j)
          ik=0
          do k=1,nat
            do l=1,3
              ik=ik+1
              geo(k,l)=geom(k,l)+Cm(ik,q1)*dqi+Cm(ik,q2)*dqj
            enddo
          enddo
          Call Comp_Fun1 (ifx,funf,parfun,parfun2,fu,nat,geo,max,mmx)
          val(i,j)=funf
          if(zcheck) then
            Call Comp_Func (ifx,fun0,fun1,fun2,parfun,parfun2,npar
     &                     ,ppp,nat,geo,max,mmx)
            valold=fun0(fu)
          endif
        enddo
      enddo
 
      do i=-2,2
        der=(val(i,-2)-8.d0*val(i,-1)+8.d0*val(i,1)-val(i,2))*dfac
        valu(i)=der
      enddo
      der=(valu(-2)-8.d0*valu(-1)+8.d0*valu(1)-valu(2))*dfac
      ader2=der
 
      return
      end


 
c----------------------------------------------------------------------c
      Subroutine D2f_dx2 (fu,ia,ix,ja,jx,d2,iflag,geom
     &                   ,parfun,parfun2,nat,npar,maxf,maxat)
c----------------------------------------------------------------------c
 
c     computes the second derivatives w.r.o.
c     cartesian coordinates (ia,ix) (ja,jx)    
c     ia,ja=atom  ix,jx=1,2,3(x,y,z)
 
      implicit real*8(a-h,o-z)
      parameter (max=1000)
      dimension iflag(maxf,0:8),parfun(*),parfun2(*),ifx(max,0:8)
     &         ,geom(maxat,*),geo(max,3)
     &         ,valu(-2:2),val(-2:2,-2:2)
     &         ,fun0(max),fun1(max,2),fun2(max,3),ppp(max)
      integer fu
      logical zvai,zcheck
      data icall/0/
      save icall
 
      zcheck=.false.
      icall=icall+1
      if(icall.eq.1 .and. zcheck) write(6,4)
    4 format(20(/5x,"WARNING: zcheck in D2f_dx2: activated"))
 
      dq=0.00010d0   !!  increment for finite differences
 
      dfac=1.d0/(12.d0*dq)
      nf=nat*3
      Call Check_Dim (nf,max,'D2F_DX2: increase max')
      mmx=max
 
c     copy iflag (for array problem)
 
      do i=1,npar
        do k=0,8
          ifx(i,k)=iflag(i,k)
        enddo
      enddo
 
c     computes the involved atoms in fun(fu)
 
      if(iflag(fu,0).ne.8) then
        zvai=.false.
        do k=1,4
          ka=ifx(fu,k)
          if(ka.eq.ia) zvai=.true.
          if(ka.eq.ja) zvai=.true.
        enddo
        d2=0.0d0
        if(.not.zvai) return
      endif
 
c     computes   d2.fu / d.Cart1 d.Cart2
 
      do i=-2,2
        dqi=dq*dble(i)
        do j=-2,2
          dqj=dq*dble(j)
          do k=1,nat
            do l=1,3
              geo(k,l)=geom(k,l)
            enddo
          enddo
          geo(ia,ix)=geo(ia,ix)+dqi
          geo(ja,jx)=geo(ja,jx)+dqj
          Call Comp_Fun1 (ifx,funf,parfun,parfun2,fu,nat,geo,max,mmx)
          val(i,j)=funf

          if(zcheck) then
            Call Comp_Func (ifx,fun0,fun1,fun2,parfun,parfun2,npar
     &                     ,ppp,nat,geo,max,mmx)
            fung=fun0(fu)
            diff=fung-funf
            if(dabs(diff).gt.1.d-10) then
              write(6,16) fu,(ifx(fu,k),k=0,8),fung,funf,diff
              stop
            endif
          endif
   16 format(/5x,'ERROR in checking d2f_dx2 routine:  func',i4
     &      ,/5x,'iflag vector ......',i4,3x,8i4
     &      ,/5x,'COMP_FUNC (general routine) .....',f16.10
     &      ,/5x,'COMP_FUN1 (specific routine) ....',f16.10
     &      ,/5x,'difference ......................',f16.10,/)

        enddo
      enddo
 
      do i=-2,2
        der=(val(i,-2)-8.d0*val(i,-1)+8.d0*val(i,1)-val(i,2))*dfac
        valu(i)=der
      enddo
      der=(valu(-2)-8.d0*valu(-1)+8.d0*valu(1)-valu(2))*dfac
      d2=der

      return
      end
