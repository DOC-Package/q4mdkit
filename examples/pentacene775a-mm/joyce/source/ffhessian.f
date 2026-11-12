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
      Subroutine FF_Hessian (par,iflag,parfun,parfun2,Hn,Ha,npar
     &                      ,geom,nat,iflagr,nric,labpar
     &                      ,maxpt,maxat,max3n,maxpar)
c----------------------------------------------------------------------c


c         compute the Hessian of the FF at the equilibruim geometry
 
c   a)  numerical mode (FF is computed on a two-dim grid)
 
c   b)  analytical mode (FF derivatives are computed by comp_func)
c       WARNING: If the gradients do not vanish, the formulae are to be
c                complemented with derivatives of the B matrix 
 
c       WARNING: distances must be in bohr !!!    
 
 
c  par     the obtained optimized parameters of the FF
c  iflag   driving array for the definition of the internal coordinate
c          to which each function is referred.
c          it is assumed that the coupling terms come after the diag. terms
c  parfun  fixed parameters values (R0,Ang0...) for model potential functions
c  parfun2 fixed parameters values (R0,Ang0...) for coupled potential functions
c  Hn      cartesian Hessian matrix: (numerical mode)
c  Ha      cartesian Hessian matrix: (analytical mode)
c  npar    n. of parameters (linear) to be optimized
c  geom    current atomic coordinates (in bohr) with nucl charge(4) and mass(5)
c  nat     n. of atoms
c  iflagr  iflag array for the internal coordinates
c  labpar  label of each parameter to be optimized
c  maxpt   max n. of points (ener+grad+hess) (in the caller)
c  maxat   max n. of atoms (in the caller)
c  max3n   maxat*3 (in the caller)
c  maxpar  max n. of FF functions to be optimized (in the caller)


      implicit real*8(a-h,o-z)
      parameter (maxf=2000)
      dimension par(*),iflag(maxpar,0:8),parfun(*),parfun2(*)
      dimension Hn(max3n,*)
     &         ,Ha(max3n,*),geom(maxat,*),iflagr(maxpar,0:8)
      dimension iflagx(maxf,0:8),U(maxf,maxf),indric(maxf)
     &         ,fun0(maxf),fun1(maxf,2),fun2(maxf,3),valic(maxf)
      character labpar(*)*(*),card*80
      common /printout/ levprt

c     copy the iflag array (for convenience) and other

      do i=1,npar
        do k=0,8
          iflagx(i,k)=iflag(i,k)
        enddo
      enddo
      ncart=nat*3
      need=Max0(ncart,npar)
      Call Check_Dim (need,maxf,'FF_Hessian: increase maxf')
      thrdo=1.0d-05
      ncol=10

c     compute the CC Hessian by direct numerical derivatives

      if (levprt .gt. 0 ) write(6,6)
    6 format(
     & ///12x,'================================================='
     &  ,/12x,'  COMPUTE FF-HESSIAN AT geom=0: numerical mode'
     &  ,/12x,'=================================================')
 
      do ia=1,nat
        do ja=1,nat
          do ix=1,3
            do jx=1,3
              aij=0.0d0
              do ip=1,npar
                Call D2f_dx2 (ip,ia,ix,ja,jx,dd2,iflagx,geom
     &                       ,parfun,parfun2,nat,npar,maxf,maxat)
                aij=aij+dd2*par(ip)
              enddo
              ii=(ia-1)*3+ix
              jj=(ja-1)*3+jx
              Hn(ii,jj)=aij
            enddo
          enddo
        enddo
      enddo
      card='FF Hessian (numer) at equil geom in cartesian (bohr,rad)'
      n=ncart
      i=n
      if(levprt.ge.3) Call Pr_Rec2 (Hn,n,i,card,ncol,max3n)

c     compute the CC Hessian by analytical derivatives

      if (levprt .gt. 0 ) write(6,12)
   12 format(
     &  //10x,'================================================='
     &  ,/10x,'  COMPUTE FF-HESSIAN at geom=0: analytical mode'
     &  ,/10x,'=================================================')
      if(levprt.ge.2) then
        write(6,*) '   '
        write(6,14) 'IFLAG ',(i,(iflag (i,k),k=0,8),i=1,npar)
        write(6,14) 'IFLAGR',(i,(iflagr(i,k),k=0,8),i=1,nric)
      endif
   14 format(/9x,'array ',a,/(2i5,2x,8i5))

c     compute the value, gradient and 2nd deriv of the FF functions
c     with respect to the internal coordinates (RICOs)

      Call Comp_Func (iflagx,fun0,fun1,fun2,parfun,parfun2,npar
     &               ,valic,nat,geom,maxat,maxf)
      do ip=1,npar
        fun0(ip)=fun0(ip)*par(ip)
        fun1(ip,1)=fun1(ip,1)*par(ip)
        fun1(ip,2)=fun1(ip,2)*par(ip)
        fun2(ip,1)=fun2(ip,1)*par(ip)
        fun2(ip,2)=fun2(ip,2)*par(ip)
        fun2(ip,3)=fun2(ip,3)*par(ip)
        ia=iflag(ip,1)
        ib=iflag(ip,2)
        ic=iflag(ip,3)
        id=iflag(ip,4)
        ichi=0
        if(iflag(ip,0).ne.8) then
          Call Gia_ic (nric,iflagr,ia,ib,ic,id,ichi,maxpar)
          if(ichi.le.0) then
            write(6,'(//" ERROR in FF_Hessian: FF has a wrong IC")')
            write(6,'(10x,"fun",i5,5x,4i5,//)') ip,ia,ib,ic,id
            stop
          endif
        endif
        indric(ip)=ichi
      enddo
c
      if(levprt.ge.1) then
        write(6,20)
        icoup=0
        do i=1,npar
          it=iflagx(i,0)
          if(it.lt.8) then
            write(6,22) i,labpar(i),fun0(i),fun1(i,1),fun2(i,1)
     &                 ,indric(i)
          else
            icoup=icoup+1
            if(icoup.eq.1) write(6,21)
            write(6,24) i,labpar(i),fun0(i),fun1(i,1),fun1(i,2)
     &                           ,fun2(i,1),fun2(i,2),fun2(i,3)
          endif
        enddo
      endif
   20 format(/8x,'Value of FF functions and deriv. with coeff.s'
     &      ,/32x,'fun     df/dIC  d2f/dIC2')
   21 format( 32x,'fun      df/da     df/db   d2f/daa'
     &                                  ,'    d2fdab   d2f/dbb')
   22 format(i4,2x,a,3f10.5,i5)
   24 format(i4,2x,a,6f10.5)

c     check the first derivatives are null

      thrder=1.d-07
      nout=0
      do i=1,npar
        it=iflagx(i,0)
        if(it.lt.8) then
          if(dabs(fun1(i,1)).gt.thrder) nout=nout+1
        else
          if(dabs(fun1(i,1)).gt.thrder) nout=nout+1
          if(dabs(fun1(i,2)).gt.thrder) nout=nout+1
        endif
      enddo
c     if(nout.gt.0) write(6,32)
      if(nout.gt.0) return
   32 format(//12x,'WARNING: Hess analyt calcul. may be inaccurate'
     &       ,/12x,'since first derivatives do not vanish'
     &       ,/12x,'the calcul is skipped',/)

c     compute the transformation matrix (RIC = U * CC)

      Call Trasf2Int (geom,nat,U,iflagr,nric,maxat,maxf,maxpar)
      ncart=nat*3
      card='Transformation B Matrix: RIC = U * CC (bohr,rad)'
      if(levprt.ge.2) Call Pr_Rec2 (U,nric,ncart,card,ncol,maxf)

c     transform the Hessian from IC to CC

      do i=1,ncart
        do j=1,ncart
          Ha(i,j)=0.0d0
        enddo
      enddo
      do ip=1,npar
        if(iflag(ip,0).ne.8) then
          jp=indric(ip)
          do j=1,ncart
            do i=1,j
              Ha(i,j)=Ha(i,j)+U(jp,i)*U(jp,j)*fun2(ip,1)
            enddo
          enddo
        else
          ia=iflag(ip,1)
          ib=iflag(ip,2)
          ja=indric(ia)
          jb=indric(ib)
          do j=1,ncart
            do i=1,j
              Ha(i,j)=Ha(i,j)+U(ja,i)*U(ja,j)*fun2(ip,1)
              Ha(i,j)=Ha(i,j)+U(ja,i)*U(jb,j)*fun2(ip,2)
              Ha(i,j)=Ha(i,j)+U(jb,i)*U(ja,j)*fun2(ip,2)
              Ha(i,j)=Ha(i,j)+U(jb,i)*U(jb,j)*fun2(ip,3)
            enddo
          enddo
        endif
      enddo
      do j=2,ncart
        do i=1,j-1
          Ha(j,i)=Ha(i,j)
        enddo
      enddo
      card='FF Hessian (analyt) at equil geom in cartesian (bohr,rad)'
      if(levprt.ge.3) Call Pr_Rec2 (Ha,n,i,card,ncol,max3n)

c     check the two Hessians computed in different way

      difmax=0.0d0
      elemax=0.0d0
      elemed=0.0d0
      ia=0
      ja=0
      do j=1,ncart
        do i=1,j
        diff=dabs(Hn(i,j)-Ha(i,j))
          if(diff.gt.difmax) then
            difmax=diff
            ia=i
            ja=j
          endif
          if(dabs(Hn(i,j)).gt.elemax) then
            elemax=dabs(Hn(i,j))
            ie=i
            je=j
          endif
          elemed=elemed+dabs(Hn(i,j))
          if(i.ne.j) elemed=elemed+dabs(Hn(i,j))
        enddo
      enddo
      elemed=elemed/dble(ncart*(ncart+1)/2)
      thrhh=1.d-05*elemed
      write(6,44) ia,ja,difmax,ie,je,elemax,elemed,thrhh
      if(difmax.gt.thrhh) write(6,46)
   44 format(/8x,'Max Diff: H(numer) , H(analyt) :'
     &          ,2i4,1pd14.3
     &      ,/8x,'Max Hessian element    ',8x,':',2i4,d14.3
     &      ,/8x,'Average |Hessian value|',8x,':',8x,d14.3
     &      ,/8x,'Max allowed difference ',8x,':',8x,d14.3)
   46 format(//18x,'WARNING: discrepancies in H are excessive',//)
 
 
      return
      end
