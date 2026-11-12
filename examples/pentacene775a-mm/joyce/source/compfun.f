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
      Subroutine Comp_Func (iflag,fun0,fun1,fun2,parfun,parfun2,npar
     &                     ,valic,nat,geom,maxat,maxpar)
c----------------------------------------------------------------------c
 
 
c----------------------------------------------------------------------c
c 1. computes the value of the functions used to describe
c    the Force Field, from the cartesian nuclear coordinates
c 2. computes the value of the involved IC 
 
c fun0(ipar)    ----->   ENERGY   function of the ipar 'th parameter
 
c fun1(ipar,k)  ----->   GRADIENT function of the ipar 'th parameter
c                        with respect to the k-coordinate
c                        k=1   for iflag(ipar,0) < 5 ; 
c                        k=1,2 for iflag(ipar,0) = 8 (coupling) 
 
c fun2(ipar,k)  ----->   2nd DER. function of the ipar 'th parameter
c                        for iflag(ipar,0) < 5   
c                        k=1  : double derivative in the first coord
c                               for iflag(ipar,0) = 8  (coupling)
c                        k=1  : double derivative in the first coord
c                        k=2  : mixed double derivative 
c                        k=3  : double derivative in the second coord
 
c    WARNING: 1- the input geometry must be in bohr
c             2- the parfun data must be in bohr or radians
c----------------------------------------------------------------------c
  
      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension parfun(*),parfun2(*),valic(*)
      dimension fun0(*),fun1(maxpar,2),fun2(maxpar,3)
      logical   z_LJassign

      common /ljassign/z_LJassign
 
      do ip=1,npar

        it=iflag(ip,0)
        ia=iflag(ip,1)
        ib=iflag(ip,2)
        ic=iflag(ip,3)
        id=iflag(ip,4)

        if(it.eq.1) then     !!! stretching
          val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &              +(geom(ia,2)-geom(ib,2))**2
     &              +(geom(ia,3)-geom(ib,3))**2 )
          valic(ip)  = val
          dval       = val-parfun(ip)
          fun0(ip)   = 0.5d0*dval*dval
          fun1(ip,1) = dval
          fun1(ip,2) = 0.d0
          fun2(ip,1) = 1.0d0
          fun2(ip,2) = 0.d0
          fun2(ip,3) = 0.d0

        elseif(it.eq.2) then  !!! bending
          Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
          valic(ip)  = val
          dval       = val-parfun(ip)
          fun0(ip)   = 0.5d0*dval*dval
          fun1(ip,1) = dval
          fun1(ip,2) = 0.d0
          fun2(ip,1) = 1.0d0
          fun2(ip,2) = 0.d0
          fun2(ip,3) = 0.d0

        elseif(it.eq.3) then  !!! improper torsion
          Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
          Call Shifta (val,parfun(ip))
          dval       = val-parfun(ip)
          valic(ip)  = val
          fun0(ip)   = 0.5d0*dval*dval
          fun1(ip,1) = dval
          fun1(ip,2) = 0.d0
          fun2(ip,1) = 1.0d0
          fun2(ip,2) = 0.d0
          fun2(ip,3) = 0.d0
          if(dabs(dval).gt.3.14159265359) then
            a0=parfun(ip)
            write(6,40) ip,ia,ib,ic,id,val,a0,fun0(ip)
            stop
          endif

        elseif(it.eq.4) then  !!! proper torsion
          Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
          valic(ip)=val
          an  = iflag(ip,5)
          arg = an*val-parfun(ip)
          cc  = dcos(arg)
          fun0(ip)  =1.d0+cc
          fun1(ip,1)=-an*dsin(arg)
          fun1(ip,2) = 0.d0
          fun2(ip,1)=-an*an*cc
          fun2(ip,2) = 0.d0
          fun2(ip,3) = 0.d0

        elseif(it.eq.5) then     !!! LJ intra
          if (z_LJassign) then 
            fun0(ip)   = 0.d0
            fun1(ip,1) = 0.d0
            fun1(ip,2) = 0.d0
            fun2(ip,1) = 0.d0
            fun2(ip,2) = 0.d0
            fun2(ip,3) = 0.d0
          else
            val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                +(geom(ia,2)-geom(ib,2))**2
     &                +(geom(ia,3)-geom(ib,3))**2 )
            valic(ip)=val
            sigma=parfun(ip)
            val1       = sigma/val
            valsq      = val1*val1
            val6       = valsq*valsq*valsq
            val12      = val6*val6
            rm1        = 1.d0/val
            rm2        = rm1*rm1
            fun0(ip)   = 4.d0*(val12-val6)
            fun1(ip,1) = 4.d0*(6.d0*val6 - 12.d0*val12)*rm1
            fun1(ip,2) = 0.d0
            fun2(ip,1) = 4.d0*(156.d0*val12 - 42.d0*val6)*rm2
            fun2(ip,2) = 0.d0
            fun2(ip,3) = 0.d0
          endif

        elseif(it.eq.6) then     !!! charge - charge interactions
          if (z_LJassign) then 
            fun0(ip)   = 0.d0
            fun1(ip,1) = 0.d0
            fun1(ip,2) = 0.d0
            fun2(ip,1) = 0.d0
            fun2(ip,2) = 0.d0
            fun2(ip,3) = 0.d0
          else
            val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                +(geom(ia,2)-geom(ib,2))**2
     &                +(geom(ia,3)-geom(ib,3))**2 )
            valic(ip)=val
            rm1        = 1.d0/val
            fun0(ip)   = rm1
            fun1(ip,1) = rm1*rm1
            fun1(ip,2) = 0.d0
            fun2(ip,1) = -2.d0*fun1(ip,1)*rm1
            fun2(ip,2) = 0.d0
            fun2(ip,3) = 0.d0
          endif
 
c       cross terms:  couplings section
 
        elseif(it.eq.8) then

          ifa=iflag(ip,1)
          ifb=iflag(ip,2)
          
c         Get geometrical values
c         Uncoupled coordinate 1 (UC1)
          it_uc1=iflag(ifa,0)
          ia    =iflag(ifa,1)
          ib    =iflag(ifa,2)
          ic    =iflag(ifa,3)
          id    =iflag(ifa,4)
          if (it_uc1 .eq.1) then       !! a = str
            val_uc1 = dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                      +(geom(ia,2)-geom(ib,2))**2
     &                      +(geom(ia,3)-geom(ib,3))**2 )
          elseif (it_uc1.eq.2) then   !! a = bend
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc1,'r')
          else                   !! a = tors(i,p)
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc1,'r')
            Call Shifta (val_uc1,parfun(ip))
          endif
c         Uncoupled coordinate 2 (UC2)
          it_uc2=iflag(ifb,0)
          ia    =iflag(ifb,1)
          ib    =iflag(ifb,2)
          ic    =iflag(ifb,3)
          id    =iflag(ifb,4)
          if (it_uc2.eq.1) then        !! b = str
            val_uc2 = dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                      +(geom(ia,2)-geom(ib,2))**2
     &                      +(geom(ia,3)-geom(ib,3))**2 )
          elseif (it_uc2.eq.2) then   !! b = bend
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc2,'r')
          else                   !! b = tors(i,p)
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc2,'r')
            Call Shifta (val_uc2,parfun2(ip))
          endif

c         Compute the functions and derivatives
          if (it_uc1.le.3)   then         !! a = str, bend, i.tors
           dval_uc1   = val_uc1-parfun(ip)
           if (it_uc2.le.3)  then         !! b = str, bend, i.tors
            dval_uc2   = val_uc2-parfun2(ip)
            fun0(ip)   = dval_uc1*dval_uc2      !! F        = (a-a0)*(b-b0)
            fun1(ip,1) = dval_uc2               !! dF/da    =  b-b0        
            fun1(ip,2) = dval_uc1               !! dF/db    =  a-a0
            fun2(ip,1) = 0.d0                   !! d2F/da2  =  0
            fun2(ip,2) = 1.d0                   !! d2F/dadb =  1
            fun2(ip,3) = 0.d0                   !! d2F/db2  =  0
           elseif(it_uc2.eq.4) then       !! b = p. torsion
            an_uc2  = iflag(ip,5)
            sdph2      = dsin(an_uc2*val_uc2-parfun2(ip))
            cdph2      = dcos(an_uc2*val_uc2-parfun2(ip))
c                      F = (a-a0)*(1+sin(n*b-g))
            fun0(ip)   = dval_uc1*(1.d0+sdph2)
c                  dF/da = (1+sin(n*b-g))
            fun1(ip,1) = (1.d0+sdph2)     
c                  dF/db = (a-a0)*n*cos(n*b-g)
            fun1(ip,2) = dval_uc1*an_uc2*cdph2
c                dF2/da2 = 0
            fun2(ip,1) = 0.d0
c                dF2/dadb  = n*cos(n*b-g)
            fun2(ip,2) = an_uc2*cdph2
c                dF2/db2   = -(a-a0)*n^2*sin(n*b-g) 
            fun2(ip,3) = -dval_uc1*an_uc2**2*sdph2
           else                      !! not defined
            write(6,'(//8x,"ERROR: wrong coupled coordinate !",/)')
            stop
           endif
          elseif(it_uc1.eq.4) then       !! a = torsion
           if(it_uc2.le.3) then          !! b = str. , bend or i.tors
            write(6,'(//8x,
     &      "ERROR: first coupled coordinate should NOT be a torsion!",
     &       /)')
            stop
           elseif(it_uc2.eq.4) then      !! b = torsion
            an_uc1     = dble(iflag(ip,5))
            an_uc2     = dble(iflag(ip,6))
            sdph1      = dsin(an_uc1*val_uc1-parfun(ip))
            sdph2      = dsin(an_uc2*val_uc2-parfun2(ip))
            cdph1      = dcos(an_uc1*val_uc1-parfun(ip))
            cdph2      = dcos(an_uc2*val_uc2-parfun2(ip))
c                      F = sin(na*a-ga)*sin(nb*b-gb)
            fun0(ip)   =   sdph1*sdph2
c                  dF/da = na*cos(n*a-ga)*sin(nb*b-g)  
            fun1(ip,1) =   an_uc1*cdph1*sdph2
c                  dF/db = nb*sin(na*a-ga)*cos(nb*b-gb)
            fun1(ip,2) =   an_uc2*sdph1*cdph2
c                dF2/da2 = -na^2*sin(na*a-ga)*sin(nb*b-gb)
            fun2(ip,1) =  -an_uc1**2*sdph1*sdph2
c                dF2/dadb  = na*nb*cos(na*a-ga)*cos(nb*b-gb)
            fun2(ip,2) =   an_uc1*an_uc2*cdph1*cdph2
c                dF2/db2   = -nb^2*sin(na*a-ga)*sin(nb*b-gb)
            fun2(ip,3) =  -an_uc2**2*sdph1*sdph2
           else                      !! not defined
            write(6,'(//8x,"ERROR: wrong coupled coordinate !",/)')
            stop
           endif

          endif

        else
          write(6,'(//8x,"ERROR: wrong function",/)')
          stop
        endif
      enddo
c
   40 format(//8x,'IRREPARABLE ERROR in subr. COMP_FUNC'
     &       ,/8x,'improper torsion',i5,5x,'atoms',4i4
     &       ,/8x,'angle , angle0',2f9.3,4x,'are differ too much'
     &       ,/8x,'function values is',f12.5,///)
      return
      end


c----------------------------------------------------------------------c
      Subroutine Comp_Fun1 (iflag,fun0,parfun,parfun2,ipar
     &                     ,nat,geom,maxat,maxpar)
c----------------------------------------------------------------------c
 
c     computes the value of the functions used to describe
c     the Force Field, from the cartesian nuclear coordinates
 
      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),geom(maxat,*),parfun(*),parfun2(*)
      character st1*1
      data st1/'r'/
      logical   z_LJassign

      common /ljassign/z_LJassign
 
      ip=ipar
      it=iflag(ip,0)
      ia=iflag(ip,1)
      ib=iflag(ip,2)
      ic=iflag(ip,3)
      id=iflag(ip,4)
      if(it.eq.1) then
        val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &            +(geom(ia,2)-geom(ib,2))**2
     &            +(geom(ia,3)-geom(ib,3))**2 )
        dval       = val-parfun(ip)
        fun0 = 0.5d0*dval*dval
      elseif(it.eq.2) then
        Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
        dval = val-parfun(ip)
        fun0 = 0.5d0*dval*dval
      elseif(it.eq.3) then
        Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
        Call Shifta (val,parfun(ip))
        dval = val-parfun(ip)
        fun0 = 0.5d0*dval*dval
        if(dabs(dval).gt.3.14159265359) then
          a0=parfun(ip)
          write(6,40) ip,ia,ib,ic,id,val,a0,fun0
          stop
        endif
      elseif(it.eq.4) then
        ic=iflag(ip,3)
        id=iflag(ip,4)
        Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
        an  = iflag(ip,5)
        fun0=1.d0+dcos(an*val-parfun(ip))
      elseif(it.eq.5) then     !!! LJ intra
        if (z_LJassign) then 
           fun0       = 0.d0
        else 
           val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &               +(geom(ia,2)-geom(ib,2))**2
     &               +(geom(ia,3)-geom(ib,3))**2 )
           sigma=parfun(ip)
           val        = sigma/val
           valsq      = val*val
           val6       = valsq*valsq*valsq
           val12      = val6*val6
           rm1        = 1.d0/val
           rm2        = rm1*rm1
c          fun0       = 4.d0*(val12-val6)+1.0d0
           fun0       = 4.d0*(val12-val6)
        endif
      elseif(it.eq.6) then     !!! charge - charge interactions
        if (z_LJassign) then 
           fun0       = 0.d0
        else 
        val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &            +(geom(ia,2)-geom(ib,2))**2
     &            +(geom(ia,3)-geom(ib,3))**2 )
        rm1        = 1.d0/val
        fun0       = rm1
        endif

c       cross terms:  couplings section
 
        elseif(it.eq.8) then

          ifa=iflag(ip,1)
          ifb=iflag(ip,2)
          
c         Get geometrical values
c         Uncoupled coordinate 1 (UC1)
          it_uc1=iflag(ifa,0)
          ia    =iflag(ifa,1)
          ib    =iflag(ifa,2)
          ic    =iflag(ifa,3)
          id    =iflag(ifa,4)
          if (it_uc1.eq.1)  then     !! a = str
            val_uc1 = dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                      +(geom(ia,2)-geom(ib,2))**2
     &                      +(geom(ia,3)-geom(ib,3))**2 )
          elseif (it_uc1.eq.2)  then !! a = bend
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc1,'r')
          else                   !! a = tors(i,p)
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc1,'r')
            Call Shifta (val_uc1,parfun(ip))
          endif
c         Uncoupled coordinate 2 (UC2)
          it_uc2=iflag(ifb,0)
          ia    =iflag(ifb,1)
          ib    =iflag(ifb,2)
          ic    =iflag(ifb,3)
          id    =iflag(ifb,4)
          if (it_uc2.eq.1)   then    !! b = str
            val_uc2 = dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                      +(geom(ia,2)-geom(ib,2))**2
     &                      +(geom(ia,3)-geom(ib,3))**2 )
          elseif (it_uc2.eq.2)  then !! b = bend
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc2,'r')
          else                   !! b = tors(i,p)
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val_uc2,'r')
            Call Shifta (val_uc2,parfun2(ip))
          endif

c         Compute the functions and derivatives
          if (it_uc1.le.3)   then        !! a = str, bend, i.tors
           dval_uc1   = val_uc1-parfun(ip)
           if (it_uc2.le.3)  then        !! b = str, bend, i.tors
            dval_uc2   = val_uc2-parfun2(ip)
            fun0       = dval_uc1*dval_uc2      !! F        = (a-a0)*(b-b0)
           elseif(it_uc2.eq.4) then   !! b = p. torsion
            an_uc2  = iflag(ip,5)
            sdph2      = dsin(an_uc2*val_uc2-parfun2(ip))
c                      F = (a-a0)*(1+sin(n*b-g))
            fun0       = dval_uc1*(1.d0+sdph2)
           else                      !! not defined
            write(6,'(//8x,"ERROR: wrong coupled coordinate !",/)')
            stop
           endif
          elseif(it_uc1.eq.4) then       !! a = torsion
           if(it_uc2.le.3) then          !! b = str. , bend or i.tors
            write(6,'(//8x,
     &      "ERROR: first coupled coordinate should NOT be a torsion!",
     &       /)')
            stop
           elseif(it_uc2.eq.4) then      !! b = torsion
            an_uc1     = dble(iflag(ip,5))
            an_uc2     = dble(iflag(ip,6))
            sdph1      = dsin(an_uc1*val_uc1-parfun(ip))
            sdph2      = dsin(an_uc2*val_uc2-parfun2(ip))
c                      F = sin(na*a-ga)*sin(nb*b-gb)
            fun0       =   sdph1*sdph2
           else                      !! not defined
            write(6,'(//8x,"ERROR: wrong coupled coordinate !",/)')
            stop
           endif

          endif

        else
          write(6,'(//8x,"ERROR: wrong function",/)')
          stop
        endif

   40 format(//8x,'IRREPARABLE ERROR in subr. COMP_FUN1'
     &       ,/8x,'improper torsion',i5,5x,'atoms',4i4
     &       ,/8x,'angle , angle0',2f9.3,4x,'are differ too much'
     &       ,/8x,'function value is',f12.5,///)
   44 format(//8x,'ERROR: iflag(',i2,',1) =',i5,5x,'wrong',/)
      return
      end

c----------------------------------------------------------------------c
      Subroutine Shifta (val,val0)
c----------------------------------------------------------------------c
 
c     shifts (if needed) the first angle (val) by 2*pi  
c     so to make it similar to the second angle (val0)
c     if the gain is relevant.
c     It should be needed only for angles near to 180°
 
      implicit double precision (a-h,o-z)
      parameter (pig=3.141592653589793d0, pig2=pig*2.0d0)
 
      dif0=dabs(val-val0)
      if(dif0.lt.pig) return
 
      val1=val+pig2
      dif1=dabs(val1-val0)
      if(dabs(dif1).lt.dif0) then
        val=val1
        return
      endif
 
      val1=val-pig2
      dif1=dabs(val1-val0)
      if(dabs(dif1).lt.dif0) then
        val=val1
        return
      endif
 
      return
      end



