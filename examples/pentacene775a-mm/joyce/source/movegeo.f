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


!----------------------------------------------------------------------!
      Subroutine Do_Scan (par,parfun,parfun2,iflag,npar,zprt
     &                   ,iflagr,labric,nric,iun
     &                   ,scamin,scamax,scainc,iscaic,nscan
     &                   ,z_plottors,iungrace
     &                   ,geom,geoms,aname,nat,kbond,maxat,maxpar)
!----------------------------------------------------------------------!
 
!   compute the FF energy for geometries changed according to ...
 
!   The new values of the IRC's in xmove must be in angstroms or degrees.
!   The same values in valic are in bohr,radians
!   geom,geoms are in bohr
 
      implicit real*8(a-h,o-z)
      parameter (maxf=2000)
      dimension par(*),parfun(*),parfun2(*) 
      dimension iflag(maxpar,0:8),iflagr(maxpar,0:8)
     &         ,scaMin(*),scaMax(*),scainc(*),iscaic(*)
     &         ,geom(maxat,*),geoms(maxat,*),kbond(maxat,*)
     &         ,xmove(40),imove(40),vc(maxf),valicf0(maxf)
      dimension iflagx(maxf,0:8),fun0(maxf),fun1(maxf,2),fun2(maxf,3)
     &         ,valic(maxf)
      character aname(*)*5,labric(*)*15,labw*15
      logical zprt,z_plottors
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/
      data har2jkmol /2625.5d0/
      thrmove=-1.d-05
      ndim=maxat
 
      need=npar
      Call Check_Dim (need,maxf,'DO_SCAN: increase maxf')
 
      do i=1,npar
        do k=0,8
          iflagx(i,k)=iflag(i,k)
        end do
        vc(i)=boh2ang
        if(iflagx(i,8)/=2) vc(i)=todeg
      end do
 
      nmove=nscan
      do i=1,nscan
        xmove(i)=scaMin(i)
        imove(i)=iscaic(i)
      end do
 
      write(iun,'(/4x,"# RESULTS OF JOYCE SCAN")')
      nsicvar=0
      if(nscan>0) goto 30
 
!     no IC to be scanned : only one calculation
 
      Call Change_Geom (iflagr,nric,nmove,xmove,imove
     &                 ,geom,geoms,aname,nat,kbond,maxat,maxpar)
      Call Comp_Func (iflagx,fun0,fun1,fun2,parfun,parfun2,npar
     &               ,valic,nat,geom,maxat,maxf)
      ene=0.0d0
      do i=1,npar
        ene=ene+par(i)*fun0(i)
      end do
      write(6,'("  energy is",f14.4,"   kJ/mol")') ene*har2jkmol
      goto 80
 
!     more than 1 IC to be scanned
 
   30 continue
      do is=1,nscan
        if(scainc(is)>Abs(thrmove)) nsicvar=nsicvar+1
        if(nsicvar>2) nsicvar=2
        xmove(is)=scaMin(is)
        imove(is)=iscaic(is)
        irc=iscaic(is)
        a1=scaMin(is)
        a2=scaMax(is)
        labw=labric(irc)
        if(is<=nsicvar)write(iun,'("#",i4,2x,a,2f9.3)')irc,labw,a1,a2
        if(is>nsicvar)write(iun,'("#",i4,2x,a,f9.3)') irc,labw,a1
        if(is<=nsicvar)write(  6,'("#",i4,2x,a,2f9.3)')irc,labw,a1,a2
        if(is>nsicvar)write(  6,'("#",i4,2x,a,f9.3)') irc,labw,a1
      end do
 
      if(nsicvar==1) goto 40
      if(nsicvar==2) goto 60
      write(6,'(//"   SORRY: up to 2 ICs can be scanned",/)')
      goto 80
 
!     only one IC to be scanned
 
   40 write(iun,'("#  first  row =", a)') labric(imove(1))
      zprt = .false.
      if(zprt) write(6,'(" x variable is  ", a)') labric(imove(1))
      dx1=scainc(1)
      x1=scaMin(1)
      dowhile ((x1-scaMin(1))*(scaMax(1)-x1)>thrmove)
        xmove(1)=x1
        Call Change_Geom (iflagr,nric,nmove,xmove,imove
     &                   ,geom,geoms,aname,nat,kbond,maxat,maxpar)
        Call Comp_Func (iflagx,fun0,fun1,fun2,parfun,parfun2,npar
     &                 ,valic,nat,geom,maxat,maxf)
        ene=0.0d0
        do i=1,npar
          ene=ene+par(i)*fun0(i)
        end do
        write(iun,'(f12.3,f14.4)') x1,ene*har2jkmol
        if (z_plottors) then 
              write(iungrace,'(f9.1,5x,f12.3)') x1,ene*har2jkmol
        endif
        if(zprt) write(6,'(f12.3,f14.4)')   x1,ene*har2jkmol
        x1=x1+scainc(1)
      end do
      goto 80
 
!     two IC
 
   60 write(iun,'("#  first  row =", a)') labric(imove(1))
      write(iun,'("#  second row =", a)') labric(imove(2))
      zprt = .false.
      if(zprt) write(6,'(" x variable is  ", a)') labric(imove(1))
      if(zprt) write(6,'(" y variable is  ", a)') labric(imove(2))
 
      dx1=scainc(1)
      x1=scaMin(1)
      dowhile ((x1-scaMin(1))*(scaMax(1)-x1)>thrmove)
        xmove(1)=x1
 
      dx2=scainc(2)
      x2=scaMin(2)
      dowhile ((x2-scaMin(2))*(scaMax(2)-x2)>thrmove)
        xmove(2)=x2
 
        Call Change_Geom (iflagr,nric,nmove,xmove,imove
     &                   ,geom,geoms,aname,nat,kbond,maxat,maxpar)
        Call Comp_Func (iflagx,fun0,fun1,fun2,parfun,parfun2,npar
     &                 ,valic,nat,geom,maxat,maxf)
        ene=0.0d0
        do i=1,npar
          ene=ene+par(i)*fun0(i)
        end do
        write(iun,'(2f12.3,f14.4)') x1,x2,ene*har2jkmol
        if(zprt) write(6,'(2f12.3,f14.4)')   x1,x2,ene*har2jkmol
        x2=x2+scainc(2)
      end do
        write(iun,'(a)')
        x1=x1+scainc(1)
      end do
      write(iun,*) ' end'
 
   80 continue
      if (z_plottors) write(iungrace,'(a)') '&'
      return
      end Subroutine Do_Scan
 
c----------------------------------------------------------------------c
      Subroutine Change_Geom (iflagr,nric,nmove,xmove,imove
     &                       ,geom,geoms,aname,nat,kbond,maxat,maxpar)
c----------------------------------------------------------------------c
 
!   change the geometry, starting from geoms, according to the
!   isn.s in the .move arrays
 
!   - The new values of the IRC's in xmove must be in angstroms or degrees.
!   - geom,geoms are in bohr
 
      implicit real*8(a-h,o-z)
      dimension iflagr(maxpar,0:8),xmove(*),imove(*)
     &         ,geom(maxat,*),geoms(maxat,*),kbond(maxat,*)
      character*5 aname(*)
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/
      common /prmovegeo/imvpr
      thrmove=1.d-06
      ndim=maxat
      imvpr=0     !! test
 
      do i=1,nat
        do k=1,3
          geom(i,k)=geoms(i,k)
        end do
      end do
 
      if(nmove<=0) return
      do 20 im=1,nmove
        irc=imove(im)   !! the internal coord to be changed
        val=xmove(im)   !! the new value of the irc-th internal coord
        if(irc<0 .or. irc>nric) then
          write(6,'(//" ERROR in Change_Geom: invalid RIC",i5,/)') irc
          stop
        end if
        nc=iflagr(irc,8)
        vvv=val/todeg
        if(nc==2) vvv=val/boh2ang
        ia=iflagr(irc,1)
        ib=iflagr(irc,2)
        ic=iflagr(irc,3)
        id=iflagr(irc,4)
        Call Value_1IC (geoms,maxat,ia,ib,ic,id,val0,'r')
        del=vvv-val0
        if(Abs(del)<thrmove) goto 20
 
!       warning: these routines require the angle /vvv/ in degrees
 
        if(nc==2) then
          Call Change_Bond (nat,geom,kbond,ia,ib,vvv,aname,ndim)
        elseif(nc==3) then
          Call Change_Angle (nat,geom,kbond,ia,ib,ic,vvv,aname,ndim)
        else
          Call Change_Died (nat,geom,kbond,ia,ib,ic,id,vvv,aname,ndim)
        end if
   20 continue
   36 format(' CHANGE GEOM  ',4i4,'  old -> new',2f9.3)
 
!     check
 
      thrt=thrmove*5.d0
      nerr=0
      do 40 im=1,nmove
        irc=imove(im)   !! the internal coord to be changed
        val=xmove(im)   !! the new value of the irc-th internal coord
        nc=iflagr(irc,8)
        vvv=val/todeg
        if(nc==2) vvv=val/boh2ang
        ia=iflagr(irc,1)
        ib=iflagr(irc,2)
        ic=iflagr(irc,3)
        id=iflagr(irc,4)
        Call Value_1IC (geom,maxat,ia,ib,ic,id,valc,'r')
        if(ic>0) Call Shifta (valc,vvv)
        del=valc-vvv
        if(Abs(del)<thrt) goto 40
        write(6,46) ia,ib,ic,id,vvv,valc,Abs(vvv-valc),thrt
        nerr=nerr+1
   40 continue
   46 format(//9x,' CHANGE GEOM: something wrong for IC',4i4
     &       ,/9x,'              requested value ......',f14.8
     &       ,/9x,'              computed  value ......',f14.8
     &       ,/9x,'              difference ...........',f14.8
     &                                             ,' >',f14.8,//)
 
      if(nerr>0) stop
      return
      end Subroutine Change_Geom
 
c----------------------------------------------------------------------c
      Subroutine Change_Bond (nat,T,kbond,ia,ib,rab,aname,ndim)
c----------------------------------------------------------------------c
 
!     move the /ib/ atom at /rab/ distance from the /ia/ atom
!     along the ia-ib direction.
!     The atoms connected to /ib/ are also moved.
 
      implicit real*8(a-h,o-z)
      parameter (max=400)
      dimension T(ndim,*),kbond(ndim,*),imove(max)
      character s*1,aname(*)*5
      common /prmovegeo/imvpr
 
      if(nat>max)             stop 'Error Change_Bond: incr. max'
      if(ia<1 .or. ia>nat) stop 'Error Change_Bond: wrong ia'
      if(ib<1 .or. ib>nat) stop 'Error Change_Bond: wrong ib'
      if(ia==ib) stop 'ERROR: coincident atoms'
      if(kbond(ia,ib)==0) then
        write(6,16) ia,ib
      end if
   16 format(/'   INPUT WARNING: atoms',2i4,'   are not bonded',/)
 
      dx0=T(ib,1)-T(ia,1)
      dy0=T(ib,2)-T(ia,2)
      dz0=T(ib,3)-T(ia,3)
 
      Call Polar3 (dx0,dy0,dz0,r0,th,fi)
      dx=(rab-r0)*dsin(th)*dcos(fi)
      dy=(rab-r0)*dsin(th)*dsin(fi)
      dz=(rab-r0)*dcos(th)
      if(imvpr>0)
     &    write(6,24) r0,dx0,dy0,dz0,dx,dy,dz,dx+dx0,dy+dy0,dz+dz0
   24 format(/8x,'Old distance .......',f9.5
     &      ,/8x,'Old R(b-a) .........',3f9.5
     &      ,/8x,'Movement of b ......',3f9.5
     &      ,/8x,'New R(b-a) .........',3f9.5)
 
      Call Get_Conn (nat,kbond,ib,ia,imove,ndim)
 
      do i=1,nat
        if(imove(i)==1) then
          T(i,1)=T(i,1)+dx
          T(i,2)=T(i,2)+dy
          T(i,3)=T(i,3)+dz
        end if
      end do
 
      if(imvpr>0) then
        write(6,44) (imove(i),i=1,nat)
        write(6,'(/6x,"GEOMETRY WITH THE NEW DIST",2i4)') ia,ib
        do i=1,nat
          s=' '
          if(i==ia) s='A'
          if(i==ib) s='B'
          write(6,54) i,aname(i),(T(i,l),l=1,3),s,imove(i)
        end do
      end if
   44 format(/8x,'MOVED ATOMS   1=yes 0=no',/(6x,4(10i1,2x)))
   54 format(i5,3x,a,3f12.6,2x,a,i5)
 
      if(imove(ia)==1) then
        write(6,48)
        stop
      end if
   48 format(/8x,'ERROR: even the first atom has been moved',//)
 
      return
      end Subroutine Change_Bond
 
c----------------------------------------------------------------------c
      Subroutine Change_Angle (nat,T,kbond,ia,ib,ic,ang,aname,ndim)
c----------------------------------------------------------------------c
 
!     move the /ic/ atom so that the ia-ib-ic angle is /ang/.
!     All the atoms connected to ic are moved.
 
!     WARNING: in this routine ang is in radians
 
      implicit real*8(a-h,o-z)
      parameter (max=400)
      dimension T(ndim,*),kbond(ndim,*),imove(max),ibari(4)
     &         ,U(3,3),Ut(3)
      character s*1,aname(*)*5
      data todeg/57.2957795130823d0/
      common /prmovegeo/imvpr
 
!     controls
 
      if(nat>max)             stop 'Error Change_Angle: incr. max'
      if(ia<1 .or. ia>nat) stop 'Error Change_Angle: wrong ia'
      if(ib<1 .or. ib>nat) stop 'Error Change_Angle: wrong ib'
      if(ic<1 .or. ic>nat) stop 'Error Change_Angle: wrong ic'
      if(ia==ib .or. ia==ic .or.ib==ic) then
        write(6,12) ia,ib,ic
        stop
      end if
   12 format(/'   INPUT ERROR: atoms',3i4,'   are repeated',/)
      if(kbond(ia,ib)==0 .or. kbond(ib,ic)==0) then
        write(6,16) ia,ib,ic
      end if
   16 format(/'   INPUT WARNING: atoms',3i4,'   are not bonded',/)
 
!     puts the molecule in a convenient ref.
 
      ibari(1)=ib
      ibari(2)=0
      iaz=ia
      iaxz=ic
      Call Disponi (T,nat,ibari,iaz,iaxz,U,Ut,ndim)
 
!     looks for connection of the ic atom
 
      Call Get_Conn (nat,kbond,ic,ib,imove,ndim)
 
!     computes the rotation angles and moves ic connected atoms accordingly
 
      id=0
      Call Value_Ang (T,ndim,ia,ib,ic,id,ang0,'d')
      alf=0.d0
      bet=ang*todeg-ang0
      gam=0.d0
      Call Ruota_S (T,nat,alf,bet,gam,imove,ndim)
 
!     performs rotation + inverse traslation
 
      do i =1,nat
        tx = T(i,1)
        ty = T(i,2)
        tz = T(i,3)
        T(i,1) = U(1,1)*tx + U(2,1)*ty + U(3,1)*tz - Ut(1)
        T(i,2) = U(1,2)*tx + U(2,2)*ty + U(3,2)*tz - Ut(2)
        T(i,3) = U(1,3)*tx + U(2,3)*ty + U(3,3)*tz - Ut(3)
      end do
 
!     control
 
      Call Value_Ang (T,ndim,ia,ib,ic,id,angf,'d')
 
      if(imvpr>0) then
      write(6,44) (imove(i),i=1,nat)
      write(6,45) ang0,angf
   44 format(/8x,'MOVED ATOMS   1=yes 0=no',/(6x,4(10i1,2x)))
   45 format(/8x,'Old angle ..........',f12.5
     #      ,/8x,'New angle ..........',f12.5,4x,'recomputed')
 
      write(6,'(/6x,"GEOMETRY WITH THE NEW ANGLE",3i4)') ia,ib,ic
      do i=1,nat
        s=' '
        if(i==ia) s='A'
        if(i==ib) s='B'
        if(i==ic) s='C'
        write(6,54) i,aname(i),(T(i,l),l=1,3),s,imove(i)
      end do
      end if
   54 format(i5,3x,a,3f12.6,2x,a,i5)
 
      return
      end Subroutine Change_Angle
 
c----------------------------------------------------------------------c
      Subroutine Change_Died (nat,T,kbond,ia,ib,ic,id,ang,aname,ndim)
c----------------------------------------------------------------------c
 
!     Rotate the /id/ atom, along the ib-ic axis, so that the
!     ia-ib-ic-id dihedral is /ang/.
!     All the atoms connected to ic and id are also moved.
 
!     WARNING: in this routine ang is in radians
 
      implicit real*8(a-h,o-z)
      parameter (max=400)
      logical z_star
      dimension T(ndim,*),kbond(ndim,*),imove(max),ibari(4),iseq(nat)
     &         ,U(3,3),Ut(3),mbond(ndim,ndim)
      character s*1,aname(*)*5
      data todeg/57.2957795130823d0/
      common /prmovegeo/imvpr
 
!     control
 
      if(nat>max)             stop 'Error Change_Died: incr. max'
      if(ia<1 .or. ia>nat) stop 'Error Change_Died: wrong ia'
      if(ib<1 .or. ib>nat) stop 'Error Change_Died: wrong ib'
      if(ic<1 .or. ic>nat) stop 'Error Change_Died: wrong ic'
      if(id<1 .or. id>nat) stop 'Error Change_Died: wrong id'
      if(ia==ib .or. ia==ic .or.ia==id .or. ib==ic
     &            .or. ib==id .or. ic==id) then
        write(6,12) ia,ib,ic,id
        stop
      end if
   12 format(/'   INPUT ERROR: atoms',4i4,'   are repeated',/)
      z_star = .false.
      if(kbond(ia,ib)==0 .or. kbond(ib,ic)==0
     &                   .or. kbond(ic,id)==0) then
        z_star = .true.
        write(6,16) ia,ib,ic,id
        write(6,18) 
      end if
   16 format(/'   INPUT WARNING: atoms',4i4,'   are not bonded',/)
   18 format( '   FIRA for improper dihedrals will be adopted',/)
 
!     puts the molecule in a convenient ref.
 
      ibari(1)=ib
      ibari(2)=0
      iaz=ic
      iaxz=0
      Call Disponi (T,nat,ibari,iaz,iaxz,U,Ut,ndim)
 
!
!   Set the (modified) connection table. This is different for
!   proper and improper dihedrals.
!
!   NOTE: the treatment of improper dihedrals is rather naif,
!         it was adapted to a 6 atom aliphatic cycle to treat
!         boat/chair distorsion. It works for dihedrals defined
!         by the ordered sequence 1 3 5 4.
!
!                     1           A    
!                  6     2     6     2 
!                  5     3     C     B 
!                     4           D    
!
!   WARNING: move_bc set the movement of the atoms bonded to
!            B and C, other than the that connecting C to A
!            (i.e. atom 6) and B to A (i.e. atom 2)
!            move_bc=0   no movement
!            move_bc=1   full movement
!            move_bc=-1  half movement
!            
      move_bc=0
c     if(move_bc==0) then
c       write(6,'("move_bc=0  atoms conn. to B,C not moved")')
c     else if(move_bc>0) then
c       write(6,'("move_bc>0  atoms conn. to B,C moved")')
c     else if(move_bc<0) then
c       write(6,'("move_bc<0  atoms conn. to B,C moved half")')
c     end if
!
      mbond(1:nat,1:nat)=kbond(1:nat,1:nat)
      if (z_star) then
        mbond(ib,ic)=0
        mbond(ic,ib)=0
!    eliminate interactions B-2 and C-6
        Call Bonded_Pair (nat,kbond,ia,ic,ninter,iseq,ndim)
        if(ninter==1 .and. kbond(ic,id)>0) mbond(ic,iseq(2))=0
        Call Bonded_Pair (nat,kbond,ia,ib,ninter,iseq,ndim)
        if(ninter==1 .and. kbond(ib,id)>0) mbond(ib,iseq(2))=0
        Call Get_Conn (nat,mbond,id,it,imove,ndim)
        do i=1,nat
          if(i==id) cycle
          if(kbond(ic,i)>0 .and. imove(i)>0) imove(i)=move_bc
          if(kbond(ib,i)>0 .and. imove(i)>0) imove(i)=move_bc
        end do
c       write(6,'(/"IMOVE",/(10i2))') imove(1:nat)
      else
        mbond(ib,ic)=0
        mbond(ic,ib)=0
        it=0
        Call Get_Conn (nat,mbond,id,it,imove,ndim)
      endif
 
ccccccccccccccccccccccccc TEST BONDED_PAIR cccccccccccccccccccccc
c     ii=3 ; jj=9
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
c     ii=5 ; jj=nat
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
c     ii=35 ; jj=9
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
c     ii=35 ; jj=10
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
c     ii= 6 ; jj=12
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
c     ii= 9 ; jj=12
c     Call Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c     Call Bonded_Pair (nat,kbond,jj,ii,ninter,iseq,ndim)
ccccccccccccccccccccccccc TEST BONDED_PAIR cccccccccccccccccccccc
!     computes the rotation angles and moves ic connected atoms accordingly
 
      Call Value_Ang (T,ndim,ia,ib,ic,id,ang0,'d')
      alf=ang0-ang*todeg
      bet=0.0d0
      gam=0.0d0
      Call Ruota_S (T,nat,alf,bet,gam,imove,ndim)
 
!     performs rotation + inverse translation
 
      do i =1,nat
        tx = T(i,1)
        ty = T(i,2)
        tz = T(i,3)
        T(i,1) = U(1,1)*tx + U(2,1)*ty + U(3,1)*tz - Ut(1)
        T(i,2) = U(1,2)*tx + U(2,2)*ty + U(3,2)*tz - Ut(2)
        T(i,3) = U(1,3)*tx + U(2,3)*ty + U(3,3)*tz - Ut(3)
      end do
 
!     control
 
      Call Value_Ang (T,ndim,ia,ib,ic,id,angf,'d')
 
      imvpr = 0
      if(imvpr>0) then
      write(6,44) (imove(i),i=1,nat)
      write(6,45) ang0,angf
   44 format(/8x,'MOVED ATOMS   1=yes 0=no',/(6x,4(10i1,2x)))
   45 format(/8x,'Old dihedral .......',f12.5
     &      ,/8x,'New dihedral .......',f12.5,4x,'recomputed')
 
      write(6,'(/6x,"GEOMETRY WITH NEW DIHEDRAL",4i4)') ia,ib,ic,id
      do i=1,nat
        s=' '
        if(i==ia) s='A'
        if(i==ib) s='B'
        if(i==ic) s='C'
        if(i==id) s='D'
        write(6,54) i,aname(i),(T(i,l),l=1,3),s,imove(i)
      end do
   54 format(i5,3x,a,3f12.6,2x,a,i5)
      end if
 
      return
      end Subroutine Change_Died
 
c----------------------------------------------------------------------c
      Subroutine Get_Conn (nat,kbond,iat,ino,iconn,ndim)
c----------------------------------------------------------------------c
 
!     get a connectivity vector with respect to the atom /iat/
!      iconn(i)=1  atom i-th is connected with atom /iat/
!      iconn(i)=0  atom i-th is not connected with atom /iat/
 
      implicit real*8(a-h,o-z)
      dimension iconn(*),kbond(ndim,*),ibond(ndim,ndim)
 
      if(iat<1 .or. iat>nat) stop 'Error in Get_Conn: wrong iat'
 
      do i=1,nat
        iconn(i)=0
        do j=1,nat
          ibond(i,j)=kbond(i,j)
        end do
      end do
      iconn(iat)=1
      if(ino>0 .and. ino<=nat) ibond(iat,ino)=0
      if(ino>0 .and. ino<=nat) ibond(ino,iat)=0
c     do i=1,nat
c     write(6,'("ibond",i3,3x,5(10i1,1x))') i,ibond(i,1:nat)
c     end do
 
      icyc=0
  100 nconn=0
      icyc=icyc+1
      do i=1,nat
        if(iconn(i)==1) then
          do j=1,nat
            if(iconn(j)==0 .and. ibond(i,j)>0) then
              nconn=nconn+1
              iconn(j)=1
            end if
          end do
        end if
      end do
ccc   write(6,'(/"ICONN  cycle",2i3/(10i2))') icyc,nconn,iconn(1:nat)
      if(nconn>0) goto 100
 
      return
      end Subroutine Get_Conn
c----------------------------------------------------------------------c
      Subroutine Add_Connect (nat,kbond,ia,ib,ic,id,iconn,ndim)
c----------------------------------------------------------------------c
 
!   this routine modify the connectivity table (iconn) of an
!   improper dihedral geometry change, by adding the atoms
!   bonded to IB and IC, but those connecting IC with IA and
!   IB with IA, which should not be moved.
 
      implicit real*8(a-h,o-z)
      dimension iconn(*),kbond(ndim,*),ibond(ndim,ndim)
 
      if(ia<1 .or. ia>nat) stop 'Error in Add_Connect: wrong ia'
      if(ib<1 .or. ib>nat) stop 'Error in Add_Connect: wrong ib'
      if(ic<1 .or. ic>nat) stop 'Error in Add_Connect: wrong ic'
      if(id<1 .or. id>nat) stop 'Error in Add_Connect: wrong id'
      icyc=0
  100 nconn=0
      icyc=icyc+1
      do i=1,nat
        if(iconn(i)==1) then
          do j=1,nat
            if(iconn(j)==0 .and. ibond(i,j)>0) then
              nconn=nconn+1
              iconn(j)=1
            end if
          end do
        end if
      end do
ccc   write(6,'(/"ICONN  cycle",2i3/(10i2))') icyc,nconn,iconn(1:nat)
      if(nconn>0) goto 100
 
      return
      end Subroutine Add_Connect
 
 
c----------------------------------------------------------------------c
      Subroutine Bonded_Pair (nat,kbond,ii,jj,ninter,iseq,ndim)
c----------------------------------------------------------------------c
 
!     Given the connectivity tabel, this routine computes the n.
!     of atoms interposed between a given pair of atoms (ii and jj).
!     In case two or more bond sequences connect atoms ii and jj, 
!     the shortest one is sorted.
!
!  nat,kbond  #atoms and connectivity table
!  ii,jj      the pair of atoms
!  ninter     n. of atoms interposed between the atoms ii and jj
!             at return: ninter=0   atoms ii and jj are bonded
!                      : ninter>0   n. of atoms interposed
!                      : ninter=-1  atoms ii and jj are not linked, they
!                                   belong to non bonded fragments
!  iseq       the sequence of bonded atoms connecting ii and jj,
!             ninter+2 elements. iseq(1)=ii .... iseq(ninter+2)=jj
!
      implicit real*8(a-h,o-z)
      dimension kbond(ndim,*),iconn(nat),iseq(*)
      logical zconn
 
      if(ii<1 .or. ii>nat) stop 'Error in Bonded_Pair: wrong ii'
      if(jj<1 .or. jj>nat) stop 'Error in Bonded_Pair: wrong jj'
c     do i=1,nat
c       write(6,'("kbond",i3,3x,5(10i1,1x))') i,kbond(i,1:nat)
c     end do
 
!   special case: ii and jj are bonded to each other
      if(kbond(ii,jj)>0) then
        ninter=0
        write(6,24) ii,jj,ninter
        goto 20
      end if
!
c     write(6,'(//"ATOMS ===========",2i3)') ii,jj
      ninter=-1
      iconn(1:nat)=0
      iconn(ii)=1
      icyc=0
   10 nconn=0
      icyc=icyc+1
      do i=1,nat
        if(iconn(i)/=icyc) cycle
        do j=1,nat
          if(iconn(j)==0 .and. kbond(i,j)>0) then
            nconn=nconn+1
            iconn(j)=icyc+1
            if(kbond(j,jj)>0) then
              ninter=icyc
              ilast=j
              iconn(jj)=icyc+2
              write(6,24) ii,jj,ninter
              goto 20
            end if
ccc         write(6,*) 'j i iconn(j) ninter',j,i,iconn(j)
          end if
        end do
      end do
ccc   write(6,'(/"ICONN  cycle",2i3/(10i2))') icyc,nconn,iconn(1:nat)
      if(nconn>0) goto 10
   20 continue
      if(ninter<=0) return
   24 format( ' Bonded_Pair: atoms',2i4,4x,'#atoms interposed',i4)
!
!   now set the sequence connecting atoms ii and jj
!
      iseq(1:nat)=0
      if(ninter==1) then
        iseq(1)=ii
        iseq(2)=ilast
        iseq(3)=jj
        goto 40
      end if
!
!   eliminate those atoms which are not connected with both
!   lower and higher connection levels
!
   30 continue
c     write(6,364) iconn(1:nat)
c 364 format(/'ICONN' ,/(10i3,3x,10i3))
      ndrop=-1
      dowhile (ndrop/=0)
        ndrop=0
        do ic=2,ninter+1
          ic1=ic-1 ; ic2=ic+1
          do i=1,nat
            if(iconn(i)/=ic) cycle
            iok1=0
            do j=1,nat
              if(iconn(j)/=ic1) cycle
              if(kbond(i,j)>0) iok1=1
            end do
            iok2=0
            do j=1,nat
              if(iconn(j)/=ic2) cycle
              if(kbond(i,j)>0) iok2=1
            end do
            if(iok1+iok2==2) cycle
            iconn(i)=0
            ndrop=ndrop+1
c           write(6,*) 'eliminato atomo',i,'    at level',ic
          end do
        end do
        n=0
        do i=1,nat
          if(iconn(i)>0) n=n+1
        end do
c       write(6,*) ' -------------- remaining atoms',n
      end do
c     if(n/=ninter+2) write(6,36) n,iconn(1:nat)
c  36 format(//'WARNING in BONDED_PAIR: wrong n. remaining atoms',i4
c    #       ,/(10i3,3x,10i3))
!   in the case n>ninter+2 two or more paths are allowed,
!   then eliminate all but one paths
      if(n>ninter+2) then
        nelim=0
        do ic=2,ninter
          nc=0
          do i=1,nat
            if(iconn(i)==ic) nc=nc+1
          end do
          if(nc>1) then
            do i=1,nat
              if(iconn(i)/=ic) cycle
              iconn(i)=0 ; nelim=nelim+1
c             write(6,*) ' forced elim atom',i
              exit
            end do
          end if
          if(nelim>0) exit
        end do
        goto 30
      end if
!   pick up the remaining atoms forming the connecting sequence
      it=0
      do ic=1,ninter+2
        do i=1,nat
          if(iconn(i)/=ic) cycle
          it=it+1
          iseq(it)=i
          exit
        end do
      end do
!
   40 write(6,42) iseq(1:ninter+2)
   42 format( ' Sequence: ',10i4:/(11x,10i4))

      return
      end Subroutine Bonded_Pair

c----------------------------------------------------------------------c
      Subroutine Polar3 (x,y,z,r,t,f)
c----------------------------------------------------------------------c
 
!     computes the polar coordinates from the cartesian ones
 
      implicit real*8(a-h,o-z)
      parameter (pig=3.141592653589793238462d0, thrr=1.d-13)
      t=0.d0
      f=0.d0
      r=Sqrt(x*x+y*y+z*z)
      if(r<thrr) return
      xy=Sqrt(x*x+y*y)
      st=xy/r
      ct=z/r
      t=Getang(st,ct)
      if(xy<thrr) return
      sf=y/xy
      cf=x/xy
      f=Getang(sf,cf)
      return
      end Subroutine Polar3
 
c----------------------------------------------------------------------c
      Subroutine Disponi (coord,nat,ibari,iaz,iaxz,T,Tr,maxat)
c----------------------------------------------------------------------c
 
!     traslates e rotates molecular nuclei so that:
!      1. origin in CM
!      2. atom /iaz/ on z axis
!      3. atom /iaxz/ on xz plane
 
      implicit real*8(a-h,o-y), logical(z)
      dimension coord(maxat,*),ibari(*),T(3,*),Tr(*),Ra(999)
      parameter (pig=3.141592653589793238d0, fang=180.d0/pig)
      logical zdeb
      zdeb=.false.
 
      if(nat<1)  stop 'ERROR in Disponi: nat <=0'
 
!     1. trasl
 
      cx=0.d0
      cy=0.d0
      cz=0.d0
      n=1
      dowhile (ibari(n)>0 .and. ibari(n)<=nat)
        ib=ibari(n)
        cx=cx+coord(ib,1)
        cy=cy+coord(ib,2)
        cz=cz+coord(ib,3)
        n=n+1
      end do
      if(n==1) goto 15
      Tr(1)=-cx/Dble(n-1)
      Tr(2)=-cy/Dble(n-1)
      Tr(3)=-cz/Dble(n-1)
      do i=1,nat
        ar=0.d0
        do k=1,3
          coord(i,k)=coord(i,k)+Tr(k)
          ar=ar+coord(i,k)**2
        end do
        Ra(i)=Sqrt(ar)
      end do
      if(zdeb) write(6,12) nat,(Tr(k),k=1,3)
     &    ,(i,coord(i,4),(coord(i,k),k=1,3),Ra(i),i=1,nat)
   12 format(/'  ------ GEOMETRY after translation',i6,' -------'
     &      ,/'  transl vec =',3f11.4
     &      ,/20x,'X',10x,'Y',10x,'Z',11x,'R'
     &    ,/(i6,f5.0,2x,3f11.5,f11.4))
 
!     2. rot atom iaz
!     3. rot atom iaxz
 
   15 Call Polar3 (coord(iaz,1),coord(iaz,2),coord(iaz,3),r,th,ph)
      alpha=0.d0            !! rotaz about  Z
      beta=th*fang          !! rotaz about  Y'
      gamma=ph*fang         !! rotaz about  Z"
      Call Gimme_RotMat (alpha,beta,gamma,T)
      if(zdeb) write(6,'(" alf bet gam",3f10.4)') alpha,beta,gamma
 
      if(iaxz<=0) goto 20
      ax=coord(iaxz,1)
      ay=coord(iaxz,2)
      az=coord(iaxz,3)
      bx=T(1,1)*ax+T(1,2)*ay+T(1,3)*az
      by=T(2,1)*ax+T(2,2)*ay+T(2,3)*az
      bz=T(3,1)*ax+T(3,2)*ay+T(3,3)*az
      Call Polar3 (bx,by,bz,r,th,ph)
      alpha=(ph+pig)*fang       !! rotaz about Z (for xzplane)
      Call Gimme_RotMat (alpha,beta,gamma,T)
 
   20 do i=1,nat
        ax=coord(i,1)
        ay=coord(i,2)
        az=coord(i,3)
        coord(i,1)=T(1,1)*ax+T(1,2)*ay+T(1,3)*az
        coord(i,2)=T(2,1)*ax+T(2,2)*ay+T(2,3)*az
        coord(i,3)=T(3,1)*ax+T(3,2)*ay+T(3,3)*az
        ar=0.d0
        do k=1,3
          ar=ar+coord(i,k)**2
        end do
        Ra(i)=Sqrt(ar)
      end do
 
      if(zdeb) write(6,'(/" ROTATION MATRIX",/(3f18.13))')
     &                          ((T(i,j),j=1,3),i=1,3)
      if(zdeb) write(6,44) nat
     &    ,(i,coord(i,4),(coord(i,k),k=1,3),Ra(i),i=1,nat)
   44 format(/'  ------ GEOMETRY after transl+rot ',i6,' -------'
     &      ,/20x,'X',10x,'Y',10x,'Z',11x,'R'
     &    ,/(i6,f5.0,2x,3f11.5,f11.4))
 
      return
      end Subroutine Disponi
 
c----------------------------------------------------------------------c
      Function Getang (ss,cc)
c----------------------------------------------------------------------c
 
!     given ss=sen(a) and cc=cos(a), computes a in [-pig ... pig]
 
      implicit real*8(a-h,o-z)
      parameter (pig=3.141592653589793238462d0,tol=1.d-14)
 
      if(Abs(cc)<tol) then
        if(ss>0.d0) f= 0.5d0*pig
        if(ss<0.d0) f=-0.5d0*pig
      else
        f=datan(ss/cc)
        if(cc<0.0d0) f=f+pig
        if(f>pig) f=f-pig-pig
      end if
      Getang=f
 
!     internal check
 
      ss1=dsin(f)
      cc1=dcos(f)
      ers=Abs(ss-ss1)
      erc=Abs(cc-cc1)
      if(ers>1.d-12 .or. erc>1.d-12) then
        fg=f*180/pig
        write(6,34) ss,cc,ss1,cc1,fg
        stop
      end if
   34 format(//9x,'SEVERE  INTERNAL  ERROR  in  GETANG'
     &    ,/4x,'input      sin cos  ',2f17.12
     &    ,/4x,'recomputed sin cos  ',2f17.12
     &    ,/4x,'computed angle      ',f14.3,//)
      return
      end Function Getang
 
c----------------------------------------------------------------------c
      Subroutine Ruota_S (ca,nat,alf,bet,gam,irot,maxat)
c----------------------------------------------------------------------c
 
!     molecular rotations (with selected atoms)
!     UPDATE jan 2015: the atoms with irot<0 are rotated by
!                      the half values of alf,bet,gam
!
      implicit real*8(a-h,o-z)
      parameter (pig=3.141592653589793238462d0)
      dimension ca(maxat,*),eu(3,3),irot(*)
 
      Call Gimme_RotMat (alf,bet,gam,eu)
      do i =1,nat
        if(irot(i)>0) then
          tx = ca(i,1)
          ty = ca(i,2)
          tz = ca(i,3)
          ca(i,1) = eu(1,1)*tx + eu(1,2)*ty + eu(1,3)*tz
          ca(i,2) = eu(2,1)*tx + eu(2,2)*ty + eu(2,3)*tz
          ca(i,3) = eu(3,1)*tx + eu(3,2)*ty + eu(3,3)*tz
        end if
      end do
 
      alfm=alf*0.5d0
      betm=bet*0.5d0
      gamm=gam*0.5d0
      Call Gimme_RotMat (alfm,betm,gamm,eu)
      do i =1,nat
        if(irot(i)<0) then
          tx = ca(i,1)
          ty = ca(i,2)
          tz = ca(i,3)
          ca(i,1) = eu(1,1)*tx + eu(1,2)*ty + eu(1,3)*tz
          ca(i,2) = eu(2,1)*tx + eu(2,2)*ty + eu(2,3)*tz
          ca(i,3) = eu(3,1)*tx + eu(3,2)*ty + eu(3,3)*tz
        end if
      end do
 
      return
      return
      end Subroutine Ruota_S
 
c----------------------------------------------------------------------c
      Subroutine Gimme_RotMat (alf,bet,gam,eu)
c----------------------------------------------------------------------c
 
!   computes rotation matrix of the three Euler angles  (deg).
 
!   WARNING: this version does not control angle beta
 
!     rotation OF the ref. system (passive mode)
!      1)  alf  = clockwise rotation about z
!      2)  bet  = clockwise rotation about y'
!      3)  gam  = clockwise rotation about z"
 
!     rotation  IN the  ref. system  (active mode)
!      1)  gam  =  clockwise rotation about z
!      2)  bet  =  clockwise rotation about y
!      3)  alf  =  clockwise rotation about z
 
!     use column matrix
 
      implicit real*8(a-h,o-z)
      parameter (fang=3.141592653589793238462d0/180.d0)
      dimension eu(3,3)
 
      ralf = alf*fang
      rbet = bet*fang
      rgam = gam*fang
      calf = dcos(ralf)
      cbet = dcos(rbet)
      cgam = dcos(rgam)
      salf = dsin(ralf)
      sbet = dsin(rbet)
      sgam = dsin(rgam)
 
      eu(1,1) =  cgam*cbet*calf - salf*sgam
      eu(1,2) =  calf*cbet*sgam + salf*cgam
      eu(1,3) = -calf*sbet
      eu(2,1) = -salf*cbet*cgam - calf*sgam
      eu(2,2) = -salf*cbet*sgam + calf*cgam
      eu(2,3) =  salf*sbet
      eu(3,1) =  sbet*cgam
      eu(3,2) =  sgam*sbet
      eu(3,3) =  cbet
 
      return
      end Subroutine Gimme_RotMat
 
c----------------------------------------------------------------------c
      Subroutine Read_Xyz (nat,geom,aname,maxat,iun)
c----------------------------------------------------------------------c
 
!     reads xyz fortmat geometry
!                           :  C   10. 4.0 -2.0
!                           :  6.  10. 4.0 -2.0
 
      implicit real*8(a-h,o-z)
      dimension geom(maxat,*)
      character aname(*)*5,card*90,word*50
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
 
      write(6,48)
      read(iun,*) nat0
      if(nat0/=nat) goto 210
      read(iun,'(a)') card
      do ia=1,nat
        read(iun,'(a)',end=200) card
        Call SelWrd (card,word,1,lw,ipos)
        read(card(ipos+lw+1:80),*) (geom(ia,i),i=1,3)
        write(6,50) ia,aname(ia),geom(ia,4),(geom(ia,i),i=1,3)
        geom(ia,1)=geom(ia,1)/boh2ang
        geom(ia,2)=geom(ia,2)/boh2ang
        geom(ia,3)=geom(ia,3)/boh2ang
      end do
      return
   48 format(/5x,'New Reference geometry for scan (warning: no check)')
   50 format(i4,3x,a,f5.0,3x,3f12.6)
 
  200 write(6,'(//" END-OF-FILE on input file",//)')
      stop
  210 write(6,'(//" ERROR: wrong /nat/ on geometry file",//)')
      stop
      end Subroutine Read_Xyz

c     ------------------------------------------------------------
      Subroutine PrintGtors(iungrace,icountsc)
c     ------------------------------------------------------------

      implicit real*8 (a-h,o-z)
      character*20 tmpfile
      character*3  num
      logical   z_LJassign
      dimension delta(1000),E(1000),Ecorr(1000)
      common /ljassign/z_LJassign

c     reading scanned QM data
      if (icountsc .lt. 10)   write(num,'(i1)') icountsc
      if (icountsc .lt. 100)  write(num,'(i2)') icountsc
      Call Compact(num,lw) 
      tmpfile = 'qmscan.'//num(1:lw)//'.dat'
      open(23,file=tmpfile,form='formatted',status='old')
      read(23,*) 
      iqm = 1 
5     if (z_LJassign) then
          read(23,*,end=10) delta(iqm), Ecorr(iqm), E(iqm)
      else
          read(23,*,end=10) delta(iqm),E(iqm)
      endif
      iqm = iqm + 1
      goto 5

10    continue
      Nqm = iqm - 1
c     getting maximum for y-scale
      emax = -1.d6
      do i = 1, Nqm
         if (E(i) .gt. emax) emax = E(i)
      enddo
      de = emax*0.35
      emax = emax  + de

      if (z_LJassign) then ! printing grace file for torsion  with  LJ assign

       write(iungrace,'(a)') '# Grace project file'
       write(iungrace,'(a)') '#'
       write(iungrace,'(a)') '@version 50125'
       write(iungrace,'(a)') '@page size 792, 612'
       write(iungrace,'(a)') '@page scroll 5%'
       write(iungrace,'(a)') '@page inout 5%'
       write(iungrace,'(a)') '@link page off'
       write(iungrace,'(a)') 
     &   '@map font 0 to "Times-Roman", "Times-Roman"'
       write(iungrace,'(a)') 
     &   '@map font 1 to "Times-Italic", "Times-Italic"'
       write(iungrace,'(a)') 
     &   '@map font 2 to "Times-Bold", "Times-Bold"'
       write(iungrace,'(a)') 
     &   '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"'
       write(iungrace,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
       write(iungrace,'(a)') 
     &   '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"'
       write(iungrace,'(a)') 
     &   '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"'
       write(iungrace,'(a)') 
     & '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"'
       write(iungrace,'(a)') '@map font 8 to "Courier", "Courier"'
       write(iungrace,'(a)') 
     &   '@map font 9 to "Courier-Oblique", "Courier-Oblique"'
       write(iungrace,'(a)') 
     &   '@map font 10 to "Courier-Bold", "Courier-Bold"'
       write(iungrace,'(a)') 
     &   '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"'
       write(iungrace,'(a)') '@map font 12 to "Symbol", "Symbol"'
       write(iungrace,'(a)') 
     &   '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
       write(iungrace,'(a)') '@map color 0 to (255, 255, 255), "white"'
       write(iungrace,'(a)') '@map color 1 to (0, 0, 0), "black"'
       write(iungrace,'(a)') '@map color 2 to (255, 0, 0), "red"'
       write(iungrace,'(a)') '@map color 3 to (0, 255, 0), "green"'
       write(iungrace,'(a)') '@map color 4 to (0, 0, 255), "blue"'
       write(iungrace,'(a)') '@map color 5 to (255, 255, 0), "yellow"'
       write(iungrace,'(a)') '@map color 6 to (188, 143, 143), "brown"'
       write(iungrace,'(a)') '@map color 7 to (220, 220, 220), "grey"'
       write(iungrace,'(a)') '@map color 8 to (148, 0, 211), "violet"'
       write(iungrace,'(a)') '@map color 9 to (0, 255, 255), "cyan"'
       write(iungrace,'(a)') '@map color 10 to (255, 0, 255), "magenta"'
       write(iungrace,'(a)') '@map color 11 to (255, 165, 0), "orange"'
       write(iungrace,'(a)') '@map color 12 to (114, 33, 188), "indigo"'
       write(iungrace,'(a)') '@map color 13 to (103, 7, 72), "maroon"'
       write(iungrace,'(a)') 
     &   '@map color 14 to (64, 224, 208), "turquoise"'
       write(iungrace,'(a)') '@map color 15 to (0, 139, 0), "green4"'
       write(iungrace,'(a)') '@reference date 0'
       write(iungrace,'(a)') '@date wrap off'
       write(iungrace,'(a)') '@date wrap year 1950'
       write(iungrace,'(a)') '@default linewidth 1.0'
       write(iungrace,'(a)') '@default linestyle 1'
       write(iungrace,'(a)') '@default color 1'
       write(iungrace,'(a)') '@default pattern 1'
       write(iungrace,'(a)') '@default font 0'
       write(iungrace,'(a)') '@default char size 1.000000'
       write(iungrace,'(a)') '@default symbol size 1.000000'
       write(iungrace,'(a)') '@default sformat "%.8g"'
       write(iungrace,'(a)') '@background color 0'
       write(iungrace,'(a)') '@page background fill on'
       write(iungrace,'(a)') '@timestamp off'
       write(iungrace,'(a)') '@timestamp 0.03, 0.03'
       write(iungrace,'(a)') '@timestamp color 1'
       write(iungrace,'(a)') '@timestamp rot 0'
       write(iungrace,'(a)') '@timestamp font 0'
       write(iungrace,'(a)') '@timestamp char size 1.000000'
       write(iungrace,'(a)') '@r0 off'
       write(iungrace,'(a)') '@link r0 to g0'
       write(iungrace,'(a)') '@r0 type above'
       write(iungrace,'(a)') '@r0 linestyle 1'
       write(iungrace,'(a)') '@r0 linewidth 1.0'
       write(iungrace,'(a)') '@r0 color 1'
       write(iungrace,'(a)') '@r0 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r1 off'
       write(iungrace,'(a)') '@link r1 to g0'
       write(iungrace,'(a)') '@r1 type above'
       write(iungrace,'(a)') '@r1 linestyle 1'
       write(iungrace,'(a)') '@r1 linewidth 1.0'
       write(iungrace,'(a)') '@r1 color 1'
       write(iungrace,'(a)') '@r1 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r2 off'
       write(iungrace,'(a)') '@link r2 to g0'
       write(iungrace,'(a)') '@r2 type above'
       write(iungrace,'(a)') '@r2 linestyle 1'
       write(iungrace,'(a)') '@r2 linewidth 1.0'
       write(iungrace,'(a)') '@r2 color 1'
       write(iungrace,'(a)') '@r2 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r3 off'
       write(iungrace,'(a)') '@link r3 to g0'
       write(iungrace,'(a)') '@r3 type above'
       write(iungrace,'(a)') '@r3 linestyle 1'
       write(iungrace,'(a)') '@r3 linewidth 1.0'
       write(iungrace,'(a)') '@r3 color 1'
       write(iungrace,'(a)') '@r3 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r4 off'
       write(iungrace,'(a)') '@link r4 to g0'
       write(iungrace,'(a)') '@r4 type above'
       write(iungrace,'(a)') '@r4 linestyle 1'
       write(iungrace,'(a)') '@r4 linewidth 1.0'
       write(iungrace,'(a)') '@r4 color 1'
       write(iungrace,'(a)') '@r4 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@g0 on'
       write(iungrace,'(a)') '@g0 hidden false'
       write(iungrace,'(a)') '@g0 type XY'
       write(iungrace,'(a)') '@g0 stacked false'
       write(iungrace,'(a)') '@g0 bar hgap 0.000000'
       write(iungrace,'(a)') '@g0 fixedpoint off'
       write(iungrace,'(a)') '@g0 fixedpoint type 0'
       write(iungrace,'(a)') '@g0 fixedpoint xy 0.000000, 0.000000'
       write(iungrace,'(a)') '@g0 fixedpoint format general general'
       write(iungrace,'(a)') '@g0 fixedpoint prec 6, 6'
       write(iungrace,'(a)') '@with g0'
       write(iungrace,'(a24,1x,f12.0)') '@    world -180, 0, 180,',emax
       write(iungrace,'(a)') '@    stack world 0, 0, 0, 0'
       write(iungrace,'(a)') '@    znorm 1'
       write(iungrace,'(a)') 
     &   '@    view 0.150000, 0.150000, 1.150000, 0.850000'
       write(iungrace,'(a)') '@    title ""'
       write(iungrace,'(a)') '@    title font 0'
       write(iungrace,'(a)') '@    title size 1.500000'
       write(iungrace,'(a)') '@    title color 1'
       write(iungrace,'(a)') '@    subtitle ""'
       write(iungrace,'(a)') '@    subtitle font 0'
       write(iungrace,'(a)') '@    subtitle size 1.000000'
       write(iungrace,'(a)') '@    subtitle color 1'
       write(iungrace,'(a)') '@    xaxes scale Normal'
       write(iungrace,'(a)') '@    yaxes scale Normal'
       write(iungrace,'(a)') '@    xaxes invert off'
       write(iungrace,'(a)') '@    yaxes invert off'
       write(iungrace,'(a)') '@    xaxis  on'
       write(iungrace,'(a)') '@    xaxis  type zero false'
       write(iungrace,'(a)') '@    xaxis  offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    xaxis  bar on'
       write(iungrace,'(a)') '@    xaxis  bar color 1'
       write(iungrace,'(a)') '@    xaxis  bar linestyle 1'
       write(iungrace,'(a)') '@    xaxis  bar linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  label "\xd\f{} (degr)"'
       write(iungrace,'(a)') '@    xaxis  label layout para'
       write(iungrace,'(a)') '@    xaxis  label place auto'
       write(iungrace,'(a)') '@    xaxis  label char size 1.000000'
       write(iungrace,'(a)') '@    xaxis  label font 4'
       write(iungrace,'(a)') '@    xaxis  label color 1'
       write(iungrace,'(a)') '@    xaxis  label place normal'
       write(iungrace,'(a)') '@    xaxis  tick on'
       write(iungrace,'(a)') '@    xaxis  tick major 30'
       write(iungrace,'(a)') '@    xaxis  tick minor ticks 1'
       write(iungrace,'(a)') '@    xaxis  tick default 6'
       write(iungrace,'(a)') '@    xaxis  tick place rounded true'
       write(iungrace,'(a)') '@    xaxis  tick in'
       write(iungrace,'(a)') '@    xaxis  tick major size 1.000000'
       write(iungrace,'(a)') '@    xaxis  tick major color 1'
       write(iungrace,'(a)') '@    xaxis  tick major linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  tick major linestyle 1'
       write(iungrace,'(a)') '@    xaxis  tick major grid off'
       write(iungrace,'(a)') '@    xaxis  tick minor color 1'
       write(iungrace,'(a)') '@    xaxis  tick minor linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  tick minor linestyle 1'
       write(iungrace,'(a)') '@    xaxis  tick minor grid off'
       write(iungrace,'(a)') '@    xaxis  tick minor size 0.500000'
       write(iungrace,'(a)') '@    xaxis  ticklabel on'
       write(iungrace,'(a)') '@    xaxis  ticklabel format general'
       write(iungrace,'(a)') '@    xaxis  ticklabel prec 5'
       write(iungrace,'(a)') '@    xaxis  ticklabel formula ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel append ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel prepend ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel angle 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel skip 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel stagger 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel place normal'
       write(iungrace,'(a)') '@    xaxis  ticklabel offset auto'
       write(iungrace,'(a)') 
     &   '@    xaxis  ticklabel offset 0.000000 , 0.010000'
       write(iungrace,'(a)') '@    xaxis  ticklabel start type auto'
       write(iungrace,'(a)') '@    xaxis  ticklabel start 0.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel stop type auto'
       write(iungrace,'(a)') '@    xaxis  ticklabel stop 0.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel char size 1.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel font 4'
       write(iungrace,'(a)') '@    xaxis  ticklabel color 1'
       write(iungrace,'(a)') '@    xaxis  tick place both'
       write(iungrace,'(a)') '@    xaxis  tick spec type none'
       write(iungrace,'(a)') '@    yaxis  on'
       write(iungrace,'(a)') '@    yaxis  type zero false'
       write(iungrace,'(a)') '@    yaxis  offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    yaxis  bar on'
       write(iungrace,'(a)') '@    yaxis  bar color 1'
       write(iungrace,'(a)') '@    yaxis  bar linestyle 1'
       write(iungrace,'(a)') '@    yaxis  bar linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  label "\xD\f{}E (kJ/mol)"'
       write(iungrace,'(a)') '@    yaxis  label layout para'
       write(iungrace,'(a)') '@    yaxis  label place auto'
       write(iungrace,'(a)') '@    yaxis  label char size 1.000000'
       write(iungrace,'(a)') '@    yaxis  label font 4'
       write(iungrace,'(a)') '@    yaxis  label color 1'
       write(iungrace,'(a)') '@    yaxis  label place normal'
       write(iungrace,'(a)') '@    yaxis  tick on'
       write(iungrace,'(a)') '@    yaxis  tick major 5'
       write(iungrace,'(a)') '@    yaxis  tick minor ticks 1'
       write(iungrace,'(a)') '@    yaxis  tick default 6'
       write(iungrace,'(a)') '@    yaxis  tick place rounded true'
       write(iungrace,'(a)') '@    yaxis  tick in'
       write(iungrace,'(a)') '@    yaxis  tick major size 1.000000'
       write(iungrace,'(a)') '@    yaxis  tick major color 1'
       write(iungrace,'(a)') '@    yaxis  tick major linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  tick major linestyle 1'
       write(iungrace,'(a)') '@    yaxis  tick major grid off'
       write(iungrace,'(a)') '@    yaxis  tick minor color 1'
       write(iungrace,'(a)') '@    yaxis  tick minor linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  tick minor linestyle 1'
       write(iungrace,'(a)') '@    yaxis  tick minor grid off'
       write(iungrace,'(a)') '@    yaxis  tick minor size 0.500000'
       write(iungrace,'(a)') '@    yaxis  ticklabel on'
       write(iungrace,'(a)') '@    yaxis  ticklabel format general'
       write(iungrace,'(a)') '@    yaxis  ticklabel prec 5'
       write(iungrace,'(a)') '@    yaxis  ticklabel formula ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel append ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel prepend ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel angle 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel skip 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel stagger 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel place normal'
       write(iungrace,'(a)') '@    yaxis  ticklabel offset auto'
       write(iungrace,'(a)') 
     &   '@    yaxis  ticklabel offset 0.000000 , 0.010000'
       write(iungrace,'(a)') '@    yaxis  ticklabel start type auto'
       write(iungrace,'(a)') '@    yaxis  ticklabel start 0.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel stop type auto'
       write(iungrace,'(a)') '@    yaxis  ticklabel stop 0.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel char size 1.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel font 4'
       write(iungrace,'(a)') '@    yaxis  ticklabel color 1'
       write(iungrace,'(a)') '@    yaxis  tick place both'
       write(iungrace,'(a)') '@    yaxis  tick spec type none'
       write(iungrace,'(a)') '@    altxaxis  off'
       write(iungrace,'(a)') '@    altyaxis  off'
       write(iungrace,'(a)') '@    legend on'
       write(iungrace,'(a)') '@    legend loctype view'
       write(iungrace,'(a)') '@    legend 0.85, 0.8'
       write(iungrace,'(a)') '@    legend box color 1'
       write(iungrace,'(a)') '@    legend box pattern 1'
       write(iungrace,'(a)') '@    legend box linewidth 1.0'
       write(iungrace,'(a)') '@    legend box linestyle 1'
       write(iungrace,'(a)') '@    legend box fill color 0'
       write(iungrace,'(a)') '@    legend box fill pattern 1'
       write(iungrace,'(a)') '@    legend font 0'
       write(iungrace,'(a)') '@    legend char size 1.000000'
       write(iungrace,'(a)') '@    legend color 1'
       write(iungrace,'(a)') '@    legend length 4'
       write(iungrace,'(a)') '@    legend vgap 1'
       write(iungrace,'(a)') '@    legend hgap 1'
       write(iungrace,'(a)') '@    legend invert false'
       write(iungrace,'(a)') '@    frame type 0'
       write(iungrace,'(a)') '@    frame linestyle 1'
       write(iungrace,'(a)') '@    frame linewidth 1.0'
       write(iungrace,'(a)') '@    frame color 1'
       write(iungrace,'(a)') '@    frame pattern 1'
       write(iungrace,'(a)') '@    frame background color 0'
       write(iungrace,'(a)') '@    frame background pattern 0'
       write(iungrace,'(a)') '@    s0 hidden false'
       write(iungrace,'(a)') '@    s0 type xy'
       write(iungrace,'(a)') '@    s0 symbol 2'
       write(iungrace,'(a)') '@    s0 symbol size 1.100000'
       write(iungrace,'(a)') '@    s0 symbol color 7'
       write(iungrace,'(a)') '@    s0 symbol pattern 1'
       write(iungrace,'(a)') '@    s0 symbol fill color 7'
       write(iungrace,'(a)') '@    s0 symbol fill pattern 1'
       write(iungrace,'(a)') '@    s0 symbol linewidth 1.0'
       write(iungrace,'(a)') '@    s0 symbol linestyle 1'
       write(iungrace,'(a)') '@    s0 symbol char 65'
       write(iungrace,'(a)') '@    s0 symbol char font 0'
       write(iungrace,'(a)') '@    s0 symbol skip 0'
       write(iungrace,'(a)') '@    s0 line type 0'
       write(iungrace,'(a)') '@    s0 line linestyle 1'
       write(iungrace,'(a)') '@    s0 line linewidth 1.0'
       write(iungrace,'(a)') '@    s0 line color 7'
       write(iungrace,'(a)') '@    s0 line pattern 1'
       write(iungrace,'(a)') '@    s0 baseline type 0'
       write(iungrace,'(a)') '@    s0 baseline off'
       write(iungrace,'(a)') '@    s0 dropline off'
       write(iungrace,'(a)') '@    s0 fill type 0'
       write(iungrace,'(a)') '@    s0 fill rule 0'
       write(iungrace,'(a)') '@    s0 fill color 1'
       write(iungrace,'(a)') '@    s0 fill pattern 1'
       write(iungrace,'(a)') '@    s0 avalue off'
       write(iungrace,'(a)') '@    s0 avalue type 2'
       write(iungrace,'(a)') '@    s0 avalue char size 1.000000'
       write(iungrace,'(a)') '@    s0 avalue font 0'
       write(iungrace,'(a)') '@    s0 avalue color 1'
       write(iungrace,'(a)') '@    s0 avalue rot 0'
       write(iungrace,'(a)') '@    s0 avalue format general'
       write(iungrace,'(a)') '@    s0 avalue prec 3'
       write(iungrace,'(a)') '@    s0 avalue prepend ""'
       write(iungrace,'(a)') '@    s0 avalue append ""'
       write(iungrace,'(a)') '@    s0 avalue offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    s0 errorbar on'
       write(iungrace,'(a)') '@    s0 errorbar place both'
       write(iungrace,'(a)') '@    s0 errorbar color 7'
       write(iungrace,'(a)') '@    s0 errorbar pattern 1'
       write(iungrace,'(a)') '@    s0 errorbar size 1.000000'
       write(iungrace,'(a)') '@    s0 errorbar linewidth 1.0'
       write(iungrace,'(a)') '@    s0 errorbar linestyle 1'
       write(iungrace,'(a)') '@    s0 errorbar riser linewidth 1.0'
       write(iungrace,'(a)') '@    s0 errorbar riser linestyle 1'
       write(iungrace,'(a)') '@    s0 errorbar riser clip off'
       write(iungrace,'(a)') 
     &   '@    s0 errorbar riser clip length 0.100000'
       write(iungrace,'(a)') '@    s0 comment "Cols 1:3"'
       write(iungrace,'(a)') '@    s0 legend  "QM"'
       write(iungrace,'(a)') '@    s1 hidden false'
       write(iungrace,'(a)') '@    s1 type xy'
       write(iungrace,'(a)') '@    s1 symbol 1'
       write(iungrace,'(a)') '@    s1 symbol size 0.850000'
       write(iungrace,'(a)') '@    s1 symbol color 15'
       write(iungrace,'(a)') '@    s1 symbol pattern 1'
       write(iungrace,'(a)') '@    s1 symbol fill color 15'
       write(iungrace,'(a)') '@    s1 symbol fill pattern 1'
       write(iungrace,'(a)') '@    s1 symbol linewidth 1.0'
       write(iungrace,'(a)') '@    s1 symbol linestyle 1'
       write(iungrace,'(a)') '@    s1 symbol char 65'
       write(iungrace,'(a)') '@    s1 symbol char font 0'
       write(iungrace,'(a)') '@    s1 symbol skip 0'
       write(iungrace,'(a)') '@    s1 line type 0'
       write(iungrace,'(a)') '@    s1 line linestyle 3'
       write(iungrace,'(a)') '@    s1 line linewidth 1.0'
       write(iungrace,'(a)') '@    s1 line color 15'
       write(iungrace,'(a)') '@    s1 line pattern 1'
       write(iungrace,'(a)') '@    s1 baseline type 0'
       write(iungrace,'(a)') '@    s1 baseline off'
       write(iungrace,'(a)') '@    s1 dropline off'
       write(iungrace,'(a)') '@    s1 fill type 0'
       write(iungrace,'(a)') '@    s1 fill rule 0'
       write(iungrace,'(a)') '@    s1 fill color 1'
       write(iungrace,'(a)') '@    s1 fill pattern 1'
       write(iungrace,'(a)') '@    s1 avalue off'
       write(iungrace,'(a)') '@    s1 avalue type 2'
       write(iungrace,'(a)') '@    s1 avalue char size 1.000000'
       write(iungrace,'(a)') '@    s1 avalue font 0'
       write(iungrace,'(a)') '@    s1 avalue color 1'
       write(iungrace,'(a)') '@    s1 avalue rot 0'
       write(iungrace,'(a)') '@    s1 avalue format general'
       write(iungrace,'(a)') '@    s1 avalue prec 3'
       write(iungrace,'(a)') '@    s1 avalue prepend ""'
       write(iungrace,'(a)') '@    s1 avalue append ""'
       write(iungrace,'(a)') 
     &    '@    s1 avalue offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    s1 errorbar on'
       write(iungrace,'(a)') '@    s1 errorbar place both'
       write(iungrace,'(a)') '@    s1 errorbar color 15'
       write(iungrace,'(a)') '@    s1 errorbar pattern 1'
       write(iungrace,'(a)') '@    s1 errorbar size 1.000000'
       write(iungrace,'(a)') '@    s1 errorbar linewidth 1.0'
       write(iungrace,'(a)') '@    s1 errorbar linestyle 1'
       write(iungrace,'(a)') '@    s1 errorbar riser linewidth 1.0'
       write(iungrace,'(a)') '@    s1 errorbar riser linestyle 1'
       write(iungrace,'(a)') '@    s1 errorbar riser clip off'
       write(iungrace,'(a)') 
     &   '@    s1 errorbar riser clip length 0.100000'
       write(iungrace,'(a)') '@    s1 comment "reference QM data"'
       write(iungrace,'(a)') '@    s1 legend  "QM\Scorr\N"'
       write(iungrace,'(a)') '@    s2 hidden false'
       write(iungrace,'(a)') '@    s2 type xy'
       write(iungrace,'(a)') '@    s2 symbol 0'
       write(iungrace,'(a)') '@    s2 symbol size 1.000000'
       write(iungrace,'(a)') '@    s2 symbol color 6'
       write(iungrace,'(a)') '@    s2 symbol pattern 1'
       write(iungrace,'(a)') '@    s2 symbol fill color 6'
       write(iungrace,'(a)') '@    s2 symbol fill pattern 0'
       write(iungrace,'(a)') '@    s2 symbol linewidth 1.0'
       write(iungrace,'(a)') '@    s2 symbol linestyle 1'
       write(iungrace,'(a)') '@    s2 symbol char 65'
       write(iungrace,'(a)') '@    s2 symbol char font 0'
       write(iungrace,'(a)') '@    s2 symbol skip 0'
       write(iungrace,'(a)') '@    s2 line type 1'
       write(iungrace,'(a)') '@    s2 line linestyle 4'
       write(iungrace,'(a)') '@    s2 line linewidth 2.0'
       write(iungrace,'(a)') '@    s2 line color 6'
       write(iungrace,'(a)') '@    s2 line pattern 1'
       write(iungrace,'(a)') '@    s2 baseline type 0'
       write(iungrace,'(a)') '@    s2 baseline off'
       write(iungrace,'(a)') '@    s2 dropline off'
       write(iungrace,'(a)') '@    s2 fill type 0'
       write(iungrace,'(a)') '@    s2 fill rule 0'
       write(iungrace,'(a)') '@    s2 fill color 1'
       write(iungrace,'(a)') '@    s2 fill pattern 1'
       write(iungrace,'(a)') '@    s2 avalue off'
       write(iungrace,'(a)') '@    s2 avalue type 2'
       write(iungrace,'(a)') '@    s2 avalue char size 1.000000'
       write(iungrace,'(a)') '@    s2 avalue font 0'
       write(iungrace,'(a)') '@    s2 avalue color 1'
       write(iungrace,'(a)') '@    s2 avalue rot 0'
       write(iungrace,'(a)') '@    s2 avalue format general'
       write(iungrace,'(a)') '@    s2 avalue prec 3'
       write(iungrace,'(a)') '@    s2 avalue prepend ""'
       write(iungrace,'(a)') '@    s2 avalue append ""'
       write(iungrace,'(a)') '@    s2 avalue offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    s2 errorbar on'
       write(iungrace,'(a)') '@    s2 errorbar place both'
       write(iungrace,'(a)') '@    s2 errorbar color 6'
       write(iungrace,'(a)') '@    s2 errorbar pattern 1'
       write(iungrace,'(a)') '@    s2 errorbar size 1.000000'
       write(iungrace,'(a)') '@    s2 errorbar linewidth 1.0'
       write(iungrace,'(a)') '@    s2 errorbar linestyle 1'
       write(iungrace,'(a)') '@    s2 errorbar riser linewidth 1.0'
       write(iungrace,'(a)') '@    s2 errorbar riser linestyle 1'
       write(iungrace,'(a)') '@    s2 errorbar riser clip off'
       write(iungrace,'(a)') 
     &   '@    s2 errorbar riser clip length 0.100000'
       write(iungrace,'(a)') '@    s2 comment "QMD-FF data"'
       write(iungrace,'(a)') '@    s2 legend  "E\sft\N (\xd\f{})"'
       write(iungrace,'(a)') '@target G0.S0'
       write(iungrace,'(a)') '@type xy'
       do i = 1, Nqm
         write(iungrace,'(f9.1,5x,f12.3)') delta(i), E(i)
       enddo
       write(iungrace,'(a)') '&'
       write(iungrace,'(a)') '@target G0.S1'
       write(iungrace,'(a)') '@type xy'
       do i = 1, Nqm
         write(iungrace,'(f9.1,5x,f12.3)') delta(i), Ecorr(i)
       enddo
       write(iungrace,'(a)') '&'
       write(iungrace,'(a)') '@target G0.S2'
       write(iungrace,'(a)') '@type xy'

      else ! printing grace file for torsion  without LJ assign

       write(iungrace,'(a)') '# Grace project file'
       write(iungrace,'(a)') '#'
       write(iungrace,'(a)') '@version 50125'
       write(iungrace,'(a)') '@page size 792, 612'
       write(iungrace,'(a)') '@page scroll 5%'
       write(iungrace,'(a)') '@page inout 5%'
       write(iungrace,'(a)') '@link page off'
       write(iungrace,'(a)') 
     &  '@map font 0 to "Times-Roman", "Times-Roman"'
       write(iungrace,'(a)') 
     &  '@map font 1 to "Times-Italic", "Times-Italic"'
       write(iungrace,'(a)') 
     &  '@map font 2 to "Times-Bold", "Times-Bold"'
       write(iungrace,'(a)') 
     &    '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"'
       write(iungrace,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
       write(iungrace,'(a)') 
     &    '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"'
       write(iungrace,'(a)') 
     &    '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"'
       write(iungrace,'(a)') 
     & '@map font 7 to "Helvetica-BoldOblique","Helvetica-BoldOblique"'
       write(iungrace,'(a)') '@map font 8 to "Courier", "Courier"'
       write(iungrace,'(a)') 
     &   '@map font 9 to "Courier-Oblique", "Courier-Oblique"'
       write(iungrace,'(a)') 
     &   '@map font 10 to "Courier-Bold", "Courier-Bold"'
       write(iungrace,'(a)') 
     &   '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"'
       write(iungrace,'(a)') '@map font 12 to "Symbol", "Symbol"'
       write(iungrace,'(a)') 
     &    '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
       write(iungrace,'(a)') '@map color 0 to (255, 255, 255), "white"'
       write(iungrace,'(a)') '@map color 1 to (0, 0, 0), "black"'
       write(iungrace,'(a)') '@map color 2 to (255, 0, 0), "red"'
       write(iungrace,'(a)') '@map color 3 to (0, 255, 0), "green"'
       write(iungrace,'(a)') '@map color 4 to (0, 0, 255), "blue"'
       write(iungrace,'(a)') '@map color 5 to (255, 255, 0), "yellow"'
       write(iungrace,'(a)') '@map color 6 to (188, 143, 143), "brown"'
       write(iungrace,'(a)') '@map color 7 to (220, 220, 220), "grey"'
       write(iungrace,'(a)') '@map color 8 to (148, 0, 211), "violet"'
       write(iungrace,'(a)') '@map color 9 to (0, 255, 255), "cyan"'
       write(iungrace,'(a)') 
     &  '@map color 10 to (255, 0, 255), "magenta"'
       write(iungrace,'(a)') '@map color 11 to (255, 165, 0), "orange"'
       write(iungrace,'(a)') 
     &  '@map color 12 to (114, 33, 188), "indigo"'
       write(iungrace,'(a)') '@map color 13 to (103, 7, 72), "maroon"'
       write(iungrace,'(a)') 
     &  '@map color 14 to (64, 224, 208), "turquoise"'
       write(iungrace,'(a)') '@map color 15 to (0, 139, 0), "green4"'
       write(iungrace,'(a)') '@reference date 0'
       write(iungrace,'(a)') '@date wrap off'
       write(iungrace,'(a)') '@date wrap year 1950'
       write(iungrace,'(a)') '@default linewidth 1.0'
       write(iungrace,'(a)') '@default linestyle 1'
       write(iungrace,'(a)') '@default color 1'
       write(iungrace,'(a)') '@default pattern 1'
       write(iungrace,'(a)') '@default font 0'
       write(iungrace,'(a)') '@default char size 1.000000'
       write(iungrace,'(a)') '@default symbol size 1.000000'
       write(iungrace,'(a)') '@default sformat "%.8g"'
       write(iungrace,'(a)') '@background color 0'
       write(iungrace,'(a)') '@page background fill on'
       write(iungrace,'(a)') '@timestamp off'
       write(iungrace,'(a)') '@timestamp 0.03, 0.03'
       write(iungrace,'(a)') '@timestamp color 1'
       write(iungrace,'(a)') '@timestamp rot 0'
       write(iungrace,'(a)') '@timestamp font 0'
       write(iungrace,'(a)') '@timestamp char size 1.000000'
       write(iungrace,'(a)') '@timestamp def "      "'
       write(iungrace,'(a)') '@r0 off'
       write(iungrace,'(a)') '@link r0 to g0'
       write(iungrace,'(a)') '@r0 type above'
       write(iungrace,'(a)') '@r0 linestyle 1'
       write(iungrace,'(a)') '@r0 linewidth 1.0'
       write(iungrace,'(a)') '@r0 color 1'
       write(iungrace,'(a)') '@r0 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r1 off'
       write(iungrace,'(a)') '@link r1 to g0'
       write(iungrace,'(a)') '@r1 type above'
       write(iungrace,'(a)') '@r1 linestyle 1'
       write(iungrace,'(a)') '@r1 linewidth 1.0'
       write(iungrace,'(a)') '@r1 color 1'
       write(iungrace,'(a)') '@r1 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r2 off'
       write(iungrace,'(a)') '@link r2 to g0'
       write(iungrace,'(a)') '@r2 type above'
       write(iungrace,'(a)') '@r2 linestyle 1'
       write(iungrace,'(a)') '@r2 linewidth 1.0'
       write(iungrace,'(a)') '@r2 color 1'
       write(iungrace,'(a)') '@r2 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r3 off'
       write(iungrace,'(a)') '@link r3 to g0'
       write(iungrace,'(a)') '@r3 type above'
       write(iungrace,'(a)') '@r3 linestyle 1'
       write(iungrace,'(a)') '@r3 linewidth 1.0'
       write(iungrace,'(a)') '@r3 color 1'
       write(iungrace,'(a)') '@r3 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@r4 off'
       write(iungrace,'(a)') '@link r4 to g0'
       write(iungrace,'(a)') '@r4 type above'
       write(iungrace,'(a)') '@r4 linestyle 1'
       write(iungrace,'(a)') '@r4 linewidth 1.0'
       write(iungrace,'(a)') '@r4 color 1'
       write(iungrace,'(a)') '@r4 line 0, 0, 0, 0'
       write(iungrace,'(a)') '@g0 on'
       write(iungrace,'(a)') '@g0 hidden false'
       write(iungrace,'(a)') '@g0 type XY'
       write(iungrace,'(a)') '@g0 stacked false'
       write(iungrace,'(a)') '@g0 bar hgap 0.000000'
       write(iungrace,'(a)') '@g0 fixedpoint off'
       write(iungrace,'(a)') '@g0 fixedpoint type 0'
       write(iungrace,'(a)') '@g0 fixedpoint xy 0.000000, 0.000000'
       write(iungrace,'(a)') '@g0 fixedpoint format general general'
       write(iungrace,'(a)') '@g0 fixedpoint prec 6, 6'
       write(iungrace,'(a)') '@with g0'
       write(iungrace,'(a24,1x,f12.0)') '@    world -180, 0, 180,',emax
       write(iungrace,'(a)') '@    stack world 0, 0, 0, 0'
       write(iungrace,'(a)') '@    znorm 1'
       write(iungrace,'(a)') 
     &   '@    view 0.150000, 0.150000, 1.150000, 0.850000'
       write(iungrace,'(a)') '@    title ""'
       write(iungrace,'(a)') '@    title font 0'
       write(iungrace,'(a)') '@    title size 1.500000'
       write(iungrace,'(a)') '@    title color 1'
       write(iungrace,'(a)') '@    subtitle ""'
       write(iungrace,'(a)') '@    subtitle font 0'
       write(iungrace,'(a)') '@    subtitle size 1.000000'
       write(iungrace,'(a)') '@    subtitle color 1'
       write(iungrace,'(a)') '@    xaxes scale Normal'
       write(iungrace,'(a)') '@    yaxes scale Normal'
       write(iungrace,'(a)') '@    xaxes invert off'
       write(iungrace,'(a)') '@    yaxes invert off'
       write(iungrace,'(a)') '@    xaxis  on'
       write(iungrace,'(a)') '@    xaxis  type zero false'
       write(iungrace,'(a)') '@    xaxis  offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    xaxis  bar on'
       write(iungrace,'(a)') '@    xaxis  bar color 1'
       write(iungrace,'(a)') '@    xaxis  bar linestyle 1'
       write(iungrace,'(a)') '@    xaxis  bar linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  label "\xd\f{} (degr)"'
       write(iungrace,'(a)') '@    xaxis  label layout para'
       write(iungrace,'(a)') '@    xaxis  label place auto'
       write(iungrace,'(a)') '@    xaxis  label char size 1.000000'
       write(iungrace,'(a)') '@    xaxis  label font 4'
       write(iungrace,'(a)') '@    xaxis  label color 1'
       write(iungrace,'(a)') '@    xaxis  label place normal'
       write(iungrace,'(a)') '@    xaxis  tick on'
       write(iungrace,'(a)') '@    xaxis  tick major 30'
       write(iungrace,'(a)') '@    xaxis  tick minor ticks 1'
       write(iungrace,'(a)') '@    xaxis  tick default 6'
       write(iungrace,'(a)') '@    xaxis  tick place rounded true'
       write(iungrace,'(a)') '@    xaxis  tick in'
       write(iungrace,'(a)') '@    xaxis  tick major size 1.000000'
       write(iungrace,'(a)') '@    xaxis  tick major color 1'
       write(iungrace,'(a)') '@    xaxis  tick major linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  tick major linestyle 1'
       write(iungrace,'(a)') '@    xaxis  tick major grid off'
       write(iungrace,'(a)') '@    xaxis  tick minor color 1'
       write(iungrace,'(a)') '@    xaxis  tick minor linewidth 1.0'
       write(iungrace,'(a)') '@    xaxis  tick minor linestyle 1'
       write(iungrace,'(a)') '@    xaxis  tick minor grid off'
       write(iungrace,'(a)') '@    xaxis  tick minor size 0.500000'
       write(iungrace,'(a)') '@    xaxis  ticklabel on'
       write(iungrace,'(a)') '@    xaxis  ticklabel format general'
       write(iungrace,'(a)') '@    xaxis  ticklabel prec 5'
       write(iungrace,'(a)') '@    xaxis  ticklabel formula ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel append ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel prepend ""'
       write(iungrace,'(a)') '@    xaxis  ticklabel angle 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel skip 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel stagger 0'
       write(iungrace,'(a)') '@    xaxis  ticklabel place normal'
       write(iungrace,'(a)') '@    xaxis  ticklabel offset auto'
       write(iungrace,'(a)') 
     &    '@    xaxis  ticklabel offset 0.000000 , 0.010000'
       write(iungrace,'(a)') '@    xaxis  ticklabel start type auto'
       write(iungrace,'(a)') '@    xaxis  ticklabel start 0.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel stop type auto'
       write(iungrace,'(a)') '@    xaxis  ticklabel stop 0.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel char size 1.000000'
       write(iungrace,'(a)') '@    xaxis  ticklabel font 4'
       write(iungrace,'(a)') '@    xaxis  ticklabel color 1'
       write(iungrace,'(a)') '@    xaxis  tick place both'
       write(iungrace,'(a)') '@    xaxis  tick spec type none'
       write(iungrace,'(a)') '@    yaxis  on'
       write(iungrace,'(a)') '@    yaxis  type zero false'
       write(iungrace,'(a)') '@    yaxis  offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    yaxis  bar on'
       write(iungrace,'(a)') '@    yaxis  bar color 1'
       write(iungrace,'(a)') '@    yaxis  bar linestyle 1'
       write(iungrace,'(a)') '@    yaxis  bar linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  label "\xD\f{}E (kJ/mol)"'
       write(iungrace,'(a)') '@    yaxis  label layout para'
       write(iungrace,'(a)') '@    yaxis  label place auto'
       write(iungrace,'(a)') '@    yaxis  label char size 1.000000'
       write(iungrace,'(a)') '@    yaxis  label font 4'
       write(iungrace,'(a)') '@    yaxis  label color 1'
       write(iungrace,'(a)') '@    yaxis  label place normal'
       write(iungrace,'(a)') '@    yaxis  tick on'
       write(iungrace,'(a)') '@    yaxis  tick major 5'
       write(iungrace,'(a)') '@    yaxis  tick minor ticks 1'
       write(iungrace,'(a)') '@    yaxis  tick default 6'
       write(iungrace,'(a)') '@    yaxis  tick place rounded true'
       write(iungrace,'(a)') '@    yaxis  tick in'
       write(iungrace,'(a)') '@    yaxis  tick major size 1.000000'
       write(iungrace,'(a)') '@    yaxis  tick major color 1'
       write(iungrace,'(a)') '@    yaxis  tick major linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  tick major linestyle 1'
       write(iungrace,'(a)') '@    yaxis  tick major grid off'
       write(iungrace,'(a)') '@    yaxis  tick minor color 1'
       write(iungrace,'(a)') '@    yaxis  tick minor linewidth 1.0'
       write(iungrace,'(a)') '@    yaxis  tick minor linestyle 1'
       write(iungrace,'(a)') '@    yaxis  tick minor grid off'
       write(iungrace,'(a)') '@    yaxis  tick minor size 0.500000'
       write(iungrace,'(a)') '@    yaxis  ticklabel on'
       write(iungrace,'(a)') '@    yaxis  ticklabel format general'
       write(iungrace,'(a)') '@    yaxis  ticklabel prec 5'
       write(iungrace,'(a)') '@    yaxis  ticklabel formula ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel append ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel prepend ""'
       write(iungrace,'(a)') '@    yaxis  ticklabel angle 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel skip 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel stagger 0'
       write(iungrace,'(a)') '@    yaxis  ticklabel place normal'
       write(iungrace,'(a)') '@    yaxis  ticklabel offset auto'
       write(iungrace,'(a)') 
     &        '@    yaxis  ticklabel offset 0.000000 , 0.010000'
       write(iungrace,'(a)') '@    yaxis  ticklabel start type auto'
       write(iungrace,'(a)') '@    yaxis  ticklabel start 0.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel stop type auto'
       write(iungrace,'(a)') '@    yaxis  ticklabel stop 0.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel char size 1.000000'
       write(iungrace,'(a)') '@    yaxis  ticklabel font 4'
       write(iungrace,'(a)') '@    yaxis  ticklabel color 1'
       write(iungrace,'(a)') '@    yaxis  tick place both'
       write(iungrace,'(a)') '@    yaxis  tick spec type none'
       write(iungrace,'(a)') '@    altxaxis  off'
       write(iungrace,'(a)') '@    altyaxis  off'
       write(iungrace,'(a)') '@    legend on'
       write(iungrace,'(a)') '@    legend loctype view'
       write(iungrace,'(a)') '@    legend 0.85, 0.8'
       write(iungrace,'(a)') '@    legend box color 1'
       write(iungrace,'(a)') '@    legend box pattern 1'
       write(iungrace,'(a)') '@    legend box linewidth 1.0'
       write(iungrace,'(a)') '@    legend box linestyle 1'
       write(iungrace,'(a)') '@    legend box fill color 0'
       write(iungrace,'(a)') '@    legend box fill pattern 1'
       write(iungrace,'(a)') '@    legend font 0'
       write(iungrace,'(a)') '@    legend char size 1.000000'
       write(iungrace,'(a)') '@    legend color 1'
       write(iungrace,'(a)') '@    legend length 4'
       write(iungrace,'(a)') '@    legend vgap 1'
       write(iungrace,'(a)') '@    legend hgap 1'
       write(iungrace,'(a)') '@    legend invert false'
       write(iungrace,'(a)') '@    frame type 0'
       write(iungrace,'(a)') '@    frame linestyle 1'
       write(iungrace,'(a)') '@    frame linewidth 1.0'
       write(iungrace,'(a)') '@    frame color 1'
       write(iungrace,'(a)') '@    frame pattern 1'
       write(iungrace,'(a)') '@    frame background color 0'
       write(iungrace,'(a)') '@    frame background pattern 0'
       write(iungrace,'(a)') '@    s0 hidden false'
       write(iungrace,'(a)') '@    s0 type xy'
       write(iungrace,'(a)') '@    s0 symbol 1'
       write(iungrace,'(a)') '@    s0 symbol size 1.000000'
       write(iungrace,'(a)') '@    s0 symbol color 15'
       write(iungrace,'(a)') '@    s0 symbol pattern 1'
       write(iungrace,'(a)') '@    s0 symbol fill color 15'
       write(iungrace,'(a)') '@    s0 symbol fill pattern 1'
       write(iungrace,'(a)') '@    s0 symbol linewidth 1.0'
       write(iungrace,'(a)') '@    s0 symbol linestyle 1'
       write(iungrace,'(a)') '@    s0 symbol char 65'
       write(iungrace,'(a)') '@    s0 symbol char font 0'
       write(iungrace,'(a)') '@    s0 symbol skip 0'
       write(iungrace,'(a)') '@    s0 line type 0'
       write(iungrace,'(a)') '@    s0 line linestyle 3'
       write(iungrace,'(a)') '@    s0 line linewidth 1.0'
       write(iungrace,'(a)') '@    s0 line color 15'
       write(iungrace,'(a)') '@    s0 line pattern 1'
       write(iungrace,'(a)') '@    s0 baseline type 0'
       write(iungrace,'(a)') '@    s0 baseline off'
       write(iungrace,'(a)') '@    s0 dropline off'
       write(iungrace,'(a)') '@    s0 fill type 0'
       write(iungrace,'(a)') '@    s0 fill rule 0'
       write(iungrace,'(a)') '@    s0 fill color 1'
       write(iungrace,'(a)') '@    s0 fill pattern 1'
       write(iungrace,'(a)') '@    s0 avalue off'
       write(iungrace,'(a)') '@    s0 avalue type 2'
       write(iungrace,'(a)') '@    s0 avalue char size 1.000000'
       write(iungrace,'(a)') '@    s0 avalue font 0'
       write(iungrace,'(a)') '@    s0 avalue color 1'
       write(iungrace,'(a)') '@    s0 avalue rot 0'
       write(iungrace,'(a)') '@    s0 avalue format general'
       write(iungrace,'(a)') '@    s0 avalue prec 3'
       write(iungrace,'(a)') '@    s0 avalue prepend ""'
       write(iungrace,'(a)') '@    s0 avalue append ""'
       write(iungrace,'(a)') 
     &   '@    s0 avalue offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    s0 errorbar on'
       write(iungrace,'(a)') '@    s0 errorbar place both'
       write(iungrace,'(a)') '@    s0 errorbar color 15'
       write(iungrace,'(a)') '@    s0 errorbar pattern 1'
       write(iungrace,'(a)') '@    s0 errorbar size 1.000000'
       write(iungrace,'(a)') '@    s0 errorbar linewidth 1.0'
       write(iungrace,'(a)') '@    s0 errorbar linestyle 1'
       write(iungrace,'(a)') '@    s0 errorbar riser linewidth 1.0'
       write(iungrace,'(a)') '@    s0 errorbar riser linestyle 1'
       write(iungrace,'(a)') '@    s0 errorbar riser clip off'
       write(iungrace,'(a)') 
     &    '@    s0 errorbar riser clip length 0.100000'
       write(iungrace,'(a)') '@    s0 comment "reference QM data"'
       write(iungrace,'(a)') '@    s0 legend  "QM"'
       write(iungrace,'(a)') '@    s1 hidden false'
       write(iungrace,'(a)') '@    s1 type xy'
       write(iungrace,'(a)') '@    s1 symbol 0'
       write(iungrace,'(a)') '@    s1 symbol size 1.000000'
       write(iungrace,'(a)') '@    s1 symbol color 4'
       write(iungrace,'(a)') '@    s1 symbol pattern 1'
       write(iungrace,'(a)') '@    s1 symbol fill color 4'
       write(iungrace,'(a)') '@    s1 symbol fill pattern 0'
       write(iungrace,'(a)') '@    s1 symbol linewidth 1.0'
       write(iungrace,'(a)') '@    s1 symbol linestyle 1'
       write(iungrace,'(a)') '@    s1 symbol char 65'
       write(iungrace,'(a)') '@    s1 symbol char font 0'
       write(iungrace,'(a)') '@    s1 symbol skip 0'
       write(iungrace,'(a)') '@    s1 line type 1'
       write(iungrace,'(a)') '@    s1 line linestyle 4'
       write(iungrace,'(a)') '@    s1 line linewidth 2.0'
       write(iungrace,'(a)') '@    s1 line color 4'
       write(iungrace,'(a)') '@    s1 line pattern 1'
       write(iungrace,'(a)') '@    s1 baseline type 0'
       write(iungrace,'(a)') '@    s1 baseline off'
       write(iungrace,'(a)') '@    s1 dropline off'
       write(iungrace,'(a)') '@    s1 fill type 0'
       write(iungrace,'(a)') '@    s1 fill rule 0'
       write(iungrace,'(a)') '@    s1 fill color 1'
       write(iungrace,'(a)') '@    s1 fill pattern 1'
       write(iungrace,'(a)') '@    s1 avalue off'
       write(iungrace,'(a)') '@    s1 avalue type 2'
       write(iungrace,'(a)') '@    s1 avalue char size 1.000000'
       write(iungrace,'(a)') '@    s1 avalue font 0'
       write(iungrace,'(a)') '@    s1 avalue color 1'
       write(iungrace,'(a)') '@    s1 avalue rot 0'
       write(iungrace,'(a)') '@    s1 avalue format general'
       write(iungrace,'(a)') '@    s1 avalue prec 3'
       write(iungrace,'(a)') '@    s1 avalue prepend ""'
       write(iungrace,'(a)') '@    s1 avalue append ""'
       write(iungrace,'(a)') 
     &  '@    s1 avalue offset 0.000000 , 0.000000'
       write(iungrace,'(a)') '@    s1 errorbar on'
       write(iungrace,'(a)') '@    s1 errorbar place both'
       write(iungrace,'(a)') '@    s1 errorbar color 4'
       write(iungrace,'(a)') '@    s1 errorbar pattern 1'
       write(iungrace,'(a)') '@    s1 errorbar size 1.000000'
       write(iungrace,'(a)') '@    s1 errorbar linewidth 1.0'
       write(iungrace,'(a)') '@    s1 errorbar linestyle 1'
       write(iungrace,'(a)') '@    s1 errorbar riser linewidth 1.0'
       write(iungrace,'(a)') '@    s1 errorbar riser linestyle 1'
       write(iungrace,'(a)') '@    s1 errorbar riser clip off'
       write(iungrace,'(a)') 
     &        '@    s1 errorbar riser clip length 0.100000'
       write(iungrace,'(a)') '@    s1 comment "QMD-FF data"'
       write(iungrace,'(a)') '@    s1 legend  "QMD-FF (FIRA)"'
       write(iungrace,'(a)') '@target G0.S0'
       write(iungrace,'(a)') '@type xy'
       do i = 1, Nqm
          write(iungrace,'(f9.1,5x,f12.3)') delta(i), E(i)
       enddo
       write(iungrace,'(a)') '&'
       write(iungrace,'(a)') '@target G0.S1'
       write(iungrace,'(a)') '@type xy'
      endif

      return
      end 
