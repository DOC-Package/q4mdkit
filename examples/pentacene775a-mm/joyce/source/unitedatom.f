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
      Subroutine United_Atom (geom,nat,nator,aname,kbond
     &                       ,UAM,manua,zstop,maxat,max3n)
c----------------------------------------------------------------------c
 
 
c     manages the United Atom input and set the transformation matrix
c     for the cartesian gradients G and the cartesian Hessian matrix H.
 
 
c                INPUT EXAMPLE:
c     $UnitedAtom
c       6 - 9 12         # put together atoms 6 7 8 9 12 in a single UA
c       stop             # tell the program to stop, after reading all UA data
c       22 15 16         # put together atoms 15 16 22 in a single UA
c     $end
 
 
c     geom    (i)    current original geometry (1..nator)
c             (o)    UA geometry (1..nat)
c     nat     (io)   # of atoms in UA geometry
c     nator   (o)    # of original atoms
c     aname   (io)   name of the atoms (changed)
c     kbond   (i)    bond order matrix in the original geometry (changed ??)
c     UAM     (o)    transformation matrix for gradient and hessian
c     manua   (o)    array to manage the UA rearrangements
c     zstop   (o)    T = in the UA data a stop has been included
c     maxat,max3n  parameter value in the caller
 
 
      implicit real*8(a-h,o-z)
      parameter (max=300)
      dimension geom(maxat,*),kbond(maxat,*),UAM(max3n,*)
     &         ,manua(maxat,*)
      dimension achar(max),isel(max),iaway(max)
      character aname(*)*(*),card*90,str5*5
      logical zbari,zrew,zstop,zord
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      common /printout/ levprt
 
      Call Check_Dim (nat,max,'UNITED_ATOM: increase max')
      nator=nat
      ncart=nator*3
c
c   save the original geometry in geor
c
      do i=1,nator
        iaway(i)=0
        achar(i)=geom(i,4)
        do k=1,8
          manua(i,k)=0
        enddo
        manua(i,1)=i
        manua(i,8)=1
        manua(i,9)=0
      enddo
 
      iun=5
      zstop=.false.
      zrew=.true.
      Call Cerca_W ('$United',iun,card,ier,zrew)
      nnua=0
      if(ier.ne.0) return
c     choose where to set UA center
      zbari=.false.
      ii=index(card,'-b')
      if(ii.gt.0) zbari=.true.
      ii=index(card,'-h')
      if(ii.gt.0) zbari=.false.
      do i=1,nator
        manua(i,9)=-1            !! UA on the heavy atom
        if(zbari) manua(i,9)=+1  !! UA in the center of mass of the group
      enddo
 
      write(6,8)
      if(.not.zbari) write(6,10) 'on the heavy atom'
      if(zbari)      write(6,10) 'in the center of mass of the group'
    8 format(//16x,'=============================================='
     &       ,/16x,'                UNITED ATOM'
     &       ,/16x,' a group of atoms is treated as a single atom'
     &       ,/16x,'     the forces are assembled in the UA'
     &       ,/16x,'==============================================')
   10 format(   2x,'  POSITION of the  UA: ',a)
 
c     input begins
 
      write(6,52)
   20 read(5,'(a)',end=50) card 
      Call LowCas (card)
      if(index(card(1:9),'$end').gt.0) goto 50
      if(index(card(1:4),'#').gt.0)    goto 20
      if(index(card(1:4),'%').gt.0)    goto 20
      Call SetGrat (card)
      if(index(card,'stop').gt.0) then
        zstop=.true.
        goto 20
      endif
      nnua=nnua+1
 
      ibeg=1
      Call Read_Atoms (nator,card,ibeg,isel,nsel)
      if(nsel.eq.1) goto 100
      if(nsel.le.0) goto 20
 
c   check: only one atom can have bonds with not involved atoms
c   set: iua
 
      iua=0
      do 30 is=1,nsel
        ia=isel(is)
        do 25 ka=1,nator
          do js=1,nsel
            if(isel(js).eq.ka) goto 25
          enddo
        if(kbond(ia,ka).gt.0) then
          if(iua.gt.0 .and. ia.ne.iua) goto 105
          iua=ia
        endif
   25 continue
   30 continue
      if(iua.eq.0) goto 110
 
c     check: all the involved atoms must be bonded to atom iua 
 
      do is=1,nsel
        ia=isel(is)
        if(ia.ne.iua .and. kbond(ia,iua).le.0) goto 120
      enddo
 
c     set the /manua/ array which drives the UA transformation
 
      qua=0.0d0
      aua=0.0d0
      js=1
      do is=1,nsel
        ia=isel(is)
        qua=qua+achar(ia)
        aua=aua+geom(ia,5)
        if(iaway(ia).gt.0) goto 125
        if(ia.eq.iua) then
          iaway(ia)=2
        else
          iaway(ia)=1
          js=js+1
          manua(iua,js)=ia
        endif
        manua(ia,1)=0
      enddo
      manua(iua,1)=iua
      manua(iua,8)=nsel
      achar(iua)=qua
      ll=Lenw(aname(iua))
      str5=aname(iua)
      str5(ll+1:ll+1)='*'
      write(6,54) iua,str5,qua,aua,(isel(k),aname(isel(k)),k=1,nsel)
      aname(iua)(ll+1:ll+1)='*'
 
      goto 20
   50 continue
   52 format(/3x,'United Atom   charge  mass    -----------'
     &      ,' atoms included ------------')
   54 format(i6,2x,a,2f8.1,t33,4(i4,'-',a5)/(t33,4(i4,'-',a5)))
 
c     end input
 
 
c     Update the arrays /manua/ and /aname/, discarding the unused atoms
c     WARNING: the cardinal number of the nuclei will be changed
 
      ja=0
      do ia=1,nator
        if(iaway(ia).ne.1) then
          ja=ja+1
          aname(ja)=aname(ia)
          do k=1,8
            manua(ja,k)=manua(ia,k)
          enddo
        endif
      enddo
      nat=ja
      write(6,66) (i,aname(i),(manua(i,k),k=1,9),i=1,nat)
   66 format(/9x,'Driving array for the UA prescriptions'
     &            ,/(i8,2x,a5,'  -',9i4))
 
c   set the transformation matrix UAM for gradient and hessian
c   in cartesian coordinates
 
      ncartua=nat*3
      do i=1,ncart
        do j=1,ncartua
          UAM(i,j)=0.0d0
        enddo
      enddo
      do ia=1,nat
        k=1
        dowhile(manua(ia,k).gt.0)
          ja=manua(ia,k)
          do ix=1,3
            iax=(ia-1)*3+ix
            jax=(ja-1)*3+ix
            UAM(jax,iax)=1.0d0
          enddo
          k=k+1
        enddo
      enddo
 
c     transform the geometry according to the UA prescriptions
c     contained in the /manua/ array
 
      Call Transf_GE_UA (geom,manua,nat,nator,aname,maxat)
 
c     modify the aname array according bto the UA prescriptions
 
      zord=.true.
      Call Atom_Count (nat,geom,aname,zord,maxat)
 
c     punch data to be used for moscito.system
 
      open (unit=12,file='UnAtom.geo',form='formatted')
      rewind 12
      write(12,64) (aname(i),(boh2ang*geom(i,k),k=1,3)
     &                     ,geom(i,4),geom(i,5),i=1,nat)
   64 format(2x,a5,3f12.7,3x,2f8.1)
 
      ncol=10
      ndim=max3n
      i=nat*3
      if(levprt.ge.2)
     &Call Pr_Rec2 (UAM,ncart,i,'UA-matrix to transf. G,H',ncol,ndim)
      return
 
c     errors were found in the United Atom input
 
  100 write(6,101) card(1:Lenw(card))
  101 format(//8x,'ERROR in the UA input: only one atom'
     &    ,' has been selected',/4x,a,//)
      stop
 
  105 write(6,106) iua,ia,(isel(k),k=1,nsel)
  106 format(//8x,'ERROR in the UA input: only one atom'
     &    ,' can have bonds with external atoms'
     &       ,/8x,'Atoms that form bonds with ext atoms',2i6 
     &       ,/8x,'Complete list of atoms in the UA',12i4,/)
      write(6,107) iua,(kbond(iua,k),k=1,nat)
      write(6,107) ia ,(kbond(ia ,k),k=1,nat)
  107 format(4x,'Bonds for atom',i5,/(6(1x,10i1)))
      stop
 
  110 write(6,112) (isel(k),k=1,nsel)
  112 format(//8x,'ERROR in the UA input: no atom'
     &    ,' is bonded with external atoms'
     &       ,/8x,'Complete list of atoms in the UA',12i4,/)
      stop
 
  120 write(6,122) iua,ia,(kbond(iua,k),k=1,nat)
  122 format(//8x,'ERROR in the UA input: one involved'
     &                  ,' atom is not bonded with the United Atom'
     &       ,/8x,'UA =',i4,5x,'Involved atom',i5
     &       ,/8x,'bond array for iua',/(5(2x,10i1)),/)
      stop
 
  125 write(6,126) ia
  126 format(//8x,'ERROR in the UA input: one'
     &    ,' atom is included in two or more UAs',i5,//)
      stop
 
      end



c----------------------------------------------------------------------c
      Subroutine Transf_GE_UA (geom,manua,nat,nator,aname,maxat)
c----------------------------------------------------------------------c
 
 
c     transforms the geometry according to the united atom settings
 
c     geom   (i)   original geometry (1...nator)
c            (o)   UA geometry (1...nat)
c     manua  (i)   array to manage the UA rearrangements
c     nat    (i)   # of atoms after the UA rearrangement
c     nator  (i)   # of original atoms
c     maxat  (i)   parameter value in the caller for /geom/
 
      implicit real*8(a-h,o-z)
      parameter (max=200)
      dimension geom(maxat,*),manua(maxat,*),geomor(max,5)
     &         ,amom(3),vmom(9),baric(3),amomu(3),vmomu(9),baricu(3)
      character*5 aname(*)
      logical zbari
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      common /printout/ levprt
 
      if(nat.eq.nator) return
      Call Check_Dim (nator,max,'Transf_GE_UA: increase max')
 
      do i=1,nator
        do k=1,5
          geomor(i,k)=geom(i,k)
        enddo
      enddo
      zbari=(manua(1,9).gt.0)
      Call Flush (6)
 
      do i=1,nat
        na=manua(i,8)
        do k=1,na        !! check
          ia=manua(i,k)
          if(ia.gt.nator .or. ia.le.0) goto 40
        enddo
        ia=manua(i,1)
        do k=1,3        !! position of UA in the main atom
          geom(i,k)=geomor(ia,k)
        enddo
        geom(i,4)=geom(ia,4)
        totw=0.0d0
        do k=1,na
          totw=totw+geomor(manua(i,k),5)
        enddo
        geom(i,5)=totw

        if(zbari) then  !! position of UA in the center of mass
          do k=1,3
            baric(k)=0.d0
          enddo
          totw=0.0d0
          do k=1,na
            ia=manua(i,k)
            wi=geomor(ia,5)
            totw=totw+wi
            do jx=1,3
              baric(jx)=baric(jx)+wi*geomor(ia,jx)
            enddo
          enddo
          do k=1,3
            geom(i,k)=baric(k)/totw
          enddo
        endif
      enddo
 
c     prints the new geometry
 
      if(levprt.ge.0) then
        write(6,22) (i,aname(i),(boh2ang*geom(i,k),k=1,3) 
     &                         ,(geom(i,k),k=4,5), i=1,nat)
        if(.not.zbari) write(6,24) 'on the heavy atom'
        if(zbari)      write(6,24) 'in the center of mass of the group'
      else
        write(6,26)
      endif
   22 format(/13x,'---- GEOMETRY DATA (ang) with UAs included ----'
     &      ,/18x,'x',11x,'y',11x,'z',6x,'UA-charge   mass'
     &      ,/(1x,i3,2x,a5,3f12.7,f8.1,f10.3))
   24 format('  POSITION of the  UA: ',a)
   26 format(/8x,'United Atom geometrical changes included')
 
c     computes the moment of inertia tensor 
 
      if(levprt.gt.1) then
        Call Inertia (geomor,nator,baric,amom,vmom,max)
        Call Inertia (geom,nat,baricu,amomu,vmomu,maxat)
        write(6,32) (amom(i),i=1,3),(amomu(i),i=1,3)
        do i=1,3
          write(6,34) (vmom(j),j=i,9,3),(vmomu(j),j=i,9,3)
        enddo
        write(6,36) (baric(i),i=1,3),(baricu(i),i=1,3)
      endif
   32 format(/6x,22('x'),' principal inertia axes (au) ',22('x')
     &  ,/17x,'Full Atomic',31x,'United Atom',/(2x,3f12.5,6x,3f12.5))
   34 format(1x,3f12.5,6x,3f12.5)
   36 format(35x,'center of mass',/1x,3f12.5,6x,3f12.5)
 
      return
 
   40 write(6,42) ia
   42 format(//8x,'ERROR in Transf_GE_UA: atom <1 or >nat'
     &       ,/8x,'                     : atom number',i8,//)
      stop
      end
 

c----------------------------------------------------------------------c
      Subroutine Transf_HG_UA (G,H,UAM,D2,nat,nator,max3n)
c----------------------------------------------------------------------c
 
 
c     transform the Gradient and Hessian matrices (in Cartesian coordinates)
c     according to the united atom settings
 
c     iflag   driving array for the definition of internal coordinates
c             it is assumed that the coupling terms 
c             come after the diag. terms
c     G       cartesian Gradients (3*nator) --> (3*nat)
c     H       cartesian Hessian matrix (3*nator,3*nator) --> (3*nat,3*nat)
c     UAM     united Atom trasformation matrix (3*nat,3*nator)
c     D2      work matrix
c     nat     number of atoms (with UA included)
c     nator   original number of atoms (normal molecule)
c     max3n   dimensions of some arrays in the caller
 
 
      implicit real*8(a-h,o-z)
      parameter (max=200)
      dimension G(*),H(max3n,*),UAM(max3n,*),D2(max3n,*)
      logical zover,zdeb
      zdeb=.false.
 
      if(nat.eq.nator) return
      nnew=nat*3
      nold=nator*3
      write(6,10) nator,nat
   10 format(/4x,'Gradient and Hessian matrices are trasformed'
     &      ,' for UA',2i5)
 
      do k=1,nold
        D2(k,1)=G(k)
      enddo
      do i=1,nnew
        gi=0.0d0
        do k=1,nold
          gi=gi+UAM(k,i)*D2(k,1)
        enddo
        G(i)=gi
      enddo
      if(zdeb) write(6,14) (k,(D2(l,1),l=k*3-2,k*3),k=1,nator)
      if(zdeb) write(6,16) (k,(G (l  ),l=k*3-2,k*3),k=1,nat)
   14 format(/' full-atom gradient vector (au)',/(i4,2x,3f10.5))
   16 format(/' unit-atom gradient vector (au)',/(i4,2x,3f10.5))
 
c     compute some data to doublecheck
 
      aaa=0.0d0
      aas=0.0d0
      do i=1,nold
        do k=1,nold
          aaa=aaa+H(i,k)
          aas=aas+dabs(H(i,k))
        enddo
      enddo
c     in the case H=0 (not read) : return
      if((aas/dble(nold*nold)).lt.1.d-18) return
 
c     transforms Hessian
 
      zover=.true.
      Call Trsml2 (H,UAM,D2,nnew,nold,max3n,zover)
 
      bbb=0.0d0
      do i=1,nnew
        do k=1,nnew
          bbb=bbb+H(i,k)
        enddo
      enddo
      tol=dabs(aaa)*dlog(dabs(aas/aaa))*1.d-06
      if(dabs(aaa-bbb).gt.tol) then
        write(6,24) aaa,bbb,dabs(aaa-bbb),tol
        stop
      endif
   24 format(//9x,'TRANSF_HG_UA: something wrong for Hessian transf'
     &       ,/14x,'Initial sum(H) .......',1pd16.8
     &       ,/14x,'Final   sum(H) .......',1pd16.8
     &       ,/14x,'difference ...........',1pd16.8,' >',D16.8,//)
 
      return
      end

 
c----------------------------------------------------------------------c
      Subroutine Read_Atoms (nat,card,ibeg,isel,nsel)
c----------------------------------------------------------------------c
 
c     reads a sequence of atoms in the i - j format
c     for instance 8 - 12  ---> 8 9 10 11 12
 
      implicit real*8(a-h,o-z)
      dimension isel(*)
      character card*(*),word*9
      Call Count_W (card,nword)
 
      iat1=1
      nsel=0
      ii=ibeg-1
   10 ii=ii+1
      if(ii.gt.nword)       goto 40
      Call SelWrd (card,word,ii,lw,ipos)
      if(word(1:2).eq.'; ') goto 40
      if(word(1:2).eq.', ') goto 40
      if(word(1:2).eq.'. ') goto 40
      if(word(1:2).eq.'  ') goto 40
 
      if(word(1:1).eq.'-') then
        ii=ii+1
        if(ii.gt.nword) goto 40
        Call SelWrd (card,word,ii,lw,ipos)
        read(word,*) iat2
        do i=Min0(iat1,iat2),Max0(iat1,iat2)
          if(i.gt.0 .and. i.le.nat .and. i.ne.isel(nsel)) then
            nsel=nsel+1
            isel(nsel)=i
          endif
        enddo
      else
        read(word,*) iat1
        if(iat1.gt.0 .and. iat1.le.nat) then
          nsel=nsel+1
          isel(nsel)=iat1
        endif
      endif
      goto 10
   40 continue
 
      return
      end
