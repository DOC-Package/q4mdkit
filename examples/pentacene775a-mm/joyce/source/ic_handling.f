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
      Subroutine Generate_IC (iprt, iflagn,ifftype,labnat,ncint,
     &                        nat,kbond,geom,aname,maxat,maxpar,cr,
     &                        sigma,epsilon,charge,
     &                        z_lj,zharmdihed,zstardihed)
c----------------------------------------------------------------------c

c        generates the internal 'natural' coordinates.

c  Warning: the natural int coord defined by Pulay are different.
c           as they consider symmetry properties


      implicit real*8(a-h,o-z)
      integer   cr
      dimension iflagn(maxpar,0:8),geom(maxat,*),kbond(maxat,*)
     &         ,valic(999),iatom(8)
      dimension sigma(*),charge(*),epsilon(*)
      character aname(*)*5,labnat(*)*21,label*23,string*66,st2*2
      character conv*4
      logical   zstardihed(maxpar),z14(maxat,maxat),zexcl(maxat,maxat)
      logical   z_lj,zli,zharmdihed(maxpar)

      IF (iprt .eq. 0) THEN
  
c        zero

         do i=1,maxpar
            do k=1,8
               iflagn(i,k) = 0
            enddo
            zstardihed(i) = .false.
            zharmdihed(i) = .false.
         enddo

c        bonds

         int=0
         nb=0
         do j=1,nat
           do i=1,j
             if(kbond(i,j).gt.0) then
               nb=nb+1
               int=int+1
               iflagn(int,0)=1
               iflagn(int,1)=i
               iflagn(int,2)=j
               iflagn(int,3)=0
               iflagn(int,4)=0
               iflagn(int,8)=2
               zstardihed(int) = .false.
             endif
           enddo
         enddo

c        angles

         iz=0
         na=0
         do ib=1,nb
         i=iflagn(ib,1)
         j=iflagn(ib,2)
           do k=1,nat
             if(kbond(i,k).gt.0 .and. k.ne.j) then
               Call Gia_ic (int,iflagn,k,i,j,iz,ichi,maxpar)
               if(ichi.eq.0) then
                 na=na+1
                 int=int+1
                 iflagn(int,0)=2
                 iflagn(int,1)=k
                 iflagn(int,2)=i
                 iflagn(int,3)=j
                 iflagn(int,4)=0
                 iflagn(int,8)=3
                 zstardihed(int) = .false.
               endif
             elseif(kbond(j,k).gt.0 .and. k.ne.i) then
               Call Gia_ic (int,iflagn,i,j,k,iz,ichi,maxpar)
               if(ichi.eq.0) then
                 na=na+1
                 int=int+1
                 iflagn(int,0)=2
                 iflagn(int,1)=i
                 iflagn(int,2)=j
                 iflagn(int,3)=k
                 iflagn(int,4)=0
                 iflagn(int,8)=3
               endif
             endif
           enddo
         enddo

c        dihedrals

         nd=0

c        stiff dihedrals first 
         do ia=nb+1,nb+na   ! on angles
            i=iflagn(ia,1)
            j=iflagn(ia,2)
            k=iflagn(ia,3)

            if (kbond(j,k) .gt. 1) then ! stiff died

              do l=1,nat
               if(l.ne.i .and. l.ne.j .and. l.ne.k) then

                if(kbond(i,l).gt.0) then
                  Call Gia_ic (int,iflagn,l,i,j,k,ichi,maxpar)
                  if(ichi.eq.0) then
                   ma=l
                   mb=i
                   mc=j
                   md=k
                   Call Check_Lin (ma,mb,mc,md,nat,geom,aname,maxat,zli)
                   if (zli) then 
                    nd=nd+1
                    int=int+1
                    iflagn(int,0)=3
                    iflagn(int,1)=ma
                    iflagn(int,2)=mb
                    iflagn(int,3)=mc
                    iflagn(int,4)=md
                    iflagn(int,8)=4
                    zstardihed(int) = .false.
                    zharmdihed(int) = .true.
                   endif
                  endif
                elseif(kbond(k,l).gt.0) then
                  Call Gia_ic (int,iflagn,i,j,k,l,ichi,maxpar)
                  if(ichi.eq.0) then
                   ma=i
                   mb=j
                   mc=k
                   md=l
                   Call Check_Lin (ma,mb,mc,md,nat,geom,aname,maxat,zli)
                   if (zli) then 
                    nd=nd+1
                    int=int+1
                    iflagn(int,0)=3
                    iflagn(int,1)=ma
                    iflagn(int,2)=mb
                    iflagn(int,3)=mc
                    iflagn(int,4)=md
                    iflagn(int,8)=4
                    zharmdihed(int) = .true.
                   endif
                  endif
                endif

               endif
              enddo ! on atoms 

            endif ! on stiff dihed 
         enddo  ! on angles



c                           29 August 2006
c        Inversion coordinates for four atoms 
c        (or star dihedrals) are added,for atoms which have 3 bonds
c        In this case the first three atoms of the sequence
c        are bonded, whereas the foruth is bonded to the second
c        and NOT to the third.
c        The logical variable zstardihed is set to true

c        Two definitions are possible: 
c        I   : org
c            The second atom is the one with three bonds
c            In the case of a double bond with the central atom(2) the
c            correspondent atom will be placed n third poistion
c        II  : FF
c            Atom 1 is the one with three bonds
c            that is atom 2 is changed with 1
c        In all cases 1-2-3 are bonded while 3-4 are not

         conv='FF'  !! definition labels: 'mosc' 'ivo '
         do i=1,nat
           nbi=0
           ibi=0
           idop=0
           do k=1,nat
             if(kbond(i,k).gt.0 .and. i.ne.k) then
               nbi=nbi+1
               if(kbond(i,k).gt.1 .and. idop.eq.0) then
                 idop=k
               else
                 ibi=ibi+1
                 iatom(ibi)=k
               endif
             endif
           enddo
           if(nbi.eq.3) then
             nd=nd+1
             int=int+1
             iflagn(int,0)=3
             iflagn(int,8)=4
             zstardihed(int) = .true.
             if(idop.gt.0) then
               iflagn(int,1)=iatom(1)
               iflagn(int,2)=i
               iflagn(int,3)=idop
               iflagn(int,4)=iatom(2)
             else
               iflagn(int,0)=3
               iflagn(int,1)=iatom(1)
               iflagn(int,2)=i
               iflagn(int,3)=iatom(2)
               iflagn(int,4)=iatom(3)
             endif
             if(conv.eq.'FF') then
               iii=iflagn(int,1)
               iflagn(int,1)=iflagn(int,2)
               iflagn(int,2)=iii
             endif
           endif
         enddo
         ncint=int

c        flex dihedrals dihedrals first 
         do ia=nb+1,nb+na  ! loop on angles 
            i=iflagn(ia,1)
            j=iflagn(ia,2)
            k=iflagn(ia,3)

            if (kbond(j,k) .eq. 1) then ! flex died

             do l=1,nat
              if(l.ne.i .and. l.ne.j .and. l.ne.k) then
               if(kbond(i,l).gt.0) then
                Call Gia_ic (int,iflagn,l,i,j,k,ichi,maxpar)
                if(ichi.eq.0) then
                 ma=l
                 mb=i
                 mc=j
                 md=k
                 Call Check_Lin (ma,mb,mc,md,nat,geom,aname,maxat,zli)
                 if (zli) then 
                   nd=nd+1
                   int=int+1
                   iflagn(int,0)=3
                   iflagn(int,1)=ma
                   iflagn(int,2)=mb
                   iflagn(int,3)=mc
                   iflagn(int,4)=md
                   iflagn(int,8)=4
                   zstardihed(int) = .false.
                  endif
                 endif
               elseif(kbond(k,l).gt.0) then
                Call Gia_ic (int,iflagn,i,j,k,l,ichi,maxpar)
                if(ichi.eq.0) then
                 ma=i
                 mb=j
                 mc=k
                 md=l
                 Call Check_Lin (ma,mb,mc,md,nat,geom,aname,maxat,zli)
                 if (zli) then 
                  nd=nd+1
                  int=int+1
                  iflagn(int,0)=3
                  iflagn(int,1)=ma
                  iflagn(int,2)=mb
                  iflagn(int,3)=mc
                  iflagn(int,4)=md
                  iflagn(int,8)=4
                 endif
                endif
               endif
              endif
             enddo

          endif ! on flex dihed  
         enddo  ! on angles 
         ncint=int


c      make a label to each IC

         do int=1,ncint
           ni=iflagn(int,0)
           string=' '
           init=1
           ia=iflagn(int,1)
           do k=1,ni+1
             ib=iflagn(int,k)
             st2='- '
             if(kbond(ia,ib).gt.1) st2='= '
             if(kbond(ia,ib).lt.1) st2='--'
             write(string(init:),'(2a)') st2,aname(ib)
             init=init+8
             ia=ib
           enddo
           string(1:2)=' '
           Call Compact (string,nl)
           labnat(int)=string
         enddo


         label='NATURAL-ICs: style='//conv
         Call Anal_IC (iflagn,valic,labnat,ncint,label
     &                ,nat,kbond,geom,aname,maxat,maxpar)


      ELSEIF (iprt .eq. 1) THEN

c        make a label to each IC according to FF at labels

         do int=1,ncint
           ni=iflagn(int,0)
           string=' '
           init=1
           ia=iflagn(int,1)
           do k=1,ni+1
             ib=iflagn(int,k)
             st2='- '
             if(kbond(ia,ib).gt.1) st2='= '
             if(kbond(ia,ib).lt.1) st2='--'
             write(string(init:),'(2a)') st2,aname(ib)
             init=init+8
             ia=ib
           enddo
           string(1:2)=' '
           Call Compact (string,nl)
           labnat(int)=string
         enddo

c        punch the generated IC's (for help in making starting FF file)


         If (ifftype .eq. 1) then    !! moscito format

          do iat=1,nat
           do jat=1,nat
             zexcl(iat,jat) = .false.
             z14(iat,jat)   = .false.
           enddo
          enddo
          open (unit=14,file='generated.IC.txt')
          istr  = 0
          ibend = 0
          itors = 0
          zero  = 0.d0
          do i=1,ncint
               iat1 = iflagn(i,1)
               iat2 = iflagn(i,2)
               iat3 = iflagn(i,3)
               iat4 = iflagn(i,4)
               if(iflagn(i,0) .eq. 1) then
                  istr = istr + 1
                  if (istr .eq. 1) then 
                      write(14,'(a)') '; Stretchings'
                      write(14,'(a)') 'begin{bonds}'
                  endif
                  write(14,105) 
     &            'harmonic',iat1,iat2,zero,zero,';',i,labnat(i)
               elseif(iflagn(i,0) .eq. 2) then
                  ibend = ibend + 1
                  if (ibend .eq. 1) then
                      write(14,'(a)') 'end{bonds}'
                      write(14,'(a)') '; Bendings'
                      write(14,'(a)') 'begin{angle}'
                  endif    
                  write(14,106) iat1,iat2,iat3,zero,zero,';',i,labnat(i)
               elseif(iflagn(i,0) .eq. 3) then
                  itors = itors + 1
                  if (itors .eq. 1) then
                      write(14,'(a)') 'end{angle}'
                      write(14,'(a)') '; Torsions'
                      write(14,'(a)') 'begin{dihedral}'
                  endif
                  if(zstardihed(i)) then
                  write(14,107) 'improper',iat1,iat2,iat3,iat4,
     &                          zero,zero,'; ',i,labnat(i),'*'
                  else
                     z14(iat1,iat4) = .true.
                     z14(iat4,iat1) = .true.
                  write(14,107) 
     &         'improper',iat1,iat2,iat3,iat4,zero,zero,'; ',i,labnat(i)
                  endif
               endif   
          enddo
          write(14,'(a)') 'end{dihedral}'

c         prints out all possible nonbonded IC

          do iat=1,nat
             do ico = 1,ncint
               if(iflagn(ico,1) .eq. iat) then
                  jat = iflagn(ico,2) 
                  zexcl(iat,jat) = .true.
                  zexcl(jat,iat) = .true.
                  kat = iflagn(ico,3) 
                  lat = iflagn(ico,4) 
                  if (kat .gt. 0) then
                      zexcl(iat,kat) = .true.
                      zexcl(kat,iat) = .true.
                  elseif (lat .gt. 0) then
                      zexcl(iat,lat) = .true.
                      zexcl(lat,iat) = .true.
                  endif
               endif 
             enddo
          enddo

c         1 - 4 interactions

          n14 = 0
          do iat=1,nat-1
           do jat=iat+1,nat
              if(.not. zexcl(iat,jat) .and. z14(iat,jat)) then
               n14= n14 + 1
               if (n14 .eq. 1) then
                   write(14,'(a)') '; Nonbonded terms'
                   write(14,'(a)') 'begin{nonbonded}'
                   write(14,'(a)') '; 1-4 interactions'
               endif
               i1 = iat
               i4 = jat
               if(iat .eq. 7 .and. jat .eq. 10) then
                       write(6,*) 'trovati 2!',iat,jat
                       stop 'ad hoc'
               endif
               if (cr .eq. 2) then 
                 sij = 0.5d0*(sigma(iat)+sigma(jat))
               elseif (cr .eq. 3) then 
                 sij = dsqrt(sigma(iat)*sigma(jat))
               endif
               eij = 0.5d0*dsqrt(epsilon(iat)*epsilon(jat)) ! scaled 1-4
               fqq = 0.833d0
               write(14,108) i1,i4,sij,eij,fqq,';',
     &         aname(i1),'-',aname(i4)
             endif
           enddo
          enddo

c         everything else

          nexcl = 0
          do iat=1,nat-1
           do jat=iat+1,nat
             if(.not. zexcl(iat,jat)) then
               nexcl= nexcl + 1
               if (nexcl .eq. 1) then
                   write(14,'(a)') '; other nb interactions'
               endif
               i1 = iat
               i4 = jat
               if (.not. z14(iat,jat) ) then
                 if (cr .eq. 2) then 
                    sij = 0.5d0*(sigma(iat)+sigma(jat))
                 elseif (cr .eq. 3) then 
                    sij = dsqrt(sigma(iat)*sigma(jat))
                 endif
                 eij = dsqrt(epsilon(iat)*epsilon(jat)) ! scaled 1-4
                 fqq = 1.000d0
                 write(14,108) i1,i4,sij,eij,fqq,';',
     &           aname(i1),'-',aname(i4)
               endif
             endif
           enddo
          enddo
          write(14,'(a)') 'end{nonbonded}'
          write(14,'(a)') 
          write(14,'(a)') 'begin{exclude}'
          write(14,'(a)')  '  all   '
          write(14,'(a)') 'end{exclude}'
          write(14,'(a)') 
          close (unit=14)
        
        elseif (ifftype .eq. 2) then    !! gromacs format
          do iat=1,nat
           do jat=1,nat
             zexcl(iat,jat) = .false.
             z14(iat,jat)   = .false.
           enddo
          enddo

          open (unit=14,file='generated.IC.txt')
          istr  = 0
          ibend = 0
          itors = 0
          istar = 0
          zero  = 0.d0
          ione  = 1
          itwo  = 2

          do i=1,ncint
             iat1 = iflagn(i,1)
             iat2 = iflagn(i,2)
             iat3 = iflagn(i,3)
             iat4 = iflagn(i,4)
             if(iflagn(i,0) .eq. 1) then
                istr = istr + 1
                if (istr .eq. 1) then 
                    write(14,'(a)') '; Stretchings'
                    write(14,'(a)') '[ bonds ]'
                endif
                write(14,115) 
     &          iat1,iat2,ione,zero,zero,';',i,labnat(i)
             elseif(iflagn(i,0) .eq. 2) then
                ibend = ibend + 1
                if (ibend .eq. 1) then
                    write(14,'(a)') 
                    write(14,'(a)') '; Bendings'
                    write(14,'(a)') '[ angles ]'
                endif    
                write(14,116) 
     &          iat1,iat2,iat3,ione,zero,zero,';',i,labnat(i)
             elseif(iflagn(i,0) .eq. 3) then
                if (zharmdihed(i)) then 
                    itors = itors + 1
                    if (itors .eq. 1) then
                      write(14,'(a)') 
                      write(14,'(a)') '; Torsions'
                      write(14,'(a)') '[ dihedrals ]'
                      write(14,'(a)') '; stiff dieds'
                    endif
                    z14(iat1,iat4) = .true.
                    z14(iat4,iat1) = .true.
                    write(14,117) iat1,iat2,iat3,iat4,itwo,
     &                            zero,zero,';',i,labnat(i)
                elseif(zstardihed(i)) then
                    itors = itors + 1
                    if (itors .eq. 1) then
                      write(14,'(a)') 
                      write(14,'(a)') '; Torsions'
                      write(14,'(a)') '[ dihedrals ]'
                    endif
                    istar = istar + 1
                    if (istar .eq. 1) then
                      write(14,'(a)') '; star-like dieds'
                    endif
                    write(14,117) iat1,iat2,iat3,iat4,itwo,
     &                            zero,zero,';',i,labnat(i),'*'
                endif                        
             endif   
          enddo
          if (itors .eq. 0) then 
             write(14,'(a)') 
             write(14,'(a)') '; Torsions'
             write(14,'(a)') '[ dihedrals ]'
          endif
          iflex = 0
          do i=1,ncint
             iat1 = iflagn(i,1)
             iat2 = iflagn(i,2)
             iat3 = iflagn(i,3)
             iat4 = iflagn(i,4)
             if(iflagn(i,0) .eq. 3) then
                if (.not. zharmdihed(i) .and. .not. zstardihed(i)) then
                     iflex          = iflex + 1 
                     if (iflex .eq. 1) then
                      write(14,'(a)') '; flexible dieds'
                     endif
                     itors          = itors + 1
                     z14(iat1,iat4) = .true.
                     z14(iat4,iat1) = .true.
                     write(14,118) 
     &               iat1,iat2,iat3,iat4,ione,zero,zero,itwo,
     &               ';',i,labnat(i)
                endif
             endif   
          enddo
          write(14,'(a)') 

c         prints out all possible nonbonded IC

          do iat=1,nat
             do ico = 1,ncint
               if(iflagn(ico,1) .eq. iat) then
                  jat = iflagn(ico,2) 
                  zexcl(iat,jat) = .true.
                  zexcl(jat,iat) = .true.
                  kat = iflagn(ico,3) 
                  lat = iflagn(ico,4) 
                  if (kat .gt. 0) then
                      zexcl(iat,kat) = .true.
                      zexcl(kat,iat) = .true.
                  elseif (lat .gt. 0) then
                      zexcl(iat,lat) = .true.
                      zexcl(lat,iat) = .true.
                  endif
               endif 
             enddo
          enddo

c         prints out all possible nonbonded IC

          do iat=1,nat
             do ico = 1,ncint
               if(iflagn(ico,1) .eq. iat) then
                  jat = iflagn(ico,2) 
                  zexcl(iat,jat) = .true.
                  zexcl(jat,iat) = .true.
                  kat = iflagn(ico,3) 
                  lat = iflagn(ico,4) 
                  if (kat .gt. 0) then
                      zexcl(iat,kat) = .true.
                      zexcl(kat,iat) = .true.
                  elseif (lat .gt. 0) then
                      zexcl(iat,lat) = .true.
                      zexcl(lat,iat) = .true.
                  endif
               endif 
             enddo
          enddo

c         1 - 4 interactions

          n14 = 0
          do iat=1,nat-1
           do jat=iat+1,nat
              if(.not. zexcl(iat,jat) .and. z14(iat,jat)) then
               n14= n14 + 1
               if (n14 .eq. 1) then
                   write(14,'(a)') '; Nonbonded terms'
                   write(14,'(a)') '[ pairs ]'
                   write(14,'(a)') '; 1-4 interactions'
               endif
               i1 = iat
               i4 = jat
               if (cr .eq. 2) then 
                  sij = 0.5d0*(sigma(iat)+sigma(jat))
               elseif (cr .eq. 3) then 
                  sij = dsqrt(sigma(iat)*sigma(jat))
               endif
               eij = 0.5d0*dsqrt(epsilon(iat)*epsilon(jat)) ! scaled 1-4
               fqq = 0.833d0
               qi  = charge(iat)
               qj  = charge(jat)
               write(14,119) i1,i4,itwo,fqq,qi,qj,sij,eij,';',
     &         aname(i1),'-',aname(i4)
             endif
           enddo
          enddo

c         everything else

          nexcl = 0
          do iat=1,nat-1
           do jat=iat+1,nat
             if(.not. zexcl(iat,jat) ) then
               nexcl= nexcl + 1
               if (nexcl .eq. 1) then
                   write(14,'(a)') '; other nb interactions'
               endif
               i1 = iat
               i4 = jat
               if (.not. z14(iat,jat) ) then
                  if (cr .eq. 2) then 
                     sij = 0.5d0*(sigma(iat)+sigma(jat))
                  elseif (cr .eq. 3) then 
                     sij = dsqrt(sigma(iat)*sigma(jat))
                  endif
                  eij = dsqrt(epsilon(iat)*epsilon(jat)) ! scaled 1-5 or more
                  fqq = 1.d0
                  qi  = charge(iat)
                  qj  = charge(jat)
                  write(14,119) i1,i4,itwo,fqq,qi,qj,sij,eij,';',
     &            aname(i1),'-',aname(i4)
                endif   
             endif
           enddo
          enddo
          write(14,'(a)') 

c         prints out complete exclusion list

          write(14,'(a)') '; Exclusions from default nonbonded '
          write(14,'(a)') '[ exclusions ]'
          write(14,'(a)') '; ai        aj'
          do iat=1,nat-1
           do jat=iat+1,nat
              write(14,'(3x,i3,3x,i3)') iat,jat
           enddo
          enddo


          close (unit=14)
        endif
      ENDIF

c     moscito .system formats

105   format(1x,a8,3x,i3,1x,i3,3x,f4.1,3x,f4.1,15x,a1,2x,i4,2x,a12)
106   format(3x,i3,1x,i3,1x,i3,3x,f4.1,3x,f4.1,9x,a1,2x,i4,2x,a16)
107   format(1x,a8,3x,i3,1x,i3,1x,i3,1x,i3,3x,f4.1,3x,f4.1,3x,a2,
     &       2x,i4,2x,a20)
108   format(4x,i3,1x,i3,3x,f6.3,3x,f9.3,3x,f6.3,5x,a1,
     &       2x,a4,1x,a1,1x,a4)


c     gromacs .top formats

115   format(3x,i3,1x,i3,3x,i1,4x,f4.1,3x,f4.1,15x,a1,2x,i4,2x,a12) 
116   format(3x,i3,1x,i3,1x,i3,3x,i1,4x,3x,f4.1,3x,f4.1,9x,a1,2x,
     &       i4,2x,a16)
117   format(3x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f4.1,3x,f4.1,3x,a1,
     &       2x,i4,2x,a21,2x,a1)
118   format(3x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f4.1,3x,f4.1,3x,i2,
     &       3x,a1,2x,i4,2x,a20,2x,a1)
119   format(4x,i3,1x,i3,3x,i1,3x,f6.3,3x,f6.3,3x,f6.3,3x,f7.4,3x,
     &       f9.5,5x,a1,2x,a4,1x,a1,1x,a4)

      return
      end

c----------------------------------------------------------------------c
      Subroutine Check_Lin (ia,ib,ic,id,nat,geom,aname,maxat,zli)
c----------------------------------------------------------------------c

c     check whether the A-B-C-D dihedral has 3 aligned atoms

      implicit real*8(a-h,o-z)
      dimension geom(maxat,*)
      character aname(*)*5
      data todeg/57.2957795130823d0/
      logical    zli
      tol=0.05d0
      th=180.d0
      kk=0

      zli = .true.
c     check: A-B-C linear ? if so, exchange C-D --> A-B-D-C

      Call Value_Ang (geom,maxat,ia,ib,ic,kk,val1,'d')
      if(dabs(val1-th).lt.tol) then
c     change:  a-b-c-d  -->  a-b-d-c
        Call Value_Ang (geom,maxat,ia,ib,id,kk,val3,'d')
        Call Value_Ang (geom,maxat,ib,id,ic,kk,val4,'d')
        write(6,26) ia,ib,ic,id,ia,ib,id,ic,val1,val3,val4
        if(dabs(val3-th).lt.tol .or. dabs(val4-th).lt.tol) then
          write(6,32)
          zli = .false. 
c         comment the following stop instruction only in specific cases
          stop   
        endif
        ii=ic
        ic=id
        id=ii
      endif

c     check: B-C-D linear ? if so, exchange A-B --> B-A-C-D

      Call Value_Ang (geom,maxat,ib,ic,id,kk,val1,'d')
      if(dabs(val1-th).lt.tol) then
c     change:  a-b-c-d  -->  b-a-c-d
        Call Value_Ang (geom,maxat,ia,ic,id,kk,val3,'d')
        Call Value_Ang (geom,maxat,ib,ia,ic,kk,val4,'d')
        write(6,28) ia,ib,ic,id,ib,ia,ic,id,val1,val3,val4
        if(dabs(val3-th).lt.tol .or. dabs(val4-th).lt.tol) then
          write(6,32)
c         comment the following stop instruction only in specific cases
c         stop   
        endif
        ii=ia
        ia=ib
        ib=ii
      endif

   26 format(/5x,'WARNING: dihedral between three aligned atoms'
     &      ,/5x,'         sequence',4i4,'  changed to',4i4
     &      ,/5x,'         123 original angle',f9.4
     &      ,/5x,'         123 new angle     ',f9.4
     &      ,/5x,'         234 new angle     ',f9.4)
   28 format(/5x,'WARNING: dihedral  between three aligned atoms'
     &      ,/5x,'         sequence',4i4,'  changed to',4i4
     &      ,/5x,'         234 original angle',f9.4
     &      ,/5x,'         234 new angle     ',f9.4
     &      ,/5x,'         123 new angle     ',f9.4)
   32 format(/5x,'WARNING: this change does not work !!!',//,
     &        5x,'Do not use this dihedral definition in the FF ')
c  32 format(/5x,'ERROR  : this change does not work !!!',//)
c
      return
      end

c----------------------------------------------------------------------c
      Subroutine Anal_IC (iflagr,valic,labric,nric,labe
     &                   ,nat,kbond,geom,aname,maxat,maxpar)
c----------------------------------------------------------------------c

c     analyse the internal coordinates

      implicit real*8(a-h,o-z)
      dimension iflagr(maxpar,0:8),geom(maxat,*),kbond(maxat,*)
     &         ,valic(*)
      character aname(*)*5,labric(*)*(*),labe*(*),type(3)*8
      data type/'distance','angle   ','dihedral'/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      ll=Lenw(labe)

      write(6,44) labe(1:ll)
      do ico=1,nric
        it=iflagr(ico,0)
        nt=iflagr(ico,8)   !!   atoms involved
        ia=iflagr(ico,1)
        ib=iflagr(ico,2)
        if(nt.eq.2) then
          val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &              +(geom(ia,2)-geom(ib,2))**2
     &              +(geom(ia,3)-geom(ib,3))**2 )*boh2ang
        else
          ic=iflagr(ico,3)
          id=iflagr(ico,4)
          Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'d')
        endif
        valic(ico)=val
        write(6,48) ico,type(iflagr(ico,0)),labric(ico),val
     &             ,(iflagr(ico,k),k=1,nt)
        do i=1,nt-1
          ia=iflagr(ico,i)
          ib=iflagr(ico,i+1)
          if(kbond(ia,ib).le.0) write(6,50) ia,ib
        enddo
        do i=2,nt
          do j=1,i-1
            ia=iflagr(ico,i)
            ib=iflagr(ico,j)
            if(ia.eq.ib) then
              write(6,'(/" ERROR: atom repeated in a RIC",//)')
              stop
            endif
          enddo
        enddo
      enddo

   44 format(//10x,'====== INTERN COORD analysis:  ',a,'  ======'
     &       ,/10x,'TYPE     NAME            EQUIL.VALUE'
     &       ,/10x,'                         Angs or deg')
   46 format(10x,i4,4x,'coupling',2i4,5x,2a)
   48 format(10x,i4,2x,a,2x,a,f11.4,4x,4i4)
   50 format(10x,'warning: atoms',2i5,'  not bonded')
      return
      end

c----------------------------------------------------------------------c
      Subroutine Generate_Label (npar,iflag,labpar,aname,maxpar
     &                          ,kbond,maxat)
c----------------------------------------------------------------------c

c     given the iflag array, this routine generate a label for
c     each parameter of the FF

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),kbond(maxat,*),aname(*)
      character labpar(*)*(*),aname*5,string*66,st2*2,st5*5

      do 80 ip=1,npar
      it=iflag(ip,0)
      nt=iflag(ip,8)
      string=' '
      if(it.eq.8) goto 50

c     diagonal term

      ia=iflag(ip,1)
      init=1
      do k=1,nt
        ib=iflag(ip,k)
        st2='- '
        if(kbond(ia,ib).gt.1) st2='= '
        if(kbond(ia,ib).le.0) st2='--'
        if(iflag(ip,0).eq.5)  st2='--'   !! intra LJ
        write(string(init:),'(2a)') st2,aname(ib)
        init=init+8
        ia=ib
      enddo
      if(it.eq.4) write(string(init:),'("_n=",i2)') iflag(ip,5)
      string(1:1)=' '
      if(string(2:2).eq.'-') string(2:2)=' '
      Call Compact (string,nl)
      labpar(ip)=string
      goto 80

c     off diagonal term

   50 continue
      write(string(1:),'("coupl-",i4,"-",i4)') (iflag(ip,k),k=1,2)
      labpar(ip)=string
   80 continue

      do ip = 1,Npar   
       do jp = 1,Npar
        if (iflag(ip,0) .eq. 1 .and. iflag(jp,0) .eq. 1) then ! stretch
            ia1 = iflag(ip,1) 
            ia2 = iflag(jp,1) 
            ib1 = iflag(ip,2) 
            ib2 = iflag(jp,2) 
            if (aname(ia1) .eq. aname(ib2) .and. 
     &          aname(ia2) .eq. aname(ib1) ) labpar(jp) = labpar(ip)
        elseif (iflag(ip,0) .eq. 2 .and. iflag(jp,0) .eq. 2) then ! bend
            ia1 = iflag(ip,1) 
            ia2 = iflag(jp,1) 
            ib1 = iflag(ip,2) 
            ib2 = iflag(jp,2) 
            ic1 = iflag(ip,3) 
            ic2 = iflag(jp,3) 
            if (aname(ib1) .eq. aname(ib2)) then
              if (aname(ia1) .eq. aname(ic2) .and. 
     &            aname(ia2) .eq. aname(ic1) ) labpar(jp) = labpar(ip)
            endif
        elseif (iflag(ip,0) .eq. 3 .and. iflag(jp,0) .eq. 3) then ! harm tors
            ia1 = iflag(ip,1) 
            ia2 = iflag(jp,1) 
            ib1 = iflag(ip,2) 
            ib2 = iflag(jp,2) 
            ic1 = iflag(ip,3) 
            ic2 = iflag(jp,3) 
            id1 = iflag(ip,4) 
            id2 = iflag(jp,4) 
            if (aname(ib1) .eq. aname(ic2) .and. 
     &          aname(ic1) .eq. aname(ib2)) then
              if (aname(ia1) .eq. aname(id2) .and. 
     &            aname(ia2) .eq. aname(id1) ) labpar(jp) = labpar(ip)
            endif
        elseif (iflag(ip,0) .eq. 4 .and. iflag(jp,0) .eq. 4) then ! cos  tors
            ia1 = iflag(ip,1) 
            ia2 = iflag(jp,1) 
            ib1 = iflag(ip,2) 
            ib2 = iflag(jp,2) 
            ic1 = iflag(ip,3) 
            ic2 = iflag(jp,3) 
            id1 = iflag(ip,4) 
            id2 = iflag(jp,4) 
            if (iflag(ip,5) .eq. iflag(jp,5))  then  ! same n value
             if (aname(ib1) .eq. aname(ic2) .and. 
     &          aname(ic1) .eq. aname(ib2)) then
              if (aname(ia1) .eq. aname(id2) .and. 
     &            aname(ia2) .eq. aname(id1) ) labpar(jp) = labpar(ip)
             endif
            endif
        else
        endif
       enddo
      enddo
      return
      end

c----------------------------------------------------------------------c
      Subroutine Value_IC (iflag,valic,npar,geom,maxat,maxpar)
c----------------------------------------------------------------------c

c         compute the value of the RIC as given in iflag

c     The data are computed for each function AND NOT for each RIC,
c     therefore there may be repetitions (above all for torsional func)

c     warning: coordinates must be in bohr

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),geom(maxat,*),valic(*)
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs

      do ip=1,npar
        it=iflag(ip,0)                  !! type 
        nt=iflag(ip,8)                  !! atoms involv.
        ia=iflag(ip,1)
        ib=iflag(ip,2)
        if(it.lt.8) then                !! couplings
          if(nt.eq.2) then
            val=dsqrt( (geom(ia,1)-geom(ib,1))**2
     &                +(geom(ia,2)-geom(ib,2))**2
     &                +(geom(ia,3)-geom(ib,3))**2 )
          else
            ic=iflag(ip,3)
            id=iflag(ip,4)
            Call Value_Ang (geom,maxat,ia,ib,ic,id,val,'r')
          endif
          valic(ip)=val
        endif
      enddo
      return
      end

c----------------------------------------------------------------------c
      Subroutine Set_ParFun (iflag,valic0,parfun,parfun2,
     &                       npar,ikeepff,maxpar)
c----------------------------------------------------------------------c

c   set the adequate fixed parameters entering the functions
c   expressed in one or more IC

c      These are:
 
c     iflag(i,0)=1  nt=2   harmonic stretching               R0
c     iflag(i,0)=2  nt=3   harmonic bending                  Ang0
c     iflag(i,0)=3  nt=4   harmonic torsion                  R0
c     iflag(i,0)=4  nt=4   torsione propria (cosines)        Gamma
c     iflag(i,0)=5  nt=2   LJ intramolecolare                sigma
c     iflag(i,0)=8         coupling                          Gamma (if torsion is present)

c     ikeepff(i,j)=1      keep the forcefield data for parmaeter i (j=1,2)

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),valic0(*),parfun(*),parfun2(*),
     &          ikeepff(maxpar,2)
      data todeg/57.2957795130823d0/



      do 40 ip=1,npar
        it=iflag(ip,0)            !! function type
        nt=iflag(ip,8)            !! num of atoms involved
        ia=iflag(ip,1)            !! atom n.1
        ib=iflag(ip,2)            !! atom n.2
        if(it .eq.5) then         !! sigma is always kept
           ikeepff(ip,1) = 1
           ikeepff(ip,2) = 1
        endif
c       sets either keepff or valic0 values
        if(it.lt.8) then  ! uncoupled coord.

          if(ikeepff(ip,1).eq.1) goto 40
          if(it.eq.4) then
             n_ip    = iflag(ip,5)
             phi0    = valic0(ip)*todeg
             gnew    = dble(n_ip)*phi0 + 540.d0 
             parfun(ip) = dmod(gnew,360.d0)/todeg
          else
             parfun(ip)=valic0(ip)
          endif

        elseif(it.eq.8) then   !! couplings

          if(ikeepff(ip,1) .ne. 1) then    ! computes the first nonlinear prm in coupling
           iuc1 = iflag(ip,1)   ! uncoupled IC corresponding to the first coupled IC
           if (iflag(iuc1,0) .eq. 4) then ! proper torsion
             n_ip = iflag(ip,5) 
             phi0 = valic0(iuc1)*todeg
             gnew = dble(n_ip)*phi0 + 540.d0 
             parfun(ip) = dmod(gnew,360.d0)/todeg ! => sin(n*phi-gam) = 0
           else                          ! harmonic  torsion
             parfun(ip) = valic0(iuc1)   ! from QM geometry
           endif
          endif

          if(ikeepff(ip,2).ne. 1) then ! computes the second nonlinear prm in coupling

           iuc2 = iflag(ip,2)   ! uncoupled IC corresponding to the first coupled IC
           if (iflag(iuc2,0) .eq. 4) then ! proper torsion
             n_ip = iflag(ip,5) 
             phi0 = valic0(iuc2)*todeg
             gnew = dble(n_ip)*phi0 + 540.d0 
             parfun2(ip) = dmod(gnew,360.d0)/todeg ! => sin(n*phi-gam) = 0
           else                          ! harmonic  torsion
             parfun2(ip) = valic0(iuc2)   ! from QM geometry
           endif
c           if (iflag(ip,8) .eq. 8) then   ! two coupled torsions 

c            iuc2 = iflag(ip,2)   ! uncoupled IC corresponding to the second coupled IC
c            if (iuc2 .eq. 3) then ! harmonic tors - harmonic tors coupling
c               parfun2(ip) = valic0(iuc2)    ! from QM geometry
c            else                  ! either harmonic/proper or proper/proper
c               n_ip = iflag(ip,6) 
c               phi0 = valic0(iuc2)*todeg
c               gnew = dble(n_ip)*phi0 + 540.d0 ! => sin(n*phi-gam) = 0
c               parfun2(ip) = dmod(gnew,360.d0)/todeg
c            endif

c           elseif (iflag(ip,8) .eq. 7) then   ! bending-torsion

c            iuc2 = iflag(ip,2)   ! uncoupled IC corresponding to the second coupled IC
c            if (iuc2 .eq. 3) then ! bending  - harmonic tors 
c               parfun2(ip) = valic0(iuc2)    ! from QM geometry
c            elseif (iuc2 .eq. 4) then  ! bending - proper tors
c               n_ip = iflag(ip,5) 
c               phi0 = valic0(iuc2)*todeg
c               gnew = dble(n_ip)*phi0 + 540.d0 ! => sin(n*phi-gam) = 0
c               parfun2(ip) = dmod(gnew,360.d0)/todeg
c            endif

c           elseif (iflag(ip,8) .eq. 6) then   ! stretching-torsion (or bend-bend)

c            iuc2 = iflag(ip,2)   ! uncoupled IC corresponding to the second coupled IC
c            if (iuc2 .eq. 3) then ! strewtching  - harmonic tors 
c               parfun2(ip) = valic0(iuc2)    ! from QM geometry
c            elseif (iuc2 .eq. 4) then  ! stretching - prtoper tors
c               n_ip = iflag(ip,5) 
c               phi0 = valic0(iuc2)*todeg
c               gnew = dble(n_ip)*phi0 + 540.d0 
c               parfun2(ip) = dmod(gnew,360.d0)/todeg
c            endif
c           else          ! str-str, str-bnd and bnd
c               parfun2(ip) = valic0(iuc2)    ! from QM geometry
c           endif
          endif

        endif  ! uncoupled or coupled
   40 continue
      return
      end

c----------------------------------------------------------------------c
      Subroutine Anal_FF (iflag,par,valic,parfun,labpar,npar
     &                   ,labe,kbond,aname,maxat,maxpar)
c----------------------------------------------------------------------c

c     analyze the internal Force Field

c     warning: coordinates must be in bohr

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),kbond(maxat,*)
     &         ,par(*),parfun(*),valic(*)
      character aname(*)*5,labpar(maxpar)*(*),labe*(*),type(5)*8
      common /printout/ levprt
      data type/'stretch ','bending ','unpr-tor','pr-tors ','LJ-nbd  '/
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      ll=Lenw(labe)

      if(levprt.gt.0) write(6,44) labe(1:ll)
      ni=0
      do ip=1,npar
        it=iflag(ip,0)   !! function type
        nt=iflag(ip,8)   !! num of atoms involved
        ia=iflag(ip,1)
        ib=iflag(ip,2)
        ni=ni+1
        if(it.eq.8) then                       !! couplings
          if(levprt.gt.0)
     &    write(6,46) ip,ia,ib,labpar(ia),labpar(ib)
          ni=ni-1
        else                                   !! others
          if(nt.eq.2) then
            val0=valic(ip)*boh2ang
            parf=parfun(ip)*boh2ang
          else
            val0=valic(ip)*todeg
            parf=parfun(ip)*todeg
          endif
          if(levprt.gt.0)
     &    write(6,48) ip,type(iflag(ip,0)),labpar(ip)
     &               ,par(ip),parf,val0,(iflag(ip,k),k=1,nt)
          if(iflag(ip,0).ne.5) then
            do i=1,nt-1
              ia=iflag(ip,i)
              ib=iflag(ip,i+1)
              if(levprt.gt.0 .and.
     &           kbond(ia,ib).le.0) write(6,50) ia,ib
            enddo
          endif
          do i=2,nt
            do j=1,i-1
              ia=iflag(ip,i)
              ib=iflag(ip,j)
              if(ia.eq.ib) then
                write(6,'(/" ERROR: atom considered two times",//)')
                write(6,'(/"        On function : ",//)')
                write(6,*) ip, labpar(ip)
                stop
              endif
            enddo
          enddo
        endif
      enddo
      if(levprt.gt.0) write(6,52) ni,npar-ni

   44 format(//14x,'====== FORCE FIELD analysis:  ',a,'  ======'
     &       ,/8x,'TYPE     NAME',17x,'Linear Par  Fixed Par.'
     &                                            ,'  IC value'
     &       ,/8x,'             ',17x,' H(bohr,r)  angs,degr '
     &                                            ,' angs degr')
   46 format(i4,2x,'coupling',2i4,6x,2a)
   48 format(i4,2x,a,2x,a,f11.4,2f11.4,4x,4i4)
   50 format(9x,'warning: atoms',2i5,'  not bonded')
   52 format(19x,'there are',i4,' func.s and',i4,' couplings',/)
      return
      end

c----------------------------------------------------------------------c
      Subroutine Recover_RIC (npar,iflag,labpar,indir
     &                       ,nric,iflagr,labric,maxpar)
c----------------------------------------------------------------------c

c     set the IFLAGR array containing information on the redundant
c     internal coordinates (RICO) from the array (iflag) specific for
c     the parameters of the force fields.
c     In practice, the terms referring to the same
c     RICO as for the dihedral potentials given as a sum of cosines,
c     are eliminated. 
c        stretching , LJ, Coulomb      ---> distance
c        bending                       ---> angle
c        torsion proper and improper   ---> dihedral

c     UPDATE  may,24,2006
c      iflagr(ip,0)   1=distance   2=angle   3=dihedral
c      iflagr(ip,8)   2=distance   3=angle   4=dihedral
c     IMPORTANT NOTE:  iflag and iflagr have the same content at loc=8
c                      i.e. the number of atoms involved.
c     IMPORTANT NOTE:  internal charge-charge interaction is supposed
c                      to be set on atom pairs which ALREADY interact
c                      with LJ 

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8),iflagr(maxpar,0:8),indir(*)
      character labpar(*)*(*),labric(*)*(*)

      nl=Len(labric(1))
      ic=0
      do ip=1,npar
        irep=0
        if(iflag(ip,0).gt.7) irep=8
        if(ip.gt.1) then
          if(iflag(ip,0) .ne. 6) then  ! not dealing with charge-charge     
             do i=0,4
               if(iflag(ip,i).eq.iflag(ip-1,i)) irep=irep+1
             enddo
          elseif(iflag(ip,0) .eq. 6) then  ! dealing with charge-charge     
             irep = 8
          endif
        endif
        if(irep.lt.5) then
          ic=ic+1
          do i=0,8
            iflagr(ic,i)=iflag(ip,i)
          enddo
          iflagr(ic,5)=0
          labric(ic)=labpar(ip)
          ii=index(labric(ic),'_')
          if(ii.gt.0) labric(ic)(ii:nl)=' '
          nt=iflag(ip,8)    !! involved atoms
          iflagr(ic,0)=nt-1
          indir(ic)=ip
        endif
      enddo
      nric=ic

      return
      end

c----------------------------------------------------------------------c
      Subroutine Gia_ic (ncint,iflag,i,j,k,l,iquale,maxpar)
c----------------------------------------------------------------------c

c  given the atoms i,j,k,l this routine gives the position of the
c  corresponding internal coordinate in the iflag array

c  iquale=0  the i,j,k,l coord has not been found in the iflag array

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8)
      iquale=0

c     bond

      if(k.eq.0) then
        do int=1,ncint
          nok=0
          if(iflag(int,1).eq.i) nok=nok+1
          if(iflag(int,2).eq.j) nok=nok+1
          if(iflag(int,8).eq.2) nok=nok+1
          if(nok.eq.3) then
            iquale=int
            goto 40
          endif
          nok=0
          if(iflag(int,1).eq.j) nok=nok+1
          if(iflag(int,2).eq.i) nok=nok+1
          if(iflag(int,8).eq.2) nok=nok+1
          if(nok.eq.3) then
            iquale=int
            goto 40
          endif
        enddo
        return
      endif

c     angle

      if(l.eq.0) then
        do int=1,ncint
          nok=0
          if(iflag(int,1).eq.i) nok=nok+1
          if(iflag(int,2).eq.j) nok=nok+1
          if(iflag(int,3).eq.k) nok=nok+1
          if(iflag(int,8).eq.3) nok=nok+1
          if(nok.eq.4) then
            iquale=int
            goto 40
          endif
          nok=0
          if(iflag(int,1).eq.k) nok=nok+1
          if(iflag(int,2).eq.j) nok=nok+1
          if(iflag(int,3).eq.i) nok=nok+1
          if(iflag(int,8).eq.3) nok=nok+1
          if(nok.eq.4) then
            iquale=int
            goto 40
          endif
        enddo
        return
      endif

c     dihedral

      do int=1,ncint
        nok=0
        if(iflag(int,1).eq.i) nok=nok+1
        if(iflag(int,2).eq.j) nok=nok+1
        if(iflag(int,3).eq.k) nok=nok+1
        if(iflag(int,4).eq.l) nok=nok+1
        if(iflag(int,8).eq.4) nok=nok+1
        if(nok.eq.5) then
          iquale=int
          goto 40
        endif
        nok=0
        if(iflag(int,1).eq.l) nok=nok+1
        if(iflag(int,2).eq.k) nok=nok+1
        if(iflag(int,3).eq.j) nok=nok+1
        if(iflag(int,4).eq.i) nok=nok+1
        if(iflag(int,8).eq.4) nok=nok+1
        if(nok.eq.5) then
          iquale=int
          goto 40
        endif
      enddo
      return

   40 return
      end

c----------------------------------------------------------------------c
      Subroutine Compare_IC (iflagr,labric,ncintr,labr
     &                      ,iflagp,labrin,ncintp,labp,max3n,maxpar)
c----------------------------------------------------------------------c

c     compare RIC and natural-IC

c      RIC (moscito) are in:  iflagr
c      NATURAL       are in:  iflagp

      implicit real*8(a-h,o-z)
      dimension iflagr(maxpar,0:8),iflagp(maxpar,0:8)
     &         ,iokr(3000),iokp(3000)
      character labric(*)*15,labrin(*)*15,labr*(*),labp*(*)

      do intp=1,ncintp
        iokp(intp)=0
      enddo
      do intr=1,ncintr
        iokr(intr)=0
      enddo

      do intp=1,ncintp
        ip=iflagp(intp,1)
        jp=iflagp(intp,2)
        kp=iflagp(intp,3)
        lp=iflagp(intp,4)
        Call Gia_ic (ncintr,iflagr,ip,jp,kp,lp,ichi,maxpar)
        if(ichi.gt.0) then
          iokp(intp)=ichi
          iokr(ichi)=intp
        endif
      enddo

      write(6,*) ' '
      ndif=0
      do i=1,ncintr
        if(iokr(i).eq.0) then
          write(6,30) labr,i,(iflagr(i,k),k=1,4),labric(i),labp
          ndif=ndif+1
        endif
      enddo
      do i=1,ncintp
        if(iokp(i).eq.0) then
          write(6,30) labp,i,(iflagp(i,k),k=1,4),labrin(i),labr
          ndif=ndif+1
        endif
      enddo
      write(6,'(9x,"found",i4,"  differences",/2x)') ndif
      write(6,*) ' '
   30 format(1x,a,' n.',i4,3x,4i3,3x,a,'  not found in ',a)

      return
      end

c----------------------------------------------------------------------c
      Subroutine Punch_RICO (iflagr,labric,nric,valric,ngeom
     &                      ,maxpar,maxgeom)
c----------------------------------------------------------------------c

c  analyse the internal coordinates at the geometries
c  considered for the force field fitting

      implicit real*8(a-h,o-z)
      dimension iflagr(maxpar,0:8),valric(0:maxgeom,*)
     &         ,indep(8),fmul(8),vmin(8),vmax(8)
      character labric(*)*(*),card*90
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      b=boh2ang

      ianal=0
      Call Cerca_W ('$analric',5,card,ier,.true.)
      if(ier.ne.0) return
      Call Count_W (card,nwrd)
      if(nwrd.gt.5) nwrd=3
      if(nwrd.le.0) return
      read(card,*) (indep(i),i=1,nwrd)
      do i=1,nwrd
        if(indep(i).lt.0 .or. indep(i).gt.nric) goto 110
      enddo
      write(6,62) (i,indep(i),labric(indep(i)),i=1,nwrd)
      write(6,64)

c     print the values of the RICO's at all the geometries

      nco=5
      lsup=0
      dowhile (lsup.lt.nric)
        linf=lsup+1
        lsup=lsup+nco
        if(lsup.gt.nric) lsup=nric
        nw=lsup-linf+1
        write(6,70)  (labric(i),i=linf,lsup)
        write(6,72)  (i,i=linf,lsup)
        j=0
        do l=linf,lsup
          j=j+1
          fmul(j)=1.0d0
          if(iflagr(l,0).gt.1) fmul(j)=todeg
          xmin=1.d9
          xmax=-1.d9
          do k=0,ngeom
            if(valric(k,l).lt.xmin) xmin=valric(k,l)
            if(valric(k,l).gt.xmax) xmax=valric(k,l)
          enddo
          vmin(j)=xmin*fmul(j)
          vmax(j)=xmax*fmul(j)
        enddo
        do k=0,ngeom
          write(6,74) k,(fmul(j-linf+1)*valric(k,j)*b,j=linf,lsup)
        enddo
        write(6,75)
        write(6,76) (vmin(j),j=1,nw)
        write(6,78) (vmax(j),j=1,nw)
        write(6,80) (vmax(j)-vmin(j),j=1,nw)
      enddo

   62 format(///,9x,54('*'),/20x,'Analysis of Redundant Int Coord.s'
     &         ,/22x,'at the considered Geometries'
     &         ,/24x,'independent variable are:',/(i26,i5,3x,a))
   64 format(9x,54('*'))
   70 format(/8x,5(2x,a15))
   72 format(i15,5i17)
   74 format(' geom',i2,f10.3,5f17.3)
   75 format(5(11x,'------'))
   76 format(' min ',2x,f10.3,5f17.3)
   78 format(' max ',2x,f10.3,5f17.3)
   80 format(' max-min' ,f9.3,5f17.3)
      return

  110 write(6,112) indep(i)
  112 format(/8x,'input error in $analric: RIC',i5,' wrong',//)
      stop
      end

c----------------------------------------------------------------------c
      Subroutine Trasf2Int (geom,nat,U,iflag,ncint,maxat,max3n,maxpar)
c----------------------------------------------------------------------c

c     compute the transformation matrix for the geometry
c     given in /geom/ for each internal coordinates in terms of
c     cartesian displacements
c     Angles (coord.) derivatives are in radians (bohr)

c                q(ncint) = U(acint,3*nat) * X(3*nat)

c     WARNING: the two type of iflag arrays used in the program i.e.
c        iflag    referred to the functions
c        iflagr   referred to the redund intern coordinates
c     can be used without distinction as input data for this routine
c     since they are identical in the location used 
c     For instance: iflag(ic,8) = 2,3,4  = n. of atoms involved in
c     the internal coordinate ic,  is the same in the two types of array.

      implicit real*8(a-h,o-z)
      parameter (max=200)
      dimension geom(maxat,*),U(max3n,*),iflag(maxpar,0:8)
     &         ,faux(max,5),gaux(max,5),ddr(4),valu(4)
      data rad2deg /57.2957795130823229d0/

      Call Check_Dim (nat,max,'TRASF2INT: increase max')
      dr=1.d-04
      dnorm=12.0d0*dr
      ddr(1)=-dr-dr
      ddr(2)=-dr
      ddr(3)=+dr
      ddr(4)=+dr+dr
      do int=1,ncint
        do j=1,3*nat
          U(int,j)=0.0d0
        enddo
      enddo
      do i=1,nat
        do k=1,4
          gaux(i,k)=geom(i,k)
          faux(i,k)=geom(i,k)
        enddo
      enddo

c     DISTANCES 

      do 20 int=1,ncint
        if(iflag(int,8).ne.2) goto 20
        ia=iflag(int,1)
        ib=iflag(int,2)
        do iat=1,2       !! loop on  atoms A e B
          jj=iflag(int,iat)
          iadd=3*(jj-1)
          do k=1,3       !! loop components x,y,z
            do ispo=1,4  !! loop on displacements (jj,k) for derivative
              do l=1,3
                gaux(jj,l)=faux(jj,l)
              enddo
              gaux(jj,k)=gaux(jj,k)+ddr(ispo)
              val=dsqrt((gaux(ia,1)-gaux(ib,1))**2
     &                 +(gaux(ia,2)-gaux(ib,2))**2
     &                 +(gaux(ia,3)-gaux(ib,3))**2)
              valu(ispo)=val
            enddo
            der=(valu(1)-8.d0*valu(2)+8.d0*valu(3)-valu(4))/dnorm
            U(int,iadd+k)=der
            do l=1,3
              gaux(jj,l)=faux(jj,l)
            enddo
          enddo
        enddo
   20 continue

c     ANGLES AND DIHEDRALS

      do 30 int=1,ncint
        if(iflag(int,8).eq.2) goto 30
        n=iflag(int,0)
        ia=iflag(int,1)
        ib=iflag(int,2)
        ic=iflag(int,3)
        id=iflag(int,4)
        do iat=1,Min0(4,n+1) !! loop on  atoms (n+1 = 3 o 4)
          jj=iflag(int,iat)
          iadd=3*(jj-1)
          do k=1,3       !! loop on x,y,z
            do ispo=1,4  !! loop on displacements (jj,k) for derivative comp.
              do l=1,3
                gaux(jj,l)=faux(jj,l)
              enddo
              gaux(jj,k)=gaux(jj,k)+ddr(ispo)
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
            U(int,iadd+k)=der/rad2deg
            do l=1,3
              gaux(jj,l)=faux(jj,l)
            enddo
          enddo
        enddo
   30 continue

      return
      end

c----------------------------------------------------------------------c
      Subroutine Redund_Int (U,F,C,ncint,nat,labpar,nrid,thr0,max3n)
c----------------------------------------------------------------------c

c     compute the redundancy test on the internal coordinates
c     by diagonalizing the U*(U+) matrix U(a,i)=(dqa/dxi)


      implicit real*8(a-h,o-z)
      dimension U(max3n,*),F(*),C(*)
      dimension ndsym(5)
      character title*50,labpar(*)*15

      n=3*nat
      nsym=1
      ncol=10

c     redundancy test

      ij=0
      do jnt=1,ncint
        do int=1,jnt
          fij=0.0d0
          do k=1,n
            fij=fij+U(jnt,k)*U(int,k)
          enddo
          ij=ij+1
          F(ij)=fij
        enddo
      enddo
      ndsym(1)=ncint
      title='U*U(+) Matrix : redundancy of internal coordinates'
      Call Pr_Tri (F,nsym,ndsym,title,ncol)
      Call Eigen (F,C,ncint)
      write(6,30) (i,F(i),i=1,ncint)
      nrid=0
      do int=1,ncint
        if(F(int).lt.1.d-03) then
          ib=(int-1)*ncint
          write(6,32) F(int),(k,labpar(k),C(k+ib),k=1,ncint)
          nrid=nrid+1
        endif
      enddo
      if(nrid.eq.0) write(6,34) ncint

   30 format(/'    TEST ON INTERNAL COORDINATES SINGULARITY'
     &      ,/'    eigenvalues of the U*U(T) matrix',/(4(i6,f12.4)))
   32 format(/'  WARNING: low eigenvalue found',f15.8
     &      ,/(i6,2x,a,2x,f10.5))
   34 format( '    OK: no small eigenvalue found:   n. NRIC =',i4)

      return
      end

c----------------------------------------------------------------------c
      Subroutine Fstop (string)
c----------------------------------------------------------------------c

      character cardf*90,string*(*)
      common /forstop/ cardf
      ls=Lenw(string)
      lc=Lenw(cardf)
      write(6,4) string(1:ls),cardf(1:lc)
    4 format(//4x,'ERROR: ',a,/4x,'last card /',a,'/',/,/)
      stop
      end

c----------------------------------------------------------------------c
      Subroutine Generate_deps (ncint,iflag,labpar,atlab,maxpar,maxat)
c----------------------------------------------------------------------c

c     given the iflag array ant the atlab labels 
c     this routine generates suggested dependencies as
c     dictated by symmetry of labels

      implicit real*8(a-h,o-z)
      dimension iflag(maxpar,0:8)
      character labpar(*)*(*),atlab(maxat)*5,string*66,st2*2,st5*5,
     &   ia1*5, ia2*5, ia3*5, ia4*5, ib1*5, ib2*5, ib3*5, ib4*5,
     &   ja1*5, ja2*5, ja3*5, ja4*5, jb1*5, jb2*5, jb3*5, jb4*5
      logical   zdep(maxpar)

      open(35,file='suggdeps.txt' ,form='formatted',status='unknown')
      write(35,'(a)') '$dependence 1.2'
      
      do ico = 1,Ncint
         zdep(ico) = .false.
      enddo

      do ico = 1,Ncint-1
        do jco = 1+ico,Ncint
           if(zdep(jco)) goto 5
           ia = iflag(ico,1)
           ib = iflag(ico,2)
           ic = iflag(ico,3)
           id = iflag(ico,4)
           ja = iflag(jco,1)
           jb = iflag(jco,2)
           jc = iflag(jco,3)
           jd = iflag(jco,4)

c          considers two IC of the same type (str, bend, etc.)

           if (iflag(ico,0) .eq. iflag(jco,0)) then

c             iflag(k,0)=1 --> stretching
              if (iflag(ico,0) .eq. 1) then    
                  if(atlab(ia) .eq. atlab(ja) .and.
     &               atlab(ib) .eq. atlab(jb)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                  elseif(atlab(ia) .eq. atlab(jb) .and.
     &                   atlab(ib) .eq. atlab(ja)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                   endif

c             iflag(k,0)=2 --> bending
              elseif (iflag(ico,0) .eq. 2) then    
                  if(atlab(ia) .eq. atlab(ja) .and.
     &               atlab(ic) .eq. atlab(jc) .and.
     &               atlab(ib) .eq. atlab(jb)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                  elseif(atlab(ia) .eq. atlab(jc) .and.
     &               atlab(ic) .eq. atlab(ja) .and.
     &               atlab(ib) .eq. atlab(jb)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                   endif

c             iflag(k,0)=3 --> harmonic torsions
              elseif (iflag(ico,0) .eq. 3) then    
                  if(atlab(ia) .eq. atlab(ja) .and.
     &               atlab(ib) .eq. atlab(jb) .and.
     &               atlab(ic) .eq. atlab(jc) .and.
     &               atlab(id) .eq. atlab(jd)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                  elseif(atlab(ia) .eq. atlab(jd) .and.
     &               atlab(ib) .eq. atlab(jc) .and.
     &               atlab(ic) .eq. atlab(jb) .and.
     &               atlab(id) .eq. atlab(ja)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                  endif

c             iflag(k,0)=4 --> proper torsions (cosines)
              elseif (iflag(ico,0) .eq. 4) then    

c                 same n on cosine term
                 if(iflag(ico,5) .eq. iflag(jco,5)) then
                     if(atlab(ia) .eq. atlab(ja) .and.
     &                  atlab(ib) .eq. atlab(jb) .and.
     &                  atlab(ic) .eq. atlab(jc) .and.
     &                  atlab(id) .eq. atlab(jd)) then
                        zdep(jco) =.true.
                        write(35,105) jco,'=',ico,'*1.d0',';',
     &                                labpar(jco),'=',labpar(ico)
                     elseif(atlab(ia) .eq. atlab(jd) .and.
     &                  atlab(ib) .eq. atlab(jc) .and.
     &                  atlab(ic) .eq. atlab(jb) .and.
     &                  atlab(id) .eq. atlab(ja)) then
                        zdep(jco) =.true.
                        write(35,105) jco,'=',ico,'*1.d0',';',
     &                                labpar(jco),'=',labpar(ico)
                     endif
                  endif

c             iflag(k,0)=5 --> non-bonded interactions
              elseif (iflag(ico,0) .eq. 5) then    
                  if(atlab(ia) .eq. atlab(ja) .and.
     &               atlab(ib) .eq. atlab(jb)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                  elseif(atlab(ia) .eq. atlab(jb) .and.
     &                   atlab(ib) .eq. atlab(ja)) then
                     zdep(jco) =.true.
                     write(35,105) jco,'=',ico,'*1.d0',';',
     &                             labpar(jco),'=',labpar(ico)
                   endif
                  


c             iflag(k,0)=8 --> cross 
              elseif (iflag(ico,0) .eq. 8) then    
c               at first assign coupling type and define atoms with
c               notation:   i, j are the couplings
c                           a, b are the ICs coupled by a coupling
c                           1, 2, (3, (4)) are the atom involved 
c               (assuming gro2joyce sorted str before bnd before dih)
              ia1     = atlab(iflag(ia,1))
              ia2     = atlab(iflag(ia,2))
              ib1     = atlab(iflag(ib,1))
              ib2     = atlab(iflag(ib,2))
              if (iflag(ia,0).eq.1.and.iflag(ib,0).eq.1) then
                  icptype = 1
              elseif (iflag(ia,0).eq.1.and.iflag(ib,0).eq.2) then
                  icptype = 2
                  ib3     = atlab(iflag(ib,3))
              elseif (iflag(ia,0).eq.2.and.iflag(ib,0).eq.2) then
                  icptype = 3
                  ia3     = atlab(iflag(ia,3))
                  ib3     = atlab(iflag(ib,3))
              elseif (iflag(ia,0).eq.1.and.iflag(ib,0).eq.4) then
                  icptype = 4
                  ib3     = atlab(iflag(ib,3))
                  ib4     = atlab(iflag(ib,4))
              elseif (iflag(ia,0).eq.2.and.iflag(ib,0).eq.4) then
                  icptype = 5
                  ia3     = atlab(iflag(ia,3))
                  ib3     = atlab(iflag(ib,3))
                  ib4     = atlab(iflag(ib,4))
              elseif (iflag(ia,0).eq.4.and.iflag(ib,0).eq.4) then
                  icptype = 6
                  ia3     = atlab(iflag(ia,3))
                  ia4     = atlab(iflag(ia,4))
                  ib3     = atlab(iflag(ib,3))
                  ib4     = atlab(iflag(ib,4))
              endif
              ja1     = atlab(iflag(ja,1))
              ja2     = atlab(iflag(ja,2))
              jb1     = atlab(iflag(jb,1))
              jb2     = atlab(iflag(jb,2))
              if (iflag(ja,0).eq.1.and.iflag(jb,0).eq.1) then
                  jcptype = 1
              elseif (iflag(ja,0).eq.1.and.iflag(jb,0).eq.2) then
                  jcptype = 2
                  jb3     = atlab(iflag(jb,3))
              elseif (iflag(ja,0).eq.2.and.iflag(jb,0).eq.2) then
                  jcptype = 3
                  ja3     = atlab(iflag(ja,3))
                  jb3     = atlab(iflag(jb,3))
              elseif (iflag(ja,0).eq.1.and.iflag(jb,0).eq.4) then
                  jcptype = 4
                  jb3     = atlab(iflag(jb,3))
                  jb4     = atlab(iflag(jb,4))
              elseif (iflag(ja,0).eq.2.and.iflag(jb,0).eq.4) then
                  jcptype = 5
                  ja3     = atlab(iflag(ja,3))
                  jb3     = atlab(iflag(jb,3))
                  jb4     = atlab(iflag(jb,4))
              elseif (iflag(ja,0).eq.4.and.iflag(jb,0).eq.4) then
                  jcptype = 6
                  ja3     = atlab(iflag(ja,3))
                  ja4     = atlab(iflag(ja,4))
                  jb3     = atlab(iflag(jb,3))
                  jb4     = atlab(iflag(jb,4))
              endif
              if (icptype.eq.jcptype) then
                  if(icptype.eq.1   .and. (      !!    strstr
     &               (((ia1.eq.ja1.and.ia2.eq.ja2).or.
     &                 (ia1.eq.ja2.and.ia2.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2).or.
     &                 (ib1.eq.jb2.and.ib2.eq.jb1)) )   .or.
     &               (((ia1.eq.jb1.and.ia2.eq.jb2).or.
     &                 (ia1.eq.jb2.and.ia2.eq.jb1)) .and.
     &                ((ib1.eq.ja1.and.ib2.eq.ja2).or.
     &                 (ib1.eq.ja2.and.ib2.eq.ja1)) )   
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

                  if(icptype.eq.2   .and. (      !!    strbnd
     &                ((ia1.eq.ja1.and.ia2.eq.ja2).or.
     &                 (ia1.eq.ja2.and.ia2.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2.and.ib3.eq.jb3).or.
     &                 (ib1.eq.jb3.and.ib2.eq.jb2.and.ib3.eq.jb1)) 
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

                  if(icptype.eq.3   .and. (      !!   bndbnd 
     &               (((ia1.eq.ja1.and.ia2.eq.ja2.and.ia3.eq.ja3).or.
     &                 (ia1.eq.ja3.and.ia2.eq.ja2.and.ia3.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2.and.ib3.eq.jb3).or.
     &                 (ib1.eq.jb3.and.ib2.eq.jb2.and.ib3.eq.jb1)) ).or.
     &               (((ia1.eq.jb1.and.ia2.eq.jb2.and.ia3.eq.jb3).or.
     &                 (ia1.eq.jb3.and.ia2.eq.jb2.and.ia3.eq.jb1)) .and.
     &                ((ib1.eq.ja1.and.ib2.eq.ja2.and.ib3.eq.ja3).or.
     &                 (ib1.eq.ja3.and.ib2.eq.ja2.and.ib3.eq.ja1)) )   
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

                  if(icptype.eq.4   .and. (      !!    strdih
     &                ((ia1.eq.ja1.and.ia2.eq.ja2).or.
     &                 (ia1.eq.ja2.and.ia2.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2.and.ib3.eq.jb3
     &                                           .and.ib4.eq.jb4).or.
     &                 (ib1.eq.jb4.and.ib2.eq.jb3.and.ib3.eq.jb2
     &                                           .and.ib4.eq.jb1)) 
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

                  if(icptype.eq.5   .and. (      !!    bnddih
     &                ((ia1.eq.ja1.and.ia2.eq.ja2.and.ia3.eq.ja3).or.
     &                 (ia1.eq.ja3.and.ia2.eq.ja2.and.ia3.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2.and.ib3.eq.jb3
     &                                           .and.ib4.eq.jb4).or.
     &                 (ib1.eq.jb4.and.ib2.eq.jb3.and.ib3.eq.jb2
     &                                           .and.ib4.eq.jb1)) 
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

                  if(icptype.eq.6   .and. (      !!   dihdih 
     &               (((ia1.eq.ja1.and.ia2.eq.ja2.and.ia3.eq.ja3
     &                                           .and.ia4.eq.ja4).or.
     &                 (ia1.eq.ja4.and.ia2.eq.ja3.and.ia3.eq.ja2
     &                                           .and.ia4.eq.ja1)) .and.
     &                ((ib1.eq.jb1.and.ib2.eq.jb2.and.ib3.eq.jb3
     &                                           .and.ib4.eq.jb4).or.
     &                 (ib1.eq.jb4.and.ib2.eq.jb3.and.ib3.eq.jb2
     &                                           .and.ib4.eq.jb1)) ).or.
     &               (((ia1.eq.jb1.and.ia2.eq.jb2.and.ia3.eq.jb3
     &                                           .and.ia4.eq.jb4).or.
     &                 (ia1.eq.jb4.and.ia2.eq.jb3.and.ia3.eq.jb2
     &                                           .and.ia4.eq.jb1)) .and.
     &                ((ib1.eq.ja1.and.ib2.eq.ja2.and.ib3.eq.ja3
     &                                           .and.ib4.eq.ja4).or.
     &                 (ib1.eq.ja4.and.ib2.eq.ja3.and.ib3.eq.ja2
     &                                           .and.ib4.eq.ja1)) )   
     &              )) then
                          zdep(jco)=.true.
                          write(35,106) jco,'=',ico,'*1.d0','; <',
     &                           labpar(ja),'*',labpar(jb), '> = <',
     &                           labpar(ia),'*',labpar(ib),'>'
                  endif

              endif

            endif
           endif        
5          continue
        enddo
      enddo
      
      write(35,'(a)') '$end '
      close(35)
105   format(3x,i4,1x,a1,1x,i4,a5,3x,a1,1x,a20,1x,a1,1x,a20)      
106   format(3x,i4,1x,a1,1x,i4,a5,3x,a3,a7,1x,a1,1x,a7,1x,
     &       a5,a7,1x,a1,1x,a7,1x,a1)      

      return
      end
c
