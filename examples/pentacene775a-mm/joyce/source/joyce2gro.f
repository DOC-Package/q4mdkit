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
      Subroutine Print_Gromacs (iun,iflag,par,parfun,parfun2,npar,geom
     &                         ,nat,maxpar,maxat,aname,labpar)
c----------------------------------------------------------------------c

c      writes a new.top file for Gromacs MD program
c      substituting Joyce computed force constants
c      and equilibrum values (R0, theta0) read from G09 .fchk files

      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      parameter       (boh2nm  = 0.0529177249d0)!! 1 bohr = 0.0529177 nm
      parameter       (H2kJ    = 2625.5d0     )!! 1 H = 2625.5d0 kJ/mol
      
      character        card*90,string*90
      character        card0*90,string0*90
      character*46     word(maxspl),subword(maxspl)
      character        aname(*)*5,labpar(*)*(*)

      integer          p,pchk,qchk
      real*8           k_s,k_b,k_bl,k_d,k_di 
      real*8           k_ss, k_sb, k_st
      real*8           k_bb, k_bt, k_tt
      logical          z_swap

      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension par(maxpar),parfun(maxpar),parfun2(maxpar)


      do ip = 1,Npar
         if(iflag(ip,0) .ne. 8) Npar0=ip
      enddo

C     Open new topol.top formatted file (joyce.new.top)

      open(25,file='joyce.new.top',form='formatted',status='unknown')
      rewind (25)

C     Reads input top file from iun unit
C     prints unchanged cards, substitutes force constants 
c     At the end of each gromacs section [ directive ]:
c     msplit = 0   --> blank line in .top files before next section
c     contains "[" --> next section starts without a blank line 
c     msplit > 0   --> reads stretchings definitions and substitutes 
c                      force constants and eq. values
 
      ntype = 0
      rewind iun
 10   READ(iun,1000,end=100) card

      if( index(card,'#') .eq. 1) goto 99
      if( index(card,';') .eq. 1) goto 99
      CALL SPLIT (card,word,nsplit,maxspl)
      if(nsplit .eq. 0)         goto 99
      if(word(1)(1:1) .eq. '#') goto 99
      if(word(1)(1:1) .eq. ';') goto 99

c     Begins writing : 1) bond stretching

      IF( word(2) .eq. 'bonds') THEN

         write(25,'(a)') card(1:Lenw(card))
 11      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
c        search for comments         
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 11
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 11
         endif
         if(msplit .eq. 0) then 
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) i
            read(subword(2),*) j
            read(subword(3),*) iftype
            read(subword(4),*) r0
            read(subword(5),*) k_s
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 1) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  if (ia .eq. i .and. ib .eq. j) goto 111
                  if (ia .eq. j .and. ib .eq. i) goto 111
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match atoms in stretching',
     &         ' while writing new.top'
            write(6,'(/a27,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j         
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
111         continue
            cost = par(ip)
c           conversion 1 H/(Bohr^2) = 2625.5/(0.052917)^2 kJ/(mol*nm^2)
            cost = cost*H2kJ/(boh2nm*boh2nm)
c           conversion 1 Bohr = 0.0529177249 nm
            r0   = parfun(ip)*boh2nm
            write(25,330) ia,ib,iftype,r0,cost
            goto 11
         endif


c     2) angle bending


      ELSEIF( word(2) .eq. 'angles') THEN

         write(25,'(a)') card(1:Lenw(card))
 12      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
            write(25,'(a)') string(1:Lenw(string))
            goto 12
         elseif(subword(1)(1:1) .eq. ';') then
            write(25,'(a)') string(1:Lenw(string))
            goto 12
         endif
         if(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) i
            read(subword(2),*) j
            read(subword(3),*) k
            read(subword(4),*) iftype
            if (iftype.eq.1) then           
               read(subword(5),*) t0_b
               read(subword(6),*) k_b
               do ip = 1,Npar
                  if(iflag(ip,0) .eq. 2) then  
                     ia = iflag(ip,1)
                     ib = iflag(ip,2)
                     ic = iflag(ip,3)
                     if(ia.eq.i.and.ib.eq.j.and.ic.eq.k)goto 113
                     if(ia.eq.k.and.ib.eq.j.and.ic.eq.i)goto 113
                  endif
               enddo
               write(6,'(/,a)')  '  **********************************'
               write(6,'(/a37,1x,a25)') 
     &            'ERROR: cannot match atoms in bending',
     &            ' while writing new.system'
               write(6,'(/a27,1x,i3,1x,i3,1x,i3)') 
     &            'Atoms in input system file:',i,j,k         
               write(6,'(/,a)')  '  **********************************'
               stop 'output error'
113            continue
               cost = par(ip)
c              conversion 1 H/rad^2 = 2625.5d0 kJ/mol/rad^2
               cost = cost*H2kJ
               t0   = parfun(ip)*rad2deg
               write(25,331) ia,ib,ic,iftype,t0,cost
               goto 12
c           OBSOLETE to be deleted, Gromacs standard str-str couplings
c           elseif (iftype .eq. 3) then !! STR - STR
c              read(subword(5),*) krr_ss
c              Call Gia_ic (Npar,iflag,i,j,0,0,i_cpld1,maxpar)
c              Call Gia_ic (Npar,iflag,j,k,0,0,i_cpld2,maxpar)
c              do ip = 1,Npar
c                 if(iflag(ip,0) .eq. 8) then
c                    ia = iflag(ip,1)
c                    ib = iflag(ip,2)
c                    if(ia .eq. i_cpld1 .and. ib .eq. i_cpld2) goto 114
c                    if(ia .eq. i_cpld2 .and. ib .eq. i_cpld1) goto 114
c                 endif
c              enddo
c              write(6,'(/,a)')  '  **********************************'
c              write(6,'(/a37,1x,a25)') 
c    &            'ERROR: cannot match pairs in coupling',
c    &            ' while writing new.system'
c              write(6,'(/a27,1x,i3,1x,i3,1x,i3)') 
c    &            'Atoms in input system file:',i,j,k         
c              write(6,'(/,a)')  '  **********************************'
c              stop 'output error'
c114           continue
c              cost = par(ip)
c              conversion 1 H/(Bohr^2) = 10^6/1.066 kJ/(mol*nm^2)
c              cost = cost*1.d6/1.066
c              conversion 1 Bohr = 0.529177249/10 nm
c              r1e  = parfun(i_cpld1)*0.529177249d0/10.d0
c              r2e  = parfun(i_cpld2)*0.529177249d0/10.d0
c              write(25,303) i,j,k,iftype,r1e,r2e,cost
c              goto 12
            endif
         endif

c     3) dihedral torsions (proper and improper)


      ELSEIF( word(2) .eq. 'dihedrals') THEN
         
         iold = 0
         jold = 0
         kold = 0
         lold = 0
         write(25,'(a)') card(1:Lenw(card))

 13      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
          write(25,'(a)') string(1:Lenw(string))
          goto 13
         elseif(subword(1)(1:1) .eq. ';') then
          write(25,'(a)') string(1:Lenw(string))
          goto 13
         endif
         if(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
          call flush(6)
          read(subword(1),*) i
          read(subword(2),*) j
          read(subword(3),*) k
          read(subword(4),*) l
          read(subword(5),*) iftype
          if(iftype .eq. 1) then   !! proper dihedrals (with cosine series)
            if (iold .eq. i .and. jold .eq. j .and. 
     &          kold .eq. k .and. lold .eq. l) then
                iterm = iterm + 1
            else    
                iterm = 0
            endif
            do ip = 1,Npar  
             if(iflag(ip,0) .eq. 4) then
                ia   = iflag(ip,1)
                ib   = iflag(ip,2)
                ic   = iflag(ip,3)
                id   = iflag(ip,4)
                if (ia .eq. i .and. ib .eq. j .and. 
     &              ic .eq. k .and. id .eq. l) then
                       iold = i
                       jold = j
                       kold = k
                       lold = l
                       goto 1121
                endif
                if (ia .eq. l .and. ib .eq. k .and. 
     &              ic .eq. j .and. id .eq. i) then
                       iold = i
                       jold = j
                       kold = k
                       lold = l
                       goto 1121
                endif
             endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a45,1x,a24)') 
     &         'ERROR: cannot match atoms in proper torsion',
     &         ' while writing new.system'
            write(6,'(a,4i5)') 'atoms i,j,k,l =', i,j,k,l
            write(6,'(a)') 'line: '//string(1:Lenw(string))
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1121        continue
            cost = par(ip+iterm)
c           conversion 1 H = 2625.5d0 kJ/mol
            cost = cost*H2kJ
            nk   = iflag(ip+iterm,5)
            gk   = parfun(ip+iterm)*rad2deg
c           GP Feb. 2023 ! deleted comment for gromacs (finally!)
c           if (nk .eq. 0 ) then ! comment for gromacs
c             write(25,301) ';',ia,ib,ic,id,iftype,gk,cost,nk
c           else
c             write(25,308)
c    &         ia,ib,ic,id,iftype,gk,cost,nk,(subword(kk),kk=9,msplit)
c           endif
            write(25,332) ia,ib,ic,id,iftype,gk,cost,nk
            goto 13

          elseif(iftype .eq. 2) then   !! harmonic dihedrals
            read(subword(6),*) p0_di
            read(subword(7),*) k_di
            do ip = 1,Npar
             if(iflag(ip,0) .eq. 3) then
                ia = iflag(ip,1)
                ib = iflag(ip,2)
                ic = iflag(ip,3)
                id = iflag(ip,4)
                if (ia .eq. i .and. ib .eq. j .and. 
     &              ic .eq. k .and. id .eq. l) goto 115
                if (ia .eq. l .and. ib .eq. k .and. 
     &              ic .eq. j .and. id .eq. i) goto 115
             endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a45,1x,a24)') 
     &         'ERROR: cannot match atoms in improper torsion',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
               stop 'output error'
115         continue
            cost = par(ip)
c           conversion 1 H = 2625.5d0 kJ/mol
            cost = cost*H2kJ
            p0   = parfun(ip)*rad2deg
            write(25,333) ia,ib,ic,id,iftype,p0,cost
            goto 13

          endif
         endif

c     4) LJ intramolecular interactions

      ELSEIF( word(2) .eq. 'pairs') THEN

         write(25,'(a)') card(1:Lenw(card))
 15      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
          write(25,'(a)') string(1:Lenw(string))
          goto 15
         elseif(subword(1)(1:1) .eq. ';') then
          write(25,'(a)') string(1:Lenw(string))
          goto 15
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         endif
         read(subword(3),*) iftype

         if (iftype .eq. 1) then   !! LJ only
            read(subword(1),*) i
            read(subword(2),*) j
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 5) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  if (ia .eq. i .and. ib .eq. j) goto 1111
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match atoms in LJ intra',
     &         ' while writing new.system'
            write(6,'(/a27,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j         
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1111        continue
            cost = par(ip)
c           conversion 1 H = 2625.5d0 KJ/mol
            cost = cost*H2kJ
            sig  = parfun(ip)*boh2nm
            write(25,330) ia,ib,iftype,sig,cost
            goto 15
         elseif (iftype .eq. 2) then   !! LJ + Coulomb
            read(subword(1),*) i
            read(subword(2),*) j
            read(subword(4),*) fq
            read(subword(5),*) qi
            read(subword(6),*) qj
            read(subword(7),*) sig
            read(subword(8),*) eps
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 5) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  if (ia .eq. i .and. ib .eq. j) goto 1112
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match atoms in LJ intra',
     &         ' while writing new.system'
            write(6,'(/a27,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j         
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1112        continue
            cost = par(ip)
c           conversion 1 H = 2625.5d0 kJ/mol
            cost = cost*H2kJ
            sig  = parfun(ip)*boh2nm
            write(25,314) ia,ib,iftype,fq,qi,qj,sig,cost,
     &                    (subword(kk),kk=9,msplit)
            goto 15
         endif

c     5) couplings (Joyce .top format)
c        new devel 9/2017 GP+JC

      ELSEIF( word(2) .eq. 'bo_bo') THEN       ! stretching-stretching

         write(25,'(a)') card(1:Lenw(card))
1116     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1116
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1116
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) icss1
            read(subword(2),*) icss2
            read(subword(3),*) icssa
            read(subword(4),*) icssb
            read(subword(5),*) iftype
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc2,1)
                  lchk = iflag(iuc2,2)
                  if (ichk .eq. icss1 .and. jchk .eq. icss2  .and.
     &                kchk .eq. icssa .and. lchk .eq. icssb) then
                    z_swap=.false.
                    goto 1130
                  elseif (ichk .eq. icssa .and. jchk .eq. icssb  .and.
     &                    kchk .eq. icss1 .and. lchk .eq. icss2) then
                    z_swap=.true.
                    goto 1130
                  endif
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a41,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,4(i3,1x))') 
     &         'Atoms declared in input system file:',icss1,icss2,
     &         icssa,icssb
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1130        continue
            cost = par(ip)
            r1e  = parfun(ip)
            r2e  = parfun2(ip)
c           conversion 1 H/(Bohr^2) = 2625.5/(0.052917)^2 kJ/(mol*nm^2)
            cost = cost*H2kJ/(boh2nm*boh2nm)
c           conversion 1 Bohr = 0.0529177249 nm
            r1e  = r1e*boh2nm
            r2e  = r2e*boh2nm
            if (z_swap) then 
              write(25,315) icssa, icssb, icss1, icss2, iftype,
     &                      r2e,r1e,cost,(subword(kk),kk=9,msplit)
            else 
              write(25,315) icss1, icss2, icssa, icssb, iftype,
     &                      r1e,r2e,cost,(subword(kk),kk=9,msplit)
            endif
            goto 1116
         endif

      ELSEIF( word(2) .eq. 'bo_an') THEN       ! stretching-bending

         write(25,'(a)') card(1:Lenw(card))
1117     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1117
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1117
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) icsb1
            read(subword(2),*) icsb2
            read(subword(3),*) icsba
            read(subword(4),*) icsbb
            read(subword(5),*) icsbc
            read(subword(6),*) iftype
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc2,1)
                  lchk = iflag(iuc2,2)
                  mchk = iflag(iuc2,3)
                  if (ichk .eq. icsb1 .and. jchk .eq. icsb2 .and. 
     &                kchk .eq. icsba .and. lchk .eq. icsbb .and.
     &                mchk .eq. icsbc) goto 1131
                  if (ichk .eq. icsb2 .and. jchk .eq. icsb1 .and. 
     &                kchk .eq. icsba .and. lchk .eq. icsbb .and.
     &                mchk .eq. icsbc) goto 1131
                  if (ichk .eq. icsb1 .and. jchk .eq. icsb2 .and. 
     &                kchk .eq. icsbc .and. lchk .eq. icsbb .and.
     &                mchk .eq. icsba) goto 1131
                  if (ichk .eq. icsb2 .and. jchk .eq. icsb1 .and. 
     &                kchk .eq. icsbc .and. lchk .eq. icsbb .and.
     &                mchk .eq. icsba) goto 1131
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a40,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,5(i3,1x))') 
     &         'Atoms declared in input system file:',icsb1,icsb2,
     &         icsba,icsbb,icsbc
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1131        continue
            cost = par(ip)
            r1e  = parfun(ip)
            t2e  = parfun2(ip)
c           conversion 1 H/(Bohr) = 2625.5/0.052917 kJ/(mol*nm)
            cost = cost*H2kJ/boh2nm
c           conversion 1 Bohr = 0.529177249/10 nm
            r1e  = r1e*boh2nm
            t2e  = t2e*rad2deg
            write(25,316) icsb1, icsb2, icsba, icsbb, icsbc, iftype,
     &                      r1e,t2e,cost,(subword(kk),kk=10,msplit)
            goto 1117
         endif

      ELSEIF( word(2) .eq. 'bo_di') THEN       ! stretching-torsion

         write(25,'(a)') card(1:Lenw(card))
1118     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1118
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1118
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) icst1
            read(subword(2),*) icst2
            read(subword(3),*) icsta
            read(subword(4),*) icstb
            read(subword(5),*) icstc
            read(subword(6),*) icstd
            read(subword(7),*) iftype
            if (iftype .eq. 1) then
               read(subword(11),*) ntt
            endif
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc2,1)
                  lchk = iflag(iuc2,2)
                  mchk = iflag(iuc2,3)
                  nchk = iflag(iuc2,4)
                  if (iftype .eq. 1) then
                     nc  = iflag(ip,5)
                     if (nc .ne. ntt) cycle
                  endif
                  if (ichk .eq. icst1 .and. jchk .eq. icst2 .and. 
     &                kchk .eq. icsta .and. lchk .eq. icstb .and.
     &                mchk .eq. icstc .and. nchk .eq. icstd) goto 1132
                  if (ichk .eq. icst2 .and. jchk .eq. icst1 .and. 
     &                kchk .eq. icsta .and. lchk .eq. icstb .and.
     &                mchk .eq. icstc .and. nchk .eq. icstd) goto 1132
                  if (ichk .eq. icst1 .and. jchk .eq. icst2 .and. 
     &                kchk .eq. icstd .and. lchk .eq. icstc .and.
     &                mchk .eq. icstb .and. nchk .eq. icsta) goto 1132
                  if (ichk .eq. icst2 .and. jchk .eq. icst2 .and. 
     &                kchk .eq. icstd .and. lchk .eq. icstc .and.
     &                mchk .eq. icstb .and. nchk .eq. icsta) goto 1132
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a41,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,6(i3,1x))') 
     &         'Atoms declared in input system file:',icst1,icst2,
     &         icsta,icstb,icstc,icstd
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1132        continue
            if (iftype .eq. 1) then   ! one proper torsions
               cost = par(ip)
               r1e  = parfun(ip)
               gam  = parfun2(ip)
               nc   = iflag(ip,5)
c              conversion 1 H/(Bohr) = 2625.5/0.052917 kJ/(mol*nm)
               cost = cost*H2kJ/boh2nm
c              conversion 1 Bohr = 0.529177249/10 nm
               r1e  = r1e*boh2nm
               gam  = gam*rad2deg
               write(25,317) icst1, icst2, icsta, icstb, icstc,icstd, 
     &                        iftype,r1e,gam,cost,nc,
     &                        (subword(kk),kk=12,msplit)
               goto 1118
            elseif (iftype .eq. 2) then   ! one harmonic torsions
               cost = par(ip)
               r1e  = parfun(ip)
               p1e  = parfun2(ip)
c              conversion 1 H/(Bohr) = 2625.5/0.052917 kJ/(mol*nm)
               cost = cost*H2kJ/boh2nm
c              conversion 1 Bohr = 0.529177249/10 nm
               r1e  = r1e*boh2nm
               p1e  = p1e*rad2deg
               write(25,318) icst1, icst2, icsta, icstb, icstc,icstd, 
     &                        iftype,r1e,p1e,cost,
     &                        (subword(kk),kk=11,msplit)
               goto 1118
            endif
         endif


      ELSEIF( word(2) .eq. 'an_an') THEN       ! bending-bending

         write(25,'(a)') card(1:Lenw(card))
1119     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1119
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1119
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) icbb1
            read(subword(2),*) icbb2
            read(subword(3),*) icbb3
            read(subword(4),*) icbba
            read(subword(5),*) icbbb
            read(subword(6),*) icbbc
            read(subword(7),*) iftype
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc1,3)
                  lchk = iflag(iuc2,1)
                  mchk = iflag(iuc2,2)
                  nchk = iflag(iuc2,3)
                  if (ichk .eq. icbb1 .and. jchk .eq. icbb2 .and. 
     &                kchk .eq. icbb3 .and. lchk .eq. icbba .and.
     &                mchk .eq. icbbb .and. nchk .eq. icbbc) goto 1133
                  if (ichk .eq. icbb3 .and. jchk .eq. icbb2 .and. 
     &                kchk .eq. icbb1 .and. lchk .eq. icbba .and.
     &                mchk .eq. icbbb .and. nchk .eq. icbbc) goto 1133
                  if (ichk .eq. icbb1 .and. jchk .eq. icbb2 .and. 
     &                kchk .eq. icbb3 .and. lchk .eq. icbbc .and.
     &                mchk .eq. icbbb .and. nchk .eq. icbba) goto 1133
                  if (ichk .eq. icbb3 .and. jchk .eq. icbb2 .and. 
     &                kchk .eq. icbb1 .and. lchk .eq. icbbc .and.
     &                mchk .eq. icbbb .and. nchk .eq. icbba) goto 1133
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a41,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,6(i3,1x))') 
     &         'Atoms declared in input system file:',icbb1,icbb2,
     &         icbb3,icbba,icbbb,icbbc
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1133        continue
            cost = par(ip)
            t1e  = parfun(ip)
            t2e  = parfun2(ip)
c           conversion 1 H = 2625.5 kJ/mol
            cost = cost*H2kJ
            t1e  = r1e*rad2deg
            t2e  = t2e*rad2deg
            write(25,319) icbb1, icbb2, icbb3, icbba, icsbb, icsbc, 
     &                    iftype,r1e,t2e,cost,(subword(kk),kk=11,msplit)
            goto 1119
         endif


      ELSEIF( word(2) .eq. 'an_di') THEN       ! bending-torsion

         write(25,'(a)') card(1:Lenw(card))
1120     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1120
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1120
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) icbt1
            read(subword(2),*) icbt2
            read(subword(3),*) icbt3
            read(subword(4),*) icbta
            read(subword(5),*) icbtb
            read(subword(6),*) icbtc
            read(subword(7),*) icbtd
            read(subword(8),*) iftype
            if (iftype .eq. 1) then
               read(subword(12),*) ntt
            endif
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc1,3)
                  lchk = iflag(iuc2,1)
                  mchk = iflag(iuc2,2)
                  nchk = iflag(iuc2,3)
                  pchk = iflag(iuc2,4)
                  if (iftype .eq. 1) then
                     nc  = iflag(ip,5)
                     if (nc .ne. ntt) cycle
                  endif
                  if (ichk .eq. icbt1 .and. jchk .eq. icbt2 .and. 
     &                kchk .eq. icbt3 .and. lchk .eq. icbta .and.
     &                mchk .eq. icbtb .and. nchk .eq. icbtc.and.
     &                pchk .eq. icbtd) goto 1134
                  if (ichk .eq. icbt3 .and. jchk .eq. icbt2 .and. 
     &                kchk .eq. icbt1 .and. lchk .eq. icbta .and.
     &                mchk .eq. icbtb .and. nchk .eq. icbtc.and.
     &                pchk .eq. icbtd) goto 1134
                  if (ichk .eq. icbt1 .and. jchk .eq. icbt2 .and. 
     &                kchk .eq. icbt3 .and. lchk .eq. icbtd .and.
     &                mchk .eq. icbtc .and. nchk .eq. icbtb.and.
     &                pchk .eq. icbta) goto 1134
                  if (ichk .eq. icbt3 .and. jchk .eq. icbt2 .and. 
     &                kchk .eq. icbt1 .and. lchk .eq. icbtd .and.
     &                mchk .eq. icbtc .and. nchk .eq. icbtb.and.
     &                pchk .eq. icbta) goto 1134
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,7(i3,1x))') 
     &         'Atoms declared in input system file:',icbt1,icbt2,
     &         icbt3,icbta,icbtb,icbtc,icbtd
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1134        continue
            if (iftype .eq. 1) then   ! one proper torsions
               cost = par(ip)
               t1e  = parfun(ip)
               gam  = parfun2(ip)
               nc   = iflag(ip,5)
c              conversion 1 H = 2625.5 kJ/mol
               cost = cost*H2kJ
               t1e  = t1e*rad2deg
               gam  = gam*rad2deg
               write(25,320) icbt1, icbt2, icbt3, 
     &                       icbta, icbtb, icbtc,icbtd, 
     &                       iftype,t1e,gam,cost,nc,
     &                        (subword(kk),kk=13,msplit)
               goto 1120
            elseif (iftype .eq. 2) then   ! one harmonic torsions
               cost = par(ip)
               t1e  = parfun(ip)
               p1e  = parfun2(ip)
c              conversion 1 H = 2625.5 kJ/mol
               cost = cost*H2kJ
               t1e  = t1e*boh2nm
               p1e  = p1e*rad2deg
               write(25,321) icbt1, icbt2, icbt3, 
     &                       icbta, icbtb, icbtc,icbtd, 
     &                       iftype,t1e,p1e,cost,
     &                        (subword(kk),kk=12,msplit)
               goto 1120
            endif
         endif

      ELSEIF( word(2) .eq. 'di_di') THEN       ! torsion-torsion

         write(25,'(a)') card(1:Lenw(card))
1122     read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1122
         elseif(subword(1)(1:1) .eq. ';') then
             write(25,'(a)') string(1:Lenw(string))
             goto 1122
         elseif(msplit .eq. 0) then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(index(string,'[') .ne. 0)   then  
            backspace(iun)
            goto 10 
         else
            read(subword(1),*) ictt1
            read(subword(2),*) ictt2
            read(subword(3),*) ictt3
            read(subword(4),*) ictt4
            read(subword(5),*) ictta
            read(subword(6),*) icttb
            read(subword(7),*) icttc
            read(subword(8),*) icttd
            read(subword(9),*) iftype
            if (iftype .eq. 1) then
               read(subword(13),*)ntt1
               read(subword(14),*)ntt2
            elseif (iftype .eq. 2) then
               read(subword(13),*)ntt1
            endif
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 8) then
                  iuc1 = iflag(ip,1)
                  iuc2 = iflag(ip,2)
                  ichk = iflag(iuc1,1)
                  jchk = iflag(iuc1,2)
                  kchk = iflag(iuc1,3)
                  lchk = iflag(iuc1,4)
                  mchk = iflag(iuc2,1)
                  nchk = iflag(iuc2,2)
                  pchk = iflag(iuc2,3)
                  qchk = iflag(iuc2,4)
                  if (iftype .eq. 1) then
                     nc1  = iflag(ip,5)
                     nc2  = iflag(ip,6)
                     if (nc1 .ne. ntt1 .or. nc2 .ne. ntt2) cycle
                  elseif (iftype .eq. 2) then
                     nc1  = iflag(ip,5)
                     if (nc1 .ne. ntt1) cycle
                  endif
                  if (ichk .eq. ictt1 .and. jchk .eq. ictt2 .and. 
     &                kchk .eq. ictt3 .and. lchk .eq. ictt4 .and.
     &                mchk .eq. ictta .and. nchk .eq. icttb.and.
     &                pchk .eq. icttc .and. qchk .eq. icttd) goto 1135
                  if (ichk .eq. ictt4 .and. jchk .eq. ictt3 .and. 
     &                kchk .eq. ictt2 .and. lchk .eq. ictt1 .and.
     &                mchk .eq. ictta .and. nchk .eq. icttb.and.
     &                pchk .eq. icttc .and. qchk .eq. icttd) goto 1135
                  if (ichk .eq. ictt1 .and. jchk .eq. ictt2 .and. 
     &                kchk .eq. ictt3 .and. lchk .eq. ictt4 .and.
     &                mchk .eq. icttd .and. nchk .eq. icttc.and.
     &                pchk .eq. icttb .and. qchk .eq. ictta) goto 1135
                  if (ichk .eq. ictt4 .and. jchk .eq. ictt3 .and. 
     &                kchk .eq. ictt2 .and. lchk .eq. ictt1 .and.
     &                mchk .eq. icttd .and. nchk .eq. icttc.and.
     &                pchk .eq. icttb .and. qchk .eq. ictta) goto 1135
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match IC atoms in couplings',
     &         ' while writing new.top'
            write(6,'(/a27,1x,8(i3,1x))') 
     &         'Atoms declared in input system file:',ictt1,ictt2,
     &         ictt3,ictt4,ictta,icttb,icttc,icttd
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1135        continue
            if (iftype .eq. 1) then   ! two proper torsions
               cost = par(ip)
               gam1 = parfun(ip)
               gam2 = parfun2(ip)
               nc1  = iflag(ip,5)
               nc2  = iflag(ip,6)
c              conversion 1 H = 2625.5 kJ/mol
               cost = cost*H2kJ
               gam1 = gam1*rad2deg
               gam2 = gam2*rad2deg
               write(25,322) ictt1, ictt2, ictt3, ictt4, 
     &                       ictta, icttb, icttc, icttd, 
     &                       iftype,gam1,gam2,cost,nc1,nc2,
     &                        (subword(kk),kk=15,msplit)
               goto 1122
            elseif (iftype .eq. 2) then   ! one harmonic and one proper
               cost = par(ip)
               p1e  = parfun(ip)
               gam  = parfun2(ip)
c              conversion 1 H = 2625.5 kJ/mol
               cost = cost*H2kJ
               p1e  = t1e*boh2nm
               gam  = p1e*rad2deg
               nc   = iflag(ip,5)
               write(25,323) ictt1, ictt2, ictt3, ictt4, 
     &                       ictta, icttb, icttc, icttd, 
     &                       iftype,p1e,gam,cost,nc,
     &                        (subword(kk),kk=14,msplit)
               goto 1122
            elseif (iftype .eq. 3) then   ! two harmonic
               cost = par(ip)
               p1e  = parfun(ip)
               p2e  = parfun2(ip)
c              conversion 1 H = 2625.5 kJ/mol
               cost = cost*H2kJ
               p1e  = t1e*boh2nm
               p2e  = p1e*rad2deg
               write(25,324) ictt1, ictt2, ictt3, ictt4, 
     &                       ictta, icttb, icttc, icttd, 
     &                       iftype,p1e,p2e,cost,
     &                        (subword(kk),kk=13,msplit)
               goto 1122
            endif
         endif
      ENDIF

c     Couplings formats
315   format(1x,i3,1x,i3,3x,i3,1x,i3,3x,i1,3x,2(f8.4,3x),2x,f9.3,2x,
     &       a1,2x,a3,2x,a5,1x,a1,1x,a5)
316   format(1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,3x,i1,3x,f8.4,3x,f7.1,2x,
     &       2x,f9.3,2x,a1,2x,a3,2x,a5,1x,a1,1x,a8)
317   format(1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f8.4,2x,f7.1,2x,f9.3,3x,i2,2x,
     &       a1,2x,a3,2x,a5,1x,a1,1x,a12)
318   format(1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f8.4,2x,f7.1,2x,f9.3,2x,
     &       a1,2x,a3,2x,a5,1x,a1,1x,a12)
319   format(1x,i3,1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,2x,
     &       a1,2x,a3,2x,a8,1x,a1,1x,a8)
320   format(1x,i3,1x,i3,1x,i3,2x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,3x,i1,2x,
     &       a1,2x,a3,2x,a8,1x,a1,1x,a12)
321   format(1x,i3,1x,i3,1x,i3,2x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,3x,i1,3x,i1,2x,
     &       a1,2x,a3,2x,a8,1x,a1,1x,a12)
322   format(1x,i3,1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,3x,i1,3x,i1,2x,
     &       a1,2x,a3,2x,a12,1x,a1,1x,a12)
323   format(1x,i3,1x,i3,1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,3x,i1,2x,
     &       a1,2x,a3,2x,a12,1x,a1,1x,a12)
324   format(1x,i3,1x,i3,1x,i3,1x,i3,3x,i3,1x,i3,1x,i3,1x,i3,
     &       3x,i1,3x,f7.1,2x,f7.1,2x,f9.3,2x,
     &       a1,2x,a3,2x,a12,1x,a1,1x,a12)

 99   write(25,'(a)') card(1:Lenw(card))
      goto 10


100   continue

      write(6,'(/,3x,30("**"))')  
      write(6,'(/,3x,a)')  
     &'Gromacs .top file correctly written in file: joyce.new.top'
      write(6,'(/,3x,30("**"))')  



c 300   format(4x,i3,3x,i3,3x,i1,3x,f8.4,3x,f12.3,5x,3(a7,2x))
c 302   format(4x,i3,2x,i3,2x,i3,3x,i1,4x,f6.2,3x,f10.4,5x,2(a7,2x),a10)
300   format(4x,i3,3x,i3,3x,i1,3x,f8.4,3x,f12.3,3x,a1,2x,a4,2x,a20)
302   format(4x,i3,2x,i3,2x,i3,3x,i1,4x,f6.2,3x,f10.4,3x,
     &       a1,2x,a4,2x,a20) 
303   format(4x,i3,2x,i3,2x,i3,3x,i1,4x,f8.4,4x,f8.4,3x,f10.4)
304   format(4x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f6.1,2x,f9.3,5x,
     &       2(a4,2x),2x,a20,1x,a1)
305   format(a10,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,1x,i2,3x,f6.2,2x,f10.4)
306   format(a13,1x,i2,1x,i2,2x,f6.2,
     &       3x,i2,1x,i2,1x,i2,1x,i2,3x,i2,3x,f8.2,2x,f10.4)
307   format(a11,1x,i2,1x,i2,1x,i2,2x,f6.2,
     &       3x,i2,1x,i2,1x,i2,1x,i2,3x,i2,3x,f8.2,2x,f10.4)
308   format(4x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f7.2,3x,f12.3,
     &       2x,i2,5x,2(a4,2x),2x,a20)
330   format(4x,i3,3x,i3,3x,i1,3x,f8.4,3x,f12.3)
331   format(4x,i3,2x,i3,2x,i3,3x,i1,4x,f6.2,3x,f10.4)
332   format(4x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f7.2,3x,f12.3,2x,i2)
333   format(4x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f6.1,2x,f9.3)
301   format(a1,4x,i3,1x,i3,1x,i3,1x,i3,3x,i1,3x,f7.2,3x,f12.3,2x,i2)
309   format(a9,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,3x,f6.2,2x,f10.4)
310   format(a11,2x,i2,1x,i2,1x,i2,2x,f7.3,
     &       3x,i2,1x,i2,1x,i2,3x,f6.2,2x,f10.4)
311   format(a13,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,1x,i2,1x,i2,
     &       3x,f10.4,2x,i2,3x,f8.2)
312   format(i3,3x,i3,3x,i1,3x,f10.4,3x,f13.7,5x,3(a7,2x))
314   format(i3,3x,i3,3x,i1,3x,f6.3,2x,f8.3,2x,f8.3,2x,f10.4,3x,f13.7,
     &       5x,a1,2x,3(a3,2x))
313   format(a11,1x,4(i2,1x),1x,i2,2x,f6.2,
     &       3x,4(i2,1x),1x,i2,2x,f6.2,2x,f10.4)
1000  format(a)

      return
      end
