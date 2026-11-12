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
      Subroutine Print_Moscito (iun,iflag,par,parfun,npar,geom
     &                         ,nat,maxpar,maxat,aname,labpar)
c----------------------------------------------------------------------c
 
c     writes a new.system file for Moscito MD program
c     substituting Joyce computed force constants
c     and equilibrum values (R0, theta0) read from G03 .fchk files

      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      
      integer          p
      real*8           k_s,k_b,k_bl,k_d,k_di,k_sb,k_bd
      real*8           k_ss, k_bb, k_sd
      character        card*90,string*90
      character*46     word(maxspl),subword(maxspl)
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
 
      character        aname(*)*5,labpar(*)*(*)


      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension par(maxpar),parfun(maxpar)


      do ip = 1,Npar
         if(iflag(ip,0) .ne. 8) Npar0=ip
      enddo

c      opens a new moscito.system formatted file (new.system)

      open(25,file='joyce.new.system',form='formatted',status='unknown')
      rewind (25)

c     Reads old moscito.system files from iun unit 
c     prints unchanged cards, substitutes force constants  
 
      ntype = 0
      rewind iun
 10   READ(iun,1000) card

      if( index(card,'#') .eq. 1) goto 99
      if( index(card,'end{molecule}') .gt. 0) then
             write(25,'(a)') card(1:Lenw(card))
         goto 100
      endif
      CALL SPLIT (card,word,nsplit,maxspl)
      if( nsplit .eq. 0)        goto 99
      if(word(1)(1:1) .eq. '#') goto 99

c     Starts writing : 1) bond stretching

      IF( word(1) .eq. 'begin{bonds}') THEN

         write(25,'(a)') card(1:Lenw(card))
 11      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 11
         endif
         if(subword(1) .eq. 'end{bonds}') then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(subword(1) .eq. 'harmonic') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k_s
            read(subword(5),*) r0
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 1) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  if (ia .eq. i .and. ib .eq. j) goto 111
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a39,1x,a24)') 
     &         'ERROR: cannot match atoms in stretching',
     &         ' while writing new.system'
            write(6,'(/a27,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j         
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
111         continue
            cost = par(ip)
c           conversion: 1 H/(Bohr^2) = 10^4/1.066 kJ/(mol*A^2)
            cost = cost*1.d4/1.066
c           conversion: 1 Bohr = 0.529177249 A
            r0   = parfun(ip)*0.529177249d0
            write(25,300) 'harmonic',ia,ib,cost,r0
            goto 11
         else
            goto 11
         endif

c     2) angle bending (linear and quadratic)
c              Linear
c        V(t) =  k_bl [1 + cos(t)]
c        dV^2/dt^2 = -k_bl*cos(t0)
c        k_bl      = -dV^2/dt^2/cos(t0)
c             Quadratic
c        V(t) =  k_b [t-t0]^2
c        k_b  =  dV^2/dt^2


      ELSEIF( word(1) .eq. 'begin{angle}') THEN

         write(25,'(a)') card(1:Lenw(card))
 12      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') then
            write(25,'(a)') string(1:Lenw(string))
            goto 12
         endif
         if(subword(1) .eq. 'end{angle}') then
            write(25,'(a)') string(1:Lenw(string))
            goto 10
         elseif(subword(1) .eq. 'linear') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k
            read(subword(5),*) k_bl
            write(6,'(/a)') 
     &      'Warning: linear angle was transformed in harmonic'
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 2) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  ic = iflag(ip,3)
                  if(ia .eq. i .and. ib .eq. j .and. ic .eq. k) goto 112
                  if(ia .eq. k .and. ib .eq. j .and. ic .eq. i) goto 112
               endif
            enddo
            write(6,'(/a37,1x,a25)') 
     &         'ERROR: cannot match atoms in bending',
     &         ' while writing new.system'
            write(6,'(/a27,1x,i3,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j,k         
            stop 'output error'
112         continue
            cost = par(ip)
c           conversion: 1 H/rad^2 = 2625.5d0 kJ/mol/rad^2
            cost = cost*2625.5d0
            t0   = parfun(ip)*rad2deg
            write(25,302) ia,ib,ic,cost,t0
            goto 12
         else
            read(subword(1),*) i
            read(subword(2),*) j
            read(subword(3),*) k
            read(subword(4),*) k_b
            read(subword(5),*) t0_b
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 2) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  ic = iflag(ip,3)
                  if(ia .eq. i .and. ib .eq. j .and. ic .eq. k)goto 113
                  if(ia .eq. k .and. ib .eq. j .and. ic .eq. i)goto 113
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(/a37,1x,a25)') 
     &         'ERROR: cannot match atoms in bending',
     &         ' while writing new.system'
            write(6,'(/a27,1x,i3,1x,i3,1x,i3)') 
     &         'Atoms in input system file:',i,j,k         
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
113         continue
            cost = par(ip)
c           conversion 1 H/rad^2 = 2625.5d0 kJ/mol/rad^2
            cost = cost*2625.5d0
            t0   = parfun(ip)*rad2deg
            write(25,302) ia,ib,ic,cost,t0
            goto 12
         endif

c     3) dihedral torsions (proper and improper)


      ELSEIF( word(1) .eq. 'begin{dihedral}') THEN

         
         iold = 0
         jold = 0
         kold = 0
         lold = 0
         write(25,'(a)') card(1:Lenw(card))
         write(25,'(a)') 
     &   '# type   at1 at2 at3 at4   k_d (kJ/mol)   n    gamma' 


 13      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') then

           if(subword(3) .ne. 'at1') then
              write(25,'(a)') string(1:Lenw(string))
           endif
           goto 13

         elseif(subword(1) .eq. 'end{dihedral}') then
            write(25,'(a)') string(1:Lenw(string))
            goto 10

         elseif(subword(1) .eq. 'proper') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k
            read(subword(5),*) l
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
     &                ic .eq. k .and. id .eq. l) then
                         iold = i
                         jold = j
                         kold = k
                         lold = l
                         goto 1121
                  endif
                  if (ia .eq. l .and. ib .eq. k .and. 
     &                ic .eq. j .and. id .eq. i) then
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
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
1121        continue
            cost = par(ip+iterm)
c           conversion: 1 H = 2625.5d0 kJ/mol
            cost = cost*2625.5d0
            nk   = iflag(ip+iterm,5)
            gk   = parfun(ip+iterm)*rad2deg
            write(25,308) 'proper',ia,ib,ic,id,cost,nk,gk
            goto 13

         elseif(subword(1) .eq. 'improper') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k
            read(subword(5),*) l
            read(subword(6),*) k_di
            read(subword(7),*) p0_di
            do ip = 1,Npar
               if(iflag(ip,0) .eq. 3) then
                  ia = iflag(ip,1)
                  ib = iflag(ip,2)
                  ic = iflag(ip,3)
                  id = iflag(ip,4)
                  if (ia .eq. i .and. ib .eq. j .and. 
     &                ic .eq. k .and. id .eq. l) goto 115
                  if (ia .eq. l .and. ib .eq. k .and. 
     &                ic .eq. j .and. id .eq. i) goto 115
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
c           conversion: 1 H = 2625.5d0 kJ/mol
            cost = cost*2625.5d0
            p0   = parfun(ip)*rad2deg
            write(25,304) 'improper',ia,ib,ic,id,cost,p0
            goto 13
         else   
            goto 13
         endif

c     4) LJ intramolecular interactions

      ELSEIF( word(1) .eq. 'begin{nonbonded}') THEN

         write(25,'(a)') card(1:Lenw(card))
 15      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') then
            write(25,'(a)') string(1:Lenw(string))
            goto 15
         endif
         if (msplit .eq. 4) then   !! LJ only
            if(subword(1) .eq. '#') then
             write(25,'(a)') string(1:Lenw(string))
             goto 15
            endif
            if(subword(1) .eq. 'end{nonbonded}') then
               write(25,'(a)') string(1:Lenw(string))
               goto 10
            else
               read(subword(1),*) i
               read(subword(2),*) j
               read(subword(3),*) sig
               read(subword(4),*) eps
               do ip = 1,Npar
                  if(iflag(ip,0) .eq. 5) then
                     ia = iflag(ip,1)
                     ib = iflag(ip,2)
                     if (ia .eq. i .and. ib .eq. j) goto 1111
                  endif
               enddo
               write(6,'(/,a)')  '  **********************************'
               write(6,'(/a39,1x,a24)') 
     &            'ERROR: cannot match atoms in LJ intra',
     &            ' while writing new.system'
               write(6,'(/a27,1x,i3,1x,i3)') 
     &            'Atoms in input system file:',i,j         
               write(6,'(/,a)')  '  **********************************'
               stop 'output error'
1111           continue
               cost = par(ip)
c              conversion: 1 H = 2625.5d0/8.31451121d-3 K
               cost = cost*2625.5d0/8.31451121d-3
               sig  = parfun(ip)*0.529177249d0   
               write(25,312) ia,ib,sig,cost
               goto 15
            endif
         elseif (msplit .ge. 5) then   !! LJ + Coulomb
            if(subword(5) .eq. '#') then  !! LJ only + comment
               read(subword(1),*) i
               read(subword(2),*) j
               read(subword(3),*) sig
               read(subword(4),*) eps
               do ip = 1,Npar
                  if(iflag(ip,0) .eq. 5) then
                     ia = iflag(ip,1)
                     ib = iflag(ip,2)
                     if (ia .eq. i .and. ib .eq. j) goto 1114
                  endif
               enddo
               write(6,'(/,a)')  '  **********************************'
               write(6,'(/a39,1x,a24)') 
     &            'ERROR: cannot match atoms in LJ intra',
     &            ' while writing new.system'
               write(6,'(/a27,1x,i3,1x,i3)') 
     &            'Atoms in input system file:',i,j         
               write(6,'(/,a)')  '  **********************************'
               stop 'output error'
1114           continue
               cost = par(ip)
c              conversion: 1 H = 2625.5d0/8.31451121d-3 K
               cost = cost*2625.5d0/8.31451121d-3
               sig  = parfun(ip)*0.529177249d0    
               write(25,312) ia,ib,sig,cost
               goto 15
            endif
            if(subword(1) .eq. '#') then
               write(25,'(a)') string(1:Lenw(string))
               goto 15
            endif
            if(subword(1) .eq. 'end{nonbonded}') then
               write(25,'(a)') string(1:Lenw(string))
               goto 10
            else
               read(subword(1),*) i
               read(subword(2),*) j
               read(subword(3),*) sig
               read(subword(4),*) eps
               read(subword(5),*) fqq
               do ip = 1,Npar
                  if(iflag(ip,0) .eq. 5) then
                     ia = iflag(ip,1)
                     ib = iflag(ip,2)
                     if (ia .eq. i .and. ib .eq. j) goto 1112
                  endif
               enddo
               write(6,'(/,a)')  '  **********************************'
               write(6,'(/a39,1x,a24)') 
     &            'ERROR: cannot match atoms in LJ intra',
     &            ' while writing new.system'
               write(6,'(/a27,1x,i3,1x,i3)') 
     &            'Atoms in input system file:',i,j         
               write(6,'(/,a)')  '  **********************************'
               stop 'output error'
1112           continue
               cost = par(ip)
c              conversion: 1 H = 2625.5d0/8.31451121d-3 K
               cost = cost*2625.5d0/8.31451121d-3
               sig  = parfun(ip)*0.529177249d0    
               write(25,312) ia,ib,sig,cost,fqq
               goto 15
            endif
         endif


c     5) coupled internal coordinates (cross terms)

      ELSEIF( word(1) .eq. 'begin{couplings}') THEN
   
         write(25,'(a)') card(1:Lenw(card))
 14      read(iun,1000) string
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') then
          write(25,'(a)') string(1:Lenw(string))
          goto 14
         endif
         if(subword(1) .eq. 'end{couplings}') then
          write(25,'(a)') string(1:Lenw(string))
          goto 10

         elseif(subword(1) .eq. 'bond_bond') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) r01_ss
            read(subword(5),*) k
            read(subword(6),*) l
            read(subword(7),*) r02_ss
            read(subword(8),*) k_ss
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 1 .and. ity2 .eq. 1 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ia2 .eq. k .and. ib2 .eq. l) goto 116
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a45,1x,a24)') 
     &         'ERROR: cannot match atoms in bond-bond coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
116         continue
            cost = par(ipc)
c           conversion: 1 H/(Bohr^2) = 10^4/1.066 kJ/(mol*A^2)
            cost = cost*1.d4/1.066
            r01  = parfun(icor1)*boh2ang
            r02  = parfun(icor2)*boh2ang
            write(25,309) 'bond_bond',ia1,ib1,r01,
     &                     ia2,ib2,r02,cost
            goto 14

         elseif(subword(1) .eq. 'bond_angle') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) r0_sb
            read(subword(5),*) k
            read(subword(6),*) l
            read(subword(7),*) m
            read(subword(8),*) t0_sb
            read(subword(9),*) k_sb
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 1 .and. ity2 .eq. 2 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     ic2= iflag(icor2,3)
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ia2 .eq. k .and. ib2 .eq. l .and. 
     &                   ic2 .eq. m) goto 117
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a46,1x,a24)') 
     &         'ERROR: cannot match atoms in bond-angle coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
117         continue
            cost = par(ipc)
c           conversion: 1 H/(Bohr*rad) = 10^4/2.01553 kJ/(mol*A*rad)
            cost = cost*1.d4/2.01553d0
            r01  = parfun(icor1)*boh2ang
            t02  = parfun(icor2)*rad2deg
            write(25,305) 'bond_angle',ia1,ib2,r0,
     &                     ia2,ib2,ic2,t02,cost
              goto 14


         elseif(subword(1) .eq. 'angle_angle') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k
            read(subword(5),*) t01_bb
            read(subword(6),*) l
            read(subword(7),*) m
            read(subword(8),*) n
            read(subword(9),*) t0_bb
            read(subword(10),*)k_bb
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 2 .and. ity2 .eq. 2 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ic1= iflag(icor1,3)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     ic2= iflag(icor2,3)
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ic1 .eq. k .and. 
     &                   ia2 .eq. l .and. ib2 .eq. m .and. 
     &                   ic2 .eq. n) goto 118
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a47,1x,a24)') 
     &         'ERROR: cannot match atoms in angle-angle coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
118         continue
            cost = par(ipc)
c           conversion: 1 H/(rad^2) = 2625.5 kJ/(mol*rad^2)
            cost = cost*2625.5d0
            t01  = parfun(icor1)*rad2deg
            t02  = parfun(icor2)*rad2deg
            write(25,310) 'angle_angle',ia1,ib1,ic1,t01,
     &                     ia2,ib2,ic2,t02,cost
              goto 14

         elseif(subword(1) .eq. 'bond_dihed') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) r0_sd
            read(subword(5),*) k
            read(subword(6),*) l
            read(subword(7),*) m
            read(subword(8),*) n
            read(subword(9),*) n0_sd
            read(subword(10),*) gamma_sd
            read(subword(10),*) k_sd
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 1 .and. ity2 .eq. 4 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     ic2= iflag(icor2,3)
                     id2= iflag(icor2,4)
                     nn = iflag(icor2,5)
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ia2 .eq. k .and. ib2 .eq. l .and. 
     &                   ic2 .eq. m .and. id2 .eq. n .and.
     &                   nn  .eq. n0_sd   ) goto 119
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   id2 .eq. k .and. ic2 .eq. l .and. 
     &                   ib2 .eq. m .and. ia2 .eq. n .and.
     &                   nn  .eq. n0_sd   ) goto 119
                     if( ib1 .eq. i .and. ia1 .eq. j .and. 
     &                   ia2 .eq. k .and. ib2 .eq. l .and. 
     &                   ic2 .eq. m .and. id2 .eq. n .and.
     &                   nn  .eq. n0_sd   ) goto 119
                     if( ib1 .eq. i .and. ia1 .eq. j .and. 
     &                   id2 .eq. k .and. ic2 .eq. l .and. 
     &                   ib2 .eq. m .and. ia2 .eq. n .and.
     &                   nn  .eq. n0_sd   ) goto 119
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a49,1x,a24)') 
     &         'ERROR: cannot match atoms in bond-dihedral coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            stop 'output error'
119         continue
            cost = par(ipc)
c           conversion: 1 H/Bohr = 10^4/2.01553 kJ/(mol*A)
            cost = cost*1.d4/2.01553d0
            r01  = parfun(icor1)*boh2ang
            nk   = iflag(icor2,5)
            gk   = parfun(icor2)*rad2deg
            write(25,306) 'bond_dihedral',ia1,ib2,r0,
     &                     ia2,ib2,ic2,id2,nk,gk,cost
            goto 14

         elseif(subword(1) .eq. 'angle_dihed') then
            read(subword(2),*) i
            read(subword(3),*) j
            read(subword(4),*) k
            read(subword(5),*) t0_bd
            read(subword(6),*) l
            read(subword(7),*) m
            read(subword(8),*) n
            read(subword(9),*) p
            read(subword(10),*) n0_bd
            read(subword(11),*) gamma_bd
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 2 .and. ity2 .eq. 4 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ic1= iflag(icor1,3)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     ic2= iflag(icor2,3)
                     id2= iflag(icor2,4)
                     nn = iflag(icor2,5)
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ic1 .eq. k .and. 
     &                   ia2 .eq. l .and. ib2 .eq. m .and. 
     &                   ic2 .eq. n .and. id2 .eq. p .and.
     &                   nn  .eq. n0_bd) goto 120
                     if( ia1 .eq. k .and. ib1 .eq. j .and. 
     &                   ic1 .eq. i .and. 
     &                   ia2 .eq. l .and. ib2 .eq. m .and. 
     &                   ic2 .eq. n .and. id2 .eq. p .and.
     &                   nn  .eq. n0_bd) goto 120
                     if( ia1 .eq. i .and. ib1 .eq. j .and. 
     &                   ic1 .eq. k .and. 
     &                   ia2 .eq. p .and. ib2 .eq. n .and. 
     &                   ic2 .eq. m .and. id2 .eq. l .and.
     &                   nn  .eq. n0_bd) goto 120
                     if( ia1 .eq. k .and. ib1 .eq. j .and. 
     &                   ic1 .eq. i .and. 
     &                   ia2 .eq. p .and. ib2 .eq. n .and. 
     &                   ic2 .eq. m .and. id2 .eq. l .and.
     &                   nn  .eq. n0_bd) goto 120
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a53,1x,a24)') 
     &         'ERROR: cannot match atoms in angle-dihedral coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            write(6,*)  'readen indexes :',i,j,k,l,m,n,p,'n = ',n0_bd
            write(6,*)  'stored indexes '
            write(6,*)  'bendings'
            do ipp = 1,Npar
               if (iflag(ipp,0) .eq. 4) then
                 write(6,*) ipp,(iflag(ipp,kk),kk=1,3)
               endif
            enddo
            write(6,*)  'torsions'
            do ipp = 1,Npar
               if (iflag(ipp,0) .eq. 4) then
                 write(6,*) ipp,(iflag(ipp,kk),kk=1,4),'with n =',
     &                      iflag(ipp,5)
               endif
            enddo
            stop 'output error'
120         continue
            cost = par(ipc)
c           conversion: 1 H/(rad^2) = 2625.5 kJ/(mol*rad^2)
            cost = cost*2625.5d0
            t01  = parfun(icor1)*rad2deg
            nk   = iflag(icor2,5)
            gk   = parfun(icor2)*rad2deg
              write(25,307) 'angle_dihed',ia1,ib1,ic1,t01,
     &                       ia2,ib2,ic2,id2,nk,gk,cost
              goto 14
         elseif(subword(1) .eq. 'dihed_dihed') then
            read(subword(2),*) i1
            read(subword(3),*) j1
            read(subword(4),*) k1
            read(subword(5),*) l1
            read(subword(6),*) n1
            read(subword(7),*) gamma1
            read(subword(8),*) i2
            read(subword(9),*) j2
            read(subword(10),*) k2
            read(subword(11),*) l2
            read(subword(12),*) n2
            read(subword(13),*) gamma2
            do ipc = Npar0+1,Npar
               if(iflag(ipc,0) .eq. 8) then
                  icor1 = iflag(ipc,1)
                  icor2 = iflag(ipc,2)
                  ity1= iflag(icor1,0)
                  ity2= iflag(icor2,0)
                  if (ity1 .eq. 4 .and. ity2 .eq. 4 ) then
                     ia1= iflag(icor1,1)
                     ib1= iflag(icor1,2)
                     ic1= iflag(icor1,3)
                     id1= iflag(icor1,4)
                     nn1= iflag(icor1,5)
                     ia2= iflag(icor2,1)
                     ib2= iflag(icor2,2)
                     ic2= iflag(icor2,3)
                     id2= iflag(icor2,4)
                     nn2= iflag(icor2,5)
                     if( ia1 .eq. i1 .and. ib1 .eq. j1 .and. 
     &                   ic1 .eq. k1 .and. id1 .eq. l1 .and.
     &                   ia2 .eq. i2 .and. ib2 .eq. j2 .and. 
     &                   ic2 .eq. k2 .and. id2 .eq. l2 .and.
     &                   nn1 .eq. n1 .and. nn2 .eq. n2 ) goto 121
                     if( ia1 .eq. l1 .and. ib1 .eq. k1 .and. 
     &                   ic1 .eq. j1 .and. id1 .eq. i1 .and.
     &                   ia2 .eq. i2 .and. ib2 .eq. j2 .and. 
     &                   ic2 .eq. k2 .and. id2 .eq. l2 .and.
     &                   nn1 .eq. n1 .and. nn2 .eq. n2 ) goto 121
                     if( ia1 .eq. i1 .and. ib1 .eq. j1 .and. 
     &                   ic1 .eq. k1 .and. id1 .eq. l1 .and.
     &                   ia2 .eq. l2 .and. ib2 .eq. k2 .and. 
     &                   ic2 .eq. j2 .and. id2 .eq. i2 .and.
     &                   nn1 .eq. n1 .and. nn2 .eq. n2 ) goto 121
                     if( ia1 .eq. l1 .and. ib1 .eq. k1 .and. 
     &                   ic1 .eq. j1 .and. id1 .eq. i1 .and.
     &                   ia2 .eq. l2 .and. ib2 .eq. k2 .and. 
     &                   ic2 .eq. j2 .and. id2 .eq. i2 .and.
     &                   nn1 .eq. n1 .and. nn2 .eq. n2 ) goto 121
                  endif
               else
                 write(6,'(/,a)')'  **********************************'
                 write(6,'(a45,1x,a24)') 
     &           'ERROR: coupling should be stored for last !! ',
     &           ' Please check definitions'
                 write(6,'(/,a)')'  **********************************'
                 stop 'internal error'
               endif
            enddo
            write(6,'(/,a)')  '  **********************************'
            write(6,'(a53,1x,a28)') 
     &         'ERROR: cannot match atoms in dihedral-dihedral coupl.',
     &         ' while writing new.system'
            write(6,'(/,a)')  '  **********************************'
            write(6,*)  'readen indexes :',i1,j1,k1,l1,i2,j2,k2,l2
            write(6,*)  'readen n       :',n1,n2
            write(6,*)  'stored indexes '
            write(6,*)  'torsions'
            do ipp = 1,Npar
               if (iflag(ipp,0) .eq. 4) then
                 write(6,*) ipp,(iflag(ipp,kk),kk=1,4),'with n =',
     &                      iflag(ipp,5)
               endif
            enddo
            stop 'output error'
121         continue
            cost = par(ipc)
c           conversion: 1 H = 2625.5 kJ/mol
            cost = cost*2625.5d0
            n1       = iflag(icor1,5)
            gamma1   = parfun(icor1)*rad2deg
            n2       = iflag(icor2,5)
            gamma2   = parfun(icor2)*rad2deg
              write(25,313) 'dihed_dihed',ia1,ib1,ic1,id1,n1,gamma1,
     &                       ia2,ib2,ic2,id2,n2,gamma2,cost
              goto 14
         endif

      ELSEIF( word(1) .eq. 'begin{configuration}') THEN   

         write(25,'(a)') card(1:Lenw(card))
 16      read(iun,1000) string
         write(25,'(a)') string(1:Lenw(string))
         if(index(string,'#') .eq. 1) goto 16
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. 'end{configuration}') goto 10
         goto 16

      ENDIF

 99   write(25,'(a)') card(1:Lenw(card))
      goto 10


 100  continue



      write(6,'(/,3x,30("**"))')  
      write(6,'(/,3x,a)')  
     &'Moscito .system file correctly written in file: new.system'
      write(6,'(/,3x,30("**"))')  



300   format(a8,3x,i2,3x,i2,3x,f10.4,3x,f7.3)
301   format(a6,5x,i2,2x,i2,2x,i2,3x,f10.4)
302   format(11x,i2,2x,i2,2x,i2,3x,f10.4,3x,f6.2)
303   format(a6,3x,i2,1x,i2,1x,i2,1x,i2,3x,4(f9.3,2x,i2,2x,f6.1))
304   format(a8,3x,i2,1x,i2,1x,i2,1x,i2,3x,f9.3,2x,f6.1)
305   format(a10,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,1x,i2,3x,f6.2,2x,f10.4)
306   format(a13,1x,i2,1x,i2,2x,f6.2,
     &       3x,i2,1x,i2,1x,i2,1x,i2,3x,i2,3x,f8.2,2x,f10.4)
307   format(a11,1x,i2,1x,i2,1x,i2,2x,f6.2,
     &       3x,i2,1x,i2,1x,i2,1x,i2,3x,i2,3x,f8.2,2x,f10.4)
308   format(a6,3x,i2,1x,i2,1x,i2,1x,i2,3x,f9.3,2x,i2,2x,f7.2)
309   format(a9,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,3x,f6.2,2x,f10.4)
310   format(a11,2x,i2,1x,i2,1x,i2,2x,f7.3,
     &       3x,i2,1x,i2,1x,i2,3x,f6.2,2x,f10.4)
311   format(a13,2x,i2,1x,i2,2x,f7.3,3x,i2,1x,i2,1x,i2,1x,i2,
     &       3x,f10.4,2x,i2,3x,f8.2)
312   format(i2,3x,i2,3x,f10.4,3x,f13.7,3x,f12.4)
313   format(a11,1x,4(i2,1x),1x,i2,2x,f6.2,
     &       3x,4(i2,1x),1x,i2,2x,f6.2,2x,f10.4)
1000  format(a)

      return
      end


