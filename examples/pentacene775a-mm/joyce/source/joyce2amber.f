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
      Subroutine Print_Amber (iflag,par,parfun,npar,
     &                        nat,maxpar,maxat,aname,labpar,
     &                        nsit,namesit,amass,charge)
c----------------------------------------------------------------------c
 
c     writes a joyce.new.par file for AMBER MD program
c     substituting Joyce computed force constants
c     and equilibrum values (R0, theta0) read from G03 .fchk files

      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      
      data two6/1.122462048309373d0/  !! 2**(1/6)  r*=two6*sigma  for LJ
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      data au2kcal /6.27509541d2/     !! 1 a.u. = 627.51 kcal/mol
 
c     global variables

      character        namesit(*)*5,aname(*)*5,labpar(*)*(*)

      logical   zprtamber(maxpar)

      dimension iflag(maxpar,0:8)
      dimension par(maxpar),parfun(maxpar)
      dimension apar(maxpar),aparfun(maxpar)
      dimension amass(*),charge(*)

c     converting to AMBER units and potential functions
c     stretch: k*(r-r0)^2     [k] = kcal/mol A^-2   ; [R0]= A
c     bending: k*(t-t0)^2     [k] = kcal/mol rad^-2 ; [t0]= deg
c     Note: no harmonic torsion supported
c     torsion: k*(1+cos(np+g))[k] = kcal/mol        ; [g] = deg
c     LJ     : [eps] = kcal  [R0] = 0.5*sigma/2^-1/6 A

      write(6,'(3x,a)') 
      write(6,'(3x,a)') 
     &     '---------------------------------------------------------'
      write(6,'(3x,a)') 
     &     'Checking consistency of final parameters for AMBER format'
      write(6,'(3x,a)') 
     &     '---------------------------------------------------------'
      write(6,'(3x,a)') 


      do ip = 1,Npar
         zprtamber(ip)=.true.
         if (iflag(ip,0) .eq. 1) then
          apar(ip)    = 0.5d0*par(ip)*1.d4/(1.066*4.184d0)
          aparfun(ip) = parfun(ip)*0.529177249d0
         elseif (iflag(ip,0) .eq. 2) then
          apar(ip)    = 0.5d0*par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
         elseif (iflag(ip,0) .eq. 3) then
          apar(ip)    = 0.5d0*par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
          write(6,'(3x,a)') 
          write(6,'(3x,a)') 
     &     'Please check employed dihedral funtions:'
          write(6,'(3x,a)') 
     &     'Harmonic dihedrals are not supported in AMBER'
             stop 'Inconsistent dihedral function'
         elseif (iflag(ip,0) .eq. 4) then
          apar(ip)    = par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
         elseif (iflag(ip,0) .eq. 5) then
          apar(ip)    = par(ip)*au2kcal
          aparfun(ip) = 0.5d0*parfun(ip)*0.529177249d0*two6
         elseif (iflag(ip,0) .eq. 6) then
          zprtamber(ip)=.false.
         endif
      enddo

      do ip = 1, Npar-1
       do jp = ip+1, Npar
        if(zprtamber(ip) ) then
         if(iflag(ip,0) .eq. iflag(jp,0)) then 
          if (labpar(jp) .eq. labpar(ip)) then
           if(apar(ip) .eq. apar(jp) ) then
            if(aparfun(ip) .eq. aparfun(jp)) then   
               zprtamber(jp) = .false.
            else
               p1 = aparfun(ip)
               p2 = aparfun(jp)
               pm = 0.5d0*(p1 + p2)
               diff        = dabs(p1 - p2)
               aparfun(ip) = pm
               aparfun(jp) = pm
               if (iflag(ip,0) .ge. 4) then 
                zprtamber(jp) = .false.
               endif
               if (iflag(ip,0) .eq. 1 .and. diff .ge. 0.1d0 ) then
                write(6,'(12x,a)') 'Warning'
                write(6,'(12x,a20,1x,a9,1x,a3,1x,a9)')
     &          'Equil value of terms',labpar(ip),'and',
     &          labpar(jp)
                write(6,'(12x,a15,1x,f12.5,1x,a4,1x,f12.5,1x,a3,f12.5)') 
     &          'was averaged to',pm,'from',p1,'and',p2
               elseif (iflag(ip,0) .eq. 2 .and. diff .ge. 2.0d0 ) then
                write(6,'(12x,a)') 'Warning'
                write(6,'(12x,a20,1x,a9,1x,a3,1x,a9)')
     &          'Equil value of terms',labpar(ip),'and',
     &          labpar(jp)
                write(6,'(12x,a15,1x,f12.5,1x,a4,1x,f12.5,1x,a3,f12.5)') 
     &          'was averaged to',pm,'from',p1,'and',p2
               elseif (iflag(ip,0) .ge. 3 .and. diff .ge. 5.0d0 ) then
                write(6,'(12x,a)') 'Warning'
                write(6,'(12x,a20,1x,a9,1x,a3,1x,a9)')
     &          'Equil value of terms',labpar(ip),'and',
     &          labpar(jp)
                write(6,'(12x,a15,1x,f12.5,1x,a4,1x,f12.5,1x,a3,f12.5)') 
     &          'was averaged to',pm,'from',p1,'and',p2
               endif
            endif
           else
            if(iflag(ip,0) .le. 4) then
             write(6,'(3x,a)') 'Warning'
             write(6,'(3x,a23,1x,a15,1x,a3,1x,a15)') 
     &       'Force constant of terms',labpar(ip),'and',
     &       labpar(jp)
             write(6,'(3x,a28,1x,i3,1x,a3,1x,i3)') 
     &        'corresponding to parameters',ip, 'and',jp        
             write(6,'(3x,a)') 'have different values:'
             write(6,'(3x,f10.5,1x,a3,1x,f10.5)') 
     &       apar(ip),'and',apar(jp)
             write(6,'(3x,a)') 
     &       'Please check consistency of attypes definition'
             stop 'Inconsistent attypes'
            endif
           endif
          endif
         endif
        endif
       enddo
      enddo

C     Open new amber.par formatted file (joyce.new.frcmod)

      open(25,file='joyce.new.frcmod',form='formatted',status='unknown')
      rewind (25)



      write(25,'(a)') '# AMBER PARAMETERS FROM JOYCE FITTING '
      write(25,'(a)') 

      write(25,'(a)') 'MASS' 
      do isit = 1, Nsit
         write(25,'(a5,3x,f6.2)') namesit(isit),amass(isit)
      enddo
      write(25,'(a)') 

      icount = 1
      do ip = 1, Npar
         if(zprtamber(ip) .and. iflag(ip,0) .eq. 1) then
            if (icount .eq. 1)  write(25,'(a)') 'BOND'
            icount = icount + 1
            write(25,120) labpar(ip),apar(ip),aparfun(ip)
         endif
      enddo
      write(25,'(a)') 
120   format(a5,4x,f9.2,4x,f6.3)

      icount = 1
      do ip = 1, Npar
         if(zprtamber(ip) .and. iflag(ip,0) .eq. 2) then
            if (icount .eq. 1)  write(25,'(a)') 'ANGLE'
            icount = icount + 1
            write(25,121) labpar(ip),apar(ip),aparfun(ip)
         endif
      enddo
      write(25,'(a)') 
121   format(a8,4x,f9.2,4x,f8.2)

      write(25,'(a)') 
      write(25,'(a)') 
      write(25,'(a)') '# WARNING WARNING WARNING WARNING WARNING  '
      write(25,'(a)') '# Please carefully check consistency'
      write(25,'(a)') '# in multiple cosines defined on the same'
      write(25,'(a)') '# torsion by adjusting labels and multiplicity'
      write(25,'(a)') '# according to AMBER formalism (negative n)'
      write(25,'(a)') '# WARNING WARNING WARNING WARNING WARNING  '
      write(25,'(a)') 
      write(25,'(a)') 

      icount = 1
      npath  = 1
      do ip = 1, Npar
         if(zprtamber(ip) .and. iflag(ip,0) .eq. 4) then
            if (icount .eq. 1)  write(25,'(a)') 'DIHEDRAL'
            icount = icount + 1
            dn = dble(iflag(ip,5))
            write(25,122) labpar(ip),npath,apar(ip),aparfun(ip),dn
         endif
      enddo
      write(25,'(a)') 
122   format(a11,3x,i1,2x,f9.2,4x,f8.2,3x,f5.1)

      write(25,'(a)') 
      write(25,'(a)') 
      write(25,'(a)') '# WARNING WARNING WARNING WARNING WARNING  '
      write(25,'(a)') '# Please carefully check consistency'
      write(25,'(a)') '# between intramolecular LJ and Charges'
      write(25,'(a)') '# parameters employed in the JOYCE fittig'
      write(25,'(a)') '# and the ones to be employed with AMBER'
      write(25,'(a)') '# WARNING WARNING WARNING WARNING WARNING  '
      write(25,'(a)') 
      write(25,'(a)') 

      icount = 1
      do ip = 1, Npar
         if(zprtamber(ip) .and. iflag(ip,0) .eq. 5) then
            if (icount .eq. 1)  write(25,'(a)') 'NONBOND'
            icount = icount + 1
            write(25,123) labpar(ip),aparfun(ip),apar(ip)
         endif
      enddo
123   format(a6,3x,f10.4,4x,f10.4)
      write(25,'(a)') 
      write(25,'(a)') 
      write(25,'(a)') 

      return
      end
