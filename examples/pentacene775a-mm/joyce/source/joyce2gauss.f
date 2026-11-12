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
      Subroutine Print_Gaussian (iflag,par,parfun,npar,
     &                        nat,maxpar,maxat,aname,labpar,
     &                        nsit,name,amass,charge,geom,
     &                        sigma,epsilon,namesit)
c----------------------------------------------------------------------c
 
 
c      writes a joyce.g09.com file for Gaussian program
c      substituting Joyce computed force constants
c      and equilibrum values (R0, theta0) read from G09 .fchk files

      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      
      data two6/1.122462048309373d0/  !! 2**(1/6)  r*=two6*sigma  for LJ
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      data au2kcal /6.27509541d2/     !! 1 a.u. = 627.51 kcal/mol
 
      character name(maxat)*5,aname(*)*5,labpar(*)*(*),namesit(maxat)*5
      character string*7,sitname*2

      logical   zprtg09(maxpar),z_zero

      dimension iflag(maxpar,0:8)
      dimension par(maxpar),parfun(maxpar)
      dimension apar(maxpar),aparfun(maxpar)
      dimension amass(*),charge(*)
      dimension geom(maxat,6)
      dimension cost(maxpar,4),ngamma(maxpar,4)
      dimension id1(maxpar),id2(maxpar),id3(maxpar),id4(maxpar)
      dimension epsilon(maxat),sigma(maxat)
      dimension eps(maxat),sig(maxat)


c     converting to AMBER-like units and potential functions
c     stretch: k*(r-r0)^2     [k] = kcal/mol A^-2   ; [R0]= A
c     bending: k*(t-t0)^2     [k] = kcal/mol rad^-2 ; [t0]= deg
c     harmonic torsion: k/2(p)^2a [k] = kcla/mol rad^-2 ; note t0 forced to 0 
c     torsion: k*(1+cos(np+g))[k] = kcal/mol        ; [g] = deg
c     LJ     : [eps] = kcal  [R0] = sigma*2^-1/6 A

      write(6,'(3x,a)') 
      write(6,'(3x,a)') 
     &     '---------------------------------------------------------'
      write(6,'(3x,a)') 
     &     'Checking consistency of final parameters for G09 format'
      write(6,'(3x,a)') 
     &     '---------------------------------------------------------'
      write(6,'(3x,a)') 



      nstr = 0
      nbnd = 0
      nhtr = 0
      nftr = 0
      nljc = 0
      nqqi = 0
      do i = 1,Npar
         do k = 1,4
            ngamma(i,k) = 0
            cost(i,k)   = 0.d0
         enddo
      enddo
      do ip = 1,Npar
         zprtg09(ip)=.true.
         if (iflag(ip,0) .eq. 1) then
          nstr = nstr + 1
          apar(ip)    = 0.5d0*par(ip)*1.d4/(1.066*4.184d0)
          aparfun(ip) = parfun(ip)*0.529177249d0
         elseif (iflag(ip,0) .eq. 2) then
          nbnd = nbnd + 1
          apar(ip)    = 0.5d0*par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
         elseif (iflag(ip,0) .eq. 3) then
          nhtr = nhtr + 1
          apar(ip)    = par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
          if (aparfun(ip) .ne. 0.d0) then
            write(6,'(3x,a)') 
            write(6,'(3x,a)') 
     &       'Please check employed dihedral functions:'
            write(6,'(3x,a)') 
     &      'Harmonic dihedrals are supported in g09 only for phi0 = 0'
            stop 'Inconsistent dihedral function'
          endif  
         elseif (iflag(ip,0) .eq. 4) then
          nftr = nftr + 1
          apar(ip)    = par(ip)*au2kcal
          aparfun(ip) = parfun(ip)*rad2deg
         elseif (iflag(ip,0) .eq. 5) then
          nljc = nljc + 1
          zprtg09(ip)=.false.
         elseif (iflag(ip,0) .eq. 6) then
          nqqi = nqqi + 1
          zprtg09(ip)=.false.
         endif
      enddo
      do isit = 1,Nsit
         sig(isit) = sigma(isit)
         eps(isit) = epsilon(isit)/4.184d0
      enddo

      do ip = 1, Npar-1
       do jp = ip+1, Npar
        if(zprtg09(ip) ) then
         if(iflag(ip,0) .eq. iflag(jp,0)) then 
          if (labpar(jp) .eq. labpar(ip)) then
           if(apar(ip) .eq. apar(jp) ) then
            if(aparfun(ip) .eq. aparfun(jp)) then   
               zprtg09(jp) = .false.
            else
               p1          = aparfun(ip)
               p2          = aparfun(jp)
               pm          = 0.5d0*(p1 + p2)
               diff        = dabs(p1 - p2)
               aparfun(ip) = pm
               aparfun(jp) = pm
               if (iflag(ip,0) .le. 4) then 
                zprtg09(jp) = .false.
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

c     opens new g09 formatted input file file (joyce.g09.com)

      open(63,file='joyce.g09.com',form='formatted',status='unknown')
      rewind (63)

      write(63,'(a)') '%chk=opt.MMfromJoyce.chk'
      write(63,'(a)') '%nprocshared=4'
      write(63,'(a)') 
     &  '#p opt=(nomicro,big) amber=SoftOnly iop(4/33=3)'
      write(63,'(a)') 
      write(63,'(a)')  'MM optimization with Joyce parameterized FF' 
      write(63,'(a)') 
      write(63,'(a)')  '0   1'
      do iat = 1, Nat
         write(string,'(f7.4)') charge(iat)
         call compact(string,lw)
         sitname = name(iat)
         call compact(sitname,lw2)
         if (lw .eq. 6) then
                 if (lw2 .eq. 2) then      
                    write(63,'(a1,a1,a2,a1,a6,5x,3(f13.6,2x))') 
     &              aname(iat),'-',name(iat),'-',string(1:lw),
     &              (geom(iat,k)*0.529177249d0,k=1,3)
                 else
                    write(63,'(a1,a1,a1,a1,a6,6x,3(f13.6,2x))') 
     &              aname(iat),'-',name(iat),'-',string(1:lw),
     &              (geom(iat,k)*0.529177249d0,k=1,3)
                 endif
         else
                 if (lw2 .eq. 2) then      
                    write(63,'(a1,a1,a2,a1,a7,4x,3(f13.6,2x))') 
     &              aname(iat),'-',name(iat),'-',string(1:lw),
     &              (geom(iat,k)*0.529177249d0,k=1,3)
                 else
                    write(63,'(a1,a1,a1,a1,a7,5x,3(f13.6,2x))') 
     &              aname(iat),'-',name(iat),'-',string(1:lw),
     &              (geom(iat,k)*0.529177249d0,k=1,3)
                 endif
         endif
      enddo
      if (nstr .gt. 0) then 
        write(63,'(a)')  
        write(63,'(a)')  '!'
        write(63,'(a)')  '! Stretches '
        write(63,'(a)')  '!'
      endif

      icount = 1
      do ip = 1, Npar
       if(zprtg09(ip) .and. iflag(ip,0) .eq. 1) then
        icount = icount + 1
        ia = iflag(ip,1)
        ib = iflag(ip,2)
        write(63,120) 'HrmStr1',name(ia),name(ib),apar(ip),aparfun(ip)
       endif
      enddo
120   format(a7,1x,a3,1x,a3,2x,f8.3,2x,f6.3)

      if (nbnd .gt. 0) then 
       write(63,'(a)')  '!'
       write(63,'(a)')  '! Angles'
       write(63,'(a)')  '!'
      endif
      icount = 1
      do ip = 1, Npar
       if(zprtg09(ip) .and. iflag(ip,0) .eq. 2) then
        icount = icount + 1
        ia = iflag(ip,1)
        ib = iflag(ip,2)
        ic = iflag(ip,3)
        write(63,121) 'HrmBnd1',name(ia),name(ib),name(ic),
     &                 apar(ip),aparfun(ip)
       endif
      enddo
121   format(a7,1x,a3,1x,a3,1x,a3,2x,f8.3,2x,f7.2)

      if (nhtr .gt. 0) then 
        write(63,'(a)')  '!'
        write(63,'(a)')  '! Dihedrals harmonic'
        write(63,'(a)')  '!'
      endif
      icount = 1
      do ip = 1, Npar
       if(zprtg09(ip) .and. iflag(ip,0) .eq. 3) then
        icount = icount + 1
        ia = iflag(ip,1)
        ib = iflag(ip,2)
        ic = iflag(ip,3)
        id = iflag(ip,4)
        write(63,122) 'HrmWil1',name(ia),name(ib),name(ic),name(id),
     &                 apar(ip)
       endif
      enddo
122   format(a7,1x,a3,1x,a3,1x,a3,1x,a3,2x,f8.3)

      id     = 0   ! dihedral number
      istart = 1   ! counter for cosine terms
10    continue      
      do ip = istart, Npar
         if(zprtg09(ip) .and. iflag(ip,0) .eq. 4) then  ! Fourier-like diheds
          Ncos = iflag(ip,6)   ! number of cosine terms 
          if (ncos .gt. 5) then 
             write(6,'(a)') 'Joyce has stopped !'
             write(6,'(a)') 'Consistency error'
             write(6,'(a)') 
     &       'Gaussian FF does not support more than 4 cos terms '
             write(6,'(a)') 'Too many terms for one dihedral!' 
             write(6,'(a)') 
             write(63,'(a)') 'Joyce has stopped !'
             write(63,'(a)') 'Consistency error'
             write(63,'(a)') 
     &       'Gaussian FF does not support more than 4 cos terms '
             write(63,'(a)') 'Too many terms for one dihedral!' 
             stop 'Format consistency error with G09'
          endif
          id     = id + 1
          z_zero = .false.    
          do icos = 1, Ncos
             id1(id)     = iflag(ip,1)
             id2(id)     = iflag(ip,2)
             id3(id)     = iflag(ip,3)
             id4(id)     = iflag(ip,4)
             jp          = ip + icos -1 ! groups cos of the same dihed
             nmult       = iflag(jp,5)  ! cos(nmult*phi)
             if (nmult .gt. 4 ) then
              write(6,'(a)') 'Joyce has stopped !'
              write(6,'(a)') 'Consistency error'
              write(6,'(a)') 
     &        'Gaussian FF does not support n > 4 for cos terms '
              write(6,'(a)') 'Check you phase factors!' 
              write(6,'(a)') 
              write(63,'(a)') 'Joyce has stopped !'
              write(63,'(a)') 'Consistency error'
              write(63,'(a)') 
     &        'Gaussian FF does not support n > 4 for cos terms '
              write(63,'(a)') 'Check you phase factors!' 
              stop 'Format consistency error with G09'
             elseif (nmult .eq. 0) then ! n=0
                write(6,'(5x,a)') 
     &          '*****************************************************'
                write(6,'(27x,a)') 'WARNING'
                write(6,'(5x,a)') 
     &          '*****************************************************'
                write(6,'(a)') 
                write(6,'(5x,a)') 
     &          'Final FF contains a n = 0 constant in some dihedral'
                write(6,'(5x,a)') 
     &          'This will be ingored by G09 output !'
                write(6,'(a)') 
                write(6,'(a)') 
                z_zero         = .true.
             else
              cost(id,nmult)  = apar(jp)
              ngamma(id,nmult)= int(aparfun(jp))
              istart         = ip + ncos 
             endif        
             enddo    ! end loop on cos terms of the same dihed
             goto 10
         endif
      enddo           ! end loop on parameters


      Ndihed = id
      if (nftr .gt. 0) then 
        write(63,'(a)')  '!'
        write(63,'(a)')  '! Torsions '
        write(63,'(a)')  '!'
      endif
      do kd = 1,Ndihed
         ia = id1(kd)
         ib = id2(kd)
         ic = id3(kd)
         id = id4(kd)
         write(63,123) 'AmbTrs',name(ia),name(ib),name(ic),name(id),
     &                 (ngamma(kd,icos),icos=1,4),
     &                 (cost(kd,icos),icos=1,4),'1.0'  
      enddo
123   format(a6,1x,a2,1x,a2,1x,a2,1x,a2,2x,4(i3,1x),1x,4(f8.3,1x),1x,a3)

      if (nljc .gt. 0) then
         icount = 1
         do isit = 1, Nsit
               if (icount .eq. 1)  then 
                 write(63,'(a)')  '!'
                 write(63,'(a)')  '! Vanderwaals parameters '
                 write(63,'(a)')  '!'
               endif 
               icount = icount + 1
               write(63,124) 'VDW',namesit(isit),sig(isit),eps(isit)
         enddo
124      format(a3,1x,a3,1x,f7.3,2x,f10.4)
         write(63,'(a)')  '!'
         write(63,'(a)')  '! Nonbonded rules  '
         write(63,'(a)')  '!'
         write(63,'(a)')  ' NonBon  6 1  0 0  0.0 0.0 0.5 0.0 0.0 0.833'
      endif   
      write(63,'(a)') 
      write(63,'(a)') 
      write(63,'(a)')  

      return
      end
