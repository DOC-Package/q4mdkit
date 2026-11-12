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
      Subroutine Read_Gromacs (iun,iflag,par,parfun,parfun2 
     &                         ,npar,geom,nat
     &                         ,maxpar,maxat,aname,labpar,name,ichrgprm
     &                         ,nsit,namesit,cr,eps_sit,sig_sit
     &                         ,amass,charge,sigma,epsilon,z_lj)
c----------------------------------------------------------------------c

c               reads a Gromacs .top file as input and
c               extract an 'iflag' matrix with the requested potential
c               function according to Joyce conventions

c         *****        JOYCE conventions on iflag  matrix  ******

c         iflag(j,k)   k-th label of the j-th parameter

c         in particular, for k = 0,        
 
c         iflag(ipar,0) identifies parameter type i.e.:
c         iflag(ipar,0)  = 1     stretching
c         iflag(ipar,0)  = 2     bending
c         iflag(ipar,0)  = 3     torsion (harmonic)
c         iflag(ipar,0)  = 4     torsion proper  (with cosines)
c         iflag(ipar,0)  = 5     LJ intramolecular
c         iflag(ipar,0)  = 6     charge-charge intramolecular
c         iflag(ipar,0)  = 8     coupling               ! changed by GP in Joyce3

c         iflag(ipar,1)  = i     atom  i which is needed for the definition
c                                of the internal coordinate 
c                                contained in iflag(ipar,0)
c         iflag(ipar,2)  = j     atom  j which is needed for the definition
c                                of the internal coordinate 
c                                contained in iflag(ipar,0)
c         iflag(ipar,3)  = k     atom  k (if present) which is needed 
c                                for the definition of the i.c.
c                                contained in iflag(ipar,0)
c         iflag(ipar,4)  = l     atom  l (if present) which is needed for 
c                                the definition of the internal coordinate 
c                                contained in iflag(ipar,0)

c         Note : i < j < k < l

c         for iflag(ipar,0) = 4  (torsions)

c         iflag(ipar,5)  = n     only if iflag(ipar,0) = 4 
c                                is the n value for the cosine term
c                                Vn = 1 + cos(n*phi - gamma)
c         iflag(ipar,6)  = N     only if iflag(ipar,0) = 4 
c                                is the number of cosine terms
c         iflag(ipar,8)  = nat   numer of atoms involved in IC definition

c         Couplings: if iflag(ipar,0)  = 8

c         iflag(ipar,1)  = i     i-th coupled function
c         iflag(ipar,2)  = j     j-th coupled function
c         iflag(ipar,5)  = n     only if a proper dihedral is coupled
c                                is the n value for the cosine term
c         iflag(ipar,6)  = n2    only if two proper dihedral are coupled
c                                is the n value for the cosine term of the second f.
c         iflag(ipar,8)  = nat   numer of atoms involved in the coupling
                                
c               *****        Potential model functions ******

c         Stretching : Vs  = Ks  * (r     - r0)^2
c         Bending    : Vb  = Kb  * (theta - theta0)^2
c         Torsion  i.: Vdi = Kdi * (phi   - phi0  )^2
c         Torsion    : V0 + V1 + V2 + .....
c                      with 
c                      Vn  = Kn * [1 + cos(n*phi - gamma)]
c
c        Str-Str       Vss = kss * (r -rc1) * (r-rc2)
c        Str-Bnd       Vsb = kbs * (r -rc) * (theta-tc)
c        Str-Tors      Vst = kst * (r -rc) * (p-pc)    <== cfype = 2
c                      or 
c                      Vst = kst * (r -rc) * [V0 + V1 + ... VNc]  <== cftype = 1
c                      with 
c                      Vn = [sin(nc*phi - pc)]
c        Bnd-Bnd       Vbb = kbb * (theta -tc1) * (theta-tc2)
c        Bnd-Tors      Vbt = kbt * (theta -tc) * (p-pc)  <== cftype = 2
c                      or 
c                      Vst = kvt * (theta -tc) * [V0 + V1 + ... VNc]  <== cftype = 1
c                      with 
c                      Vn = [sin(nc*phi - pc)]
c        Tors-Tors     Vtt = ktt * (phi -pc1) * (phi-pc2)  <== cftype = 3
c                      or 
c                      Vtt = ktt * (phi -pc1) * [V0 + V1 + ... VNc]  <== cftype = 2
c                      with 
c                      Vn = [sin(nc2*phi - pc2)]
c                      or
c                      Vtt = ktt * [V0 + V1 + ... VNc] * [V'0 + V'1 + ... V'Nc]  <== cftype = 1
c                      with 
c                      Vn = [sin(nc1*phi - pc1)] and V'n = [sin(nc2*phi - pc2)]

 
      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      parameter       (boh2nm  = 0.0529177249d0)!! 1 bohr = 0.0529177 nm
      parameter       (H2kJ    = 2625.5d0     )!! 1 H = 2625.5d0 kJ/mol
      
c     characters
      character        card*90,string*90,card0*90,cardm*90
      character        name(maxat)*5,namesit(*)*5,col*1
      character        labchk1*5,labchk2*5
      character        aname(*)*5,labpar(*)*(*)
      character        labchk*5
      character*46     word(maxspl),subword(maxspl)
c     logicals
      logical          z_harmonic,z_linear,z_improper,z_lj,z_coupl
      logical          z_intcharge
      logical          zchk(2)
    
      integer   cftype,cr
c     diagonal force constants
      real*8    k_s,k_b,k_d,k_di
c     coupled force constants
      real*8    k_ss, k_sb, k_st
      real*8    k_bb,k_bt
      real*8    k_tt
c     conversion  factor
      real*8    j2kjmol

c     dimensioning
      dimension x(maxint),y(maxint),z(maxint)
      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension par(*),parfun(*),parfun2(*)
      dimension iref(maxpar)
c     atoms defining simple ICs
      dimension iat_s(maxint),jat_s(maxint),k_s(maxint),r0_s(maxint)
      dimension iat_b(maxint),jat_b(maxint),kat_b(maxint),k_b(maxint),
     &          t0_b(maxint)
      dimension iat_d(maxint),jat_d(maxint),kat_d(maxint),lat_d(maxint),
     &          k_d(maxint,maxcos),n0_d(maxint,maxcos),
     &          p0_d(maxint,maxcos),ic(maxint)
      dimension iat_di(maxint),jat_di(maxint),kat_di(maxint),
     &          lat_di(maxint),k_di(maxint),p0_di(maxint)
c     couplings
      dimension z_coupl(6)
      dimension cftype(6,maxint)
c     atoms defining coupled ICs
      dimension iat_ss(2),jat_ss(2)
      dimension iat_sb(2),jat_sb(2),kat_sb(2)
      dimension iat_st(2),jat_st(2),kat_st(2),lat_st(2)
      dimension iat_bb(2),jat_bb(2),kat_bb(2)
      dimension iat_bt(2),jat_bt(2),kat_bt(2),lat_bt(2)
      dimension iat_tt(2),jat_tt(2),kat_tt(2),lat_tt(2)
      dimension iic_ss(maxint),jic_ss(maxint)
      dimension iic_sb(maxint),jic_sb(maxint)
      dimension iic_bb(maxint),jic_bb(maxint)
      dimension iic_st(maxint),jic_st(maxint)
      dimension iic_bt(maxint),jic_bt(maxint)
      dimension iic_tt(maxint),jic_tt(maxint)
c     constants defining coupled ICs (with keepFF)
      dimension k_ss(maxint),k_sb(maxint),k_st(maxint,maxcos)
      dimension k_bb(maxint),k_bt(maxint,maxcos)
      dimension k_tt(maxint,maxcos)
      dimension rc_ss(maxint,2),tc_bb(maxint,2)
      dimension rc_sb(maxint),tc_sb(maxint)
      dimension rc_st(maxint),tc_bt(maxint)
      dimension pc_st(maxint,maxcos),nc_st(maxint,maxcos)
      dimension pc_bt(maxint,maxcos),nc_bt(maxint,maxcos)
      dimension pc_tt(maxint,maxcos,2),nc_tt(maxint,maxcos,2)
      dimension Ncos_st(maxint)
      dimension Ncos_bt(maxint)
      dimension Ncos_tt(maxint)
c     LJ and charge interactions
      dimension iat_lj(maxint),jat_lj(maxint),qi(maxint),qj(maxint),
     &          eps_lj(maxint),sig_lj(maxint),fqq(maxint),qq(maxint)
      dimension ichrgprm(2)
      dimension amass(*),charge(*),sigma(*),epsilon(*)
      dimension sig_sit(*),eps_sit(*)
      dimension stemp(maxint),etemp(maxint)


      do i=1,maxpar
        par(i)    = 0.d0
        parfun(i) = 0.d0
        parfun2(i)= 0.d0
        do j=0,8   
          iflag(i,j)=0
        enddo
      enddo
      do i = 1,maxint
         iat_s(i)   = 0
         jat_s(i)   = 0
         iat_b(i)   = 0
         jat_b(i)   = 0
         kat_b(i)   = 0
         iat_d(i)   = 0
         jat_d(i)   = 0
         kat_d(i)   = 0
         lat_d(i)   = 0
         iat_di(i)  = 0
         jat_di(i)  = 0
         kat_di(i)  = 0
         lat_di(i)  = 0
         Ncos_st(i) = 1
         Ncos_bt(i) = 1
         Ncos_tt(i) = 1
         do j = 1,6
            cftype(j,i) = 0
         enddo
         do j = 1,2
            nc_st(i,j) = 0
            nc_bt(i,j) = 0
            do icos=1,maxcos
             nc_tt(i,icos,j) = 0
            enddo
         enddo
      enddo
      do i = 1,2
         ichrgprm(i) = 0
         iat_ss(i)   = 0
         jat_ss(i)   = 0
         iat_sb(i)   = 0
         jat_sb(i)   = 0
         kat_sb(i)   = 0
         iat_st(i)   = 0
         jat_st(i)   = 0
         kat_st(i)   = 0
         lat_st(i)   = 0
         iat_bb(i)   = 0
         jat_bb(i)   = 0
         kat_bb(i)   = 0
         iat_bt(i)   = 0
         jat_bt(i)   = 0
         kat_bt(i)   = 0
         lat_bt(i)   = 0
         iat_tt(i)   = 0
         jat_tt(i)   = 0
         kat_tt(i)   = 0
         lat_tt(i)   = 0
      enddo

      z_harmonic = .false.
      z_linear   = .false.
      z_improper = .false.
      z_lj       = .false.
      z_intcharge= .false.
      do k = 1,Nat
         charge(k) = 0.d0
      enddo
c     Couplings 
      do k = 1,6
         z_coupl(k) = .false.
      enddo

      Nstretch     = 0
      Nbend        = 0
      Nlj          = 0
      Nbl          = 0
      Ndihedral    = 0
      Ndihedral_i  = 0
      Ncoupl       = 0
      Nstrstr      = 0
      Nstrbend     = 0
      Nstrtors     = 0
      Nbendbend    = 0
      Nbendtors    = 0
      Ntorstors    = 0

      string = 'blank'

c     Reads topol.top files from standard input 
 
      call flush (6)
      iline = 0
 10   READ(iun,'(a)',err=2232,end=100) card
      lw11 = lenw(card) 
      iline = iline + 1
      card0 = card
      if( index(card,'#') .eq. 1) goto 10
      if( index(card,';') .eq. 1) goto 10
      if( Lenw(card).lt.2)        goto 10
      call setcomm(card)          ! eliminates comments after ";"
      CALL SPLIT (card,word,nsplit,maxspl)

      if (word(1) .eq. '#') goto 10
      if (word(1) .eq. ';') goto 10
      if (nsplit .eq. 0) goto 10

      IF( word(2) .eq. 'defaults') THEN

 18      read(iun,1000) string
         string = adjustl(string)   
         if( index(string,'#')  .eq. 1) goto 18
         if( index(string,';')  .eq. 1) goto 18
         call setcomm(string)     ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. ';') goto 18
         if(msplit .eq. 0)   then  !  blank line
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            backspace(iun)
            goto 10
         else
            read(subword(1),*,err=2233) ft  ! function type
            if (ft .ne. 1) then
              write(6,'(a)') 
     &        'A function different form LJ was requested'
              write(6,'(a)') 
     &        'in the Gromacs topology in the [defaults] section'
              write(6,'(a)') 
     &        'Unfortunately Joyce does not support it, yet'
              write(6,'(a)') 
     &        'Please check your .top file in string:'
              write(6,'(a)') string
              stop
            endif
            read(subword(2),*,err=2233) cr  ! LJ combination rule
            if (cr .ne. 2 .and. cr .ne. 3) then 
              write(6,'(a)') 
     &        'An unknown combination rule was requested'
              write(6,'(a)') 
     &        'in the Gromacs topology in the [defaults] section'
              write(6,'(a)') 
     &        'Joyce only supports Lorentz-Berthelot (2) or '
              write(6,'(a)') 
     &        'OPLS (3) combination rules.'
              write(6,'(a)') 
     &        'Please check your .top file in string:'
              write(6,'(a)') string
              stop
            endif
            continue
         endif
         goto 18
      ENDIF

      IF( word(2) .eq. 'bonds') THEN

         ns = 0
 11      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 11
         if( index(string,';') .eq. 1) goto 11
         call setcomm(string)     ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. ';') goto 11
         if(msplit .eq. 0)      then  ! blank line
            nstretch = ns
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nstretch = ns
            backspace(iun)
            goto 10
         else
            read(subword(3),*,err=2233) itype
            if (itype .eq. 1) then 
              z_harmonic = .true.
              ns = ns + 1
              read(subword(1),*,err=2233) iat_s(ns)
              read(subword(2),*,err=2233) jat_s(ns)
              if (iat_s(ns) .gt. jat_s(ns)) then  !! ctrl i < j
                 jnew      = iat_s(ns)
                 iat_s(ns) = jat_s(ns)
                 jat_s(ns) = jnew
              endif
              read(subword(4),*,err=2233) r0_s(ns)
              read(subword(5),*,err=2233) k_s(ns)
              goto 11
            else
              write(6,'(a)') 
     &        'Not harmonic stretching found in gromacs topology'
              write(6,'(a)') 
     &        'Please check your .top file in string:'
              write(6,'(a)') string
              stop
            endif
            goto 11
         endif


      ELSEIF( word(2) .eq. 'angles') THEN

         nb  = 0
c        nbl = 0            ! obsolete
c        nss  = 0           ! obsolete
 12      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 12
         if( index(string,';') .eq. 1) goto 12
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 12
         if(subword(1)(1:1) .eq. ';') goto 12
         if(msplit .eq. 0)    then  ! blank line
            nbend    = nb
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
c           card     = string
c           obsolete ! GP 9/2017
c           nbend_l  = nbl
c           nstrstr  = nss
c           ncoupl   = nstrstr !! the only allowed in std. gromacs
c           obsolete ! GP 9/2017
            nbend    = nb
            backspace(iun)
            goto 10
         else
            read(subword(4),*,err=2233) itype
            if (itype .eq. 1) then 
               nb = nb + 1
               Call Read_1atom (iat_b(nb),subword(1),nat,string)
               Call Read_1atom (jat_b(nb),subword(2),nat,string)
               Call Read_1atom (kat_b(nb),subword(3),nat,string)
               if (iat_b(nb) .gt. kat_b(nb)) then  !! ctrl i < k
                  knew      = iat_b(nb)
                  iat_b(nb) = kat_b(nb)
                  kat_b(nb) = knew
               endif
               read(subword(5),*,err=2233) t0_b(nb)
               read(subword(6),*,err=2233) k_b(nb)
               goto 12
c           obsolete ! GP 9/2017
c           elseif(itype .eq. 3) then
c           standard gromacs 3-sites str-str cross term
c              z_coupl(7)=.true.
c              nss = nss + 1
c              Call Read_1atom (iat_ss(nss),subword(1),nat,string)
c              Call Read_1atom (jat_ss(nss),subword(2),nat,string)
c              Call Read_1atom (kat_ss(nss),subword(3),nat,string)
c              read(subword(5),*,err=2233) r1e_ss(nss)
c              read(subword(6),*,err=2233) r2e_ss(nss)
c              read(subword(7),*,err=2233) krr_ss(nss)
c              goto 12
c           obsolete ! GP 9/2017
            else
              write(6,'(a)') 
     &        'Not harmonic bending found in gromacs topology'
              write(6,'(a)') 
     &        'Please check your .top file in string:'
              write(6,'(a)') string
              stop
            endif

         endif


      ELSEIF( word(2) .eq. 'dihedrals') THEN
         nd  = 0
         ndi = 0
 13      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 13
         if( index(string,';') .eq. 1) goto 13
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') goto 13
         if(subword(1)(1:1) .eq. ';') goto 13
         if(msplit .eq. 0)  then   ! blank line
            ndihedral   = nd
            ndihedral_i = ndi
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            ndihedral   = nd
            ndihedral_i = ndi
            backspace(iun)
            goto 10
         else
            read(subword(5),*,err=2233) itype
            if (itype .eq. 1) then   ! cosines
             Ndihedral   = 0
             nd = nd + 1
             Call Read_1atom (iat_d(nd),subword(1),nat,string)
             Call Read_1atom (jat_d(nd),subword(2),nat,string)
             Call Read_1atom (kat_d(nd),subword(3),nat,string)
             Call Read_1atom (lat_d(nd),subword(4),nat,string)
             if (jat_d(nd) .gt. kat_d(nd)) then  !! ctrl j < k
                lnew      = iat_d(nd)
                knew      = jat_d(nd)
                iat_d(nd) = lat_d(nd)
                jat_d(nd) = kat_d(nd)
                kat_d(nd) = knew
                lat_d(nd) = lnew
             endif
             if (index(subword(msplit),'#') .eq. 1) 
     &           msplit = msplit -1  
             if (subword(msplit-1) .eq. '#') 
     &           msplit = msplit - 2
             if (subword(msplit-2) .eq. '#') 
     &           msplit = msplit - 3
             if (subword(9) .eq. '#') msplit = 8
             ic(nd) = 0
             do iterm=6,msplit,3
               ic(nd) = ic(nd) + 1
               read(subword(iterm),*,err=2233)   p0_d(nd,ic(nd)) 
               read(subword(iterm+1),*,err=2233) k_d(nd,ic(nd))
               read(subword(iterm+2),*,err=2233) n0_d(nd,ic(nd)) 
             enddo
             goto 13
            elseif (itype .eq. 2) then  ! harmonic
              z_improper = .true.
              ndi = ndi + 1
              Call Read_1atom (iat_di(ndi),subword(1),nat,string)
              Call Read_1atom (jat_di(ndi),subword(2),nat,string)
              Call Read_1atom (kat_di(ndi),subword(3),nat,string)
              Call Read_1atom (lat_di(ndi),subword(4),nat,string)
              if (jat_di(ndi) .gt. kat_di(ndi)) then  !! ctrl j < k
                 lnew        = iat_di(ndi)
                 knew        = jat_di(ndi)
                 iat_di(ndi) = lat_di(ndi)
                 jat_di(ndi) = kat_di(ndi)
                 kat_di(ndi) = knew
                 lat_di(ndi) = lnew
              endif
              read(subword(6),*,err=2233) p0_di(ndi)
              read(subword(7),*,err=2233) k_di(ndi)
              goto 13
            endif
         endif

      ELSEIF( word(2) .eq. 'pairs') THEN

         z_lj = .true.
         lj  = 0
 15      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 15
         if( index(string,';') .eq. 1) goto 15
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 15
         if(subword(1)(1:1) .eq. ';') goto 15
         if (msplit .eq. 0 ) then    ! blank line
            nlj  = lj
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nlj  = lj
            backspace(iun)
            goto 10
         else
            read(subword(3),*,err=2233) itype
            if (itype .eq. 1) then 
               z_intcharge = .false.
               lj = lj + 1
               read(subword(1),*,err=2233) iat_lj(lj)
               read(subword(2),*,err=2233) jat_lj(lj)
               read(subword(4),*,err=2233) sig_lj(lj)
               read(subword(5),*,err=2233) eps_lj(lj)
               if (msplit .gt. 5) then
                write(6,'(a)')'Please check pairs section in .top file'
                write(6,'(a)')'Only 5 columns are expected if Ftype=1'
                stop 
               endif
               goto 15
            elseif (itype .eq. 2) then 
               z_intcharge = .true.
               lj = lj + 1
               read(subword(1),*,err=2233) iat_lj(lj)
               read(subword(2),*,err=2233) jat_lj(lj)
               read(subword(4),*,err=2233) fqq(lj)
               read(subword(5),*,err=2233) qi(lj)
               read(subword(6),*,err=2233) qj(lj)
               read(subword(7),*,err=2233) sig_lj(lj)
               read(subword(8),*,err=2233) eps_lj(lj)
               qq(lj) = fqq(lj)*qi(lj)*qj(lj)
               goto 15
            endif
         endif


c     Couplings read in the Joyce non-standard .top format
c     New in Joyce3.0    GP&JC 9/2017

      ELSEIF( word(2) .eq. 'bo_bo') THEN       ! stretching-stretching

         nss  = 0
 190     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 190
         if( index(string,';') .eq. 1) goto 190
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 190
         if(subword(1)(1:1) .eq. ';') goto 190
         if(msplit .eq. 0) then   ! blank line
            nstrstr     = nss
            ncoupl      = nss 
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nstrstr     = nss
            ncoupl      = nss 
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 8) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'8 columns expected in str-str section'
             write(6,'(a)')'Last read line:'
             write(6,*) string
             stop 
            endif
            z_coupl(1) = .true.
            nss        = nss + 1    ! str-str IC counter
            read(subword(1),*,err=2233) iat_ss(1)
            read(subword(2),*,err=2233) jat_ss(1)
            read(subword(3),*,err=2233) iat_ss(2)
            read(subword(4),*,err=2233) jat_ss(2)
            do kk = 1,2 
              if (iat_ss(kk) .gt. jat_ss(kk)) then  !! ctrl i < j
                jnew       = iat_ss(kk)
                iat_ss(kk) = jat_ss(kk)
                jat_ss(kk) = jnew
              endif
            enddo
c           search for the coupled ICs among stretchings and stores
            do is = 1, nstretch
               ichk = abs(iat_ss(1) - iat_s(is))
               jchk = abs(jat_ss(1) - jat_s(is))
               nchk = ichk + jchk 
               if (nchk .eq. 0) then 
                iic_ss(nss) = is
                zchk(1)     = .true.
               endif
               ichk = abs(iat_ss(2) - iat_s(is))
               jchk = abs(jat_ss(2) - jat_s(is))
               nchk = ichk + jchk 
               if (nchk .eq. 0) then 
                jic_ss(nss) = is
                zchk(2)     = .true.
               endif
            enddo
c           reads other parameters
            read(subword(5),*,err=2233) cftype(1,nss)
            if (cftype(1,nss) .eq. 1) then 
               read(subword(6),*,err=2233) rc_ss(nss,1)
               read(subword(7),*,err=2233) rc_ss(nss,2)
               read(subword(8),*,err=2233) k_ss(nss)
               do kk = 1,2
                  if(.not. zchk(kk)) then 
                   write(6,'(a)') 'Topology error: in Str-str couplings'
                   write(6,'(a)') 'Coupled coordinate defined by atoms'
                   write(6,'(2(i4,2x))') iat_ss(k), jat_ss(k) 
                   write(6,'(a)') 'Not found in uncoupled list'
                   stop
                  endif
               enddo
               goto 190
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'function type undefined in str-str section'
             stop 
            endif
         endif

      ELSEIF( word(2) .eq. 'bo_an') THEN       ! stretching-bending

         nsb  = 0
 191     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 191
         if( index(string,';') .eq. 1) goto 191
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 191
         if(subword(1)(1:1) .eq. ';') goto 191
         if(msplit .eq. 0 )  then   ! blank line
            nstrbend    = nsb
            ncoupl      = nstrstr + nsb
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nstrbend    = nsb
            ncoupl      = nstrstr + nsb
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 9) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'9 columns expected in str-bnd section'
             stop 
            endif
            z_coupl(2) = .true.
            nsb        = nsb + 1  ! stretching-bending IC counter
            read(subword(1),*,err=2233) iat_sb(1)
            read(subword(2),*,err=2233) jat_sb(1)
            read(subword(3),*,err=2233) iat_sb(2)
            read(subword(4),*,err=2233) jat_sb(2)
            read(subword(5),*,err=2233) kat_sb(2)
            if (iat_sb(1) .gt. jat_sb(1)) then !! ctrl i < j str
                jnew          = iat_sb(1)
                iat_sb(1) = jat_sb(1)
                jat_sb(1) = jnew
            endif
            if (iat_sb(2) .gt. kat_sb(2)) then  !! ctrl i < k bnd
               knew      = iat_sb(2)
               iat_sb(2) = kat_sb(2)
               kat_sb(2) = knew
            endif
c           search for the coupled ICs among stretch. and bend.  and stores
            do is = 1, nstretch
               ichk = abs(iat_sb(1) - iat_s(is))
               jchk = abs(jat_sb(1) - jat_s(is))
               nchk = ichk + jchk 
               if (nchk .eq. 0) then 
                iic_sb(nsb) = is
                zchk(1)     = .true.
               endif
            enddo
            do ib = 1, nbend
               ichk = abs(iat_sb(2) - iat_b(ib))
               jchk = abs(jat_sb(2) - jat_b(ib))
               kchk = abs(kat_sb(2) - kat_b(ib))
               nchk = ichk + jchk  + kchk
               if (nchk .eq. 0) then 
                jic_sb(nsb) = ib
                zchk(2)     = .true.
               endif
            enddo
c           reads other parameters
            read(subword(6),*,err=2233) cftype(2,nsb)
            if (cftype(2,nsb) .eq. 1) then 
               read(subword(7),*,err=2233) rc_sb(nsb)
               read(subword(8),*,err=2233) tc_sb(nsb)
               read(subword(9),*,err=2233) k_sb(nsb)
               do kk = 1,2
                 if(.not. zchk(kk)) then 
                  write(6,'(a)') 'Topology error: in str-bend couplings'
                  write(6,'(a)') 'Coupled coordinate defined by atoms'
                  write(6,'(3(i4,2x))') 
     &            iat_sb(kk), jat_sb(kk), kat_sb(kk) 
                  write(6,'(a)') 'Not found in uncoupled list'
                  stop
                 endif
               enddo
               goto 191
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'function type undefined in str-bnd section'
             stop 
            endif
         endif

      ELSEIF( word(2) .eq. 'bo_di') THEN       ! stretching-torsion

         nst  = 0
 192     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 192
         if( index(string,';') .eq. 1) goto 192
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 192
         if(subword(1)(1:1) .eq. ';') goto 192
         if(msplit .eq. 0) then   ! blank line
            nstrtors    = nst
            ncoupl      = nstrstr + nstrbend + nst
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nstrtors    = nst
            ncoupl      = nstrstr + nstrbend + nst
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 10) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')
     &         'at least 10 columns expected in str-tors section'
             stop 
            endif
            z_coupl(3) = .true.
            nst = nst + 1   ! stretching-torsion IC counter
            read(subword(1),*,err=2233)  iat_st(1)
            read(subword(2),*,err=2233)  jat_st(1)
            read(subword(3),*,err=2233)  iat_st(2)
            read(subword(4),*,err=2233)  jat_st(2)
            read(subword(5),*,err=2233)  kat_st(2)
            read(subword(6),*,err=2233)  lat_st(2)
            if (iat_st(1) .gt. jat_st(1)) then !! ctrl i < j str
                jnew      = iat_st(1)
                iat_st(1) = jat_st(1)
                jat_st(1) = jnew
            endif
            if (jat_st(2) .gt. kat_st(2)) then  !! ctrl j < k  tors
                lnew      = iat_st(2)
                knew      = jat_st(2)
                iat_st(2) = lat_st(2)
                jat_st(2) = kat_st(2)
                kat_st(2) = knew
                lat_st(2) = lnew
            endif
c           search for the coupled ICs among stretch. and torsions and stores
            do is = 1, nstretch
               ichk = abs(iat_st(1) - iat_s(is))
               jchk = abs(jat_st(1) - jat_s(is))
               nchk = ichk + jchk 
               if (nchk .eq. 0) then 
                iic_st(nst) = is
                zchk(1)     = .true.
               endif
            enddo
            read(subword(7),*,err=2233) cftype(3,nst)
            if (cftype(3,nst) .eq. 1) then   ! torsion proper (cosines)
             do it = 1, ndihedral
               ichk = abs(iat_st(2) - iat_d(it))
               jchk = abs(jat_st(2) - jat_d(it))
               kchk = abs(kat_st(2) - kat_d(it))
               lchk = abs(lat_st(2) - lat_d(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                jic_st(nst) = it
                zchk(2)     = .true.
               endif
             enddo
c            reads other parameters
             read(subword(8),*,err=2233)  rc_st(nst)
             read(subword(9),*,err=2233)  pc_st(nst,1) ! 1 = first cosine, later grouped
             read(subword(10),*,err=2233) k_st(nst,1)  ! 1 = first cosine, later grouped
             read(subword(11),*,err=2233) nc_st(nst,1) ! 1 = first cosine, later grouped
             do kk = 1,2
               if(.not. zchk(kk)) then 
                write(6,'(a)') 'Topology error: in str-tors couplings'
                write(6,'(a)') 'Coupled coordinate defined by atoms'
                write(6,'(4(i4,2x))') iat_st(kk), jat_st(kk), 
     &                                kat_st(kk), lat_st(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
               endif
             enddo
             goto 192
            elseif (cftype(3,nst) .eq. 2) then   ! torsion improper (harmonic)
             do it = 1, ndihedral_i
               ichk = abs(iat_st(2) - iat_di(it))
               jchk = abs(jat_st(2) - jat_di(it))
               kchk = abs(kat_st(2) - kat_di(it))
               lchk = abs(lat_st(2) - lat_di(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                jic_st(nst) = it
                zchk(2)     = .true.
               endif
             enddo
             read(subword(8),*,err=2233)  rc_st(nst)
             read(subword(9),*,err=2233)  pc_st(nst,1) 
             read(subword(10),*,err=2233) k_st(nst,1)
             do kk = 1,2
               if(.not. zchk(kk)) then 
                write(6,'(a)') 'Topology error: in str-tors couplings'
                write(6,'(a)') 'Coupled coordinate defined by atoms'
                write(6,'(4(i4,2x))') iat_st(kk), jat_st(kk), 
     &                                kat_st(kk), lat_st(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
               endif
             enddo
             goto 192
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'function type undefined in str-tors section'
             stop 
            endif
         endif


      ELSEIF( word(2) .eq. 'an_an') THEN       ! bending-bending

         nbb  = 0
 193     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 193
         if( index(string,';') .eq. 1) goto 193
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 193
         if(subword(1)(1:1) .eq. ';') goto 193
         if(msplit .eq. 0) then ! blank line
            nbendbend   = nbb
            ncoupl      = nstrstr + nstrbend + nstrtors + nbb
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nbendbend   = nbb
            ncoupl      = nstrstr + nstrbend + nstrtors + nbb
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 10) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'10 columns expected in bnd-bnd section'
             stop 
            endif
            z_coupl(4) = .true.
            nbb        = nbb + 1  ! bending-bending IC counter
            read(subword(1),*,err=2233)  iat_bb(1)
            read(subword(2),*,err=2233)  jat_bb(1)
            read(subword(3),*,err=2233)  kat_bb(1)
            read(subword(4),*,err=2233)  iat_bb(2)
            read(subword(5),*,err=2233)  jat_bb(2)
            read(subword(6),*,err=2233)  kat_bb(2)
            do kk = 1,2
              if (iat_bb(kk) .gt. kat_bb(kk)) then  !! ctrl i < k bnd
                 knew       = iat_bb(kk)
                 iat_bb(kk) = kat_bb(kk)
                 kat_bb(kk) = knew
              endif
            enddo
c           search for the coupled ICs among bendings and stores
            do ib = 1, nbend
               ichk = abs(iat_bb(1) - iat_b(ib))
               jchk = abs(jat_bb(1) - jat_b(ib))
               kchk = abs(kat_bb(1) - kat_b(ib))
               nchk = ichk + jchk  + kchk
               if (nchk .eq. 0) then 
                  iic_bb(nbb) = ib
                  zchk(1)     = .true.
               endif
               ichk = abs(iat_bb(2) - iat_b(ib))
               jchk = abs(jat_bb(2) - jat_b(ib))
               kchk = abs(kat_bb(2) - kat_b(ib))
               nchk = ichk + jchk  + kchk
               if (nchk .eq. 0) then 
                  jic_bb(nbb) = ib
                  zchk(2)     = .true.
               endif
            enddo
c           reads other parameters
            read(subword(7),*,err=2233) cftype(4,nbb)
            if (cftype(4,nbb) .eq. 1) then 
               read(subword(8),*,err=2233)  tc_bb(nbb,1)
               read(subword(9),*,err=2233)  tc_bb(nbb,2)
               read(subword(10),*,err=2233) k_bb(nbb)
               do kk = 1,2
                if(.not. zchk(kk)) then 
                 write(6,'(a)') 'Topology error: in bend-bend couplings'
                 write(6,'(a)') 'Coupled coordinate defined by atoms'
                 write(6,'(3(i4,2x))')iat_bb(kk), jat_bb(kk), kat_bb(kk)
                 write(6,'(a)') 'Not found in uncoupled list'
                stop
                endif
               enddo 
               goto 193
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'function type undefined in bnd-bnd section'
             stop 
            endif
         endif

      ELSEIF( word(2) .eq. 'an_di') THEN       ! bending-torsion

         nbt  = 0
 194     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 194
         if( index(string,';') .eq. 1) goto 194
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 194
         if(subword(1)(1:1) .eq. ';') goto 194
         if(msplit .eq. 0) then   ! blank line
            nbendtors   = nbt
            ncoupl      = nstrstr + nstrbend + nstrtors + 
     &                    nbendbend + nbt
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            nbendtors   = nbt
            ncoupl      = nstrstr + nstrbend + nstrtors + 
     &                    nbendbend + nbt
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 11) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')
     &         'at least 11 columns expected in bnd-tors section'
             stop 
            endif
            z_coupl(5) = .true.
            nbt        = nbt + 1  ! bending-torsion IC counter
            read(subword(1),*,err=2233)  iat_bt(1)
            read(subword(2),*,err=2233)  jat_bt(1)
            read(subword(3),*,err=2233)  kat_bt(1)
            read(subword(4),*,err=2233)  iat_bt(2)
            read(subword(5),*,err=2233)  jat_bt(2)
            read(subword(6),*,err=2233)  kat_bt(2)
            read(subword(7),*,err=2233)  lat_bt(2)
            if (iat_bt(1) .gt. kat_bt(1)) then  !! ctrl i < k bnd
               knew      = iat_bt(1)
               iat_bt(1) = kat_bt(1)
               kat_bt(1) = knew
            endif
            if (jat_bt(2) .gt. kat_bt(2)) then  !! ctrl j < k  tors
                lnew      = iat_bt(2)
                knew      = jat_bt(2)
                iat_bt(2) = lat_bt(2)
                jat_bt(2) = kat_bt(2)
                kat_bt(2) = knew
                lat_bt(2) = lnew
            endif
c           search for the coupled ICs among bendings and stores
            do ib = 1, nbend
               ichk = abs(iat_bt(1) - iat_b(ib))
               jchk = abs(jat_bt(1) - jat_b(ib))
               kchk = abs(kat_bt(1) - kat_b(ib))
               nchk = ichk + jchk  + kchk
               if (nchk .eq. 0) then 
                 iic_bt(nbt) = ib
                 zchk(1)     = .true.
               endif
            enddo
            read(subword(8),*,err=2233) cftype(5,nbt)
            if (cftype(5,nbt) .eq. 1) then   ! torsion proper (cosines)
             do it = 1, ndihedral
               ichk = abs(iat_bt(2) - iat_d(it))
               jchk = abs(jat_bt(2) - jat_d(it))
               kchk = abs(kat_bt(2) - kat_d(it))
               lchk = abs(lat_bt(2) - lat_d(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                jic_bt(nbt) = it
                zchk(2)     = .true.
               endif
             enddo
c            reads other parameters
             read(subword(9),*,err=2233)  tc_bt(nbt)
             read(subword(10),*,err=2233) pc_bt(nbt,1) ! 1 = first cosine, groups later
             read(subword(11),*,err=2233) k_bt(nbt,1)  ! 1 = first cosine, groups later
             read(subword(12),*,err=2233) nc_bt(nbt,1) ! 1 = first cosine, groups later
               do kk = 1,2
                if(.not. zchk(kk)) then 
                 write(6,'(a)') 'Topology error: in bend-tors couplings'
                 write(6,'(a)') 'Coupled coordinate defined by atoms'
                 write(6,'(4(i4,2x))') iat_bt(kk), jat_bt(kk), 
     &                                kat_bt(kk), lat_bt(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
                endif
               enddo
             goto 194
            elseif (cftype(5,nbt) .eq. 2) then   ! torsion improper (harmonic)
              do it = 1, ndihedral_i
              ichk = abs(iat_bt(2) - iat_di(it))
              jchk = abs(jat_bt(2) - jat_di(it))
              kchk = abs(kat_bt(2) - kat_di(it))
              lchk = abs(lat_bt(2) - lat_di(it))
              nchk = ichk + jchk  + kchk + lchk
              if (nchk .eq. 0) then 
                jic_bt(nbt) = it
                zchk(2)     = .true.
              endif
             enddo
c            reads other parameters
             read(subword(9),*,err=2233)  tc_bt(nbt)
             read(subword(10),*,err=2233) pc_bt(nbt,1) 
             read(subword(11),*,err=2233) k_bt(nbt,1)
             goto 194
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')'function type undefined in bnd-tors section'
             stop 
            endif
         endif


      ELSEIF( word(2) .eq. 'di_di') THEN       ! torsion-torsion

         ntt  = 0
 195     read(iun,1000) string
         do kk = 1,2
            zchk(kk) = .false.
         enddo
         if( index(string,'#') .eq. 1) goto 195
         if( index(string,';') .eq. 1) goto 195
         call setcomm(string)          ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 195
         if(subword(1)(1:1) .eq. ';') goto 195
         if(msplit .eq. 0) then  ! blank line
            ntorstors   = ntt
            ncoupl      = nstrstr + nstrbend + nstrtors + 
     &                    nbendbend + nbendtors + ntt
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            ntorstors   = ntt
            ncoupl      = nstrstr + nstrbend + nstrtors + 
     &                    nbendbend + nbendtors + ntt
            backspace(iun)
            goto 10
         else
            if (msplit .lt. 12) then
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')
     &         'at least 12 columns expected in tors-tors section'
             stop 
            endif
            z_coupl(6) = .true.
            ntt        = ntt + 1    ! torsion-torsion IC counter
            read(subword(1),*,err=2233)  iat_tt(1)
            read(subword(2),*,err=2233)  jat_tt(1)
            read(subword(3),*,err=2233)  kat_tt(1)
            read(subword(4),*,err=2233)  lat_tt(1)
            read(subword(5),*,err=2233)  iat_tt(2)
            read(subword(6),*,err=2233)  jat_tt(2)
            read(subword(7),*,err=2233)  kat_tt(2)
            read(subword(8),*,err=2233)  lat_tt(2)
            do kk = 1,2
               if (jat_tt(kk) .gt. kat_tt(kk)) then  !! ctrl j < k  tors
                   lnew       = iat_tt(kk)
                   knew       = jat_tt(kk)
                   iat_tt(kk) = lat_tt(kk)
                   jat_tt(kk) = kat_tt(kk)
                   kat_tt(kk) = knew
                   lat_tt(kk) = lnew
               endif
            enddo
c           search for the coupled ICs among torsions and stores
            read(subword(9),*,err=2233) cftype(6,ntt)
            if (cftype(6,ntt) .eq. 1) then   ! both torsion proper (cosines)
             do it = 1, ndihedral
               ichk = abs(iat_tt(1) - iat_d(it))
               jchk = abs(jat_tt(1) - jat_d(it))
               kchk = abs(kat_tt(1) - kat_d(it))
               lchk = abs(lat_tt(1) - lat_d(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                   iic_tt(ntt) = it
                   zchk(1)     = .true.
               endif
               ichk = abs(iat_tt(2) - iat_d(it))
               jchk = abs(jat_tt(2) - jat_d(it))
               kchk = abs(kat_tt(2) - kat_d(it))
               lchk = abs(lat_tt(2) - lat_d(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                   jic_tt(ntt) = it
                   zchk(2)     = .true.
               endif
             enddo
             read(subword(10),*,err=2233) pc_tt(ntt,1,1)! 1 = first cosine, groups later
             read(subword(11),*,err=2233) pc_tt(ntt,1,2)! 1 = first cosine, groups later
             read(subword(12),*,err=2233) k_tt(ntt,1)! 1 = first cosine, groups later
             read(subword(13),*,err=2233) nc_tt(ntt,1,1)! 1 = first cosine, groups later
             read(subword(14),*,err=2233) nc_tt(ntt,1,2)! 1 = first cosine, groups later
             do kk = 1,2
              if(.not. zchk(kk)) then 
                write(6,'(a)') 'Topology error: in tors-tors couplings'
                write(6,'(a)') 'Coupled coordinate defined by atoms'
                write(6,'(4(i4,2x))') iat_tt(kk), jat_tt(kk), 
     &                                kat_tt(kk), lat_tt(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
              endif
             enddo
             goto 195
            elseif (cftype(6,ntt) .eq. 2) then   ! one torsion proper and one improper 
             do it = 1, ndihedral_i
               ichk = abs(iat_tt(1) - iat_di(it))
               jchk = abs(jat_tt(1) - jat_di(it))
               kchk = abs(kat_tt(1) - kat_di(it))
               lchk = abs(lat_tt(1) - lat_di(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                iic_tt(ntt) = it
                zchk(1)     = .true.
               endif
             enddo
             do it = 1, ndihedral
               ichk = abs(iat_tt(2) - iat_d(it))
               jchk = abs(jat_tt(2) - jat_d(it))
               kchk = abs(kat_tt(2) - kat_d(it))
               lchk = abs(lat_tt(2) - lat_d(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                jic_tt(ntt) = it
                zchk(2)     = .true.
               endif
             enddo
             read(subword(10),*,err=2233) pc_tt(ntt,1,1)
             read(subword(11),*,err=2233) pc_tt(ntt,1,2)
             read(subword(12),*,err=2233) k_tt(ntt,1)! 1 = first cosine, groups later
             read(subword(13),*,err=2233) nc_tt(ntt,1,2)
             do kk = 1,2
              if(.not. zchk(kk)) then 
                write(6,'(a)') 'Topology error: in tors-tors couplings'
                write(6,'(a)') 'Coupled coordinate defined by atoms'
                write(6,'(4(i4,2x))') iat_tt(kk), jat_tt(kk), 
     &                                kat_tt(kk), lat_tt(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
              endif
             enddo
             goto 195
            elseif (cftype(6,ntt) .eq. 3) then   ! both harmonic torsion 
             do it = 1, ndihedral_i
               ichk = abs(iat_tt(1) - iat_di(it))
               jchk = abs(jat_tt(1) - jat_di(it))
               kchk = abs(kat_tt(1) - kat_di(it))
               lchk = abs(lat_tt(1) - lat_di(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                iic_tt(ntt) = it
                zchk(1)     = .true.
               endif
               ichk = abs(iat_tt(2) - iat_di(it))
               jchk = abs(jat_tt(2) - jat_di(it))
               kchk = abs(kat_tt(2) - kat_di(it))
               lchk = abs(lat_tt(2) - lat_di(it))
               nchk = ichk + jchk  + kchk + lchk
               if (nchk .eq. 0) then 
                jic_tt(ntt) = it
                zchk(2)     = .true.
               endif
             enddo
             read(subword(10),*,err=2233) pc_tt(ntt,1,1)
             read(subword(11),*,err=2233) pc_tt(ntt,1,2)
             read(subword(12),*,err=2233) k_tt(ntt,1)
             do kk = 1,2
              if(.not. zchk(kk)) then 
                write(6,'(a)') 'Topology error: in tors-tors couplings'
                write(6,'(a)') 'Coupled coordinate defined by atoms'
                write(6,'(4(i4,2x))') iat_tt(kk), jat_tt(kk), 
     &                                kat_tt(kk), lat_tt(kk) 
                write(6,'(a)') 'Not found in uncoupled list'
                stop
              endif
             enddo
             goto 195
            else
             write(6,'(a)')'Please check the .top file'
             write(6,'(a)')
     &        'function type undefined in tors-tors section'
             stop 
            endif
         endif


      ELSEIF( word(2) .eq. 'atoms') THEN   

         iat  = 0
 16      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 16
         if( index(string,';') .eq. 1) goto 16
         call setcomm(string)     ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 16
         if(subword(1)(1:1) .eq. ';') goto 16
         if(msplit .eq. 0)         then  ! blank line
            Nato = iat
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            Nato = iat
            backspace(iun)
            goto 10
         else
            iat = iat + 1
            read(subword(2),*,err=2233) name(iat)
            read(subword(7),*,err=2233) charge(iat)
            goto 16
         endif            

      ELSEIF( word(2) .eq. 'atomtypes') THEN   
              

         isit  = 0
 17      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 17
         if( index(string,';') .eq. 1) goto 17
         call setcomm(string)     ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 17
         if(subword(1)(1:1) .eq. ';') goto 17
         if(msplit .eq. 0)  then
            Nsit = isit
            goto 10
         elseif(index(string,'[') .ne. 0)   then  ! new section
            Nsit = isit
            backspace(iun)
            goto 10
         else
            ! Determine the format based on the position of PType entry
            ! which can take the value A, S, V or D
            if(subword(4) .eq. 'A' .or. 
     &         subword(4) .eq. 'S' .or.
     &         subword(4) .eq. 'V' .or.
     &         subword(4) .eq. 'D') then
               iname = 1 
               imass = 2
               isig  = 5
               ieps  = 6
            elseif (subword(5) .eq. 'A' .or.  
     &              subword(5) .eq. 'S' .or.
     &              subword(5) .eq. 'V' .or.
     &              subword(5) .eq. 'D') then
               iname = 1 
               imass = 3
               isig  = 6
               ieps  = 7
            elseif (subword(6) .eq. 'A' .or.  
     &              subword(6) .eq. 'S' .or.
     &              subword(6) .eq. 'V' .or.
     &              subword(6) .eq. 'D') then
               iname = 1 ! or 2 (BondType)
               imass = 4
               isig  = 7
               ieps  = 8
            else
               goto 2233
            endif

            isit = isit + 1
            read(subword(iname),*,err=2233) namesit(isit)
            read(subword(imass),*,err=2233) amass(isit)
            read(subword(isig ),*,err=2233) sig_sit(isit)
            read(subword(ieps ),*,err=2233) eps_sit(isit)

            goto 17
         endif            

      ELSEIF( word(2) .eq. 'moleculetype') THEN   

 117      read(iun,1000) string
          if (index(string,';') .eq. 1)  goto 117
          goto 10 

      ELSEIF( word(2) .eq. 'exclusions') THEN   

 118     read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 118
         if( index(string,';') .eq. 1) goto 118
         call setcomm(string)     ! eliminates comments after ";"
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#')      goto 118
         if(subword(1)(1:1) .eq. ';') goto 118
         if(msplit .eq. 0)  then
            goto 118
         elseif(index(string,'[') .ne. 0)   then  ! new section
            backspace(iun)
            goto 10
         elseif (msplit .eq. 2) then 
            read(subword(1),*,iostat=nchk1)  ichk
            read(subword(2),*,iostat=nchk2)  jchk
            nchk = nchk1+nchk2
            if (nchk .eq. 0) then 
                goto 118
            else
               card = string
               goto 119
            endif
         else
            card = string
            goto 119
         endif

      ELSEIF( word(2) .eq. 'system') THEN   

        goto 100

      ELSE

        write(6,'(a)') 'Syntax error found in starting .top file'
        write(6,'(a)') 'Please check the following line:'
        write(6,'(a)') card
        stop 'syntax error'

      ENDIF
      goto 10

119   continue
      write(6,'(a)') 'Syntax error found in starting .top file'
      write(6,'(a)') 'Please check the following line:'
      write(6,'(a)') card
      stop 'syntax error'
 100  continue

c     obsolete to be eventually deleted (old Gromacs couplings)  ! GP 9/2017
c     if( z_coupl(1) .or. z_coupl(2) .or. z_coupl(3)  
c    &  .or. z_coupl(4) .or. z_coupl(5) .or. z_coupl(6) ) then 
c        if ( z_coupl(7) ) then
c          write(6,'(a)')'Please check the .top file.'
c          write(6,'(a)')'Coupling must be either in the standard'
c          write(6,'(a)')'or in the Joyce format, not both.'
c          stop 
c        endif
c     endif

C     Groups cosine terms in proper torsions 
      
      do it = 1,Ndihedral-1
       do jt = it+1,Ndihedral
        ichk = abs(iat_d(it) - iat_d(jt))
        jchk = abs(jat_d(it) - jat_d(jt))
        kchk = abs(kat_d(it) - kat_d(jt))
        lchk = abs(lat_d(it) - lat_d(jt))
        nchk = ichk + jchk + kchk + lchk
        if( nchk .eq. 0) then
           do jcos = 1,ic(jt)
              icnew = ic(it) + jcos
              k_d(it,icnew)  = k_d(jt,jcos)
              n0_d(it,icnew) = n0_d(jt,jcos)
              p0_d(it,icnew) = p0_d(jt,jcos)
           enddo   
           ic(it) = ic(it) + ic(jt)
           do lt = jt,Ndihedral-1
              iat_d(lt) = iat_d(lt+1)
              jat_d(lt) = jat_d(lt+1)
              kat_d(lt) = kat_d(lt+1)
              lat_d(lt) = lat_d(lt+1)
              ic(lt) = ic(lt+1)
              do lcos = 1,ic(lt)
                 k_d(lt,lcos)  = k_d(lt+1,lcos)
                 n0_d(lt,lcos) = n0_d(lt+1,lcos)
                 p0_d(lt,lcos) = p0_d(lt+1,lcos)
              enddo
           enddo
           Ndihedral = Ndihedral -1
           goto 100
        endif
       enddo
      enddo

C     Groups cosine terms in proper torsions  within couplings 

c     stretching-torsion
      if (z_coupl(3)) then
       icos = 1
101    continue
       do it = 1,Nstrtors-1
        if (cftype(3,it) .eq. 1) then 
         do jt = it+1,Nstrtors
          if (iic_st(it) .eq. iic_st(jt)) then  ! same bond coupled IC
            if(jic_st(it) .eq. jic_st(jt)) then ! same dihed coupled IC => removes from count
               if (rc_st(it) .ne. rc_st(jt)) then 
                write(6,'(a)') 'Topology error: in str-tors couplings'
                write(6,'(a)') 'Rc was found to have different values'
                write(6,'(2(f8.2,2x))') rc_st(it), rc_st(jt)
                write(6,'(a)') 
     &                   'when the bond is coupled to the same dihedral'
                stop
               endif
               if (ncos_st(it) .gt. 1) then
                   icos =  ncos_st(it) + 1
               else
                   icos =  icos + 1
               endif
               nc_st(it,icos) = nc_st(jt,1)
               pc_st(it,icos) = pc_st(jt,1)
               k_st(it,icos)  = k_st(jt,1)
               do lt = jt,Nstrtors-1
                  iic_st(lt)  = iic_st(lt+1)
                  jic_st(lt)  = jic_st(lt+1)
                  pc_st(lt,1) = pc_st(lt+1,1)
                  nc_st(lt,1) = nc_st(lt+1,1)
                  rc_st(lt)   = rc_st(lt+1)
                  k_st(lt,1)  = k_st(lt+1,1)
                  cftype(3,lt)= cftype(3,lt+1)
               enddo
               Nstrtors     = Nstrtors -1
               Ncos_st(it)  = icos
               goto 101
            endif
          endif 
         enddo ! jt
         icos = 1
        endif  ! if proper
       enddo   ! it
      endif

c     bending-torsion
      if (z_coupl(5)) then
       icos = 1
102    continue
       do it = 1,Nbendtors-1
        if (cftype(5,it) .eq. 1) then 
         do jt = it+1,Nbendtors
          if (iic_bt(it) .eq. iic_bt(jt)) then  ! same angle coupled IC
           if( jic_bt(it) .eq. jic_bt(jt)) then ! same dih. coupled IC => removes form count
               if (tc_bt(jt) .ne. tc_bt(it)) then 
                write(6,'(a)') 'Topology error: in bend-tors couplings'
                write(6,'(a)') 'Tc was found to have different values'
                write(6,'(2(f8.2,1x))') tc_bt(it), tc_bt(jt)
                write(6,'(a)') 
     &                 'when the angle is coupled to the same dihedral'
                stop
               endif
               if (Ncos_bt(it) .gt. 1) then
                   icos =  Ncos_bt(it) + 1
               else
                   icos =  icos + 1
               endif
               nc_bt(it,icos) = nc_bt(jt,1)
               pc_bt(it,icos) = pc_bt(jt,1)
               k_bt(it,icos)  = k_bt(jt,1)
               do lt = jt,Nbendtors-1
                  iic_bt(lt)  = iic_bt(lt+1)
                  jic_bt(lt)  = jic_bt(lt+1)
                  pc_bt(lt,1) =  pc_bt(lt+1,1)
                  nc_bt(lt,1) =  nc_bt(lt+1,1)
                  tc_bt(lt)   =  tc_bt(lt+1)
                  k_bt(lt,1)  =  k_bt(lt+1,1)
                  cftype(5,lt)=  cftype(5,lt+1)
               enddo
               Nbendtors   = Nbendtors -1
               Ncos_bt(it) = icos
               goto 102
            endif
          endif
         enddo ! jt
         icos = 1
        endif  ! if proper
       enddo ! it
      endif

c     torsion-torsion
      if (z_coupl(6)) then
       icos = 1
103    continue
       do it = 1,Ntorstors-1
        if (cftype(6,it) .eq. 1) then   ! two torsion proper (cosines)
         do jt = it+1,Ntorstors
          if(iic_tt(it) .eq. iic_tt(jt)) then  ! same 1st dihedral coupled IC
           if(jic_tt(it) .eq. jic_tt(jt)) then ! same 2nd dih. coupl. IC => remove from count
              if (Ncos_tt(it) .gt. 1) then
                 icos =  Ncos_tt(it) + 1
              else
                 icos =  icos + 1
              endif
              do kk = 1,2
                   nc_tt(it,icos,kk) = nc_tt(jt,1,kk)
                   pc_tt(it,icos,kk) = pc_tt(jt,1,kk)
              enddo 
              k_tt(it,icos)  = k_tt(jt,1)
              do lt = jt,Ntorstors-1
                 iic_tt(lt)  = iic_tt(lt+1)
                 jic_tt(lt)  = jic_tt(lt+1)
                 do kk= 1,2
                   pc_tt(lt,1,kk) =  pc_tt(lt+1,1,kk)
                   nc_tt(lt,1,kk) =  nc_tt(lt+1,1,kk)
                 enddo 
                 k_tt(lt,1)  =  k_tt(lt+1,1)
                 cftype(6,lt)=  cftype(6,lt+1)
              enddo
                Ntorstors   = Ntorstors -1
                Ncos_tt(it) = icos
                goto 103
           endif
          endif
         enddo ! jt
         icos = 1
        elseif (cftype(6,it) .eq. 2) then   ! 1 torsion proper (cosines)
         do jt = it+1,Ntorstors
          if(iic_tt(it) .eq. iic_tt(jt)) then  ! same first harmonic dihedral coupled IC
           if(jic_tt(it) .eq. jic_tt(jt)) then ! same 2nd dih. coupl. IC => remove from count
            if (pc_tt(it,1,1) .ne. pc_tt(jt,1,1)) then 
              write(6,'(a)') 'Topology error: in tors-tors couplings'
              write(6,'(a)') 'Phic was found to have different values'
              write(6,'(2(f8.2,1x))') pc_tt(it,1,1), pc_tt(jt,1,1)
              write(6,'(a)') 'when the harmonic dihedral is coupled'
              write(6,'(a)') 'to the same proper dihedral'
              stop
            endif
            if (Ncos_tt(it) .gt. 1) then
                icos =  ncos_tt(it) + 1
            else
                icos =  icos + 1
            endif
            do kk = 1,2
               pc_tt(it,icos,kk) = pc_tt(jt,1,kk)
               nc_tt(it,icos,kk) = nc_tt(jt,1,kk)
            enddo 
            k_tt(it,icos)    = k_tt(jt,1)
            do lt = jt,Ntorstors-1
               iic_tt(lt)  = iic_tt(lt+1)
               jic_tt(lt)  = jic_tt(lt+1)
               do kk= 1,2
                  pc_tt(lt,1,kk) =  pc_tt(lt+1,1,kk)
                  nc_tt(lt,1,kk) =  nc_tt(lt+1,1,kk)
               enddo 
               k_tt(lt,1)    =  k_tt(lt+1,1)
               cftype(6,lt)  =  cftype(6,lt+1)
            enddo
            Ntorstors = Ntorstors -1
            Ncos_tt(it) = icos
            goto 103
           endif
          endif
         enddo! jt
         icos = 1
        endif ! harmonic/proper
       enddo  ! it
      endif


C     sets defaults for intramolecular nonbonded


      do iat = 1,Nato
         labchk = name(iat)
         call compact(labchk,lw)
         do isit =1,Nsit
            call compact(namesit(isit),lw3)
            if(namesit(isit)(1:lw3) .eq. labchk(1:lw)) then
               stemp(iat) = sig_sit(isit)
               etemp(iat) = eps_sit(isit)
            endif
         enddo  
      enddo  
      do iat = 1,Nato
         sigma(iat)   = stemp(iat)
         epsilon(iat) = etemp(iat)
      enddo  

C     Unit Conversion to Joyce   
C     1 kJ/mol       = 1/2625.5d0 Hartree        
C     1 nm           = 1/0.052917726d0 Bohr    
C     1 kJ/(mol*nm^2)= 1/(2625.d0*(0.0529117726d0)^2) = 1.066*10^(-6) H/(Bohr^2)
C     1 kJ/(mol*nm)  = 1/(0.052917726d0*2625.d0) H/Bohr
C     e^2/4*pi*eps0  = 1 H bohr --> no conversion needed      


      do is=1,Nstretch
         k_s(is)  = k_s(is)*boh2nm *boh2nm /H2kJ
         r0_s(is) = r0_s(is)/boh2nm
      enddo

      do ib=1,Nbend
         k_b(ib)  = k_b(ib)/H2kJ
         t0_b(ib) = t0_b(ib)/rad2deg
      enddo
      
      do id = 1,Ndihedral
         do it = 1,ic(id)
            k_d(id,it)  = k_d(id,it)/H2kJ
            p0_d(id,it) = p0_d(id,it)/rad2deg
         enddo
      enddo

      do it = 1,Ndihedral_i
         k_di(it)  = k_di(it)/H2kJ
         p0_di(it) = p0_di(it)/rad2deg
      enddo

      do lj = 1,Nlj
         eps_lj(lj)= eps_lj(lj)/H2kJ
         sig_lj(lj)= sig_lj(lj)/boh2nm 
      enddo

      ntt = Nstretch + Nbend + Ndihedral + Ndihedral_i
      ntl = ntt + Nlj

      if (z_coupl(1)) then 
         do icc = 1,Nstrstr
            k_ss(icc)  = k_ss(icc)*boh2nm *boh2nm /H2kJ
            do kk = 1,2
               rc_ss(icc,kk) = rc_ss(icc,kk)/boh2nm
            enddo
         enddo
      endif

      if (z_coupl(2)) then 
         do icc = 1,Nstrbend
            k_sb(icc)  = k_sb(icc)*boh2nm /H2kJ
            rc_sb(icc) = rc_sb(icc)/boh2nm
            tc_sb(icc) = tc_sb(icc)/rad2deg
         enddo
      endif

      if (z_coupl(3)) then 
         do icc = 1,Nstrtors
          rc_st(icc)      = rc_st(icc)/boh2nm
          do icos = 1,Ncos_st(icc)
            k_st(icc,icos)  = k_st(icc,icos)*boh2nm /H2kJ
            pc_st(icc,icos) = pc_st(icc,icos)/rad2deg
          enddo
         enddo
      endif

      if (z_coupl(4)) then 
         do icc = 1,Nbendbend
            k_bb(icc)  = k_bb(icc)/H2kJ
            do kk = 1,2
               tc_bb(icc,kk) = tc_bb(icc,kk)/rad2deg
            enddo
         enddo
      endif

      if (z_coupl(5)) then 
         do icc = 1,Nbendtors
          tc_bt(icc)  = tc_bt(icc)/rad2deg
          do icos = 1,Ncos_bt(icc)
            k_bt(icc,icos)  = k_bt(icc,icos)/H2kJ
            pc_bt(icc,icos) = pc_bt(icc,icos)/rad2deg
          enddo
         enddo
      endif

      if (z_coupl(6)) then 
         do icc = 1,Ntorstors
          do icos = 1,Ncos_tt(icc)
            k_tt(icc,icos)  = k_tt(icc,icos)/H2kJ
            do kk = 1,2
              pc_tt(icc,icos,kk) = pc_tt(icc,icos,kk)/rad2deg
            enddo
          enddo
         enddo
      endif

c     puts data on iflag array and parameters in par vector      
c     the order counts! 
c     1) Diagonal terms 2) LJ  3) Couplings into iflag vector
      
      none = 0
      col  = ';'

      int=0                 !! parameter index
c     iflag(int,0) = 1  ----> stretching
      do is=1,Nstretch
        int=int+1
        iflag(int,0)=1
        par(int)    = k_s(is)
        iflag(int,1)= iat_s(is)
        iflag(int,2)= jat_s(is)
        iflag(int,7)= is
        iflag(int,8)= 2
        parfun(int) = r0_s(is)
        parfun2(int)= 0.d0
        card=' '
        write(card,'(a,"-",a)') aname(iat_s(is)),aname(jat_s(is))
        Call Compact (card,lw)
        labpar(int)=card
      enddo

c     iflag(int,0) = 2  ----> bending
      do ib=1,Nbend
        int=int+1
        iflag(int,0)=2
        par(int)    = k_b(ib)
        iflag(int,1)= iat_b(ib)
        iflag(int,2)= jat_b(ib)
        iflag(int,3)= kat_b(ib)
        iflag(int,7)= ib
        iflag(int,8)= 3
        parfun(int) = t0_b(ib)
        parfun2(int)= 0.d0
        card=' '
        write(card,'(a,"-",a,"-",a)') aname(iat_b(ib))
     &               ,aname(jat_b(ib)),aname(kat_b(ib))
        Call Compact (card,lw)
        labpar(int)=card
      enddo

c     iflag(int,0) = 3  ----> harmonic torsion
      do it = 1,Ndihedral_i
        int=int+1
        iflag(int,0)=3
        par(int)    =k_di(it)
        iflag(int,1)=iat_di(it)
        iflag(int,2)=jat_di(it)
        iflag(int,3)=kat_di(it)
        iflag(int,4)=lat_di(it)
        iflag(int,8)= 4
        parfun(int) = p0_di(it)
        parfun2(int)= 0.d0
        card=' '
        write(card,'(a,"-",a,"-",a,"-",a)') aname(iat_di(it))
     &   ,aname(jat_di(it)),aname(kat_di(it)),aname(lat_di(it))
        Call Compact (card,lw)
        labpar(int)=card
      enddo

c     iflag(int,0) = 4  ----> proper torsion
      do id = 1,Ndihedral
         do it = 1,ic(id)
            int=int+1
            iflag(int,0)=4
            par(int)    = k_d(id,it)
            iflag(int,1)= iat_d(id)
            iflag(int,2)= jat_d(id)
            iflag(int,3)= kat_d(id)
            iflag(int,4)= lat_d(id)
            iflag(int,5)= n0_d(id,it)
            iflag(int,6)=ic(id)
            iflag(int,8)= 4
            parfun(int) = p0_d(id,it)
            parfun2(int)= 0.d0
         enddo
         card=' '
         write(card,'(a,"-",a,"-",a,"-",a)') aname(iat_d(id))
     &     ,aname(jat_d(id)),aname(kat_d(id)),aname(lat_d(id))
         Call Compact (card,lw)
         labpar(int)=card
      enddo


c     iflag(int,0) = 5  ----> intramolecular LJ 
      do lj = 1,Nlj
         int=int+1
         if (int .gt. maxint)  goto 2255
         iflag(int,0)=5
         par(int)    = eps_lj(lj)
         iflag(int,1)= iat_lj(lj)
         iflag(int,2)= jat_lj(lj)
         iflag(int,8)= 2
         parfun(int) = sig_lj(lj)
         parfun2(int) =sig_lj(lj)
         card=' '
         write(card,'(a,"-",a)') aname(iat_lj(lj)),aname(jat_lj(lj))
         Call Compact (card,lw)
         labpar(int)=card
      enddo

c     iflag(int,0) = 6  ----> charge-charge
      if (z_intcharge) then
         Nchch = Nlj
         ichrgprm(1) = int + 1
         do lj = 1,Nchch
            int=int+1
            if (int .gt. maxint)  goto 2255
            iflag(int,0)=6
            par(int)    = qq(lj)
            iflag(int,1)= iat_lj(lj)
            iflag(int,2)= jat_lj(lj)
            iflag(int,8)= 2
            parfun(int) = 0.d0
            parfun2(int)= 0.d0
            card=' '
            write(card,'(a,"-",a)') aname(iat_lj(lj)),aname(jat_lj(lj))
            Call Compact (card,lw)
            labpar(int)=card
         enddo
         ichrgprm(2)= int
      endif


c     iflag(int,0) = 8  ----> couplings
c     couplings in Joyce modified .top format   ! udpdated 9/2017 GP

      if (z_coupl(1)) then   !! stretch-stretch
          do iss=1,Nstrstr   
             int = int + 1
             iflag(int,0) = 8
             if(iflag(iic_ss(iss),0).eq.1.and.
     &          iflag(jic_ss(iss),0).eq.1) then
                par(int)     = k_ss(iss)
                parfun(int)  = rc_ss(iss,1) 
                parfun2(int) = rc_ss(iss,2) 
                iflag(int,1) = iic_ss(iss)
                iflag(int,2) = jic_ss(iss)
                iflag(int,8) = 4
              else
                write(6,'(a)') 'Topology error: check  str-str pairs'
                stop
              endif
              card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,'(a,"-",a,")(",a,"-",a)')
     &          aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &          aname(iflag(ic2,1)), aname(iflag(ic2,2))  
              Call Compact (card,lw)
              labpar(int)=card
          enddo
      endif

      if (z_coupl(2)) then            !! stretch-bend
          do isb=1,Nstrbend   
             int = int + 1
             iflag(int,0) = 8
             par(int)     = k_sb(isb)
             istr         = iic_sb(isb)
             ibnd         = jic_sb(isb) + Nstretch
             if(iflag(istr,0).eq.1.and.iflag(ibnd,0).eq.2) then
                parfun(int)  = rc_sb(isb)
                parfun2(int) = tc_sb(isb)
                iflag(int,1) = istr
                iflag(int,2) = ibnd 
                iflag(int,8) = 5
             else
              write(6,'(a)') 'Topology error: check str-bnd pairs'
              stop
             endif
             card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,'(a,"-",a,")(",a,"-",a,"-",a)')
     &           aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &           aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &           aname(iflag(ic2,3)) 
             Call Compact (card,lw)
             labpar(int)=card
          enddo
      endif

      if (z_coupl(3)) then
          do ist=1,Nstrtors   !! stretch-tors
           if (cftype(3,ist) .eq. 2) then  ! harmonic torsion
             int          = int + 1
             par(int)     = k_st(ist,1)
             iflag(int,0) = 8
             istr         = iic_st(ist)
             itor         = jic_st(ist) + Nstretch + Nbend
             if(iflag(istr,0).eq.1.and.iflag(itor,0).eq.3) then
                parfun(int)  = rc_st(ist)
                parfun2(int) = pc_st(ist,1)
                iflag(int,1) = istr
                iflag(int,2) = itor
                iflag(int,8) = 6
             else
              write(6,'(a)')
     &            'Topology error: check str-harmonic tors pairs'
              stop
             endif
             card=' '
             ic1 = iflag(int,1)   !! ic 1 to be coupled
             ic2 = iflag(int,2)   !! ic 2 to be coupled
             write(card,'(a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &         aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &         aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &         aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
             Call Compact (card,lw)
             labpar(int)=card
           elseif (cftype(3,ist) .eq. 1) then  ! proper torsion
             do icos=1,Ncos_st(ist) 
              int          = int + 1
              iflag(int,0) = 8
              par(int)     = k_st(ist,icos)
              istr         = iic_st(ist)
              itor         = jic_st(ist) +Nstretch+Nbend+Ndihedral_i
              iflag(int,6) = Ncos_st(ist)
              if(iflag(istr,0).eq.1.and.iflag(itor,0).eq.4) then
                 parfun(int)  = rc_st(ist)
                 parfun2(int) = pc_st(ist,icos)
                 iflag(int,1) = istr
                 iflag(int,2) = itor
                 iflag(int,5) = nc_st(ist,icos)
                 iflag(int,8) = 6
              else
               write(6,'(a)')
     &             'Topology error: check str-proper tors pairs'
               stop
              endif
              card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,'(a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &          aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &          aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &          aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
              Call Compact (card,lw)
              labpar(int)=card
             enddo
           endif
          enddo
      endif

      if (z_coupl(4)) then     !! bend-bend
          do ibb=1,Nbendbend   
              int = int + 1
              iflag(int,0) = 8
              par(int)     = k_bb(ibb)
              ibnd1        = iic_bb(ibb) + Nstretch
              ibnd2        = jic_bb(ibb) + Nstretch
              if(iflag(ibnd1,0).eq.2.and.iflag(ibnd2,0).eq.2) then
                 parfun(int)  = tc_bb(ibb,1)
                 parfun2(int) = tc_bb(ibb,2)
                 iflag(int,1) = ibnd1
                 iflag(int,2) = ibnd2
                 iflag(int,8) = 6
              else
              write(6,'(a)') 'Topology error: check Bend-bend pairs'
              stop
              endif
              card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,'(a,"-",a,"-",a,")(",a,"-",a,"-",a)')
     &          aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &          aname(iflag(ic1,3)), 
     &          aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &          aname(iflag(ic2,3)) 
              Call Compact (card,lw)
              labpar(int)=card
          enddo
      endif


      if (z_coupl(5)) then      !! bend-tors
          do ibt=1,Nbendtors   
           if (cftype(5,ibt) .eq. 2) then  ! harmonic torsion
              int          = int + 1
              iflag(int,0) = 8
              par(int)     = k_bt(ibt,1)
              ibnd         = iic_bt(ibt) + Nstretch
              itor         = jic_bt(ibt) + Nstretch + Nbend
              if(iflag(ibnd,0).eq.2.and.iflag(itor,0).eq.3) then
                 parfun(int)  = tc_bt(ibt)
                 parfun2(int) = pc_bt(ibt,1)
                 iflag(int,1) = ibnd
                 iflag(int,2) = itor
                 iflag(int,8) = 7
              else
               write(6,'(a)') 'Topology error: check bend-tors pairs'
               stop
              endif
              card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,'(a,"-",a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &          aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &          aname(iflag(ic1,3)), 
     &          aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &          aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
              Call Compact (card,lw)
              labpar(int)=card
           elseif (cftype(5,ibt) .eq. 1) then  ! proper torsion
              do icos=1,Ncos_bt(ibt)
                 int          = int + 1
                 iflag(int,0) = 8
                 par(int)     = k_bt(ibt,icos)
                 ibnd         = iic_bt(ibt) + Nstretch
                 itor         = jic_bt(ibt) +Nstretch+Nbend+Ndihedral_i
                 if(iflag(ibnd,0).eq.2.and.iflag(itor,0).eq.4) then
                    parfun(int)  = tc_bt(ibt)
                    parfun2(int) = pc_bt(ibt,icos)
                    iflag(int,1) = ibnd
                    iflag(int,2) = itor
                    iflag(int,5) = nc_bt(ibt,icos)
                    iflag(int,8) = 7
                 else
                  write(6,'(a)') 'Topology error: check bend-tors pairs'
                  stop
                 endif
                 card=' '
                 ic1 = iflag(int,1)   !! ic 1 to be coupled
                 ic2 = iflag(int,2)   !! ic 2 to be coupled
                 write(card,'(a,"-",a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &             aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &             aname(iflag(ic1,3)), 
     &             aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &             aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
                 Call Compact (card,lw)
                 labpar(int)=card
              enddo
           endif
          enddo
      endif

      if (z_coupl(6)) then
          do itt=1,Ntorstors   !! tors-tors coupling
           if (cftype(6,itt) .eq. 3) then  ! two harmonic torsions
              int          = int + 1
              iflag(int,0) = 8
              par(int)     = k_tt(itt,1)
              itor1        = iic_tt(itt) + Nstretch + Nbend
              itor2        = jic_tt(itt) + Nstretch + Nbend
              if(iflag(itor1,0).eq.3.and.iflag(itor2,0).eq.3) then
                   parfun(int)  = pc_tt(itt,1,1)
                   parfun2(int) = pc_tt(itt,1,2)
                   iflag(int,1) = itor1
                   iflag(int,2) = itor2
                   iflag(int,8) = 8
              else
               write(6,'(a)') 
     &         'Topology error: check harmonic tors-tors pairs'
               stop
              endif
              card=' '
              ic1 = iflag(int,1)   !! ic 1 to be coupled
              ic2 = iflag(int,2)   !! ic 2 to be coupled
              write(card,
     &         '(a,"-",a,"-",a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &         aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &         aname(iflag(ic1,3)), aname(iflag(ic1,4)),
     &         aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &         aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
              Call Compact (card,lw)
              labpar(int)=card
           elseif (cftype(6,itt) .eq. 2) then  ! one harmonic and one proper torsion
              do icos=1,Ncos_tt(itt)
                 int          = int + 1
                 iflag(int,0) = 8
                 par(int)     = k_tt(itt,icos)
                 itor1        = iic_tt(itt) +Nstretch+Nbend
                 itor2        = jic_tt(itt) +Nstretch+Nbend+Ndihedral_i
                 if(iflag(itor1,0).eq.3.and.iflag(itor2,0).eq.4) then
                      parfun(int)  = pc_tt(itt,1,1)
                      parfun2(int) = pc_tt(itt,icos,2)
                      iflag(int,1) = itor1
                      iflag(int,2) = itor2
                      iflag(int,5) = nc_tt(itt,icos,2)
                      iflag(int,8) = 8
                 else
                  write(6,'(a)') 
     &            'Topology error: check harm.tors- proper tors pairs'
                  stop
                 endif
                 card=' '
                 ic1 = iflag(int,1)   !! ic 1 to be coupled
                 ic2 = iflag(int,2)   !! ic 2 to be coupled
                 write(card,
     &            '(a,"-",a,"-",a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &            aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &            aname(iflag(ic1,3)), aname(iflag(ic1,4)),
     &            aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &            aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
                 Call Compact (card,lw)
                 labpar(int)=card
              enddo
           elseif (cftype(6,itt) .eq. 1) then  ! two proper torsions
              do icos=1,Ncos_tt(itt)
                 int          = int + 1
                 iflag(int,0) = 8
                 par(int)     = k_tt(itt,icos)
                 itor1        = iic_tt(itt) +Nstretch+Nbend+Ndihedral_i
                 itor2        = jic_tt(itt) +Nstretch+Nbend+Ndihedral_i
                 if(iflag(itor1,0).eq.4.and.iflag(itor2,0).eq.4) then
                      parfun(int)  = pc_tt(itt,icos,1)
                      parfun2(int) = pc_tt(itt,icos,2)
                      iflag(int,1) = itor1
                      iflag(int,2) = itor2
                      iflag(int,5) = nc_tt(itt,icos,1)
                      iflag(int,6) = nc_tt(itt,icos,2)
                      iflag(int,8) = 8
                 else
                  write(6,'(a)') 
     &            'Topology error: check proper tors-tors pairs'
                  stop
                 endif
                 card=' '
                 ic1 = iflag(int,1)   !! ic 1 to be coupled
                 ic2 = iflag(int,2)   !! ic 2 to be coupled
                 write(card,
     &            '(a,"-",a,"-",a,"-",a,")(",a,"-",a,"-",a,"-",a)')
     &            aname(iflag(ic1,1)), aname(iflag(ic1,2)), 
     &            aname(iflag(ic1,3)), aname(iflag(ic1,4)),
     &            aname(iflag(ic2,1)), aname(iflag(ic2,2)), 
     &            aname(iflag(ic2,3)), aname(iflag(ic2,4)) 
                 Call Compact (card,lw)
                 labpar(int)=card
              enddo
           endif
          enddo
      endif


      Npar = int  !! Npar including charges and couplings

c     Prints out read potential terms  

        write(6,*)
        write(6,'(22x,a)') '1) Atom Types'
        write(6,*)
        write(6,1001) 'Site','Name','Charge','Mass','Sigma','Epsilon'
        do isit = 1,Nsit
           write(6,1002) isit,namesit(isit),charge(isit),amass(isit),
     &                       sig_sit(isit),eps_sit(isit)
        enddo
        write(6,1000) 

c     ! now prints directly from iflag vector ! GP+JC 9/2017

      iprt = 1
      if (Nstretch .gt. 0) then
       do  ip = 1, Npar
         if (iprt .eq. 1) then 
          write(6,1000) 
          write(6,2000) '2) Stretching parameters'
          write(6,1000) 
          write(6,1003) 'Bond','Atoms','k_s','r0','FF term','Atoms'
          iprt = 0
          ns   = 0
         endif
         if (iflag(ip,0) .eq. 1) then 
            ns       = ns + 1
            iref(ip) = ns
            i1       = iflag(ip,1) 
            i2       = iflag(ip,2) 
            r0       = parfun(ip)*boh2nm*10.d0 
            c0       = par(ip)*H2kJ/(boh2nm *boh2nm*100.d0)
            write(6,1004) ns,name(i1),name(i2),c0,r0,ip,i1,i2
         endif
       enddo
       if (ns .ne. Nstretch) then 
         write(6,*) 'Expected number of stretchings:',Nstretch
         write(6,*) 'Internal counter:',ns
         stop 'internal error: check stretchings'
       endif
       write(6,1000) 
      endif

      iprt = 1
      if (Nbend .gt. 0) then
       write(6,1000) 
       do  ip = 1, Npar
         if (iprt .eq. 1) then 
          write(6,2000) '3) Bending parameters'
          write(6,1000) 
          write(6,1005) 'Angle','Atoms','k_b','theta0','FF term','Atoms'
          iprt = 0
          nb   = 0
         endif
         if (iflag(ip,0) .eq. 2) then 
            nb       = nb + 1
            iref(ip) = nb
            i1       = iflag(ip,1) 
            i2       = iflag(ip,2) 
            i3       = iflag(ip,3) 
            t0       = parfun(ip)*rad2deg
            c0       = par(ip)*H2kJ
            write(6,1006)nb,name(i1),name(i2),name(i3),c0,t0,ip,i1,i2,i3
         endif
       enddo
       if (nb .ne. Nbend) then 
         write(6,*) 'Expected number of bendings:',Nbend
         write(6,*) 'Internal counter:',nb
         stop 'internal error: check bendings'
       endif
       write(6,1000) 
      endif

      iprt = 1
      if (Ndihedral_i .gt. 0) then
       write(6,1000) 
       do  ip = 1, Npar
         if (iprt .eq. 1) then 
          write(6,2000) '4) Torsional parameters'
          write(6,1000) 
          write(6,2000) '4.1 Harmonic torsions'
          write(6,1000) 
          write(6,1007) 'Dihedral','Atoms','    k_d','    gamma',
     &                  'FF term', 'Atoms'
          iprt = 0
          nt   = 0
         endif
         if (iflag(ip,0) .eq. 3) then 
            nt = nt + 1
            iref(ip) = nt
            i1       = iflag(ip,1) 
            i2       = iflag(ip,2) 
            i3       = iflag(ip,3) 
            i4       = iflag(ip,4) 
            p0       = parfun(ip)*rad2deg
            c0       = par(ip)*H2kJ
            write(6,1008) nt,name(i1),name(i2),name(i3),name(i4),
     &                    c0,p0,ip,i1,i2,i3,i4
         endif
       enddo
       if (nt .ne. Ndihedral_i) then 
         write(6,*) 'Expected number of harmonic torsions:',Ndihedral_i
         write(6,*) 'Internal counter:',nt
         stop 'internal error: check harmonic torsions'
       endif
       write(6,1000) 
      endif

      iprt = 1
      if (Ndihedral .gt. 0) then
       do  ip = 1, Npar
         if (iprt .eq. 1) then 
          write(6,1000) 
          write(6,2000) '4.2 Fourier torsions'
          write(6,1000) 
          write(6,1009)
          iprt = 0
          nt   = 0
          i10  = 0
          i20  = 0
          i30  = 0
          i40  = 0
         endif
         if (iflag(ip,0) .eq. 4) then 
            i1   = iflag(ip,1) 
            i2   = iflag(ip,2) 
            i3   = iflag(ip,3) 
            i4   = iflag(ip,4) 
            nc   = iflag(ip,5) 
            nct  = iflag(ip,6) 
            ga   = parfun(ip)*rad2deg
            c0   = par(ip)*H2kJ
            if (i1 .eq. i10 .and. i2 .eq. i20 .and.
     &          i3 .eq. i30 .and. i4 .eq. i40) then 
              icos = icos + 1
              iref(ip) = nt
              write(6,1011) c0,nc,ga,ip
            else
              icos     = 1
              nt       = nt + 1
              iref(ip) = nt
              write(6,1010) nt,name(i1),name(i2),name(i3),name(i4),
     &                      nct,c0,nc,ga,ip,i1,i2,i3,i4
              i10  = i1
              i20  = i2
              i30  = i3
              i40  = i4
            endif
         endif
       enddo
       if (nt .ne. Ndihedral) then 
         write(6,*) 'Expected number of proper torsions:',Ndihedral
         write(6,*) 'Internal counter:',nt
         stop 'internal error: check proper torsions'
       endif
       write(6,1000) 
      endif


c     LJ terms if any
      ipr  = 1
      if (z_lj) then
         write(6,2000) '5) Nonbonded parameters'
         write(6,1000) 
         do ip = 1, Npar
          if (iflag(ip,0) .eq. 5) then   ! lj
            ilj1 = iflag(ip,1)
            ilj2 = iflag(ip,2)
            ep   = par(ip)*H2kJ
            si   = parfun(ip)*boh2nm*10.d0
            if (z_intcharge) then
              do ip2 = 1, Npar
                 if (iflag(ip2,0) .eq. 6) then 
                   iq1 = iflag(ip2,1)
                   iq2 = iflag(ip2,2)
                   ch2 = par(ip2)
                   if (iq1 .eq. ilj1 .and. iq2 .eq. ilj2 ) then 
                    if (ipr .eq. 1) then 
                        write(6,1014) 'Couple','Atoms','sigma (A)',
     &                 'epsilon (kJ/mol)','f*q*q (e^2)','FF terms'
                        write(6,1015) ilj1,ilj2,name(ilj1),name(ilj2),
     &                               si,ep,ch2,ip,',',ip+1
                          ipr = 0
                    else
                        write(6,1015) ilj1,ilj2,name(ilj1),name(ilj2),
     &                                si,ep,ch2,ip,',',ip+1
                    endif
                   endif
                 endif
              enddo
            else
              if (ipr .eq. 1) then 
                 write(6,1014) 'Couple','Atoms','sigma (A)',
     &                         'epsilon (kJ/mol)','FF term'
                 write(6,1016) ilj1,ilj2,name(ilj1),name(ilj2),
     &                         si,ep,ip
                 ipr = 0
              else
                        write(6,1016) ilj1,ilj2,name(ilj1),name(ilj2),
     &                               si,ep,ip
              endif
            endif
          endif
         enddo
       endif 
 
c     couplings if any 
      do kc = 1,6
         if (z_coupl(kc)) then 
           write(6,1000) 
           write(6,1000) 
           write(6,'(37x,a)') '6) Couplings'
           goto 206
         endif
      enddo
 206  continue
c 
      if (z_coupl(1)) then
        write(6,1000) 
        write(6,2000) '6.1 Stretching-stretching coupling'
        write(6,1000) 
        write(6,1032) 'Bon','Bon','Atoms','k_ss','rc1','rc2',
     &                'Coupled terms','FF term' 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 1 .and. iflag(ic2,0) .eq. 1) then 
              ia1 = iflag(ic1,1)
              ia2 = iflag(ic1,2)
              ib1 = iflag(ic2,1)
              ib2 = iflag(ic2,2)
              c0  = par(ip)*H2kJ/(boh2nm*boh2nm*100.d0) 
              r1e = parfun(ip)*boh2nm*10.d0
              r2e = parfun2(ip)*boh2nm*10.d0
              ir1 = iref(ic1)
              ir2 = iref(ic2)
              write(6,1033) ir1,ir2,ia1,ia2,ib1,ib2,c0,r1e,r2e,
     &                      ic1,ic2,ip
           endif
          endif
        enddo
      endif
 

      if (z_coupl(2)) then
        write(6,1000) 
        write(6,1000) 
        write(6,2000) '6.2 Stretching-bending coupling'
        write(6,1000) 
        write(6,1034) 'Bon','Ang','Atoms','k_sb','rc1','tc2', 
     &                'Coupled terms','FF term' 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 1 .and. iflag(ic2,0) .eq. 2) then 
              ia1 = iflag(ic1,1)
              ia2 = iflag(ic1,2)
              ib1 = iflag(ic2,1)
              ib2 = iflag(ic2,2)
              ib3 = iflag(ic2,3)
              c0  = par(ip)*H2kJ/(boh2nm*10.d0) 
              r1e = parfun(ip)*boh2nm*10.d0
              t2e = parfun2(ip)*rad2deg
              ir1 = iref(ic1)
              ir2 = iref(ic2)
              write(6,1035) ir1,ir2,ia1,ia2,ib1,ib2,ib3,c0,r1e,t2e,
     &                      ic1,ic2,ip
           endif
          endif
        enddo
      endif
  
      
      if (z_coupl(3)) then
        write(6,1000) 
        write(6,1000) 
        write(6,2000) '6.3 Stretching-torsion coupling'
        write(6,1000) 
        ipr  = 1
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 1 .and. iflag(ic2,0) .eq. 3) then 
            if (ipr .eq. 1) then 
             write(6,2000) '6.3a Harmonic torsions'
             write(6,1000) 
             write(6,1036) 'Bon','Dih','Atoms','k_st','rc1','phi_c',
     &                     'Coupled terms','FF term'
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ/(boh2nm*10.d0) 
            r1  = parfun(ip)*boh2nm*10.d0
            p1  = parfun2(ip)*rad2deg
            ir1 = iref(ic1)
            ir2 = iref(ic2)
            write(6,1037)ir1,ir2,ia1,ia2,ib1,ib2,ib3,ib4,c0,r1,p1,
     &                   ic1,ic2,ip
           endif
         endif
        enddo
        ipr  = 1
        write(6,1000) 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 1 .and. iflag(ic2,0) .eq. 4) then 
            if (ipr .eq. 1) then 
             write(6,1000) '6.3b Proper torsions'
             write(6,1000) 
             write(6,1038)'Bon','Dih','Atoms','k_st','rc1','gamma','n',
     &                    'Coupled terms','FF term'
             ia10 = 0
             ia20 = 0
             ib10 = 0
             ib20 = 0
             ib30 = 0
             ib40 = 0
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ/(boh2nm*10.d0) 
            r1  = parfun(ip)*boh2nm*10.d0
            ga  = parfun2(ip)*rad2deg
            mm  = iflag(ip,5)
            if(ia1 .eq. ia10 .and. ia2 .eq. ia20 .and.
     &         ib1 .eq. ib10 .and. ib2 .eq. ib20 .and.
     &         ib3 .eq. ib30 .and. ib4 .eq. ib40 ) then
             write(6,1039) c0, r1, ga, mm, ip
            else
             ia10 = ia1
             ia20 = ia2
             ib10 = ib1
             ib20 = ib2
             ib30 = ib3
             ib40 = ib4
             ir1  = iref(ic1)
             ir2  = iref(ic2)
             write(6,1040)ir1,ir2,ia1,ia2,ib1,ib2,ib3,ib4,c0,r1,ga,mm,
     &                    ic1,ic2,ip
            endif
           endif
          endif
        enddo
      endif
 
      if (z_coupl(4)) then
        write(6,1000) 
        write(6,1000) '6.4 Bending-bending coupling'
        write(6,1000) 
        write(6,1041) 'Ang','Ang','Atoms','k_bb','tc1','tc2', 
     &                'Coupled terms','FF term' 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 2 .and. iflag(ic2,0) .eq. 2) then 
               ia1 = iflag(ic1,1)
               ia2 = iflag(ic1,2)
               ia3 = iflag(ic1,3)
               ib1 = iflag(ic2,1)
               ib2 = iflag(ic2,2)
               ib3 = iflag(ic2,3)
               c0  = par(ip)*H2kJ
               t1e = parfun(ip)*rad2deg
               t2e = parfun2(ip)*rad2deg
               ir1 = iref(ic1)
               ir2 = iref(ic2)
               write(6,1042) ir1,ir2,ia1,ia2,ia3,ib1,ib2,ib3,c0,
     &                       t1e,t2e,ic1,ic2,ip
           endif
          endif
        enddo
      endif
 
      if (z_coupl(5)) then
        write(6,1000) 
        write(6,1000) 
        write(6,1000) '6.5 Bending-torsion coupling'
        write(6,1000) 
        ipr  = 1
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 2 .and. iflag(ic2,0) .eq. 3) then 
            if (ipr .eq. 1) then 
             write(6,1000) '5.5a Harmonic torsions'
             write(6,1000) 
             write(6,1043) 'Ang','Dih','Atoms','k_bt','tc1','phi_c',
     &                     'Coupled terms','FF term'
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ia3 = iflag(ic1,3)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ
            t1  = parfun(ip)*rad2deg
            p1  = parfun2(ip)*rad2deg
            ir1 = iref(ic1)
            ir2 = iref(ic2)
            write(6,1044)ir1,ir2,ia1,ia2,ia3,ib1,ib2,ib3,ib4,c0,t1,p1,
     &                   ic1,ic2,ip
           endif
         endif
        enddo
        ipr  = 1
        write(6,1000) 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 2 .and. iflag(ic2,0) .eq. 4) then 
            if (ipr .eq. 1) then 
             write(6,1000) '5.5b Proper torsions'
             write(6,1000) 
             write(6,1045)'Ang','Dih','Atoms','k_bt','tc1','gamma','n',
     &                    'Coupled terms','FF term'
             ia10 = 0
             ia20 = 0
             ia30 = 0
             ib10 = 0
             ib20 = 0
             ib30 = 0
             ib40 = 0
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ia3 = iflag(ic1,3)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ
            t1  = parfun(ip)*rad2deg
            ga  = parfun2(ip)*rad2deg
            mm  = iflag(ip,5)
            if(ia1 .eq. ia10 .and. ia2 .eq. ia20 .and.
     &         ia3 .eq. ia30 .and. ib1 .eq. ib10 .and. 
     &         ib2 .eq. ib20 .and. ib3 .eq. ib30 .and. 
     &         ib4 .eq. ib40 ) then
             write(6,1047) c0, t1, ga, mm, ip
            else
             ia10 = ia1
             ia20 = ia2
             ia30 = ia3
             ib10 = ib1
             ib20 = ib2
             ib30 = ib3
             ib40 = ib4
             ir1  = iref(ic1)
             ir2  = iref(ic2)
             write(6,1046)ir1,ir2,ia1,ia2,ia3,ib1,ib2,ib3,ib4,c0,t1,ga,
     &                    mm,ic1,ic2, ip
            endif
           endif
          endif
        enddo
      endif

      if (z_coupl(6)) then
        write(6,1000) 
        write(6,1000) 
        write(6,1000) '6.6 Torsion-torsion coupling'
        write(6,1000) 
        ipr  = 1
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 3 .and. iflag(ic2,0) .eq. 3) then 
            if (ipr .eq. 1) then 
             write(6,1000) '6.5a Two harmonic torsions'
             write(6,1000) 
             write(6,1048) 'Dih','Dih','Atoms','k_tt','phi1','phi2',
     &                     'Coupled terms','FF term'
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ia3 = iflag(ic1,3)
            ia4 = iflag(ic1,4)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ
            p1  = parfun(ip)*rad2deg
            p2  = parfun2(ip)*rad2deg
            ir1 = iref(ic1)
            ir2 = iref(ic2)
            write(6,1049)ir1,ir2,ia1,ia2,ia3,ia4,ib1,ib2,ib3,ib4,
     &                   c0,p1,p2,ic1,ic2,ip
           endif
         endif
        enddo
        ipr  = 1
        write(6,1000) 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 3 .and. iflag(ic2,0) .eq. 4) then 
            if (ipr .eq. 1) then 
             write(6,1000) '6.5b One harmonic and one proper torsion'
             write(6,1000) 
             write(6,1050)'Dih','Dih','Atoms','k_tt','phi','gamma',
     &                    'n','Coupled terms','FF term'
             ia10 = 0
             ia20 = 0
             ia30 = 0
             ia40 = 0
             ib10 = 0
             ib20 = 0
             ib30 = 0
             ib40 = 0
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ia3 = iflag(ic1,3)
            ia4 = iflag(ic1,4)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ
            p1  = parfun(ip)*rad2deg
            ga  = parfun2(ip)*rad2deg
            mm  = iflag(ip,5)
            if(ia1 .eq. ia10 .and. ia2 .eq. ia20 .and.
     &         ia3 .eq. ia30 .and. ia4 .eq. ia40 .and.
     &         ib1 .eq. ib10 .and. ib2 .eq. ib20 .and. 
     &         ib3 .eq. ib30 .and. ib4 .eq. ib40 ) then
             write(6,1051) c0, p1, ga, mm, ip
            else
             ia10 = ia1
             ia20 = ia2
             ia30 = ia3
             ia40 = ia4
             ib10 = ib1
             ib20 = ib2
             ib30 = ib3
             ib40 = ib4
             ir1  = iref(ic1)
             ir2  = iref(ic2)
             write(6,1052)ir1,ir2,ia1,ia2,ia3,ia4,ib1,ib2,ib3,ib4,
     &                    c0,p1,ga,mm,ic1,ic2,ip 
            endif
           endif
          endif
        enddo
        ipr  = 1
        write(6,1000) 
        do  ip = 1, Npar
          if (iflag(ip,0) .eq. 8) then 
           ic1 = iflag(ip,1)
           ic2 = iflag(ip,2)
           if (iflag(ic1,0) .eq. 4 .and. iflag(ic2,0) .eq. 4) then 
            if (ipr .eq. 1) then 
             write(6,1000) '6.5c Two proper torsions'
             write(6,1000) 
             write(6,1053)'Dih','Dih','Atoms','k_tt','ga1','ga2',
     &                    'n1','n2','Coupled terms', 'FF term'
             ia10 = 0
             ia20 = 0
             ia30 = 0
             ia40 = 0
             ib10 = 0
             ib20 = 0
             ib30 = 0
             ib40 = 0
             ipr  = 0
            endif
            ia1 = iflag(ic1,1)
            ia2 = iflag(ic1,2)
            ia3 = iflag(ic1,3)
            ia4 = iflag(ic1,4)
            ib1 = iflag(ic2,1)
            ib2 = iflag(ic2,2)
            ib3 = iflag(ic2,3)
            ib4 = iflag(ic2,4)
            c0  = par(ip)*H2kJ
            ga1 = parfun(ip)*rad2deg
            ga2 = parfun2(ip)*rad2deg
            m1  = iflag(ip,5)
            m2  = iflag(ip,6)
            if(ia1 .eq. ia10 .and. ia2 .eq. ia20 .and.
     &         ia3 .eq. ia30 .and. ia4 .eq. ia40 .and.
     &         ib1 .eq. ib10 .and. ib2 .eq. ib20 .and. 
     &         ib3 .eq. ib30 .and. ib4 .eq. ib40 ) then
             write(6,1054) c0, ga1, ga2, m1, m2, ip
            else
             ia10 = ia1
             ia20 = ia2
             ia30 = ia3
             ia40 = ia4
             ib10 = ib1
             ib20 = ib2
             ib30 = ib3
             ib40 = ib4
             ir1  = iref(ic1)
             ir2  = iref(ic2)
             write(6,1055)ir1,ir2,ia1,ia2,ia3,ia4,ib1,ib2,ib3,ib4,
     &                    c0,ga1,ga2,m1,m2,ic1,ic2, ip
            endif
           endif
          endif
        enddo
      endif

 
 
c     write(6,1000) 
c     write(6,1000) 
c     write(6,1000)  
c    & '       ********************************************'
c     write(6,1000)  
c    & '          E N D  G R O M A C S  I N P U T'
c     write(6,1000)  
c    & '       ********************************************'
c     write(6,1000) 
c     write(6,1000) 

c     formats

  555 format(i3,2x,7i3,3f9.3)
 1000 format(a)
 2000 format(10x,a)
c     LJ formats
 1001 format(10x,a4,2x,a4,3x,2(4x,a6,2x),6x,a5,9x,a7)
 1002 format(10x,i2,5x,a3,3x,5(f10.3,3x))
c     stretching formats
 1003 format(11x,a4,6x,a5,12x,a3,8x,a2,3x,a7,5x,a5)
 1004 format(11x,i3,5x,2(a3,3x),4x,f9.2,2x,f7.3,1x,i4,6x,2(i3,2x))
c     bending formats
 1005 format(11x,a5,7x,a5,10x,a3,5x,a6,3x,a7,7x,a5)
 1006 format(12x,i3,4x,3(a3,2x),f9.2,2x,f7.3,3x,i4,6x,3(i3,2x))
c     harmonic torsion formats
 1007 format(10x,a8,7x,a5,7x,a7,4x,a9,2x,a7,7x,a5)
 1008 format(13x,i2,6x,4(a3,1x),f9.2,3x,f8.2,3x,i4,3x,4(i3,1x))
c     proper torsion formats
 1009 format(10x,' Dihedral',8x,'Atoms',6x,'Ncos   K_d'
     &              ,'     n    gamma',3x,'FF term',5x,'Atoms')
 1010 format(13x,i3,6x,4(a3,1x),i2,f10.4,i4,f9.2,3x,i4,3x,4(i3,1x))
 1011 format(40x,f10.4,i4,f9.2,3x,i4,3x,4(i3,1x))
c     non bonded parameters
 1014 format(14x,a6,7x,a5,11x,a9,3x,a16,3x,a12,3x,a8)
 1015 format(12x,i3,1x,'-',i3,5x,2(a3,3x),4x,f8.2,5x,f9.3,10x,f9.3,6x,
     &       i4,a1,i4)
 1016 format(12x,i3,1x,'-',i3,5x,2(a3,3x),4x,f8.2,5x,f9.3,10x,i4)
c     coupling related formats
c     1) stretching-stretching
 1032 format(11x,a3,1x,'-',1x,a3,6x,a6,11x,a4,2x,
     &       2(3x,a3,3x),1x,a13,3x,a7)
 1033 format(11x,i3,1x,'-',i3,2x,4(i3,1x),3x,f9.2,2(f7.2,1x),
     &       4x,2(i3,1x),6x,i4)
c     2) stretching-bending
 1034 format(11x,a3,1x,'-',1x,a3,12x,a6,14x,a4,2x,2(3x,a3,1x),5x,a13,
     &       3x,a7)
 1035 format(11x,i3,1x,'-',i3,6x,5(i3,1x),3x,f9.2,2(f7.2,2x),
     &       4x,2(i3,1x),6x,i4)
c     3) stretching-torsion
 1036 format(11x,a3,1x,'-',1x,a3,11x,a5,11x,a6,4x,a4,4x,a5,4x,a13,3x,
     &       a7)
 1037 format(11x,i3,1x,'-',i3,2x,6(i3,1x),f8.2,1x,f8.2,1x,f7.1,2x,4x,
     &       2(i3,1x),6x,i4)
 1038 format(11x,a3,1x,'-',1x,a3,11x,a5,11x,a6,4x,a4,4x,a5,4x,
     &       a1,5x,a13,3x,a7)
 1039 format(44x,f9.2,1x,f8.2,f8.1,3x,i2,21x,i4)
 1040 format(11x,i3,1x,'-',i3,2x,6(i3,1x),f8.2,1x,f8.2,1x,f7.1,3x,
     &       i2,4x,2(i3,1x),6x,i7)
c     4) bending-bending
 1041 format(11x,a3,1x,'-',1x,a3,7x,a6,16x,a4,4x,2(3x,a3,5x),a13,
     &       3x,a7)
 1042 format(11x,i3,1x,'-',i3,1x,6(i3,1x),1x,f9.2,3x,2(f8.2,3x),
     &       1x,2(i3,1x),6x,i4)
c     5) bending-torsion
 1043 format(11x,a3,1x,'-',1x,a3,11x,a5,15x,a6,4x,a4,4x,a5,4x,a13,
     &       3x,a7)
 1044 format(11x,i3,1x,'-',i3,2x,7(i3,1x),f8.2,1x,f8.2,1x,f7.1,
     &       7x,2(i3,1x),6x,i4)
 1045 format(11x,a3,1x,'-',1x,a3,11x,a5,15x,a6,4x,a4,4x,a5,3x,a1,4x,a13,
     &       3x,a7)
 1046 format(11x,i3,1x,'-',i3,2x,7(i3,1x),f8.2,1x,f8.2,1x,f7.1,
     &       2x,i2,5x,2(i3,1x),6x,i4)
 1047 format(48x,f9.2,1x,f8.2,f8.1,2x,i2,19x,i4)
c     6) torsion--torsion
 1048 format(11x,a3,1x,'-',1x,a3,13x,a5,15x,a6,4x,a4,4x,a5,4x,a13,
     &       3x,a7)
 1049 format(11x,i3,1x,'-',i3,1x,8(i3,1x),f8.2,1x,f8.2,1x,f7.1,
     &       5x,2(i3,1x),6x,i4)
 1050 format(11x,a3,1x,'-',1x,a3,13x,a5,15x,a6,4x,a4,4x,a5,4x,a1,3x,a13,
     &       3x,a7)
 1051 format(51x,f9.2,2x,f7.2,f7.1,4x,i3,21x,i4)
 1052 format(11x,i3,1x,'-',i3,1x,8(i3,1x),f7.2,1x,f8.1,f7.1,
     &       4x,i2,4x,2(i3,1x),6x,i4)
 1053 format(11x,a3,1x,'-',1x,a3,13x,a5,15x,a6,5x,a3,4x,a3,3x,a2,
     &       2x,a2,3x,a13,3x,a7)
 1054 format(52x,f8.2,1x,f7.1,f7.1,3x,i2,2x,i2,19x,i4)
 1055 format(11x,i3,1x,'-',i3,1x,8(i3,1x),f8.2,1x,f7.1,f7.1,
     &       3x,i2,2x,i2,5x,2(i3,1x),6x,i4)



      return

c     I/O errors

2232  write(6,1000) 'I/O error: check your .system file'      
      write(6,1000) 'Last correct line:'
      write(6,2238) iline,':',card0
2238  format(i4,1x,a1,1x,a90)      
      stop

2233  write(6,1000) 'I/O error: check your .top file'      
      write(6,1000) 'Syntax error on line:'
      write(6,1000) string
      stop

2255  write(6,1000) 'Dimension error: too many int. coords'
      write(6,2256) 'Number of coordinates requested:  ',int
      write(6,2256) 'Number of max coordinates allowed:',maxint
2256  format(a37,1x,i7) 
      stop
      
      
      end

