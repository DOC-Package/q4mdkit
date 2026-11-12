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
      Subroutine Read_Moscito (iun,iflag,par,parfun,npar,geom,nat
     &                        ,maxpar,maxat,aname,labpar,name,ichrgprm
     &                        ,nsit,namesit,cr,eps_sit,sig_sit
     &                        ,amass,charge,sigma,epsilon,z_lj)
c----------------------------------------------------------------------c

 
c         reads a Moscito .system file as input and
c         extract an 'iflag' matrix with the requested potential
c         function according to Joyce conventions
 

c----------------------------------------------------------------------c
c         *****        JOYCE conventions on iflag  matrix  ******
c----------------------------------------------------------------------c
c                                                                      c
c         iflag(j,k)      k-th label of the j-th parameter             c
c                                                                      c
c         in particular, for k = 0,                                    c
c                                                                      c
c         iflag(ipar,0) identifies parameter type i.e.:                c
c         iflag(ipar,0)  = 1     stretching                            c
c         iflag(ipar,0)  = 2     bending                               c
c         iflag(ipar,0)  = 3     torsione (harmonic)                   c
c         iflag(ipar,0)  = 4     torsione proper  (with cosines)       c
c         iflag(ipar,0)  = 5     LJ intramolecular                     c
c         iflag(ipar,0)  = 6     charge-charge intramolecular          c
c         iflag(ipar,0)  = 8     coupling                              c
c                                                                      c
c         for k = 1,2,...                                              c
c                                                                      c
c         iflag(ipar,1)  = i     atom  i which is needed               c
c                                for the definition                    c
c                                of the internal coordinate            c
c                                contained in iflag(ipar,0)            c
c         iflag(ipar,2)  = j     atom  j which is needed               c
c                                for the definition                    c
c                                of the internal coordinate            c
c                                contained in iflag(ipar,0)            c
c         iflag(ipar,3)  = k     atom  k (if present) which is needed  c
c                                for the definition of the i.c.        c
c                                contained in iflag(ipar,0)            c
c         iflag(ipar,4)  = l     atom  l (if present) which is needed  c
c                                for the definition of the i.c.        c
c                                contained in iflag(ipar,0)            c
c                                                                      c
c         Note : i < j < k < l                                         c
c                                                                      c
c         for k = 5                                                    c
c                                                                      c
c         iflag(ipar,5)  = n     only if iflag(ipar,0) = 4             c
c                                is the n value for the cosine term    c
c                                Vn = 1 + cos(n*phi - gamma)           c
c         iflag(ipar,6)  = N     only if iflag(ipar,0) = 4             c
c                                is the number of cosines functions    c
c                                                                      c
c         if iflag(ipar,0)  = 8                                        c
c                                                                      c
c         iflag(ipar,1)  = i     i-th coupled function                 c
c         iflag(ipar,2)  = j     j-th coupled function                 c
c                                                                      c
c                                                                      c
c                                                                      c
c         *****        Potential model functions ******                c
c                                                                      c
c         Stretching : Vs  = Ks  * (r     - r0)^2                      c
c         Bending    : Vb  = Kb  * (theta - theta0)^2                  c
c         Torsione i.: Vdi = Kdi * (phi   - phi0  )^2                  c
c         Torsione   : V0 + V1 + V2 + .....                            c
c                      con                                             c
c                      Vn  = Kn * [1 + cos(n*phi - gamma)]             c
c                                                                      c
c         Couplings                                                    c
c                                                                      c
c         Str - str  : Vss = Kss * (r     - r01)     * (r     - r02)   c
c         Bnd - bnd  : Vbb = Kbb * (theta - theta01)*(theta - theta02) c
c         Str - bnd  : Vsb = Ksb * (r     - r0) * (theta - theta0)     c
c         Str - tors : Vst = Kst * (r     - r0) * [ Vc0 + Vc1 + ...]   c
c                      con                                             c
c                      Vcn = Kcn * [sin(n*phi + gamma)]                c
c         Bnd - tors : Vbt = Kbt * (theta - theta0)*[Vc0 + Vc1 + ...]  c
c                      con                                             c
c                      Vcn = Kcn * [sin(n*phi + gamma)]                c
c         Tors- tors : Vtt0 + Vtt1 + Vtt2 + .....                      c
c                      con                                             c
c                      Vttn = kttn * {sin(n1*phi1 - gamma1)*           c
c                             sin(n2*phi2 - gamma2)}                   c
c                                                                      c
c     PARFUN                                                           c
c                                                                      c
c     iflag(i,0)=1  nt=2   harmonic stretching            R0           c
c     iflag(i,0)=2  nt=3   harmonic bending               Ang0         c
c     iflag(i,0)=3  nt=4   harmonic torsion               Phi0         c
c     iflag(i,0)=4  nt=4   torsione proper  (cosines)     Gamma        c
c     WARNING: parfun is in radians or bohr  !!!                       c
c     iflag(i,0)=5  nt=2   LJ intramolecular              R0           c
c     iflag(i,0)=6  nt=2   QQ intramolecular              R0           c
c     WARNING: parfun is in radians or bohr  !!!                       c
c     iflag(i,0)=8         coupling                        --          c
c                                                                      c
c----------------------------------------------------------------------cc

 
      implicit real*8 (a-h,o-z)

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      
      integer          pat_bd,pp,pnew,pat_dd,rat_dd
      integer          cr
      real*8           k_s,k_b,k_bl,k_d,k_di,k_sb,k_bd
      real*8           k_ss, k_bb, k_sd, k_dd
      real*8           j2kjmol
      character        card*90,string*90,card0*90
      character*46     word(maxspl),subword(maxspl)
      character        name(maxat)*5,namesit(*)*5,col*1
      character        namesit2(maxat)*5
      character        labchk1*5,labchk2*5
      logical          z_harmonic,z_linear,z_improper,z_lj,z_coupl
      logical          z_c,z_test,z_intcharge
 
      data todeg/57.2957795130823d0/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs



      dimension x(maxint),y(maxint),z(maxint)
      dimension iat_s(maxint),jat_s(maxint),k_s(maxint),r0_s(maxint)
      dimension iat_b(maxint),jat_b(maxint),kat_b(maxint),k_b(maxint),
     &          t0_b(maxint)
      dimension iat_bl(maxint),jat_bl(maxint),kat_bl(maxint),
     &          k_bl(maxint),t0_bl(maxint)
      dimension iat_d(maxint),jat_d(maxint),kat_d(maxint),lat_d(maxint),
     &          k_d(maxint,maxcos),n0_d(maxint,maxcos),
     &          p0_d(maxint,maxcos),ic(maxint)
      dimension iat_di(maxint),jat_di(maxint),kat_di(maxint),
     &          lat_di(maxint),k_di(maxint),p0_di(maxint)
      dimension z_coupl(6)
      dimension iat_ss(maxint),jat_ss(maxint),kat_ss(maxint),
     &          lat_ss(maxint),
     &          k_ss(maxint),r1_ss(maxint),r2_ss(maxint)
      dimension iat_bb(maxint),jat_bb(maxint),kat_bb(maxint),
     &          lat_bb(maxint),mat_bb(maxint),nat_bb(maxint),
     &          k_bb(maxint),t1_bb(maxint),t2_bb(maxint)
      dimension iat_sb(maxint),jat_sb(maxint),kat_sb(maxint),
     &          lat_sb(maxint),mat_sb(maxint),
     &          k_sb(maxint),r0_sb(maxint),t0_sb(maxint)
      dimension iat_sd(maxint),jat_sd(maxint),kat_sd(maxint),
     &          lat_sd(maxint),mat_sd(maxint),nat_sd(maxint),
     &          k_sd(maxint),r0_sd(maxint),n0_sd(maxint),g_sd(maxint)
      dimension iat_bd(maxint),jat_bd(maxint),kat_bd(maxint),
     &          lat_bd(maxint),mat_bd(maxint),nat_bd(maxint),
     &          pat_bd(maxint),
     &          k_bd(maxint),t0_bd(maxint),n0_bd(maxint),g_bd(maxint)
      dimension iat_dd(maxint),jat_dd(maxint),kat_dd(maxint),
     &          lat_dd(maxint),mat_dd(maxint),nat_dd(maxint),
     &          pat_dd(maxint),rat_dd(maxint),
     &          k_dd(maxint),n1_dd(maxint),n2_dd(maxint),
     &          g1_dd(maxint),g2_dd(maxint)
      dimension iat_lj(maxint),jat_lj(maxint),
     &          eps_lj(maxint),sig_lj(maxint),fqq(maxint),qq(maxint)
      dimension ichrgprm(2)
      dimension etemp(maxat),stemp(maxat)

      character aname(*)*5,labpar(*)*(*)
      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension par(*),parfun(*),charge(*),amass(*)
      dimension sigma(*),epsilon(*)
      dimension sig_sit(*),eps_sit(*)

      do i=1,maxpar
        par(i)    = 0.d0
        parfun(i) = 0.d0
        do j=0,8
          iflag(i,j)=0
        enddo
      enddo
      do i = 1,2
         ichrgprm(i) = 0
      enddo

      z_test     = .true.           !! test on geometry

      z_harmonic = .false.
      z_linear   = .false.
      z_improper = .false.
      z_lj       = .false.
      z_intcharge= .false.
      z_c        = .false.
      do k = 1,6
         z_coupl(k) = .false.
      enddo

      Nlj          = 0
      Nbl          = 0
      Ndihedral_i  = 0
      Ncoupl       = 0
      Nstrstr      = 0
      Nstrbend     = 0
      Nbendbend    = 0
      Nstrdihed    = 0
      Nbenddihed   = 0
      Ndihedihed   = 0



c     Reads moscito .system file 
 
      iline = 0
 10   READ(iun,'(a)',err=2232) card
      iline = iline + 1
      card0 = card
      if( index(card,'#') .eq. 1) goto 10
      if( Lenw(card).lt.2)        goto 10
      if( index(card,'end{molecule}') .gt. 0) goto 100
      CALL SPLIT (card,word,nsplit,maxspl)

      if (word(1) .eq. '#') goto 10
      IF( word(1) .eq. 'begin{bonds}') THEN

         ns = 0
 11      read(iun,1000) string
         if( index(card,'#') .eq. 1) goto 11
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 11
         if(subword(1) .eq. 'end{bonds}') then
            nstretch = ns
            goto 10
         elseif(subword(1) .eq. 'harmonic') then
            z_harmonic = .true.
            ns = ns + 1
            read(subword(2),*,err=2233) iat_s(ns)
            read(subword(3),*,err=2233) jat_s(ns)
            if (iat_s(ns) .gt. jat_s(ns)) then  !! ctrl i < j
               jnew      = iat_s(ns)
               iat_s(ns) = jat_s(ns)
               jat_s(ns) = jnew
            endif
            read(subword(4),*,err=2233) k_s(ns)
            read(subword(5),*,err=2233) r0_s(ns)
            goto 11
         else
            goto 11
         endif


      ELSEIF( word(1) .eq. 'begin{angle}') THEN

         nb  = 0
         nbl = 0
 12      read(iun,1000) string
         if( index(card,'#') .eq. 1) goto 12
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 12
         if(subword(1) .eq. 'end{angle}') then
            nbend   = nb
            nbend_l = nbl
            goto 10
         elseif(subword(1) .eq. 'linear') then
            z_linear = .true.
            nbl = nbl + 1
            Call Read_1atom (iat_bl(nbl),subword(2),nat,string)
            Call Read_1atom (jat_bl(nbl),subword(3),nat,string)
            Call Read_1atom (kat_bl(nbl),subword(4),nat,string)
            if (iat_bl(nbl) .gt. kat_bl(nbl)) then  !! ctrl i < k
               knew        = iat_bl(nbl)
               iat_bl(nbl) = kat_bl(nbl)
               kat_bl(nbl) = knew
            endif
            read(subword(5),*,err=2233) k_bl(nbl)
            read(subword(6),*,err=2233) t0_bl(nbl)   
            goto 12
         else
            nb = nb + 1
            Call Read_1atom (iat_b(nb),subword(1),nat,string)
            Call Read_1atom (jat_b(nb),subword(2),nat,string)
            Call Read_1atom (kat_b(nb),subword(3),nat,string)
            if (iat_b(nb) .gt. kat_b(nb)) then  !! ctrl i < k
               knew      = iat_b(nb)
               iat_b(nb) = kat_b(nb)
               kat_b(nb) = knew
            endif
            read(subword(4),*,err=2233) k_b(nb)
            read(subword(5),*,err=2233) t0_b(nb)
            goto 12
         endif


      ELSEIF( word(1) .eq. 'begin{dihedral}') THEN

         nd  = 0
         ndi = 0
 13      read(iun,1000) string
         if( index(card,'#') .eq. 1) goto 13
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1)(1:1) .eq. '#') goto 13
         if(subword(1) .eq. 'end{dihedral}') then
            ndihedral   = nd
            ndihedral_i = ndi
            goto 10
         elseif(subword(1) .eq. 'proper') then
            Ndihedral   = 0
            nd = nd + 1
            Call Read_1atom (iat_d(nd),subword(2),nat,string)
            Call Read_1atom (jat_d(nd),subword(3),nat,string)
            Call Read_1atom (kat_d(nd),subword(4),nat,string)
            Call Read_1atom (lat_d(nd),subword(5),nat,string)
            if (jat_d(nd) .gt. kat_d(nd)) then  !! ctrl j < k
               lnew      = iat_d(nd)
               knew      = jat_d(nd)
               iat_d(nd) = lat_d(nd)
               jat_d(nd) = kat_d(nd)
               kat_d(nd) = knew
               lat_d(nd) = lnew
            endif
            if (index(subword(msplit),'#') .eq. 1) 
     &          msplit = msplit -1  
            if (subword(msplit-1) .eq. '#') 
     &          msplit = msplit - 2
            if (subword(msplit-2) .eq. '#') 
     &          msplit = msplit - 3
            ic(nd) = 0
            do iterm=6,msplit,3
              ic(nd) = ic(nd) + 1
              read(subword(iterm),*,err=2233)   k_d(nd,ic(nd))
              read(subword(iterm+1),*,err=2233) n0_d(nd,ic(nd))
              read(subword(iterm+2),*,err=2233) p0_d(nd,ic(nd))
            enddo
            goto 13
         elseif(subword(1) .eq. 'improper') then
            z_improper = .true.
            ndi = ndi + 1
            Call Read_1atom (iat_di(ndi),subword(2),nat,string)
            Call Read_1atom (jat_di(ndi),subword(3),nat,string)
            Call Read_1atom (kat_di(ndi),subword(4),nat,string)
            Call Read_1atom (lat_di(ndi),subword(5),nat,string)
            if (jat_di(ndi) .gt. kat_di(ndi)) then  !! ctrl j < k
               lnew        = iat_di(ndi)
               knew        = jat_di(ndi)
               iat_di(ndi) = lat_di(ndi)
               jat_di(ndi) = kat_di(ndi)
               kat_di(ndi) = knew
               lat_di(ndi) = lnew
            endif
            read(subword(6),*,err=2233) k_di(ndi)
            read(subword(7),*,err=2233) p0_di(ndi)
            goto 13
         endif

      ELSEIF( word(1) .eq. 'begin{couplings}') THEN

         z_c  = .true.
         nss  = 0
         nbb  = 0
         nsb  = 0
         nsd  = 0
         nbd  = 0
         ndd  = 0
 14      read(iun,1000) string
         if( index(card,'#') .eq. 1) goto 14
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 14
         if(subword(1) .eq. 'end{couplings}') then
            Nstrstr     = nss
            Nbendbend   = nbb
            Nstrbend    = nsb
            Nstrdihed   = nsd
            Nbenddihed  = nbd
            Ndihedihed  = ndd
            goto 10
         elseif(subword(1) .eq. 'bond_bond') then
            z_coupl(1) = .true.
            nss = nss + 1
            read(subword(2),*,err=2233)  iat_ss(nss)
            read(subword(3),*,err=2233)  jat_ss(nss)
            if (iat_ss(nss) .gt. jat_ss(nss)) then  !! ctrl i < j
               jnew        = iat_ss(nss)
               iat_ss(nss) = jat_ss(nss)
               jat_ss(nss) = jnew
            endif
            read(subword(4),*,err=2233)  r1_ss(nss)
            read(subword(5),*,err=2233)  kat_ss(nss)
            read(subword(6),*,err=2233)  lat_ss(nss)
            if (kat_ss(nss) .gt. lat_ss(nss)) then  !! ctrl k < l
               lnew        = kat_ss(nss)
               kat_ss(nss) = lat_ss(nss)
               lat_ss(nss) = lnew
            endif
            read(subword(7),*,err=2233)  r2_ss(nss)
            read(subword(8),*,err=2233)  k_ss(nss)
            goto 14
         elseif(subword(1) .eq. 'angle_angle') then
            z_coupl(2) = .true.
            nbb = nbb + 1
            read(subword(2),*,err=2233)  iat_bb(nbb)
            read(subword(3),*,err=2233)  jat_bb(nbb)
            read(subword(4),*,err=2233)  kat_bb(nbb)
            if (iat_bb(nbb) .gt. kat_bb(nbb)) then  !! ctrl i < k
               knew        = iat_bb(nbb)
               iat_bb(nbb) = kat_bb(nbb)
               kat_bb(nbb) = knew
            endif
            read(subword(5),*,err=2233)  t1_bb(nbb)
            read(subword(6),*,err=2233)  lat_bb(nbb)
            read(subword(7),*,err=2233)  mat_bb(nbb)
            read(subword(8),*,err=2233)  nat_bb(nbb)
            if (lat_bb(nbb) .gt. nat_bb(nbb)) then  !! ctrl l < n
               nnew        = lat_bb(nbb)
               lat_bb(nbb) = nat_bb(nbb)
               nat_bb(nbb) = nnew
            endif
            read(subword(9),*,err=2233)  t2_bb(nbb)
            read(subword(10),*,err=2233) k_bb(nbb)
            goto 14
         elseif(subword(1) .eq. 'bond_angle') then
            z_coupl(3) = .true.
            nsb = nsb + 1
            read(subword(2),*,err=2233)  iat_sb(nsb)
            read(subword(3),*,err=2233)  jat_sb(nsb)
            if (iat_sb(nsb) .gt. jat_sb(nsb)) then  !! ctrl i < j
               jnew        = iat_sb(nsb)
               iat_sb(nsb) = jat_sb(nsb)
               jat_sb(nsb) = jnew
            endif
            read(subword(4),*,err=2233)  r0_sb(nsb)
            read(subword(5),*,err=2233)  kat_sb(nsb)
            read(subword(6),*,err=2233)  lat_sb(nsb)
            read(subword(7),*,err=2233)  mat_sb(nsb)
            if (kat_sb(nsb) .gt. mat_sb(nsb)) then  !! ctrl k < m
               mnew        = kat_sb(nsb)
               kat_sb(nsb) = mat_sb(nsb)
               mat_sb(nsb) = mnew
            endif
            read(subword(8),*,err=2233)  t0_sb(nsb)
            read(subword(9),*,err=2233)  k_sb(nsb)
            goto 14
         elseif(subword(1) .eq. 'bond_dihed') then
            z_coupl(4) = .true.
            nsd = nsd + 1
            read(subword(2),*,err=2233)  iat_sd(nsd)
            read(subword(3),*,err=2233)  jat_sd(nsd)
            if (iat_sd(nsd) .gt. jat_sd(nsd)) then  !! ctrl i < j
               jnew        = iat_sd(nsd)
               iat_sd(nsd) = jat_sd(nsd)
               jat_sd(nsd) = jnew
            endif
            read(subword(4),*,err=2233)  r0_sd(nsd)
            read(subword(5),*,err=2233)  kat_sd(nsd)
            read(subword(6),*,err=2233)  lat_sd(nsd)
            read(subword(7),*,err=2233)  mat_sd(nsd)
            read(subword(8),*,err=2233)  nat_sd(nsd)
            if (lat_sd(nsd) .gt. mat_sd(nsd)) then  !! ctrl l < m
               nnew        = kat_sd(nsd)
               mnew        = lat_sd(nsd)
               kat_sd(nsd) = nat_sd(nsd)
               lat_sd(nsd) = mat_sd(nsd)
               mat_sd(nsd) = mnew
               nat_sd(nsd) = nnew
            endif
            read(subword(9),*,err=2233)  n0_sd(nsd)
            if(n0_sd(nsd) .eq. 0) then
              write(6,'(//8x,"ERROR in definition of couplings",//)')
              write(6,'(//8x,"CHECK moscito.system file       ",//)')
              write(6,'(//8x," in bond - dihedral coupling ",//)')
              write(6,*) 'coupled torsional term cannot have n = 0'
              stop 'input error'
            endif
            read(subword(10),*,err=2233) g_sd(nsd)
            read(subword(11),*,err=2233) k_sd(nsd)
            goto 14
         elseif(subword(1) .eq. 'angle_dihed') then
            z_coupl(5) = .true.
            nbd = nbd + 1
            read(subword(2),*,err=2233)  iat_bd(nbd)
            read(subword(3),*,err=2233)  jat_bd(nbd)
            read(subword(4),*,err=2233)  kat_bd(nbd)
            if (iat_bd(nbd) .gt. kat_bd(nbd)) then  !! ctrl i < k
               knew        = iat_bd(nbd)
               iat_bd(nbd) = kat_bd(nbd)
               kat_bd(nbd) = knew
            endif
            read(subword(5),*,err=2233)  t0_bd(nbd)
            read(subword(6),*,err=2233)  lat_bd(nbd)
            read(subword(7),*,err=2233)  mat_bd(nbd)
            read(subword(8),*,err=2233)  nat_bd(nbd)
            read(subword(9),*,err=2233)  pat_bd(nbd)
            if (mat_bd(nbd) .gt. nat_bd(nbd)) then  !! ctrl m < n
               pnew        = lat_bd(nbd)
               nnew        = mat_bd(nbd)
               lat_bd(nbd) = pat_bd(nbd)
               mat_bd(nbd) = nat_bd(nbd)
               nat_bd(nbd) = nnew
               pat_bd(nbd) = pnew
            endif
            read(subword(10),*,err=2233) n0_bd(nbd)
            if(n0_bd(nbd) .eq. 0) then
              write(6,'(//8x,"ERROR in definition of couplings",//)')
              write(6,'(//8x,"CHECK moscito.system file       ",//)')
              write(6,'(//8x," in angle - dihedral coupling ",//)')
              write(6,*) 'coupled torsional term cannot have n = 0'
              stop 'input error'
            endif
            read(subword(11),*,err=2233) g_bd(nbd)
            read(subword(12),*,err=2233) k_bd(nbd)
            goto 14
         elseif(subword(1) .eq. 'dihed_dihed') then
            z_coupl(6) = .true.
            ndd = ndd + 1
            read(subword(2),*,err=2233)  iat_dd(ndd)
            read(subword(3),*,err=2233)  jat_dd(ndd)
            read(subword(4),*,err=2233)  kat_dd(ndd)
            read(subword(5),*,err=2233)  lat_dd(ndd)
            read(subword(6),*,err=2233)  n1_dd(ndd)
            if(n1_dd(ndd) .eq. 0) then
              write(6,'(//8x,"ERROR in definition of couplings",//)')
              write(6,'(//8x,"CHECK moscito.system file       ",//)')
              write(6,'(//8x," in dihedral - dihedral coupling ",//)')
              write(6,*) 'coupled torsional term cannot have n = 0'
              stop 'input error'
            endif
            read(subword(7),*,err=2233)  g1_dd(ndd)
            if (jat_dd(ndd) .gt. kat_dd(ndd)) then  !! ctrl j < k
               lnew        = iat_dd(ndd)
               knew        = jat_dd(ndd)
               iat_dd(ndd) = lat_dd(ndd)
               jat_dd(ndd) = kat_dd(ndd)
               kat_dd(ndd) = knew
               lat_dd(ndd) = lnew
            endif
            read(subword(8),*,err=2233)  mat_dd(ndd)
            read(subword(9),*,err=2233)  nat_dd(ndd)
            read(subword(10),*,err=2233) pat_dd(ndd)
            read(subword(11),*,err=2233) rat_dd(ndd)
            read(subword(12),*,err=2233) n2_dd(ndd)
            if(n2_dd(ndd) .eq. 0) then
              write(6,'(//8x,"ERROR in definition of couplings",//)')
              write(6,'(//8x,"CHECK moscito.system file       ",//)')
              write(6,'(//8x," in dihedral - dihedral coupling ",//)')
              write(6,*) 'coupled torsional term cannot have n = 0'
              stop 'input error'
            endif
            read(subword(13),*,err=2233) g2_dd(ndd)
            if (mat_dd(ndd) .gt. pat_dd(ndd)) then  !! ctrl m < p
               rnew        = mat_dd(ndd)
               pnew        = nat_dd(ndd)
               mat_dd(ndd) = rat_dd(ndd)
               nat_dd(ndd) = pat_dd(ndd)
               pat_dd(ndd) = pnew
               rat_dd(ndd) = rnew
            endif
            read(subword(14),*,err=2233) k_dd(ndd)
            goto 14
         else 
            write(6,'(//8x,"ERROR in definition of couplings",//)')
            write(6,'(//8x,"CHECK moscito.system file       ",//)')
            stop 'SYNTAX ERROR'
         endif

      ELSEIF( word(1) .eq. 'begin{nonbonded}') THEN

         z_lj = .true.
         lj  = 0
 15      read(iun,1000) string
         if( index(card,'#') .eq. 1) goto 15
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 15
         if(subword(1) .eq. 'end{nonbonded}') then
            nlj = lj
            goto 10
         else
            lj = lj + 1
            read(subword(1),*,err=2233) iat_lj(lj)
            read(subword(2),*,err=2233) jat_lj(lj)
            read(subword(3),*,err=2233) sig_lj(lj)
            read(subword(4),*,err=2233) eps_lj(lj)
            if (msplit .ge. 5) then
                if(subword(5) .ne. '#') then    
                  z_intcharge = .true.
                  read(subword(5),*,err=2233) fqq(lj)
                endif
            endif
            goto 15
         endif

      ELSEIF( word(1) .eq. 'begin{configuration}') THEN   

         iat  = 0
 16      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 16
         if( Lenw(string).lt. 3)       goto 16
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 16
         if(subword(1) .eq. 'end{configuration}') then
             Nato = iat
             goto 10
         else
             iat = iat + 1
             read(subword(1),*,err=2233) name(iat)
             read(subword(2),*,err=2233) x(iat)
             read(subword(3),*,err=2233) y(iat)
             read(subword(4),*,err=2233) z(iat)
             goto 16
         endif            

      ELSEIF( word(1) .eq. 'begin{sites}') THEN   

         isit  = 0
 17      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 17
         if( Lenw(string).lt. 3)       goto 17
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 17
         if(subword(1) .eq. 'end{sites}') then
             Nsit = isit
             goto 10
         else
             isit = isit + 1
             read(subword(1),*,err=2233) namesit(isit)
             read(subword(2),*,err=2233) charge(isit)
             read(subword(3),*,err=2233) amass(isit)
             goto 17
         endif            

      ELSEIF( word(1) .eq. 'begin{lorentz_lj}') THEN   

         isit  = 0
 18      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 18
         if( Lenw(string).lt. 3)       goto 18
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 18
         if(subword(1) .eq. 'end{lorentz_lj}') then
             cr   = 2
             Nsit = isit
             goto 10
         else
             isit = isit + 1
             read(subword(1),*,err=2233) namesit2(isit)
             read(subword(2),*,err=2233) sig_sit(isit)
             read(subword(3),*,err=2233) eps_sit(isit)
             goto 18
         endif            

      ELSEIF( word(1) .eq. 'begin{opls_lj}') THEN   

         isit  = 0
 19      read(iun,1000) string
         if( index(string,'#') .eq. 1) goto 19
         if( Lenw(string).lt. 3)       goto 19
         CALL SPLIT (string,subword,msplit,maxspl)
         if(subword(1) .eq. '#') goto 19
         if(subword(1) .eq. 'end{opls_lj}') then
             cr   = 3
             Nsit = isit
             goto 10
         else
             isit = isit + 1
             read(subword(1),*,err=2233) namesit2(isit)
             read(subword(2),*,err=2233) sig_sit(isit)
             read(subword(3),*,err=2233) eps_sit(isit)
             goto 19
         endif            


      ENDIF
      goto 10


 100  continue

c     Assings charge values to charge-charge

      if (z_intcharge) then
              do lj = 1,Nlj     
                 labchk1 = name(iat_lj(lj)) 
                 labchk2 = name(jat_lj(lj)) 
                 call compact(labchk1,lw1)
                 call compact(labchk2,lw2)
                 do isit =1,Nsit
                    call compact(namesit(isit),lw3)
                    if(namesit(isit)(1:lw3) .eq. labchk1(1:lw1))
     &              qi = charge(isit) 
                    if(namesit(isit)(1:lw3) .eq. labchk2(1:lw2))
     &              qj = charge(isit) 
                 enddo
                 qq(lj) = fqq(lj)*qi*qj
              enddo
      endif
      

c    sets defaults for intramolecular nonbonded

      do iat = 1,Nato
         labchk1 = name(iat)
         call compact(labchk1,lw)
         do isit =1,Nsit
            call compact(namesit(isit),lw3)
            if(namesit(isit)(1:lw3) .eq. labchk1(1:lw)) then
               stemp(iat) = sig_sit(isit)
               etemp(iat) = eps_sit(isit)
            endif
         enddo
      enddo
      do iat = 1,Nato
         sigma(iat)   = stemp(iat)
         epsilon(iat) = etemp(iat)
      enddo

c     groups cosine terms in proper torsions 
      
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


c     Prints out potential terms

      write(6,'(/a)') '0) Molecular Geometry'
      write(6,*)
      write(6,1001) 'Site','Name','x','y','z'
      do iat = 1,Nato
         write(6,1002) iat,name(iat),x(iat),y(iat),z(iat)
      enddo
      write(6,1000) 

      write(6,'(/a)') '1) Intermolecular parameters'
      write(6,*)
      write(6,1031) 'Site','Name','charge','sigma','epsilon'
      do iat = 1,Nato
         write(6,1032) iat,name(iat),charge(iat),sigma(iat),epsilon(iat)
      enddo
      write(6,1000) 

      if (z_harmonic) then
        write(6,1000) '2) Stretching parameters'
        write(6,1000) 
        write(6,1003) 'Bond','Atoms','k_s','r0'
        do is = 1,Nstretch
           write(6,1004) 
     &           is,name(iat_s(is)),name(jat_s(is)),k_s(is),r0_s(is)
     &          ,iat_s(is),jat_s(is)
        enddo
        write(6,1000) 
      endif

      write(6,1000) '3) Bending parameters'
      write(6,1000) 
      if (z_linear) then
        write(6,1000) '3.1 Harmonic bending'
        write(6,1000) 
      endif
      write(6,1005) 'Angle','Atoms','k_b','theta0'
      do ib = 1,Nbend
         write(6,1006)ib,name(iat_b(ib)),name(jat_b(ib)),
     &                 name(kat_b(ib)),k_b(ib),t0_b(ib)
     &                ,iat_b(ib),jat_b(ib),kat_b(ib)
      enddo
      if (z_linear) then
         write(6,1000) 
         write(6,1000) '3.2 Linear parameters'
         write(6,1000)
         write(6,1007) 'Angle','Atoms','k_b linear'
         do ib = 1,Nbend_l
         write(6,1006)ib,name(iat_bl(ib)),name(jat_bl(ib)),
     &                 name(kat_bl(ib)),k_bl(ib)
     &                ,iat_bl(ib),jat_bl(ib),kat_bl(ib)
         enddo
         write(6,1000) 
      endif

      write(6,1000) 
      write(6,1000) '4) Torsional parameters'
      write(6,1000) 
      if (Ndihedral .gt. 0) then 
        write(6,1000) '4.1 Fourier torsions'
        write(6,1000) 
        write(6,1008)
      endif
      do id = 1,Ndihedral
         write(6,1009)id,name(iat_d(id)),name(jat_d(id)),
     &        name(kat_d(id)),name(lat_d(id)),ic(id)
     &                   ,1,k_d(id,1),n0_d(id,1),p0_d(id,1)
     &                   ,iat_d(id),jat_d(id),kat_d(id),lat_d(id)
         do it = 2,ic(id)
            write(6,1011) it,k_d(id,it),n0_d(id,it),p0_d(id,it)
         enddo
      enddo
      if (z_improper) then 
         write(6,1000) 
         write(6,1000) '4.2 Harmonic torsions'
         write(6,1000) 
         write(6,1012) 'Dihedral','Atoms','k_d','gamma0'
         do id = 1,Ndihedral_i
         write(6,1013)id,name(iat_di(id)),name(jat_di(id)),
     &        name(kat_di(id)),name(lat_di(id)),k_di(id),
     &        p0_di(id),iat_di(id),jat_di(id),kat_di(id),lat_di(id)
         enddo
      endif
      write(6,1000)
 
      if (z_c) then
        write(6,1000) '5) Coupled parameters'
        write(6,1000) 
      endif
      if (z_coupl(1)) then
        write(6,1000) '5.1  Stretching-stretching coupling'
        write(6,1000) 
        write(6,1016) 'Term','Atoms','k_ss','r01','r02'
        do iss = 1,Nstrstr
            write(6,1017)iss,name(iat_ss(iss)),name(jat_ss(iss)),
     &                   '<->',name(kat_ss(iss)),name(lat_ss(iss)),
     &                   k_ss(iss),r1_ss(iss),r2_ss(iss)
        enddo
      endif
      if (z_coupl(2)) then
        write(6,1000) 
        write(6,1000) '5.2  Bending-bending coupling'
        write(6,1000) 
        write(6,1018) 'Term','Atoms','k_bb','theta01','theta02'
        do ibb = 1,Nbendbend
            write(6,1019)ibb,name(iat_bb(ibb)),name(jat_bb(ibb)),
     &                   name(kat_bb(ibb)),'<->',name(lat_bb(ibb)),
     &                   name(mat_bb(ibb)),name(nat_bb(ibb)),
     &                   k_bb(ibb),t1_bb(ibb),t2_bb(ibb)
        enddo
      endif
      if (z_coupl(3)) then
        write(6,1000) 
        write(6,1000) '5.3  Stretching-bending coupling'
        write(6,1000) 
        write(6,1021) 'Term','Atoms','k_sb','r0','theta0'
        do isb = 1,Nstrbend
            write(6,1022)isb,name(iat_sb(isb)),name(jat_sb(isb)),
     &                   '<->',name(kat_sb(isb)),name(lat_sb(isb)),
     &                   name(mat_sb(isb)),
     &                   k_sb(isb),r0_sb(isb),t0_sb(isb)
        enddo
      endif
      if (z_coupl(4)) then
        write(6,1000) 
        write(6,1000) '5.4 Stretching-torsion coupling'
        write(6,1000)
        write(6,1023) 'Term','Atoms','k_sd','r0','n','gamma'
        do isd = 1,Nstrdihed
            write(6,1024)isd,name(iat_sd(isd)),name(jat_sd(isd)),
     &                   '<->',name(kat_sd(isd)),name(lat_sd(isd)),
     &                   name(mat_sd(isd)),name(nat_sd(isd)),
     &                   k_sd(isd),r0_sd(isd),n0_sd(isd),g_sd(isd)
         enddo
         write(6,1000) 
      endif
      if (z_coupl(5)) then
        write(6,1000) 
        write(6,1000) '5.5 Bending-torsion coupling'
        write(6,1000)
        write(6,1025) 'Term','Atoms','k_bd','theta0','n','gamma'
        do ibd = 1,Nbenddihed
            write(6,1026)ibd,name(iat_bd(ibd)),name(jat_bd(ibd)),
     &                   name(kat_bd(ibd)),'<->',name(lat_bd(ibd)),
     &                   name(mat_bd(ibd)),name(nat_bd(ibd)),
     &                   name(pat_bd(ibd)),
     &                   k_bd(ibd),t0_bd(ibd),n0_bd(ibd),g_bd(ibd)
         enddo
         write(6,1000) 
      endif
      if (z_coupl(6)) then
        write(6,1000) 
        write(6,1000) '5.6 Torsion-torsion coupling'
        write(6,1000)
        write(6,1027) 'Term','Atoms','k_dd','n1','n2','gamma1','gamma2'
        do idd = 1,Ndihedihed
            write(6,1028)idd,name(iat_dd(idd)),name(jat_dd(idd)),
     &                   name(kat_dd(idd)),name(lat_dd(idd)),
     &                   '<->',name(mat_dd(idd)),name(nat_dd(idd)),
     &                   name(pat_dd(idd)),name(rat_dd(idd)),
     &                   k_dd(idd),
     &                   n1_dd(idd),n2_dd(idd),g1_dd(idd),g2_dd(idd)
         enddo
         write(6,1000) 
      endif


      if (z_lj) then
         write(6,1000) '6) Nonbonded parameters'
         write(6,1000) 
         if (z_intcharge) then
           write(6,1000)  'Atomic charges'
           do isit =1,Nsit
             write(6,1030) namesit(isit), charge(isit)
           enddo
           write(6,1000) 
           write(6,1014) 'Couple','Atoms','sigma (A)',
     &                   'epsilon (K)','f','f*q*q (e^2)'
           do  lj=1,Nlj
              write(6,1015)  
     &        lj,name(iat_lj(lj)),name(jat_lj(lj)),
     &        sig_lj(lj),eps_lj(lj),fqq(lj),qq(lj)
           enddo
         else        
           write(6,1014) 'Couple','Atoms','sigma (A)','epsilon (K)'
           do  lj=1,Nlj
              write(6,1015)  
     &        lj,name(iat_lj(lj)),name(jat_lj(lj)),sig_lj(lj),eps_lj(lj)
           enddo
         endif
      endif 


 1000  format(a)
 1001  format(a4,2x,a4,3x,3(7x,a1,7x))
 1002  format(i2,5x,a3,3x,3(f12.6,3x))
 1003  format(a4,6x,a5,12x,a3,9x,a2)
 1004  format(1x,i2,5x,2(a3,3x),4x,f9.4,3x,f8.3,i8,i4)
 1005  format(a5,8x,a5,7x,a3,9x,a6)
 1006  format(2x,i2,4x,3(a3,2x),f9.4,2x,f8.3,i8,2i4)
 1007  format(a5,8x,a5,7x,a12)
 1008  format(' Dihedral',6x,'Atoms',8x,'Nterms iterm   K_d'
     &              ,'     n   gamma')
 1009  format(3x,i2,6x,4(a3,1x),i5,i6,f10.4,i4,f9.2,i8,3i4)
 1010  format(10x,a5,3x,a3,6x,a1,4x,a5)
 1011  format(33x,i5,f10.4,i4,f9.2)
 1012  format(a8,6x,a5,15x,a3,5x,a6)
 1013  format(3x,i2,6x,4(a3,1x),7x,f9.4,3x,f8.2,i8,3i4)
 1014  format(a6,4x,a5,11x,a9,3x,a12,4x,a2,2x,a15)
 1015  format(2x,i3,5x,2(a3,3x),4x,f8.3,5x,f9.4,3x,f9.4,1x,f9.4)
 1016  format(a5,9x,a5,22x,a4,6x,a3,6x,a3)
 1017  format(1x,i2,6x,2(a3,1x),a3,1x,2(a3,1x),8x,f9.4,3x,f6.2,3x,f6.2)
 1018  format(a5,9x,a5,22x,a4,5x,a7,2x,a7)
 1019  format(1x,i2,6x,3(a3,1x),a3,1x,3(a3,1x),f9.4,3x,f7.1,3x,f7.1)
 1021  format(a5,9x,a5,22x,a4,6x,a3,6x,a6)
 1022  format(1x,i2,6x,2(a3,1x),a3,1x,3(a3,1x),4x,f9.4,4x,f5.2,3x,f7.1)
 1023  format(a5,9x,a5,22x,a4,8x,a2,4x,a1,3x,a5)
 1024  format(1x,i2,6x,3(a3,1x),a3,1x,3(a3,1x),f9.4,4x,f5.2,3x,i2,1x,
     &        f7.1)
 1025  format(a5,9x,a5,22x,a4,4x,a6,4x,a1,3x,a5)
 1026  format(1x,i2,2x,3(a3,1x),a3,1x,4(a3,1x),f9.4,3x,f6.2,3x,i2,2x,
     &        f7.1)
 1027  format(a5,9x,a5,22x,a4,4x,2(a3,2x),2(a6,2x))
 1028  format(1x,i2,1x,4(a3),a3,1x,4(a3,1x),f9.4,3x,2(i3,2x),
     &        2(f6.1,2x))
 1030  format(2x,a4,3x,f8.4)           
 1031  format(a4,2x,a4,3x,3(3x,a7,3x))
 1032  format(i2,5x,a3,3x,5(f12.6,3x))


c     before going on, check geometry !
 
      IF (z_test) THEN
          write(6,'(/" CHECK GEOMETRY: begins")')
          if(nat.ne.Nato) then
            write(6,'(/4x,"Error: wrong n. of atoms",i6,/)') Nato
            stop
          endif
          ermax=0.0d0
          do i=1,nat
            dx=geom(i,1)-x(i)
            dy=geom(i,2)-y(i)
            dz=geom(i,3)-z(i)
            dist=dsqrt(dx*dx+dy*dy+dz*dz)
            if(dist.gt.0.100d0) then
              write(6,68) i,dist
              stop
            endif
          enddo
          write(6,'( " CHECK GEOMETRY:    OK!")')
      ELSE
          write(6,'(//" WARNING: TEST ON GEOMETRY IS SKIPPED !!!",/)')
          write(6,'( /" WARNING: TEST ON GEOMETRY IS SKIPPED !!!",/)')
          write(6,'( /" WARNING: TEST ON GEOMETRY IS SKIPPED !!!",/)')
      ENDIF
   68 format(/4x,'Error: atom',i4,' in g03 and in moscito are'
     &                 ,f12.7,' distant',//)


c      Unit Conversions to Joyce 
c        1 kJ/mol       = 1/2625.5 Hartree  
c        1 K            = 8.31451121d-3/2625.5 Hartree 
c        1 A            = 1/0.52917726d0 Bohr    
c        1 kJ/(mol*A^2) = 1.066*10^(-4) H/(Bohr^2) 
c        e^2/4*pi*eps0  = 1 H bohr --> no conversion needed      

      do is=1,Nstretch
         k_s(is) = k_s(is)*1.066/10000.d0
      enddo

      do ib=1,Nbend
         k_b(ib) = k_b(ib)/2625.5d0
      enddo
      
      do ibl=1,Nbend_l
         k_bl(ibl) = k_bl(ibl)/2625.5d0
      enddo

      do id = 1,Ndihedral
         do it = 1,ic(id)
            k_d(id,it) = k_d(id,it)/2625.5d0
         enddo
      enddo

      do it = 1,Ndihedral_i
         k_di(it) = k_di(it)/2625.5d0
      enddo

      do iss=1,Nstrstr
         k_ss(iss) = k_ss(iss)*1.066/10000.d0
      enddo

      do ibb=1,Nbendbend
         k_sb(ibb) = k_bb(ibb)/2625.5d0
      enddo

      do isb=1,Nstrbend
         k_sb(isb) = k_sb(isb)*2.01553/10000.d0
      enddo

      do isd=1,Nstrdihed
         k_sd(isd) = k_sd(isd)*2.01553/10000.d0
      enddo

      do ibd=1,Nbenddihed
         k_bd(ibd) = k_bd(ibd)/2625.5d0
      enddo

      do idd=1,Ndihedihed
         k_dd(idd) = k_dd(idd)/2625.5d0
      enddo

      do lj = 1,Nlj
         eps_lj(lj)= eps_lj(lj)*8.31451121d-3/2625.5d0
         sig_lj(lj)= sig_lj(lj)/0.52917726d0
      enddo


c---------------------------------------------------------------c
c            PUTS DATA ON iflag ARRAY                    
c                       AND                              
c            PARAMETERS IN PAR VECTOR                    
c---------------------------------------------------------------c
      
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
        iflag(int,8)= 2
        parfun(int) =r0_s(is)/boh2ang
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
        iflag(int,1)=iat_b(ib)
        iflag(int,2)=jat_b(ib)
        iflag(int,3)=kat_b(ib)
        iflag(int,8)= 3
        parfun(int) =t0_b(ib)/todeg 
        card=' '
        write(card,'(a,"-",a,"-",a)') aname(iat_b(ib))
     &               ,aname(jat_b(ib)),aname(kat_b(ib))
        Call Compact (card,lw)
        labpar(int)=card
      enddo

c     iflag(int,0) = 2  ----> bending (linear)
      do ibl=1,Nbend_l
        int=int+1
        iflag(int,0)=2
        par(int)    = k_bl(ibl)
        iflag(int,1)=iat_bl(ibl)
        iflag(int,2)=jat_bl(ibl)
        iflag(int,3)=kat_bl(ibl)
        iflag(int,8)= 3
        parfun(int) =t0_bl(ibl)/todeg 
        card=' '
        write(card,'(a,"-",a,"-",a)') aname(iat_bl(ibl))
     &             ,aname(jat_bl(ibl)),aname(kat_bl(ibl))
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
        parfun(int) =p0_di(it)/todeg 
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
            par(int)    =k_d(id,it)
            iflag(int,1)=iat_d(id)
            iflag(int,2)=jat_d(id)
            iflag(int,3)=kat_d(id)
            iflag(int,4)=lat_d(id)
            iflag(int,5)=n0_d(id,it)
            iflag(int,6)=ic(id)
            iflag(int,8)= 4
            parfun(int)    =p0_d(id,it)/todeg
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
         iflag(int,0)=5
         par(int)    =eps_lj(lj)
         iflag(int,1)=iat_lj(lj)
         iflag(int,2)=jat_lj(lj)
         iflag(int,8)=2
         parfun(int)    =sig_lj(lj)
         card=' '
         write(card,'(a,"-",a)') aname(iat_lj(lj)),aname(jat_lj(lj))
         Call Compact (card,lw)
         labpar(int)=card
      enddo


c     iflag(int,0) = 6  ----> charge-charge
      if (z_intcharge) then
         ichrgprm(1)= int + 1
         do lj = 1,Nlj
            int=int+1
            iflag(int,0)=6
            par(int)    = qq(lj)
            iflag(int,1)=iat_lj(lj)
            iflag(int,2)=jat_lj(lj)
            iflag(int,8)=2
            parfun(int) = 0.d0
            card=' '
            write(card,'(a,"-",a)') aname(iat_lj(lj)),aname(jat_lj(lj))
            Call Compact (card,lw)
            labpar(int)=card
         enddo
         ichrgprm(2)= int
      endif

      Npar = int  !! Npar including charges

c     COUPLINGS

c     1) gets atoms ivolved in bond-bond couplings        

      IF (z_coupl(1)) THEN
       do iss=1,Nstrstr
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_ss(iss)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 1) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            if     (ii .eq. iat_ss(iss) .and. jj .eq. jat_ss(iss)) then
                iflag(int,1) = ichk
                ictrl = ictrl + 1
            elseif (ii .eq. kat_ss(iss) .and. jj .eq. lat_ss(iss)) then
                iflag(int,2) = ichk
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3000
         enddo
         write(6,'(//8x,"ERROR in definition of str-str coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j   +  k  l")')
         write(6,'(8x,i2,4x,4(i2,4x))') 
     &        iss,iat_ss(iss),jat_ss(iss),kat_ss(iss),lat_ss(iss)
         stop 'ERROR'
3000     continue
       enddo
      ENDIF

c     2) gets atoms involved in angle-angle couplings        

      IF (z_coupl(2)) THEN
       do ibb=1,Nbendbend
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_bb(ibb)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 2) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            kk = iflag(ichk,3)
            if     (ii .eq. iat_bb(ibb) .and. 
     &              jj .eq. jat_bb(ibb) .and.
     &              kk .eq. kat_bb(ibb)) then
                iflag(int,1) = ichk
                ictrl = ictrl + 1
            elseif (ii .eq. lat_bb(ibb) .and. 
     &              jj .eq. mat_bb(ibb) .and.
     &              kk .eq. nat_bb(ibb)) then
                iflag(int,2) = ichk
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3100
         enddo
         write(6,'(//8x,"ERROR in definition of bnd-bnd coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j   k  +  l  m   n")')
         write(6,'(8x,i2,4x,6(i2,4x))') 
     &        ibb,iat_bb(ibb),jat_bb(ibb),kat_bb(ibb),lat_bb(ibb),
     &        mat_bb(ibb),nat_bb(ibb)
         stop 'ERROR'
3100     continue
       enddo
      ENDIF

c     3) gets atoms involved in bond-angle couplings        

      IF (z_coupl(3)) THEN
       do isb=1,Nstrbend
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_sb(isb)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 1) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            if     (ii .eq. iat_sb(isb) .and. jj .eq. jat_sb(isb)) then
                iflag(int,1) = ichk
                ictrl        = ictrl + 1
            endif
          elseif (iflag(ichk,0) .eq. 2) then
            kk = iflag(ichk,1)
            ll = iflag(ichk,2)
            mm = iflag(ichk,3)
            if (kk .eq. kat_sb(isb) .and. 
     &          ll .eq. lat_sb(isb) .and.
     &          mm .eq. mat_sb(isb)) then
                iflag(int,2) = ichk
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3200
         enddo
         write(6,'(//8x,"ERROR in definition of str-bnd coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j    +  k  l   m")')
         write(6,'(8x,i2,4x,5(i2,4x))') 
     &        isb,iat_sb(isb),jat_sb(isb),kat_sb(isb),lat_sb(isb),
     &        mat_sb(isb)
         stop 'ERROR'
3200     continue
       enddo
      ENDIF

c     4) gets atoms involved in stretching-dihedral couplings        

      IF (z_coupl(4)) THEN
       do isd=1,Nstrdihed
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_sd(isd)
         parfun(int)  = g_sd(isd)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 1) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            if     (ii .eq. iat_sd(isd) .and.
     &              jj .eq. jat_sd(isd)) then
                iflag(int,1) = ichk
                ictrl        = ictrl + 1
            endif
          elseif (iflag(ichk,0) .eq. 4) then
            kk = iflag(ichk,1)
            ll = iflag(ichk,2)
            mm = iflag(ichk,3)
            nn = iflag(ichk,4)
            np = iflag(ichk,5)
            if     (kk .eq. kat_sd(isd) .and. 
     &              ll .eq. lat_sd(isd) .and.
     &              mm .eq. mat_sd(isd) .and.
     &              nn .eq. nat_sd(isd) .and.
     &              np .eq.  n0_sd(isd)) then
                iflag(int,2) = ichk
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3300
         enddo
         write(6,'(//8x,"ERROR in definition of str-tors coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j   +  k  l   m  n ")')
         write(6,'(8x,i2,4x,6(i2,4x))') 
     &        isd,iat_sd(isd),jat_sd(isd),kat_sd(isd),lat_sd(isd),
     &        mat_sd(isd),nat_sd(isd)
         stop 'ERROR'
3300     continue
       enddo
      ENDIF

c     4) gets atoms involved in bend-dihedral couplings        

      IF (z_coupl(5)) THEN
       do ibd=1,Nbenddihed
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_bd(ibd)
         parfun(int)  = g_bd(ibd)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 2) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            kk = iflag(ichk,3)
            if     (ii .eq. iat_bd(ibd) .and.
     &              jj .eq. jat_bd(ibd) .and.
     &              kk .eq. kat_bd(ibd)) then
                iflag(int,1) = ichk
                ictrl        = ictrl + 1
            endif
          elseif (iflag(ichk,0) .eq. 4) then
            ll = iflag(ichk,1)
            mm = iflag(ichk,2)
            nn = iflag(ichk,3)
            pp = iflag(ichk,4)
            np = iflag(ichk,5)
            if     (ll .eq. lat_bd(ibd) .and. 
     &              mm .eq. mat_bd(ibd) .and.
     &              nn .eq. nat_bd(ibd) .and.
     &              pp .eq. pat_bd(ibd) .and.
     &              np .eq.  n0_bd(ibd)) then
                iflag(int,2) = ichk
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3400
         enddo
         write(6,'(//8x,"ERROR in definition of bnd-tors coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j  k  +  l  m   n  p ")')
         write(6,'(8x,i2,4x,6(i2,4x))') 
     &        ibd,iat_bd(ibd),jat_bd(ibd),kat_bd(ibd),lat_bd(ibd),
     &        mat_bd(ibd),nat_bd(ibd),pat_bd(ibd)
         stop 'ERROR'
3400     continue
       enddo
      ENDIF

c     5) gets atoms involved in diheral-dihedral couplings        

      IF (z_coupl(6)) THEN
       do idd=1,Ndihedihed
         int          = int+1
         iflag(int,0) = 8
         par(int)     = k_dd(idd)
         ictrl        = 0
         do ichk = 1, Npar
          if (iflag(ichk,0) .eq. 4) then
            ii = iflag(ichk,1)
            jj = iflag(ichk,2)
            kk = iflag(ichk,3)
            ll = iflag(ichk,4)
            nn = iflag(ichk,5)
            if     (ii .eq. iat_dd(idd) .and. 
     &              jj .eq. jat_dd(idd) .and.
     &              kk .eq. kat_dd(idd) .and.
     &              ll .eq. lat_dd(idd) .and. 
     &              nn .eq. n1_dd(idd) ) then
                iflag(int,1) = ichk
                parfun(int)  = parfun(ichk)
                ictrl        = ictrl + 1
            elseif (ii .eq. mat_dd(idd) .and. 
     &              jj .eq. nat_dd(idd) .and.
     &              kk .eq. pat_dd(idd) .and.
     &              ll .eq. rat_dd(idd) .and.
     &              nn .eq. n2_dd(idd) ) then
                iflag(int,2) = ichk
                parfun(int)  = parfun(ichk)
                ictrl        = ictrl + 1
            endif
          endif
          if (ictrl .eq. 2) goto 3500
         enddo
         write(6,'(//8x,"ERROR in definition of tors-tors coupl.",//)')
         write(6,'(/8x,"Coupling #     i  j  k  l +  m   n  p  r ")')
         write(6,'(8x,i2,4x,8(i2,4x))') 
     &        idd,iat_dd(idd),jat_dd(idd),kat_dd(idd),lat_dd(idd),
     &        mat_dd(idd),nat_dd(idd),pat_dd(idd),rat_dd(idd)
         stop 'ERROR'
3500     continue
       enddo
      ENDIF

      Npar = int

  555 format(i3,2x,7i3,3f9.3)


 1020  format(a6,4(2x,i2),2x,a1,1x,f9.6,3x,a1,4(a3,1x),a1)


      write(6,1000) 
      write(6,1000) 
      write(6,1000)  
     & '       ********************************************'
      write(6,1000)  
     & '          E N D  M O S C I T O  I N P U T'
      write(6,1000)  
     & '       ********************************************'
      write(6,1000) 
      write(6,1000) 



      return

c     I/O errors
2232  write(6,1000) 'I/O error: check your .system file'      
      write(6,1000) 'Last correct line:'
      write(6,2238) iline,':',card0
2238  format(i4,1x,a1,1x,a90)      
      stop
c
2233  write(6,1000) 'I/O error: check your .system file'      
      write(6,1000) 'Syntax error on line:'
      write(6,1000) string
      stop
c      
      end


c-------------------------------------------------------------c
      Subroutine Read_1atom (ia,word,nat,string)
c-------------------------------------------------------------c
      character*(*) word,string
      read(word,*) ia
      if(ia.lt.1) then
        write(6,22) ia,string(1:60)
        stop
      elseif(ia.gt.nat) then
        write(6,24) ia,nat,string(1:60)
        stop
      endif
      return
   22 format(//'  FF file INPUT ERROR: atom n. <1',i6,/4x,a,//)
   24 format(//'  FF file INPUT ERROR: atom n. >NAT',2i6,/,4x,a,//)
      end



c-------------------------------------------------------------c
      Subroutine Split (A,B,nsplit,maxspl)
c-------------------------------------------------------------c


      character*(*) A,B(maxspl)
 
      lena=Len(A)
      nsplit=0
      i=0
   10 i=i+1
      if (i.gt.lena) return
      !Transform tab into blank space
      if (iachar(A(i:i)).eq.9) A(i:i)=' '
      if (A(i:i).eq.' ') goto 10
      if(nsplit.eq.maxspl) return
      nsplit=nsplit+1
      ibeg=i
   20 i=i+1
      if (i.gt.lena) then
      iend=i-1
      B(nsplit)=' '
      do k=ibeg,iend
        mm=k-ibeg+1
        B(nsplit)(mm:mm)=A(k:k)
      enddo
      return
      endif
      if (A(i:i).ne.' ') goto 20
      iend=i-1
      B(nsplit)=' '
      do k=ibeg,iend
        mm=k-ibeg+1
        B(nsplit)(mm:mm)=A(k:k)
      enddo
      goto 10
      end
