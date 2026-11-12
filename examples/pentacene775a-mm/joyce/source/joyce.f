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
      program joyce
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      integer  cr
      parameter (maxat=300, maxpt=20000, maxgeom=20
     &          ,max3n=maxat*3, maxpar=2000, maxspl = 100)
      parameter (maxint=2000)
      dimension geom(maxat,6),amass(maxat),manua(maxat,9)
     &         ,gang(maxat,6),dist(maxat,maxat),kbond(maxat,maxat)
     &         ,wmass(max3n),T(max3n,max3n),FCC(max3n,max3n)
     &         ,FCC_el(max3n,max3n),FCC_whole(max3n,max3n)
     &         ,F(max3n*max3n),D1(max3n*max3n),D2(max3n,max3n)
     &         ,C(max3n,max3n),E(max3n),geom_el(maxat,6)
     &         ,Utrot(6,max3n),valric(0:maxgeom,maxpar),indir(maxpar)
     &         ,iflag(maxpar,0:8),par(maxpar) 
     &         ,parfun(maxpar),parfun2(maxpar)
     &         ,Un(max3n,max3n),iclass(maxpar),ikeepff(maxpar,2) 
     &         ,Fg03(max3n,max3n),Grad(max3n),Grad_el(max3n)
     &         ,UAM(max3n,max3n),Grad_whole(max3n)
     &         ,iflagn(maxpar,0:8),iflagr(maxpar,0:8),vc(maxpar)
     &         ,valic(maxpar),valic0(maxpar),ricf(maxpar),rics(maxpar)
     &         ,alpha(maxpt,maxpar),beta(maxpt),weig(maxpt)
     &         ,work(maxpar),anda(maxpar),inda(maxpar)
     &         ,geom0(maxat,6),geom00(maxat,6),geoms(maxat,6)
     &         ,scamin(maxpar),scamax(maxpar),scainc(maxpar)
     &         ,iscaic(maxpar),iflagrf(maxpar,0:8)
     &         ,xfroz(99),ifroz(99),Freq(max3n)
     &         ,iscan(maxgeom),DEQM(maxgeom)
     &         ,charge(maxat),sitmass(maxat)
     &         ,sigma(maxat),epsilon(maxat)
     &         ,sig_sit(maxat),eps_sit(maxat)
      dimension A(max3n,max3n),B(max3n,1),X(max3n,1)
      dimension ichrgprm(2)
      character card*300,card0*140,word*90,fileinp*90,FFstartfile*90
     &         ,file0*90,string*90,filesc*90,filegg*90,title*120
     &         ,fileqmnrg*90,gnum*30,gnum0*30,filetype*30
     &         ,labpar(maxint)*24,labric(maxpar)*15,labnat(maxpar)*15
     &         ,aname(maxat)*5,bname(maxat)*5,cname(maxat)*5
     &         ,labpt(maxpt)*11,lab1*10,lab2*10,cardf*90
     &         ,scelta*40,stype(8)*6,atlab(maxat)*5
     &         ,namesit(maxat)*5,file0el*90,subword(maxspl)*30
     &         ,ftype*4, freqgracefile*25, torsgracefile*25
     &         ,torsgracefile0*25
      logical zrew,zord,zg03,zover,zmosc,zoptim,zaddtr,zredund
     &       ,zcheck,zorico,zhess,zpic,zreduok,zstop,znohess,znograd
     &       ,zmass,zfreqff,zexist,zgrah,zprt,z_amber,z_lj,z_fitlj
     &       ,z_g09,z_sepel,z_ljassign,z_scan,z_plotfreq,z_plottors
      logical zstardihed(maxpar)
      logical zharmdihed(maxpar)
      common /forstop/ cardf
      common /printout/ levprt
      common /ljassign/z_LJassign
      common /gracefiles/freqgracefile, torsgracefile
      data stype/'stretc','bend  ','tors.i','tors.p'
     &          ,'LJ-nbd','free  ','free  ','coupl '/
      data amu     /1822.8880d0/  !!   1 amu = 1822.888 electron mass
      data todeg   /57.2957795130823d0/
      data boh2ang /0.529177249d0/    !! 1 bohr = 0.529177 Angs
      data boltzev /8.617343d-05/     !! boltzmann const eV/K
      boltzau=boltzev/27.2114d0       !! boltzmann const Hartee/K
 
 
      ncol=10
      write(6,*)
      write(6,*) 
     & '==============================================================='
      write(6,88)
88    format(/14x,'\__    _/(  ___  )|\     /|(  ____ \(  ____ \'
     &      ,/14x,'   )  (  | (   ) |( \   / )| (    \/| (    \/'
     &      ,/14x,'   |  |  | |   | | \ (_) / | |      | (__    '
     &      ,/14x,'   |  |  | |   | |  \   /  | |      |  __)   '
     &      ,/14x,'   |  |  | |   | |   ) (   | |      | (      '
     &      ,/14x,'|\_)  )  | (___) |   | |   | (____/\| (____/\'
     &      ,/14x,'(____/   (_______)   \_/   (_______/(_______/',//
     &      ,/14x,'a program to optimize internal FF from QM data',/)
      write(6,*)
      write(6,886)
  886 format(/14x,'  version 3.0        Copyright (C) 2017    ' 
     &      ,/14x,'                                           '
     &      ,/14x,'            written by                     ' 
     &      ,/14x,'                                           '
     &      ,/14x,'Ivo Cacelli, Javier Cerezo, Samuele Giannini,'
     &      ,/14x,'Nicola De Mitri and Giacomo Prampolini     '
     &      ,/14x,'                                           '
     &      ,/14x,'                                           '
     &      ,/14x,'                                           '
     &      ,/14x,'This program comes with ABSOLUTELY NO WARRANTY; '
     &      ,/14x,'for details type : <go.joyce -lic>              '
     &      ,/14x,'This is free software, and you are welcome to   '
     &      ,/14x,'redistribute it under certain conditions;       '
     &      ,/14x,'type <go.joyce -lic> for details'
     &      ,/14x,'                                           '
     &      ,/14x,'                                           ',/)
      call print_version
      write(6,*) 
     & '==============================================================='


      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')
     &  '==============================================='
      write(6,'(10x,a)')  'JOYCE  Parameterization starts'
      write(6,'(10x,a)')
     &  '==============================================='
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '
c     checks input files for syntax errors
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  'A)  Scanning input file       '
      write(6,'(10x,a)')  '------------------------------'
      iun= 5
      iline = 0
      rewind iun
5     read(iun,'(a)', end=585) card      
      iline=iline+1
      CALL SPLIT (card,subword,msplit,maxspl)
      l0 = Lenw(subword(1))
      if (subword(1)(1:l0) .eq. '$title' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$print' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$gracefreq' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$gracetors' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$sep_el' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$equil' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$mass' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$forcefield' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$fitLJ' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$amber' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$gaussian' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$generate' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$keepff' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$Ckeepff' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$normal' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$wfre' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$whess' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$boltz' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$zero' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$LJassign' ) then 
          goto 5
      elseif (subword(1)(1:l0) .eq. '$dependence' ) then 
6         read(iun,'(a)') word           
          iline=iline+1
          CALL SPLIT (word,subword,msplit,maxspl)
          lw1 = Lenw(subword(1)) 
          if (subword(1)(1:lw1) .eq. '$end' ) then  
             goto 5
          elseif (subword(1)(1:1) .eq. ';' ) then 
             goto 6 
          elseif (subword(1)(1:1) .eq. '#' ) then 
             goto 6 
          else
            goto 6 
          endif
      elseif (subword(1)(1:l0) .eq. '$assign' ) then 
7         read(iun,'(a)') word           
          iline=iline+1
          CALL SPLIT (word,subword,msplit,maxspl)
          lw1 = Lenw(subword(1)) 
          if (subword(1)(1:lw1) .eq. '$end' ) then  
             goto 5
          elseif (subword(1)(1:1) .eq. ';' ) then 
             goto 7 
          elseif (subword(1)(1:1) .eq. '#' ) then 
             goto 7 
          else
            goto 7 
          endif
      elseif (subword(1)(1:l0) .eq. '$geom' ) then 
8         read(iun,'(a)') word           
          iline=iline+1
          CALL SPLIT (word,subword,msplit,maxspl)
          lw1 = Lenw(subword(1)) 
          if (subword(1)(1:lw1) .eq. '$end' ) then  
             goto 5
          elseif (subword(1)(1:1) .eq. ';' ) then 
             goto 8 
          elseif (subword(1)(1:1) .eq. '#' ) then 
             goto 8 
          else
            goto 8 
          endif
      elseif (subword(1)(1:l0) .eq. '$scan' ) then 
9         read(iun,'(a)') word           
          iline=iline+1
          CALL SPLIT (word,subword,msplit,maxspl)
          lw1 = Lenw(subword(1)) 
          if (subword(1)(1:lw1) .eq. '$end' ) then  
             goto 5
          elseif (subword(1)(1:1) .eq. ';' ) then 
             goto 9 
          elseif (subword(1)(1:1) .eq. '#' ) then 
             goto 9 
          else
            goto 9 
          endif
      elseif (subword(1)(1:1) .eq. ';' ) then 
          goto 5
      elseif (subword(1)(1:1) .eq. '#' ) then 
          goto 5
      elseif (msplit .eq. 0 ) then 
          goto 5
      elseif (subword(1)(1:1) .eq.'$') then 
c       Detect an unknown keyword (starting with $)
        write(6,'(15x,a)')'Syntax error detected on input.'
        write(6,'(15x,a)')' Unknown keyword: '//subword(1)
        write(6,'(15x,a18,1x,i3)')'Check line number:', iline       
        write(6,'(15x,a13)') 'Line content:'
        write(6,'(15x,a)') card
        stop 'SYNTAX ERROR: check input'  
      else
c       This detects comments out of any keyword block
        goto 5
      endif


585   continue

      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')
     &  '              input file ok                    '
      write(6,'(10x,a)')' '
      write(6,'(10x,a)')' '

c     printout options (0:3)
      Call Cerca_W ('$print',iun,card,ier,.true.)
      levprt=0
      if(ier.eq.0) read(card,*) levprt
      if(levprt.lt.0) levprt=-1
      if(levprt.gt.4) levprt=4
      if (levprt .gt. 0) write(6,'(10x,"print level:",i2)') levprt
      if (levprt .gt. 0) write(6,'(10x,a)')' '

c     options for grace plots
      z_plotfreq = .false.
      Call Cerca_W ('$gracefreq',iun,card,ier,.true.)
      if (ier. eq. 0) z_plotfreq = .true.
      if (z_plotfreq) then 
        freqgracefile = 'joyce.freq_check.agr'
        CALL SPLIT (card,subword,msplit,maxspl)
        if (msplit .gt. 0 ) read(card,*) freqgracefile
        lwt=Lenw(freqgracefile)
        write(6,'(10x,"Frequency plot required:",a25)')
     &        freqgracefile(1:lwt)
        write(6,'(10x,a)')' '
      endif
      z_plottors = .false.
      Call Cerca_W ('$gracetors',iun,card,ier,.true.)
      if (ier. eq. 0) z_plottors = .true.
      if (z_plottors) then 
        torsgracefile  = 'joyce.tors_check'
        CALL SPLIT (card,subword,msplit,maxspl)
        if (msplit .gt. 0 ) read(card,*) torsgracefile
        torsgracefile0 = torsgracefile
        lwt=Lenw(torsgracefile)
        write(6,'(10x,"Torsions plot required: ",a25)')
     &        torsgracefile(1:lwt)
        write(6,'(10x,a)')' '
      endif

c     read title
      iun=5
      zrew=.true.
      Call Cerca_W ('$tit',iun,title,ier,zrew)
      lwt=Lenw(title)
      if(ier.eq.0) write(6,'(10x,"Title: | ",a," |")') title(1:lwt)
 
 
c     option for separate electrostatic (e.g. MMpol)
c     Obsolete & deprecated : GP October 2023
c     When electrostatic are separately handled, Joyce is 
c     expecting two fchk files (e.g. QM and MMpol)
      z_sepel = .false.
      Call Cerca_W ('$sep_el',iun,card,ier,.true.)
      if (ier. eq. 0) z_sepel = .true.
      if (z_sepel) then 
         write(6,'(/14x,"---------------------------------")') 
         write(6,'(14x,"            WARNING              ")') 
         write(6,'(14x,"---------------------------------")') 
         write(6,'(/14x,"Electrostatics will be subtracted")') 
         write(6,'(14x,"from QM data, and should be added")') 
         write(6,'(14x,"again to correctly complement ")') 
         write(6,'(14x,"Joyce FF.")') 
      endif
 
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  'B) Reading QM training data   '
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  
 
      Call Cerca_W ('$equil',5,card,ier,.true.)
      if(ier.ne.0) then
        write(6,'(//9x,"ERROR: keyword $equil not found",//)')
        stop
      endif
      Call Selwrd (card,file0,1,lwi,ipos)
      word=file0
      Call LowCAs (word)
      lw=Lenw(word)
      if (word(lw-3:lw) .eq. '.fcc') then
        ftype = 'fcc '
      else
        ftype = 'fchk'
        ii=index(word,'.fchk')
        if(ii.le.0) file0=file0(1:lw)//'.fchk'
      endif
      write(6,14) file0(1:Lenw(file0))
   14 format(/10x,'QM INPUT FILE for geom-0 (QM):  ',a,/)

      do j=1,max3n
        do i=1,max3n
          FCC(i,j) = 0.d0
        enddo
      enddo
      iun=31
      open(iun,file=file0,form='formatted',status='old',err=2335)
      goto 2336
2335  write(6,'(a)') 'I/O error: file was not found'      
      write(6,'(a8,2x,a1,1x,a30)') 'Filename',':',file0
      stop
2336  continue      
c     reads QM optimized geometry, energy, gradients and H
c     input QM data for the equil. geometry 
c     (absolute minimum)
c     it must be a fchk file from g03/g09/g16 
c     (freq,opt=verytight)
c     or fcc state file (from FCclasses3)
      if (ftype .eq. 'fcc ') then
        Call Read_fcc  (nat,geom,aname,Etot0,FCC,Grad,iun,maxat,max3n)
      else
        Call Read_fchk (nat,geom,aname,Etot0,FCC,Grad,iun,maxat,max3n)
      endif
      ncart=nat*3
      Etot0_whole = Etot0
      do j=1,ncart
        Grad_whole(j) = Grad(j)
        do i=1,ncart
          Fg03(i,j)=FCC(i,j)
          FCC_whole(i,j)=FCC(i,j)
        enddo
      enddo
      close (iun)
c     checks if Hessian have been read
      znohess=.true.
      do i=1,nat*3
        do j=1,nat*3
          if(dabs(FCC(i,j)).gt.1.d-08) znohess=.false.
        enddo
      enddo
      if(znohess) then 
            write(6,'(a)')
     &     '$equil error: no Hessian elements found: check fchk file'
      endif
 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c              Obsolete & depreacated to be removed . GP october 2023
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     reads gaussian files containing information about
c     electrostatics (e.g. MMpol) checking they are referring
c     to the same geometry
      if (z_sepel) then                         ! option for separate elctrostatics (obsolete)
        Call Selwrd (card,file0el,2,lwi,ipos)
        if (lwi .eq. 0) then
           write(6,'(a)') 
     &     '$equil error: please specify the second filename'
           stop
        endif
        word=file0el
        Call LowCAs (word)
        lw=Lenw(word)
        ii=index(word,'.fchk')
        if(ii.le.0) file0el=file0el(1:lw)//'.fchk'
        write(6,664) file0el(1:Lenw(file0el))
        iun=41
        open(iun,file=file0el,form='formatted',status='old',err=2337)
        goto 2338
2337    write(6,'(a)') 'I/O error: file was not found'      
        write(6,'(a8,2x,a1,1x,a30)') 'Filename',':',file0el
        stop
2338    continue      
        Call Read_fchk (nat_el,geom_el,bname,Etot0_el,
     &                   FCC_el,Grad_el,iun,maxat,max3n)
c       checks basics
        if (nat_el .ne. nat) then
           write(6,'(a)') 'I/O error: QM and low level files have'
           write(6,'(a)') ' different number of atoms!'
        endif
c       checks geometry between files
        ncart=nat*3
        do i=1,nat
          do k=1,3
            check = dabs(geom(i,k) - geom_el(i,k))
            if (check .gt. 1.d-5) then
             write(6,'(a)') 'I/O error: QM and low level geometry'
             write(6,'(a)') 'differ by more than treshold (1.d-5)'
             write(6,'(a15,1x,f12.6)') 'current error:',check
             stop
            endif
          enddo
        enddo
        close (iun)
c       subtract the electrostatic term from the whole QM data
        Etot0 = Etot0_whole - Etot0_el
        do j=1,ncart
           Grad(j) = Grad_whole(j) - Grad_el(j)
           do i=1,ncart
             Fg03(i,j) = FCC_whole(i,j) - FCC_el(i,j)
           enddo
        enddo
c       printouts
        write(6,2012) file0(1:lenw(file0)),file0el(1:lenw(file0el))
     &           ,Etot0_whole, Etot0_el, Etot0
664     format(/14x,'FCHK elctrostatics INPUT FILE for geom-0   : ',a,/)
2012    format(/3x,'============ final Gaussian data ==========='
     &      ,/3x,'QM data read from:        ',a
     &      ,/3x,'Electrostatics read from: ',a
     &      ,/3x,'E(QM):                    ',f14.8
     &      ,/3x,'E(Electr.):               ',f14.8
     &      ,/3x,'E(QM-Electr.):            ',f14.8)
      endif  ! endif separate electrostatic (obsolete) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c     computes atoms bonds (in kbond)
      zord=.true.
      Call Atom_Count (nat,geom,aname,zord,maxat)
      do i=1,nat
        bname(i)=aname(i)
        Call Via_Numer(bname(i))
      enddo
      iout=6
      Call Get_Bonds (nat,geom,aname,dist,kbond,maxat,iout)

c     Set the array mass for each atom. 
c     The same mass can be used for all atoms.
      fmass=-1.0d0
      Call Cerca_W ('$mass',5,card,ier,.true.)
      if(ier.eq.0) read(card,*) fmass
      zmass=.false.
      if(fmass.gt.0.4999d0) zmass=.true.
      do i=1,nat
        if(zmass) then
          amass(i)=fmass
        else
          iz=geom(i,4)+1.d-09
          aa=Atom_Mass(iz)
          amass(i)=aa
        endif
        geom(i,5)=amass(i)
      enddo
      if (levprt .gt. 0 ) then 
         if(zmass) then
           write(6,18) fmass
         else
           write(6,19)
         endif
      endif
   18 format(/8x,'All nuclear masses are set to',f9.3)
   19 format(/8x,'Normal nuclear masses are used')

c     transform cartesian coord.s from Angstrom to bohr
c     and save (before UA changes) geom in /geom00/
      do i=1,nat
        do k=1,3
          geom(i,k)=geom(i,k)/boh2ang    !! bohr
          geom00(i,k)=geom(i,k)          !! bohr
        enddo
        geom00(i,4)=geom(i,4)            !! charge
        geom00(i,5)=geom(i,5)            !! mass
      enddo

c     read data for United Atoms in the molecule
      nator=nat
      Call United_Atom (geom,nat,nator,aname,kbond
     &                 ,UAM,manua,zstop,maxat,max3n)
c     save current reference geometry after UA changes in /geom0/
      ncart=nat*3
      do i=1,nat
        do k=1,3
          geom0(i,k)=geom(i,k)
          gang (i,k)=geom(i,k)*boh2ang
        enddo
        gang (i,4)=geom(i,4)
        gang (i,5)=geom(i,5)
        geom0(i,4)=geom(i,4)
        geom0(i,5)=geom(i,5)
      enddo
c     bond analysis (if UA is active)
      if(nat.lt.nator)
     &  Call Get_Bonds (nat,gang,aname,dist,kbond,maxat,iout)

c     punch file for visualization
      open (45,file='joyce.g0.xyz')
      rewind 45
      write(45,'(i8,"    from last joyce",/2x)') nat
      do i=1,nat
        write(45,'(a,3x,3f14.8)') bname(i),(gang(i,k),k=1,3)
      enddo
      close (unit=45)

c     compute array: sqrt(mass)
      m=0
      do i=1,nat
        aa=geom(i,5)
        wmass(m+1)=dsqrt(aa*amu)
        wmass(m+2)=wmass(m+1)
        wmass(m+3)=wmass(m+1)
        m=m+3
      enddo

      iprt = 0
c     generate the NATURAL internal coordinates
      Call Generate_IC (iprt,iflagn,ifftype,labnat,ncintn,
     &                  nat,kbond,geom,aname,maxat,maxpar,cr,
     &                  sigma,epsilon,charge,z_lj,
     &                  zharmdihed,zstardihed)




      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  'C) Reading FF & IC definition '
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  

c     read data from selected FF file
c     the parfun data are in bohr,radians
c     the LJ param is Sigma not R0
c     possible formats are:
c     1   moscito
c     2   gromacs
      ifftype = 0 
      Call Cerca_W ('$forcefield',5,card,ier,.true.)
c
      if(ier.ne.0) then
        write(6,24)
        stop
      endif
24    format(//' ERROR: keyword $forcefield not found'//)
c
      Call Selwrd (card,filetype,1,lwm,ipos)
      if (filetype(1:lwm) .eq. 'moscito') then
          ifftype = 1 
      elseif (filetype(1:lwm) .eq. 'gromacs') then
          ifftype = 2 
c     GP 3/10/2023 insert here new formats as lammps or amber
      else 
          write(6,244)
          stop
      endif
244   format(//' ERROR: format not specified in $forcefield '//)
      Call Selwrd (card,FFstartfile,2,lwm,ipos)
c
      if (ifftype .eq. 1) then    !! moscito format   
c        write(6,25)
         write(6,26) FFstartfile(1:lwm)
c5       format(//8x,44('*')/17x,'M O S C I T O    F F'/8x,44('*'))
26       format(/10x,'Moscito input file ..........: ',a,/)
         if(lwm.le.1) Call Fstop ('card $forcefield: no file given')
         iumosc=12
         open(iumosc,file=FFstartfile,form='formatted',
     &        status='old',err=2339)
         goto 2340
c
2339     write(6,'(a)') 'I/O error: file was not found'      
         write(6,'(a8,2x,a1,1x,a30)') 'Filename',':',FFstartfile
         stop
c
2340     continue      
         Call Read_Moscito (iumosc,iflag,par,parfun,npar,gang,nat
     &                     ,maxpar,maxat,aname,labpar,atlab,ichrgprm
     &                     ,nsit,namesit,cr,eps_sit,sig_sit
     &                     ,sitmass,charge,sigma,epsilon,z_lj)
c
      elseif (ifftype .eq. 2) then     !! gromacs format 
c        write(6,27)
         write(6,28) FFstartfile(1:lwm)
c7       format(//8x,44('*')/17x,'G R O M A C S    F F'/8x,44('*'))
28       format(/10x,'Gromacs input file ..........: ',a,/)
         if(lwm.le.1) Call Fstop ('card $forcefield: no file given')
         iugro=12
         open(iugro,file=FFstartfile,form='formatted',
     &        status='old',err=2341)
         goto 2342
c
2341     write(6,'(a)') 'I/O error: file was not found'      
         write(6,'(a8,2x,a1,1x,a30)') 'Filename',':',FFstartfile
         stop
c
2342     continue      
         Call Read_Gromacs (iugro,iflag,par,parfun,parfun2,npar,gang,nat
     &                     ,maxpar,maxat,aname,labpar,atlab,ichrgprm
     &                     ,nsit,namesit,cr,eps_sit,sig_sit
     &                     ,sitmass,charge,sigma,epsilon,z_lj)
c
      else
         write(6,'(5x,a)') 'Joyce input syntax error:'
         write(6,'(5x,a)') 'Force field format not specified'
         write(6,'(5x,a)') '$forcefield key contains:'
         write(6,'(5x,a)') card
         stop
      endif


c     fits LJ parameters also (assigned by default)

      z_fitlj = .false. 
      Call Cerca_W ('$fitLJ',5,card,ier,.true.)
      if(ier.eq.0) then
        z_fitlj = .true.
      endif

c     generates final parameters in AMBER format

      z_amber = .false. 
      Call Cerca_W ('$amber',5,card,ier,.true.)
      if(ier.eq.0) then
        z_amber = .true.
      endif

c     generates final parameters in G09 format

      z_g09 = .false. 
      Call Cerca_W ('$gaussian',5,card,ier,.true.)
      if(ier.eq.0) then
        z_g09 = .true.
      endif

c     punches NATURAL internal coordinates with FF atom labels

      iprt = 1
      Call Generate_IC (iprt,iflagn,ifftype,labnat,ncintn,
     &                  nat,kbond,geom,atlab,maxat,maxpar,cr,
     &                  sigma,epsilon,charge,z_lj,
     &                  zharmdihed,zstardihed)

c     generate the labels for each FF parameter

      Call Generate_Label (npar,iflag,labpar,atlab,maxpar
     &                    ,kbond,maxat)

c     generate file with dependencies suggested 
c     by symmetry of atomic labels

      Call Generate_deps (npar,iflag,labpar,atlab,maxpar
     &                    ,maxat)




c     stops the program after generation of 
c     suggested IC and deps
c     activated by the $generate key

      Call Cerca_W ('$generate',5,card,ier,.true.)
      if(ier.eq.0) then
          write(6,'(a)') 
          write(6,'(a)') 
     &     '*********************************************************'
          write(6,'(a)') 
          write(6,'(a)') 
     &     '$generate key has been activated in the input section:'
          write(6,'(a)') 
     &     ' no fitting will be performed'
          write(6,'(a)') 
          write(6,'(a)') 
     &     'Suggested Internal coordinates have been generated'
          write(6,'(a)') 
     &     'according to the provided site labels in file:' 
          write(6,'(a)') 
     &     '                    generated.IC.txt'
          write(6,'(a)') 
          write(6,'(a)') 
     &     '*********************************************************'
          write(6,'(a)') 
          write(6,'(a)') 
          write(6,'(a)') 
          write(6,'(a)') 
          write(6,'(//15x,"NORMAL EXIT",/)')
          write(6,'(a)') 
          stop 
      endif


c     keeps the R0 (or Ang0, or gamma) value
c     of FF file, even if this does not match the equil geom
c     since 9/2017  keepff has become a double vector to 
c     allow for coupling constants to be kept, as gamma1 and gamma2

      do i=1,npar
       iclass(i)=0
       do kk = 1,2
         ikeepff(i,kk)=0     !! 0=compute Par Fun  1=keep from FF file
       enddo
      enddo
      Call Cerca_W ('$keepff',5,card,ier,.true.)
      if(ier.eq.0) then
        Call Read_Seq (card,iclass,nk,1,npar)
        if(nk.gt.0) then
          do i=1,nk
            k=iclass(i)
            ikeepff(k,1)=1
          enddo
          write(6,42) (iclass(i),labpar(iclass(i)),i=1,nk)
        endif
      endif
   42 format(/5x,'The following functions keep the R0/Ang0/Gamma'
     &          ,' values as given in FF file',/(5x,i4,3x,a))

c     search for parameters to be kept in coupling (ikeepff(i,2))
      Call Cerca_W ('$Ckeepff',5,card,ier,.true.)
      if(ier.eq.0) then
        Call Read_Seq (card,iclass,nk,1,npar)
        if(nk.gt.0) then
          do i=1,nk
            k=iclass(i)
            ikeepff(k,2)=1
          enddo
          write(6,43) (iclass(i),labpar(iclass(i)),i=1,nk)
        endif
      endif
   43 format(/5x,'The following coupled functions will also keep', 
     &           'the R0/Ang0/Gamma',
     &           ' values as given in FF file',/(5x,i4,3x,a))



c     compute the values of RIC's at the reference geometry
c     valic0 in bohr,rad

      Call Value_IC (iflag,valic0,npar,geom,maxat,maxpar)

c     set the adequate R0,Phi0,Gamma values to each function,
c     in order to get a minimum of FF at the g03 refer. geometry
c     In the case ikeepff(i,j)=1 the force-field parameter is kept
c     and the FF may not have a minimum at the refer. geometry

      Call Set_ParFun (iflag,valic0,parfun,parfun2,npar,ikeepff,maxpar)


c     analyze the RICO's of the Force Field

      Call Anal_FF (iflag,par,valic0,parfun,labpar,npar
     &             ,'AT ABSOL-MIN',kbond,aname,maxat,maxpar)

c     set the Redund Int Coord.s (RICOs) from the FF parameters

      Call Recover_RIC (npar,iflag,labpar,indir
     &                 ,nric,iflagr,labric,maxpar)
      do ir=1,nric
        vc(ir)=boh2ang
        if(iflagr(ir,8).ne.2) vc(ir)=todeg
      enddo
      Call Anal_IC (iflagr,valic,labric,nric,'AT ABSOL-MIN'
     &             ,nat,kbond,geom,aname,maxat,maxpar)

c     compare the NATURAL ICs with the RICOs

      write(6,32) nric,ncintn
   32 format(//9x,'COMPARE:',i6,' RICOs   with',i6,' natural-ICs')
      lab1='Redund-IC '
      lab2='Natural-IC'
      Call Compare_IC (iflagr,labric,nric  ,lab1
     &                ,iflagn,labnat,ncintn,lab2,max3n,maxpar)

c     read option for the normal mode calculation

      zgrah=.true.
      Call Cerca_W ('$normal',5,card,ier,.true.)
      Call LowCas (card)
      if(index(card,'zero').gt.0) zgrah=.false.

c     read scaling factor for frequencies

      wfreq=-1.0d0
      Call Cerca_W ('$wfre',5,card,ier,.true.)
      if(ier.eq.0) read(card,*) wfreq

c     CONSIDER THE ABSOLUTE MINIMUM
      string='ABSOLUTE MINIMUM GEOMETRY'

      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  'D) Working with the QM Hessian'
      write(6,'(10x,a)')  '------------------------------'
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  

c     if requested rescales QM enrgy for LJ intra contrib. 
      z_LJassign = .false.
      Call Cerca_W ('$LJassign',5,card,ier,.true.)
      if(ier.eq.0) then
          write(6,'(10x,a)') 
     &     '$LJassign key has been activated in the input section:'
          write(6,'(10x,a35,1x,a41)') 
     &     'intramolecular LJ prms are assigned',
     &     'and their energy  subtracted from QM data'
          write(6,'(a)') 
          write(6,'(a)') 
          z_LJassign = .true.
      endif
c     Rescales for the intramolecular LJ term if requested ! GP & SG October 2023
      if (z_LJassign) then 
              z_scan = .false.
              euncorr0 = Etot0
              Call CorrectQMdata(Etot0, nator, geom, iflag, par,
     &             parfun, npar, z_scan, maxat, max3n, maxpar)
              write(6,2013) euncorr0, Etot0, (euncorr0 - Etot0)*2625.5d0
2013          format(/10x,'=========== final Gaussian data ===========',
     &        /10x,'Original QM energy (a.u.):     ', f15.8,
     &        /10x,'LJ-corrected QM energy (a.u.): ', f15.8,
     &        /10x,'LJ-correction (kJ/mol):        ', f15.8)
      endif

c     first of all, trasform matrices according to UA data
c     If there aren't UA, nator=nat and no action

      igeom=0
      whd=5000.0d0
      Call Cerca_W ('$whess',5,card,ier,.true.)
      if(ier.eq.0) read(card,*) whd
      who=whd/2.d0    
      Call Count_W (card,nw)
      if(nw.gt.1) read(card,*) whd,who
      whd=dabs(whd)
      who=dabs(who)
      if(znohess .and. (whd+who).gt.1.0d-12) then
        whd=0.0d0
        who=0.0d0
        write(6,'(//8x,"WARNING: Hess+Grad not considered because")')
        write(6,'(  8x,"       : no data were found in file fchk ",/)')
      endif
c
      do i=1,nric
        ip=indir(i)
        valric(igeom,i)=valic0(ip)
      enddo
c
      npt=0
      we=0.0d0
      wg=0.0d0

c     transform Grad,Hess according to the United Atom instructions

      Call Transf_HG_UA (Grad,Fg03,UAM,D2,nat,nator,max3n)

c     check: U(t+r)* * Hessian * U(t+r) = 0

      Call U_for_Trot (geom,nat,Utrot,ncart,maxat)
      do j=1,6
        do i=1,ncart
          D2(i,j)=Utrot(j,i)/wmass(i)
        enddo
      enddo
      zover=.false.
      n6=6
      Call Trsml2 (Fg03,D2,T,n6,ncart,max3n,zover)
      if (nator .ne. nat) then 
          write(6,62)
          do i=1,n6
            write(6,64) (T(i,j),j=1,n6)
          enddo
          write(6,*) ' '
   62     format(//9x,'=== IMPORTANT CHECK FOR THE UA TREATMENT ==='
     &       ,/9x,'=== the follow matrix should be null     ==='
     &       ,/9x,'Print:  Ut(t+r) * Hess * U(t+r) = 0 !'
     &       ,/9x,'units = (H/bohr2)'
     &       ,/9x,'1=tr.x 2=tr,y 3=tr.z 4=rot.1 5=rot.2 6=rot.3')
   64     format(3x,6f12.8)
      endif

c     compute the alpha,beta matrices for the reference geometry

      Call Comp_Syst (iflag,parfun,parfun2,Etot,Grad,Fg03,npar,string
     &               ,alpha,beta,weig,we,wg,whd,who,labpt,npt
     &               ,wfreq,igeom,valic0,Etot0,zgrah
     &               ,geom,nat,labpar
     &               ,iflagr,labric,nric
     &               ,maxpt,maxat,max3n,maxpar)

c     read option for Boltzmann weights for the energy points
c     $boltzmann 1.0 300.

      fboltz=-1.0d0
      boltzmin=0.01d0
      Call Cerca_W ('$boltz',5,card,ier,.true.)
      if(ier.eq.0) then
        fboltz=1.0d0
        Call Selwrd (card,word,1,lwi,ipos)
        if(word(1:1).ne.'*') read(word,*) fboltz
        tempe=300.0d0
        Call Selwrd (card,word,2,lwi,ipos)
        if(word(1:1).ne.'*') read(word,*) tempe
        if(tempe.lt.0.01d0) then
          write(6,'(/5x,"INPUT ERROR: negative temperature",/)')
          stop
        endif
      endif
      if(fboltz.gt.boltzmin) write(6,72) fboltz,tempe
   72 format(/8x,'The weight of the energies are given by the'
     &   ,' Boltzmann formula'
     &      ,/8x,'Multiplicative factor ...',f9.4
     &      ,/8x,'Temperature (K) .........',f9.4)


c     reads geometry files other then equilibrium (if any) and
c     computes Alpha and Beta matrices each
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  
      write(6,'(10x,a)')  'E) Working with relaxed QM scans'
      write(6,'(10x,a)')  '--------------------------------'

      Call Cerca_W ('$geom',5,card,ier,.true.)
      iscan(1)=1   !! scanned coordinate
      if(ier.ne.0) goto 200
      nfroz=0
      tolr=0.03d0  !! tolerance for frozen geometry changes: distances
      tola=3.00d0  !! tolerance for frozen geometry changes: angles
c
      iline=0
      iun5=5
  100 card=' '
      Call Read_One_Card (iun5,card,iline) 
      card0=card
      Call LowCas (card0)
      if(index(card0(1:4),'#').gt.0)    goto 100
      if(index(card0(1:4),'%').gt.0)    goto 100
      if(Lenw(card0).lt.3)              goto 100
      if(index(card0(1:9),'$end').gt.0) goto 200
      Call SetGrat (card)
c
      igeom=igeom+1          !! geometry index
      if(igeom .gt. 1) iscan(igeom)=iscan(igeom-1)    !! scanned coordinate
c
      Call Check_Dim 
     &     (igeom,maxgeom,'Too many geometries: increase maxgeom')
      Call Selwrd (card,fileinp,1,lw,ipos)
      word=fileinp
      Call LowCAs (word)
      lw=Lenw(word)
      if (word(lw-3:lw) .eq. '.fcc') then
        ftype = 'fcc '
      else
        ftype = 'fchk'
        ii=index(word,'.fchk')
        if(ii.le.0) fileinp(lw+1:lw+5)='.fchk'
      endif
c
      Call Selwrd (card,word,2,lw,ipos)
      we=dble(nat*(nat+1)/2)
      if(word(1:1).ne.'*') read(word,*) we
      Call Selwrd (card,word,3,lw,ipos)
      wg=0.0d0    
      if(word(1:1).ne.'*') read(word,*) wg
      Call Selwrd (card,word,4,lw,ipos)
      whd=0.0d0    
      if(word(1:1).ne.'*') read(word,*) whd
      who=0.0d0    
      Call Selwrd (card,word,5,lw,ipos)
      if(lw.gt.0 .and. word(1:1).ne.'*') read(word,*) who
      we=dabs(we)
      wg=dabs(wg)
      whd=dabs(whd)
      who=dabs(who)
c
      inquire (file=fileinp,exist=zexist)
      if(.not.zexist) then
        write(6,106) fileinp(1:Lenw(fileinp))
        stop
      endif
  106 format(/,/4x,'Error: QM file does not exist:  ',a,///)
      iun=32
      open(iun,file=fileinp,form='formatted',status='old')
      write(6,174) igeom,fileinp(1:Lenw(fileinp))
      levprt=levprt-1
      if (ftype .eq. 'fcc ') then
        Call Read_fcc  (nator,geom,cname,Etot,FCC,Grad,iun,maxat,max3n)
      else
        Call Read_fchk (nator,geom,cname,Etot,FCC,Grad,iun,maxat,max3n)
      endif


      if (z_LJassign) then   ! added by GP & SG october 2023
              z_scan = .true.
              euncorr = Etot
              Call CorrectQMdata(Etot, nator, geom, iflag, par,
     &             parfun, npar, z_scan, maxat, max3n, maxpar)
              write(6,2013) euncorr , Etot , (euncorr  - Etot )*2625.5d0
      endif

      deuncorr    =  (euncorr - euncorr0)*2625.5d0
      DEQM(igeom) =  (Etot-Etot0)*2625.5d0

      write(6,107) DEQM(igeom)
  107 format(10x,'E(tot) - E(reference geom) ........',f15.8
     &                                               ,f10.2, 'kJ/mol')
      levprt=levprt+1
      do i=1,nator
        geom(i,5)=geom00(i,5)   !! masses
      enddo

c     check the atom sequence

      do i=1,nator
        ad=dabs(geom(i,4)-geom00(i,4))
        if(ad.gt.1.d-08) then
          write(6,108) i,ad,(k,geom00(k,4),geom(k,4),k=1,nator)
          write(6,'(//9x)')
          stop
        endif
      enddo
  108 format(//9x,'ERROR: current atomic sequence is different'
     &       ,/9x,'     : from the reference sequence'
     &       ,/9x,'     : changed atom is',i5,d14.3
     &       ,/9x,'    ref       actual'
     &       ,/(i5,2f12.4))

c     change angstrom --> bohr

      do i=1,nat
        do k=1,3
          geom(i,k)=geom(i,k)/boh2ang    !! bohr
          geoms(i,k)=geom(i,k)
        enddo
        geoms(i,4)=geom(i,4)
      enddo

c     transform geometry according to the United Atom instructions

      Call Transf_GE_UA (geom,manua,nat,nator,aname,maxat)


c     read data for frozen geometrical changes of the reference geom,
c     on the basis of one (or more) internal coordinate read from g03
c     WARNING:  the geometry to be changed is the UA geometry
c
c          filename 1.0 1.0 0.0 0.0  ; 4 6 8  ; 3 4 ;


      if((wg+whd+who).gt.1.d-04) goto 150

c     case: igeom>1   use the previusly read data

      if(igeom.gt.1 .and. nfroz.gt.0) then
        irep=index(card,'-f')
        if(irep.gt.0) goto 130
        irep=index(card,'-F')
        if(irep.gt.0) goto 130
      endif
      if(igeom .gt. 1) iscan(igeom) = iscan(igeom) + 1
c
      lwc=Lenw(card)
      nfroz=0
  120 ii=index(card,';')
      if(ii.le.1)               goto 130
      card0=' '
      card0=card(ii+1:lwc)
      jj=index(card0,';')
      if(jj.gt.1) word=card0(1:jj-1)
      if(jj.le.1) word=card0(1:lwc)
      card(ii:ii)=' '
      lww=Lenw(word)
      Call Count_W (word,nw)
      if(nw.gt.4) nw=4
      if(nw.eq.1) then
        write(6,'(//" ERROR in Froz Geom Change: too few atoms")')
        write(6,'(12x,a,//)') word(1:lww)
        stop
      elseif(nw.eq.0) then
        goto 130
      endif
      nfroz=nfroz+1
      do k=0,8
        iflagrf(nfroz,k)=0
      enddo
      read(word,*) (iflagrf(nfroz,k),k=1,nw)
      iflagrf(nfroz,0)=nw
      iflagrf(nfroz,8)=nw
      do k=1,nw
        ia=iflagrf(nfroz,k)
        if(ia.lt.1 .or. ia.gt.nat) then
          write(6,122) ia
          stop
        endif
      enddo
      goto 120
  122 format(//8x,'ERROR in frozen geom input: atom <1 or >nat'
     &       ,/8x,'atom number',i8,//)

c     end read,  make geometrical changes

  130 continue
      if(nfroz.le.0) goto 150
      do i=1,nric
        inda(i)=0
      enddo
      write(6,134)
      write(6,135)
      do ifr=1,nfroz
        ia=iflagrf(ifr,1)
        ib=iflagrf(ifr,2)
        ic=iflagrf(ifr,3)
        id=iflagrf(ifr,4)
        nc=iflagrf(ifr,8)
        ifroz(ifr)=ifr
        Call Value_1IC (geom ,maxat,ia,ib,ic,id,val,'r')  ! @relaxed  QM geom
        Call Value_1IC (geom0,maxat,ia,ib,ic,id,val0,'r') ! @full opt QM geom
        if(nc.gt.2) Call Shifta (val,val0)
        fco=boh2ang
        if(nc.gt.2) fco=todeg
        xfroz(ifr)=val*fco
        x1=val*fco
        x0=val0*fco
        write(6,136) ifr,ia,ib,ic,id,x1,x0,x1-x0
        Call Gia_ic (nric,iflagr,ia,ib,ic,id,ichi,maxpar)
        if(ichi.gt.0) inda(ichi)=ifr
      enddo


c     Prints outs QM energy scans

      if (igeom .eq. 1 ) then
             write(gnum0,'(i3)') iscan(igeom)
             Call Selwrd (gnum0,gnum,1,lw,ipos)
             fileqmnrg = 'qmscan.'//gnum(1:lw)//'.dat'
             open(74,file=fileqmnrg,form='formatted',status='unknown')
             if (z_LJassign) then 
                 write(74,'(a)') 
     &           '#  coordinate   DEcorr (kj/mol)   DE (kj/mol) '
             else
                 write(74,'(a)') '#  coordinate   DE (kj/mol)'
             endif
      elseif (iscan(igeom) .ne. iscan(igeom-1)) then
             close(74)
             write(gnum0,'(i3)') iscan(igeom)
             Call Selwrd (gnum0,gnum,1,lw,ipos)
             fileqmnrg = 'qmscan.'//gnum(1:lw)//'.dat'
             open(74,file=fileqmnrg,form='formatted',status='unknown')
             if (z_LJassign) then 
                 write(74,'(a)') 
     &           '#  coordinate   DEcorr (kj/mol)   DE (kj/mol) '
             else
                 write(74,'(a)') '#  coordinate   DE (kj/mol)'
             endif
      endif
      rc_scanned = xfroz(1)
      if (iflagrf(1,0) .eq. 4) then !! dihedral
         if ( xfroz(1) .le. -180.d0) then
              rc_scanned = rc_scanned + 360.d0
         elseif ( xfroz(1) .gt. 180.d0) then
              rc_scanned = rc_scanned - 360.d0
         endif     
      endif     
      if (z_LJassign) then 
       write(74,'(2x,f9.3,4x,f12.4,f12.4)') 
     &          rc_scanned,DEQM(igeom),deuncorr
      else
       write(74,'(2x,f9.3,4x,f12.4)')  rc_scanned,DEQM(igeom)
      endif

c     now change the geometry on the basis of the few RICs extracted
c     from the relaxed geometry read from the last .fchk/.fcc file

      Call Change_Geom (iflagrf,nric,nfroz,xfroz,ifroz
     &                 ,geom,geom0,aname,nat,kbond,maxat,maxpar)

c     test: are the other RIC practically unchanged in the real
c           geometry and in the geometry obtained by frozen changes ?

      Call Value_IC (iflagr,ricf,nric,geom ,maxat,maxpar)
      Call Value_IC (iflagr,rics,nric,geoms,maxat,maxpar)
      do i=1,nric
        if(iflagr(i,8).gt.2) Call Shifta (ricf(i),rics(i))
      enddo
      if(levprt.gt.1) then
        write(6,140)
        do i=1,nric
          fc=boh2ang
          if(iflagr(i,8).gt.2) fc=todeg
          dif=ricf(i)-rics(i)
          write(6,142) i,fc*ricf(i),fc*rics(i),fc*dif
     &                ,(iflagr(i,k),k=1,4),inda(i)
        enddo
      endif
      nwrong=0
      do i=1,nric
        if(inda(i).eq.0) then
          dif=ricf(i)-rics(i)
          if(iflagr(i,8).gt.2)then
            fc=todeg
            tol=tola
          else
            fc=boh2ang
            tol=tolr
          endif
          dif=dif*fc
          if(dabs(dif).gt.tol) then
            nwrong=nwrong+1
            if(nwrong.eq.1) write(6,144)
            err=dabs(dif)/tol
            write(6,146) i,fc*ricf(i),fc*rics(i),dif
     &                  ,(iflagr(i,k),k=1,4),tol,err
          endif
        endif
      enddo
c
c
      if(levprt.ge.0)
     &write(6,138) (i,aname(i),(geom(i,k)*boh2ang,k=1,4),i=1,nator)
c
  134 format(/10x,'The geometry is changed according to the FROZEN'
     &                                       ,' options'
     &      ,/10x,'1. the follow. RICs are obtained by the current'
     &                                       ,' QM geom '
     &      ,/10x,'2. the reference geom is changed accordingly'
     &      ,/10x,'   and the obtained geom is used in the follow')
  135 format(/10x,'  RIC    --- atoms ----'
     &               ,5x,'current    refer   change')
  136 format(10x,i4,3x,4i4,3x,3f9.3)
  138 format(/10x,'----- GEOMETRY from frozen changes (angstrom) -----'
     &      ,/28x,'x',11x,'y',11x,'z',7x,'Nucl.ch'
     &      ,/(10x,i3,2x,a5,3f12.7,f8.1))
  140 format(/10x,'Value of the RICs in the geometry obtained'
     &        , ' by frozen changes'
     &      ,/10x,'FrozGeo   RelaxGe    change')
  142 format(10x,i4,3f10.3,3x,4i4,i6)
  144 format(/10x,'WARNING: the following RICs are changed too much' 
     &      ,/10x,'FrozGeo   RelaxGe    change',21x,'allowed')
  146 format(10x,i6,3f10.3,3x,4i4,2f8.3)


c     end of frozen geometrical changes section


  150 continue
      do i=1,nator
        do k=1,3
          gang(i,k)=geom(i,k)*boh2ang    !! angstrom
        enddo
      enddo
      close (unit=32)
  174 format(//10x,40('='),/21x,'GEOMETRY n.',i3,/10x,40('=')
     &      ,/10x,'QM input file',3x,a,/)
      write(string,'(a,i3)') 'GEOMETRY =',igeom
c     checks Hessian have been read
      znohess=.true.
      do i=1,nat*3
        do j=1,nat*3
          if(dabs(FCC(i,j)).gt.1.d-08) znohess=.false.
        enddo
      enddo
c     checks Gradients have been read
      znograd=.true.
      do i=1,nat*3
        if(dabs(Grad(i)).ne.-99999.d0) znograd=.false.
      enddo
c
      Call Value_IC (iflag,valic,npar,geom,maxat,maxpar)
      do i=1,nric
        ip=indir(i)
        valric(igeom,i)=valic(ip)
      enddo
c
      Call Anal_FF (iflag,par,valic,parfun,labpar,npar
     &             ,string,kbond,aname,maxat,maxpar)
c
      if(znohess .and. (whd+who).gt.1.0d-12) then
        whd=0.0d0
        who=0.0d0
        write(6,'(//8x,"ERROR  : Hess     not considered because")')
        write(6,'(  8x,"       : no data were found in file fchk ",/)')
        stop 'check g09 files'
      endif

      if(znograd .and. (wg).gt.1.0d-12) then
        write(6,'(//8x,"ERROR  : Gradient not considered because")')
        write(6,'(  8x,"       : no data were found in file fchk ",/)')
        stop 'check g09 files'
      endif

c     transform Grad,Hess according to the United Atom instructions

      Call Transf_HG_UA (Grad,FCC,UAM,D2,nat,nator,max3n)
c
      if(fboltz.gt.boltzmin) then
        we=dexp(-(etot-etot0)/boltzau/tempe)
      endif

c     compute the new part of the Alpha,Beta matrices

      Call Comp_Syst (iflag,parfun,parfun2,etot,Grad,FCC,npar,string
     &               ,alpha,beta,weig,we,wg,whd,who,labpt,npt
     &               ,wfreq,igeom,valic0,etot0,zgrah
     &               ,geom,nat,labpar
     &               ,iflagr,labric,nric
     &               ,maxpt,maxat,max3n,maxpar)
      goto 100
  200 continue
      write(6,'(//20x,"END OF ALPHA,BETA CALCULATION AT GEOM.s",//)')


c     punch and analyze the RICO at the considered geometries

      ngeom=igeom
      if(ngeom.gt.1) then
        Call Punch_RICO (iflagr,labric,nric,valric,ngeom
     &                  ,maxpar,maxgeom)
      endif


c     printout: points

      if(levprt.ge.1) then
        write(6,205)
        nk=4
        iop=0
        do i=1,npt
          do k=1,npar
            work(k)=Alpha(i,k)
          enddo
          Call Sort_Max (work,anda,inda,npar,nk)
          write(6,206) i,labpt(i),weig(i),beta(i)*1000.
     &                ,(inda(k),1000.d0*anda(k),k=1,nk)
        enddo
      endif
  205 format(//9x,'Summary of the points to be fitted with Force Field'
     &       ,/7x,'--label--   weight   value(m.au)'
     &       ,5x,'most important functions (in m.au)')
  206 format(i4,3x,a,2f10.5,2x,5(i3,1x,f9.3))


c     solve the linear system (with dependencies)

      do ip=1,npar
        iclass(ip)=iflag(ip,0)
      enddo
      if(levprt.ge.1) then
        write(6,208)
        do ip=1,npar
          it=iflag(ip,0)
          write(6,'(i4,2x,a,6i4)') ip,stype(it),(iflag(ip,k),k=0,5)
        enddo
      endif
  208 format(//9x,'Summary of the considered parameters')
c
      Call Solve_Syst (par,npar,iflag,iclass,labpar,ichrgprm
     &                ,z_fitlj,alpha,beta,weig,labpt,npt
     &                ,maxpar,maxpt)



c     computes the vibrational normal modes at the equil.
c     geometry using the optimized FF

 
c     compute the CC Hessian by direct numerical derivatives
c     and by analytical method

      Call FF_Hessian (par,iflag,parfun,parfun2,T,FCC,npar
     &                ,geom0,nat,iflagr,nric,labpar
     &                ,maxpt,maxat,max3n,maxpar)

c     analyze FF Hessian

      null=0
      do j=1,ncart
        aa=0.0d0
        do i=1,ncart
          aa=aa+dabs(T(i,j))
        enddo
        if(aa.lt.1.d-10) null=null+1
      enddo
      if(null.gt.9) then
        write(6,240) null
        goto 300
      endif
  240 format(/12x,'Too many null rows are found:',i4
     &      ,/12x,'The FF is unable to describe normal modes at g0'
     &      ,/12x,'Frequencies are not computed',/)

c     compute the transformation matrix (RIC = U * CC)

      mx=maxpar
      Call Trasf2Int (geom0,nat,Un,iflagr,nric,maxat,max3n,maxpar)
      card='Transformation U Matrix: RIC = U * CC (bohr,rad)'
      if(levprt.ge.3) Call Pr_Rec2 (Un,nric,ncart,card,ncol,max3n)

c     compute normal modes from the optimized FF

      scelta='by the optim. FF (EQUIL GEOMETRY)'
      zpic=.false.
      Call Comp_Vibr (T,Un,F,E,Freq,C,geom0,scelta,zpic,zgrah
     &               ,nric,nat,nvib,labric,max3n,maxat)

c     read the normal modes at equil geom with the g09 Hessian
c     and compare them with those arising from the Hessian-FF
c     note that when electrostatic term is handled separately, 
c     real frequencies are still compared (i.e. corrected 
c     for the electrosatic term)

      if (z_sepel) then
        write(6,269)
        scelta='by the optim. FF + elctr. term '
        do j=1,ncart
         do i=1,ncart
           T(i,j) = T(i,j) + FCC_el(i,j)
         enddo
        enddo
        Call Comp_Vibr (T,Un,F,E,Freq,C,geom0,scelta,zpic,zgrah
     &                 ,nric,nat,nvib,labric,max3n,maxat)
        scelta='by the whole QM (EQUIL GEOMETRY)'
        zpic=.false.
        Call Comp_Vibr (FCC_whole,Un,F,work,Freq,T,geom0,scelta,zpic
     &               ,zgrah,nric,nat,nvib,labric,max3n,maxat)
        mx=max3n
        Call Compare_NM (z_plotfreq,E,work,C,T,D2,nvib,
     &                   ncart,iclass,max3n)
        zpic=.false.
      else
         rewind 48
         read(48) kn
         if(kn.gt.0)  then
           rewind 48
           read(48) kn,km,(work(i),i=1,nvib)
           read(48) ((T(i,nm),i=1,ncart),nm=1,nvib)
           close (unit=48)
           mx=max3n
           Call Compare_NM (z_plotfreq, E,work,C,T,D2,
     &                      nvib,ncart,iclass,max3n)
         else
           write(6,270)
         endif
      endif
  269 format(/9x,'=============================================='
     &      ,/9x,'COMPARISON BETWEEN G09 AND FF NORMAL MODES'
     &      ,/9x,'electrostatic has been re-added to both terms'
     &      ,/9x,'==============================================')
  270 format(/9x,'=============================================='
     &      ,/9x,'NO COMPARISON BETWEEN G09 AND FF NORMAL MODES'
     &      ,/9x,'reason: no G09 data were calculated at geom0'
     &      ,/9x,'==============================================')

c     punch results for MD progs (Moscito, Gromacs)

300   continue

      if (ifftype .eq. 1) then    !! moscito format   
        Call Print_Moscito (iumosc,iflag,par,parfun,npar,geom0
     &                      ,nat,maxpar,maxat,aname,labpar)
      elseif (ifftype .eq. 2) then    !! gromacs format   
        Call Print_Gromacs (iugro,iflag,par,parfun,parfun2,npar,geom0
     &                      ,nat,maxpar,maxat,aname,labpar)
      endif
c     punch results for Amber

      if(z_amber) then
         Call Print_Amber (iflag,par,parfun,npar
     &                     ,nat,maxpar,maxat,aname,labpar,
     &                      nsit,namesit,sitmass,charge)

      endif

c     punch results for G09 MM optimizations

      if(z_g09) then
         Call Print_Gaussian (iflag,par,parfun,npar,
     &                        nat,maxpar,maxat,aname,labpar,
     &                        nsit,atlab,sitmass,charge,geom0,
     &                        sig_sit,eps_sit,namesit)

      endif

c     Read data for energy scan along internal coord.s

      thrdx=1.d-8
      zrew=.true.
      icountsc=0
      iunsc=46
      iungg=47
  600 Call Cerca_W ('$scan',5,card,ier,zrew)
      zrew=.false.
      if(ier.ne.0) goto 700
      icountsc=icountsc+1
      Call Read_Opt (card,'-n',zprt)
      zprt=.not.zprt

c     prepare the grace files if requested

      if (z_plottors) then 
         if (icountsc .lt. 10) write(gnum0,'(i1)') icountsc
         if (icountsc .lt. 100) write(gnum0,'(i2)') icountsc
         call Compact(gnum0,lw)
         lwg= Lenw(torsgracefile0)
         torsgracefile = 
     &   'joyce.'//torsgracefile0(1:lwg)//'_'//gnum0(1:lw)//'.agr'
         iungrace = 179 + icountsc
         open(iungrace,file=torsgracefile,form='formatted',
     &        status='unknown')
         Call PrintGtors(iungrace,icountsc) 
      endif

c     read the filename where scan results will be punched

      Call Selwrd (card,filesc,1,lw,ipos)
      if(icountsc.eq.1) then
        if(lw.le.1) then
          Call Getarg (1,word)
          filesc='joyce.'//word(1:Lenw(word))//'.scan.pun'
        endif
        open(iunsc,file=filesc,form='formatted',status='unknown')
      else
        if(lw.gt.1) then
          open(iunsc,file=filesc,form='formatted',status='unknown')
        endif
      endif

c     read the filename for the reference geometry

      Call Selwrd (card,filegg,2,lw,ipos)
      if(lw.gt.1) then
        iungg=47
        open(iungg,file=filegg,form='formatted',status='unknown')
      else
        iungg=0
        filegg='standard geometry  '
        do i=1,nat
          do k=1,5
            geoms(i,k)=geom0(i,k)
          enddo
        enddo
      endif
c
      write(6,604) filesc(1:Lenw(filesc)),filegg(1:Lenw(filegg))
  604 format(/
     & ,/9x,'======================================================'
     & ,/9x,'===     ENERGY SCAN ALONG INTERNAL COORDINATES     ==='
     & ,/9x,'===  Punch results: ',a                      ,t61,'==='
     & ,/9x,'===  Pefer. geom  : ',a                      ,t61,'===')
      if (z_plottors) then 
         write(6,605) torsgracefile(1:Lenw(torsgracefile))
  605    format(
     &   9x,'===  Grace file   : ',a,t61,'===')
      endif
      write(6,606)
  606 format(
     &   9x,'======================================================')



      if(iungg.gt.0) Call Read_xyz (nat,geoms,aname,maxat,iungg)
      nscan=0
  610 read(5,'(a)',end=620) card
      card0=card
      Call LowCas (card0)
      if(index(card0(1:4),'#').gt.0)    goto 610
      if(index(card0(1:4),'%').gt.0)    goto 610
      if(Lenw(card0).lt.3)              goto 610
      if(index(card0(1:9),'$end').gt.0) goto 620
      Call SetGrat (card)
      lc=Lenw(card)

c     analyze the card
c     example:  6 12 7 ;  30. 60. 10

      ii=index(card,';')
      jj=index(card,',')
      if(ii.le.0 .and. jj.gt.1) ii=jj
      if(ii.le.0) goto 610
      Call Count_W (card(1:ii-1),nw)
      if(nw.le.1) then
        write(6,'(//" ERROR in SCAN: too few atoms")')
        write(6,'(12x,a,//)') word(1:lww)
        stop
      endif
      ic=0
      id=0
      if(nw.eq.2) read(card,*) ia,ib
      if(nw.eq.3) read(card,*) ia,ib,ic
      if(nw.eq.4) read(card,*) ia,ib,ic,id
      Call Gia_ic (nric,iflagr,ia,ib,ic,id,ichi,maxpar)
      if(ichi.eq.0) then
        write(6,'(//" ERROR in SCAN: wrong int coord")')
        write(6,'(12x,4i5,//)') ia,ib,ic,id
        stop
      endif
      Call Count_W (card(ii+1:lc),nw)
      if(nw.le.0) then
        write(6,'(/" ERROR: after ; thre are no data")')
        write(6,'(8x,"|",a,"|",/,/)') card(1:lc)
        stop
      endif

c     the IC is fixed

      if(nw.lt.3) then
        read(card(ii+1:lc),*) amin
        amax=0.0d0
        dx=0.0d0
      else
        read(card(ii+1:lc),*) amin,amax,dx
      endif
      nscan=nscan+1
      scamin(nscan)=amin
      scamax(nscan)=amax
      scainc(nscan)=dabs(dx)
      if(amin.gt.amax .and. dabs(dx).gt.thrdx) scainc(nscan)=-dabs(dx)
      iscaic(nscan)=ichi
      goto 610
  620 continue
c
      Call Do_Scan (par,parfun,parfun2,iflag,npar,zprt
     &             ,iflagr,labric,nric,iunsc
     &             ,scamin,scamax,scainc,iscaic,nscan
     &             ,z_plottors,iungrace
     &             ,geom,geoms,aname,nat,kbond,maxat,maxpar)
c
c


      goto 600
  700 continue


      write(6,'(//9x,"NORMAL EXIT",/)')
      stop
      end



