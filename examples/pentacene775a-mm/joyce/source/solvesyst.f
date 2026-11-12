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
      Subroutine Solve_Syst (par,npar,iflag,iclass,labpar,ichrgprm
     &                      ,z_fitlj,alpha,beta,weig,labpt,npt
     &                      ,maxpar,maxpt)
c----------------------------------------------------------------------c
 
 
c----------------------------------------------------------------------c
c     solve the linear system from the alpha,beta,weig arrays          c
c                                                                      c
c  par     o value of the parameters                                   c
c  npar    i # of parameters                                           c
c  iclass  i an array to identify the class of each parameter          c
c            To avoid splitting in classes,                            c
c            give iclass(1)=...=iclass(npar)=0                         c
c            Usually stretching, bending and torsional parameters      c
c            belong to different classes.  (iclass(i) == iflag(i,0))   c
c  labpar  i label of the parameters                                   c
c  ichrgprm i starting (ichrgprm[1]) and final (ichrgprm[2])           c
c             number identifying charge-charge interactions            c
c             which are ALWAYS assigned and NOT optimized              c
c  alpha   i alpha(mu,j) = funz. j at the point mu.    dim=(npt,npar)  c
c  beta    i beta (mu)   = value to be fitted in the point mu          c
c  weig    i a weight of each point                                    c
c  labpt   i a label for each point                                    c
c  npt     i # of points                                               c
c  maxpar  i dimension of some two-dim arrays in the caller            c
c  maxpt   i dimension of same two-dim arrays in the caller            c
c                                                                      c
c   This routine try to read some data from unit=5. These data are     c
c----------------------------------------------------------------------c

      implicit real*8(a-h,o-z)
      parameter (max1=2000, max2=20000)
      dimension iclass(*),par(*)
     &         ,alpha(maxpt,*),beta(*),weig(*)
      dimension ichrgprm(2)
      dimension iflag(maxpar,0:8)
      dimension idep(max1),rdep(max1),ifree(max1)
     &         ,idepa(max1),rdepa(max1)
     &         ,iass(max1),pass(max1),ileave(max1),kleave(max1)
     &         ,A(max1,max1),B(max1,1),X(max1,1),betas(max2)
     &         ,G(max1,max1),Q(max1,1),Sens(max1)
     &         ,T1(max1,max1),T2(max2,max1),alphas(max2,max1)
     &         ,iwra(max1),iwrb(max1),wwab(max1)
      character labpar(*)*(*),labpt(*)*(*)
     &         ,title*70,card*140,word*90,labi*9
      logical zrew,zmosc,z_fitlj
      common /printout/ levprt
      equivalence (T1(1,1),T2(1,1))
 
c     conversion factors
c            1 H/rad^2      = 2625.5d0 kJ/mol/rad^2
c            1 H/(Bohr*rad) = 10^4/2.01553 kJ/(mol*A*rad)
c            1 H/Bohr       = 10^4/2.01553 kJ/(mol*A)
c            1 H/(Bohr^2)   = 10^4/1.066   kJ/(mol*A^2)
c    fene    1 H            = 2625.5               kJ/mol
c    fgrd    1 H/(Bohr)     = 2625.5/0.529177      kJ/mol/A
c    fhes    1 H/(Bohr^2)   = 2625.5/(0.529177)^2  kJ/mol/A^2
 
      fene   = 2625.5
      fgrd   = 2625.5/0.529177
      fhes   = 2625.5/0.529177/0.529177

 
      ncol=10
      npar1=npar+1
      Call Check_Dim (npar1,max1,'SOLVE_SYST: increase max1')
      Call Check_Dim (npt  ,max2,'SOLVE_SYST: increase max2')
 
      write(6,'(//3x)')
      write(6,'(10x,a)') ' F) QMD-FF parameterization' 
      write(6,'(10x,a)') '----------------------------------'
      write(6,15) npar,npt
   15 format(/10x,'================================================='
     &      ,/10x,'  S O L V E   T H E   L I N E A R   S Y S T E M'
     &      ,/10x,'  F O R   T H E   B E S T   P A R A M E T E R S'
     &      ,/10x,'================================================='
     &    ,/,/10x,'    n. of parameters ..............',i5
     &      ,/10x,'    n. of points ..................',i5,/)
 
c     sets a factor for alpha,beta (for printout)
 
      aijmax=0.0d0
      do mu=1,npt
        do ip=1,npar
          alij=dabs(alpha(mu,ip))
          if(alij.gt.aijmax) aijmax=alij
        enddo
      enddo
      ff=10.d0/(aijmax+1.d-15)
      if(ff.lt.5.0d0) then
        ff=1.0d0
      else
        aa=1.0d0
        dowhile (aa.lt.ff)
          aa=aa*10.d0
        enddo
        ff=aa
      endif
      alphafac=ff
 
c     prints the alpha and beta matrices
 
      if(levprt.ge.3) then
        do mu=1,npt
          do i=1,npar
            T2(mu,i)=alpha(mu,i)*alphafac
          enddo
        enddo
        do mu=1,npt
          T2(mu,npar+1)=beta(mu)*alphafac
        enddo
        nab=npar+1
        write(word,'("  (*",1pd10.2,")")') alphafac
        title='ALPHA matrix + BETA in the last column'//word(1:15)
        Call Pr_Rec2 (T2,npt,nab,title,ncol,max2)
      endif

c ---------------------------------------------------------------
c                       DEPENDENCIES
c ---------------------------------------------------------------
 
 
c     reads the dependencies: some param.s are constrained 
c     to take a value dependent from the value of other param.s
c                           or
c     two o more param.s must be in a definte ratio
c
c     |$dependencies   0.999  |
c     | 4 = 18 * 1.20         |   P4=P18*1.2
c     | 6 = 7  * 1.00         |   P6=P7
c     |$end                   |
c
c     The first card includes a number /t/ used for automatic dependencies
c     between redundant parameters. If the dot product (fa.fb)>t then
c     the param. Pb is set as dependent from Pa: Pb = +/- Pa depending on the
c     sign of the dot product. The default value is 0.99.
c     If t>1 no dependency is automatically included.
c     If the $dependence card is not present, the default value is assumed.
c
c     In the second card P18 is free, whereas P4 is constrained
c     Really: both parameters are subjected to a constraint
c
c     This section forms the idep array:
c     if idep(I)=J then par(I)=par(J)*rapp(I) , par(i) is dependent
 
      do ip=1,npar
        idep(ip)=0
        ifree(ip)=0
      enddo
      thrdep=0.995d0
      write(6,16)
   16 format(/10x,'----------- read dependences -----------')
      iun=5
      zrew=.true.
      Call Cerca_W ('$depen',iun,card,ier,zrew)
      ndep=0
      if(ier.ne.0) goto 30
      Call Selwrd (card,word,1,lw,ipos)
      if(lw.gt.1) read(word,*) thrdep
      if(thrdep.lt.0.5d0) thrdep=0.5d0
   20 read(5,'(a)',end=30) card 
      Call LowCas (card)
      if(Lenw(card).lt.2)              goto 30
      if(index(card(1:9),'$end').gt.0) goto 30
      if(index(card(1:4),'#').gt.0)    goto 20
      if(index(card(1:4),'%').gt.0)    goto 20
      Call SetGrat (card)
      ndep=ndep+1
      do m=1,Lenw(card)
        if(card(m:m).eq.'*') card(m:m)=' '
        if(card(m:m).eq.'=') card(m:m)=' '
      enddo
 
      Call Selwrd (card,word,3,lw,ipos)
      if(lw.eq.0) then
        write(6,32)
        stop 
      endif
      read(card,*) ia,ib,rapp
      if(ia.gt.npar .or. ia.lt.0) then
        write(6,33) ia
        stop 
      endif
      if(idep(ia).gt.0) then
        write(6,34) ia
        stop 
      endif
      idep(ia)=ib
      rdep(ia)=rapp
      if(iclass(ia).ne.iclass(ib)) then
        write(6,35) card(1:Lenw(card))
        stop
      endif

c     test: dependencies on already dependent param.

      if(idep(ib).eq.0) then 
        write(6,36) ndep,ia,ib,rapp
        ifree(ib)=ia
      else
        ic=idep(ib)
        idep(ia)=ic
        rdep(ia)=rdep(ia)*rdep(ib)
        write(6,37) ndep,ia,ib,rapp,ic,rdep(ia)
      endif

c     test: conflicting dependences

      do ip=1,npar
        if(ia.eq.idep(ip)) then
          write(6,38) ia
          stop
        endif
      enddo

c     test: self dependence

      if(ia.eq.ib) then
        write(6,39)
        stop
      endif
 
      goto 20
   30 continue
      if(ndep.eq.0) write(6,'(12x,"no given dependences")')
 
   32 format(//10x,'DEPENDENCE ERROR: too few data',/)
   33 format(//10x,'DEPENDENCE ERROR: index out of range',i4,/)
   34 format(//10x,'DEPENDENCE ERROR: two dependencies of param.',i5,/)
   35 format(//10x,'DEPENDENCE ERROR: the involved param.s'
     &     ,' belong to different classes',/3x,a//)
   36 format(i4,10x,'param',i4,3x,'=   param',i4,'  *',f10.4)
   37 format(i4,10x,'param',i4,3x,'=   param',i4,'  *',f10.4,' ='
     &         ,10x,'param',i4,'  *',f10.4)
   38 format(//10x,'DEPENDENCE ERROR: depend. par was already used as'
     &           ,' a refer. param',i5,/)
   39 format(//10x,'DEPENDENCE ERROR: one param is dep. by itself ?',/)
 
c ---------------------------------------------------------------
c        ASSIGNED PARAMETERS
c ---------------------------------------------------------------
 
c     WARNING: charge-charge parameters are always assigned to
c              values read in the system/topology file
c     WARNING: intramolecular LJ parameters are assigned for default
c              unless the $fitLJ keyword is activated in the input section
 
      nass=0
      do ip=1,npar
        iass(ip)=0
        if (ip .ge. ichrgprm(1) .and. ip .le. ichrgprm(2))  then
             nass = nass + 1
             iass(ip) = 1
             pass(ip) = par(ip)  ! note par(ip) is already in H*bohr
        endif        
        if (.not. z_fitlj) then
            if(iclass(ip) .eq. 5) then  !LJ
             nass = nass + 1
             iass(ip) = 1
             pass(ip) = par(ip)  ! note par(ip) is already in H*bohr
            endif
        endif
      enddo
      iun=5
      zrew=.true.
      Call Cerca_W ('$assign',iun,card,ier,zrew)
      if(ier .ne. 0 ) goto 60
      zmosc=.true.
      Call LowCas (card)
      write(6,46) 
   46 format(/10x,'------- read assigned parameters -------' )
      write(6,47) 
   47 format(10x,'------- Expected input units:    -------' )
      write(6,48) 
   48 format(10x,'------- [L] =A ; [E] = kJ/mol    -------' )
      if(ier.ne.0) goto 60
   50 read(5,'(a)',end=60) card 
      Call LowCas (card)
      if(Lenw(card).lt.2)              goto 60
      if(index(card(1:9),'$end').gt.0) goto 60
      if(index(card(1:4),'#').gt.0)    goto 50
      if(index(card(1:4),'%').gt.0)    goto 50
      Call SetGrat (card)
      do m=1,Lenw(card)
        if(card(m:m).eq.'=') card(m:m)=' '
      enddo
 
      Call Selwrd (card,word,2,lw,ipos)
      if(lw.eq.0) then
        write(6,61)
        stop 
      endif
      read(card,*) ia,valu
      if(ia.gt.npar .or. ia.le.0) then
        write(6,63)
        stop 
      endif
      nass=nass+1
      iass(ia)=1
      pass(ia)=valu
      if(zmosc) then
        n0=iclass(ia)
        pass(ia)=valu/fene
        if(n0.eq.1) pass(ia)=valu/fhes
        if(n0.eq. 8) then  ! coupling requested
          ic1 = iflag(ia,1)
          ic2 = iflag(ia,2)
          it1 = iflag(ic1,0) 
          it2 = iflag(ic2,0) 
          if (it1 .eq. 1 ) then  ! stretching coupled with 
             if (it2 .eq. 1 ) then  ! stretching 
                 pass(ia) = valu/fhes
             else                   ! bending or torsion
                 pass(ia) = valu/fgrd
             endif
          else                   ! bending or torsion with
             if (it2 .eq. 1 ) then  ! stretching 
                 pass(ia) = valu/fgrd
             else                   ! bending or torsion
                 pass(ia) = valu/fene
             endif
          endif
        endif
      endif

c     test: conflicting assignments

      icass = 0
      do ip=1,npar
        if(idep(ip).eq.ia) then 
          write(6,67)
          stop
        endif
        if(idep(ia).ne.0) then
          write(6,68)
          stop
        endif
      enddo
 
      goto 50
   60 continue

      if(z_fitlj .and. nass.eq.0) 
     &   write(6,'(10x,"---------- no assigned parameters -------")')
      if(.not. z_fitlj)
     &   write(6,'(10x,"- all nonbonded prms have been assigned -")')

c     prints assigned parameters 

      icass = 0
      do ip=1,npar
        if (iass(ip) .eq. 1) then 
           icass = icass + 1     
           if (iclass(ip) .eq. 1) then 
              valu = pass(ip)*fhes
           elseif (ip .ge. ichrgprm(1) .and. ip .le. ichrgprm(2)) then
              valu = pass(ip)
           elseif (iclass(ip) .eq. 8) then 
              ic1 = iflag(ip,1)
              ic2 = iflag(ip,2)
              it1 = iflag(ic1,0) 
              it2 = iflag(ic2,0) 
              if (it1 .eq. 1 ) then  ! stretching coupled with 
                 if (it2 .eq. 1 ) then  ! stretching 
                     valu = pass(ip)*fhes
                 else                   ! bending or torsion
                     valu = pass(ip)*fgrd
                 endif
              else                   ! bending or torsion with
                 if (it2 .eq. 1 ) then  ! stretching 
                     valu = pass(ip)*fgrd
                 else                   ! bending or torsion
                     valu = pass(ip)*fene
                 endif
              endif
           else
              valu = pass(ip)*fene
           endif        
           write(6,64) icass,ip,pass(ip),valu
        endif  
      enddo
 
   61 format(//10x,'ASSIGN ERROR: too few data',/)
   62 format(//10x,'ASSIGN ERROR: coupling par.s not allowed',/)
   63 format(//10x,'ASSIGN ERROR: index out of range',/)
   64 format(10x,i4,3x,'param',i4,3x,'=  ',f11.5,4x,'input = ',f10.4)
   66 format(//10x,'ASSIGN ERROR: one param is dep. by itself ?',/)
   67 format(//10x,'ASSIGN ERROR: assigned param is a dependent param')
   68 format(//10x,'ASSIGN ERROR: assigned param is a param'
     &           ,' used as ref for depend',/)
 

c     Saves the alpha and beta arrays
 
      do mu=1,npt
        betas(mu)=beta(mu)
        do ip=1,npar
          alphas(mu,ip)=alpha(mu,ip)
        enddo
      enddo

c     modifies the beta array in the case some parameters have been
c     assigned. Some elements of the alpha array are set to zero.

      if(nass.gt.0) then
        do ip=1,npar
          if(iass(ip).gt.0) then
            do mu=1,npt
              beta(mu)=beta(mu)-pass(ip)*alpha(mu,ip)
              alpha(mu,ip)=0.0d0
            enddo
          endif
        enddo
      endif

c     modifies the Alpha array according to the dependencies
c     among parameters.

      do ia=1,npar
        ib=idep(ia)
        if(ib.gt.0) then
          r=rdep(ia)
          do mu=1,npt
            alpha(mu,ib)=alpha(mu,ib)+r*alpha(mu,ia)
            alpha(mu,ia)=0.0d0
          enddo
        endif
      enddo

c     normalizes the sum of the weigths

      iun=5
      zrew=.true.
      Call Cerca_W ('$sumw',iun,card,ier,zrew)
      sumweig=10.d0
      if(ier.eq.0) read(card,*) sumweig
      if(sumweig.lt.1.d0) sumweig=1.0d0
      if(sumweig.gt.1.d5) sumweig=1.0d5
      www=0.0d0
      do i=1,npt
        www=www+weig(i)
      enddo
      ff=sumweig/www
      do i=1,npt
        weig(i)=weig(i)*ff
      enddo

c     analyze and print the redundancy of the functions

      ndepa=0
      do jp=1,npar
       idepa(jp)=0
       T1(jp,jp)=1.0d0
      enddo
      nparfree=npar-nass-ndep
      if(nparfree.gt.2) then
        sumnew=1000.d0
        facred=sumnew/sumweig
        do jp=1,npar
         do ip=1,npar
            aij=0.0d0
            do mu=1,npt
              aij=aij+alpha(mu,ip)*alpha(mu,jp)*weig(mu)*facred
            enddo
            G(ip,jp)=aij
         enddo
        enddo
 
        nred=0
        if (levprt .gt. 0 ) write(6,80) thrdep
        redmax=0.0d0
        thrnew=sumnew*1.d-24
        do jp=2,npar
          T1(jp,jp)=1.0d0
          anormj=G(jp,jp)
          do ip=1,jp-1
            anormi=G(ip,ip)
            T1(ip,jp)=0.0d0
            T1(jp,ip)=0.0d0
            if(anormi*anormj.gt.thrnew) then
              red=G(ip,jp)/dsqrt(anormi*anormj)
              T1(ip,jp)=red
              T1(jp,ip)=red
              if(red.gt.redmax) then
                redmax=red
                ki=ip
                kj=jp
              endif
c             prepares data for print
              if(red.gt.0.95d0) then
                if(levprt.ge.2) write(6,82) ip,jp,red
                if(nred.lt.max1) then
                  nred=nred+1
                  iwra(nred)=jp
                  iwrb(nred)=ip
                  wwab(nred)=red
                endif
              endif
c             sets an automatic dependency between redundant parameters
              if(red.gt.thrdep) then
                if(idepa(jp).eq.0 .and. ifree(jp).eq.0) then
                  ndepa=ndepa+1
                  idepa(jp)=ip
                  rdepa(jp)=1.0d0
                  if(red.lt.0.0d0) rdepa(jp)=-1.0d0
                  if (levprt .gt. 0 ) 
     &             write(6,84) ndepa,jp,ip,rdepa(jp),red
                endif
              endif
            endif
          enddo
        enddo
        if (levprt .gt. 0) 
     &     write(6,86) ki,kj,redmax,labpar(ki),labpar(kj)
        if(nred.gt.0 .and. levprt .gt. 0 ) 
     &      write(6,88) (iwra(i),iwrb(i),wwab(i),i=1,nred)
        if(ndepa.gt.0) write(6,90) ndepa
        title='Redundancy Matrix of the Functions  '
        nab=npar
        if(levprt.ge.2) Call Pr_Rec2 (T1,npar,nab,title,ncol,max1)
      endif

   80 format(/10x,'========= REDUNDANCY OF THE FUNCTIONS ========='
     &      ,/10x,'threshold for automatic dependencies',f9.4)
   82 format(10x,'functions',2i5,'  are very redundant:',f9.4)
   84 format(i4,'-',4x,'par',i5,' is made dependent from par',i4
     &     ,' with const',f7.3,'     overlap',f8.4)
   86 format(10x,'max redundancy',2i4,f9.4,4x,a,3x,a)
   88 format(/10x,'Very high redundancies',/(4(i8,i4,f7.3)))
   90 format(10x,i8,'  automatic dependencies are added')


c     modifies the Alpha array according to the automatic
c     dependencies just added to those given by input file


      if(ndepa.gt.0) then
        do ia=1,npar
          ib=idepa(ia)
          if(ib.gt.0) then
            r=rdepa(ia)
            do mu=1,npt
              alpha(mu,ib)=alpha(mu,ib)+r*alpha(mu,ia)
              alpha(mu,ia)=0.0d0
            enddo
          endif
        enddo
      endif

c    computes the n. of free parameters

      nfree = 0
      do ip=1,npar
        ifree(ip)=1
        if(idep(ip) .ne.0) ifree(ip)=0
        if(idepa(ip).ne.0) ifree(ip)=0
        if(iass(ip) .ne.0) ifree(ip)=0
        nfree=nfree+ifree(ip)
      enddo
      write(6,102) nfree
  102 format(/16x,'Number of free parameters .........',i5)
  103 format( 16x,'skip the system resolution')
      if(nfree.eq.0) then
        write(6,103)
        goto 140
      endif

c     prints the alpha and beta matrices (scaled)

      levp=-1
      if(ndep.gt.0 .or. ndepa.gt.0) levp=levprt
      if(nass.gt.0)                 levp=levprt
      if(levp.ge.3) then
        do mu=1,npt
          do i=1,npar
            T2(mu,i)=alpha(mu,i)*alphafac
          enddo
        enddo
        do mu=1,npt
          T2(mu,npar+1)=beta(mu)*alphafac
        enddo
        nab=npar+1
        write(word,'("  (*",1pd10.2,")")') alphafac
        title='ALPHA , BETA (with depend+assign)'//word(1:15)
        Call Pr_Rec2 (T2,npt,nab,title,ncol,max2)
      endif
 
c     computes the A,B arrays
 
      do jp=1,npar
        do ip=1,npar
          aij=0.0d0
          do mu=1,npt
            aij=aij+alpha(mu,ip)*alpha(mu,jp)*weig(mu)
          enddo
          A(ip,jp)=aij
        enddo
      enddo
      do ip=1,npar
        bi=0.0d0
        do mu=1,npt
          bi=bi+alpha(mu,ip)*beta(mu)*weig(mu)
        enddo
        B(ip,1)=bi
      enddo

c     reads the threshold for null eigenvalues
c     and multipy A,B matrices by a factor
 
      iun=5
      zrew=.true.
      thr0=1.d-10
      Call Cerca_W ('$zero',iun,card,ier,zrew)
      if(ier.eq.0) read(card,*) thr0
      thr0=dabs(thr0)
 
      aijmax=0.0d0
      do jp=1,npar
        do ip=1,jp
          if(dabs(A(ip,jp)).gt.aijmax) aijmax=dabs(A(ip,jp)) 
        enddo
      enddo
      fmult=10.d0/aijmax
      if(fmult.lt.5.0d0) then
        fmult=1.0d0
      else
        aa=1.0d0
        dowhile (aa.lt.fmult)
          aa=aa*10.d0
        enddo
        fmult=aa
      endif
 
      write(6,114) thr0,aijmax,fmult,sumweig
  114 format(/10x,'given threshold for null eigenvalues ....',1d12.3
     &      ,' * Max(eigenvalue)'
     &      ,/10x,'maximum value of the A matrix ...........',d12.3
     &      ,/10x,'both A and B are multiplied by ..........',d12.3
     &      ,/10x,'sum (over points) of the weights ........',d12.3)

      do i=1,npar
        B(i,1)=B(i,1)*fmult
        do j=1,npar
          A(i,j)=A(i,j)*fmult
        enddo
      enddo

c     saves the A matrix (needed for leave-one-out)

      do i=1,npar
        do j=1,npar
          G(i,j)=A(i,j)
        enddo
      enddo
 
c     prints the A and B matrices
 
      do i=1,npar
        A(i,npar+1)=B(i,1)
      enddo
      nab=npar+1
      title='Lin. Syst. Matrix A with vect. B in the last column'
      if(levprt.ge.2) Call Pr_Rec2 (A,npar,nab,title,ncol,max1)


c ---------------------------------------------------------------
c                  SOLVES THE LINEAR SYSTEM
c    (possibly) modified for dependences and assignments
c ---------------------------------------------------------------

      m=1
      Call Lin_Eq_SVD (A,B,X,npar,m,thr0,Sens,max1,ierr)
      if(ierr.ne.0)  then
        write(6,'(//8x,"ABNORMAL TERMINATION OF JOYCE",//)')
        stop
      endif

c     recovers the parameters

  140 continue  !! (previous inst.s are skipped if nfree=0)
      neff = 0
      sums = 0
      do ip=1,npar
        kp=idep(ip)
        if(kp.gt.0) X(ip,1)=X(kp,1)*rdep(ip)
        if(kp.gt.0) Sens(ip)=Sens(kp)
        kp=idepa(ip)
        if(kp.gt.0) X(ip,1)=X(kp,1)*rdepa(ip)
        if(kp.gt.0) Sens(ip)=Sens(kp)
        mp=iass(ip)
        if(mp.gt.0) X(ip,1)=pass(ip)
        if(mp.gt.0) Sens(ip)=0.0d0
        par(ip)=X(ip,1)
        if(kp.eq.0 .and. mp.eq.0) then
          neff=neff+1
          sums=sums+Sens(ip)
        endif
      enddo

c     prints the obtained parameters

      write(6,142) npar,neff
      do ip=1,npar
        labi='        '
        if(idep(ip).gt.0)  write(labi,'("depend",i3)') idep(ip)
        if(idepa(ip).gt.0) write(labi,'("depend",i3)') idepa(ip)
        if(iass(ip).gt.0) write(labi,'("assigned ")')
        ll=Lenw(labpar(ip))
        sip=Sens(ip)
        write(6,144) ip,par(ip),labi,sip    ,labpar(ip)(1:ll)
      enddo
      if (levprt .gt. 0 ) write(6,146) sums
  142 format(//10x,'==== VALUE OF THE BEST LINEAR PARAMETERS ===='
     &       ,/10x,'# param.s',i5,8x,'#. of free param.s',i5
     &       ,/10x,'sensitivity index:    =1 OK,  =0 useless function'
     &       ,/10x,'---------------------------------------------'
     &   ,/10x,'#     value',15x,'Sens Index')
  144 format(10x,i4,f11.6,2x,a,t27,f10.3,3x,a)
  146 format(10x,'Sum of sens indeces',t27,f10.3)

c     prints the global residues and variance (au)

      wtot=0.0d0
      do mu=1,npt
        wtot=wtot+weig(mu)
      enddo
      fx=1.0d5
      if (levprt .gt. 0) write(6,162) fx
      chisq=0.0d0
      do mu=1,npt
        val=0.0d0
        do ip=1,npar
          val=val+par(ip)*alphas(mu,ip)
        enddo
        res=betas(mu)-val
        chisq=chisq+weig(mu)*res*res/wtot
        if (levprt .gt. 0 ) 
     &      write(6,164) mu,fx*betas(mu),fx*val,fx*res,fx*fx*chisq,
     &                   labpt(mu),weig(mu)
      enddo
      stdev=dsqrt(chisq)
      write(6,166) stdev
  162 format(//10x,'RESULTS OF THE LINEAR FITTING (values*',1d10.1,')'
     &    ,/10x,'exact  computed   residue    Chisq')
  164 format(10x,i6,3f10.4,f10.5,3x,a,2x,f8.4)
  166 format(10x,'TOTAL STANDARD DEVIATION',1d10.3,'   atomic units')
 
c     prints the global residues and variance (kJ/mol)
 
      write(6,172)
      chisq=0.0d0
      do mu=1,npt
      fx=fene
      if(index(labpt(mu),'gr').gt.0) fx=fgrd
      if(index(labpt(mu),'he').gt.0) fx=fhes
        val=0.0d0
        do ip=1,npar
          val=val+par(ip)*alphas(mu,ip)
        enddo
        res=(betas(mu)-val)*fx
        chisq=chisq+weig(mu)*res*res/wtot
        write(6,174) mu,fx*betas(mu),fx*val,res,chisq,labpt(mu)
     #              ,weig(mu)
      enddo
      stdev=dsqrt(chisq)
      write(6,176) stdev
  172 format(//10x,'RESULTS OF THE LINEAR FITTING (kJ/mol L=angst)'
     #    ,/10x,'exact  computed   residue    Chisq')
  174 format(10x,i6,3f10.4,f10.5,3x,a,2x,f8.4)
  176 format(10x,'TOTAL STANDARD DEVIATION',1pd10.3,'  kJ/mol L=angst')
 

c     leave-one-out method

      iun=5
      zrew=.true.
      Call Cerca_W ('$leave',iun,card,ier,zrew)
      if(ier.ne.0) goto 240
      Call Setgrat (card)
      Call Count_w (card,nw)
      write(6,224)
      do 220 il=1,nw
        Call SelWrd (card,word,il,lw,ipos)
        read(word,*) ipar
        if(ipar.lt.0 .or. ipar.gt.npar) goto 220
        if(iass(ipar).gt.0) then
          write(6,'("   param",i5,"  has been assigned")') ipar
          goto 220
        endif
c       prepare the new A matrix
        do i=1,npar
          ileave(i)=0
        enddo
        ileave(ipar)=1
        do i=1,npar
          if(idep(i).eq.ipar) ileave(i)=1
        enddo
        do i=1,npar
          do j=1,npar
            A(i,j)=G(i,j)
            if(ileave(i)+ileave(j).gt.0) A(i,j)=0.0d0
          enddo
          if(ileave(i).gt.0) A(i,i)=1.0d0
        enddo
c       solve the linear system
        m=1
        Call Lin_Eq_SVD (A,B,X,npar,m,thr0,Sens,max1,ierr)
        if(ierr.ne.0) stop
c       store the new parameterss
        neff=0
        do ip=1,npar
          kp=idep(ip)
          if(kp.gt.0) X(ip,1)=X(kp,1)*rdep(ip)
          mp=iass(ip)
          if(mp.gt.0) X(ip,1)=pass(ip)
          if(kp.eq.0 .and. mp.eq.0) neff=neff+1
          if(ileave(ip).gt.0) X(ip,1)=0.0d0
        enddo
c       compute stand. dev. and print
        chisq=0.0d0
        do mu=1,npt
          val=0.0d0
          do ip=1,npar
            val=val+X(ip,1)*alpha(mu,ip)
          enddo
          res=beta(mu)-val
          chisq=chisq+weig(mu)*res*res/wtot
        enddo
        sd=dsqrt(chisq)
        ii=0
        do ip=1,npar
          if(ileave(ip).gt.0) then
            ii=ii+1
            kleave(ii)=ip
          endif
        enddo
        write(6,226) sd,sd/stdev,(kleave(m),m=1,ii)
  220 continue
c
  224 format(///7x,'------- L E A V E - O N E - O U T ---------'
     #  ,/5x,'St.Dev   SD/SD0    null param.s')
  226 format(i4,2f12.6,4x,9i4)
  240 continue
 
c     punch parameters in a form suitable to be used in the
c     '$assign' input data (for successive runs)
 
      iun=33
      open(iun,file='assign.dat',form='formatted')
      write(iun,'("$assign")')
      do i=1,npar
          n0=iclass(i)
          valu=par(i)*fene
          if(n0.eq.1) valu=par(i)*fhes
          if (n0 .eq. 8) then  ! couplings
              ic1 = iflag(i,1)
              ic2 = iflag(i,2)
              it1 = iflag(ic1,0) 
              it2 = iflag(ic2,0) 
              if (it1 .eq. 1 ) then  ! stretching coupled with 
                 if (it2 .eq. 1 ) then  ! stretching 
                     valu = par(i)*fhes
                 else                   ! bending or torsion
                     valu = par(i)*fgrd
                 endif
              else                   ! bending or torsion with
                 if (it2 .eq. 1 ) then  ! stretching 
                     valu = par(i)*fgrd
                 else                   ! bending or torsion
                     valu = par(i)*fene
                 endif
              endif
          endif        
          if(n0.eq.5) then   !! intra LJ
             write(iun,'(i3," = ",f14.7,5x,a,5x," non-bonded LJ")') 
     &             i,valu,labpar(i)  !!  nonbonded is assigned by default
             elseif(n0.eq.6) then  !! intra charge-charge
             valu=par(i)
             write(iun,'(i3," = ",f14.7,5x,a,5x," non-bonded Coul")') 
     &             i,valu,labpar(i)  !!  nonbonded is assigned by default
          else   
             write(iun,'(i3," = ",f14.7,5x,a)') i,valu,labpar(i)
          endif
      enddo
      write(iun,'("$end")')
      close (iun)
 
      return
      end
