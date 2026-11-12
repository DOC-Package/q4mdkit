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
      Subroutine Compare_NM (z_plotfreq,Eff,Eok,Cff,Cok,W,
     &                       nvib,ncart,iw,ndim)
c----------------------------------------------------------------------c
 
c  compares frequencies and normal modes of QM and MM calculations
c              ff = from Force Field     ok = from g09
 
      implicit real*8(a-h,o-z)
      dimension Eff(*),Eok(*),Cff(ndim,*),Cok(ndim,*),W(ndim,*),iw(*)
      dimension overlap(ndim)
      character*80 title
      character*25 freqgracefile, torsgracefile
      logical z_plotfreq
      common /printout/ levprt
      common /gracefiles/freqgracefile, torsgracefile
      data har2cm1 /219474.6156d0/    !! 1 hartree = 221474 1/cm
      ilead=1        !! 0=arrange freq acconding to:  0=FF  1=G09
c
      ncol=10
      write(6,15)
   15 format(/9x,'================================================='
     &      ,/9x,'     Compare Norm Modes from QM and FF'
     &      ,/9x,'=================================================')
      do i=1,nvib
        iw(i)=0
      enddo
 
c   overlap between the two normal modes
 
      do i=1,nvib
        do j=1,nvib
        aij=0.0d0
        do k=1,ncart
          aij=aij+Cff(k,i)*Cok(k,j)
        enddo
        W(i,j)=aij
        enddo
      enddo
      title='Overlap between Norm Vibr: <FFoptim)|G09>'
      ni=nvib
      if(levprt.ge.2) Call Pr_Rec2 (W,nvib,ni,title,ncol,ndim)
 
c   leading sequence of frequencies: from FF 
 
      if(ilead.gt.0) goto 20
      write(6,32)
      stdv=0.0d0
      do i=nvib,1,-1
        omax=0.0d0
        do j=1,nvib
          omij=dabs(W(i,j))
          if(iw(j).eq.0 .and. omij.gt.omax) then
            omax=omij
            jmax=j
          endif
        enddo
 
        xx=dabs(Eff(i))
        fi=dsqrt(xx)*har2cm1
        if(Eff(i).lt.0.0d0) fi=-fi
 
        j=jmax
        iw(j)=i
        xx=dabs(Eok(j))
        fj=dsqrt(xx)*har2cm1
        if(Eok(j).lt.0.0d0) fj=-fj
 
        err=fi-fj
        stdv=stdv+err*err
        write(6,34) i,fi,j,fj,omax,err
      enddo
      stdv=dsqrt(stdv/dble(nvib))
      write(6,38) stdv
      return

c   leading sequence of frequencies: from G09
 
   20 write(6,32)
      stdv=0.0d0
      do j=nvib,1,-1
        omax=0.0d0
        do i=1,nvib
          omij=dabs(W(i,j))
          if(iw(i).eq.0 .and. omij.gt.omax) then
            omax=omij
            imax=i
          endif
        enddo
 
        i=imax
        xx=dabs(Eff(i))
        fi=dsqrt(xx)*har2cm1
        if(Eff(i).lt.0.0d0) fi=-fi
        iw(i)=j
 
        xx=dabs(Eok(j))
        fj=dsqrt(xx)*har2cm1
        if(Eok(j).lt.0.0d0) fj=-fj
 
        err=fi-fj
        stdv=stdv+err*err
        write(6,34) i,fi,j,fj,omax,err
        overlap(j) = omax
      enddo
 
      stdv=dsqrt(stdv/dble(nvib))
      write(6,38) stdv

      if (z_plotfreq) then 
       write(6,'(10x,a)') 'Plotting data on grace file'
       open(173,file=freqgracefile,form='formatted',status='unknown')
       write(173,'(a)') '# Joyce made Grace project file'
       write(173,'(a)') '# '
c-------------------

       write(173,'(a)') '@version 50125'
       write(173,'(a)') '@page size 792, 612'
       write(173,'(a)') '@page scroll 5%'
       write(173,'(a)') '@page inout 5%'
       write(173,'(a)') '@link page off'
       write(173,'(a)') '@map font 0 to "Times-Roman", "Times-Roman"'
       write(173,'(a)') '@map font 1 to "Times-Italic", "Times-Italic"'
       write(173,'(a)') '@map font 2 to "Times-Bold", "Times-Bold"'
       write(173,'(a)') 
     &       '@map font 3 to "Times-BoldItalic", "Times-BoldItalic"'
       write(173,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
       write(173,'(a)') 
     &       '@map font 5 to "Helvetica-Oblique", "Helvetica-Oblique"'
       write(173,'(a)') 
     &       '@map font 6 to "Helvetica-Bold", "Helvetica-Bold"'
       write(173,'(a)') 
     & '@map font 7 to "Helvetica-BoldOblique", "Helvetica-BoldOblique"'
       write(173,'(a)') '@map font 8 to "Courier", "Courier"'
       write(173,'(a)') 
     &       '@map font 9 to "Courier-Oblique", "Courier-Oblique"'
       write(173,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
       write(173,'(a)') 
     & '@map font 11 to "Courier-BoldOblique", "Courier-BoldOblique"'
       write(173,'(a)') '@map font 12 to "Symbol", "Symbol"'
       write(173,'(a)') '@map font 13 to "ZapfDingbats", "ZapfDingbats"'
       write(173,'(a)') '@map color 0 to (255, 255, 255), "white"'
       write(173,'(a)') '@map color 1 to (0, 0, 0), "black"'
       write(173,'(a)') '@map color 2 to (255, 0, 0), "red"'
       write(173,'(a)') '@map color 3 to (0, 255, 0), "green"'
       write(173,'(a)') '@map color 4 to (0, 0, 255), "blue"'
       write(173,'(a)') '@map color 5 to (255, 255, 0), "yellow"'
       write(173,'(a)') '@map color 6 to (188, 143, 143), "brown"'
       write(173,'(a)') '@map color 7 to (220, 220, 220), "grey"'
       write(173,'(a)') '@map color 8 to (148, 0, 211), "violet"'
       write(173,'(a)') '@map color 9 to (0, 255, 255), "cyan"'
       write(173,'(a)') '@map color 10 to (255, 0, 255), "magenta"'
       write(173,'(a)') '@map color 11 to (255, 165, 0), "orange"'
       write(173,'(a)') '@map color 12 to (114, 33, 188), "indigo"'
       write(173,'(a)') '@map color 13 to (103, 7, 72), "maroon"'
       write(173,'(a)') '@map color 14 to (64, 224, 208), "turquoise"'
       write(173,'(a)') '@map color 15 to (0, 139, 0), "green4"'
       write(173,'(a)') '@reference date 0'
       write(173,'(a)') '@date wrap off'
       write(173,'(a)') '@date wrap year 1950'
       write(173,'(a)') '@default linewidth 1.0'
       write(173,'(a)') '@default linestyle 1'
       write(173,'(a)') '@default color 1'
       write(173,'(a)') '@default pattern 1'
       write(173,'(a)') '@default font 0'
       write(173,'(a)') '@default char size 1.000000'
       write(173,'(a)') '@default symbol size 1.000000'
       write(173,'(a)') '@default sformat "%.8g"'
       write(173,'(a)') '@background color 0'
       write(173,'(a)') '@page background fill on'
       write(173,'(a)') '@timestamp off'
       write(173,'(a)') '@timestamp 0.03, 0.03'
       write(173,'(a)') '@timestamp color 1'
       write(173,'(a)') '@timestamp rot 0'
       write(173,'(a)') '@timestamp font 0'
       write(173,'(a)') '@timestamp char size 1.000000'
       write(173,'(a)') '@timestamp def "          "'
       write(173,'(a)') '@r0 off'
       write(173,'(a)') '@link r0 to g0'
       write(173,'(a)') '@r0 type above'
       write(173,'(a)') '@r0 linestyle 1'
       write(173,'(a)') '@r0 linewidth 1.0'
       write(173,'(a)') '@r0 color 1'
       write(173,'(a)') '@r0 line 0, 0, 0, 0'
       write(173,'(a)') '@r1 off'
       write(173,'(a)') '@link r1 to g0'
       write(173,'(a)') '@r1 type above'
       write(173,'(a)') '@r1 linestyle 1'
       write(173,'(a)') '@r1 linewidth 1.0'
       write(173,'(a)') '@r1 color 1'
       write(173,'(a)') '@r1 line 0, 0, 0, 0'
       write(173,'(a)') '@r2 off'
       write(173,'(a)') '@link r2 to g0'
       write(173,'(a)') '@r2 type above'
       write(173,'(a)') '@r2 linestyle 1'
       write(173,'(a)') '@r2 linewidth 1.0'
       write(173,'(a)') '@r2 color 1'
       write(173,'(a)') '@r2 line 0, 0, 0, 0'
       write(173,'(a)') '@r3 off'
       write(173,'(a)') '@link r3 to g0'
       write(173,'(a)') '@r3 type above'
       write(173,'(a)') '@r3 linestyle 1'
       write(173,'(a)') '@r3 linewidth 1.0'
       write(173,'(a)') '@r3 color 1'
       write(173,'(a)') '@r3 line 0, 0, 0, 0'
       write(173,'(a)') '@r4 off'
       write(173,'(a)') '@link r4 to g0'
       write(173,'(a)') '@r4 type above'
       write(173,'(a)') '@r4 linestyle 1'
       write(173,'(a)') '@r4 linewidth 1.0'
       write(173,'(a)') '@r4 color 1'
       write(173,'(a)') '@r4 line 0, 0, 0, 0'
       write(173,'(a)') '@g0 on'
       write(173,'(a)') '@g0 hidden false'
       write(173,'(a)') '@g0 type XY'
       write(173,'(a)') '@g0 stacked false'
       write(173,'(a)') '@g0 bar hgap 0.000000'
       write(173,'(a)') '@g0 fixedpoint off'
       write(173,'(a)') '@g0 fixedpoint type 0'
       write(173,'(a)') '@g0 fixedpoint xy 0.000000, 0.000000'
       write(173,'(a)') '@g0 fixedpoint format general general'
       write(173,'(a)') '@g0 fixedpoint prec 6, 6'
       write(173,'(a)') '@with g0'
       write(173,'(a17,1x,i4,1x,a3)')'@    world 0, 0,',Nvib, ', 1'
       write(173,'(a)') '@    stack world 0, 0, 0, 0'
       write(173,'(a)') '@    znorm 1'
       write(173,'(a)') 
     &      '@    view 0.150000, 0.670000, 1.000000, 0.850000'
       write(173,'(a)') '@    title ""'
       write(173,'(a)') '@    title font 0'
       write(173,'(a)') '@    title size 1.500000'
       write(173,'(a)') '@    title color 1'
       write(173,'(a)') '@    subtitle ""'
       write(173,'(a)') '@    subtitle font 0'
       write(173,'(a)') '@    subtitle size 1.000000'
       write(173,'(a)') '@    subtitle color 1'
       write(173,'(a)') '@    xaxes scale Normal'
       write(173,'(a)') '@    yaxes scale Normal'
       write(173,'(a)') '@    xaxes invert off'
       write(173,'(a)') '@    yaxes invert off'
       write(173,'(a)') '@    xaxis  on'
       write(173,'(a)') '@    xaxis  type zero false'
       write(173,'(a)') '@    xaxis  offset 0.000000 , 0.000000'
       write(173,'(a)') '@    xaxis  bar on'
       write(173,'(a)') '@    xaxis  bar color 1'
       write(173,'(a)') '@    xaxis  bar linestyle 1'
       write(173,'(a)') '@    xaxis  bar linewidth 1.0'
       write(173,'(a)') '@    xaxis  label "Q\sK\N"'
       write(173,'(a)') '@    xaxis  label layout para'
       write(173,'(a)') '@    xaxis  label place auto'
       write(173,'(a)') '@    xaxis  label char size 1.000000'
       write(173,'(a)') '@    xaxis  label font 6'
       write(173,'(a)') '@    xaxis  label color 1'
       write(173,'(a)') '@    xaxis  label place opposite'
       write(173,'(a)') '@    xaxis  tick on'
       write(173,'(a)') '@    xaxis  tick major 25'
       write(173,'(a)') '@    xaxis  tick minor ticks 4'
       write(173,'(a)') '@    xaxis  tick default 6'
       write(173,'(a)') '@    xaxis  tick place rounded true'
       write(173,'(a)') '@    xaxis  tick in'
       write(173,'(a)') '@    xaxis  tick major size 1.000000'
       write(173,'(a)') '@    xaxis  tick major color 1'
       write(173,'(a)') '@    xaxis  tick major linewidth 1.0'
       write(173,'(a)') '@    xaxis  tick major linestyle 1'
       write(173,'(a)') '@    xaxis  tick major grid off'
       write(173,'(a)') '@    xaxis  tick minor color 1'
       write(173,'(a)') '@    xaxis  tick minor linewidth 1.0'
       write(173,'(a)') '@    xaxis  tick minor linestyle 1'
       write(173,'(a)') '@    xaxis  tick minor grid off'
       write(173,'(a)') '@    xaxis  tick minor size 0.500000'
       write(173,'(a)') '@    xaxis  ticklabel on'
       write(173,'(a)') '@    xaxis  ticklabel format general'
       write(173,'(a)') '@    xaxis  ticklabel prec 5'
       write(173,'(a)') '@    xaxis  ticklabel formula ""'
       write(173,'(a)') '@    xaxis  ticklabel append ""'
       write(173,'(a)') '@    xaxis  ticklabel prepend ""'
       write(173,'(a)') '@    xaxis  ticklabel angle 0'
       write(173,'(a)') '@    xaxis  ticklabel skip 0'
       write(173,'(a)') '@    xaxis  ticklabel stagger 0'
       write(173,'(a)') '@    xaxis  ticklabel place opposite'
       write(173,'(a)') '@    xaxis  ticklabel offset auto'
       write(173,'(a)') 
     &   '@    xaxis  ticklabel offset 0.000000 , 0.010000'
       write(173,'(a)') '@    xaxis  ticklabel start type spec'
       write(173,'(a)') '@    xaxis  ticklabel start 1.000000'
       write(173,'(a)') '@    xaxis  ticklabel stop type auto'
       write(173,'(a)') '@    xaxis  ticklabel stop 0.000000'
       write(173,'(a)') '@    xaxis  ticklabel char size 1.000000'
       write(173,'(a)') '@    xaxis  ticklabel font 4'
       write(173,'(a)') '@    xaxis  ticklabel color 1'
       write(173,'(a)') '@    xaxis  tick place both'
       write(173,'(a)') '@    xaxis  tick spec type none'
       write(173,'(a)') '@    yaxis  on'
       write(173,'(a)') '@    yaxis  type zero false'
       write(173,'(a)') '@    yaxis  offset 0.000000 , 0.000000'
       write(173,'(a)') '@    yaxis  bar on'
       write(173,'(a)') '@    yaxis  bar color 1'
       write(173,'(a)') '@    yaxis  bar linestyle 1'
       write(173,'(a)') '@    yaxis  bar linewidth 1.0'
       write(173,'(a)') '@    yaxis  label "<Q\SQM\N|Q\SQMD-FF\N>"'
       write(173,'(a)') '@    yaxis  label layout para'
       write(173,'(a)') '@    yaxis  label place auto'
       write(173,'(a)') '@    yaxis  label char size 1.000000'
       write(173,'(a)') '@    yaxis  label font 6'
       write(173,'(a)') '@    yaxis  label color 1'
       write(173,'(a)') '@    yaxis  label place normal'
       write(173,'(a)') '@    yaxis  tick on'
       write(173,'(a)') '@    yaxis  tick major 0.5'
       write(173,'(a)') '@    yaxis  tick minor ticks 1'
       write(173,'(a)') '@    yaxis  tick default 6'
       write(173,'(a)') '@    yaxis  tick place rounded true'
       write(173,'(a)') '@    yaxis  tick in'
       write(173,'(a)') '@    yaxis  tick major size 1.000000'
       write(173,'(a)') '@    yaxis  tick major color 1'
       write(173,'(a)') '@    yaxis  tick major linewidth 1.0'
       write(173,'(a)') '@    yaxis  tick major linestyle 1'
       write(173,'(a)') '@    yaxis  tick major grid off'
       write(173,'(a)') '@    yaxis  tick minor color 1'
       write(173,'(a)') '@    yaxis  tick minor linewidth 1.0'
       write(173,'(a)') '@    yaxis  tick minor linestyle 1'
       write(173,'(a)') '@    yaxis  tick minor grid off'
       write(173,'(a)') '@    yaxis  tick minor size 0.500000'
       write(173,'(a)') '@    yaxis  ticklabel on'
       write(173,'(a)') '@    yaxis  ticklabel format general'
       write(173,'(a)') '@    yaxis  ticklabel prec 5'
       write(173,'(a)') '@    yaxis  ticklabel formula ""'
       write(173,'(a)') '@    yaxis  ticklabel append ""'
       write(173,'(a)') '@    yaxis  ticklabel prepend ""'
       write(173,'(a)') '@    yaxis  ticklabel angle 0'
       write(173,'(a)') '@    yaxis  ticklabel skip 0'
       write(173,'(a)') '@    yaxis  ticklabel stagger 0'
       write(173,'(a)') '@    yaxis  ticklabel place normal'
       write(173,'(a)') '@    yaxis  ticklabel offset auto'
       write(173,'(a)') 
     &    '@    yaxis  ticklabel offset 0.000000 , 0.010000'
       write(173,'(a)') '@    yaxis  ticklabel start type auto'
       write(173,'(a)') '@    yaxis  ticklabel start 0.000000'
       write(173,'(a)') '@    yaxis  ticklabel stop type auto'
       write(173,'(a)') '@    yaxis  ticklabel stop 0.000000'
       write(173,'(a)') '@    yaxis  ticklabel char size 1.000000'
       write(173,'(a)') '@    yaxis  ticklabel font 4'
       write(173,'(a)') '@    yaxis  ticklabel color 1'
       write(173,'(a)') '@    yaxis  tick place normal'
       write(173,'(a)') '@    yaxis  tick spec type none'
       write(173,'(a)') '@    altxaxis  off'
       write(173,'(a)') '@    altyaxis  off'
       write(173,'(a)') '@    legend on'
       write(173,'(a)') '@    legend loctype view'
       write(173,'(a)') '@    legend 0.85, 0.8'
       write(173,'(a)') '@    legend box color 1'
       write(173,'(a)') '@    legend box pattern 1'
       write(173,'(a)') '@    legend box linewidth 1.0'
       write(173,'(a)') '@    legend box linestyle 1'
       write(173,'(a)') '@    legend box fill color 0'
       write(173,'(a)') '@    legend box fill pattern 1'
       write(173,'(a)') '@    legend font 0'
       write(173,'(a)') '@    legend char size 1.000000'
       write(173,'(a)') '@    legend color 1'
       write(173,'(a)') '@    legend length 4'
       write(173,'(a)') '@    legend vgap 1'
       write(173,'(a)') '@    legend hgap 1'
       write(173,'(a)') '@    legend invert false'
       write(173,'(a)') '@    frame type 0'
       write(173,'(a)') '@    frame linestyle 1'
       write(173,'(a)') '@    frame linewidth 1.0'
       write(173,'(a)') '@    frame color 1'
       write(173,'(a)') '@    frame pattern 1'
       write(173,'(a)') '@    frame background color 0'
       write(173,'(a)') '@    frame background pattern 0'
       write(173,'(a)') '@    s0 hidden false'
       write(173,'(a)') '@    s0 type bar'
       write(173,'(a)') '@    s0 symbol 0'
       write(173,'(a)') '@    s0 symbol size 0.300000'
       write(173,'(a)') '@    s0 symbol color 15'
       write(173,'(a)') '@    s0 symbol pattern 1'
       write(173,'(a)') '@    s0 symbol fill color 15'
       write(173,'(a)') '@    s0 symbol fill pattern 1'
       write(173,'(a)') '@    s0 symbol linewidth 3.0'
       write(173,'(a)') '@    s0 symbol linestyle 1'
       write(173,'(a)') '@    s0 symbol char 65'
       write(173,'(a)') '@    s0 symbol char font 0'
       write(173,'(a)') '@    s0 symbol skip 0'
       write(173,'(a)') '@    s0 line type 0'
       write(173,'(a)') '@    s0 line linestyle 1'
       write(173,'(a)') '@    s0 line linewidth 1.0'
       write(173,'(a)') '@    s0 line color 1'
       write(173,'(a)') '@    s0 line pattern 1'
       write(173,'(a)') '@    s0 baseline type 0'
       write(173,'(a)') '@    s0 baseline off'
       write(173,'(a)') '@    s0 dropline off'
       write(173,'(a)') '@    s0 fill type 0'
       write(173,'(a)') '@    s0 fill rule 0'
       write(173,'(a)') '@    s0 fill color 1'
       write(173,'(a)') '@    s0 fill pattern 1'
       write(173,'(a)') '@    s0 avalue off'
       write(173,'(a)') '@    s0 avalue type 2'
       write(173,'(a)') '@    s0 avalue char size 1.000000'
       write(173,'(a)') '@    s0 avalue font 0'
       write(173,'(a)') '@    s0 avalue color 1'
       write(173,'(a)') '@    s0 avalue rot 0'
       write(173,'(a)') '@    s0 avalue format general'
       write(173,'(a)') '@    s0 avalue prec 3'
       write(173,'(a)') '@    s0 avalue prepend ""'
       write(173,'(a)') '@    s0 avalue append ""'
       write(173,'(a)') '@    s0 avalue offset 0.000000 , 0.000000'
       write(173,'(a)') '@    s0 errorbar on'
       write(173,'(a)') '@    s0 errorbar place both'
       write(173,'(a)') '@    s0 errorbar color 15'
       write(173,'(a)') '@    s0 errorbar pattern 1'
       write(173,'(a)') '@    s0 errorbar size 1.000000'
       write(173,'(a)') '@    s0 errorbar linewidth 1.0'
       write(173,'(a)') '@    s0 errorbar linestyle 1'
       write(173,'(a)') '@    s0 errorbar riser linewidth 1.0'
       write(173,'(a)') '@    s0 errorbar riser linestyle 1'
       write(173,'(a)') '@    s0 errorbar riser clip off'
       write(173,'(a)') '@    s0 errorbar riser clip length 0.100000'
       write(173,'(a)') '@    s0 comment "Cols 3:5"'
       write(173,'(a)') '@    s0 legend  ""'
       write(173,'(a)') '@g1 on'
       write(173,'(a)') '@g1 hidden false'
       write(173,'(a)') '@g1 type XY'
       write(173,'(a)') '@g1 stacked false'
       write(173,'(a)') '@g1 bar hgap 0.000000'
       write(173,'(a)') '@g1 fixedpoint off'
       write(173,'(a)') '@g1 fixedpoint type 0'
       write(173,'(a)') '@g1 fixedpoint xy 0.000000, 0.000000'
       write(173,'(a)') '@g1 fixedpoint format general general'
       write(173,'(a)') '@g1 fixedpoint prec 6, 6'
       write(173,'(a)') '@with g1'
       write(173,'(a)') '@    world 0, 0, 4000, 4000'
       write(173,'(a)') '@    stack world 0, 0, 0, 0'
       write(173,'(a)') '@    znorm 1'
       write(173,'(a)') 
     &  '@    view 0.150000, 0.150000, 1.000000, 0.650000'
       write(173,'(a)') '@    title ""'
       write(173,'(a)') '@    title font 0'
       write(173,'(a)') '@    title size 1.500000'
       write(173,'(a)') '@    title color 1'
       write(173,'(a)') '@    subtitle ""'
       write(173,'(a)') '@    subtitle font 0'
       write(173,'(a)') '@    subtitle size 1.000000'
       write(173,'(a)') '@    subtitle color 1'
       write(173,'(a)') '@    xaxes scale Normal'
       write(173,'(a)') '@    yaxes scale Normal'
       write(173,'(a)') '@    xaxes invert off'
       write(173,'(a)') '@    yaxes invert off'
       write(173,'(a)') '@    xaxis  on'
       write(173,'(a)') '@    xaxis  type zero false'
       write(173,'(a)') '@    xaxis  offset 0.000000 , 0.000000'
       write(173,'(a)') '@    xaxis  bar on'
       write(173,'(a)') '@    xaxis  bar color 1'
       write(173,'(a)') '@    xaxis  bar linestyle 1'
       write(173,'(a)') '@    xaxis  bar linewidth 1.0'
       write(173,'(a)') '@    xaxis  label "\xn\f{}\SQM\N (cm\S-1\N)"'
       write(173,'(a)') '@    xaxis  label layout para'
       write(173,'(a)') '@    xaxis  label place auto'
       write(173,'(a)') '@    xaxis  label char size 1.000000'
       write(173,'(a)') '@    xaxis  label font 6'
       write(173,'(a)') '@    xaxis  label color 1'
       write(173,'(a)') '@    xaxis  label place normal'
       write(173,'(a)') '@    xaxis  tick on'
       write(173,'(a)') '@    xaxis  tick major 500'
       write(173,'(a)') '@    xaxis  tick minor ticks 1'
       write(173,'(a)') '@    xaxis  tick default 6'
       write(173,'(a)') '@    xaxis  tick place rounded true'
       write(173,'(a)') '@    xaxis  tick in'
       write(173,'(a)') '@    xaxis  tick major size 1.000000'
       write(173,'(a)') '@    xaxis  tick major color 1'
       write(173,'(a)') '@    xaxis  tick major linewidth 1.0'
       write(173,'(a)') '@    xaxis  tick major linestyle 1'
       write(173,'(a)') '@    xaxis  tick major grid off'
       write(173,'(a)') '@    xaxis  tick minor color 1'
       write(173,'(a)') '@    xaxis  tick minor linewidth 1.0'
       write(173,'(a)') '@    xaxis  tick minor linestyle 1'
       write(173,'(a)') '@    xaxis  tick minor grid off'
       write(173,'(a)') '@    xaxis  tick minor size 0.500000'
       write(173,'(a)') '@    xaxis  ticklabel on'
       write(173,'(a)') '@    xaxis  ticklabel format general'
       write(173,'(a)') '@    xaxis  ticklabel prec 5'
       write(173,'(a)') '@    xaxis  ticklabel formula ""'
       write(173,'(a)') '@    xaxis  ticklabel append ""'
       write(173,'(a)') '@    xaxis  ticklabel prepend ""'
       write(173,'(a)') '@    xaxis  ticklabel angle 0'
       write(173,'(a)') '@    xaxis  ticklabel skip 0'
       write(173,'(a)') '@    xaxis  ticklabel stagger 0'
       write(173,'(a)') '@    xaxis  ticklabel place normal'
       write(173,'(a)') '@    xaxis  ticklabel offset auto'
       write(173,'(a)') 
     &    '@    xaxis  ticklabel offset 0.000000 , 0.010000'
       write(173,'(a)') '@    xaxis  ticklabel start type auto'
       write(173,'(a)') '@    xaxis  ticklabel start 0.000000'
       write(173,'(a)') '@    xaxis  ticklabel stop type auto'
       write(173,'(a)') '@    xaxis  ticklabel stop 0.000000'
       write(173,'(a)') '@    xaxis  ticklabel char size 1.000000'
       write(173,'(a)') '@    xaxis  ticklabel font 4'
       write(173,'(a)') '@    xaxis  ticklabel color 1'
       write(173,'(a)') '@    xaxis  tick place both'
       write(173,'(a)') '@    xaxis  tick spec type none'
       write(173,'(a)') '@    yaxis  on'
       write(173,'(a)') '@    yaxis  type zero false'
       write(173,'(a)') '@    yaxis  offset 0.000000 , 0.000000'
       write(173,'(a)') '@    yaxis  bar on'
       write(173,'(a)') '@    yaxis  bar color 1'
       write(173,'(a)') '@    yaxis  bar linestyle 1'
       write(173,'(a)') '@    yaxis  bar linewidth 1.0'
       write(173,'(a)') 
     &       '@    yaxis  label "\xn\f{}\SQMD-FF\N (cm\S-1\N)"'
       write(173,'(a)') '@    yaxis  label layout para'
       write(173,'(a)') '@    yaxis  label place auto'
       write(173,'(a)') '@    yaxis  label char size 1.000000'
       write(173,'(a)') '@    yaxis  label font 6'
       write(173,'(a)') '@    yaxis  label color 1'
       write(173,'(a)') '@    yaxis  label place normal'
       write(173,'(a)') '@    yaxis  tick on'
       write(173,'(a)') '@    yaxis  tick major 1000'
       write(173,'(a)') '@    yaxis  tick minor ticks 1'
       write(173,'(a)') '@    yaxis  tick default 6'
       write(173,'(a)') '@    yaxis  tick place rounded true'
       write(173,'(a)') '@    yaxis  tick in'
       write(173,'(a)') '@    yaxis  tick major size 1.000000'
       write(173,'(a)') '@    yaxis  tick major color 1'
       write(173,'(a)') '@    yaxis  tick major linewidth 1.0'
       write(173,'(a)') '@    yaxis  tick major linestyle 1'
       write(173,'(a)') '@    yaxis  tick major grid off'
       write(173,'(a)') '@    yaxis  tick minor color 1'
       write(173,'(a)') '@    yaxis  tick minor linewidth 1.0'
       write(173,'(a)') '@    yaxis  tick minor linestyle 1'
       write(173,'(a)') '@    yaxis  tick minor grid off'
       write(173,'(a)') '@    yaxis  tick minor size 0.500000'
       write(173,'(a)') '@    yaxis  ticklabel on'
       write(173,'(a)') '@    yaxis  ticklabel format general'
       write(173,'(a)') '@    yaxis  ticklabel prec 5'
       write(173,'(a)') '@    yaxis  ticklabel formula ""'
       write(173,'(a)') '@    yaxis  ticklabel append ""'
       write(173,'(a)') '@    yaxis  ticklabel prepend ""'
       write(173,'(a)') '@    yaxis  ticklabel angle 0'
       write(173,'(a)') '@    yaxis  ticklabel skip 0'
       write(173,'(a)') '@    yaxis  ticklabel stagger 0'
       write(173,'(a)') '@    yaxis  ticklabel place normal'
       write(173,'(a)') '@    yaxis  ticklabel offset auto'
       write(173,'(a)') 
     &       '@    yaxis  ticklabel offset 0.000000 , 0.010000'
       write(173,'(a)') '@    yaxis  ticklabel start type auto'
       write(173,'(a)') '@    yaxis  ticklabel start 0.000000'
       write(173,'(a)') '@    yaxis  ticklabel stop type auto'
       write(173,'(a)') '@    yaxis  ticklabel stop 0.000000'
       write(173,'(a)') '@    yaxis  ticklabel char size 1.000000'
       write(173,'(a)') '@    yaxis  ticklabel font 4'
       write(173,'(a)') '@    yaxis  ticklabel color 1'
       write(173,'(a)') '@    yaxis  tick place both'
       write(173,'(a)') '@    yaxis  tick spec type none'
       write(173,'(a)') '@    altxaxis  off'
       write(173,'(a)') '@    altyaxis  off'
       write(173,'(a)') '@    legend on'
       write(173,'(a)') '@    legend loctype view'
       write(173,'(a)') '@    legend 0.5, 0.8'
       write(173,'(a)') '@    legend box color 1'
       write(173,'(a)') '@    legend box pattern 1'
       write(173,'(a)') '@    legend box linewidth 1.0'
       write(173,'(a)') '@    legend box linestyle 1'
       write(173,'(a)') '@    legend box fill color 0'
       write(173,'(a)') '@    legend box fill pattern 1'
       write(173,'(a)') '@    legend font 0'
       write(173,'(a)') '@    legend char size 1.000000'
       write(173,'(a)') '@    legend color 1'
       write(173,'(a)') '@    legend length 4'
       write(173,'(a)') '@    legend vgap 1'
       write(173,'(a)') '@    legend hgap 1'
       write(173,'(a)') '@    legend invert false'
       write(173,'(a)') '@    frame type 0'
       write(173,'(a)') '@    frame linestyle 1'
       write(173,'(a)') '@    frame linewidth 1.0'
       write(173,'(a)') '@    frame color 1'
       write(173,'(a)') '@    frame pattern 1'
       write(173,'(a)') '@    frame background color 0'
       write(173,'(a)') '@    frame background pattern 0'
       write(173,'(a)') '@    s0 hidden false'
       write(173,'(a)') '@    s0 type xy'
       write(173,'(a)') '@    s0 symbol 0'
       write(173,'(a)') '@    s0 symbol size 0.400000'
       write(173,'(a)') '@    s0 symbol color 1'
       write(173,'(a)') '@    s0 symbol pattern 1'
       write(173,'(a)') '@    s0 symbol fill color 1'
       write(173,'(a)') '@    s0 symbol fill pattern 1'
       write(173,'(a)') '@    s0 symbol linewidth 1.0'
       write(173,'(a)') '@    s0 symbol linestyle 1'
       write(173,'(a)') '@    s0 symbol char 65'
       write(173,'(a)') '@    s0 symbol char font 0'
       write(173,'(a)') '@    s0 symbol skip 0'
       write(173,'(a)') '@    s0 line type 1'
       write(173,'(a)') '@    s0 line linestyle 3'
       write(173,'(a)') '@    s0 line linewidth 1.0'
       write(173,'(a)') '@    s0 line color 1'
       write(173,'(a)') '@    s0 line pattern 1'
       write(173,'(a)') '@    s0 baseline type 0'
       write(173,'(a)') '@    s0 baseline off'
       write(173,'(a)') '@    s0 dropline off'
       write(173,'(a)') '@    s0 fill type 0'
       write(173,'(a)') '@    s0 fill rule 0'
       write(173,'(a)') '@    s0 fill color 1'
       write(173,'(a)') '@    s0 fill pattern 1'
       write(173,'(a)') '@    s0 avalue off'
       write(173,'(a)') '@    s0 avalue type 2'
       write(173,'(a)') '@    s0 avalue char size 1.000000'
       write(173,'(a)') '@    s0 avalue font 0'
       write(173,'(a)') '@    s0 avalue color 1'
       write(173,'(a)') '@    s0 avalue rot 0'
       write(173,'(a)') '@    s0 avalue format general'
       write(173,'(a)') '@    s0 avalue prec 3'
       write(173,'(a)') '@    s0 avalue prepend ""'
       write(173,'(a)') '@    s0 avalue append ""'
       write(173,'(a)') '@    s0 avalue offset 0.000000 , 0.000000'
       write(173,'(a)') '@    s0 errorbar on'
       write(173,'(a)') '@    s0 errorbar place both'
       write(173,'(a)') '@    s0 errorbar color 1'
       write(173,'(a)') '@    s0 errorbar pattern 1'
       write(173,'(a)') '@    s0 errorbar size 1.000000'
       write(173,'(a)') '@    s0 errorbar linewidth 1.0'
       write(173,'(a)') '@    s0 errorbar linestyle 1'
       write(173,'(a)') '@    s0 errorbar riser linewidth 1.0'
       write(173,'(a)') '@    s0 errorbar riser linestyle 1'
       write(173,'(a)') '@    s0 errorbar riser clip off'
       write(173,'(a)') '@    s0 errorbar riser clip length 0.100000'
       write(173,'(a)') '@    s0 comment "correlation line"'
       write(173,'(a)') '@    s0 legend  ""'
       write(173,'(a)') '@    s1 hidden false'
       write(173,'(a)') '@    s1 type xy'
       write(173,'(a)') '@    s1 symbol 1'
       write(173,'(a)') '@    s1 symbol size 0.500000'
       write(173,'(a)') '@    s1 symbol color 15'
       write(173,'(a)') '@    s1 symbol pattern 1'
       write(173,'(a)') '@    s1 symbol fill color 15'
       write(173,'(a)') '@    s1 symbol fill pattern 1'
       write(173,'(a)') '@    s1 symbol linewidth 1.0'
       write(173,'(a)') '@    s1 symbol linestyle 1'
       write(173,'(a)') '@    s1 symbol char 65'
       write(173,'(a)') '@    s1 symbol char font 0'
       write(173,'(a)') '@    s1 symbol skip 0'
       write(173,'(a)') '@    s1 line type 0'
       write(173,'(a)') '@    s1 line linestyle 1'
       write(173,'(a)') '@    s1 line linewidth 1.0'
       write(173,'(a)') '@    s1 line color 15'
       write(173,'(a)') '@    s1 line pattern 1'
       write(173,'(a)') '@    s1 baseline type 0'
       write(173,'(a)') '@    s1 baseline off'
       write(173,'(a)') '@    s1 dropline off'
       write(173,'(a)') '@    s1 fill type 0'
       write(173,'(a)') '@    s1 fill rule 0'
       write(173,'(a)') '@    s1 fill color 1'
       write(173,'(a)') '@    s1 fill pattern 1'
       write(173,'(a)') '@    s1 avalue off'
       write(173,'(a)') '@    s1 avalue type 2'
       write(173,'(a)') '@    s1 avalue char size 1.000000'
       write(173,'(a)') '@    s1 avalue font 0'
       write(173,'(a)') '@    s1 avalue color 1'
       write(173,'(a)') '@    s1 avalue rot 0'
       write(173,'(a)') '@    s1 avalue format general'
       write(173,'(a)') '@    s1 avalue prec 3'
       write(173,'(a)') '@    s1 avalue prepend ""'
       write(173,'(a)') '@    s1 avalue append ""'
       write(173,'(a)') '@    s1 avalue offset 0.000000 , 0.000000'
       write(173,'(a)') '@    s1 errorbar on'
       write(173,'(a)') '@    s1 errorbar place both'
       write(173,'(a)') '@    s1 errorbar color 15'
       write(173,'(a)') '@    s1 errorbar pattern 1'
       write(173,'(a)') '@    s1 errorbar size 1.000000'
       write(173,'(a)') '@    s1 errorbar linewidth 1.0'
       write(173,'(a)') '@    s1 errorbar linestyle 1'
       write(173,'(a)') '@    s1 errorbar riser linewidth 1.0'
       write(173,'(a)') '@    s1 errorbar riser linestyle 1'
       write(173,'(a)') '@    s1 errorbar riser clip off'
       write(173,'(a)') '@    s1 errorbar riser clip length 0.100000'
       write(173,'(a)') '@    s1 comment "Cols 3:4"'
       write(173,'(a)') '@    s1 legend  ""'
       write(173,'(a)') '@target G0.S0'
       write(173,'(a)') '@type bar'
       do ivib = 1,Nvib
         write(173,'(i4,4x,f10.3)') ivib, overlap(ivib)
       enddo
       write(173,'(a)') '&'
       write(173,'(a)') '@target G1.S0'
       write(173,'(a)') '@type xy'
       write(173,'(a)') '0       0'
       write(173,'(a)') '500   500'
       write(173,'(a)') '1000 1000'
       write(173,'(a)') '2000 2000'
       write(173,'(a)') '3000 3000'
       write(173,'(a)') '4000 4000'
       write(173,'(a)') '&'
       write(173,'(a)') '@target G1.S1'
       write(173,'(a)') '@type xy'
       do ivib = 1,Nvib
        xx=dabs(Eok(ivib))
        fi=dsqrt(xx)*har2cm1
        if(Eok(ivib).lt.0.0d0) fi=-fi
        xx=dabs(Eff(ivib))
        fj=dsqrt(xx)*har2cm1
        if(Eff(ivib).lt.0.0d0) fj=-fj
         write(173,'(2(f12.1,3x))') fi, fj
       enddo
       write(173,'(a)') '&'
       close(173)
      endif

      return
 
   32 format(/8x,'iFF Freq/FF   ig09 freq/g09      overl    err')
   34 format(i10,f9.1,i6,f9.1,f12.3,f8.1)
   38 format(18x,'Standard deviation (cm-1)',f12.2,///)

      end
