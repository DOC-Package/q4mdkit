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
      Subroutine CorrectQMdata(Etot, nator, geom, iflag, par,
     &             parfun, npar, z_scan, maxat, max3n, maxpar) 
c----------------------------------------------------------------------c

c      corrects training QM data (E,G,H) by subtracting the LJ (+charge)
c      contribution

      implicit real*8 (a-h,o-z)
      logical  z_scan

      parameter       (maxint=5000,maxcos=10,maxspl=300)
      parameter       (pi = 3.141592653589793238d0)
      parameter       (rad2deg = 180.d0/pi)
      parameter       (boh2nm  = 0.0529177249d0)!! 1 bohr = 0.0529177 nm
      parameter       (H2kJ    = 2625.5d0     )!! 1 H = 2625.5d0 kJ/mol

      
      dimension iflag(maxpar,0:8),geom(maxat,*)
      dimension par(maxpar),parfun(maxpar)

c     computes LJ and charge-charge contribution to current geometry

      elj = 0.d0 
      ech = 0.d0
      do ip=1,npar

        it=iflag(ip,0)
        ia=iflag(ip,1)
        ib=iflag(ip,2)
        if(it.eq.5) then     !!! intra non-bonded LJ intramolecular
          xa = geom(ia,1)
          xb = geom(ib,1)
          ya = geom(ia,2)
          yb = geom(ib,2)
          za = geom(ia,3)
          zb = geom(ib,3)
          val=dsqrt( (xa-xb)**2 +(ya-yb)**2 +(za-zb)**2 )
          if (z_scan) then  ! geom in A 
             val=val/10.d0
          else              ! geom in bohr
             val=val*boh2nm
          endif
          sigma=parfun(ip)*boh2nm
          eps  =par(ip)*H2kJ
          val1       = sigma/val
          valsq      = val1*val1
          val6       = valsq*valsq*valsq
          val12      = val6*val6
          rm1        = 1.d0/val
          rm2        = rm1*rm1
          eparz      = 4.d0*eps*(val12-val6)
          elj = elj + eparz

        elseif(it.eq.6) then     !!! intramolecular charge - charge interactions

          xa = geom(ia,1)
          xb = geom(ib,1)
          ya = geom(ia,2)
          yb = geom(ib,2)
          za = geom(ia,3)
          zb = geom(ib,3)
          val=dsqrt( (xa-xb)**2 +(ya-yb)**2 +(za-zb)**2 )
          if (z_scan) then  ! geom in A 
             val=val/10.d0
          else              ! geom in bohr
             val=val*boh2nm
          endif
          rm1        = 1.d0/val
          ecparz     = par(ip)*rm1
          ech = ech + ecparz
         endif

      enddo    
      enbtot = elj + ech 
      enbtot = enbtot/H2kJ

      eold = etot
      etot = etot - enbtot


      return 
      end
