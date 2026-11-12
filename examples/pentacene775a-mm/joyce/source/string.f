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
c       Package of Subroutines to Manage Strings
c       written by I. Cacelli                   
c----------------------------------------------------------------------c
 
c----------------------------------------------------------------------c
      Subroutine Firstw (strin,word,lword)
c----------------------------------------------------------------------c
 
 
c     Selects the first word in *strin* and returns it in *word*.
c     The unvisible characters are set to blanks.
 
c     strin   input string
c     word    selected first word in *string*
c     lword   length of the selected word; if no word is found 
c          a 0 value is returned and *word* is a blank string
 
      character*(*) strin,word
      nlw=len(word)
      nls=len(strin)
      if(nls*nlw.le.0) stop 'Wrong string-length in FIRSTW'
      do 5 i=1,nlw
    5 word(i:i)=' '
      do 6 i=1,nls
      ic=ichar(strin(i:i))
    6 if(ic.lt.32.or.ic.gt.126) strin(i:i)=' '
 
      lword=0
      j1=0
   10 j1=j1+1
      if(j1.gt.nls)           return
      if(strin(j1:j1).eq.' ') goto 10
      j2=j1
   20 j2=j2+1
      if(strin(j2:j2).ne.' '.and.j2.lt.nls) goto 20
 
      lword=j2-j1
      if(j2.eq.nls .and. strin(j2:j2).ne.' ') lword=lword+1 
      if(nlw.lt.lword) goto 50
      do 40 i=1,lword
      j=j1+i-1
   40 word(i:i)=strin(j:j)
      return
  
   50 lword=0
      write(6,51) strin(1:Lenw(strin))
   51 format(/6x,'Internal error in FIRSTW: word too small'
     &         ,/1x,a,/)
      stop 'Internal Error in FIRSTW: word too small'
      end
 
c-----------------------------------------------------------------
      Subroutine Selwrd (strin,word,iword,lword,ipos)
c-----------------------------------------------------------------
 
c     selects the iword-th word in *strin* and returns it in *word*.
c     a word corresponds to a sequence of non-blank characters
 
c     strin   input string
c     word    selected word from *string*
c     iword   ordinal number of the wanted word
c     lword   length of the selected word; if the word is not found
c          a 0 value is returned and *word* is a blank string
c     ipos    position of *word* in *string* i.e. *word* matches
c          *string* from the ipos-th to the (ipos+lword-1) char.s
 
      character*(*) strin,word
      nlw=len(word)
      nls=len(strin)
      if(nls*nlw.le.0) stop 'Wrong string-length in SELWRD'
      ipos=0
      do 5 i=1,nlw
    5 word(i:i)=' '
      if(iword.le.0) return
      do 6 i=1,nls
      ic=ichar(strin(i:i))
    6 if(ic.lt.32.or.ic.gt.126) strin(i:i)=' '
 
      lword=0
      j2=0
      do 30 j=1,iword
      j1=j2
   10 j1=j1+1
      if(j1.gt.nls)           return
      if(strin(j1:j1).eq.' ') goto 10
      j2=j1
   20 j2=j2+1
      if(strin(j2:j2).ne.' '.and.j2.lt.nls) goto 20
   30 continue
 
      if(j2.eq.nls .and. strin(j2:j2).ne.' ') j2=j2+1
      lword=j2-j1
      if(nlw.lt.lword) goto 50
 
      i=0
      do 40 j=j1,j2-1
      i=i+1
   40 word(i:i)=strin(j:j)
      ipos=j1
      return
 
   50 lword=0
      write(6,51) strin(1:Lenw(strin))
   51 format(/6x,'Internal error in SELWRD: word too small'
     &         ,/1x,a,/)
      stop 'Internal Error in SELWRD: word too small'
      end


c-----------------------------------------------------------------
      Function Nword (string)
c-----------------------------------------------------------------
 
c     Returns the number of words in string
 
      character string*(*),word*88
      lw=Lenw0(string)
      Nword=0
      if(lw.lt.1) return
      Nword=1
      if(lw.eq.1) return
 
      iw=0
      do i=2,lw
        if(string(i:i).eq.' ') then
          if(string(i-1:i-1).ne.' ') iw=iw+1
        endif
      enddo
      if(string(lw:lw).ne.' ') iw=iw+1
      Nword=iw
      return
      end
 
 
c-----------------------------------------------------------------
      Function Lenw (string)
c-----------------------------------------------------------------
 
c     Lenw = n. of char. of *string* without the final blanks
 
      character*(*) string
      nc=len(string)
      dowhile (string(nc:nc).eq.' ' .and. nc.gt.1 )
      nc=nc-1
      enddo
      Lenw=nc
      return
      end


c-----------------------------------------------------------------
      Function Lenw0 (string)
c-----------------------------------------------------------------
 
c      lenw = n. of char. of *string* without the final blanks
 
      character*(*) string
      nc=len(string)
      dowhile (string(nc:nc).eq.' ' .and. nc.gt.1 )
        nc=nc-1
      enddo
      Lenw0=nc
      if(nc.eq.1 .and. string(1:1).eq.' ') Lenw0=0
      return
      end


c-----------------------------------------------------------------
      Subroutine Delbad (strin)
c-----------------------------------------------------------------
 
 
c     Substitute the unvisible characters in *string* with blanks.
c     The visible characters are from 32 to 126.
 
      character*(*) strin
      nl=len(strin)
      if(nl.le.0) stop 'Wrong string-length in DELBAD'
      do i=1,nl
        ic=ichar(strin(i:i))
        if(ic.lt.32.or.ic.gt.126) strin(i:i)=' '
      enddo
      return
      end

c-----------------------------------------------------------------
      Subroutine PrDate (strdate,nch)
c-----------------------------------------------------------------
 
c     returns a string containing the date informations
 
      character strdate*(*),datefile*10
      datefile='filefortim'
      Call System ('date +"%a %d%b%y %H.%M" >'//datefile)
 
      open (unit=99,file=datefile,status='unknown')
      rewind 99
      Call Setblk (strdate)
      read (99,'(a)') strdate
      nch=Lenw(strdate)
      close (unit=99,status='delete')
 
      return
      end


c-----------------------------------------------------------------
      Subroutine Lowcas (strin)
c-----------------------------------------------------------------
 
c     converts a string to lower case characters
 
      character*(*) strin
      n=len(strin)
      do j=1,n
        is=ichar(strin(j:j))
        if(is.ge.65.and.is.le.90) then
          is=is+32
          strin(j:j)=char(is)
        endif
      enddo
      return
      end


c-----------------------------------------------------------------
      Subroutine Uppcas (strin)
c-----------------------------------------------------------------
 
c     converts a string to upper case characters
 
      character*(*) strin
      n=Len(strin)
      do j=1,n
        is=ichar(strin(j:j))
        if(is.ge.97.and.is.le.122) then
          is=is-32
          strin(j:j)=char(is)
        endif
      enddo
      return
      end



c-----------------------------------------------------------------
      Subroutine Getnum (word,iret,fnum,inum)
c-----------------------------------------------------------------
 
c     Checks if the first sequence of contiguous characters
c     in *word* may be read in a numerical fortran format
c     The results are stored in *fnum* and *inum*
 
c        iret =-2  wrong input string or blank string
c        iret =-1  data not readable (may be alphanumeric)
c        iret = 0  integer resulting number
c        iret = 1  floating point resulting number
 
      implicit logical (z)
      real*8 fnum
      character*1 allow(16),chj
      character*(*) word
      data allow /'1','2','3','4','5','6','7','8','9','0'
     &           ,'+','-','d','D','e','E'/
c     to prevent accidental case (upper,lower) changes
      allow(13)=char(68)
      allow(14)=char(69)
      allow(15)=char(100)
      allow(16)=char(101)
 
      inum=0
      fnum=0.d0
      lw=len(word)
      iret=-2
      if(lw.eq.0) return
 
      j1=1
      dowhile (j1.lt.lw .and. word(j1:j1).eq.' ')
      j1=j1+1
      enddo
      j2=j1+1
      dowhile (j2.le.lw .and. word(j2:j2).ne.' ')
      j2=j2+1
      enddo
      j2=j2-1
 
c     check integer number
 
      iret=0
      nok=12
      do j=j1,j2
      ii=0
        do i=1,nok
        if(word(j:j).eq.allow(i)) ii=i
        enddo
      if(ii.eq.0) iret=-1
      nok=10
      enddo
      if(iret.eq.0) then
      read(word,*) inum
      fnum=inum
      return
      endif
 
c     check floating point number
 
      if((j2-j1+1).lt.2) return
      zs=.true.
      zp=.true.
      ndig=0
      ndex=0
      j=j1-1
   10 j=j+1
      if(j.gt.j2) goto 40
      chj=word(j:j)
 
      if(zp.and.chj.eq.'.') then
      zs=.false.
      zp=.false.
      goto 10
      endif
 
      if(zs.and.(chj.eq.'+'.or.chj.eq.'-')) then
      zs=.false.
      goto 10
      endif
c     digits
      ii=0
      do i=1,10
      if(chj.eq.allow(i)) ii=i
      enddo
      if(ii.gt.0) then
      ndig=ndig+1
      zs=.false.
      goto 10
      endif
c     check  E,D,e,d
      ii=0
      do i=13,16
      if(chj.eq.allow(i)) ii=i
      enddo
      if(ii.eq.0) return
 
c     exponent
 
      if(ndig.lt.1) return
      zs=.true.
   20 j=j+1
      chj=word(j:j)
      if(j.gt.j2) goto 30
 
      if(zs.and.(chj.eq.'+'.or.chj.eq.'-')) then
      zs=.false.
      goto 20
      endif
c     digits
      ii=0
      do i=1,10
      if(chj.eq.allow(i)) ii=i
      enddo
      if(ii.gt.0) then
      ndex=ndex+1
      zs=.false.
      goto 20
      endif
      return
 
   30 continue
      if(ndex.eq.0.or.ndex.gt.2) return
   40 read(word,*) fnum
      inum=fnum
      iret=1
      return
      end


c--------------------------------------------------------------------      
      Subroutine Setblk (string)
c--------------------------------------------------------------------      
 
c     fill string with blanks
 
      character*(*) string
      do j=1,len(string)
        string(j:j)=' '
      enddo
      return
      end



c--------------------------------------------------------------------      
      Subroutine Setgrat (string)
c--------------------------------------------------------------------      
 
c     substitute blanks from the '#' onward in string
 
      character*(*) string
      nc=Len(string)
      ik=Index(string,'#')
      if(ik.le.0) return
      do k=ik,nc
        string(k:k)=' '
      enddo
      return
      end

c--------------------------------------------------------------------      
      Subroutine Setbl1 (string)
c--------------------------------------------------------------------      
 
c     substitute blanks to the first word of *string*
 
      character*(*) string
      nc=len(string)
      ic=1
      dowhile (string(ic:ic).eq.' ' .and. ic.lt.nc )
      ic=ic+1
      enddo
 
      dowhile (string(ic:ic).ne.' ' .and. ic.lt.nc )
      string(ic:ic)=' '
      ic=ic+1
      enddo
      string(ic:ic)=' '
      return
      end



c--------------------------------------------------------------------      
      Subroutine Compact (word,lw)
c--------------------------------------------------------------------      
 
c     Compact the non-blank characters of *word* 
 
      character*(*) word
      llw=len(word)
      lw=0
      if(llw.eq.0) return
      do j=1,llw
        if(word(j:j).ne.' ') then
          lw=lw+1
          if(lw.lt.j) then
            word(lw:lw)=word(j:j)
            word(j:j)=' '
          endif
        endif
      enddo
      return
      end



c--------------------------------------------------------------------      
      Subroutine FmtAdj (card,nch)
c--------------------------------------------------------------------      
 
c   adjust the string *card* in the following way:
c   1. the contiguous blanks are reduced to 1 blank
c   2. the right adjusted zeros of real numbers are dropped
 
c   card  (inp)  input string (with literal words and numbers)
c         (out)  adjusted-format string 
c   nch   (out)  total number of final char. in *card*
 
      parameter (max=200)
      character*(*) card
      character*200 word
      logical zkeep(max)
 
      lwc=len(card)
      nch=0
      if(lwc.gt.max) goto 200
      iw=0
   10 iw=iw+1
      Call Selwrd (card,word,iw,lw,ipos)
      if(lw.le.0) goto 20
      if(index(word(1:lw),'.').le.0) goto 10
      if(lw.le.4)                    goto 10
      do j=ipos+lw-1,ipos+2,-1
      if(card(j:j).ne.'0') goto 10
      card(j:j)=' '
      enddo
      goto 10
 
   20 if(iw.eq.1) return
      zkeep(1)=.true.
      do j=1,lwc-1
      zkeep(j+1)=.true.
      if(card(j:j+1).eq.'  ') zkeep(j+1)=.false.
      enddo
      card(lwc:lwc)=' '
      j=0
      do k=1,lwc
      if(zkeep(k)) then
      j=j+1
      card(j:j)=card(k:k)
      if(k.gt.j) card(k:k)=' '
      endif
      enddo
      nch=j
      return
 
  200 write(6,202) lwc
  202 format(/,/6x,'Program Error in routine FMTADJ'
     &        ,/6x,'insufficient dimensions',i6,/)
      stop
      end



c--------------------------------------------------------------------      
      Subroutine Skip_Blk (iun)
c--------------------------------------------------------------------      
 
c     Skips blank lines until a non blank line is found
 
      character*60 card,cblk
      do i=1,Len(card)
        card(i:i)=' '
        cblk(i:i)=' '
      enddo
      dowhile (card.eq.cblk)
        read(iun,'(a)',end=20) card
      enddo
      backspace iun
      return
   20 write(6,'(/4x,"Error: end-of-file on unit",i4,/)') iun
      stop 'Error: end-of-file on unit'
      end



c--------------------------------------------------------------------      
      Subroutine Delblk (strin)
c--------------------------------------------------------------------      
 
c   remove the initial blank characters of string and
c   move backwards the remaining non-blank charaters 
c    deleting the successive not-blank characters
 
      character*(*) strin
      nlen=Len(strin)
 
      if(strin(1:1).ne.' ') return
      i=1
      dowhile (i.le.nlen)
c       Placing the test here avoids evaluating strin(nlen+1), which is out of bounds
        if (strin(i:i).ne.' ') exit
        i=i+1
      enddo
      if(i.gt.nlen) return
      ik=0
      do k=i,nlen
        ik=ik+1
        strin(ik:ik)=strin(k:k)
      enddo
      do i=ik+1,nlen
      strin(i:i)=' '
      enddo
      return
      end


c--------------------------------------------------------------------      
      Subroutine Get_Path (file,n)
c--------------------------------------------------------------------      
 
c     Get the absolute path of the present file and remove blanks
 
      character*(*) file
      Call Get_First (file)
      Call Move_to_SCR (file)
      n=Lenw(file)
      return
      end


c--------------------------------------------------------------------      
      Subroutine Get_First (strin)
c--------------------------------------------------------------------      
 
c   sort the first word in strin deleting the initial blanks
c   WARNING: string is changed
 
      character strin*(*),strin1*300
      nlen=Len(strin)
      if(nlen.gt.Len(strin1)) stop 'ERROR in Get_First: increase dim'
      do 5 j=1,nlen
    5 strin1(j:j)=' '
 
      k=0
      j=0
   10 j=j+1
      if(j.le.nlen) goto 14
      write(6,12)
   12 format(//5x,'error in Get_First routine : input blank string'/)
      stop
   14 if(strin(j:j).eq.' ') goto 16
      k=k+1
      strin1(k:k)=strin(j:j)
      goto 10
   16 if(k.eq.0) goto 10
      nl=k
 
      do 20 i=1,nlen
   20 strin(i:i)=strin1(i:i)
      return
      end


c--------------------------------------------------------------------      
      Subroutine Read_Opt (strin,opt,zfound)
c--------------------------------------------------------------------      
 
c   check if opt is contained in string: in that case remove
c   opt from string and set zfound=T
c   WARNING: opt must contain at least two characters
 
      character strin*(*),opt*(*)
      logical zfound
      zfound=.false.
      ll=Lenw(opt)
      if(ll.le.1) return
      ii=index(strin,opt)
      if(ii.le.0) return
      zfound=.true.
      do i=ii,ii+ll-1
        strin(i:i)=' '
      enddo
      return
      end

c--------------------------------------------------------------------      
      Subroutine SortS (card,word,nc,ch1,n1,n2)
c--------------------------------------------------------------------      
 
c    Sorts a string enclosed into quotes or double-quotes
c    The choice between ' and " (as bounds of the string) is
c    determined by the first occurrence of one of them.
 
c  card   (i)   chr    input card
c  word   (o)   chr    selected string (without the bounding quotes)
c  nc     (o)   int    number of characters of the selected string
c                      0 is returned if no string was found
c  ch1    (i)   chr*1  char. to replace the string in *card*
c  n1,n2  (o)   int    position of the bounding quotes in *card*
 
      character card*(*),word*(*),ch1*1,q1*1,q2*1
      common /units/ inp,iou,ipu,ims,iun,iterm,inull
      nc=0
      do i=1,Len(word)
        word(i:i)=' '
      enddo
      q1=char(39)
      q2=char(34)
      n1=index(card,q1)    ! searches for '
      n2=index(card,q2)    ! searches for "
      if(n1.le.0) then
        if(n2.le.0) return
        n1=n2
        q1=q2
      else
        if(n2.gt.0 .and. n2.lt.n1) then
          n1=n2
          q1=q2
        endif
      endif
 
      nl=Lenw(card)
      do i=nl,n1+1,-1
        n2=i
        if(card(i:i).eq.q1) goto 20
      enddo
      write(ims,18) q1,card(1:nl)
   18 format(/' WARNING: only one begin-string /',a,'/ is found in:'
     &      ,/'  /',a,'/',/)
      nc=0
      n1=0
      return
 
   20 nc=n2-n1-1
      if(nc.gt.0) then
        word(1:nc)=card(n1+1:n2-1)
      else
        nc=1
        word(1:1)=' '
      endif
      do i=n1,n2
        card(i:i)=ch1
      enddo
      return
      end


c--------------------------------------------------------------------
      Subroutine Cerca_W (key,iun,card,ier,zrew)
c--------------------------------------------------------------------
 
c     Search for the first line containing the string /key/
c     (from the current position) and returns the remaining part
c     of the read line in /card/
c     WARNING: this version look for the string /key/ preceded by a blank
 
      implicit real*8(a-h,o-y), logical(z)
      character key*(*),skey*30,card*(*),card0*99,card1*99,back*1
      logical zdeb
      back=char(92)
      zdeb=.false.
 
      iword=0  !! 0 = complete string , 1 = first word of /key/ 
      icase=0  !! 0 = case independent , 1 = case dependent
 
      ier=1
      if(iword.eq.1) then
        Call Selwrd (key,skey,1,lws,ipos)
        if(lws.le.1) goto 30
      else
        skey=' '
        lws=Lenw(key)
        skey=key(1:lws)
      endif
      if(icase.eq.0) Call Lowcas (skey)
 
      if(zrew) rewind iun
      ic=0
   10 read(iun,'(a)',end=20) card
      card1=card
      Call Lowcas(card1)
      if(index(card1,'nd-of-fil').gt.0) goto 20
      if(index(card1,'nd of fil').gt.0) goto 20
      ic=ic+1
      Call SetGrat (card)
      if(Lenw(card).le.1)  goto 10
      card0=card
      if(icase.eq.0) Call Lowcas (card0)
      if(zdeb) write(6,'(i5,2x,a)') ic,card0(i:Lenw(card))
      ii=Index(card0,skey(1:lws))
      if(ii.gt.0) then
        do i=1,ii
          card(i:i)=' '
        enddo
        ii=ii+1
        ier=0
        lw=Lenw(card)
        dowhile (card(ii:ii).ne.' ')
          card(ii:ii)=' '
          ii=ii+1
        enddo
        if(card(lw:lw).eq.back) then
          read(iun,'(a)') card1
          lw1=Lenw(card1)
          card(lw+1:lw+lw1)=card1(1:lw1)
          lw=Lenw(card)
        endif
        Call DelBlk (card)
        return
      endif
      goto 10
   20 return
   30 stop 'ERROR in Cerca_W: null key string'
      end


c--------------------------------------------------------------------
      Subroutine Cerca_L (key,iun,card,ier,zrew,iline)
c--------------------------------------------------------------------
 
c   Search for the first line containing the string /key/
c   (from the current position) and returns the remaining part
c   of the read line in /card/. If the selected card ends with backslash '\'
c   the next card is added to the current card, and so on.
c
      implicit real*8(a-h,o-y), logical(z)
      character key*(*),skey*30,card*(*),card1*90,back*1
      logical zdeb
      back=char(92)
      zdeb=.false.
 
      iword=0  !! 0 = complete string , 1 = first word of /key/ 
      icase=1  !! 0 = case independent , 1 = case dependent
 
      ier=1
      if(iword.eq.1) then
        Call Selwrd (key,skey,1,lws,ipos)
        if(lws.le.1) goto 30
      else
        skey=' '
        lws=Lenw(key)
        skey=key(1:lws)
      endif
      if(icase.eq.0) Call Lowcas (skey)
 
      if(zrew) rewind iun
      if(zrew) iline=0
      ic=0
   10 read(iun,'(a)',end=20) card
      card1=card
      Call Lowcas(card1)
      if(index(card1,'nd-of-fil').gt.0) goto 20
      if(index(card1,'nd of fil').gt.0) goto 20
      if(index(card1,'nd.of.fil').gt.0) goto 20
      iline=iline+1
      ic=ic+1
      if(Lenw(card).le.1)  goto 10
      if(icase.eq.0) Call Lowcas (card)
      if(zdeb) write(6,'(i5,2x,a)') ic,card(1:Lenw(card))
      ii=Index(card,skey(1:lws))
      if(ii.gt.0) then
        ier=0
        do i=1,ii+lws-1
          card(i:i)=' '
        enddo
        i=ii+lws
        dowhile (card(i:i).ne.' ')
          card(i:i)=' '
          i=i+1
        enddo
 
        lw=Lenw(card)
        dowhile(card(lw:lw).eq.back)
          read(iun,'(a)') card1
          iline=iline+1
          lw1=Lenw(card1)
          card(lw:lw+lw1-1)=card1(1:lw1)
          lw=Lenw(card)
        enddo
        if(zdeb) write(6,'(i5,2x,a)') ic,card(1:Lenw(card))
        Call DelBlk (card)
        return
      endif
      goto 10
   20 return
   30 stop 'ERROR in Cerca_L: null key string'
      end


c--------------------------------------------------------------------      
      Subroutine Read_one_card (iun,card,iline)
c--------------------------------------------------------------------

c     read the current card with the following specificities:
c     1 - if the last char of card is '\' read the next card and
c         add the content to the current card, and so on
c     2 - if the card begins with '#' or '%' then skip to next card
c     3 - if the card is blank, then skip to next card
 
      implicit real*8(a-h,o-y), logical(z)
      character card*(*),card1*190,back*1
      back=char(92)
      lwc=Len(card)
 
   10 read(iun,'(a)',end=20) card
      iline=iline+1
      if(index(card(1:2),'#').gt.0) goto 10
      if(index(card(1:2),'%').gt.0) goto 10
      lw=Lenw(card)
      dowhile(card(lw:lw).eq.back)
        read(iun,'(a)',end=20) card1
        iline=iline+1
        lw1=Lenw(card1)
        if(lw+lw1-1.gt.lwc) goto 30
        card(lw:lw+lw1-1)=card1(1:lw1)
        lw=Lenw(card)
      enddo
      Call DelBlk (card)
      return
   20 write(6,21) 
   21 format(//8x,'ERROR in READ_ONE_CARD: end of file',//)
      stop 
   30 write(6,31) lw+lw1-1
   31 format(//8x,'ERROR in READ_ONE_CARD: increase length of card',//)
      stop 
      end


c--------------------------------------------------------------------
      Subroutine Skip_File (iun,nline)
c--------------------------------------------------------------------
 
c     skip the file /iun/ of /nline/ lines
 
      character*2 card
      rewind iun
      do i=1,nline
        read(iun,'(a)',end=4) card
      enddo
      return
    4 write(6,8) iun,nline
    8 format(//5x,'ERROR in SKIP_FILE unit',i4,5x,'skip =',i7,///)
      stop
      end



c--------------------------------------------------------------------
      Subroutine Read_Seq (card,isel,nsel,nmin,nmax)
c--------------------------------------------------------------------
 
c     reads a sequence of integer numbers 
c     with the following convention
c     8 - 12  ---> 8 9 10 11 12
c     extracting only numers between min and max
 
      implicit real*8(a-h,o-z)
      dimension isel(*)
      character card*(*),word*8
 
      Call Count_W (card,nword)
      i1=1
      nsel=0
      ii=0
   10 ii=ii+1
      if(ii.gt.nword)       goto 40
      Call SelWrd (card,word,ii,lw,ipos)

      if(word(1:2).eq.'; ') goto 40
      if(word(1:1).eq.'#')  goto 40
      if(word(1:2).eq.'. ') goto 40
      if(word(1:2).eq.'  ') goto 40
 
      if(word(1:2).eq.'- ') then
        ii=ii+1
        Call SelWrd (card,word,ii,lw,ipos)
        if(ii.gt.nword) goto 40
        read(word,*) i2
        do i=Min0(i1,i2)+1,Max0(i1,i2)
          if(i.ge.nmin .and. i.le.nmax) then
            nsel=nsel+1
            isel(nsel)=i
          endif
        enddo
      elseif (index(word,'-') .gt. 1) then 
        write(6,*) 
        write(6,*) 
        write(6,'(10x,a)') '*********************' 
        write(6,'(14x,a)') 'ERROR: syntax error !'
        write(6,'(7x,a)')
     &     'Please specifiy number sequence in the correct format:'
        write(6,'(10x,a)') ' 1 - 19 or  3 7 8 -20'
        write(6,'(10x,a)') '*********************' 
        write(6,*) 
        write(6,*) 
        stop 'I/O error' 
      else
        read(word,*) i1
        if(i1.ge.nmin .and. i1.le.nmax) then
          nsel=nsel+1
          isel(nsel)=i1
        endif
      endif
      goto 10
   40 continue
      if(nsel.le.1) return
 
c  rearrange isel in increasing order and eliminates the doubles
 
      do i=2,nsel
        do k=1,i-1
          if(isel(i).eq.isel(k)) isel(i)=0
        enddo
      enddo
      j=0
      do i=1,nsel
        if(isel(i).gt.0) then
          j=j+1
          isel(j)=isel(i)
        endif
      enddo
      if(j.lt.nsel) then
        do k=j+1,nsel
          isel(k)=0
        enddo
      endif
      nsel=j
 
      return
      end
 

c--------------------------------------------------------------------
      Subroutine Parse (A,B,nwrd,maxwrd)
c--------------------------------------------------------------------
 
c     copy the words in /A/ in the /B/ array
 
      character*(*) A,B(*)
      do i=1,maxwrd
        B(i)=' '
      enddo
      lena=Len(A)
      lenb=Len(B(1))
      nwrd=0
      i=0
   10 i=i+1
      if (i.gt.lena) return
      if (A(i:i).eq.' ') goto 10
      if(nwrd.eq.maxwrd) return
      nwrd=nwrd+1
      ibeg=i
   20 i=i+1
      if (i.gt.lena) then
      iend=i-1
      B(nwrd)=' '
      do k=ibeg,iend
        kb=k-ibeg+1
        B(nwrd)(kb:kb)=A(k:k)
        if(kb.le.lenb) B(nwrd)(kb:kb)=A(k:k)
      enddo
      return
      endif
      if (A(i:i).ne.' ') goto 20
      iend=i-1
      B(nwrd)=' '
      do k=ibeg,iend
        kb=k-ibeg+1
        B(nwrd)(kb:kb)=A(k:k)
        if(kb.le.lenb) B(nwrd)(kb:kb)=A(k:k)
      enddo
      goto 10
      end


c--------------------------------------------------------------------
      Subroutine Count_W (string,nwrd)
c--------------------------------------------------------------------
 
c     count the words of string
 
      character*(*) string
      lens=Len(string)
      nwrd=0
      i=0
   10 i=i+1
      if(i.gt.lens) goto 30
      if(string(i:i).eq.' ') goto 10
      nwrd=nwrd+1
      ibeg=i
   20 i=i+1
      if(i.gt.lens) goto 30
      if(string(i:i).ne.' ') goto 20
      goto 10
   30 return
      end


 
c--------------------------------------------------------------------
      Subroutine Via_Numer (word)
c--------------------------------------------------------------------
 
c     discard the numerical characters in word and compact
 
      character*(*) word
      character*1 num(11)
      data num/'1','2','3','4','5','6','7','8','9','0','.'/
      ll=Lenw(word)
      do i=1,ll
        do k=1,11
          if(word(i:i).eq.num(k)) word(i:i)=' '
        enddo
      enddo
      k=0
      do i=1,ll
        if(word(i:i).ne.' ') then
          k=k+1
          word(k:k)=word(i:i)
          if(i.gt.k) word(i:i)=' '
        endif
      enddo
      return
      end



c--------------------------------------------------------------------
      Subroutine Count_Numer (word,inum,ialfa)
c--------------------------------------------------------------------
c
c     counts character numbers in word and gives back
c         inum =  number of numerical characters (0...9)
c         ialf =  other characters
c     warning: /word/ should be without blanks
 
      character*(*) word
      character*1 num(10)
      data num/'1','2','3','4','5','6','7','8','9','0'/
      ll=Lenw(word)
      ii=0
      do i=1,ll
        do k=1,10
          if(word(i:i).eq.num(k)) ii=ii+1
        enddo
      enddo
      inum=ii
      ialfa=ll-inum
      return
      end


c--------------------------------------------------------------------
      Subroutine Check_Dim (need,nreal,string)
c--------------------------------------------------------------------
 
 
      character*(*) string
      if(need.gt.nreal) then
        ll=Lenw(string)
        write(6,8) string,need,nreal
        stop
      endif
      return
    8 format(//8x,'DIMENSION ERROR in ',a
     &       ,/8x,'requested dimension .............',i6
     &       ,/8x,'effective dimension .............',i6,//)
      end
 
c--------------------------------------------------------------------
      Subroutine Sort_Max (a,b,inda,na,nb)
c--------------------------------------------------------------------
 
c     extraxt max, min nb values from navalues in a array 
 
      implicit real*8(a-h,o-y), logical(z)
      dimension a(*),b(*),inda(*)
 
      if(nb.ge.na) then
        write(6,'(//" ERROR in SORT_MAX   nb>=na",2i6)') na,nb
        stop
      endif
      if(na.le.0) then
        write(6,'(//" ERROR in SORT_MAX   na<=0",i6)') na
        stop
      endif
 
      am=0.0d0
      do i=1,na
        if(dabs(a(i)).gt.am) then
          am=dabs(a(i))
          im=i
        endif
      enddo
      b(1)=a(im)
      inda(1)=im
 
      if(nb.le.1) return
      do ib=2,nb
        am=0.0d0
        do i=1,na
          zok=.true.
          do j=1,ib-1
            if(inda(j).eq.i) zok=.false.
          enddo
          if(zok .and. dabs(a(i)).gt.am) then
            am=dabs(a(i))
            im=i
          endif
        enddo
        b(ib)=a(im)
        inda(ib)=im
      enddo
 
      return
      end


c--------------------------------------------------------------------
      Subroutine Move_to_SCR (file)
c--------------------------------------------------------------------
 
c     subtitutes the conventional directory name 'PEGSCR'
c     with the $PEGSCR environment variable
 
      character file*(*),file1*120,word*50
      file1=file
      Call LowCas (file1)
      if(Index(file1(1:9),'pegscr').gt.0) then
        ii=index(file1(1:9),'scr')
        Call Getenv ('PEGSCR',word)
        if(Lenw(word).lt.2) stop 'ERROR: env. var. PEGSCR not set'
        file1=' '
        file1=word(1:Lenw(word))//'/'//file(ii+3:Lenw(file))
        do i=2,Lenw(file1)
          if(file1(i-1:i).eq.'//') file1(i-1:i)=' /'
        enddo
        file=file1
        Call Compact (file,lw)
      endif
      return
      end

c--------------------------------------------------------------------
      Subroutine Setcomm (string)
c--------------------------------------------------------------------
 
c   substitutes blanks from the ';' onward in string
c
      character*(*) string
      nc=Len(string)
      ik=Index(string,';')
      if(ik.le.0) return
      do k=ik,nc
        string(k:k)=' '
      enddo
      return
      end
c

