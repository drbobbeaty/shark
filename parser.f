c
c     What follows are some routines that will make the use of the
c     simulator MUCH easier. These routines for the parser et. al.
c     were originally written by Jeff Gray of Purdue University.
c
c     The motivating factor is that the simulator does NOT need to
c     be recompiled each time it is run if there is a reasonably
c     advanced parser to get the imformation specific to this run.
c
c

      subroutine parser(idev,ierr)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	reserved common block name:
c		/parsd1/
c		/parsd2/
c
c	reserved subprogram names (called by parser):
c		gnzzz
c		getzzz
c		izzzz
c		zzzr
c		gdzz
c		zzzip
c		zzzpk
c
c
c	this subroutine reads from the file idev 'cards'
c	which are of the form:
c
c	*titlecard anything at all
c
c	or
c
c	cardname var1=value1,var2=string2 var3=value3
c	+    array1=va1/va2/va3, array2=sa1/sa2/sa3/sa4
c
c	line beginning with a blank or a $ are ignored
c
c	a title card can only be one line. its name must
c	start in column 1 with a *. The contents of the
c	title line is put into the array ititle in 8a8
c	format. also, kflag=1.
c
c
c	a card can be up to five lines long. The card name
c	must start in column 1. the continuation symbol, +,
c	must also appear in column 1. for this type of card
c	kflag=0.
c
c	commas or blanks are assumed to be separaters. any
c	number of separators may appear between assignments.
c	an assignment cannot contain any blanks, i.e.
c
c		block = - 12.0
c
c	is not valid. It must be written as
c
c		block=-12.0
c
c	instead.
c
c	multiple values can be assigned to the variable
c	by separating the values by /'s. For example
c
c		wl=.3/.4/.5/.6/.7/.8
c
c	and
c
c		ohmic=yes/no/no/no/yes
c
c	The subroutine parser examines the card and puts
c	the results in the common blocks /parsd1/ and /parsd2/.
c	these are as follows:
c
c	ncard		the cardname (a8 format)
c	kflag		relates the type of card read
c			=-1, end of file reached
c			=0, comment card, ignored
c			=1, card contains assignments
c			=2, card contains a title
c	numr		if kflag=1, numr is the number
c			of real values found on the card (<=70)
c	nreal		an array containing the list of
c			real variable names found on the
c			card (a8 format)
c	vreal		an array containing the values
c			corresponding to the names in nreal
c	numc		if kflag=1, numc is the number of
c			character constants found on the card (<=30)
c	nchar		an array containing the list of
c			character constant names found on
c			the card (a8 format)
c	vchar		an array containing the strings
c			corresponding to the names in nchar
c			(a8 format)
c	ititle		if kflag=2, ititle contains the comment
c			on the title card in 10a8 format
c	ibuf		contains the entire card in 400a1 format
c
c
c	subroutine arguments:
c
c	idev	the unit number of the file containing the
c		input deck
c
c	ierr	=0 if no syntax error is detected
c		error codes-
c
c		101	the continuation character is the
c			first character on the card. Possible
c			problems:
c			  a title card cannot be continued
c			  a cad cannot be more than 5 lines
c
c		102	unrecognized value
c
c		103	character string greater than 8 characters
c			cardnames, variable names, and strings
c			cannot be more than 8 characters
c
c		104	value contains more than 20 characters.
c			real values containing more than 20
c			characters are illegal.
c
c		105	unrecognizable real value.
c
c
c
c	example-
c
c	column
c	1
c
c	*mess  hello world!
c
c		ncard='*mess   '
c		kflag=2
c		numr=0
c		numc=0
c		ititle(1)=' hello w'
c		ititle(2)='world!  '
c
c
c	newcard  x=3.0,rat=cat  alpha=1.0e12
c	+  wl=.3/.4/.5 string=ab/cde/fghi/j
c
c		ncard='newcard '
c		kflag=1
c		numr=5
c		nreal(1)='x       '
c		nreal(2)='alpha   '
c		nreal(3)='wl      '
c		nreal(4)='wl      '
c		nreal(5)='wl      '
c		vreal(1)=3.0
c		vreal(2)=1.0e12
c		vreal(3)=.3
c		vreal(4)=.4
c		vreal(5)=.5
c		numc=5
c		nchar(1)='rat     '
c		nchar(2)='string  '
c		nchar(3)='string  '
c		nchar(4)='string  '
c		nchar(5)='string  '
c		vchar(1)='cat     '
c		vchar(2)='ab      '
c		vchar(3)='cde     '
c		vchar(4)='fghi    '
c		vchar(5)='j       '
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      common /parsd1/ kflag, numr, numc, vreal(70)
      integer kflag, numr, numc
      real vreal
      common /parsd2/ ncard, nreal(70), nchar(30),
     +                vchar(30), ititle(10), ibuf(400)
      character*8 ncard, nreal, nchar, vchar, ititle
      character*1 ibuf
c
c     ncard		card name
c     kflag		kind of card
c			   = -1  end of file reached
c			   =  0  comment
c			   =  1  normal card
c			   =  2  title card
c
c
c
c
c
      dimension nbuf(8)
      character*1 nbuf
      character*8 ntemp, ipack, cval
      integer ierr, idev, i, j, k1, k2, ki, ic, k, itype, imax
      real rval
c
      ierr=0
      ncard='        '
      kflag=0
      numr=0
      numc=0
      do 901 i=1,400
  901   ibuf(i)=' '
      do 902 i=1,8
  902   ititle(i)='        '
c
c     read all lines of card (maximum of 5)
c
c     if first character is a $ or a space, the line is a
c     comment line, and is ignored.
c
c     if the first character is a *, the card name is assumed
c     to be of the form *cccccc, where c is any character.
c     the entire line after *cccccc is put into the array ititle
c     and the line is not further parsed.
c
c     if the first character of the line is a +, an error is
c     assumed since + is reserved as a continuation character
c     between lines of the same card.
c
c

 1    continue
      read(idev,100,end=999)(ibuf(i),i=1,80)
 100  format(80a1)
      i = 1
 110  if ( .not. ( (ibuf(i+1) .eq. ' ') .and.
     +           (ibuf(i+2) .eq. ' ') .and. (ibuf(i+3) .eq. ' ') .and.
     +           (ibuf(i+4) .eq. ' ') .and. (ibuf(i+5) .eq. ' ') .and.
     +           (ibuf(i+6) .eq. ' ') .and. (ibuf(i+7) .eq. ' ') .and.
     +           (ibuf(i+8) .eq. ' ') .and. (ibuf(i+9) .eq. ' ') .and.
     +           (ibuf(i+10) .eq. ' ') ) ) then
        i=i+1
        goto 110
      endif
      imax=i
      write(6,200) (ibuf(i),i=1,imax)
 200  format(1x,80a1)
      if (ibuf(1) .eq. '$') then
        go to 1
      elseif (ibuf(1) .eq. ' ') then
        go to 1
      elseif (ibuf(1) .eq. '+') then
        ierr=101
        return
      endif
c
      do 101 j=1,4
        k1=j*80+1
        k2=k1+79
        read(idev,100,end=102) (ibuf(i),i=k1,k2)
        if (ibuf(k1) .eq. '+') then
          write(6,200) (ibuf(i),i=k1,k2)
          ibuf(k1)=' '
        else
          backspace idev
          go to 102
        endif
 101  continue
      go to 104
 102  continue
      k2=k2-80
      do 103 ki=k2+1,400
 103    ibuf(ki)=' '
 104  continue
c
c     decode the card name
c
      ic=1
      call gnzzz(ic,ibuf,' ',ncard,ierr)
c
c     check for *______ form of card name
c
      if (ibuf(1) .eq. '*') then
        kflag=2
        do 23 i=1,10
          call zzzip(nbuf,8)
          k=8*(i-1)
          do 22 j=1,8
 22         nbuf(j)=ibuf(ic+k+j-1)
          call zzzpk(ipack,nbuf)
 23       ititle(i)=ipack
        return
      endif
c
c
c     find start of next variable name
c
      kflag=1
      ic=ic-1
 2    continue
      if (ic .ge. k2) go to 3
      ic=ic+1
      if (ibuf(ic) .eq. ' ') go to 2
      if (ibuf(ic) .eq. ',') go to 2
c
c     decode the variable name
c
      call gnzzz(ic,ibuf,'=',ntemp,ierr)
c
c     decode variable value, given type (real=1, char=2)
c
 12   continue
      call getzzz(ic,ibuf,',',rval,cval,itype,ierr)
      if (iabs(itype) .eq. 1) then
        numr=numr+1
        nreal(numr)=ntemp
        vreal(numr)=rval
      elseif (iabs(itype) .eq. 2) then
        numc=numc+1
        nchar(numc)=ntemp
        vchar(numc)=cval
      else
c       unregocnizable constant
        ierr=102
        return
      endif
      if (itype .lt. 0) go to 12
      ic=ic-1
      go to 2
c
 3    continue
      return
 999  kflag=-1
      return
      end
c
c
c
c
      subroutine zzzip(iz,n)
      integer i,n
      character*1 iz(n)
      do 1 i=1,n
 1      iz(i)=' '
      return
      end
c
c
c
c
c
      subroutine gnzzz(ic,ibuf,istop,ntemp,ierr)
      dimension ibuf(400),nbuf(8)
      character*1 nbuf,ibuf,istop
      character*8 ntemp,ipack
      integer ic,ierr,nc
c
      call zzzip(nbuf,8)
      nc=1
      nbuf(nc)=ibuf(ic)
 1    ic=ic+1
      if (ibuf(ic) .eq. istop) go to 2
      nc=nc+1
      if (nc .gt. 8) then
        ierr=103
        go to 1
      endif
      nbuf(nc)=ibuf(ic)
      go to 1
 2    continue
      call zzzpk(ipack,nbuf)
      ntemp=ipack
      ic=ic+1
      return
      end
c
c
c
c
c
      subroutine getzzz(ic,ibuf,istop,rval,cval,ityp,ierr)
      dimension ibuf(400),nbuf(20)
      character*1 nbuf,ibuf,istop
      character*8 ipack,cval
      integer ic,ityp,ierr,nc,izzzz
      real rval,zzzr
      logical invalid
c
      call zzzip(nbuf,20)
      nc=1
      nbuf(nc)=ibuf(ic)
 1    ic=ic+1
      if (ibuf(ic) .eq. ' ') go to 2
      if (ibuf(ic) .eq. istop) go to 2
      if (ibuf(ic) .eq. '/') go to 2
      nc=nc+1
      if (nc .gt. 20) then
        ierr=104
        go to 1
      endif
      nbuf(nc)=ibuf(ic)
      go to 1
 2    continue
      ityp=izzzz(nbuf,20)
      if (ibuf(ic) .eq. '/') ityp=-ityp
      if (iabs(ityp) .eq. 1) then
        rval=zzzr(nbuf,nc,invalid)
        if (invalid) ierr=105
      elseif (iabs(ityp) .eq. 2) then
        call zzzpk(ipack,nbuf)
        cval=ipack
        if (nc .gt. 8) ierr=103
      else
c       inrecognizable constant
        ierr=102
        return
      endif
      ic=ic+1
      return
      end
c
c
c
c
c
      subroutine zzzpk(ipack,c)
      character*1 ipack(8)
      character*1 c(8)
      integer i
c
c     concatenate characters
c
      do 1 i=1,8
 1      ipack(i)=c(i)
      return
      end
c
c
c
      integer function izzzz(nbuf,nc)
      integer i,nc
      dimension nbuf(nc),idig(13)
      character*1 nbuf,idig
      data idig/'+','-','.','0','1','2','3','4','5','6','7','8','9'/
      izzzz=2
      do 1 i=1,13
        if (nbuf(1) .eq. idig(i)) izzzz=1
 1    continue
      return
      end
c
c
      real function zzzr(nbuf,nc,invalid)
      integer nc,ndec,i,idec
      dimension nbuf(nc)
      character*1 nbuf
      logical invalid
      real value,sign,x,gdzz,xpon,se,e
c
      ndec=0
      zzzr=0.0
      invalid=.false.
      if (nc .eq. 0) then
        invalid=.true.
        return
      endif
c
c     decode mantissa
c
      value=-1.0
      sign=1.0
      if (nbuf(1) .eq. '-') then
        sign=-1.0
        i=2
      elseif (nbuf(1) .eq. '+') then
        i=2
      elseif (nbuf(1) .eq. '.') then
        ndec=1
        idec=1
        i=2
      else
        i=1
      endif
c
      if (i .gt. nc) then
        invalid=.true.
        zzzr=0.0
        return
      endif
      x=gdzz(nbuf(i))
      if (x .eq. -3) then
        idec=i
        ndec=ndec+1
        if (ndec .gt. 1) then
          invalid=.true.
          zzzr=0.0
          return
        endif
      elseif (x .lt. 0.0) then
        invalid=.true.
        zzzr=0.0
        return
      else
        value=x
      endif
c
 5    i=i+1
c
      if (i .gt. nc) then
        if (ndec .eq. 0) idec=nc
        zzzr=sign*value/10**(nc-idec)
        return
      endif
c
      x=gdzz(nbuf(i))
      if (x .eq. -3) then
        idec=i
        ndec=ndec+1
        if (ndec .gt. 1) then
          invalid=.true.
          zzzr=0.0
          return
        endif
      elseif (x .eq. -4.0 .and. i .eq. nc) then
        invalid=.true.
        zzzr=0.0
        return
      elseif (x .eq. -4.0 .and. i .lt. nc) then
        if (value .lt. 0.0) then
          invalid=.true.
          zzzr=0.0
          return
        endif
        if (ndec .eq. 0) idec=i-1
        value=value/10**(i-1-idec)
        go to 10
      elseif (x .lt. 0.0) then
        invalid=.true.
        zzzr=0.0
        return
      else
        if (value .lt. 0.0) value=0.0
        value=10.0*value+x
      endif
      go to 5
c
c     decode exponent
c
 10   xpon=-1.0
      i=i+1
      if (i .gt. nc) then
        invalid=.true.
        zzzr=0.0
        return
      endif
      se=1.0
      e=gdzz(nbuf(i))
      if (e .eq. -1) then
        se=1.0
        i=i+1
      elseif (e .eq. -2) then
        se=-1.0
        i=i+1
      elseif (e .lt. 0.0) then
        invalid=.true.
        zzzr=0.0
        return
      endif
      if (i .gt. nc) then
        invalid=.true.
        zzzr=0.0
        return
      endif
c
 15   continue
      e=gdzz(nbuf(i))
      if (e .ge. 0.0) then
        if (xpon .lt. 0.0) xpon=0.0
        xpon=xpon*10.0+e
      else
        invalid=.true.
        zzzr=0.0
        return
      endif
      if (i .lt. nc) then
        i=i+1
        go to 15
      else
        if (xpon .lt. 0.0) then
          invalid=.true.
          zzzr=0.0
          return
        else
          if (xpon .gt. 100) then
            invalid=.true.
            zzzr=0.0
            return
          endif
          zzzr=sign*value*10.0**(se*xpon)
          return
        endif
      endif
      invalid=.true.
      zzzr=0.0
      return
      end
c
c
c
      real function gdzz(ichar)
      integer i
      dimension idig(14)
      character*1 idig,ichar
      data idig/'e','.','-','+','0','1','2','3','4','5','6','7','8','9'/
      do 1 i=1,14
        if (ichar .eq. idig(i)) then
          gdzz=float(i) - 5.0
          return
        endif
 1    continue
      gdzz=-5
      return
      end
c
c
c
c
c
      logical function logic(ichar,invalid)
      character*8 ichar
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine tests the string ichar for a logical
c     value. valid strings are t,f,true,false,on,off,yes,no
c     invalid=.true. and logic=.false. if not a valid logical
c     term
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      logical invalid
      if (ichar .eq. 't') then
        logic=.true.
      elseif (ichar .eq. 'true') then
        logic=.true.
      elseif (ichar .eq. 'on') then
        logic=.true.
      elseif (ichar .eq. 'yes') then
        logic=.true.
      elseif (ichar .eq. 'f') then
        logic=.false.
      elseif (ichar .eq. 'false') then
        logic=.false.
      elseif (ichar .eq. 'off') then
        logic=.false.
      elseif (ichar .eq. 'no') then
        logic=.false.
      else
        invalid=.true.
        logic=.false.
        return
      endif
      invalid=.false.
      return
      end
c
c
c
      integer function ntger(xx,invalid)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this function tests the real value, xx, to see if it
c     can be converted to an integer without changing it's
c     value. If not, invalid is set to true and ntger is
c     the truncated value of xx
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer i
      real xx,diff,x
      logical invalid
      data diff/1.0e-13/
      x=xx
      i=ifix(x)
      if (x .ne. 0.0) then
        if (abs((x-float(i))/x) .gt. diff) then
          invalid=.true.
          ntger=i
          return
        endif
      else
        if (x .ne. float(i)) then
          invalid=.true.
          ntger=i
          return
        endif
      endif
      invalid=.false.
      ntger=i
      return
      end
