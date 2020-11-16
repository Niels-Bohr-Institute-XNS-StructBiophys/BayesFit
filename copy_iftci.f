c************************************************************************************
c             Download of source code iftci.f
c************************************************************************************
c
c   1) Compilation (linux):  gfortran -m32 -O2 iftci.f -o iftci
c                            gfortran -march=native -O3 source/iftci.f -o iftci
c
c   2) Run:                  iftci < inputfile.d                    
c
c   3) The file: inputfile.d has to contain the 14 lines:
c
c   line 1:  the name of the data file  - compulsory -
c   line 2:  value for q_min            or a blank line
c   line 3:  value for q_max            or a blank line
c   line 4:  value for d_max            or a blank line
c   line 5:  value for eta              or a blank line
c   line 6:  value for alpha            or a blank line
c   line 7:  value for smearing const.  or a blank line
c   line 8:  value for ratio            or a blank line
c   line 9:  value for method           or a blank line
c   line 10: value for no of points     or a blank line
c   line 11: value for no of extra calc or a blank line
c   line 12: value for transformation   or a blank line
c   line 13: value for background       or a blank line
c   line 14: value for screensize       or a blank line
c
c   NB NB The input values are the same as at the web site 
c   Only the first line has to be a non blank line
c   For further information use the ?-marks at the web-site
c
c   PLS ALSO NOTE: The graphics routine was written for the 
c   web site and will not work without several modifications.
c   The gfortran compiler may give problems with Windows -
c   it is probably better to use another compiler.
c
c   The program was tested with gfortran 4.4 and it may 
c   be necessary to include the "-static" flag. 
c************************************************************************************


C************************************************************************************
C  Short outline of the order of the functions and subroutines below:
C************************************************************************************
C  Program IFTc          calculates the optimal (p_1, p_excl, p_struct, alpha, eta, d, ratio) 
C  subroutine powell     calculates the optimal (alpha, eta, d, ratio) (used in IFTc)
C    +linmin             (used for powell)
C    +brent              (used for powell)
C    +f1dim              (used for powell)
C    +mnbrak             (used for powell)
C  function func         calculates the optimal (p_1, p_excl, p_struct) (used in powell)
C  subroutine prior      calculates the prior estimae of (p_1, p_excl, p_struct) (used in func)
C  subroutine excl       calculates the excluded volume distribution p_excl(r) (used in prior)
C  subroutine sphere     calculates the correlation function for a sphere (used in prior)
C  subroutine ellips     calculates the dist. distribution function for an ellipsoid
C  subroutine trans      construction of the Fourier transformation matrix (used in func)
C  subroutine svdcmp     singular value decompositon (used in func)
C  function ran1         random number generator for error estimation (used in IFTc)
C  subroutine plot       makes a plot file to be used by wgnuplot (used in IFTc)
C  subroutine alpha_e  calculates the evidence for several values of alpha (used in IFTc)
C*********************************************************************************
      PROGRAM IFTc
C********************************************************************************
C     Various parameters may be changed to improve the convergence.
C
C     XPREC (precision of FUNC) 
C     OMEGA (steplenght in FUNC) 
C     MAXIT (maximum number of iterations in FUNC)
C     RPREC (precision of POWELL)         
C     MAXCPU (maximum cputime allowed for POWELL)
C
C     The constraints for S(q) may be changed to obtain S(q) -> 1 for large q.
C     Search for the line: "Change the constraints for S(q)"
C********************************************************************************
      PARAMETER (NMAX=1000)
      PARAMETER (NDIST=1000)
      parameter (ndim=4)

      external func
      REAL A(NMAX,0:NMAX),B(0:NMAX,0:NMAX),Adec(NMAX,0:NMAX)
      real FT(NMAX,NMAX),SMEAR(NMAX,NMAX),xhmin(ndim)
      REAL X(NMAX),Y(NMAX),SD(NMAX),XF(0:NMAX),F(0:NMAX),u(nmax,nmax)
      REAL x_original(NMAX),sd_original(NMAX)
      REAL M(0:NMAX),YSUM(0:NMAX),FM(NMAX),w(nmax),dd(1000)
      real ddx(1000),al(1000),etaa(1000),p(ndim),xi2(ndim,ndim)
      real sigma(0:nmax),sigf(0:nmax),ntott(ndist),I0(ndist),int(1000)
      real ftot(0:nmax,ndist),prob(ndist),xmean(0:nmax),qmin,qmax
      real factor(100),fmdec(nmax)
      integer ptot,ir,clock,clockold
      CHARACTER*48 ANAME
      character*12 redaname
      chARACTER*50 BNAME,aaname,plname
      CHARACTER*53 hxname,gsname,gxname,sname
      CHARACTER*54 KNAME
      CHARACTER*55 ProbNAME,QNAME,in1name,in2name,dataname,diamname
      character*60 xin1name,xin2name
      CHARACTER*56 FNAME,GNAME,HNAME
      CHARACTER*6 CHAR
      CHARACTER*1 ANSWER,answer2,answer5,answer8,answer7,answerm
      character*30 dummy 
      character*80 string
      real pi

      common /data1/ nextra,dmax,eta,diam,nfactor,nadd
      common /data2/ c,cnst,ntot2,rlogdet
      common /data4/ mtot
      common /data3/ x,sd,y,xf,ntot,prob,ftot,nof,sigma
      common /data5/ ratio,ratioold
      common /data6/ alphaest,dest,etaest
      common /data7/ f,sigf
      common /data8/ nfit,nerror,answer,answer2
      common /data9/ niter
      common /data11/ itetot,itetotmax
      common /data15/ answer5
      common /data25/ csmear,cexp
      common /screen/ answer8,answer7
      common /back/ background
      common /data16/ xhmin,evidencesave,ngrid
      common /cputime/ clockold,cpumax,answerm
      common /test/ ndtest

C********************************************************************
c     Initialization
c********************************************************************
      PI=acos(-1.)
c max cpu time in seconds
      cpumax=900
      nof=0
      itetot=0
      itetotmax=8000000

      do 10 i=0,nmax
      do 11 j=1,ndist
      ftot(i,j)=0
  11  continue
  10  continue
c********************************************************************
c     Name of data file and input of data
c********************************************************************
      call system_clock(clock)
      clockold=clock
      open(unit=555,file='inputfile.d',status='unknown')
    1 WRITE(6,3)
    3 FORMAT(1X,'Write name of the inputfile - format (x,y,sd)   => ',$)
    4 FORMAT(A)
      read(555,4)aname

      if(aname.eq.'      ') aname='data.d'

c      bname='me'//ANAME
      bname='pr.d'

c      FNAME='trans_'//BNAME        
      fname='fit.d'
      GNAME='decon_'//BNAME        
      HNAME='PRIOR_'//BNAME   
      KNAME='out_'//BNAME   
      gsNAME='gs_'//BNAME   
      gxNAME='gx_'//BNAME   
c      hxname='in_'//bname
      hxname='data.d'
      sname='st_'//bname
      plname='pl'//aname
      diamname='diam_'//bname
      redaname=aname

      qmin=0.0001
      qmax=1000
      WRITE(6,131)
  131 FORMAT(1X,'qmin (min 0.0001) or <enter>                    => ',$)
      read(555,4)dummy
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)qmin
      close(50)
      endif
      if(qmin.le.0.0001) qmin=0.0001

      WRITE(6,132)
  132 FORMAT(1X,'qmax or <enter>                                 => ',$)
      read(555,4)dummy
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)qmax
      close(50)
      endif

      ndata0=1
      WRITE(6,1329)
 1329 FORMAT(1X,'Rebinning points -> number of points used       => ',$)
      nrebin=20000
      OPEN(1,FILE=aname,STATUS='UNKNOWN')
c************************************************************
c skipping text lines in data file and writing data to dummy
c*************************************************************
      open(111,file='dummy.d',status='unknown')
      open(112,file='datatot.d',status='unknown')
      ndata=0
      k=0
   2  if(ndata.ge.1) goto 8
  22  continue
      read(1,*,err=2,end=8)x1,y1,sd1
        k=k+1
        if(k.ge.100000) goto 8
        if((x1.gt.50).and.(k.le.10)) backspace 1
        if((x1.gt.50).and.(k.le.10)) goto 22
      if(sd1.le.0) ndata0=0
      write(112,*)x1,y1,sd1
      if((x1.lt.qmin).or.(x1.gt.qmax)) goto 22
      ndata=ndata+1
      write(111,*)x1,y1,sd1
        if(sd1.le.0) then 
        ndata=-2
        goto 6769
        endif
        if(abs(y1/sd1).gt.1.e6) then 
        ndata=-1
        goto 6769
        endif
      goto 22
   8  close(111)
      close(112)
      close(1)
 6769 write(6,*)'ndata = ',ndata,k,qmax,x1

cccccccccccccccccccccccccccc changed feb 2012
      if(ndata.le.1) then
      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig2.gif"'
       write(88,*)
     -'set label " Wrong input format or q-limitation " at -4.8,0.55'
       write(88,*)
     -'set label "    No points read from data file   " at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot5.pl',status='unknown')
 	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig5.gif"'
       write(88,*)
     -'set label " Wrong input format or q-limitation " at -4.8,0.55'
       write(88,*)
     -'set label "    No points read from data file   " at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
        close(88)
      write(6,*)'error'

      STOP
      endif
cccccccccccccccccccccccccccc changed feb 2012

cccccccccccccccccccccccccccc changed feb 2012
 1234 if(ndata0.eq.0) then
      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig2.gif"'
       write(88,*)
     -'set label "  Wrong input format in data file   " at -4.8,0.6'
       write(88,*)
     -'set label "   Zero standard deviation read     " at -4.8,0.5'
       write(88,*)
     -'set label "      Pls. check the data file      " at -4.5,0.4'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot5.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig5.gif"'
       write(88,*)
     -'set label "  Wrong input format in data file   " at -4.8,0.6'
       write(88,*)
     -'set label "   Zero standard deviation read     " at -4.8,0.5'
       write(88,*)
     -'set label "      Pls. check the data file      " at -4.5,0.4'       
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
        close(88)
      write(6,*)'error2'

      STOP
      endif
cccccccccccccccccccccccccccc changed feb 2012

 6767 continue


      if(ndata.lt.nrebin) then
      OPEN(1,FILE='dummy.d',STATUS='UNKNOWN')
      x1sum=x1
      goto 7777
      endif
c****************************************************************
c if ndata > xx reading data from dummy and rebinning into dummy2
c****************************************************************
      nr=ndata/nrebin
      nr2=ndata/nr
      nrest=ndata-nr2*nr

c     nr: rebin no.
c     nr2(+1): number of data points used (after rebinning)

      open(111,file='dummy.d',status='old')
      open(44,file='dummy2.d',status='unknown')

      if(nrest.eq.0) goto 1310
      x1sum=0
      y1sum=0
      sd1sum=0
      do 1311 i=1,nrest
      read(111,*)x1,y1,sd1
      x1sum=x1sum+x1
      y1sum=y1sum+y1
      sd1sum=sd1sum+sd1**2
 1311 continue
      x1sum=x1sum/nrest
      y1sum=y1sum/nrest
      sd1sum=sqrt(sd1sum)/nrest
      write(44,*)x1sum,y1sum,sd1sum

 1310 do 1313 j=1,nr2
      x1sum=0
      y1sum=0
      sd1sum=0
      do 1312 i=1,nr
      read(111,*)x1,y1,sd1
      x1sum=x1sum+x1
      y1sum=y1sum+y1
      sd1sum=sd1sum+sd1**2
 1312 continue
      x1sum=x1sum/nr
      y1sum=y1sum/nr
      sd1sum=sqrt(sd1sum)/nr
      write(44,*)x1sum,y1sum,sd1sum
 1313 continue

      close(44)
c****************************************************************
c Reading data from dummy3 using only points within [q_min;q_max]
c Adding points to constrain S(q) -> 1 for q > qmax/1.6
c The number 1.6 should probably be changed to 2.0 or even larger.
c The error of the added points may also be changed to ensure that
c S(q) -> 1 
c Change the constraints for S(q) here:
c****************************************************************
      OPEN(1,FILE='dummy2.d',STATUS='UNKNOWN')
 7777 close(111)
      mtot=0
      nadd=0
      qmax2=qmax
      if(qmax.eq.10) qmax2=x1sum
      if(qmax.ge.x1sum) qmax2=x1sum
      do 7 i=1,10000
    5 READ(1,*,end=798)Xi,Yi,SDi
        if((xi.ge.qmin).and.(xi.le.qmax)) then
        mtot=mtot+1
        x(mtot)=xi
        x_original(i)=xi
        sd_original(i)=sdi
c****************************************************************
        Y(mtot)=yi-0.0
        sd(mtot)=sdi

Change this line???
        if(xi.ge.qmax2/1.6) nadd=nadd+1
        endif
    7 continue
  798 CONTINUE
      close(1)

      rangeq=x(mtot)-x(1)
      xmtotxx=x(mtot)
      x1xx=x(1)
      write(6,*)'ndata2 = ',ndata,rangeq,x(mtot),x(1)
c********************************************************************
c     End of input of data
c********************************************************************
      mtotxx=mtot
c********************************************************************
c     Start input of parameters
c******************************************************************** 
      WRITE(6,1335)
 1335 FORMAT(1X,'Estimate of maximum diameter d or <enter>       => ',$)
c      WRITE(6,*)'(A negative value entered will fix the diameter)'
      ndtest=0
      read(555,4)dummy
      dest=0.5*3.14159/x(1)
      write(6,*)'dest = ',dest
      write(6,*)'x(mtot)/x(1) = ',x(mtot)/x(1)
      if(x(mtot)/x(1).gt.30) destnew=0.5*3.14159/(x(mtot)/30.) 
      write(6,*)'destnew = ',destnew
      xdest=1.0
      if(dummy.ne.'       ') then
        open(unit=50,file='dummy.d',status='unknown')
          if(dummy(1:1).eq.'f') dummy(1:1)='F'
          if(dummy(1:1).eq.'F') then
          xdest=0.0000001
c          xdest=0
          ndest=1.0
          write(50,4)dummy(2:)
          endif
          if(dummy(1:1).ne.'F') write(50,4)dummy
        close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)dest
      ndtest=1.0
      close(50)
      endif

  
      WRITE(6,1321)
 1321 FORMAT(1X,'Estimate of vol frac (start max 0.1) or <enter> => ',$)
      read(555,4)dummy
      etaest=0
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)etaest
      close(50)
      endif

      WRITE(6,1322)
 1322 FORMAT(1X,'Estimate Lagr. multiplier log(alpha) or <enter> => ',$)
      read(555,4)dummy
      alphaest=-100.
      xalphaest=1.0
      if(dummy.ne.'       ') then
        open(unit=50,file='dummy.d',status='unknown')
          if(dummy(1:1).eq.'f') dummy(1:1)='F'
          if(dummy(1:1).eq.'F') then
          xalphaest=0
          nalphaest=1.0
          write(50,4)dummy(2:)
          endif
          if(dummy(1:1).ne.'F') write(50,4)dummy
        close(50)
        open(unit=50,file='dummy.d',status='unknown')
        read(50,*)alphaest
        close(50)
      endif
     

      WRITE(6,1139)
 1139 FORMAT(1X,'Enter value for c_smear     => ',$)
      read(555,4)dummy
      csmear=0
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)csmear
      close(50)
      cexp=2
        if(csmear.lt.0) then
        csmear=-csmear
        cexp=6.
        endif
      endif


      WRITE(6,135)
  135 FORMAT(1X,'Estimate of axial ratio or <enter>              => ',$)
      read(555,4)dummy
      ratio=1.0
      nratio=1
      xratio=1.
      if(dummy.ne.'       ') then

        open(unit=50,file='dummy.d',status='unknown')
          if(dummy(1:1).eq.'f') dummy(1:1)='F'
          if(dummy(1:1).eq.'F') then
          nratio=0
          write(50,4)dummy(2:)
          endif
          if(dummy(1:1).ne.'F') write(50,4)dummy
        close(50)

c      open(unit=50,file='dummy.d',status='unknown')
c      write(50,4)dummy
c      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)ratio
      close(50)
      endif

      WRITE(6,136)
  136 FORMAT(1X,'Fit ratio by evidence (E) moments (M) or no (N) => ',$)
      read(555,4)dummy
      answer='M'
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,4)answer
      close(50)
      endif


      if(answer.eq.'e') answer='E'
      if((etaest.eq.0).and.(answer.eq.'E')) answer='M'
c     ((the axial ratio cannot be fitted by evidence without an excluded volume ddf.))

      if(answer.eq.'n') answer='N'
      if(answer.eq.'m') answer='M'
c 
c     Fitting alpha(1),d_max(2),eta(3),ratio(4)
c
c     ratio(4) is only fitted when eta > 0 i.e.
c
c     for eta=0 only alpha(1),d_max(2) are fitted.
c
c     UNLESS answerm='M' then alpha,d_max AND ratio are fitted.
c 
      if(answer.eq.'N') nfit=3
      if(answer.eq.'M') nfit=3
      if(answer.eq.'E') nfit=4

      WRITE(6,139)
  139 FORMAT(1X,'Number of points used for p(r)    (50-100)     => ',$)
      read(555,4)dummy
      ntot=50
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)ntot
      close(50)
      endif
      if(ntot.gt.100) ntot=100
      if(ntot.lt.20)   ntot=20

      if((etaest.gt.0).and.(ntot.gt.80)) cpumax=900

c      WRITE(6,133)
c  133 FORMAT(1X,'Enter relative precision         (0.1-0.001)    => ',$)
c      read(555,4)dummy
c      rprec=0.05
      rprec=0.1
c      rprec=0.0001
c      if(dummy.ne.'       ') then
c      open(unit=50,file='dummy.d',status='unknown')
c      write(50,4)dummy
c      close(50)
c      open(unit=50,file='dummy.d',status='unknown')
c      read(50,*)rprec
c      close(50)
c      endif


      WRITE(6,134)
  134 FORMAT(1X,'Number of extra solutions for calc. of errors   => ',$)
      read(555,4)dummy
      answer2='Y'
      nerror=2
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,*)nerror
      close(50)
      if(nerror.eq.0) answer2='N'
      if(nerror.ne.0) answer2='Y'
      endif
      if(nerror.gt.1000) nerror=1000
      if(nerror.lt.2) nerror=2


      WRITE(6,138)
  138 FORMAT(1X,'Transformation 
     -[D]ebye (def), [B]essel, [C]osine, [S]ize Neg.Debye[N] => ',$)
      read(555,4)dummy
      answer5='D'
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,4)answer5
      close(50)
      endif

      answerm='O'
      if(answer5.eq.'s') answer5='S'
      if(answer5.eq.'c') answer5='C'
      if(answer5.eq.'d') answer5='D'
      if(answer5.eq.'b') answer5='B'
      if(answer5.eq.'n') answer5='N'
      if(answer5.eq.'k') answer5='K'
      if(answer5.eq.'m') answer5='M'
      if(answer5.eq.'q') answer5='Q'
      if(answer5.eq.'M') then
      answer5='D'
      answerm='M'
      endif
      if(answer5.eq.'Q') then
      answer5='D'
      answerm='Q'
      endif
      if(nfit.eq.3.and.answerm.eq.'M') rprec=0.0001
      if(nfit.eq.3.and.answerm.eq.'Q') rprec=0.0001


      if(etaest.ne.0.and.answerm.eq.'M') then
      write(6,*)'error3'
      stop
      endif
      if(etaest.eq.0.and.answerm.eq.'M') then
      cpumax=999
      answer='E'
c ratio is fitted by evidence
      nfit=4
c 1 is subtracted from nfit below
      endif

      WRITE(6,1138)
 1138 FORMAT(1X,'Fitting flat background [N]o (default) or [Y]es => ',$)
      read(555,4)dummy
      answer7='N'
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,4)answer7
      close(50)
      endif
      if(answer7.eq.'y') answer7='Y'
      if(answer7.eq.'n') answer7='N'
      if(answer7.eq.'Y') sigf(0)=10
      if(answer7.eq.'N') sigf(0)=0
      jmin=1
      if(answer7.eq.'Y') jmin=0


      WRITE(6,1388)
 1388 FORMAT(1X,'Screenplot (S)mall (default) or (L)arge         => ',$)
      read(555,4)dummy
      answer8='S'
      if(dummy.ne.'       ') then
      open(unit=50,file='dummy.d',status='unknown')
      write(50,4)dummy
      close(50)
      open(unit=50,file='dummy.d',status='unknown')
      read(50,4)answer8
      close(50)
      endif

      if(answer8.eq.'l') answer8='L'
      if(answer8.eq.'s') answer8='S'      
c********************************************************************
c     End input of parameters
c******************************************************************** 

c********************************************************************
c     Addition of points for structure factor
c******************************************************************** 

      if(etaest.eq.0) nadd=0
      if(etaest.eq.0) nextra=0

      write(6,*)'Number of data points in data file = ',ndata
      write(6,*)'Number of data points after rebin  = ',mtot
      write(6,*)'Number of points added (S(q)->1)   = ',nadd

      mtot=mtot+nadd
      do 210 i=mtot,mtot-nadd+1,-1
      x(i)=x(i-nadd)
      y(i)=0

c     Change this line???
      sd(i)=1.0*sd(i-nadd)
 210  continue

      OPEN(24,FILE=HxNAME,STATUS='UNKNOWN')
      DO 201 I=1,MTOT-nadd
      WRITE(24,*)X(I),Y(I),sd(i)
  201 CONTINUE
      close(24)

c     for some reason we square sigma to get the variance here...
      DO 600 I=1,MTOT
      SD(I)=1.*SD(I)**2
  600 CONTINUE


c********************************************************************
c     End addition of points for structure factor
c******************************************************************** 

      OPEN(36,FILE=KNAME,STATUS='unknown',access='append')
      nextra=ntot
      write(36,6670)'-----------------------------------------
     ------------------------------------'
      write(36,6668)'Input values:   ',qmin,qmax,dest,
     -etaest,alphaest,ratio,answer,ntot,rprec,nerror
 6668 format(a,f6.4,x,f6.4,x,f6.2,x,f6.4,x,f7.2,x,f6.3,x,a,x,i3,
     -x,f7.4,x,i3)
      write(36,*)
 6669 write(6,*)'                                          '
      write(6,*)'                                          '
      if(answer7.eq.'Y') then
      write(6,*)' ite   diam  log(alpha)    chi    axratio    check  
     -backg   -evidence'
      else
      write(6,*)' ite   diam  log(alpha)    chi    axratio    check  
     - eta   -evidence'
      endif
      write(6,*)
 6670 format(a)
c********************************************************************
c     Calculation of best set of hyperparameters -
c     uses the powell algorithm from Numerical Recipes
c********************************************************************
      do 6661 i=1,ndim
      do 6662 j=1,ndim
      xi2(i,j)=0
 6662 continue
 6661 continue

      if(etaest.eq.0) nfit=nfit-1

      write(6,*)
      write(6,*) 'nfit = ',nfit
      write(6,*)
      if(nfit.eq.4) then
      p(1)=etaest
      p(2)=dest
      p(3)=ratio
      p(4)=alphaest
      ndiam=2
      nalpha=4
     
      xi2(1,1)=etaest*0.1
      xi2(2,2)=dest*0.05*xdest
      xi2(3,3)=ratio*0.1
      xi2(4,4)=0.5*xalphaest
      endif
c maxent
       if(nfit.eq.3) then
       write(6,*)' maxent *************'
        if(etaest.eq.0) then
        p(1)=dest
        p(2)=ratio
        p(3)=alphaest
      ndiam=1
      nalpha=3

        xi2(1,1)=dest*0.05*xdest
        xi2(2,2)=ratio*0.1
        xi2(3,3)=0.1*xalphaest
        endif
      
        if(etaest.ne.0) then
        p(1)=etaest
        p(2)=dest
        p(3)=alphaest
        xi2(1,1)=etaest*0.1
        xi2(2,2)=dest*0.05*xdest
        xi2(3,3)=0.1*xalphaest
        ndiam=2
        nalpha=3

        endif
      endif

      if(nfit.eq.2) then
        p(1)=dest
        p(2)=alphaest 
        xi2(1,1)=dest*0.05*xdest
        xi2(2,2)=0.1*xalphaest
      ndiam=1
      nalpha=2
      endif

      nof=1
      ftol=rprec
      ngrid=1
      if(alphaest.eq.-100) then
      call alpha_e(p,ndtest)
      if(answerm.eq.'M') p(3)=p(3)+1.
      endif

c      goto 9900
      ndtest=2
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(nratio.eq.0)   then 
         write(6,*)'******************************'       
         xi2(2,2)=ratio*0.000
         xratio=0.
         goto 9900
         endif
      if(answerm.eq.'M'.and.ratio.eq.1.0.and.nfit.eq.3) then
      write(6,*)'--------'
      evidence=func(p)
      rg=0
      sum=0
      sum1=0
      sum3=0
      do 5139 j=1,ntot
      rg=rg+xf(j)**2*abs(f(j))
      sum1=sum1+xf(j)*abs(f(j))
      sum3=sum3+xf(j)**3*abs(f(j))
      sum=sum+f(j)
 5139 continue
      rg=sqrt(abs(0.5*rg/sum))
      t=(sum3/sum)/(sum1/sum)**3

c     prolate:
c      if(ratio.ge.1.00) then
      if(t.le.1.4003198) t=1.40031981
      if(t.le.1.617)
     -ax=1+1.318*(t-1.4003198)**0.47+1.633*(t-1.4005)
      if(t.ge.1.617)
     -ax=1+1.603*(t-1.4)**4.594+3.1702*(t-1.30244)
c      endif

c     oblate
c      if(ratio.le.1.00) then
      if(t.le.1.4003198) t=1.40031981
      ax2=1-0.73649*(t-1.4)**0.190812-1.6129*(t-1.5017)
      if(t.ge.1.697) ax2=0.1
c      endif
       write(6,*)'Ratio estimated from p(r) = ',ax,ax2
c       if(ax.le.1.2) ax=1.0
       if(ax.ge.2) ax=1.8
       ax=0.5+0.5*ax
       if(ax.ge.1.1) ax=1.1
       if(ax.le.0.9) ax=0.9
       ratio=ax
       p(2)=ratio
       xi2(2,2)=ratio*0.1
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 9900 ngrid=0
      call powell(p,xi2,nfit,ndim,ftol,iter,fret)
      evidence=func(p)
      if(itetot.ge.itetotmax) then
      write(6,*)
      write(6,*)
      write(6,*) ' NO CONVERGENCE FOR MAXIMUM ITERATIONS = ',	itetotmax
      write(6,*)
      write(6,*)
      endif

      write(6,*)'                                          '
      if(answer7.eq.'Y') then
      write(6,*)' ite   diam  log(alpha)    chi    axratio    check  
     -backg   -evidence'
      else
      write(6,*)' ite   diam  log(alpha)    chi    axratio    check  
     - eta   -evidence'
      endif
      write(6,*)'                                          '
 6667 nofmax=nof-1
      nerrorold=nofmax
      if(answer2.eq.'N') goto 6634
      if(answer2.eq.'N') write(6,*)
c********************************************************************
c
c     Calculation of errors via extra distributions via Monte Carlo method
c
c     NB - not really suited for the web version due to the time limitation !!
c           (too few points found within the time limit)
c
c********************************************************************      
      write(6,*)' Calculating errors in',nerror,' points'
      write(6,*)' using',nerror-nofmax, ' points'
      niter=niter+5
      dum=func(p)
      dumold=dum
      if(nofmax+nerror.gt.990) nerror=990-nofmax
      nofmax=nerror+nofmax
      nerrorold=nofmax+10
c optimal values of parameters saved:
      sp1=p(1)
      sp2=p(2)
      sp3=p(3)
      sp4=p(4)
      spn=p(nfit)
      IR=-eta*2**30
      fac=0.03
 

      p(1)=(1+fac)*sp1
      factor(1)=fac/(abs(func(p)-dumold)+0.0001)
      p(1)=sp1

      p(2)=(1+fac)*sp2
      factor(2)=fac/abs(func(p)-dumold+0.00001)
      p(2)=sp2

      p(3)=(1+fac)*sp3
      factor(3)=fac/abs(func(p)-dumold+0.000001)
      p(3)=sp3

      p(4)=(1+fac)*sp4
      factor(4)=fac/abs(func(p)-dumold+0.000001)
      p(4)=sp4

c     p(nfit)=(1+fac)*spn
c     factor(5)=fac/abs(func(p)-dumold+0.000001)
c     p(nfit)=spn

c try this instead - above is unstable...


      factor(1)=0.1
      factor(2)=0.1
      factor(3)=0.1
      factor(4)=0.1

      if(xdest.le.0.5) factor(ndiam)=0
      if(xalphaest.le.0.5) factor(nalpha)=0
      if(xratio.le.0.5) factor(2)=0

      xmult=3.0
      do 6333 i=1,nerror
      p(1)=sp1+factor(1)*sp1*(ran1(ir)-0.5)*xmult
      p(2)=sp2+factor(2)*sp2*(ran1(ir)-0.5)*xmult
      p(3)=sp3+factor(3)*sp3*(ran1(ir)-0.5)*xmult
      p(4)=sp4+factor(4)*sp4*(ran1(ir)-0.5)*xmult
      p(nfit)=spn+factor(5)*spn*(ran1(ir)-0.5)*xmult
      dum=func(p)
      if(dum.ge.9.e8) then
      xmult=xmult*0.5
      write(6,*)'xmult reduced to ',xmult
c      if(xmult.le.0.001) call error
      if(xmult.le.0.001) then
      write(6,*)
      write(6,*)'Error in calculation of errors!'
      write(6,*) goto 6634
      endif
      endif
 6333 continue     
c********************************************************************
c     Calculation of probabilities
c********************************************************************
 6634 evimax=-1.e20
      do 9991 nof=1,nofmax
      if(prob(nof).gt.evimax) then
      evimax=prob(nof)
c 3 lines added
      cmin=ftot(nmax-2,nof)
      if(cmin.le.0.01) cmin=0.01
      if(cmin.ge.100) cmin=100
c      cmin=1.
      endif
 9991 continue

      do 9992 nof=1,nofmax
c      prob(nof)=exp(prob(nof)-evimax)
      prob(nof)=(exp(prob(nof)-evimax))**(1./cmin)
 9992 continue

      totprob=0.
      do 9118 nof=1,nofmax
      totprob=totprob+prob(nof)
 9118 continue

      do 9119 nof=1,nofmax
      prob(nof)=prob(nof)/totprob
 9119 continue
c********************************************************************
c     Calculation of the set of solutions
c********************************************************************
      dmax=0
      probx=0
      do 1002 nof=1,nofmax
      if(prob(nof).ge.probx) then
      probx=prob(nof)
      dmax=ftot(nmax-1,nof)
      endif
 1002 continue
      dmax=1.2*dmax
c
c  maximum diameter of average is set at 1.2*diam(max prob)
c
      ntotold=ntot
c ptot has to be even
      ptot=(ntot*1.2)*0.5
      ptot=ptot*2

      do 1011 j=1,ptot
      xf(j)=(1.0*j*dmax)/ptot
 1011 continue
      jpmax=ptot*0.5

      do 1012 j=1,jpmax
      xf(j+ptot)=(6.*j*dmax)/ptot
      j2=j+1.5*ptot
      xf(j2)=(6.*j*dmax)/ptot
 1012 continue

ccccccccc------------------------------ 
      ri0sum=0
      ri0sum1=0
      do 1010 nof=1,nofmax
      ri0=ftot(0,nof)*0.1*ftot(nmax-1,nof)/ntot
      ri0sum=ri0sum+ri0
      ri0sum1=ri0sum1+ri0*prob(nof)

      do 1013 j=1,ptot*2
      f(j)=0
 1013 continue   

      dfx=ftot(nmax-1,nof)/ntot
      dfx6=6*dfx

      backgroundnof=ftot(0,nof)*0.1*ftot(nmax-1,nof)/ntot
      ftot(0,nof)=backgroundnof 

c**** p_1---------------------------
c changed
c      do 1008 k=1,ntot
      do 1008 k=1,ntot+1
      do 1009 j=1,ptot     
      if((xf(j).le.dfx*k).and.(xf(j).gt.dfx*(k-1))) then
      frac=(xf(j)-dfx*(k-1))/dfx
      if(k.eq.1) f(j)=frac*ftot(k,nof)
      if(k.gt.1) f(j)=frac*ftot(k,nof)+(1-frac)*ftot(k-1,nof)
      if(k.eq.ntot+1) f(j)=(1-frac)*ftot(ntot,nof)
      endif
 1009 continue
 1008 continue
      
c**** p_excl--------------------------

      nkmin=ntot+1
      nkmax=1.5*ntot
      njmin=ptot+1
      njmax=1.5*ptot
      do 1028 k=nkmin,nkmax
      do 1029 j=njmin,njmax
      if((xf(j).le.dfx6*(k-ntot)).and.(xf(j).gt.dfx6*(k-ntot-1))) then
      frac=(xf(j)-dfx6*(k-ntot-1))/dfx6
      if(k.eq.ntot+1) f(j)=frac*ftot(k,nof)
      if(k.gt.ntot+1) f(j)=frac*ftot(k,nof)+(1-frac)*ftot(k-1,nof)
      endif
 1029 continue
 1028 continue
      
c**** p_struct------------------------

      nkmin=1.5*ntot+1
      nkmax=2*ntot
      njmin=1.5*ptot+1
      njmax=2*ptot
      do 1048 k=nkmin,nkmax
      do 1049 j=njmin,njmax
      j2=1.5*ntot
      if((xf(j).le.dfx6*(k-j2)).and.
     -(xf(j).gt.dfx6*(k-j2-1))) then
      frac=(xf(j)-dfx6*(k-j2-1))/dfx6
      if(k.eq.j2+1) f(j)=frac*ftot(k,nof)
      if(k.gt.j2+1) f(j)=frac*ftot(k,nof)+(1-frac)*ftot(k-1,nof)
      endif
 1049 continue
 1048 continue
   
c***** store ---------------------------

      do 1145 j=1,2*ptot
      ftot(j,nof)=f(j)
 1145 continue

 1010 continue
C*****************************************************************************
C     Output estimate of average p_1, p_excl, p_struct - including errors
C*****************************************************************************
      do 1116 j=0,2*ptot
      f(j)=0
      do 1117 nof=1,nofmax
      f(j)=f(j)+ftot(j,nof)*prob(nof)
 1117 continue
 1116 continue

      do 9651 j=0,2*ptot
      sigf(j)=0
 9651 continue

      do 9195 j=0,2*ptot
      do 9321 nof=1,nofmax
      sigf(j)=sigf(j)+prob(nof)*(ftot(j,nof)-f(j))**2
 9321 continue
 9195 continue
c--------------------------------------------
      OPEN(2,FILE=gxNAME,STATUS='UNKNOWN')
      write(2,*)0,0,0
      jpmax=1.5*ptot
      do 3755 j=ptot+1,jpmax
      write(2,*)xf(j),f(j),sqrt(abs(sigf(j)))
 3755 continue
      close(2)

      OPEN(2,FILE=gsNAME,STATUS='UNKNOWN')
      write(2,*)0,0,0
      jpmin=1.5*ptot+1
      do 3456 j=jpmin,2*ptot
      write(2,*)xf(j),f(j),sqrt(abs(sigf(j)))
 3456 continue
      close(2)

      OPEN(2,FILE=BNAME,STATUS='UNKNOWN')
      if(answer5.ne.'C') write(2,*)0,0,0
      sumf=0
      do 3478 j=1,ptot
      if(answer5.eq.'K') then
      write(2,*) xf(j),f(j)*xf(j)**1,sqrt(abs(sigf(j)))*xf(j)**1
      else
      write(2,*) xf(j),f(j),sqrt(abs(sigf(j)))
      sumf=sumf+f(j)
      endif
 3478 continue
      sumf=sumf*(xf(3)-xf(2))
      close(2)
c**********************************************************     
c     Calculation of I(0)  Rg  Eta  for each estimate (nof)
c**********************************************************
      ntot=ptot
      nextra=ntot
      dfx=dmax/ntot
      rmax1=0
      do 3130 nof=1,nofmax
      I0(nof)=0
      DO 8153 J=0,ptot
      I0(nof)=I0(nof)+Ftot(J,nof)*dfx
 8153 CONTINUE

      rg=0
      sum=0
      sum1=0
      sum3=0
      do 3139 j=1,ptot
      rg=rg+xf(j)**2*ftot(j,nof)
      sum1=sum1+xf(j)*ftot(j,nof)
      sum3=sum3+xf(j)**3*ftot(j,nof)
      sum=sum+ftot(j,nof)
 3139 continue
      rg=sqrt(abs(0.5*rg/sum))
      t=(sum3/sum)/(sum1/sum)**3

c     prolate:
c      if(ratio.ge.1.00) then
      if(t.le.1.4003198) t=1.40031981
      if(t.le.1.617)
     -ax=1+1.318*(t-1.4003198)**0.47+1.633*(t-1.4005)
      if(t.ge.1.617)
     -ax=1+1.603*(t-1.4)**4.594+3.1702*(t-1.30244)
c      endif

c     oblate
c      if(ratio.le.1.00) then
      if(t.le.1.4003198) t=1.40031981
      ax2=1-0.73649*(t-1.4)**0.190812-1.6129*(t-1.5017)
      if(t.ge.1.697) ax2=0.1
c      endif

      sump1=0
      do 8003 j=1,ptot
      sump1=sump1+ftot(j,nof)*(xf(2)-xf(1))
 8003 continue

      sumpx=0
      jpmax=1.5*ptot
      do 8004 j=ptot+1,jpmax
      sumpx=sumpx+ftot(j,nof)*(xf(ptot+2)-xf(ptot+1))
 8004 continue

      ftot(nmax-5,nof)=ax
      ftot(nmax-13,nof)=ax2
      ftot(nmax-8,nof)=-(1./8)*sumpx/sump1
      ftot(nmax-9,nof)=rg
      ftot(nmax-10,nof)=1
 3130 continue
c***********************************************************************
c     Output fit of data for average estimate
c***********************************************************************

C**** addition of extra points to q=0

      nn=20
      dr=1.0*x(1)/(nn+0.1)
      mtot=mtot+nn
      do 8150 i=mtot,nn+1,-1
      x(i)=x(i-nn)
 8150 continue
      do 8151 i=1,nn
      x(i)=x(nn+1)-dr*(nn-i+1)
 8151 continue
      if(x(1).le.0) x(1)=1.e-6
 8152 continue

C**** end addition

c     f(0) is changed from background to f(0) to be used with matrix     
      f(0)=f(0)*10/(xf(2)-xf(1))

      CALL TRANS(A,FT,SMEAR,X,XF,ntot,DFX,NMAX)
      csmear=0.0
      CALL TRANS(Adec,FT,SMEAR,X,XF,ntot,DFX,NMAX)

      OPEN(213,FILE=SNAME,STATUS='UNKNOWN')
 8142 OPEN(11,FILE=FNAME,STATUS='UNKNOWN')
     
      DO 8141 I=1,MTOT-nadd
          
          FM(I)=0
          fmdec(i)=0
          DO 8140 J=0,2*ptot
          fm(i)=FM(I)+A(I,J)*F(J)
          fmdec(i)=FMdec(I)+Adec(I,J)*F(J)
 8140     CONTINUE

          fmi=0
          DO 8143 J=0,ptot
          FMI=FMI+A(I,J)*F(J)

 8143     CONTINUE

c        WRITE(11,*)X(I),fmdec(i),fm(i),fmi
        WRITE(11,*)X(I),fmdec(i),fm(i)
        write(213,*)x(i),fm(i)/fmi,0
 8141   CONTINUE
      close(11)
      close(213)
c******************************************************
c     Output of parameters
c******************************************************
      dmax=0
      cav=0
      fm1=0
      etaav=0
      eta=0
      alp=0
      evi=0
      s=0
      rg=0
      ax=0
      ax2=0
      dot=0
      rtio=0
      sumng=0
      do 1132 nof=1,nofmax
      fm1=fm1+i0(nof)*prob(nof)
      dmax=dmax+ftot(nmax-1,nof)*prob(nof)
      cav=cav+ftot(nmax-2,nof)*prob(nof)
      alp=alp+ftot(nmax-3,nof)*prob(nof)
      dot=dot+ftot(nmax-4,nof)*prob(nof)
      ax=ax+ftot(nmax-5,nof)*prob(nof)
      ax2=ax2+ftot(nmax-13,nof)*prob(nof)
      s=s+ftot(nmax-6,nof)*prob(nof)
      evi=evi+ftot(nmax-7,nof)*prob(nof)
      etaav=etaav+ftot(nmax-8,nof)*prob(nof)
      rg=rg+ftot(nmax-9,nof)*prob(nof)
      eta=eta+ftot(nmax-10,nof)*prob(nof)
      rtio=rtio+ftot(nmax-11,nof)*prob(nof)
      sumng=sumng+ftot(nmax-12,nof)*prob(nof)
      write(155,*)nof,ftot(nmax-12,nof),prob(nof)
 1132 continue

      sumsh=(dmax/3.14159)*rangeq
     
c     error estimates

      sddmax=0
      sdcav=0
      sdfm1=0
      sdetaav=0
      sds=0
      sdeta=0
      sdalp=0
      sdevi=0
      sdrg=0
      sdax=0
      sdax2=0
      srtio=0
      ssumng=0
      do 5132 nof=1,nofmax
      sdfm1=sdfm1+(fm1-i0(nof))**2*prob(nof)
      sddmax=sddmax+(dmax-ftot(nmax-1,nof))**2*prob(nof)
      sdcav=sdcav+(cav-ftot(nmax-2,nof))**2*prob(nof)
      sdalp=sdalp+(alp-ftot(nmax-3,nof))**2*prob(nof)
      sdax=sdax+(ax-ftot(nmax-5,nof))**2*prob(nof)
      sdax2=sdax2+(ax2-ftot(nmax-13,nof))**2*prob(nof)
      sds=sds+(s-ftot(nmax-6,nof))**2*prob(nof)
      sdevi=sdevi+(evi-ftot(nmax-7,nof))**2*prob(nof)
      sdetaav=sdetaav+(etaav-ftot(nmax-8,nof))**2*prob(nof)
      sdrg=sdrg+(rg-ftot(nmax-9,nof))**2*prob(nof)
      sdeta=sdeta+(eta-ftot(nmax-10,nof))**2*prob(nof)
      srtio=srtio+(rtio-ftot(nmax-11,nof))**2*prob(nof)
      ssumng=ssumng+(sumng-ftot(nmax-12,nof))**2*prob(nof)
 5132 continue

      sdmax=sqrt(sddmax)
      scav=sqrt(sdcav)
      sfm1=sqrt(sdfm1)
      setaav=sqrt(sdetaav)
      sax=sqrt(sdax)
      sax2=sqrt(sdax2)
      seta=sqrt(sdeta)
      sevi=sqrt(sdevi)
      salp=sqrt(sdalp)
      ss=sqrt(sds)
      srg=sqrt(sdrg)
      srtio=sqrt(srtio)
      ssumng=sqrt(ssumng)

c
c     output of prior for maxent
c 
      if(answerm.eq.'M'.or.answerm.eq.'Q') then
      dfx=dmax/ptot
      ratio=rtio
      CALL PRIOR(M,ptot,y,XF,DFX,NMAX)
      OPEN(313,FILE='prior.d',STATUS='UNKNOWN')
      summ=0
      write(313,*)0.,0.,0.
      do 3132 j=1,ptot
      summ=summ+m(j)
 3132 continue
      summ=summ*(xf(3)-xf(2))
      do 3133 j=1,ptot
      write(313,*)xf(j),m(j)*sumf/summ,m(j)
 3133 continue
      endif

c      call ellips(m,ptot,nmax,dmax,ratio,ppmin)
c      do 4556 J=1,ptot
c      write(313,*)xf(j),m(j)
c 4556 continue      
c      write(313,*)'### axial = ',ratio,dfx,ptot,dmax
c      close(313)


c      if(rg.ge.1000) srg=srg/1000
c      if(rg.ge.1000) rg=rg/1000
c      if(dmax.ge.1000) sdmax=sdmax/1000
c      if(dmax.ge.1000) dmax=dmax/1000

c     correction of excl.volume using cylinder w. hemispherical caps
      az=ax
      if(ax.lt.1) az=1./ax
      correction=0.125*(6*az**2+12*az-2)/(3*az-1)
      etaav=etaav/correction
      setaav=setaav/correction
c      write(6,*)'correction of eta = ',correction


c      if(answer5.eq.'D') 
      if(answer5.eq.'S') then
      ax=1 
      sax=0
      ax2=1 
      sax2=0
      rtio=1
      srtio=0
      endif

      if(answer5.eq.'C') then
      ax=0
      sax=0
      ax2=0
      sax2=0
      rtio=0
      srtio=0
      endif
 
      if(answer5.eq.'B') then
      ax=0
      sax=0
      ax2=0
      sax2=0
      rtio=0
      srtio=0
      endif

      sold=s
      ssold=ss
      background=0
      sbackground=0 
      write(6,*)
 
      if(f(0).eq.0) then
      write(36,*)' Rg  ax ratio  -K     chi^2  log(alpha)  eta   dmax 
     - I(0)  evidence'
      write(6,*)' Rg  ax ratio  -K     chi^2  log(alpha)  eta   dmax 
     - I(0)  evidence'
      endif
      if(f(0).ne.0) then
      write(36,*)' Rg  ax ratio backg   Chi  log(alpha)   eta      d 
     - I(0)  evidence'
      write(6,*)' Rg  ax ratio  backg   Chi  log(alpha)   eta      d 
     - I(0)  evidence'
      S=f(0)*(xf(2)-xf(1))*0.1
      background=s
      sS=sigf(0)*(xf(2)-xf(1))*0.1
      sbackground=ss
      endif
      write(6,2135)rg,ax,-S,CAV,alp,etaav,dmax,fm1,-evimax
      write(6,2134)srg,sax,sS,sCAV,salp,setaav,sdmax,sfm1,sevi
      write(6,*)
      write(36,2135)rg,ax,-S,CAV,alp,eTAav,dmax,fm1,-evimax
      write(36,2134)srg,sax,sS,sCAV,salp,seTAav,sdmax,sfm1
     -,sevi
 2133 format
     -(f5.2,1x,F6.2,1X,f8.5,1x,
     -f6.3,1x,f8.3,1x,f6.3,1x,f6.2,1x,f7.3,1X,f9.2)
 2134 format
     -(f6.2,1x,F5.2,1X,e8.2,1x,f7.3,1x,
     -f8.3,1x,f6.3,1x,f6.2,1x,e8.3,1X,f9.2)
 2135 format
     -(f6.2,1x,F5.2,1X,e8.2,1x,f7.3,1x,
     -f8.3,1x,f6.3,1x,f6.2,1x,e8.3,1X,f9.2)
c*********************************************
c     Calculate p-value
c*********************************************
      DoF=mtotxx-sumng-2. 
      chi2=cav*mtotxx
      pval1=gammp(DoF*.5,chi2*.5)
      pval2=gammq(DoF*.5,chi2*.5)
      if(pval1.le.pval2)then
        pval=2.*pval1
      else
        pval=2.*pval2
      endif
c*********************************************
c     Correct errors
c*********************************************
      chi2r=chi2/(mtotxx-sumng-2.)
      beta=sqrt(chi2r)
c      REAL sd_corr(0:mtotxx)
      
      
      OPEN(21,FILE='rescale.d',STATUS='UNKNOWN')
      do 111,i=1,mtotxx
c        write(6,*)x_original(i),y(i),sd_original(i)*beta
        write(21,*)x_original(i),y(i),sd_original(i)*beta
111   continue
c*********************************************
c     Start OUTPUT for download for p(r) only
c*********************************************
      open(unit=166,file='parameters.d',status='unknown')
      if((answer5.eq.'D').or.(answer5.eq.'N')) then
      write(166,954)'Input file name  = ',aname
  954 format(x,a,x,a)
      write(166,953)'qmin and qmax    =  ', x1xx,xmtotxx
  953 format(x,a,f8.4,x,f8.4)
      write(166,*)'Number of data points               = ',mtotxx
      write(166,*)'Number of points in p(r)            = ',ntot
      write(166,*)'Number of error calculations        = ',nerrorold    
      write(166,*)'Errors are too low when few error calc. are used'
      write(166,*)
      if((dot.ge.0.9).and.(dot.le.0.95))
     -write(166,*)'Convergence of algorithm almost OK...'
      if(dot.gt.0.95)
     -write(166,*)'Convergence of algorithm OK'
      if(dot.lt.0.9) then
      write(166,*)'Convergence of algorithm NOT OK'
      write(166,*)'Try fewer points in p(r)'
      write(166,*)'and/or provide estimate of dmax and/or log(alpha)'
      endif

      write(166,*)
      xlimit=log(abs(fm1))
      xdlimit=log(abs(dmax))
      xrlimit=log(abs(rg))
      if((xlimit.le.11.52).and.(xlimit.ge.0)) then
      write(166,956)fm1,sfm1
  956 format(1x,'I(0) estimated             : ',f16.2,
     -'  +- ',f8.2,' ')
      else
      write(166,986)fm1,sfm1
  986 format(1x,'I(0) estimated             : ',f16.5,
     -'  +- ',f10.6,' ')
      endif
      if((xdlimit.le.11.52).and.(xdlimit.ge.0)) then
      write(166,958)dmax,sdmax
  958 format(1x,'Maximum diameter           : ',f16.2,
     -'  +- ',f10.2,' ')
      else
      write(166,1958)dmax,sdmax
 1958 format(1x,'Maximum diameter           : ',f16.2,
     -'  +- ',e10.2,' ')
      endif
      if((xrlimit.le.11.52).and.(xrlimit.ge.0)) then
      write(166,957)rg,srg
      else
      write(166,9957)rg,srg
      endif
  957 format(1x,'Radius of gyration         : ',f16.2,
     -'  +- ',f10.2,' ')
 9957 format(1x,'Radius of gyration         : ',f16.2,
     -'  +- ',e10.2,' ')

      write(166,1959)ax,sax
 1959 format(1x,'Axial ratio from p(r) (pro): ',f16.2,
     -'  +- ',f10.2,' ')

      write(166,1960)ax2,sax2
 1960 format(1x,'Axial ratio from p(r) (obl): ',f16.2,
     -'  +- ',f10.2,' ')

      if(answerm.eq.'M') then
      write(166,1965)rtio,srtio
 1965 format(1x,'Axial ratio from prior     : ',f16.2,
     -'  +- ',f10.2,' ')
      endif

      write(166,961)chi2r,scav*chi2r_c
  961 format(1x,'Reduced Chi-square         : ',f16.2,
     -'  +- ',f10.2,' ')
      
      if(background.ne.0) then
      write(166,962)background
  962 format(1x,'Background estimated       : ',f16.7)
      endif
      
      if(background.eq.0) then
      nbackground=background
      write(166,1962)nbackground
 1962 format(1x,'Background estimated       : ',13x,i9)
      endif
      
      if(etaav.ge.0.0001) then
      write(166,963)etaav,setaav
  963 format(1x,'Volume fraction            : ',f16.2,
     -'  +- ',f10.2,' ')
      endif
      write(166,964)alp,salp
  964 format(1x,'Log(alpha) (smoothness)    : ',f16.2,
     -'  +- ',f10.2,' ')
      write(166,976)sumng,ssumng
  976 format(1x,'Number of good parameters  : ',f16.2,
     -'  +- ',f10.2,' ')
      write(166,978)sumsh
  978 format(1x,'Number of Shannon channels : ',f16.2)
      write(166,965)evimax,sevi
  965 format(1x,'Evidence at maximum        : ',f16.2,
     -'  +- ',f10.2,' ')
      write(166,979)pval
  979 format(1x,'Probability of chi-square  : ',f16.7)
      if(pval.gt.0.003) then
        write(166,980)'Yes, probably ok'
      else
        if (chi2r.le.1) then
          write(166,980)'Probably overestimated'
        else
          write(166,980)'Probably underestimated'
        endif
      endif
c      write(166,980)answer 
  980 format(1x,'Are the exp errors correct : ',a)
      call system_clock(clock)
      cpu=(clock-clockold)*0.001
      write(166,955)cpu,cpumax
  955 format(1x,'Cpu time used: ',f8.1,
     -' seconds of maximum ',f8.1,' seconds')
      else
      write(166,*)'This logfile is calculated when p(r) is estimated'
      endif
      close(166)
    
c****************************************
c     End OUTPUT for download
c****************************************
      smear1=0.
      smear2=0.
      open(unit=137,file='deconout.d',status='unknown')
      write(137,*)dmax
      write(137,*)cav
      write(137,*)rtio
      write(137,*)background
      write(137,*)qmin
      write(137,*)qmax 
      write(137,*)smear1
      write(137,*)smear2
      write(137,4)answer5
      close(137)
     
      WRITE(6,900)BNAME
  900 FORMAT(1X,'p(r) in               :   ',A)
      if(etaest.ne.0)  WRITE(6,906)GsNAME
  906 FORMAT(1X,'p_struct(r) in        :   ',A)
      if(etaest.ne.0)  WRITE(6,904)gxNAME
  904 FORMAT(1X,'p_excl(r) in          :   ',A)
      WRITE(6,901)FNAME
  901 FORMAT(1X,'Fit of data in        :   ',A)
      write(6,902)'rescaled.d'
  902 format(1x,'Rescaled data in      :   ',A)
      WRITE(6,912)hxNAME
  912 FORMAT(1X,'Data used in          :   ',A)
      if(etaest.ne.0) WRITE(6,910)sNAME
  910 FORMAT(1X,'Structure factor in   :   ',A)
      WRITE(6,915)KNAME
  915 FORMAT(1X,'Rg, I(0) and Dmax in  :   ',A)
      WRITE(6,919)'parameters.d'
  919 FORMAT(1X,'Parameters in         :   ',A)
      write(6,*)
      write(6,947)cpu
  947 format(1x,'Cpu time used         :   ',f4.1,
     -' seconds')
      CLOSE(14)
      CLOSE(20) 
      CLOSE(36) 
      CLOSE(37) 

      x(1)=x(nn+1)
      call plot(xf(ptot),f,ptot,fm,mtot,nadd,aname,x,cav,alp,dot,m)
      open(unit=37,file='p.bat',status='unknown')
      write(37,9997)plname
 9997 format(1x,'wgnuplot ',a,' -')
      goto 9999

 9999 CLOSE(11)
      write(6,*)'end'
      stop
      END

c*************************************************************************
c     Powell algorithm from Numerical Recipes
c     finding the optimum set of parameters
c     via the evidence, which is calculated in the function FUNC
c*************************************************************************
      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
      INTEGER iter,n,np,NMAX,ITMAX
      REAL fret,ftol,p(np),xi(np,np),func
      EXTERNAL func
      PARAMETER (NMAX=20,ITMAX=200)
CU    USES func,linmin
      INTEGER i,ibig,j
      REAL del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
      common /data9/ niter
      common /data11/ itetot,itetotmax
call 1
      fret=func(p)
      do 11 j=1,n
        pt(j)=p(j)
11    continue
      iter=0
1     iter=iter+1

      niter=iter
      fp=fret
      ibig=0
      del=0.
      do 13 i=1,n
        do 12 j=1,n
          xit(j)=xi(j,i)
12      continue
        fptt=fret
c        write(6,*)'***',p(1),p(2),p(3),xit(1),xit(2),xit(3),n
        call linmin(p,xit,n,fret)
c        write(6,*)'**',p(1),p(2),p(3),xit(1),xit(2),xit(3),n
        if(abs(fptt-fret).gt.del)then
          del=abs(fptt-fret)
          ibig=i
        endif
13    continue
      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))then
c      write(6,*)'end here',fp,fret,ftol
      return
      endif
ccc      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
      do 14 j=1,n
        ptt(j)=2.*p(j)-pt(j)
        xit(j)=p(j)-pt(j)
        pt(j)=p(j)
14    continue
call 2
c      write(6,*)'fptt'
      fptt=func(ptt)
      if(fptt.ge.fp)goto 1
      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2

      if(t.ge.0.)goto 1
      write(6,*)'call linmin'
      call linmin(p,xit,n,fret)
      do 15 j=1,n
        xi(j,ibig)=xi(j,n)
        xi(j,n)=xit(j)
15    continue
      goto 1
99    END

C     Used for the Powell algorithm - also from Numerical Recipes
      SUBROUTINE linmin(p,xi,n,fret)
      INTEGER n,NMAX
      REAL fret,p(n),xi(n),TOL
c      PARAMETER (NMAX=50,TOL=5.e-1) ???????
      PARAMETER (NMAX=50,TOL=1.e-3)
CU    USES brent,f1dim,mnbrak
      INTEGER j,ncom
      REAL ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com/ pcom,xicom,ncom
      EXTERNAL f1dim
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
      ax=0.
      xx=1.
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=brent(ax,xx,bx,f1dim,TOL,xmin)
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
      return
      END

C     Used for the Powell algorithm - also from Numerical Recipes
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
      INTEGER ITMAX
      REAL brent,ax,bx,cx,tol,xmin,f,CGOLD,ZEPS
      EXTERNAL f
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0e-4)
      INTEGER iter
      REAL a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x)) 
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
        fu=f(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
ccc      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END

C     Used for the Powell algorithm - also from Numerical Recipes
      FUNCTION f1dim(x)
      INTEGER NMAX
      REAL f1dim,func,x
      PARAMETER (NMAX=50)
CU    USES func
      INTEGER j,ncom
      REAL pcom(NMAX),xicom(NMAX),xt(NMAX)
      COMMON /f1com/ pcom,xicom,ncom
c      write(6,*)'x = ',x,xicom(1),xicom(2),xicom(3)
      test=0
      do 10 j=1,ncom
      test=test+xicom(j)
 10   continue
      if(test.eq.0) xicom(1)=.5
      do 11 j=1,ncom
c        write(6,*)j,xt(j),pcom(j),x,xicom(j)
        xt(j)=pcom(j)+x*xicom(j)
11    continue
call 3
      f1dim=func(xt)
c      write(6,*)'call 3 ok',f1dim

      return
      END

C     Used for the Powell algorithm - also from Numerical Recipes
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      REAL ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
      EXTERNAL func
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
      common /data11/ itetot,itetotmax

call 4+5
c      write(6,*)'call 4+5'
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
call 6
c      write(6,*)'call 6'
      fc=func(cx)
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
call 7
          fu=func(u)
c      write(6,*)'call 7'
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
call 8
c      write(6,*)'call 8'
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0.)then
call 9
c      write(6,*)'call 9'
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
call 10
c      write(6,*)'call 10'
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
call 11
c      write(6,*)'call 11'
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
call 12
c      write(6,*)'call 12'
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END

c **********************************************************************
c    Function gammp og gammq - to calculate chi2 distribution
c    chi-square is a special case of the gamma functions
c    From numerical recipes for fortran 77, section 6.2
c **********************************************************************

      FUNCTION gammp(a,x)
      REAL a,gammp,x
c     USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)write(6,*)'bad arguments in gammp' 
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf 
      endif
      return 
      END

      FUNCTION gammq(a,x)
      REAL a,gammq,x
c     USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)write(6,*) 'bad arguments in gammq' 
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf 
      endif
      return 
      END

      SUBROUTINE gser(gamser,a,x,gln) 
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS 
      PARAMETER (ITMAX=100,EPS=3.e-7)
c     USES gammln
      INTEGER n
      REAL ap,del,sum,gammln 
      gln=gammln(a) 
      if(x.le.0.)then
        if(x.lt.0.)write(6,*)'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a 
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del 
        if(abs(del).lt.abs(sum)*EPS)goto 666
11    continue
      write(6,*)'a too large, ITMAX too small in gser'
 666  gamser=sum*exp(-x+a*log(x)-gln)
      return 
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
c     USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln 
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b 
        if(abs(d).lt.FPMIN)then
          d=FPMIN 
          c=b+an/c 
        endif
        if(abs(c).lt.FPMIN)then
          c=FPMIN  
          d=1./d
          del=d*c
          h=h*del
          if(abs(del-1.).lt.EPS)goto 667 
        endif
11    continue
      write(6,*)'a too large, ITMAX too small in gcf'
 667  gammcf=exp(-x+a*log(x)-gln)*h
      return 
      END

      FUNCTION gammln(xx) 
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     - 24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     --.5395239384953d-5,2.5066282746310005d0/ 
      x=xx
      y=x
      tmp=x+5.5d0 
      tmp=(x+0.5d0)*log(tmp)-tmp 
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y 
11    continue
      gammln=tmp+log(stp*ser/x) 
      return
      END

c*****************************************************************************
c     Function FUNC calculates the optimal solution (p_1, p_excl, p_struct)
c     for a given set of parameters (alpha, d, eta, ratio)
c     the function returns the value for the evidence for the optimal solution
c     (the evidence is used by the Powell algorithm to calculate the optimal set
c     of parameters)
c*****************************************************************************
      function func(xh)
     
      PARAMETER (NMAX=1000)
      PARAMETER (NDIST=1000)
      parameter (ndim=4)

      REAL A(NMAX,0:NMAX),B(0:NMAX,0:NMAX)
      real FT(NMAX,NMAX),SMEAR(NMAX,NMAX),xhmin(ndim)
      REAL X(NMAX),Y(NMAX),SD(NMAX),XF(0:NMAX),F(0:NMAX),u(nmax,nmax)
      REAL M(0:NMAX),YSUM(0:NMAX),FM(NMAX),w(nmax),xh(ndim),func
      real sigma(0:nmax),sigf(0:nmax),qmax,qmin,omega
      real ftot(0:nmax,ndist),prob(ndist)
      real ng,w1(nmax),w2(nmax)
      character*1 answer,answer2,answer5,answerm
      integer clock,clockold

      common /data1/ nextra,dmax,eta,diam,nfactor,nadd
      common /data2/ c,cnst,ntot2,rlogdet
      common /data3/ x,sd,y,xf,ntot,prob,ftot,nof,sigma
      common /data4/ mtot
      common /data5/ ratio,ratioold
      common /data6/ alphaest,dest,etaest
      common /data7/ f,sigf
      common /data8/ nfit,nerror,answer,answer2
      common /data9/ niter
      common /data10/ xprec,dotsp
      common /data11/ itetot,itetotmax
      common /data15/ answer5
      common /data25/ csmear,cexp
      common /data16/ xhmin,evidencesave,ngrid
      common /cputime/ clockold,cpumax,answerm
      common /test/ ndtest

      if(sigf(0).eq.0) jmin=1
      if(sigf(0).ne.0) jmin=0

      if(etaest.eq.0) nextra=0

      if(nfit.eq.4) then
      eta=xh(1)
      diam=xh(2)
      ratio=xh(3)
      alpha=exp(xh(4))
      endif

      if(nfit.eq.3) then
        if(etaest.lt.1.e-6) then
        diam=xh(1)
          if(xh(2).le.0.05) xh(2)=0.05
          if(xh(2).ge.20) xh(2)=20
        ratio=xh(2)
        alpha=exp(xh(3))
        eta=0.00000001
        endif
      
        if(etaest.ge.1.e-6) then
        eta=xh(1)
        diam=xh(2)
        alpha=exp(xh(3))
        ratio=ratio
        endif
      endif

      if(nfit.eq.2) then
        diam=xh(1)
        alpha=exp(xh(2))
        ratio=ratio
        eta=0.00000001
      endif

      if(diam.gt.1.e-20.and.diam.lt.1.e20) goto 822
      write(6,*)'Problem!! Diameter = ',diam,nfit,xh
      diam=-diam
      evidence=1.e9
      func=evidence
      prob(nof)=0
      goto 9999
  822 continue

      if(eta.le.0) eta=-eta
      if(eta.ge.0.5) eta=0.45+eta/1000
      if(ratio.le.0.05) ratio=0.05+ratio/1000
      if(ratio.ge.20) ratio=20-ratio/1000

      dfx=diam/ntot
      dmax=diam
      omega=0.5
      xprec=0.999

      maxit=1000+niter*200
      if(ntot.ge.60) maxit=500+niter*200
      if(ntot.ge.80) maxit=1000+niter*200
      if(ntot.ge.60) omega=0.5
      if(maxit.ge.8000) maxit=2000
      if(answerm.eq.'Q'.or.answerm.eq.'M') maxit=maxit*0.7
      if(answerm.eq.'Q'.or.answerm.eq.'M') xprec=0.99
      if(answerm.eq.'Q'.or.answerm.eq.'M') omega=0.5
c      if(maxit.ge.8000) maxit=3000

      CALL PRIOR(M,NTOT,y,XF,DFX,NMAX)

      CALL TRANS(A,FT,SMEAR,X,XF,NTOT,DFX,NMAX)
      m(0)=0
      f(0)=0

      DO 48 I=jmin,NTOT+nextra
      ysum(i)=0.
      DO 47 K=1,MTOT
      YSUM(I)=YSUM(I)+Y(K)*A(K,I)/SD(K)
   47 CONTINUE
   48 CONTINUE
 
      DO 14 I=jmin,NTOT+nextra
      DO 13 J=jmin,NTOT+nextra
        BIJ=0
          DO 12 K=1,MTOT
          BIJ=BIJ+A(K,I)*A(K,J)/SD(K)
   12     CONTINUE
        B(I,J)=BIJ
   13 CONTINUE
   14 CONTINUE

 1399 C1=0
      C2=0
      SUMM=0
      DO 1411 I=1,3 
      FM(I)=0
      DO 1400 J=1,NTOT+nextra
      FM(I)=FM(I)+A(I,J)*M(J)
 1400 CONTINUE
      C1=C1+FM(I)/SD(I)
      C2=C2+Y(I)/SD(I)
 1411 CONTINUE
      SUMM=C2/C1
      sumchec=0
      DO 35 j=1,NTOT+nextra
      M(j)=SUMM*M(j)
      F(j)=M(j)*1.001                      
      sumchec=sumchec+f(j)
   35 CONTINUE

      do 3811 j=ntot+1,ntot+nextra/2
      f(j)=1.001*M(j)
      sigma(j)=0.0001
 3811 continue
      do 384 j=ntot+nextra/2+1,ntot+nextra
      f(j)=m(j)
  384 continue
C*************************************************************************
C     Calculation of p(r) for given diameter and alpha
C*************************************************************************

      ITE=1
   50 ITE=ITE+1
ccccccccccccccccccccccccccccccccccccccccccccccccc
cc    0  < J < Ntot
ccccccccccccccccccccccccccccccccccccccccccccccccc
      DO 63 J=jmin,NTOT
      sigma(j)=abs(m(j)+1.e-10)*1
      if((answer5.eq.'K').or.(answer5.eq.'N').or.
     -answer5.eq.'B'.or.answer5.eq.'C') then
      sigma(j)=1
      goto 63
      endif
      if(answerm.eq.'M') goto 63
      if(m(j).le.0) m(j)=-m(j)+1.e-10
      if(f(j).le.0) f(j)=-f(j)+1.e-10
 6633 continue
   63 CONTINUE

      SUMM=0.                              
      SUMF=0.
      do 6631 j=1,ntot
      SUMM=SUMM+abs(M(J))
      SUMF=SUMF+abs(F(J))
 6631 continue
      SUMFACT=SUMF/SUMM
c*****************************
c     Maxent prior
c*****************************
      if(answerm.eq.'M'.or.answerm.eq.'Q') then
      do 6664 j=1,ntot
      M(J)=abs(SUMFACT*M(J))
c      sigma(j)=1.
c      sigma(j)=ntot/summ
      sigma(j)=abs(m(j))
 6664 continue
      goto 6465
      endif
c******************************
c     Smoothness prior
c******************************
      DO 64 J=2,NTOT-1
      M(J)=SUMFACT*M(J)
      m(j)=(f(j-1)+f(j+1))/2 
   64 CONTINUE
      m(1)=f(2)/2
      m(ntot)=m(ntot-1)/2
      if(answer5.eq.'C') m(1)=f(2)
c*******************************
c     smoothness end
c*******************************
 6465 m(0)=f(0)
      sigma(0)=10
cccccccccccccccccccccccccccccccccccccccccccccccccc
cc   Ntot1+1 < J < Ntot+nextra/2
cccccccccccccccccccccccccccccccccccccccccccccccccc
      SUMM=0.                              
      SUMF=0.
      DO 633 J=ntot+1,ntot+nextra/2
      if(f(j).le.0) SUMM=SUMM+M(J)
      if(f(j).le.0) SUMF=SUMF+F(J)
      SUMM=SUMM+M(J)
      SUMF=SUMF+F(J)
  633 CONTINUE
      SUMFACT=abs(SUMF/SUMM)
      sumav=summ/(nextra/2.)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc    ntot + ntextra/2+1  < J  < ntot + nextra
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      cut=.6
      DO 6435 J=ntot+nextra/2+1,ntot+nextra-1
      if(xf(j).le.cut*dmax) M(J)=0.0
      if(xf(j).gt.cut*dmax) M(J)=(f(j+1)+f(j-1))/2.
      sigma(j)=-sumav
      if(xf(j).le.0.5*cut*dmax) sigma(j)=-0.1*sumav
      if(nfit.eq.2) m(j)=0
      if(nfit.eq.2) sigma(j)=0.00001
 6435 continue
      m(ntot+nextra)=m(ntot+nextra-1)/2
      sigma(ntot+nextra)=sigma(ntot+nextra-1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc   END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   71 DO 100 I=jmin,NTOT+nextra
      FSUMI=0
      DO 75 J=jmin,NTOT+nextra
      FSUMI=FSUMI+B(I,J)*F(J)
   75 CONTINUE
      FSUMI=FSUMI-B(I,I)*F(I)  
        if(answerm.eq.'Q'.and.i.ge.1) then
        EE=2*(YSUM(I)-FSUMI)/ALPHA
        BB=M(I)*2*B(I,I)/ALPHA
        FX=ALPHA/(2*B(I,I))
        fx=fx*xx(bb,ee)
        else
      fx=(alpha*m(i)/sigma(i)+ysum(i)-fsumi)/(alpha/sigma(i)+b(i,i))
        endif
      F(I)=(1.-OMEGA)*F(I)+OMEGA*FX
      if(i.eq.0) f(0)=fx
  101 continue
  100 CONTINUE
c**************************************************************************
  110 S=0
      C=0
      GRADSI=0
      WGRADS=0
      GRADCI=0
      WGRADC=0
      DOTSP=0

      DO 150 I=jmin,NTOT+nextra
      SADD=-(f(i)-m(i))**2/sigma(i)
      if(answerm.eq.'Q'.and.i.ge.1) SADD=-F(I)*LOG(F(I)/M(I))+F(I)-M(I)
      S=S+SADD
  133 GRADSI=-2*(f(i)-m(i))/sigma(i)
      if(answerm.eq.'Q'.and.i.ge.1) GRADSI=-LOG(F(I)/M(I))
      WGRADS=WGRADS+GRADSI**2
  134 FSUMI=0
      DO 130 J=jmin,NTOT+nextra
      FSUMI=FSUMI+B(I,J)*F(J)
  130 CONTINUE
      GRADCI=2*(FSUMI-YSUM(I))
      WGRADC=WGRADC+GRADCI**2
      DOTSP=DOTSP+GRADSI*GRADCI
  150 CONTINUE

      DO 141 I=1,MTOT
      FM(I)=0
      DO 140 J=jmin,NTOT+nextra
      FM(I)=FM(I)+A(I,J)*F(J)
  140 CONTINUE
      CADD=(Y(I)-FM(I))**2/SD(I)
      C=C+CADD
  141 CONTINUE
      C=C/MTOT

      if((wgradc.ge.1.e30).or.(wgrads.ge.1.e30)) then
      evidence=1.e9
      evidence=sqrt(abs(wgradc))+sqrt(abs(wgrads))
      func=evidence
      prob(nof)=0
      do 799 j=1,ntot+nextra
 799  continue
      goto 9999
      endif

      if((wgradc.le.1.e-30).or.(wgrads.le.1.e-30)) then
      evidence=1.e9
      evidence=1./sqrt(abs(wgradc))+1./sqrt(abs(wgrads))
      func=evidence
      prob(nof)=0
      goto 9999
      endif

 565  WGRADS=SQRT(abs(WGRADS))
      WGRADC=SQRT(abs(WGRADC))
      egrad=wgrads*alpha/wgradc
      DOTSP=DOTSP/(WGRADS*WGRADC)

 1621 IF(ITE.EQ.MAXIT) GOTO 1799
      IF((ITE.ge.100).and.(dotsp.ge.xprec)) GOTO 1799
      goto 50
C************************************************************************
C     Estimate of p(r) written to file
C     and the probability for this particular solution is calculated
C************************************************************************
 1799 continue
      do 1801 i=1,ntot
      do 1797 j=1,ntot
      if(answer5.eq.'N') then
      u(i,j)=b(i,j)/alpha
      else if(answerm.eq.'M'.or.answerm.eq.'Q') then
      u(i,j)=sqrt(abs(m(i)*m(j)))*b(i,j)/alpha
      else
      u(i,j)=sqrt(abs(f(i)*f(j)))*b(i,j)/alpha
      endif
      if(abs(i-j).eq.0) u(i,j)=u(i,j)+1.
 1797 continue
 1801 continue 
      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W,smear)
 1813 rlogdet=0
      do 1802 i=1,ntot
      rlogdet=rlogdet+log(abs(w(i)))  
 1802 continue
c******************************************************
c START of calculation of Ng for volume fraction = 0
c******************************************************
      if(eta.gt.0.0001) goto 5000
c  numerator:
      do 4801 i=1,ntot
      do 4797 j=1,ntot
      if(answer5.eq.'N') then
      u(i,j)=b(i,j)
      else if(answerm.eq.'M'.or.answerm.eq.'Q') then
      u(i,j)=sqrt(abs(m(i)*m(j)))*b(i,j)
      else
      u(i,j)=sqrt(abs(f(i)*f(j)))*b(i,j)
      endif
 4797 continue
 4801 continue 

      call SVDCMP(u,ntot,Ntot,nmax,Nmax,W1,smear)

      ng=0
      do 4810 j=1,ntot
      ng=ng+w1(j)/(alpha*w(j))
 4810 continue
      ftot(nmax-12,nof)=ng
c********************************************************
c END of calculation of Ng
c********************************************************
 5000 continue 
      rnorm=0
c     saves 1/3 of the cpu-time:
      if(eta.le.0.00001) goto 2228
ccc-----------------------------------------------------
ccc   to ntot above   - extra below
ccc   new sub matrix of dimensions nextra/2 x nextra/2
ccc-----------------------------------------------------

      ntot2=nextra/2
      ntotx=ntot-ntot2

      do 1794 k1=ntot+ntot2+1,ntot+nextra
      do 1811 k2=ntot+ntot2+1,ntot+nextra
      i=k1-ntot-ntot2
      j=k2-ntot-ntot2
      u(i,j)=sqrt(abs(f(i)*f(j)))*b(i,j)/alpha
      if(abs(i-j).eq.0) u(i,j)=u(i,j)+1
 1811 continue
 1794 continue

      call SVDCMP(u,ntotx,Ntotx,nmax,Nmax,W,smear)

      do 1822 i=1,ntotx
      rlogdet=rlogdet+log(abs(w(i)))  
 1822 continue

 2228  evidence=0.5*rnorm-log(abs(dmax))+
     -(alpha*S-0.5*c*mtot)-0.5*rlogdet-log(abs(alpha))
c for lamellar symmetry smoothness needs to be weighted more.. ???
        if(answer5.eq.'C') evidence=evidence+alpha*S*20
        if(eta.gt.0.00001) then
        evidence=evidence-log(abs(eta))
        write(6,*)'eta > 0.00001'
        endif
      if(ratio.ge.1.) axial=ratio
      if(ratio.lt.1.) axial=1./ratio
      if(answer.eq.'E') evidence=evidence+log(axial)
c******************************************************************
c******************************************************************
      corm=30.
      if(evidence.le.0) then
      if((dotsp.lt.xprec).and.(nof.le.20)) evidence=evidence*corm
      if((dotsp.lt.xprec).and.(nof.gt.20)) evidence=evidence*corm
      endif
      if(evidence.gt.0) then
      if((dotsp.lt.xprec).and.(nof.le.20)) evidence=evidence/corm
      if((dotsp.lt.xprec).and.(nof.gt.20)) evidence=evidence/corm
      endif
c      if(dotsp.lt.xprec.and.answerm.eq.'M') evidence=-1.e10-ite*10000.
      evidence=-evidence
      func=evidence
  161 FORMAT(1X,I4,1x,f6.2,1X,f10.4,
     -1X,F10.4,1x,f7.3,2x,f7.4,1x,e8.1,1x,f10.1)
         
      WRITE(6,161)ite,dmax,log(abs(ALPHA))
     -,c,ratio,DOTSP,alpha*s,evidence
   
c  180 WRITE(6,161)ite,dmax,log(abs(ALPHA))
c     -,c,ratio,DOTSP,f(20),evidence

      itetot=itetot+ite
      call system_clock(clock)
      cpu=(clock-clockold)*0.001

      if((itetot.ge.itetotmax).or.(cpu.ge.cpumax)) then

      open(unit=88,file='plot1.pl',status='unknown')
	write(88,*)'set term gif'
	write(88,*)'set nokey'
	write(88,*)'set output "fig1.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88)

      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig2.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot3.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig3.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot4.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig4.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88)

      open(unit=88,file='plot5.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig5.gif"'
       write(88,*)
     -'set label "Too many iterations or max cpu time." at -6,0.55'
       write(88,*)
     -'set label "Try reducing the number of points." at -4.5,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
        close(88)

        close(88)
      write(6,*)'error4',itetot,cpu,' > ',itetotmax,cpumax

      STOP
      endif

      goto 8142

      DO 8141 I=1,MTOT-nadd
      FM(I)=0
      DO 8140 J=jmin,2*ntot
      fm(i)=FM(I)+A(I,J)*F(J)
 8140 CONTINUE
      fmi=0
      DO 8143 J=jmin,ntot
      FMI=FMI+A(I,J)*F(J)
 8143 CONTINUE
      write(213,*)x(i),fm(i)/fmi,0
 8141 CONTINUE

 8142 continue
      if(nof.ge.1) then

      do 997 j=jmin,ntot+nextra
      ftot(j,nof)=f(j)
  997 continue
      prob(nof)=-evidence
      ftot(nmax-1,nof)=dmax
      ftot(nmax-2,nof)=c
      ftot(nmax-3,nof)=log(alpha)
      ftot(nmax-4,nof)=dotsp     
      ftot(nmax-6,nof)=S
      ftot(nmax-7,nof)=evidence
      ftot(nmax-11,nof)=ratio

      nof=nof+1
      if(nof.gt.ndist) then
      write(6,*)'Error - too slow convergence => too many files'
      write(6,*)'Erasing first estimates  => max 1000 solutions'
      nof=1
      endif
 
      endif

      if(nof.eq.2) evimin=1.e6
      if(nof.eq.2) evimin2=1.e6

      if(nof.gt.2) then
      if(evidence.lt.evimin2) then
      evimin2=evidence
      xhmin(1)=xh(1)
      xhmin(2)=xh(2)
      xhmin(3)=xh(3)
      xhmin(4)=xh(4)
      evidencesave=evidence
      endif
      endif

      if(nfit.eq.4) goto 9998
c     if(nerror.ne.0) goto 9999
c     for test of alpha,dmax do not calculate ratio
      if(ngrid.eq.1) goto 9999
      if(dotsp.le.xprec) goto 9999

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if((nof.ge.5).and.(nof.le.250).and.(answer.eq.'M')) then
      if(evidence.ge.evimin) goto 9999
      evimin=evidence
      sum=0
      sum1=0
      sum3=0
      do 3139 j=2,ntot
      sum=sum+f(j)
      sum1=sum1+xf(j)*f(j)
      sum3=sum3+xf(j)**3*f(j)
 3139 continue
      t=(sum3/sum)/(sum1/sum)**3

c     prolate:
cccccccccccccxxxxxxxxxxx
      if(ratio.ge.1.00) then
      if(t.le.1.4003198) t=1.4003198
      if(t.le.1.617)
     -ax=1+1.318*(t-1.4003198)**0.47+1.633*(t-1.4005)
      if(t.ge.1.617)
     -ax=1+1.603*(t-1.4)**4.594+3.1702*(t-1.30244)
      ax=ax+0.0003
      endif
ccccccccccccccxxxxxxxxxxx
c     oblate
ccccccccccccccyyyyyyyyyy
      if(ratio.lt.1.00) then
      if(t.le.1.4003198) t=1.4003198
      ax=1-0.73649*(t-1.4)**0.190812-1.6129*(t-1.5017)
      endif
ccccccccccccccyyyyyyyyyy
      ratioold=ratio
      ratio=0.25*ax+0.75*ratio
      endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 9998 if(ratio.le.0.05) ratio = 0.05
      if(ratio.ge.20) ratio = 20

 9999 continue
c      write(6,*)'exit from func'
      return
      end

**************************************************************************
*     Definition of the prior
**************************************************************************
      SUBROUTINE PRIOR(M,NTOT,y,XF,DFX,NMAX)
      REAL M(0:NMAX),XF(0:NMAX),PPMIN,y(nmax)
      real pi
      integer clockold
      character*1 answerm
      common /data1/ nextra,dmax,eta,diam,nfactor,nadd
      common /data2/ c,cnst,ntot2,rlogdet
      common /data4/ mtot
      common /cputime/ clockold,cpumax,answerm
      common /data5/ ratio,ratioold
      
      PI=acos(-1.)
      diam=dmax

      DO 12 J=1,NTOT
      XF(J)=DFX*J
   12 CONTINUE
      xf(0)=0

c dmaxtotal = nfactor*dmax

      nfactor=3
c****************************************************
c     prior for p_1(r)
c****************************************************
      PPMIN=0.005
      DIAM=dfx*ntot

      CALL sphere(M,NTOT,DFX,Y(1),PPMIN,DIAM)
      do 301 j=0,ntot
      m(j)=m(j)*4*pi*xf(j)**2
 301  continue
      if(answerm.eq.'M'.or.answerm.eq.'Q')
     -call ellips(m,ntot,nmax,dmax,ratio,ppmin)
c****************************************************
c     prior for p_excl(r)
c****************************************************
 88   dfx=dmax/ntot
      DIAM=dfx*ntot
      xmin=0
      do 89 j=ntot+1,ntot+nextra/2      
      xf(j)=(j-ntot)*2*nfactor*dfx
  89  continue

c      ratio 1.0
c****************************************************
      if(ratio.eq.1000) then
      do 90 j=ntot+1,ntot+nextra/2      
c rescale of xf, to get max r = 4 (max radius = 1)
      r=2.*xf(j)/diam
        if(r.le.2) then
        c=1-(3./16)*r**3+(9./160)*r**4-r**6/2240
        goto 10
        endif
        if(r.gt.2) then
        c=8-(144/(35*r))-(18*r/5)+(5*r**3/16)-(9*r**4/160)+(r**6/2240)
        endif
        if(r.gt.4) then
        c=1.e-15
        endif
c from above gamma(0)=1 og gamma(2*diam)=0
  10  continue
      m(j)=-c*xf(j)**2
      if((m(j).lt.xmin).and.(j.ge.1)) xmin=m(j)
  90  continue
      endif
c-------------ratio 1.0 end----------------------------------

c      ratio general
c****************************************************

      call excl(xf,m,diam,ratio,ntot,nmax)

c****************************************************
cc  prior for p_struct(r)
c****************************************************
      xxmax=0
      fa = 5.37407
	fb = 0.67996
	fc = 0.0898459
	fd = 26.129
	fe = 41.9965
	ff = 56.6519
	fg = 9.79245
	fh = 6.92455
	fi  = 4.86883          

      do 100 j=ntot+nextra/2+1,ntot+nextra
      xf(j)=(j-ntot-nextra/2)*nfactor*2*dfx
      xf2=(20/diam)*xf(j)
      fx=fa*exp(-((xf2-fd)/fg)**2)-
     -fb*exp(-((xf2-fe)/fh)**2)+fc*exp(-((xf2-ff)/fi)**2) 
      m(j)=fx*(20/diam)**2
      if(m(j).gt.xxmax) xxmax=m(j)
 100  continue

c     For hard spheres: max(p_struct)/max(p_excl)=0.341
c     For eta = 0.1: max(p_1)=2.010*max(p_excl)

      fact=0.341*(-xmin/xxmax)

      do 200 j=ntot+nextra/2+1,ntot+nextra
      m(j)=fact*m(j)
 200  continue
c****************************************************
cc    The integral of p_excl during fitting.
c****************************************************
c     Integral(p_excl) is proportional to the volume of particle.
c     By multiplication below with ratio**n the area of p_excl is conserved
c     when changing ratio, which improves the convergence of the program
c     (otherwise changing the axial ratio will mainly effect I(0) etc.)

c     Regarding the interpretation of eta after the renormalization:
c     From the calculated p_excl it appears that p_excl/(8*p_1) = 1
c     for a large interval around the axial ratio 1.0
c     Consequently the calculated p_excl is simply divided by 8*p_1
c     to give an estimate of the volume fraction.
c
c     NB this should be corrected for the increased excluded volume
c     for deviation for spherical symmetry
c   
      do 300 j=ntot+1,ntot+nextra
      if(ratio.ge.1) m(j)=eta*m(j)*RATIO**2
      if(ratio.lt.1) m(j)=eta*m(j)/RATIO**1
 300  continue

 999  RETURN
      END
c*******************************************************************



c********************************************************************
c  Ellipsoidal prior start
c********************************************************************
      subroutine ellips(m,ntot,nmax,dmax,v,ppmin)
      REAL M(0:nmax),pi4,ppmin

      pi4=4*3.14159
      dr=dmax/ntot

      if(ax.lt.0.05) ax=0.05
      if(ax.gt.20) ax=20.

ccccccccccccccccccccccccccc
c     Prolate
ccccccccccccccccccccccccccc
      if(v.gt.1.001) then

      A=dmax/v

      VM1=SQRT(V**2-1)
      C11=1./(4*V**3)+3./(8*V)+3.*V/(8*VM1)*atan(VM1)
      C12=1./(2*V)+V/(2*VM1)*aTAN(VM1)

      DO 20 J=0,NTOT
      R=DR*J
        IF(R.GT.Dmax) THEN
        M(J)=0
        GOTO 20
        ENDIF
      c1=0
      c2=0
      C1=0.5*C11*(R/A)**3-3./2*C12*(R/A)+1
        IF(R.GT.A) THEN
        SQ=SQRT((R/A)**2-1)
        C2=3*V/(16*VM1)*(SQ*(2*A/R+R/A)+((R/A)**3-4*R/A)*aTAN(SQ))
        ENDIF
      M(J)=abs(pi4*R**2*(c1-c2))
   20 CONTINUE
      endif
ccccccccccccccccccccccccccc
c     Oblate
ccccccccccccccccccccccccccc
      if(v.lt.0.999) then
      A=dmax
      VM1=SQRT(-(V**2-1))
      C11=1./(4*V**3)+3./(8*V)+3.*V/(8*VM1)*atanh(VM1)
      C12=1./(2*V)+V/(2*VM1)*aTANh(VM1)
      DO 10 J=0,NTOT
      R=DR*J
        IF(R.GT.Dmax) THEN
        M(J)=0
        GOTO 10
        endif
      C=0.5*C11*(R/A)**3-3./2*C12*(R/A)+1
        IF(R.GT.v*A) THEN
        SQ=SQRT(-((R/A)**2-1))
        C=3*V/(16*VM1)*(SQ*(2*A/R+R/A)+((R/A)**3-4*R/A)*aTANh(SQ))
        ENDIF
      M(J)=abs(pi4*R**2*C)
   10 CONTINUE
      endif
ccccccccccccccccccccccccccc
c     Sphere
ccccccccccccccccccccccccccc
      if(v.ge.0.999.and.v.le.1.001) then
      do 30 j=1,ntot
      r=dr*j
        if(r.gt.dmax) then
        m(j)=0
        goto 30
        endif
      m(j)=pi4*R**2*(1-1.5*(R/Dmax)+.5*(R/Dmax)**3)
  30  continue
      endif
    
      sum=0.
      do 50 j=1,ntot
      sum=sum+m(j)
   50 continue

c      xmin=0.006*(sum/ntot)
      xmin=0.002*(sum/ntot)
      do 51 j=1,ntot
      if(m(j).lt.xmin) m(j)=xmin
   51 continue
      

   99 return 
      END

c********************************************
      function atanh(x)
      atanh=0.5*log((1+x)/(1-x))
      return
      end
c********************************************************************
c  Ellipsoidal prior end
c********************************************************************

********************************************************************* 
*     Calculation of excluded volume p_excl(r)
*********************************************************************
      subroutine excl(xf,m,diam,ax,ntot,nmax)
      real xf(0:nmax),m(0:nmax)
      real pi
      PI=acos(-1.)
c 
c     p_excl(r) for prolate and oblate ellipsoids of axial ratios between 0.1 and 10.
c     Maximum diameter of scatterer = 200.
c     Maximum dimension for p_excl(r) is 400. 
c
      x=ax
      if(x.le.3.33) then
	a               = -0.00640898
	b               = 0.00729263
	c               = -0.00256838
	d               = -0.000589684
	e               = 0.00204167
	f               = -0.000217912
      endif
      if((x.le.1.115).and.(x.ge.0.8)) then
	a               = 0.000472864
	b               = -0.00407438
	c               = 0.00645582
	d               = -0.00328017       
      e=0
      f=0
      endif
      atot=f*x**5+e*x**4+a*x**3+b*x**2+c*x+d
      if(x.gt.3.33) then
	a               = -0.000248348
	b               = 0.00450629
	c               = -0.0153692
      atot=a*x**2+b*x+c
      endif
  
 
      if(x.le.2.5) btot=1.98567
      if(x.gt.2.5) then
	a = -0.00122915
	b = 0.00507403
	c = 0.00197944
	d = 1.96828
      btot=a*x**3+b*x**2+c*x+d
      endif


      if(x.le.3.33) then
	a               = 0.00603866
	b               = -0.00682305
	c               = 0.00234489
	d               = 0.000597595
	e               = -0.00192782
	f               = 0.000206001
      endif
      if((x.le.1.115).and.(x.ge.0.8)) then
	a               = 1.72976e-005
	b               = 0.00258405
	c               = -0.00496145
	d               = 0.00277141
      e=0
      f=0
      endif
      ctot=f*x**5+e*x**4+a*x**3+b*x**2+c*x+d
      if(x.gt.3.33) then
	a               = 0.000346191
	b               = -0.00506421
	c               = 0.0160596
      ctot=a*x**2+b*x+c
      endif
  
      if(x.le.1.0) then
	a               = 4.77001e-012
	b               = 4.04391e-010
	c               = -4.85954e-010
	d               = -1.29272e-010
      endif
      if(x.gt.1.0) then
	a               = -3.1314e-010
	b               = 2.29874e-009
	c               = -4.71816e-009
	d               = 2.55382e-009
      endif
      dtot=a*x**3+b*x**2+c*x+d
      if(x.ge.3.4) dtot=0
      
      a=atot
      b=btot
      c=ctot
      d=dtot
      ax=x

c if max diameter of scatterer = 200 no change in x

      do 100 j=ntot+1,ntot+ntot/2
      x=xf(j)*(200/diam)
      c0=exp(-a*x**b-c*x**2+d*x**4) 
      if((x.gt.300).and.(ax.ge.2.5)) c0=0      
      m(j)=-c0*4*pi*xf(j)**2
  100 continue
  999 return
      end
********************************************************************* 
*     Correlation function for Sphere
*********************************************************************
      SUBROUTINE sphere(M,NTOT,DR,RX,PMIN,D)
      REAL M(0:NTOT)
      common /data1/ nextra,dmax,eta,diam,nfactor,nadd
      common /data2/ c,cnst,ntot2,rlogdet
c      common /data4/ mtot
      DO 20 J=0,NTOT
      R=DR*j
      IF(R.GT.D) THEN
      M(J)=0.
      GOTO 20
      ENDIF
      M(J)=(1-1.5*(R/D)+.5*(R/D)**3)
   20 CONTINUE
      RETURN
      END

C*********************************************************************
C     Construction of transformation matrix A   
C     If instrument resolution is to be included it should be
C     defined in the matrix SMEAR and A ->  A*SMEAR
C*********************************************************************
      SUBROUTINE TRANS(A,FT,SMEAR,x,xf,NTOT,DFX,NMAX)
      REAL A(NMAX,0:NMAX),SMEAR(NMAX,NMAX),FT(NMAX,NMAX)
      REAL x(nmax),xf(0:nmax),m(0:nmax)
      character*1 answer5
      integer clockold
      common /data1/ nextra,dmax,eta,diam,nfactor,nadd
      common /data2/ c,cnst,ntot2,rlogdet
      common /data4/ mtot
      common /data15/ answer5
      common /data25/ csmear,cexp
 
      dr=xf(2)-xf(1)
      DO 4 I=1,mtot-nadd
      DO 3 J=0,NTOT
      
      if((x(i).eq.0).or.(xf(j).eq.0)) then
      a(i,j)=dr*0.1
      goto 3
      endif
      
      A(I,J)=dr*SIN(X(I)*XF(J))/(X(I)*XF(J))

      if(answer5.eq.'C') A(I,J)=dr*cos(X(I)*XF(J))/X(I)**2
      if(answer5.eq.'B')  A(I,J)=dr*bessj0(X(I)*XF(J))/X(I)
      if(answer5.eq.'S')  then
      rq=x(i)*xf(j)
      a(i,j)=dr*((3*(sin(rq)-rq*cos(rq)))/rq**3)**2
c     a(i,j)=a(i,j)*(4*3.14159*xf(j)**3/3)**2
      endif
      if(answer5.eq.'K') then
      a(i,j)=0
      drk=dmax/400
      nrmax=xf(j)/drk
      do 8 k=1,nrmax
      rk=k*drk
      a(i,j)=a(i,j)+
     -drk*dr*(xf(j)-rk)*rk*sin(x(i)*rk)/(x(i)*xf(j)**3)
   8  continue
      endif
   3  CONTINUE   
   4  CONTINUE

      aold=a(2,2)

c smearing************
      if(csmear.ne.0) then

      nt=100
      tmax=4/csmear
      if(cexp.eq.6) tmax=tmax/15.
      dt=tmax/nt
      renorm=2*csmear/sqrt(3.14159)*dt
  778 c2=csmear*csmear

      testx=0
      do 600 k=0,nt
      t=k*dt
      t2=t*t
      if(cexp.eq.6) t2=t2**3
      testx=testx+exp(-c2*t2)
  600 continue
      testx=testx-0.5
c      write(6,*)renorm, 1./testx
      if(cexp.eq.6) renorm=1./testx   

      DO 803 I=1,MTOT
      DO 802 J=1,NTOT
      aij=0
      DO 801 k=0,nt
      t=k*dt
      t2=t*t
      rq=sqrt(x(i)*x(i)+t2)*xf(j)
      if(rq.le.1e-15) rq=1.e-15
      if(cexp.eq.6.) t2=t2**3
      aij=aij+exp(-c2*t2)*sin(rq)/rq
  801 CONTINUE
      rq=x(i)*xf(j)
      if(rq.le.1e-15) rq=1.e-15
      aij=aij-0.5*sin(rq)/rq
      a(i,j)=renorm*aij*dr
  802 CONTINUE
  803 CONTINUE
      
      endif

c smearing************

      DO 14 I=mtot-nadd+1,mtot
      DO 13 J=0,NTOT
      if((x(i).eq.0).or.(xf(j).eq.0)) then
      a(i,j)=0
      goto 13
      endif
      A(I,J)=0
  13  CONTINUE
  14  CONTINUE

      dr=xf(ntot+2)-xf(ntot+1)

      DO 44 I=1,mtot
      DO 33 J=ntot+1,ntot+nextra
 
      A(I,J)=dr*SIN(X(I)*xf(j))/(X(I)*xf(j))

      if(answer5.eq.'C') A(I,J)=dr*cos(X(I)*XF(J))/X(I)**2
      if(answer5.eq.'B') A(I,J)=dr*bessj0(X(I)*XF(J))/X(I)
      if(answer5.eq.'S')  then
      rq=x(i)*xf(j)
      a(i,j)=dr*((3*(sin(rq)-rq*cos(rq)))/rq**3)**2
c      a(i,j)=a(i,j)*(4*3.14159*xf(j)**3/3)**2
      endif

  33  CONTINUE
  44  CONTINUE
      
  999 return
      end

c*******************************************************************
c     singular value decompostion from Numerical Recipes
c******************************************************************* 
      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=1000)
      DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
      G=0.0
      SCALE=0.0
      ANORM=0.0
      DO 25 I=1,N
        L=I+1
        RV1(I)=SCALE*G
        G=0.0
        S=0.0
        SCALE=0.0
        IF (I.LE.M) THEN
          DO 11 K=I,M
            SCALE=SCALE+ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,M
              A(K,I)=A(K,I)/SCALE
              S=S+A(K,I)*A(K,I)
12          CONTINUE
            F=A(I,I)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,I)=F-G
            IF (I.NE.N) THEN
              DO 15 J=L,N
                S=0.0
                DO 13 K=I,M
                  S=S+A(K,I)*A(K,J)
13              CONTINUE
                F=S/H
                DO 14 K=I,M
                  A(K,J)=A(K,J)+F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K= I,M
              A(K,I)=SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I)=SCALE *G
        G=0.0
        S=0.0
        SCALE=0.0
        IF ((I.LE.M).AND.(I.NE.N)) THEN
          DO 17 K=L,N
            SCALE=SCALE+ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,N
              A(I,K)=A(I,K)/SCALE
              S=S+A(I,K)*A(I,K)
18          CONTINUE
            F=A(I,L)
            G=-SIGN(SQRT(S),F)
            H=F*G-S
            A(I,L)=F-G
            DO 19 K=L,N
              RV1(K)=A(I,K)/H
19          CONTINUE
            IF (I.NE.M) THEN
              DO 23 J=L,M
                S=0.0
                DO 21 K=L,N
                  S=S+A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,N
                  A(J,K)=A(J,K)+S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,N
              A(I,K)=SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE
      DO 32 I=N,1,-1
        IF (I.LT.N) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,N
              V(J,I)=(A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,N
              S=0.0
              DO 27 K=L,N
                S=S+A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,N
                V(K,J)=V(K,J)+S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,N
            V(I,J)=0.0
            V(J,I)=0.0
31        CONTINUE
        ENDIF
        V(I,I)=1.0
        G=RV1(I)
        L=I
32    CONTINUE
      DO 39 I=N,1,-1
        L=I+1
        G=W(I)
        IF (I.LT.N) THEN
          DO 33 J=L,N
            A(I,J)=0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G=1.0/G
          IF (I.NE.N) THEN
            DO 36 J=L,N
              S=0.0
              DO 34 K=L,M
                S=S+A(K,I)*A(K,J)
34            CONTINUE
              F=(S/A(I,I))*G
              DO 35 K=I,M
                A(K,J)=A(K,J)+F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,M
            A(J,I)=A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J= I,M
            A(J,I)=0.0
38        CONTINUE
        ENDIF
        A(I,I)=A(I,I)+1.0
39    CONTINUE
      DO 49 K=N,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM=L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z=W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
ccc          IF (ITS.EQ.300) PAUSE 'No convergence in 300 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
            Z=SQRT(F*F+H*H)
            W(J)=Z
            IF (Z.NE.0.0) THEN
              Z=1.0/Z
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3       CONTINUE
49    CONTINUE
      RETURN
      END
      
***********************************************************************
*     Random number generator from "Numerical Recipes"
***********************************************************************
      FUNCTION ran1(idum)
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
C     REAL MBIG,MSEED,MZ
      REAL ran1,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
C     PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
      INTEGER i,iff,ii,inext,inextp,k
      INTEGER mj,mk,ma(55)
C     REAL mj,mk,ma(55)
      SAVE iff,inext,inextp,ma
      DATA iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(i)=ma(i)-ma(1+mod(i+30,55))
            if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
        inext=0
        inextp=31
        idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
      inextp=inextp+1
      if(inextp.eq.56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.MZ)mj=mj+MBIG
      ma(inext)=mj
      ran1=mj*FAC
      return
      END
c**********************************************************************
c     This makes a plot file to be used by wgnuplot
c**********************************************************************
      subroutine plot(dmax,f,ptot,fm,mtot,nadd,oname,x,chi,alpha,dot,m)
      parameter (nmax=1000)
      CHARACTER*18 oname
      character*10 aname
      chARACTER*20 plname,comment
c      CHARACTER*23 hxname,gsname,gxname,sname
      CHARACTER*26 HNAME
      character*40 screen
      character*10 bname
      character*5 fname,gname
      character*4 gsname,gxname,sname
      character*6 hxname
      character*1 answer8,answer5,answer7,answerm
      integer ptot,clockold

      common /data11/ itetot,itetotmax
      common /data15/ answer5
      common /data25/ csmear,cexp
      common /back/ background
      common /screen/ answer8,answer7
      common /cputime/ clockold,cpumax,answerm

      real xf(0:nmax),f(0:nmax),fm(nmax),x(nmax),m(0:nmax)
      plname='pl'//aname
      aname=oname

c      bname='me'//ANAME
c
c      FNAME='trans_'//BNAME        
c      GNAME='decon_'//BNAME        
c      HNAME='PRIOR_'//BNAME   
c      gsNAME='gs_'//BNAME   
c      gxNAME='gx_'//BNAME   
c      hxname='in_'//bname
c      sname='st_'//bname

      bname='estimate.d'
      FNAME='fit.d'        
      GNAME='dec.d'        
      HNAME='PRIOR_'//BNAME   
      gsNAME='gs.d'   
      gxNAME='gx.d'   
      hxname='data.d'
      sname='st.d'
      if(answer8.eq.'S') screen='set term gif medium size 500,360'
      if(answer8.eq.'L') screen='set term gif large size 1000,720'

      xmax=0
      do 9 i=1,mtot
      if(x(i).gt.xmax) xmax=x(i)
   9  continue
      fmax=0
      fmin=1.e20
      do 10 j=1,ptot
      if(f(j).gt.fmax) fmax=f(j)
      if(f(j).le.fmin) fmin=f(j)
  10  continue
      exclmin=0
      do 11 j=ptot+1,2*ptot
      if(f(j).lt.exclmin) exclmin=f(j)
  11  continue
      fmmin=1.e6
      fmmax=0
      do 12 i=1,mtot-nadd
      if(fm(i).lt.fmmin) fmmin=fm(i)
      if(fm(i).gt.fmmax) fmmax=fm(i)
  12  continue


      if(dot.ge.0.9) comment='Conv. OK'

      if(dot.le.0.9) comment='NOT OK !'

      mt=mtot-nadd
      fratiox=(fm(1)+fm(2)+fm(3))/(fm(mt-2)+fm(mt-1)+fm(mt))
      fratiox=abs(fratiox)
      write(6,*)
 
c**********************'' plot 1 ******************************************

        open(unit=88,file='plot1.pl',status='unknown')

	write(88,*)'set noclip points'
	write(88,*)'set clip one'
	write(88,*)'set noclip two'
	write(88,*)'set bar 1.000000'
	write(88,*)'set xdata'
	write(88,*)'set ydata'
	write(88,*)'set zdata'
	write(88,*)'set x2data'
	write(88,*)'set y2data'
	write(88,*)'set boxwidth'
	write(88,*)'set dummy x,y'
	write(88,*)'set format x "%g"'
	write(88,*)'set format y "%g"'
	write(88,*)'set format x2 "%g"'
	write(88,*)'set format y2 "%g"'
	write(88,*)'set format z "%g"'
	write(88,*)'set angles radians'
	write(88,*)'set nogrid'
	write(88,*)'set key title ""'
	write(88,*)'set nolabel'
	write(88,*)'set noarrow'
	write(88,*)'set nolinestyle'
	write(88,*)'set nologscale'
	write(88,*)'set offsets 0, 0, 0, 0'
	write(88,*)'set pointsize 1'
	write(88,*)'set encoding default'
	write(88,*)'set nopolar'
	write(88,*)'set noparametric'
	write(88,*)'set view 60, 30, 1, 1'
	write(88,*)'set samples 100, 100'
	write(88,*)'set isosamples 10, 10'
	write(88,*)'set surface'
	write(88,*)'set nocontour'
	write(88,*)'set clabel ','''%8.3g'''
	write(88,*)'set mapping cartesian'
	write(88,*)'set nohidden3d'
	write(88,*)'set cntrparam order 4'
	write(88,*)'set cntrparam linear'
	write(88,*)'set cntrparam levels auto 5'
	write(88,*)'set cntrparam points 5'
	write(88,*)'set data style points'
	write(88,*)'set function style lines'
	write(88,*)'set xzeroaxis lt -2 lw 1.000'
	write(88,*)'set x2zeroaxis lt -2 lw 1.000'
	write(88,*)'set yzeroaxis lt -2 lw 1.000'
	write(88,*)'set y2zeroaxis lt -2 lw 1.000'
	write(88,*)'set tics in'
	write(88,*)'set ticslevel 0.5'
	write(88,*)'set ticscale 1 0.5'
	write(88,*)'set mxtics default'
	write(88,*)'set mytics default'
	write(88,*)'set mx2tics default'
	write(88,*)'set my2tics default'
	write(88,*)'set xtics border mirror norotate autofreq' 
	write(88,*)'set ytics border mirror norotate autofreq' 
	write(88,*)'set ztics border nomirror norotate autofreq' 
	write(88,*)'set nox2tics'
	write(88,*)'set noy2tics'
	write(88,*)'set rrange [ * : * ] noreverse nowriteback'  
	write(88,*)'set trange [ * : * ] noreverse nowriteback'  
	write(88,*)'set urange [ * : * ] noreverse nowriteback'  
	write(88,*)'set vrange [ * : * ] noreverse nowriteback'
	write(88,*)'set xrange [ * : * ] noreverse nowriteback'  
	write(88,*)'set x2range [ * : * ] noreverse nowriteback'  
	write(88,*)'set yrange [ * : * ] noreverse nowriteback' 
	write(88,*)'set y2range [ * : * ] noreverse nowriteback'
	write(88,*)'set zrange [ * : * ] noreverse nowriteback' 
	write(88,*)'set zero 1e-008'
	write(88,*)'set lmargin 0'
	write(88,*)'set bmargin 0'
	write(88,*)'set rmargin 0'
	write(88,*)'set tmargin 0'
	write(88,*)'set nolog'
	write(88,*)'set nolabel'
	write(88,*)'set nokey'
	write(88,*)'reset'

	write(88,*)'set nokey'
	write(88,*)'set nolabel'

	write(88,*)'set yrange [0:*]'
	write(88,*)'set xrange [*:',xmax,']'


       xmin=x(1)
       xmax=x(mtot)
       if(x(1).ge.x(mtot))then
       xmax=x(1)
       xmin=x(mtot)
       endif
        write(88,744)aname,xmin,xmax
 744  format(1x,'set title "S(q) for ',a,' from ',f8.4,' to ',f8.4,'"')

	write(88,*)'set xlabel "q"  '
	write(88,*)'set ylabel "S(q)" '

	write(88,*)'set xtics 1'
      if(xmax.le.4) write(88,*)'set xtics 0.5'
	if(xmax.le.1.6) write(88,*)'set xtics 0.2'
	if(xmax.le.0.5) write(88,*)'set xtics 0.1'
	if(xmax.le.0.26) write(88,*)'set xtics 0.05'
	if(xmax.le.0.16) write(88,*)'set xtics 0.02'
	if(xmax.le.0.08) write(88,*)'set xtics 0.01'
	if(xmax.le.0.04) write(88,*)'set xtics 0.005'
	if(xmax.le.0.016) write(88,*)'set xtics 0.002'
	if(xmax.le.0.005) write(88,*)'set xtics 0.001'
	
  818 format(a)
      write(88,818)screen
	write(88,*)'set output "fig1.gif"'
	write(88,881)sname
  881 format(1x,'plot 1 w l 3,''',a,''' w l 1') 
        write(88,*)'exit gnuplot'
        close(88)

c**********************'' plot 2 ******************************************
       xl=80
       if(fratiox.ge.xl) open(unit=88,file='plot2.pl',status='unknown')
       if(fratiox.lt.xl) open(unit=88,file='plot3.pl',status='unknown')

	write(88,*)'set yrange [*:*]'
	write(88,*)'set xrange [*:*]'
	write(88,*)'set logsc y'
	write(88,*)'set xlabel "q"  0.0,0.0'
	write(88,*)'set ylabel ','''I(q)  [a.u.]''',' 0.0'

      if(itetot.ge.itetotmax) then
      write(88,*)'set label "No convergence" at 0.5,0.5'
      endif

	write(88,*)'set nokey'
      if(fmmin.le.0) fmmin=1.e-5
	if(xmax.le.1) write(88,*)'set xtics 0.1'
	if((xmax.ge.1).and.(xmax.le.4)) write(88,*)'set xtics 0.5'
	if(xmax.gt.4) write(88,*)'set xtics 1'
	write(88,*)'set  yrange [',0.8*fmmin,':',2*fmmax,']'
	write(88,*)'set  xrange [0:',xmax,']'

	write(88,*)'set xtics 1'
      if(xmax.le.4) write(88,*)'set xtics 0.5'
	if(xmax.le.1.6) write(88,*)'set xtics 0.2'
	if(xmax.le.0.5) write(88,*)'set xtics 0.1'
	if(xmax.le.0.3) write(88,*)'set xtics 0.05'
	if(xmax.le.0.16) write(88,*)'set xtics 0.02'
	if(xmax.le.0.08) write(88,*)'set xtics 0.01'
	if(xmax.le.0.04) write(88,*)'set xtics 0.005'
	if(xmax.le.0.016) write(88,*)'set xtics 0.002'
	if(xmax.le.0.005) write(88,*)'set xtics 0.001'

       xmin=x(1)
       xmax=x(mtot)
       if(x(1).ge.x(mtot))then
       xmax=x(1)
       xmin=x(mtot)
       endif
        write(88,444)aname,xmin,xmax
 444  format(1x,'set title "data + fit ',a,' from',f8.4,' to',f8.4,'"')

      if((ptot.lt.90).and.(fratiox.ge.xl*5.)) then
      write(88,476)"Try to increase the number",xmax*0.4,1.05*fmmax*0.8
      write(88,476)"of points used for p(r) or",xmax*0.4,1.05*fmmax*0.4
      write(88,476)"to reduce the used q-range",xmax*0.4,1.05*fmmax*0.205
 476  format(1x'set label "',a,'" at ',e8.2,',',e8.2)
      endif

      write(88,818)screen

       if(fratiox.ge.xl) write(88,*)'set output "fig2.gif"'
       if(fratiox.lt.xl) write(88,*)'set output "fig3.gif"'
	write(88,882)hxname,fname,fname
  882 format(1x,'plot ''',a,''' w e 1,''',a,''' w l 2,''',a,''' 
     -using ','''%lf %*lf %lf''',' w l 3')
        write(88,*)'exit gnuplot'
        close(88)

c**********************'' plot 3 ******************************************

       if(fratiox.ge.xl) open(unit=88,file='plot3.pl',status='unknown')
       if(fratiox.lt.xl) open(unit=88,file='plot2.pl',status='unknown')

	write(88,*)'set yrange [*:*]'
	write(88,*)'set xrange [*:*]'
	write(88,*)'set nolog'
	write(88,*)'set nolabel'
	write(88,*)'set xlabel "q"  0.0,0.0'
	write(88,*)'set ylabel ','''I(q)  [a.u.]''',' 0.0'
	write(88,*)'set nokey'
	write(88,*)'set  yrange [1.e-8:',fmmax*1.05,']'
	write(88,*)'set  xrange [0:',xmax,']'
      write(88,818)screen

	write(88,*)'set xtics 1'
      if(xmax.le.4) write(88,*)'set xtics 0.5'
	if(xmax.le.1.6) write(88,*)'set xtics 0.2'
	if(xmax.le.0.5) write(88,*)'set xtics 0.1'
	if(xmax.le.0.3) write(88,*)'set xtics 0.05'
	if(xmax.le.0.16) write(88,*)'set xtics 0.02'
	if(xmax.le.0.08) write(88,*)'set xtics 0.01'
	if(xmax.le.0.04) write(88,*)'set xtics 0.005'
	if(xmax.le.0.016) write(88,*)'set xtics 0.002'
	if(xmax.le.0.005) write(88,*)'set xtics 0.001'


       xmin=x(1)
       xmax=x(mtot)
       if(x(1).ge.x(mtot))then
       xmax=x(1)
       xmin=x(mtot)
       endif
       write(88,444)aname,xmin,xmax

       if(fratiox.ge.xl) write(88,*)'set output "fig3.gif"'
       if(fratiox.lt.xl) write(88,*)'set output "fig2.gif"'
      	write(88,882)hxname,fname,fname
        write(88,*)'exit gnuplot'
        close(88)

c**********************'' plot 4 ******************************************

        open(unit=88,file='plot4.pl',status='unknown')

	write(88,*)'set yrange [*:*]'
	write(88,*)'set xrange [*:*]'
	write(88,*)'set nokey'
       dmax=3*dmax
	if(dmax.le.10) write(88,*)'set xtics 4'
	if((dmax.gt.10).and.(dmax.le.20)) write(88,*)'set xtics 5'
	if((dmax.gt.20).and.(dmax.le.50)) write(88,*)'set xtics 10'
	if((dmax.gt.50).and.(dmax.le.100)) write(88,*)'set xtics 20'
	if((dmax.gt.100).and.(dmax.le.500)) write(88,*)'set xtics 100'
	if(dmax.gt.500) write(88,*)'set xtics 100'
	if(dmax.gt.1000) write(88,*)'set xtics 200'
	if(dmax.gt.2000) write(88,*)'set xtics 500'
	if(dmax.gt.5000) write(88,*)'set xtics 1000'
       dmax=0.3333*dmax

	write(88,*)'set nolabel'
	write(88,*)'set nologsc'
	write(88,*)'set ylabel "p(r)  [a.u.]" 0.0,0.0'
	write(88,*)'set xlabel ','''r''',' 0.0,0.0'
	write(88,*)'set xrange [0:',2.5*dmax,']'
	write(88,*)'set yrange [',1.05*exclmin,':',1.05*fmax,']'
	write(88,*)'set nokey'

c*******************
      write(88,574)aname,chi,alpha

      write(88,575)comment,dmax*1.4,1.05*fmax*0.9
c*****************

      write(88,818)screen
	write(88,*)'set output "fig4.gif"'
c	write(88,885)'pr.d',gsname,gxname,'pr.d',gsname,gxname
c 885  format(1x,'plot 0 w l 3,''',a,''' w e 1,''',a,''' w e 2 ,''',a,''' 
c     -  w e 4,
c     -  ''',a,''' w l 1,''',a,''' w l 2 ,''',a,''' 
c     -  w l 4')
        write(88,*)'plot "pr.d" with errorbars'
        write(88,*)'set yrange [*:*]'
	write(88,*)'set xrange [*:*]'
	write(88,*)'set nokey'

        write(88,*)'exit gnuplot'

        close(88)

c**********************'' plot 5 ******************************************

        open(unit=88,file='plot5.pl',status='unknown')
	write(88,*)'set nolabel'
	if(dmax.le.10) write(88,*)'set xtics 1'
	if((dmax.gt.10).and.(dmax.le.20)) write(88,*)'set xtics 2'
	if((dmax.gt.20).and.(dmax.le.50)) write(88,*)'set xtics 5'
	if((dmax.gt.50).and.(dmax.le.100)) write(88,*)'set xtics 10'
	if((dmax.gt.100).and.(dmax.le.500)) write(88,*)'set xtics 50'
	if(dmax.gt.500) write(88,*)'set xtics 100'
	if(dmax.gt.1000) write(88,*)'set xtics 200'
	if(dmax.gt.2000) write(88,*)'set xtics 1000'
	if(dmax.gt.5000) write(88,*)'set xtics 2000'
	write(88,*)'set nologsc'
      
      if(answer5.eq.'D') then
	write(88,*)'set ylabel "p(r)  [a.u.]" 0.0,0.0'
	write(88,*)'set xlabel ','''r''',' 0.0,0.0'
      endif 
      if(answer5.eq.'K') then
	write(88,*)'set ylabel "G(l)*l**4  [a.u.]" 0.0,0.0'
	write(88,*)'set xlabel ','''l''',' 0.0,0.0'
      endif
 
      if(answer5.eq.'B') then
	write(88,*)'set ylabel "Cross Section distr.  [a.u.]" 0.0,0.0'
	write(88,*)'set xlabel ','''r (cross section)''',' 0.0,0.0'
      endif
 
      if(answer5.eq.'C') then
	write(88,*)'set ylabel "Thickness distr.  [a.u.]" 0.0,0.0'
	write(88,*)'set xlabel ','''r (thickness)''',' 0.0,0.0'
      endif

      if(answer5.eq.'S') then
	write(88,*)'set ylabel "Sizedistribution  [a.u.]" 0.0,0.0'
      write(88,*)'set xlabel ','''R (number weighted) ''',' 0.0,0.0'
      endif       

         if(answerm.eq.'M'.or.answerm.eq.'Q') then
       	write(88,*)'set xrange [0:',1.05*dmax,']'
              else
       	write(88,*)'set xrange [0:',dmax,']'
         endif
	write(88,*)'set yrange [',0,':',1.05*fmax,']'
       if((answer5.eq.'K').or.(answer5.eq.'N').or.answer5.eq.'C'.
     -or.answerm.eq.'M') then
         if(fmin.ge.0.or.abs(fmin).le.0.005*fmax) then
         write(88,*)'set yrange [',0.,':*]'
         write(88,*)'set yrange [',0.,':',1.05*fmax,']'
         else
	write(88,*)'set yrange [*:*]'
       endif
       endif
	write(88,*)'set nokey'

      if(answer5.ne.'S') then
      write(88,574)aname,chi,alpha
 574  format(1x,'set title "p(r) ',
     -'',a,' chi = ',f7.3,' log(alpha) =',f6.1,'"')
      else
      write(88,5744)aname,chi,alpha
 5744 format(1x,'set title "Size ',
     -'',a,' chi = ',f7.3,' log(alpha) =',f6.1,'"')     
      endif
      if(answer5.eq.'K') then
      write(88,597)aname,chi,alpha
 597  format(1x,'set title "G(l)*l**4  ',
     -'',a,' chi = ',f7.3,' log(alpha) =',f6.1,'"')
      endif

      if(answer7.eq.'Y') 
     -write(88,576)"backgr = ",background,dmax*0.54,1.05*fmax*0.8
 576  format(1x'set label "',a,e9.3,'" at ',e8.2,',',e8.2)
      write(88,575)comment,dmax*0.72,1.05*fmax*0.9
 575  format(1x'set label "',a,'" at ',e8.2,',',e8.2)
      if((f(1).gt.0.85*f(2)).and.(answer7.ne.'Y')) 
     -write(88,577)"Fit backgr ?",dmax*0.72,1.05*fmax*0.8
 577  format(1x'set label "',a,'" at ',e8.2,',',e8.2)

c****************************

c	write(88,*)'set term gif medium size 640,480'
c	write(88,*)'set term gif medium size 500,360'
      write(88,818)screen

	write(88,*)'set output "fig5.gif"'
       if((answer5.eq.'K').or.(answer5.eq.'N').or.answer5.eq.'C') then
	write(88,898)bname,bname
 898  format
     -(1x,'plot ''',a,''' w e 1,''',a,''' w l 1, 0 w l 3')
       else
         if(answerm.eq.'M'.or.answerm.eq.'Q') then
	write(88,887)bname,bname,'prior.d'
 887  format(1x,'plot ''',a,''' w e 1,''',a,''' w l 1,''',a,''' w l 2,
     -0 w l 3')
         goto 8888
         endif
	write(88,888)bname,bname
 888  format(1x,'plot ''',a,''' w e 1,''',a,''' w l 1')
 8888 continue
       endif


      write(88,*)'exit gnuplot'

      close(88)

      return
      end
c****************************************************************************
c     This subroutine calculates the evidence for several values of alpha
c     (to provide the powell algorithm with a proper stating point)
c****************************************************************************
      subroutine alpha_e(p,ndtest)
      parameter (ndim=4)
      real p(ndim)
      character*1 answer,answer2
      common /data6/ alphaest,dest,etaest
      common /data8/ nfit,nerror,answer,answer2
      common /data10/ xprec,dotsp

      d1=0.5
      d2=2.0
      jmin=1
      jmax=12
      
      if(ndtest.eq.1) then
      d1=1.0
      d2=1.0
      jmin=8
      jmax=8
      endif

      fmin=1.e10
           
      do 10 nalphaestx=-6,28,2
      do 9 j=jmin,jmax
      destx=j*dest*0.125
      if(nfit.eq.4) p(nfit-2)=destx
      if((nfit.eq.3).and.(etaest.ne.0)) p(nfit-1)=destx
      if((nfit.eq.3).and.(etaest.eq.0)) p(nfit-2)=destx
      if(nfit.eq.2) p(nfit-1)=destx
      p(nfit)=nalphaestx
      evidence=func(p)
      xerror=1-5*(1-xprec)
       if(evidence.le.fmin.and.dotsp.ge.xprec) then
       alphamin=nalphaestx
       dmax=destx
       fmin=evidence
       endif
      nof=0

   9  continue

  10  continue

      ndtest=2
      

      write(6,*)
      if(fmin.lt.1.e10) then
      write(6,*)'Using alpha_est = ',alphamin+0.5,' and d_max = ',dmax
      write(6,*)
      else
      write(6,*)'Problem - no usable alpha found. 
     -Try entering other parameter values'
      write(6,*)

cccccccccccccc
      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig2.gif"'
       write(88,*)'set label "NO USABLE ALPHA/DIAMETER" at -3.5,0.55'
       write(88,*)
     -'set label "Try entering estimates" at -3.2,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot3.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig3.gif"'
       write(88,*)'set label "NO USABLE ALPHA/DIAMETER" at -3.5,0.55'
       write(88,*)
     -'set label "Try entering estimates" at -3.2,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot5.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig5.gif"'
       write(88,*)'set label "NO USABLE ALPHA/DIAMETER" at -3.5,0.45'
       write(88,*)
     -'set label "Try entering estimates" at -3.2,0.45'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
        close(88)
       write(6,*)'error5'
ccccccccccccccccccccccccccccc

      STOP
      
      endif
      p(nfit)=alphamin+0.5
      if(nfit.eq.4) p(nfit-2)=dmax
      if((nfit.eq.3).and.(etaest.ne.0)) p(nfit-1)=dmax
      if((nfit.eq.3).and.(etaest.eq.0)) p(nfit-2)=dmax
      if(nfit.eq.2) p(nfit-1)=dmax

      return
      end

      FUNCTION bessj0(x)
      REAL bessj0,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *-.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-1,
     *.1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      DATA r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,
     *651619640.7d0,-11214424.18d0,77392.33017d0,-184.9052456d0/,s1,s2,
     *s3,s4,s5,s6/57568490411.d0,1029532985.d0,9494680.718d0,
     *59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*
     *(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software H'U,]..
c************************************************************
c     In case of an error this subroutine is called:
c************************************************************
      subroutine error
      open(unit=88,file='plot2.pl',status='unknown')
	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig2.gif"'
       write(88,*)
     -'set label " Sorry - something went wrong ...     " at -4.8,0.5'
       write(88,*)
     -'set label "Perhaps empty spaces in the file name?" at -4.8,0.4'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
      close(88) 

      open(unit=88,file='plot5.pl',status='unknown')
 	write(88,*)'set nokey'
	write(88,*)'set term gif'
	write(88,*)'set output "fig5.gif"'
       write(88,*)
     -'set label " Sorry - something went wrong ...     " at -4.8,0.5'
       write(88,*)
     -'set label "Perhaps empty spaces in the file name?" at -4.8,0.4'
	write(88,*)'set yrange [0:1]'
	write(88,*)'plot 0'
       write(88,*)'exit gnuplot'
        close(88)
      write(6,*)'error6'
      stop
      return
      end
C*************************************************************************
C     Solution of the equation XX*EXP(XX)=BB*EXP(EE)=AA
C*************************************************************************
      FUNCTION XX(BB,EE)
      IF(EE.GT.50) GOTO 85
      if(ee.lt.-50) then
      xx=0
      goto 90
      endif
      AA=BB*EXP(EE)
      IF(AA.LE.0.004) XX=-.5+SQRT(.25+AA)
      IF(AA.LE.0.004) GOTO 90
      XO=LOG(AA+1.)
   80 FXO=XO*EXP(XO)
      XX=XO*(1-(FXO-AA)/(FXO*XO+1))
      IF(XX.LT.0) XX=0
      FX=XX*EXP(XX)
      DF=ABS(FX-AA)/AA
      IF(DF.LT.0.0001) GOTO 90
      XO=XX
      GOTO 80

   85 CORR=EE+LOG(BB)
      XX=CORR
   86 CORX=XX+LOG(XX)
      TEST=ABS((CORR-CORX)/CORR)
      IF(TEST.LT.0.0001) GOTO 90
      XO=XX
      XX=XO+(CORR-CORX)*(1+1./XO)
      GOTO 86

   90 RETURN
      END

