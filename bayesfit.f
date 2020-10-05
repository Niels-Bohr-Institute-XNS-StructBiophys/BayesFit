c
c     Program Bayesfit
c     Version to upload at GitHub
c     Program used in Larsen, Arleth, Hansen (2018). J Appl Cryst, 51, 1151-1161.
c

      program bayesfit
      INTEGER ma,nca,ndata,ia(20),MMAX
      REAL alamda,chisq,fit(20),alpha(20,20),covar(20,20),
     *sd(5000),x(5000),y(5000),smear(20,20),w(20)
      real Ng,hessian2(20,20),nge,w2(20),dyda(20),ng2,ngmax
      real dummy1(20,20),dummy2(20,20),ch(20),al(20)
      character*13 par(20),model
      integer clock,clockold
      character*2 numberx
      character*30 fname
cc  comment out fit_*.d, f2name, unit 9
cc      character*34 f2name,f2file
      character*34 f2name
      character*37 f3name
      character*36 f4name
      character*41 f5name
      character*11 iname
      character*9 riname
      character*30 aname
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      common /parameter/par,nite,ite
      common /hess/ hessian(20,20),alpha2(20,20)
cc      common /diff2/diff2(20)
      common /model/ model
      common /start/istart
      call system_clock(clock)
      clockold=clock
      call getarg(1,aname)
      istart=0
cc comment out Qcheck.d, unit 77
cc      open(unit=77,file='Qcheck.d',status='unknown')
cc comment out global_par.d, unit 68
cc     open(unit=68,file='global_par.d',status='unknown')
cc comment out errors.d, unit 16
cc    open(unit=16,file='errors.d',status='unknown')
CC      write(68,*)
CC      close(68)
      open(unit=8,file=aname,status='old')
c     read from the inputfile (8)
      read(8,3)fname
      read(8,3)model
cc do not read in alphamin,alphamax and itot (number of alphas)
cc type in values instead
cc      read(8,*)alphamin,alphamax,itot,nite,maxit
      alphamin=-10
      alphamax=10
      itot=1
      read(8,*)nite,maxit
      close(8)
c     generate outputfilenames
cc  comment out fit_*.d, f2name, unit 9
cc      f2name='fit_'//fname
      f3name='output_'//fname
cc change name start_*.dat to prior.d
cc      f4name='start_'//fname
      f4name='prior.d'
cc change name parameters_*.d to parameteres.d
cc      f5name='parameters_'//fname
      f5name='parameters.d'
      open(unit=15,file=f5name,status='unknown')
      write(15,*)'# paramters from BayesFit'
cc      write(15,*)'# alpha     evidence   chi^2_r      N_g      alpha*S
cc     -  sum        chi^2     Iterations  mtot   Parameters --> '
cc      write(16,*)'# alpha     evidence   chi^2_r      N_g      alpha*S
cc     -  sum        chi^2     Iterations  mtot   Parameters --> '

cc  make new output file
cc      open(unit=88,file='params.d',status='unknown')
cc      write(88,*)'# Parameters from fit'

      itot=itot-1
      j1=0
  10  continue
      if(j1.eq.3) goto 16
      j1=j1+1
      tol=0.01
      ax=-10
      bx=0
      cx=10
      R=.61803399
      C=1.-R
      x0=ax
      x3=cx
      if(abs(cx-bx).gt.abs(bx-ax))then
        x1=bx
        x2=bx+C*(cx-bx)
      else
        x2=bx
        x1=bx-C*(bx-ax)
      endif
         alf=exp(x1)
         goto 2
 12      f1=ev
         j1=j1+1
         alf=exp(x2)
         goto 2
 13      f2=ev
         if(j1.ne.3) goto 16
 14       f2=ev   
          goto 16     
 15       f1=ev
 16      j1=3
1     if(abs(x3-x0).gt.tol*(abs(x1)+abs(x2)))then
        if(f2.lt.f1)then
          j2=1
          x0=x1
          x1=x2
          x2=R*x1+C*x3
          f1=f2
          alf=exp(x2)
          goto 2
        else
          j2=2
          x3=x2
          x2=x1
          x1=R*x2+C*x0
          f2=f1
          alf=exp(x1)
          goto 2
        endif
      goto 1
      endif
      if(f1.lt.f2)then
        golden=f1
        xmin=x1
      else
        golden=f2
        xmin=x2
      endif
        write(6,*)
        write(6,66)'Best evidence = ',ev,' found for alpha = ',alf
  66  format(1x,a,f8.3,x,a,e8.3) 
        write(6,*)
        goto 991
  2   mmax=100
      sumchi=0
  3   format(a)
      open(unit=8,file=aname,status='old')
      read(8,3)fname
      read(8,*)
      read(8,*)
      read(8,*)xmin,xmax
      mafit=0
      do 200 i=1,100
      read(8,*,end=202)fit(i),sigm(i),ia(i)
      fit(i)=fit(i)/sigm(i)
      if(ia(i).ne.0) mafit=mafit+1
      par(i)='Parameter'
  200 continue
  202 ma=i-1
      do 201 i=1,ma
      rp(i)=fit(i)
      sigm(i)=sigm(i)**2
  201 continue
      close(8)
c
c     ma :    Number of parameters
c
c     mafit : Number of parameters fitted
c
c     extra prior if mafit < ma
c
      j=0
      do 980 i=1,ma
      ch(i)=0.
      if(ia(i).ne.0) then
      j=j+1
      sigm2(j)=sigm(i)
      rp2(j)=rp(i)
      endif
 980  continue

      write(6,*)
      write(6,31)ma,' parameters read from ',aname
      write(6,*)
      write(6,32)'Input data file  :  ',fname
      write(6,32)'Output start I(q):  ',f4name
      write(6,32)'Output fit I(q)  :  ','fit.d'
      write(6,32)'Output parameters:  ',f5name
   31 format(x,i2,x,a,x,a)
   32 format(x,a,x,a)
      write(6,*)
cc      open(unit=9,file=f2name,status='unknown')
   22 continue
      open(unit=10,file=fname,status='old')
  20  continue
      open(unit=10,file=fname,status='old')
      do 120 i=1,5000
 110  continue
 125  read(10,*,end=111,err=110)x(i),y(i),sd(i)
      y(i)=y(i)
      sd(i)=sd(i)
      if(x(i).lt.xmin) goto 125
      if(y(i).eq.0) goto 125
      if(x(i).gt.xmax) goto 111
      if(sd(i).eq.0) sd(i)=1
  120 continue
  111 close(10)
      mtot=i-1
      open(unit=81,file=f4name,status='unknown')
      do 800 i=1,mtot
      call funcs2(x(i),fit,ymod,dyda,ma)
      write(81,*)x(i),ymod,0
  800 continue
      close(81)

  850 chisqo=1.e30
      chisq=1.e30
      nca=20
      alamda=-2
      write(6,*)'   chi^2/M           lambda             alpha 
     -   Q              dQ              S'
      do 400 kk=1,maxit
      call mrqmin(x,y,sd,mtot,fit,ia,ma,covar,alpha,nca,chisq,
     *alamda)

      kkk=mod(kk,4)+1
      ch(kkk)=chisq+s
      al(kkk)=alamda
      cht=(ch(1)+ch(2)+ch(3)+ch(4))/4.
      dq=abs(1-(chisq+s)/cht)

      if(dq.lt.1.e-4.and.kk.gt.5.and.al(kkk).eq.al(kkk-2)) then
      call random_number(rr)
      kkk=rr*ma+1
      write(6,*)'fit',kkk,' changed'
      fit(kkk)=fit(kkk)*1.00001
      endif

      write(6,*)chisq/mtot,alamda,alf,chisq+s,dq,s/alf
      if(alamda.ge.1.e6) goto 401
  400 continue
  401 write(6,*)
      write(6,402)
     -'End lambda = ',alamda,'File no =',jj+1,'Iterations =',kk-1
  402 format(1x,a,e10.2,5x,a,2x,i3,4x,a,2x,i3)
      call mrqmin(x,y,sd,mtot,fit,ia,ma,covar,alpha,nca,chisq,0.)

cc write to stdout     
      write(6,*)
      write(6,*)'                            FIT                       
     -     PRIOR      '
      do 600 k=1,ma
      sig=sqrt(sigm(k))
      write(6,601)par(k),fit(k)*sig,' +/- ',
     -sqrt(covar(k,k))*sig,rp(k)*sig,
     -' +/- ',sqrt(sigm(k))
  600 continue
  601 format(1x,a,3x,e10.4,x,a,x,e10.4,6x,e10.4,x,a,x,e10.4)
      write(6,*)

  605 chi2=0
      open(unit=21,file='fit.d',status='unknown')
      write(21,*)'#  fit for alpha = ',alf,' below: '
cc  comment out unit 9, fit_*.d, f2name
cc      write(9,*)'#  fit for alpha = ',alf,' below: '
      do 1000 i=1,mtot
      call funcs2(x(i),fit,ymod,dyda,ma)
      write(21,*)x(i),ymod,0
cc      write(9,*)x(i),ymod,0
      ymod2=exp(-(x(i)*fit(8))**2/2)+fit(2)
      chi2=chi2+(ymod-y(i))**2/sd(i)**2
 1000 continue
      close(21)
      do 1100 k=1,ma
 1100 continue
      sumchi=sumchi+chisq/mtot
  900 continue
c
c 3x inserts 3 blank spaces, x inserts 1 blank space etc. i5 for integer
c
  901 format(7(x,e10.4),3x,i5,3x,i5,3x,12(x,e10.4))
  902 format(12(x,e10.4))
  903 format(9(x,e10.4))
  904 format(10(x,e10.4))

c
c     alpha2 is (half the) curvature of Q
c
c     alpha2 is diagonalized below and the eigenvalues are w(ii)
c
      call SVDCMP(alpha2,mafit,mafit,20,20,W,smear)
c
c     The hessian is diagonalized below and the eigenvalues are w2(ii)
c
c     Rescaling of the Hessian is done here, because the matrix is destroyed in svdcmp
c
      do 8884 inew=1,mafit
      do 8883 jnew=1,mafit
      hessian2(inew,jnew)=hessian(inew,jnew)
 8883 continue
 8884 continue

      call SVDCMP(hessian,mafit,mafit,20,20,W2,smear)
c
c     als=alf/sigm2(ii) are the elements in the (diagonal) curvature of alf*S
c
c     In Hansen (2000) Eq. (17) 
c
c     log(evidence)=0.5*log det [curv(S)] + Q - 0.5*log det [curv(Q/alf)]
c
c     Below sum=log det [curv(alpha*S)/curv(Q)]
c
c  Below alf/sigm2(ii) are the diagonal elements of the curvature of alpha*S
c  w(ii) are the diagonal elements of (the diagonalized) curvature of Q = (chi^2/2 + alpha*S)
c  sum is (-log of) the occam factor =  volume (posterior parameter space) / volume(prior parameter space)
c  = prod eigenvaluse curv(alpha*S) / prod eigenvalues curv(Q) 

c
c     Calculation of the Occam facator (sum)
c
      sum=0
      ng=0
      do 985 ii=1,mafit
      w(ii)=abs(w(ii))
      w2(ii)=abs(w2(ii))
cxxx      sum=sum-log((alf/sigm2(ii))/w(ii))
      sum=sum-log(abs(alf)/w(ii))
      dng=w2(ii)/w(ii)
      ng=ng+dng
 985  continue

c************************************************************************************
c     Start - Calculation of number of good parameters Ng after rescaling of Hessian
c************************************************************************************
      call SVDCMP(hessian2,mafit,mafit,20,20,W2,smear)

      ng2=0
      do 986 ii=1,mafit
      dng=w2(ii)/(alf+w2(ii))
      ng2=ng2+dng
 986  continue
c**********************************************************************************
c     End - Calculation of number of good parameters Ng after rescaling of Hessian
c**********************************************************************************

c
c     ev is -log(evidence)
c
      Q=0.5*sumchi*mtot+s
      ev=0.5*sum + Q + log(alf)
    
cc change from outputting Chi2/M to outputting Chi2/(M-Ng) 
cc      write(6,987)
cc     -'Chi^2_r, -log(Evidence), Ng = ',sumchi,ev,ng2
cc 987  format(1x,a,f8.3,3x,f10.3,3x,f8.3)
cc      write(6,987)
cc     -'Chi^2_r, -log(Evidence), Ng, M ='    
cc     -,sumchi*mtot/(mtot-ng2),ev,ng2,mtot
 987  format(1x,a,f8.3,3x,f10.3,3x,f8.3,3x,i5)
      write(6,988) 'Chi2...........= ',sumchi*mtot
 988  format(1x,a,f15.4)
      write(6,988) '-log(Evidence).= ',ev
      write(6,988) 'Ng.............= ',ng2
      write(6,989) 'M..............= ',mtot
 989  format(1x,a,i15)
      write(6,988) 'Chi2/M.........= ',sumchi
      write(6,988) 'Chi2/(M-Ng)....= ',sumchi*mtot/(mtot-Ng)
cc write to new outputfile, file 88, params.d
       
c
c check the shape of the curves for Q and the Occam factor as a function of alpha
c  
cc comment out Qcheck.d, unit 77
cc      write(77,*)alf,Q,0.5*sum,0.5*sumchi*mtot,s,log(alf),s/alf

      ite=-1
      call funcs2(x(i),fit,ymod,dyda,ma)
      if(j1.eq.1) goto 12
      if(j1.eq.2) goto 13
      if(j1.eq.3.and.j2.eq.1) goto 14
      if(j1.eq.3.and.j2.eq.2) goto 15
 991  continue
      write(6,*)
cc      write(15,901)alf,ev,sumchi,ng2,s,sum,sumchi*mtot,
cc     - kk-1,mtot,(fit(k)*sqrt(sigm2(k)),k=1,ma)
cc      write(15,901)alf,ev,sumchi,ng2,s,sum,sumchi*mtot,
cc     - kk-1,mtot,(sqrt(covar(k,k))*sqrt(sigm2(k)),k=1,ma)
cc      write(15,901)alf,ev,sumchi,ng2,s,sum,sumchi*mtot,
cc     - kk-1,mtot,(sqrt(covar(k,k))*sqrt(sigm2(k)),k=1,ma)

      write(15,988) '# Chi2...........=',sumchi*mtot
      write(15,988) '# alpha......... =',alf
      write(15,988) '# S..............=',s
      write(15,988) '# alpha*S........=',alf*S
      write(15,988) '# -log(Evidence).=',ev
      write(15,989) '# M..............=',mtot
      write(15,988) '# Ng.............=',ng2
      write(15,988) '# Chi2/M.........=',sumchi
      write(15,988) '# Chi2/(M-Ng)....=',sumchi*mtot/(mtot-Ng)

cc write parameters to file parameters.d
      write(15,*)'                            FIT
     -     PRIOR      '
      do 602 k=1,ma
      sig=sqrt(sigm(k))
      write(15,601)par(k),fit(k)*sig,' +/- ',
     -sqrt(covar(k,k))*sig,rp(k)*sig,
     -' +/- ',sqrt(sigm(k))
  602 continue

      call system_clock(clock)
      cpu=(clock-clockold)*0.001
      write(*,955)cpu
  955 format(1x,'Cpu time used: ',f7.1,' seconds')
      close(15)
cc close new output file 88, params.d
cc      close(88)
cc  comment out errors.d, unit 16 (write errors to parameters_Isim.d
cc      close(16)
cc comment out Qcheck.d, unit 77
cc      close(77)
cc  comment out fit_*.d, unit 9, f2name
cc     close(9)
c      stop
      end
c****************************************************************************
      SUBROUTINE funcs2(x,a,y,dyda,na)
      INTEGER na
      REAL x,y,a(20),dyda(20)
      call fct(x,y,a,na)
c     stepsize
      h=0.0005
     
      do 11 i=1,na
      aold=a(i)
 5    a(i)=aold*(1+h)
      call fct(x,yp,a,na)
      a(i)=aold*(1-h)
      call fct(x,ym,a,na)
      a(i)=aold*(1+2*h)
      call fct(x,yp2,a,na)
      a(i)=aold*(1-2*h)
      call fct(x,ym2,a,na)
      a(i)=aold
      dyda(i)=(ym2-8*ym+8*yp-yp2)/(12.*h*aold)
 11   continue
      return
      end

      subroutine fct(x,y,a,na)
      real x,y,fit(20),head,N,a(20)
      character*13 par(20),model
      real psi,l,h,pi,I_q,q,h_t,h_h,h_p
      common /parameter/ par,nite,ite
      common /model/ model
      y=0
      if(model.eq.'micelle') call fct_mic(x,y,a,na)
      if(model.eq.'nanodisc') call fct_nano(x,y,a,na)
      if(model.eq.'coreshell') call fct_coreshell(x,y,a,na)
      if(model.eq.'test') call fct_test(x,y,a,na)
      if(y.eq.0) then
      write(6,*)'No model found'
      stop
      endif
      return
      end

      SUBROUTINE mrqmin(x,y,sig,ndata,a,ia,ma,covar,alpha,nca,chisq,
     *alamda)
      INTEGER ma,nca,ndata,ia(ma),MMAX
      REAL alamda,chisq,a(ma),alpha(nca,nca),covar(nca,nca),
     *sig(ndata),x(ndata),y(ndata)
      PARAMETER (MMAX=20)
CU    USES covsrt,gaussj,mrqcof
      INTEGER j,k,l,m,mfit
      REAL ochisq,atry(MMAX),beta(MMAX),da(MMAX)
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      SAVE ochisq,atry,beta,da,mfit,otchisq

      if(alamda.lt.0.)then
        mfit=0
        do 11 j=1,ma
          if (ia(j).ne.0) mfit=mfit+1
11      continue
        alamda=1.
        call mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq)
        ochisq=chisq
        otchisq=chisq+s
        do 12 j=1,ma
          atry(j)=a(j)
12      continue
      endif
      j=0
      do 14 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          k=0
          do 13 m=1,ma
            if(ia(m).ne.0) then
              k=k+1
              covar(j,k)=alpha(j,k)
            endif
13        continue
          covar(j,j)=alpha(j,j)*(1.+alamda)
          da(j)=beta(j)
        endif
14    continue
      call gaussj(covar,mfit,nca,da,1,1)
      if(alamda.eq.0.)then
        call covsrt(covar,nca,ma,ia,mfit)
        return
      endif
      j=0
      do 15 l=1,ma
        if(ia(l).ne.0) then
          j=j+1
          atry(l)=a(l)+da(j)*1.0
c ***************************************
c  START Limits for parameters
c  no more than 5 sigma from prior mean
c ***************************************
          rlimit=5.
          rmax=rp(l)+rlimit
          rmin=rp(l)-rlimit
          if(atry(l).gt.rmax) atry(l)=rmax
          if(atry(l).lt.rmin) atry(l)=rmin
        endif
        if(ia(l).eq.2) then
          if(atry(l).lt.0) atry(l)=rp(l)*0.1
        endif
15    continue
c ***************************************
c  END Limits for parameters
c ***************************************
      call mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq)
      tchisq=chisq+s*2
      if(tchisq.lt.otchisq)then
        alamda=0.1*alamda
        ochisq=chisq
        otchisq=tchisq
        j=0
        do 17 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
            k=0
            do 16 m=1,ma
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=covar(j,k)
              endif
16          continue
            beta(j)=da(j)
            a(l)=atry(l)
          endif
17      continue
      else
        alamda=10.*alamda
        chisq=ochisq
        tchisq=otchisq
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]4&93.
      SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
      INTEGER ma,mfit,npc,ia(ma)
      REAL covar(npc,npc)
      INTEGER i,j,k
      REAL swap
      do 12 i=mfit+1,ma
        do 11 j=1,i
          covar(i,j)=0.
          covar(j,i)=0.
11      continue
12    continue
      k=mfit
      do 15 j=ma,1,-1
        if(ia(j).ne.0)then
          do 13 i=1,ma
            swap=covar(i,k)
            covar(i,k)=covar(i,j)
            covar(i,j)=swap
13        continue
          do 14 i=1,ma
            swap=covar(k,i)
            covar(k,i)=covar(j,i)
            covar(j,i)=swap
14        continue
          k=k-1
        endif
15    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]4&93.
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=20)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                write(6,*) 'singular matrix in gaussj 1'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) write(6,*) 'singular matrix in gaussj 2'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END

      SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
      PARAMETER (NMAX=20)
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

C  (C) Copr. 1986-92 Numerical Recipes Software ]4&93.
      SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq)
      INTEGER ma,nalp,ndata,ia(ma),MMAX
      REAL chisq,a(ma),alpha(nalp,nalp),beta(ma),sig(ndata),x(ndata),
     *y(ndata)
      PARAMETER (MMAX=20)
      INTEGER mfit,i,j,k,l,m
      REAL dy,sig2i,wt,ymod,dyda(MMAX),a2(20)
      common /hess/ hessian(20,20),alpha2(20,20)
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      mfit=0
      do 11 j=1,ma
        if (ia(j).ne.0) mfit=mfit+1
11    continue
      do 13 j=1,mfit
        do 12 k=1,j
          alpha(j,k)=0.
          hessian(j,k)=0.
12      continue
        beta(j)=0.
13    continue
      chisq=0.
      do 16 i=1,ndata
        call funcs2(x(i),a,ymod,dyda,ma)
        sig2i=1./(sig(i)*sig(i))
        dy=y(i)-ymod
        j=0
        do 15 l=1,ma
          if(ia(l).ne.0) then
            j=j+1
c*******************************
            a2(j)=a(l)
c*******************************    
            wt=dyda(l)*sig2i
            k=0
            do 14 m=1,l
              if(ia(m).ne.0) then
                k=k+1
                alpha(j,k)=alpha(j,k)+wt*dyda(m)
                hessian(j,k)=hessian(j,k)+wt*dyda(m)
                alpha2(j,k)=alpha(j,k)
              endif
14          continue
            beta(j)=beta(j)+dy*wt
        if(i.eq.1) then
cxxx        alpha(j,j)=alpha(j,j)+alf*1./sigm2(j)
        alpha(j,j)=alpha(j,j)+alf
        alpha2(j,j)=alpha(j,j)
        beta(j)=beta(j)+alf*(rp2(j)-a2(j))
cxxx        beta(j)=beta(j)+alf*(rp2(j)-a2(j))/sigm2(j)
        endif
          endif
15      continue
        chisq=chisq+dy*dy*sig2i
16    continue  
      k=0
      s=0
      do 166 m=1,l
      if(ia(m).ne.0) then
      k=k+1
cxxx      s=s+alf*(rp2(k)-a2(k))**2/(2*sigm2(k))
      s=s+alf*(rp2(k)-a2(k))**2/(2.)
cccc      write(6,*)k,rp2(k),a2(k)

      endif
 166  continue
      do 18 j=2,mfit
        do 17 k=1,j-1
          alpha(k,j)=alpha(j,k)
          alpha2(k,j)=alpha(j,k)
          hessian(k,j)=hessian(j,k)
17      continue
18    continue

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ]4&93.

      function V(r,h)
      real pi
      pi = acos(-1.)
      V=pi*r**2*h
      return
      end

      function psi(q,alpha,r,h)
      x=q*r*sin(alpha)
      y=q*h*cos(alpha)/2.
      psi=2*(bessj1(x)/x)*sin(y)/y
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


      FUNCTION bessj1(x)
      REAL bessj1,x
      REAL ax,xx,z
      DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *s1,s2,s3,s4,s5,s6,y
      SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
     *s5,s6
      DATA r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,
     *242396853.1d0,-2972611.439d0,15704.48260d0,-30.16036606d0/,s1,s2,
     *s3,s4,s5,s6/144725228442.d0,2300535178.d0,18583304.74d0,
     *99447.43394d0,376.9991397d0,1.d0/
      DATA p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,
     *.2457520174d-5,-.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,
     *-.2002690873d-3,.8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+
     *y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*
     *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.,x)
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software H'U,]..

c******************************************************************************
c     Micelle model - elliptical - input_micelle.d 
c******************************************************************************  
      subroutine fct_mic(x,y,fit,na)
      character*13 par(20)
      real x,y,fit(20),N,N_A
      real psi,pi,I_q,q
      common /parameter/ par,nite,ite
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      pi=acos(-1.)
      N_A=6.022e23

      do 5 i=1,na
      fit(i)=fit(i)*sqrt(sigm(i))
  5   continue
c     
c     electron scattering length
c
      b_e=2.82e-13
      p_solvent=1./3
      q=x
      par(1)='Concentration'
      par(2)='Background   '
      par(3)='N            '
      par(4)='Contrast she  '
      par(5)='Contrast cor '
      par(6)='Sigma        '
      par(7)='Epsilon      '
      c=fit(1)
      B=fit(2)
      N=fit(3)
      rho_s=fit(4)*1.1
      rho_c=fit(5)*1.54
      sigma=fit(6)
      eps=fit(7)
      V_h=350
      V_tail=350.2
      CMC=0.17
      V_c=N*V_tail
      V_s=N*V_h
      R_c=(3*V_c/(eps*4*pi))**(1./3)
      cx=-(2+eps)*R_c/3.
      bx=(2*eps+1)*R_c**2/3.
      ax=cx**3-3*bx*cx/2.+3*V_s/(8.*pi)
      T=(ax+sqrt(abs(ax**2+(bx-cx**2)**3)))**(1./3.)+
     -  (ax-sqrt(abs(ax**2+(bx-cx**2)**3)))**(1./3.)+cx      
      V_t=4*pi*(eps*R_c+T)*(R_c+T)**2/3
      R_t=R_c+T
      eps_t=(eps*R_c+T)/(R_c+T)
      y=0
      imax=200
      dalpha=(pi/2.)/imax

         do 1 i=1,imax
         alpha=i*dalpha
         psi_s=(V_t*psi_ellip(q,R_t,eps_t,alpha)-
     -   V_c*psi_ellip(q,R_c,eps,alpha))/V_s
         psi_c=psi_ellip(q,R_c,eps,alpha)
         I_q=(rho_s*V_s*psi_s+rho_c*V_c*psi_c)**2*sin(alpha)
         y=y+I_q
    1    continue

      y=y*dalpha*b_e**2
      y=(c-CMC)/N*N_A*1.e-6*y*exp(-(q*sigma)**2/2)+B


      do 15 i=1,na
      fit(i)=fit(i)/sqrt(sigm(i))
  15  continue

      return
      end

      function psi_ellip(q,R_c,eps,alpha)
      r=R_c*sqrt(abs((sin(alpha))**2+(eps*cos(alpha))**2))     
      x=q*r
      psi_ellip=psi_sph(x)
      return
      end

      function psi_sph(x)
      psi_sph=3*(sin(x)-x*cos(x))/x**3
      return
      end

c******************************************************************************
c
c     Nano disc model - elliptical - input_nano.d
c     Nanodisc with MSP1D1 with His tags and DLPC as lipids
c
c     abbreviations:
c     _p = protein, _h = head, _t = tail, _l = lipid, _c = chain, _s = solvent
c
c******************************************************************************  
      subroutine fct_nano(x,y,fit,na)
      real x,y,fit(20),head,N,pro,N_A,global(20),n_w
      character*13 par(20)
      character*11 glob(20)
      real psi,l,h,pi,I_q,q,h_t,h_h,h_p
      common /parameter/ par,nite,ite
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      common /start/istart

      pi=acos(-1.)

      do 5 i=1,na
      fit(i)=fit(i)*sqrt(sigm(i))
  5   continue

      istart=istart+1
      if(istart.eq.1) then
      glob(1)="   alpha    "
      glob(2)="     a      "
      glob(3)="     b      "
      glob(4)="   V_core   "
      glob(5)="   V_shell  "
      glob(6)="   V_belt   "
      glob(7)="    h_p     "
      glob(8)="     H      "
      glob(9)="    h_t     "
      glob(10)="   rho_h    "
      glob(11)="   rho_t    "
      glob(12)="   rho_p    "
      glob(13)="   rho_c    "
cc  comment out global_par.d, unit 68
cc      open(unit=68,file='global_par.d',status='old',access='append')
cc      write(68,900)(glob(k),k=1,13)
  900 format(13a)
      endif

c     Volumes (V) and scattering lengths (b)
c     specific for the DLPC/MSP1D1 nanodisc with tags
      V_p = 54293.0
      V_c = 3143.0
      V_s = 30.0

      b_p = 23473.0
      b_h = 164.0723
      b_t = 178.0440
      b_c = 1.4250e+03
      b_s = 10.0

c     constants(avogadros number,electron scattering length)
      N_A = 6.022e23
      b_e = 2.82e-13

c     import fitting parameters
      Bg = fit(1)
      c = fit(2)
      V_l = fit(3)
      V_t = fit(4)
      CV_p = fit(5)
      N = fit(6)
      T = fit(7)
      sigma = fit(8)
      Ar = fit(9)
      eps = fit(10)
      n_w = fit(11)
      Rg = fit(12)

c     derived parameters
      V_h = V_l - V_t
      V_p = CV_p*V_p
      V_c = CV_p*V_c

c     add 7.9 water molecules per lipid headgroup (Kucerka et al., 2005)
      b_h = b_h + n_w * b_s
      V_h = V_h + n_w * V_s

c     reparametrization (from vol to scattering contrasts)
      p_s = b_s/V_s ! scattering length density of solvent
      dp_p = b_p/V_p - p_s
c      dp_c = dp_p
      dp_c = b_c/V_c - p_s
      dp_t = b_t/V_t - p_s
      dp_h = b_h/V_h - p_s

c     parameter labels
      par(1)='Background   '
      par(2)='Concentration'
      par(3)='V_lip     '
      par(4)='V_tail     '
      par(5)='CV_p       '
      par(6)='N            '
      par(7)='T            '
      par(8)='Sigma        '
      par(9)='Area         '
      par(10)='epsilon      '
      par(11)='n_w        '
      par(12)='Rg_tag       '

      ite=ite+1
      q=x
      xx=(q*Rg)**2
      P_c=(exp(-xx)+xx-1)/(xx/2.)
      P_tot=0
      F_tot=0

      b=sqrt(abs(N*Ar/(2*pi*eps)))
      a=eps*b

      jmax=nite
      kmax=nite
      dalpha=pi/(2*jmax)
      dfi=pi/(2*kmax)

      do 100 j=1,jmax
      alpha=j*dalpha
      do  50 k=1,kmax
      fi=k*dfi

      r_t=sqrt((a*sin(fi))**2+(b*cos(fi))**2)
      R=sqrt( ((a+T)*sin(fi))**2 +((b+T)*cos(fi))**2)
      h_p=V_p/(pi*((a+T)*(b+T)-a*b))

      h_t=2*V_t/Ar
      h_h=V_h/Ar
      H=h_t+2*h_h

      Reff=R+abs(Rg)
      yy=q*Reff*sin(alpha)
      ya=q*h_p*cos(alpha)/2

      psi_cc=(1-exp(-xx))/xx*bessj0(yy)*sin(ya)/ya

      tail=psi(q,alpha,r_t,h_t)

      pro=V(R,h_p)*psi(q,alpha,R,h_p)-V(r_t,h_p)*psi(q,alpha,r_t,h_p)
      pro=pro/(V(R,h_p)-V(r_t,h_p))

      head=(H*psi(q,alpha,r_t,H)-h_t*psi(q,alpha,r_t,h_t))/(2*h_h)

      V_nd=N*(V_t+V_h)+V_p
      dp_nd=(dp_t*N*V_t+dp_h*N*V_h+dp_p*V_p)/V_nd

      psi_nd=dp_t*(N*V_t)*tail+dp_h*(N*V_h)*head+dp_p*V_p*pro

      psi_nd=psi_nd/(dp_nd*V_nd)

      S_nd=psi_nd**2
      S_nd_c=psi_cc*psi_nd
      S_cc=psi_cc**2

      F=(dp_nd*V_nd)**2*S_nd+4*dp_c*V_c*dp_nd*V_nd*S_nd_c+
     -2*(dp_c*V_c)**2*(S_cc+P_c)
      F_tot=F_tot+F*sin(alpha)
   50 continue
  100 continue
      F_tot=F_tot*dalpha*dfi*(2./pi)

      if(ite.eq.0) then 
cc  comment out global_par.d, unit 68
cc      open(unit=68,file='global_par.d',status='old',access='append')
      global(1)=alf
      global(2)=a
      global(3)=b
      global(4)=V_t*N
      global(5)=N*V_h
      global(6)=V_p
      global(7)=h_p
      global(8)=H
      global(9)=h_t
      global(10)=dp_h
      global(11)=dp_t
      global(12)=dp_p
      global(13)=dp_c
cc      write(68,901)(global(k),k=1,13)
  901 format(13(x,e10.4))
cc      close(68)
      endif

      V_tot=V_nd+2*V_c
      dp_tot=(V_nd*dp_nd+2*V_c*dp_c)/V_tot

      P_tot=F_tot/(dp_tot*V_tot)**2

      y=c*1.e-9*N_A*F_tot*exp(-q**2*sigma**2)*b_e**2+Bg

      do 15 i=1,na
      fit(i)=fit(i)/sqrt(sigm(i))
  15  continue


 999  return
      end
c******************************************************************************
c     Core Shell model - input_coreshell.d
c******************************************************************************
      subroutine fct_coreshell(x,y,fit,na)
      character*13 par(20)
      real x,y,fit(20)
      real pi
      common /parameter/ par,nite,ite
      common /prior/ rp(20),sigm(20),alf,s,sigm2(20),rp2(20)
      pi=acos(-1.)
      
c     change from unitless params to params with units      
      do 5 i=1,na
      fit(i)=fit(i)*sqrt(sigm(i))
  5   continue


c     parameter labels
      par(1)='I0'
      par(2)='Background'
      par(3)='r'
      par(4)='R'
      par(5)='rho_r'
      par(6)='rho_R'

c     import fitting parameters
      c=fit(1)
      b=fit(2)
      r1=fit(3)
      r2=fit(4)
      p1=fit(5)
      p2=fit(6)
      
c     debugging  
c      print *, 'c', c, 'b', b, 'r1', r1, 'r2', r2, 'p1', p1, 'p2', p2

c     derived parameters
      qr1=x*r1
      qr2=x*r2
      v1=4.*pi*r1**3/3.
      v2=4.*pi*r2**3/3.
      s1=psi_sph(qr1)
      s2=psi_sph(qr2)

c     form factor sphere
c      y=c*(3*(sin(qr)-qr*cos(qr))/qr**3)**2+b
c      y=c*psi_sph(qr1)**2+b
      sum_psi=p1*v1*s1+p2*v2*s2-p2*v1*s1
      norm=p1*v1+p2*v2-p2*v1
      psi_coreshell=sum_psi/norm
      p_coreshell=psi_coreshell**2
c      y=p_coreshell
c      y=c*p_coreshell
      y=c*p_coreshell+b

c     make parameters unitless for the rest of the program
      do 15 i=1,na
      fit(i)=fit(i)/sqrt(sigm(i))
  15  continue

      return 
      end 

c******************************************************************************
c     Test model - simple sphere with background 
c******************************************************************************  
      subroutine fct_test(x,y,fit,na)
      character*13 par(20)
      real x,y,fit(20)
      real pi
      common /parameter/ par,nite,ite
      pi=acos(-1.)

      par(1)='Concentration'
      par(2)='Background   '
      par(3)='Radius            '

      c=fit(1)
      b=fit(2)
      r=fit(3)

      qr=x*r
      y=c*(3*(sin(qr)-qr*cos(qr))/qr**3)**2+b

      return
      end
