c%----------------------------------------------%	
c			           	        |
c        Extreme Value Distribution  Fitting    |
c           Feb. 2005 Dengming Ming             |
c			           	        |
c%----------------------------------------------%
	Program NONLINEAR FITTING
	implicit none
C 
C      Driver for DNLS1E example.
C 
	CHARACTER*200 DATAFILE,file_fitting,flag_showfit*100,infocalc*100
	INTEGER       calcinfo,NDATA,NSEC,M,N           !M: number of data  N: number of variables
	INTEGER       I,IOPT,M0,N0,NPRINT,JNFO,INFO,INFO2,LWA
	INTEGER       IW(3)
	parameter (M0=10000,N0=10,LWA=N0*(M0+5)+M0,NPRINT = 0)
	parameter (IOPT = 1)
	DOUBLE PRECISION TOL,FNORM,X(N0),FVEC(M0),WA(LWA)
	DOUBLE PRECISION ENORM,D1MACH,SERR,STANDERROR,RERR,RCORR
	real*8 cutperct,cutdpa,datamin,delta
	EXTERNAL FCN
C
	calcinfo=0
C
	TOL = SQRT(D1MACH(4))
	N = 2	
C	
	read(*,*) datafile,nsec,cutperct,flag_showfit
	open (unit=20,file=datafile,status='old',err=5000)
	close(20)

	if(flag_showfit(1:1) .eq. 's' .or. flag_showfit(1:1) .eq. 'S' .or. 
     +      flag_showfit(1:1) .eq. 'y' .or. flag_showfit(1:1) .eq. 'Y')
     +       read(*,*) file_fitting
	call ini_data(calcinfo,infocalc,datafile,ndata)
	if(calcinfo.lt.0) then
	   write(*,'(A)') infocalc
	   stop
	endif

	call preproccdata(ndata,nsec,M,X,datamin,delta)

	CALL DNLS1E(IOPT,M,N,X,FVEC,TOL,NPRINT,INFO,IW,WA,LWA)
	if(INFO .le. 0) then 
	   write(*,1002) -2
	   stop
	endif
	call findcutdpa(cutperct,datamin,delta,X,cutdpa,INFO2)
	if(INFO2 .le. 0) then 
	   write(*,1003) -1
	   stop
	endif
	call corrl(M,X,rcorr)

	FNORM = ENORM(M,FVEC)	!FVEC(I) = X(I) - FIT(I)
	SERR = STANDERROR(M,FVEC)

!	write(*,999) cutdpa,SERR,RCORR(M,X,FVEC),(X(i),i=1,2)
	write(*,999) cutdpa,SERR,RCORR,(X(i),i=1,2)
	write(*,*) "==========================================="
	write(*,*) "          FITTING RESULTS                  "
	write(*,*) "==========================================="
	WRITE (*,1000) INFO,(X(I),I=1,N),FNORM,SERR,delta
	WRITE (*,1001) 

	if(flag_showfit(1:1) .eq. 's' .or. flag_showfit(1:1) .eq. 'S' .or. 
     +      flag_showfit(1:1) .eq. 'y' .or. flag_showfit(1:1) .eq. 'Y')
     +       call printresults(file_fitting,X,M,FVEC)

	STOP
 999	FORMAT (2x,F10.5,1x,F12.7,1x,F10.7,3x,'DPA_CUT  SERROR  RERROR',/
     +            2x,F15.5,2x,F15.5,3x, 'FITTING PARAMETERs ',/)
 
 1000	FORMAT (2x,I15,20x, 'EXIT PATAMETER (POSITIVE: NORMAL EXIT) ',/,
     +      2x,2F15.7,5x,'FINAL APPROXIMATE SOLUTION',/,2x,F15.7,20x,
     +      'FINAL L2 NORM OF THE RESIDUALS',/,2x,F15.7,20x,
     +      'FINAL STAND ERROR',/,2x,F15.7,20x,'DELTA X')
 1001	FORMAT (//,'Extreme value fitting: ',//,
     +      10x, 'f(x) = exp(-(x-a)/b)*exp(-(exp(-(x-a)/b))/b',/)
 1002	FORMAT (2x,I15,20x, 'EXIT PATAMETER < 0, ERROR FITTING  ',//)
 1003	FORMAT (2x,I15,20x, 'ERROR INPUT FOR PERCENTAGE',//)

 5000	calcinfo=-9900
	infocalc='-9900 CAN NOT OPEN FITTING DATA'
	write(*,'(A)') infocalc
	return
	END

C
	SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
	implicit none
C
C       Fitting Extreme value distribution
C       f(x) = exp(-(x-x0)/b)*exp(-exp(-(x-x0)/b))/b
C
C     This is the form of the FCN routine if IOPT=3,
C     that is, if the user calculates the Jacobian row by row.
	INTEGER I,M,N,IFLAG,LDFJAC
	DOUBLE PRECISION X(N),FVEC(M),FJAC(N)
	DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
	integer n1,n2,nsec,id
	parameter (n1=10000)
	real*8  X1(n1),Y1(n1),extremeval
	common /DATA1/ X1,Y1
	external extremeval 
C       
	IF (IFLAG .NE. 0) GO TO 5
C       
C       Insert print statements here when NPRINT is positive.
C       
	RETURN
 5	CONTINUE

	IF( IFLAG.NE.1) GO TO 20
	do i=1,m
	   !fvec(i) = y1(i)  + (x1(i)-x(1))*x(2) + dexp(-(x1(i)-x(1))*x(2)) - dlog(x(2))
	   !fvec(i)=y1(i)- x(1) * x1(i) -x(2)

	   !fvec(i) = y1(i)-( exp(-(x1(i)-x(1))/x(2))*exp(-exp(-(x1(i)-x(1))/x(2)))/x(2)  )
	   fvec(i) = y1(i) - extremeval(x1(i),x(1),x(2))

	enddo

	RETURN
C       
C       Below, calculate the LDFJAC-th row of the Jacobian.
C       
 20	CONTINUE

	write(*,*) '-9999 DO_not_use_JACOBI'
	stop
	
c$$$	I = LDFJAC
c$$$	TMP1 = I
c$$$	TMP2 = 16 - I
c$$$	TMP3 = TMP1
c$$$	IF (I .GT. 8) TMP3 = TMP2
c$$$	TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
c$$$	FJAC(1) = -1.E0
c$$$	FJAC(2) = TMP1*TMP2/TMP4
c$$$	FJAC(3) = TMP1*TMP3/TMP4
	RETURN
	END

	real*8 Function extremeval(x,a,b)
	implicit none
	real*8 x,a,b,c
	
	extremeval=dexp(-(x-a)/b)*dexp(-dexp(-(x-a)/b))/b
!	c=1./b
!	extremeval=dexp(-(x-a)*c)*dexp(-dexp(-(x-a)*c)*c

	return
	end

	
	subroutine preproccdata(n,nsec,imax,xini,ymin0,delta)
	implicit none
	integer n,n0,n1,n2,nsec,imax,id,i,k
	parameter (n0=10000,n1=10000)
	real*8  y(n0),X1(n1),Y1(n1),yt(n1),xini(*),tmp
	real*8  ymin0,ymax0,ymin,ymax,delta
	real*8  dexp,dlog
	integer index(n1)
	common /DATA0/ Y
	common /DATA1/ X1,Y1

	
	do i=1,n1
	   y1(i)=0.
	   x1(i)=0.
	   yt(i)=0.
	enddo

	ymin=Y(1)
	ymax=y(1)
	do i=1,n
	   if(ymin .gt. y(i)) ymin=y(i)
	   if(ymax .lt. y(i)) ymax=y(i)
	enddo

	delta=(ymax-ymin)/nsec
	imax=1+int((ymax-ymin)/delta)
	ymin0=ymin
	ymax0=ymax
!assuming the second parameter as the distribution width
	xini(2)=(ymax-ymin)/5.
!calc. number of data in a specific section
	do i=1,n
	   id = 1+int((y(i)-ymin)/delta)
	   yt(id)=yt(id)+1	
!	   write(*,*) 'PRC1> ',id,yt(id),y(i),delta,int(y(i)-ymin)/delta
	enddo
!interpolate value for sections no sample y(id) values are found
!a average value from nearest up and down sections (have smple values) are used
	n2=n
 	do i=1,imax
	   if(yt(i) .le. 0) then
	      ymin=0.
	      ymax=0.
	      if(i .eq. 1 .or. i .eq. imax) then
		 goto 300
	      endif
	      do k=i-1,1,-1
		 ymin=yt(k)
		 if(ymin.gt.0)    goto 100
	      enddo
 100	      continue
	      do k=i+1,imax
		 ymax=yt(k)
		 if(ymax.gt.0) goto 200
	      enddo
 200	   continue
 	   yt(i)=0.5*(ymin+ymax)
	   n2=n2+yt(i)
	   endif
 300	   continue 
	enddo
!Normalization: sum{y(i)*delta_x(i)} = 1 
	do i=1,imax
	   y1(i)=yt(i)/n2/delta
	   x1(i)=ymin0 + (i-0.5)*delta
	enddo

!assuming the first parameter xini(1) as x_i, where y_i has an extreme
	ymax=y1(1)
	do i=1,n2
	   if(ymax < y1(i)) then 
	      ymax=y1(i)
	      xini(1)=x1(i)
	   endif
	enddo
c$$$
c$$$	tmp=0.
c$$$	do i=1,imax
c$$$	   tmp=tmp+delta*y1(i)
c$$$	   write(23,500) i,x1(i),y1(i),tmp
c$$$	enddo
 500	format(1x,I5,1x,F10.6,1x,F10.6,2x,F15.6,10x,"I  X(I)  Y(I)  SUM")
	return
	end
	subroutine ini_data(calcinfo,infocalc,datafile,ndata)
	implicit none
	character*(*) datafile,infocalc
	integer       calcinfo,ndata,ndata0,i
	parameter (ndata0=10000)
	real*8  y(ndata0),tmp
	common /DATA0/ Y

	open (unit=20,file=datafile,status='old',err=1000)
	ndata=0
	do while (ndata .ge. 0) 
	   read(20,*,end=100,err=500)  tmp
	   ndata=ndata+1
	   if(ndata.gt.10000) goto 2000
	   Y(ndata)=tmp
	enddo
 100	close(20) 
	return
c$$$	do i=1,ndata
c$$$	   write(11,'(A,2x,I5,5x,F10.5)') 'READ> ', i,Y(i) 
c$$$	enddo
 500	calcinfo=-9500
	infocalc='-9500 READING ERROR in FITTING DATA'
	return

 1000	calcinfo=-9900
	infocalc='-9900 CAN NOT OPEN FITTING DATA'
	return

 2000	calcinfo=-9000
	infocalc='-9000 NUMBER OF FITTING DATA > 10000 '
	return
C
	end

	subroutine findcutdpa(cutperct,datamin,delta,X,cutx1,INFO)
	implicit none
	integer INFO
	real*8 cutperct,cutx1,x(*),datamin,delta
	integer ndata,ndata0
	parameter (ndata0=10000)
	real*8  y(ndata0),tmp
	common /DATA0/ Y

	if (cutperct .le. 0. .or. cutperct .gt. 0.99999999) then
	   INFO = -9999
	   return
	endif

	cutx1 =  - x(2)*dlog(-(dlog(cutperct))) + x(1)

	
!	write(*,*) 'CUTX1>  ', cutx1,datamin,delta

	INFO = 1

	return
	end

	subroutine printresults(file_output,X,ntot,FVEC)
	implicit none
	character*(*) file_output
	integer n,n0,n1,ntot,i
	parameter (n0=10000,n1=10000)
	real*8  y(n0),X1(n1),Y1(n1),FVEC(*),fitted_value
	real*8  ymin,ymax,delta,X(*),extremeval
	real*8  dexp,dlog
	common /DATA0/ Y
	common /DATA1/ X1,Y1
	external extremeval
	
	open(unit=21,file=file_output,status='unknown',err=2000)
	write(21,2001)"#","X(i)","Y(i)","Fitted_Y(i)","DELTA"

	do i=1,ntot
	   fitted_value=extremeval(x1(i),x(1),x(2))
	   write(21,2002) x1(i),y1(i),fitted_value,fvec(i)
	enddo
	write(*,2003) 'Fitting results data see: extremefittingresult.dat'
	return

 2000	write(*,*) ' Can not open file for WRITE'
 2001	format(A1,A11,2x,3(A12,2x))
 2002   format(4(F12.7,2x)) 
 2003   format(A,//)
	return
	end


	real*8 Function rcorr1(n,X,FVEC)
	implicit none
	integer n,i,j,k
	real*8  X(*),FVEC(*),Y(n)
	real*8  avgx,avgy,sxx,syy,sxy

	do i=1,n
	   y(i)=x(i)-fvec(i)
	enddo
	avgx=0
	avgy=0
	do i=1,n
	   avgx=avgx+x(i)
	   avgy=avgy+y(i)
	enddo

	sxx=0
	sxy=0
	syy=0
	do i=1,n
	   sxx=sxx+(x(i)-avgx)**2
	   sxy=sxy+(x(i)-avgx)*(y(i)-avgy)
	   sxx=sxx+(y(i)-avgy)**2
	enddo

!	rseq=sxy*sxy/sxx/syy
	rcorr1=sqrt(sxy*sxy/sxx/syy)


	return
	end


	subroutine corrl(N,X,rcorr)
	implicit none
	INTEGER N,I,J,K
	real*8 x(*),avgx,avgy,sxx,sxy,syy,rcorr
	integer n1
	parameter (n1=10000)
	real*8  X1(n1),Y1(n1),FY1(n1),extremeval
	common /DATA1/ X1,Y1

	external extremeval 
	
	do i=1,n
	   fy1(i)=extremeval(x1(i),x(1),x(2))
	   !WRITE(19,*) X1(I),Y1(I),Y1(I)-FY1(I),FY1(I)
	enddo

	avgx=0.
	avgy=0.
	do i=1,n
	   avgx=avgx+y1(i)
	   avgy=avgy+fy1(i)
	enddo
	avgx=avgx/dfloat(n)
	avgy=avgy/dfloat(n)

	sxx=0.
	sxy=0.
	syy=0.
	do i=1,n
	   sxx=sxx+(y1(i)-avgx)**2
	   sxy=sxy+(y1(i)-avgx)*(fy1(i)-avgy)
	   syy=syy+(fy1(i)-avgy)**2
	enddo

!	rseq=sxy*sxy/sxx/syy
	rcorr=sqrt(sxy*sxy/sxx/syy)
!	write(*,*) 'RCORR> ',avgx,avgy
!	write(*,*) 'RCORR> ',rcorr,sxy,sxx,syy
c
	return
	end
