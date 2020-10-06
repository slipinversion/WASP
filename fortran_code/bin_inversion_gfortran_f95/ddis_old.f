      subroutine distaz(lat_sta,lon_sta,lat_e,lon_e,dis,azz,baz)
      parameter(const=6371*2*3.1415926/360)
      COMMON/SAVE/CO,C2,C4
      real lat_sta,lon_sta,lat_e,lon_e,dis,azz,baz
      dimension tH(2),PHI(2),XDEG(2),AZINV(2),DIST(2),az(2) 
      th(1)=lat_sta
      th(2)=lat_e
      phi(1)=lon_sta
      phi(2)=lon_e
      call disaz1(th, phi, 1, 2, xdeg, az, azinv, dist )
      azz=az(1)
      baz=azinv(1)
      dis=dist(1)
      gcarc=xdeg(1)
      return
      end	

      subroutine disaz1(th, phi, n, k, xdeg, az, azinv, dist )
      COMMON/SAVE/CO,C2,C4
      dimension tH(2),PHI(2),XDEG(2),AZINV(2),DIST(2)
         dimension  az(2)
      iERROR=0
      iERR=0
      RAD=6378.155
      cONSt=57.2957795
      nk = n+1
      do 99  i=1, nk
      th(i) = th(i)/CONST
  99  phi(i)= phi(i)/cONST
      GL = .99327733
      THG = atan(GL * tan(th(k)))
      c = sin(THG)
      XK= cos(tHG)
      a = xK * cos(phi(k))
      b = xK * sin(phi(k))
      ieRR=IERR+1
      if (IERROR.EQ. 1) write(*, 9000)IERR,N
9000  format( 2x, 'disaz to point number',2i5  )
      do 25 i=1, n
      dl = phi(i)-phi(k)
	THG=ATAN(GL*tan(th(i)))
	C1=SIN(THG)
	A1= cos(thg) * cos(phi(i))
	B1= cos(thg) * sin(phi(i))
      sc =  a * a1 + b * b1 + c * c1
      sd1= (a - a1)**2 + (b-b1)**2 + (c-c1)**2
      sd = (a + a1)**2 + (b+b1)**2 + (c+c1)**2
      sd = sqrt(sd1 * sd  / 4.0)
      xdeg(i) = atan(sd / sc) * CONST
	IERR=IERR+1
      if(IERROR.eq.1) write(*,9000)IERR,I
      if(sc) 1,2,2
  1   xdeg(i)=xdeg(i)+180.0
  2   ec=0.672267002e-2
      el=ec/(1.0-ec)
      E11=1.+el
	XX = (E11 + tan(th(i)) ** 2) / (E11 + tan(th(k)) ** 2)
	AL12 = tan(th(i)) / (E11 * tan(th(k))) + EC * sqrt(xx)
	DUM= sin(th(k)) * (AL12-cos(dl))
        A12=ATAN(sin(dl)/DUM)
	AL21 = tan(th(k)) /(E11 * tan(th(i))) + EC * sqrt(1.0 / xx)
	DUM1= sin(th(i)) * (AL21-cos(dl))
	A21=ATAN(-sin(dl)/DUM1)
       az(i) =A12*CONST
       azinv(i)=A21*CONST
	IERR=IERR+1
	if(IERROR.eq.1) write(*,9000) iERR,I
	if(sin(dl))  1005,1006,1006
 1005   if(duM)  1007,1008,1008
 1007  az(i) = az(i) - 180.
	goto 1008
 1006  if(duM) 1009,1008,1008
 1009  az(i)=180.+az(i)
 1008  if(-sin(dl)) 1010, 1011, 1011
 1010  if(dUM1) 1012,1013,1013
 1012  azinv(i)=azinv(i)-180.
     	goto 1013
 1011  if(DUM1 ) 1014,1013,1013
 1014  azinv(i)=azinv(i)+180.
 1013  cs12=cos(a12)
	IERR=IERR+1
       if(IERROR.eq.1) write(*,9000) iERR,I
	Eo= EL * ((cos(th(k)) * CS12) ** 2 + sin(th(k)) ** 2)
	co=1.+eo/4.0-3.*(eo**2)/64.+5.*(eo**3) /256.
	c2=-eo/8.0+(eo**2)/32.-15.*(eo**3) /1024.
	c4 = -(eo**2)/256.+3.*(eo**3) /1024.
	SE21=EC* sin(th(k)) ** 2
	V1=RAD/SQRT(1.-SE21)
	SE22 = EC * sin(th(i)) ** 2
	V2=RAD/SQRT(1.-SE22)
	Z1=V1*(1.-EC) * sin(th(k))
	Z2 = V2 * (1.0 - EC) * sin(th(i))
	X2 = V2 * cos(th(i)) * cos(dl)
	Y2 = V2 * cos(th(i)) * sin(dl)
	ARG= tan(th(k)) / CS12 / SQRT(1.+EO)
	U1=ATAN(ARG)
	X=(X2*CS12 - Y2 * sin(th(k)) * SIN(A12))*SQRT(1.+EO)
	ARG=(V1 * sin(th(k)) + (1.+EO)*(Z2-Z1))/X
	U2=ATAN(ARG)
	X = 1.0 + EL * (cos(th(k)) * CS12) ** 2
	BO=V1*SQRT(X)/(1.+EO)
	x=dd(BO,U1,U2)
        dist(i)=abs(X)
        CHCK=xdeg(i)*111.32
	TEST=CHCK-DIST(I)
	TEST=ABS(TEST)
            if(tEST.lt.30.0)  goto  302
  301  u2=u2+3.1415626535
	X=DD(BO,U1,U2)
            dist(i)=abs(X)
	TEST=CHCK-DIST(I)
	TEST=ABS(TEST)
	IF(TEST.LT.30.) GOTO 302
            u2=u2-(3.1415626535*2.0)
	X=DD(BO,U1,U2)
	dist(i)=abs(X)
	TEST=CHCK-DIST(I)
	TEST=ABS(TEST)
	IF(TEST.LT.30.) GOTO 302
	DIST(I)=CHCK
  302   th(i)=  th(i)*cONST
       phi(i)=phi(i)*cONST
         if(az(i).lt.0.0)   az(i)=360.0+az(i)
         if(azinv(i).lt.0.0)   azinv(i)=360.0+azinv(i)
   25   continue
      return
       end

              subroutine    vintrp(x,y,n,arg,val,ipoint)
       dimension    x(200),Y(200)
         do 10 i=1,n
        ipoint=i
       d=x(i)-arg
       if(d) 11,12,10
   12  val=y(i)
       goto 76
   11  slope=(y(i)-Y(i-1))/(x(i)-x(i-1))
       b=y(i)-slope*x(i)
       val=slope*arg+b
       goto 76
   10  continue
   76  return
       end
       subroutine amin(h,m,depth,imin)
       dimension h(1)
       arg1=depth-h(1)
       arg1=abs(arg1)
       imin=1
       if(arg1.eq.0.0) goto 76
       do 8 i=2,m
       arg2=depth-h(i)
       if(arg2.eq.0.0) goto 10
       arg2=abs(arg2)
       if(arg2.le.arg1) goto 9
       goto 8
   9   arg1=arg2
       imin=i
       goto 8
   10  imin=i
       goto 76
    8  continue
  76   return
       end
	FUNCTION DD(B,U,S)
	COMMON/SAVE/C0,C2,C4
	U2=U*2.
	SINU2=SIN(U2)
	U4=U*4.
	SINU4=SIN(U4)
	S2=2.*S
	SINS2=SIN(S2)
	S4=4.*S
	SINS4=SIN(S4)
	X=C0*(U-S)+C2*(SINU2-SINS2)+C4*(SINU4-SINS4)
	DD=B*X
	RETURN
	END
c       subroutine tictoc(tt,Ir,mn,ssec,ihr,imin,sec)
c       tt=tt/3600.0
c       tb=aint(tt)
c       Ir=tb+ihr
c       ta=(tt-tb)*60.0+imin
c       mn=aint(ta)
c       ssec=(ta-mn)*60.0+sec
c       if(ssec.lt.60.0) goto 307
c       ssec=ssec-60.0
c       mn=mn+1
c  307  if(mn.lt.60) goto 305
c       mn=mn-60
c       Ir=Ir+1
c  305  return
c       end
