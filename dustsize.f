
        subroutine dustsize(alpha,sigma,tgas,amean,opred)

        implicit none
        
c#include "diskpars.h"
        doubleprecision sigma,tvf,dgr,alpha,tgas

C sigma = gas surface density in g cm-2
C tvf = largest turbulent velocity in cm s-1
C dgr = density in g cm-3 of dust material
C alpha = turbulent alpha
C tgas = gas temperature

        doubleprecision cs,Re,uf,V,J,term1,term2,opred
        doubleprecision amin,amax,a12,ugas,eps,da
        integer i,na,n,iL,iP,iR
        doubleprecision x,y,pi,asett,abt,aL,aP,aR,amean
        doubleprecision  a(200),u(200),ue(200),st(200)


        pi=acos(-1.0d0)
        tvf=10.0d0
        dgr=3.6d0

C Sound speed

        cs=sqrt(1.38d-16*tgas/(2.3*1.67d-24))

C Fragmentation threshold or critical collision velocity in cm/s
c  cannot exceed the largest turbulent relative velocity as growth
c is considered to be limited by fragmentation.

        uf=min(50.0,sqrt(alpha)*cs)
    

C Relative gas veocity as felt by grain
        ugas=cs*(1.5d0*alpha)**0.5

C Re number
        Re=alpha*sigma*2.0d-15/(2.0*2.3*1.67d-24)

C Settling size
       asett=2.0*alpha*sigma/(pi*dgr)

C Brownian motion - turbulence critical size
       aBT=((8.0*sigma*Re**(-0.25)/(pi*dgr))*(2.3*1.67d-24/
     +    (3.0*alpha*pi))**0.5*(4.0*pi*dgr/3.0)**(-0.5))**0.4

C Calculation of amax
        amax=2.0*sigma*uf**2.0/(3.1416*alpha*dgr*cs*cs)
        a12=2.0*sigma/(1.6d0*3.1416*dgr*sqrt(Re))
        a12=max(a12,5.0d-4)
        amax=max(amax,5.0*a12)
        

c Calculation of the grain size distribution, relative
c velocities 
        amin=0.025d-4 ! Set by Til's simulations

c  Simpler version to see what is happening
c  02/28/2013
        amean=sqrt(amin*amax)
        opred=1.0
        return

C       For a_i+1/a_i = 1.12 or less....
        na=1+(dlog10(amax)-dlog10(amin))/dlog10(1.11d0)
        do n=1,na
             a(n)=amin*1.11**(n-1)
C            Stokes parameter 
             st(n)=a(n)*dgr*pi/(2.0*sigma)
C            Eps factor
             eps=(a(n)-a12)/(4.0*a12)
C            Monomer collision velocity
             if(a(n).lt.a12) then
                 u(n)=ugas*Re**0.25*(st(n)-st(1))
                 ue(n)=0.0d0
             elseif(a(n).lt.(5.0*a12)) then
                 u(n)=(1.0-eps)*ugas*Re**0.25*(st(n)-st(1))
     +            + eps*(st(n)-st(1))*(3.0*st(n))**0.5
                 ue(n)=sqrt(2./3.)*u(n)
             else
                 u(n)=ugas*(3.0*st(n))**0.5
                 ue(n)=sqrt(2./3.)*u(n)
             endif
             if(u(n).le.(0.8*uf)) iL=n
             if(ue(n).le.(0.8*uf)) iP=n
             if(ue(n).le.uf) iR=n
        enddo

C J and V factors
        V=cs*(8.0*2.3*1.67d-24*sigma/(alpha*2.0d-15))**0.25*
     +      (3.0*alpha/(4.0*sigma*1.6))**0.5
        J=(2.5**(-9.0)+(1.1**9.0+(1.0+2.0*3.0**0.5*V/uf)**9)
     +        **(-1.0))**(-1./9.)          

C Compute mean grain cross-sectional size  
        term1=0.0d0
        term2=0.0d0
        aP=a(iP)
        do i=1,na-1
          da=a(i+1)-a(i)
          if(a(i).lt.asett) then
           if(a(i).lt.aBT) then
            term1=term1+a(i)**(1.5-2.0)*da
            term2=term2+a(i)**(1.5-4.0)*da
           elseif(a(i).lt.a12) then
            term1=term1+a(i)**(0.25-2.0)*da
            term2=term2+a(i)**(0.25-4.0)*da
           elseif(a(i).lt.aP) then
            term1=term1+a(i)**(0.5-2.0)*da
            term2=term2+a(i)**(0.5-4.0)*da
           else
              term1=term1+0.0
              term2=term2+0.0
           endif
          else
           if(a(i).lt.aBT) then
            term1=term1+a(i)**(1.25-2.0)*da
            term2=term2+a(i)**(1.25-4.0)*da
           elseif(a(i).lt.a12) then
            term1=term1+a(i)**(0.0-2.0)*da
            term2=term2+a(i)**(0.0-4.0)*da
           elseif(a(i).lt.aP) then
            term1=term1+a(i)**(0.25-2.0)*da
            term2=term2+a(i)**(0.25-4.0)*da
           else
              term1=term1+0.0
              term2=term2+0.0
           endif
          endif
         enddo
         amean= dsqrt(term1/term2)
         opred=1.0

 
        return 
        end
        

