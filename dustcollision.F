c Calculation of dust size distribution as in Birnstiel+2011

       subroutine dustcollision(alphat,sigma,sigmad,
     %                   tmgas,amean,rad,rsnow)

        use dust

c Inputs:
c alpha - turbulence parameter alpha
c sigma - local surface density of gas
c sigmad - local surface density of dust
c tmgas - local gas temperature

c Outputs:
c amean - mean dust size 
c dmfrac - mass fraction in a given dust size, defined in module 
c         DUST and recomputed here when size distribution changes.

        implicit none
        doubleprecision sigma,tvf,tmgas,dista,alpha,dummy,sigmad
        
C tvf = largest turbulent velocity in cm s-1
C dgr = density in g cm-3 of dust material

        doubleprecision cs,Re,uf,V,J,term1,term2,dgr,c1,rad,rsnow
        doubleprecision acfmin,acfmax,a12,ugas,eps,dlna,alphat
        integer i,l,na,n,iL,iP,iR
        doubleprecision x,y,asett,abt,aL,aP,aR,amean,adrift
        doubleprecision abounce
        doubleprecision  a(200),u(200),ue(200),st(200)
        doubleprecision delta(nbin),const(nbin)

         

        alpha=alphat ! Viscous alpha = input alphat from alpha.dat
        alpha=1.0e-4 !  new turbulent alpha != viscous alpha

        dgr=drho(1)
C Sound speed
        cs=sqrt(1.38d-16*tmgas/(2.3*1.67d-24))
c Turbulent velocity
        tvf=sqrt(alpha)*cs
C Fragmentation threshold or critical collision velocity in cm/s
c  cannot exceed the largest turbulent relative velocity as growth
c is considered to be limited by fragmentation.
        uf=min(1000.0,tvf)
        if(rad.le.rsnow) uf=min(100.0,tvf) 
C Relative gas veocity as felt by grain
        ugas=tvf*sqrt(1.5d0) 
C Re number
        Re=alpha*sigma*2.0d-15/(2.0*2.3*1.67d-24)
C Settling size
       asett=2.0*alpha*sigma/(pi*dgr)
C Brownian motion - turbulence critical size
       aBT=((8.0*sigma*Re**(-0.25)/(pi*dgr))*(2.3*1.67d-24/
     +    (3.0*alpha*pi))**0.5*(4.0*pi*dgr/3.0)**(-0.5))**0.4
C Calculation of acfmax
        acfmax=2.0*sigma*uf**2.0/(3.1416*alpha*dgr*cs*cs)
        a12=2.0*sigma/(1.6d0*3.1416*dgr*sqrt(Re))
        a12=max(a12,5.0d-4)
        acfmax=max(acfmax,5.0*a12)
c Calculation of the grain size distribution, relative
c velocities 
        acfmin=0.025d-4 ! Set by Til's simulations


c Make a fine size distribution to calculate thresholds for various
c regimes

C       For a_i+1/a_i = 1.12 or less....
        na=1+(dlog10(acfmax)-dlog10(acfmin))/dlog10(1.11d0)
        na=min(na,200)
        do n=1,na
             a(n)=acfmin*1.11**(n-1)
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
             if(u(n).lt.(3.69e-4/a(n))**(5./3.)) abounce=a(n)
             if(u(n).le.(0.8*uf)) iL=n
             if(ue(n).le.(0.8*uf)) iP=n
             if(ue(n).le.uf) iR=n
        enddo

        aP=a(iP)

C J and V factors
        V=cs*(8.0*2.3*1.67d-24*sigma/(alpha*2.0d-15))**0.25*
     +      (3.0*alpha/(4.0*sigma*1.6))**0.5
        J=(2.5**(-9.0)+(1.1**9.0+(1.0+2.0*3.0**0.5*V/uf)**9)
     +        **(-1.0))**(-1./9.)          

c  Change all dust sizes to microns
         aBT=aBT*1.0d4
         a12=a12*1.0d4
         asett=asett*1.0d4
         aP=aP*1.0d4
         abounce=abounce*1.0d4


C Compute mass fraction in a grain size bin
C Compute mean grain cross-sectional size  
       c1=1.0d0

      do l=1,nds
        term1=0.0d0
        do i=1,nbin
c          dlna=dlog(adust(l,i+1))-dlog(adust(l,i))
c          dlna=1.0
          if(asett<aBT) then
           if(adust(l,i).lt.asett) then
            delta(i)=1.5d0
            const(i)=c1
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.aBT) then
            delta(i)=1.25d0
            const(i)=c1*asett**0.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.a12) then
            delta(i)=0.0d0
            const(i)=c1*asett**0.25d0*aBT**1.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.aP) then
            delta(i)=0.25d0
            const(i)=c1*asett**0.25d0*aBT**1.25d0*a12**(-0.25d0)/J
            term1=term1+const(i)*adust(l,i)**delta(i)
           else
            const(i)=0.0d0
            delta(i)=-100.0d0
            term1=term1+0.0d0
           endif
         elseif(asett<a12) then
           if(adust(l,i).lt.aBT) then
            delta(i)=1.5d0
            const(i)=c1
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.asett) then
            delta(i)=0.25d0
            const(i)=c1*aBT**1.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.a12) then
            delta(i)=0.0d0
            const(i)=c1*asett**0.25d0*aBT**1.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.aP) then
            delta(i)=0.25d0
            const(i)=c1*asett**0.25d0*aBT**1.25d0*a12**(-0.25d0)/J
            term1=term1+const(i)*adust(l,i)**delta(i)
           else
            const(i)=0.0d0
            delta(i)=-100.0d0
            term1=term1+0.0d0
           endif
         elseif(asett<aP) then
           if(adust(l,i).lt.aBT) then
            delta(i)=1.5d0
            const(i)=c1
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.a12) then
            delta(i)=0.25d0
            const(i)=c1*aBT**1.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.asett) then
            delta(i)=0.5d0
            const(i)=c1*a12**(-0.25d0)*aBT**1.25d0/J
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.aP) then
            delta(i)=0.25d0
            const(i)=c1*asett**0.25d0*aBT**1.25d0*a12**(-0.25d0)/J
            term1=term1+const(i)*adust(l,i)**delta(i)
           else
            const(i)=0.0d0
            delta(i)=-100.0d0
            term1=term1+0.0d0
           endif
        else
           if(adust(l,i).lt.aBT) then
            delta(i)=1.5d0
            const(i)=c1
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.a12) then
            delta(i)=0.25d0
            const(i)=c1*aBT**1.25d0
            term1=term1+const(i)*adust(l,i)**delta(i)
           elseif(adust(l,i).lt.aP) then
            delta(i)=0.5d0
            const(i)=c1*a12**(-0.25d0)*aBT**1.25d0/J
            term1=term1+const(i)*adust(l,i)**delta(i)
           else
            const(i)=0.0d0
            delta(i)=-100.0d0
            term1=term1+0.0d0
           endif
         endif
         enddo
c         delta(nbin)=delta(nbin-1)
       enddo
 
c      Renormalize dmfrac for accuracy
       amean=0.0d0
       term2=0.0d0
       dummy=0.0d0
       do l=1,nds
         do i=1,nbin
           if(delta(i).gt.-99) then
            dmfrac(l,i)=const(i)*adust(l,i)**delta(i)/term1
c            dmfrac(l,i)=const(i)*adust(l,i)**delta(i)
           else
            dmfrac(l,i)=0.0d0
           endif
           dummy=dummy+dmfrac(l,i)
           amean=amean+dmfrac(l,i)
           term2=term2+dmfrac(l,i)/adust(l,i)
         enddo
       enddo
       amean=1.0d-4/term2 ! Back to cm for use in GasT6

c     divide by dummy which should have been 1
        do l=1,nds
          do i=1,nbin
            dmfrac(l,i)=dmfrac(l,i)*1.0d0/dummy
          enddo
        enddo

c      Remove all grains larger than adrift for drift-limited size.
c    adrift from Birnstiel K E 2012, and using gamma=3 and N=0.5
c     Eq(20) for Stokes used and converted to a size 
c     mstar not passed here, so need to update if changed!
c      vkep=sqrt(GM*/r) = 1.149e13/sqrt(r)
c      convert adrift to microns. 
       adrift=1.0d4*2.0*sigma*uf*(1.149e13/sqrt(rad))/
     %  (pi*dgr*cs*cs*3.0*0.5) 
c     Another criterion from Eq.18 of BKE12 using gamma=3
        dummy=1.0d4*2.0*sigmad*(1.149e13*1.149e13/rad)/
     %  (pi*dgr*cs*cs*3.0)
        adrift=min(adrift,dummy)
c       With bouncing criterion, choose the smaller of adrift and abounce
c       only bouncing inside snowline then use for rad<rsnow
        !if(rad.lt.rsnow) adrift=min(adrift,abounce)


       do i=nbin,2,-1
c         limit dust max size to adrift if drift important
         if(adust(1,i).gt.adrift) then
            dmfrac(1,i-1)=dmfrac(1,i-1)+dmfrac(1,i)
            dmfrac(1,i)=0.0
         else
           exit
         endif
       enddo

c      Find amean all over again
       term2=0.0d0
       do i=1,nbin
           term2=term2+dmfrac(1,i)/adust(1,i)
       enddo
       amean=1.0d-4/term2

        if(dmfrac(1,1).eq.0.0.and.dmfrac(1,nbin).eq.0.0) then
          print*,'amean,alpha,sigma,tmgas'
          print*,amean,alpha,sigma,tmgas
         read*
        endif 

c      Some checks to conserve mass
c      Is all mass in grains smaller than amin? 
         if(amean*1.0d4.lt.adust(1,1)) dmfrac(1,1)=1.0 
c      Is all mass in grains larger than amax?
         if(amean*1.0d4.gt.adust(1,nbin)) dmfrac(1,nbin)=1.0 
      

        return 
        end
        

