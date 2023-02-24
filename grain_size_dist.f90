!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        program grainsizes
        use dust
        implicit none
        integer i
        double precision sigma,tgas,t,alpha,a,r

        alpha=1.0e-2        
        sigma=1000.0
        tgas=80.0
        r=10.0*1.5e13 

        call readcomp
        call dustsize

        call dustcollision(alpha,sigma,tgas,a,r)

        do i=1,nbin
!          write(*,'(2ES12.3)') adust(1,i),dmfrac(1,i)
        enddo

        end
        
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Calculation of dust size distribution as in Birnstiel+2011

        subroutine dustcollision(alphat,sigma,tmgas,amean,rad)

        use dust

! Inputs:
! alpha - turbulence parameter alpha
! sigma - local surface density of gas
! tmgas - local gas temperature

! Outputs:
! amean - mean dust size 
! dmfrac - mass fraction in a given dust size, defined in module 
!         DUST and recomputed here when size distribution changes.

        implicit none
        doubleprecision sigma,tvf,tmgas,dista,alpha,dummy
        
! tvf = largest turbulent velocity in cm s-1
! dgr = density in g cm-3 of dust material

        doubleprecision cs,Re,uf,V,J,term1,term2,dgr,c1,rad
        doubleprecision acfmin,acfmax,a12,ugas,eps,dlna,alphat
        integer i,l,na,n,iL,iP,iR
        doubleprecision x,y,asett,abt,aL,aP,aR,amean,adrift
        doubleprecision  a(200),u(200),ue(200),st(200)
        doubleprecision delta(nbin),const(nbin)

         

         alpha=alphat ! Viscous alpha = turbulent alpha
         !alpha=1.0e-3 !  new turbulent alpha != viscous alpha

         dgr=3.0! drho(1)
! Sound speed
        cs=sqrt(1.38d-16*tmgas/(2.3*1.67d-24))
! Turbulent velocity
        tvf=sqrt(alpha)*cs
! Fragmentation threshold or critical collision velocity in cm/s
!  cannot exceed the largest turbulent relative velocity as growth
! is considered to be limited by fragmentation.
        uf=min(100.0,tvf)
! Relative gas veocity as felt by grain
        ugas=tvf*sqrt(1.5d0) 
! Re number
        Re=alpha*sigma*2.0d-15/(2.0*2.3*1.67d-24)
! Settling size
       asett=2.0*alpha*sigma/(pi*dgr)
! Brownian motion - turbulence critical size
       aBT=((8.0*sigma*Re**(-0.25)/(pi*dgr))*(2.3*1.67d-24/ &
         (3.0*alpha*pi))**0.5*(4.0*pi*dgr/3.0)**(-0.5))**0.4
! Calculation of acfmax
        acfmax=2.0*sigma*uf**2.0/(3.1416*alpha*dgr*cs*cs)
        a12=2.0*sigma/(1.6d0*3.1416*dgr*sqrt(Re))
        a12=max(a12,5.0d-4)
        acfmax=max(acfmax,5.0*a12)
! Calculation of the grain size distribution, relative
! velocities 
        acfmin=0.025d-4 ! Set by Til's simulations

! Make a fine size distribution to calculate thresholds for various
! regimes

!       For a_i+1/a_i = 1.12 or less....
        na=1+(dlog10(acfmax)-dlog10(acfmin))/dlog10(1.11d0)
        do n=1,na
             a(n)=acfmin*1.11**(n-1)
!            Stokes parameter 
             st(n)=a(n)*dgr*pi/(2.0*sigma)
!            Eps factor
             eps=(a(n)-a12)/(4.0*a12)
!            Monomer collision velocity
             if(a(n).lt.a12) then
                 u(n)=ugas*Re**0.25*(st(n)-st(1))
                 ue(n)=0.0d0
             elseif(a(n).lt.(5.0*a12)) then
                 u(n)=(1.0-eps)*ugas*Re**0.25*(st(n)-st(1)) &
                 + eps*(st(n)-st(1))*(3.0*st(n))**0.5
                 ue(n)=sqrt(2./3.)*u(n)
             else
                 u(n)=ugas*(3.0*st(n))**0.5
                 ue(n)=sqrt(2./3.)*u(n)
             endif
             if(u(n).le.(0.8*uf)) iL=n
             if(ue(n).le.(0.8*uf)) iP=n
             if(ue(n).le.uf) iR=n
             V=(3.69e-4/a(n))**(1.0/0.6)
             print*,n,a(n)*1.0e4,u(n),V
        enddo
        aP=a(iP)

! J and V factors
        V=cs*(8.0*2.3*1.67d-24*sigma/(alpha*2.0d-15))**0.25* &
           (3.0*alpha/(4.0*sigma*1.6))**0.5
        J=(2.5**(-9.0)+(1.1**9.0+(1.0+2.0*3.0**0.5*V/uf)**9) &
             **(-1.0))**(-1./9.)          

!  Change all dust sizes to microns
         aBT=aBT*1.0d4
         a12=a12*1.0d4
         asett=asett*1.0d4
         aP=aP*1.0d4


! Compute mass fraction in a grain size bin
! Compute mean grain cross-sectional size  
       c1=1.0d0

      do l=1,nds
        term1=0.0d0
        do i=1,nbin
!          dlna=dlog(adust(l,i+1))-dlog(adust(l,i))
!          dlna=1.0
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
!         delta(nbin)=delta(nbin-1)
       enddo
 
!      Renormalize dmfrac for accuracy
       amean=0.0d0
       term2=0.0d0
       dummy=0.0d0
       do l=1,nds
         do i=1,nbin
           if(delta(i).gt.-99) then
            dmfrac(l,i)=const(i)*adust(l,i)**delta(i)/term1
!            dmfrac(l,i)=const(i)*adust(l,i)**delta(i)
           else
            dmfrac(l,i)=0.0d0
           endif
           dummy=dummy+dmfrac(l,i)
           amean=amean+dmfrac(l,i)
           term2=term2+dmfrac(l,i)/adust(l,i)
         enddo
       enddo
       amean=1.0d-4/term2 ! Back to cm for use in GasT6

!     divide by dummy which should have been 1
        do l=1,nds
          do i=1,nbin
            dmfrac(l,i)=dmfrac(l,i)*1.0d0/dummy
          enddo
        enddo

!      Remove all grains larger than adrift for drift-limited size.
!    adrift from Birnstiel K E 2012, and using gamma=3 and N=0.5
!     mstar not passed here, so need to update if changed!
!      vkep=sqrt(GM*/r) = 1.149e13/sqrt(r)
       adrift=1.0d4*2.0*sigma*uf*(1.149e13/sqrt(rad))/ &
       (pi*dgr*cs*cs*3.0*0.5)



       do i=nbin,2,-1
         if(adust(1,i).gt.adrift) then
            dmfrac(1,i-1)=dmfrac(1,i-1)+dmfrac(1,i)
            dmfrac(1,i)=0.0
         else
           exit
         endif
       enddo

!      Find amean all over again
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

!      Some checks to conserve mass
!      Is all mass in grains smaller than amin? 
         if(amean*1.0d4.lt.adust(1,1)) dmfrac(1,1)=1.0 
!      Is all mass in grains larger than amax?
         if(amean*1.0d4.gt.adust(1,nbin)) dmfrac(1,nbin)=1.0 
      

        return 
        end
        

MODULE constants      
  implicit none
  integer, parameter :: DP = SELECTED_REAL_KIND(14)
  real(KIND=DP), parameter :: pi  = 3.14159265358979324
  real(KIND=DP), parameter :: kB  = 1.3807d-16 ! Boltz. constant [erg/K]
  real(KIND=DP), parameter :: mp  = 1.6726d-24 ! Mass of proton [g]
  real(KIND=DP), parameter :: GG  = 6.672d-8 ! Gravitational constant
  real(KIND=DP), parameter :: ss = 5.6703d-5 ! Stefan-Boltz const  [erg/cm^2/K^4/s]
  real(KIND=DP), parameter :: au=1.496d13    ! Astronomical unit [cm]
  real(KIND=DP), parameter :: solm=1.99d33   ! Solar mass [g]
  real(KIND=DP), parameter :: solr=6.99d10  ! Solar radius [cm]
END MODULE constants
MODULE DUST
  USE constants
  IMPLICIT NONE
  INTEGER nds,iredo
  integer, parameter :: nbin=10,ndsmax=5,nwv=200,nbinmax=32
  REAL(KIND=DP) tf,amin,amax,as,wv(nwv),da,tstar,adustmean
  REAL(KIND=DP) delh,dmfracism
  REAL(KIND=DP), dimension(ndsmax) :: drho,dfrac,tsub
  REAL(KIND=DP), dimension(ndsmax,nbin) :: tdust  
  REAL(KIND=DP), dimension(ndsmax,nbin) :: adust  
  real(kind=dp), dimension (ndsmax,nbin) :: dmfrac
  REAL(KIND=DP), dimension(ndsmax,nbinmax,nwv) :: q ,qa
CONTAINS
!----------------------------------------------------------------------
   SUBROUTINE DUSTINIT(tstar,tf) ! Initializes all dust variables
      real(kind=dp), intent(in) :: tstar,tf
      call readcomp
      call dustsize
      call qval(tf)
      if(iredo==1) call tdusttables(tstar)
   END SUBROUTINE DUSTINIT
!----------------------------------------------------------------------
  SUBROUTINE readcomp ! Reads in dust composition data
    integer i
   
    drho=0.0
    dfrac=0.0
    tsub=0.0
    open(1,file='dust.inp')
    read(1,*) nds
    if(nds>ndsmax) then
       stop ' Dust species has to be <= 5 '
    endif
    do i=1,nds   
       read(1,*) drho(i),dfrac(i),tsub(i)
    end do
    read(1,*) amin
    read(1,*) amax
    if((amax/amin).lt.1.0001) amax=amin*1.001
    read(1,*) as
    read(1,*) iredo
    close(1) 
    return
    

  END SUBROUTINE readcomp
 
!----------------------------------------------------------------------
  SUBROUTINE dustsize
!    Generate dust size bins in microns 
     integer l,i
     real(kind=dp)  dummy

!     First find dmfrac in smallest bin for ism dust to normalize
!     PAH abundance; amin-0.005 and amax=-.2 both in um
      do l=1,nds
       da=(log10(0.2)-log10(0.005))/(nbin-1)
       adust(l,1)=0.005
       dummy=sqrt(0.005)
       do i=2,nbin
        adust(l,i)=adust(l,i-1)*10.0**da
        dummy=dummy+(adust(l,i))**(4.0-as)
       enddo
      enddo
      dmfracism=0.005**(4-as)/dummy 

      do l=1,nds
       da=(log10(amax)-log10(amin))/(nbin-1)
       adust(l,1)=amin
       dummy=sqrt(amin)
       do i=2,nbin
        adust(l,i)=adust(l,i-1)*10.0**da
        dummy=dummy+(adust(l,i))**(4.0-as)
       enddo
      enddo
      dmfrac=adust**(4-as)/dummy 
      adustmean=sqrt(amin*amax)*1.0d-4
      return
  END SUBROUTINE dustsize

!----------------------------------------------------------------------
      SUBROUTINE qval(taufac)
      real(kind=dp), dimension (ndsmax,nbinmax) :: rad
      integer, dimension (ndsmax,nbinmax) :: nsize
      integer ng,l,i,j,ivis
      real(kind=dp) taufac,frac,diff1,diff2,mfac

!     Q_abs values for dust of given composition from S.Wolf
!     for different ng grain sizes from amin to amax  and 
!     nwv wavelengths from 0.001-2000 um

      open(1,file='dustq.dat')
      read(1,*) l,ng
!      if(l/=nds.or.ng/=nbin) STOP 
      do l=1,nds
        do i=1,ng
         read(1,*) rad(l,i)
         do j=1,nwv
           read(1,*) wv(j),q(l,i,j)
         enddo
      enddo
      enddo
      close(1)

!     Save Qs for the grain size bins used, use nearest radius

comp:      do l=1,nds

gsize:      do i=1,nbin
           nsize(l,i)=0
           do j=2,ng
            if(adust(l,i)>=rad(l,j-1).and.adust(l,i)<=rad(l,j)) then
              diff1=abs(adust(l,i)-rad(l,j-1))/adust(l,i)
              diff2=abs(adust(l,i)-rad(l,j))/adust(l,i)
              if(diff1<diff2) then
               nsize(l,i)=j-1
              else
               nsize(l,i)=j
             endif
           endif
          enddo
          if(nsize(l,i)==0) then
             if(adust(l,i)<rad(l,1)) nsize(l,i)=1
             if(adust(l,i)>rad(l,ng)) nsize(l,i)=ng
          endif
          end do gsize
         end do comp
!   Rewrite relevant Q data in qabs.dat

!     Find index for visual wavelength bang in wv(i)
       ivis=0
       do i=1,nwv
          if(wv(i)>=0.55) then
             ivis=i
             exit
          endif
       enddo
       if(ivis==0) pause ' Did not find Visual wavelength in qabs !'


       mfac=0.0d0
       do j=1,nbin
        mfac=mfac+(adust(1,j)*1.0d-4)**(3-as)*  &
             1.0d-4*(adust(1,j)*(10**da-1.0d0))
       enddo


       taufac=0.0d0
       open(2,file='qabs.dat')
       do l=1,nds
       frac=dfrac(l)/drho(l)
       do i=1,nwv
         if(l==1) wv(i)=1.0e4*wv(i)  ! um to Angstroms
         write(2,'(E9.3,10F8.5)') wv(i),(q(l,nsize(l,j),i),j=1,nbin)
         do j=1,nbin
           qa(l,j,i)=q(l,nsize(l,j),i)
           if(i==ivis) then
           taufac=taufac+qa(l,j,i)*(adust(l,j)*1.0d-4)  &
            **(2.0-as)*frac*(1.0d-4*(adust(1,j)*(10**da-1.0d0)))/mfac
           endif
         enddo
       enddo
       enddo
! Do for only one composition as frac is included in taufac
 
 
       taufac=(3.0/4.0)*taufac ! Multiply this by Dz*n*m_H*eta to get tau

       close(2)
       return
       END SUBROUTINE qval
!----------------------------------------------------------------------
   SUBROUTINE tdusttables(tstar)
      real(kind=dp), intent(in) :: tstar
      integer i,j,k,l,m,nl 
      real(kind=dp) q,w1,w2,f1,f2,w,f,dw,tdd,a,qsum
      real(kind=dp), dimension(nbin) :: sint,flux,tdint
      real(kind=dp), dimension(ndsmax,nbin) :: qavg
  

      open(1,file="dustints.dat")

!     Evaluate the integral over the stellar flux 
       sint=0.0d0
       flux=0.0d0
       qavg=0.0

comp: do l=1,nds

      open(2,file="starspectrum")
      read(2,*) nl
      read(2,*) w1,f1
spec:    do i=2,nl
         read(2,*) w2,f2
         dw=w2-w1
         w=0.5*(w1+w2)
         f=0.5*(f1+f2)
gsize:   do j=1,nbin
           m=0
           if(w.le.wv(1)) m=1
           if(w.ge.wv(nwv)) m=nwv
           if(m.eq.0) then
wloop:          do k=1,nwv
                 if(w1.gt.wv(k).and.w1.le.wv(k+1)) then
                  q=qa(l,j,k)+(w1-wv(k))*(qa(l,j,k)-qa(l,j,k+1))/&
                 (wv(k)-wv(k+1))
                 exit wloop
                 endif
                end do wloop
            else
              q=qa(l,j,m)
            endif
            sint(j)=sint(j)+ f*q* &
                3.1415*(adust(l,j)*1.0e-4)**2.0*dw
            qavg(l,j)=qavg(l,j)+q*f*dw
            flux(j)=flux(j)+f*dw
        end do gsize                
        w1=w2
        f1=f2
       end do spec

       do i=1,nbin
         qavg(l,i)=qavg(l,i)/flux(i)
       enddo

       close(2)  ! This file is re-opened for each dust composition?!

          write(1,'(65E12.5)') tstar,sint
!         Evaluate integral for different dust temperatures
          do i=1,2000
            tdd=i
            do j=1,nbin
               a=adust(l,j)
               call area1(9.0d-2,1.0d3,a,tdd,qsum,j,l)
               tdint(j)=qsum
            enddo
            write(1,'(33E12.5)') tdd,tdint
!            write(*,'(33E12.5)') tdd,tdint
          enddo

       end do comp
       close(1)
        print*,'---> Done with dust integral tabulations'
       return
       stop

   END SUBROUTINE tdusttables 

!----------------------------------------------------------------------
      SUBROUTINE area1(xl,xu,ad,t,fsum,ja,l)
      real(kind=dp) xl,xu,ad,t,fsum,agr,tdd,a,b,y,up
      integer ja,l,jgr,lsp
      common /particle/ agr,tdd,jgr,lsp

      agr=ad
      jgr=ja
      tdd=t
      lsp=l
      a=xl
      b=xu
      call qromb1(a,b,y,up)
      fsum=y
      return
      END SUBROUTINE area1 

!----------------------------------------------------------------------
      function planck(w) 
      implicit double precision (A-H)
      implicit double precision (O-Z) 
      integer, parameter :: nwv=200        
      integer l,k,lsp,jgr,m
      real(kind=dp) w,w1,dqabs,planck,agr
      common /particle/ agr,tdd,jgr,lsp

      planck=(1.1927d11/(w**5.0))*3.1416/(exp(1.44d4/(w*tdd))-1.0)
      if(planck<1.0d-100) then
         planck=0.0d0
         return
      endif
      l=lsp
      w1=w*1.0e4
      m=0
      if(w1.le.wv(1)) m=1
      if(w1.ge.wv(nwv)) m=nwv
      if(m.eq.0) then
        do k=1,nwv
         if(w1.gt.wv(k).and.w1.le.wv(k+1)) then
          dqabs=qa(l,jgr,k)+(w1-wv(k))*(qa(l,jgr,k)-qa(l,jgr,k+1))/&
                (wv(k)-wv(k+1))
           exit
         endif
       enddo
      else
        dqabs=qa(l,jgr,m)
      endif
      planck=planck*dqabs*12.56*(agr*1.0e-4)**2.0
      return
      end function planck

!----------------------------------------------------------------------
      SUBROUTINE qromb1(a,b,ss,up)
      implicit double precision (A-H)
      implicit double precision (O-Z)  
      implicit integer (I-N)

!     Returns as S the integral of the funtion integ from a to b
!     Integration is performed by Romberg's method of order 2k,
!     where, e.g., k=2 is Simpson's rule
      integer, parameter :: jmax=50, jmaxp=jmax+1, k=5, km=k-1
      real(KIND=dp), parameter :: eps=1.0d-4
      real(KIND=dp), dimension(jmaxp) :: s, h

!     Here eps is the fractional accuracy desired, as determined
!     by the extrapolation error estimate; jmax limits the total
!     number of steps; k is the number of points used in the
!     extrapolation.
!     These store the successive trapezoidal approximations and
!     their relative step-sizes.

      up=0.0
      h(1)=1.
      do j=1,jmax
         call trapzd1(a,b,s(j),j)
         if(j.ge.k) then
            call polint1(h(j-km),s(j-km),k,0.0d0,ss,dss)
            if(abs(dss).le.eps*abs(ss))return
         endif
         s(j+1)=s(j)
         h(j+1)=0.25*h(j)
      enddo
      up=1.0
      return
      END SUBROUTINE qromb1


      SUBROUTINE trapzd1(a,b,s,n)
      implicit double precision (A-H)
      implicit double precision (O-Z)  
      implicit integer (I-N)
      integer n,it
      save it

!     This routine computes the nth stage of refinement of an
!     extended trapezoidal rule. f_int is input as the name of
!     the function to be integrated between limits a and b, also
!     input. When called with n=1, the routine returns as s the
!     crudest estimate of the integral. Subsequent calls with n=2,3..
!     will improve the accuracy of s by adding 2^(n-2) additional
!     interior points. s should not be modified between sequential
!     calls.

      if(n.eq.1) then
          s=0.5*(b-a)*(planck(a)+planck(b))
          it=1  
!        it is the number of points to be added on the next call
      else
         tnm=it
         del=(b-a)/tnm
!        This is the spacing of the points to be added
         x=a+0.5*del
         sum=0.
         do j=1,it
           sum=sum+planck(x)
           x=x+del
         enddo
         s=0.5*(s+(b-a)*sum/tnm)
!        This replaces s by its refined value
         it=2*it
      endif
      return
      END SUBROUTINE trapzd1
           

!----------------------------------------------------------------------
      SUBROUTINE polint1(xa,ya,n,x,y,dy)
       implicit double precision (A-H)
       implicit double precision (O-Z)  
      implicit integer (I-N)
        integer n

!     Given arrays xa and ya, each of length n, and given a 
!     value x, this routine returns a value y, and an error
!     estimate dy. If P(x) is the polynomial of degree n-1 such
!     that P(xa_i)=ya_i, i=1,...n, then the returned value y=P(x)

      integer, parameter :: nmax=10
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif) then
           ns=i
           dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo
      y=ya(ns)
      ns=ns-1
      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
                print*, ' 2 input x s are equal in polint-dust1.f'
                stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if(2*ns.lt.n-m) then
           dy=c(ns+1)
        else
           dy=d(ns)
           ns=ns-1
        endif
        y=y+dy
      enddo
      return
      END SUBROUTINE polint1



!-----------------------------------------------------------------
      SUBROUTINE gettdust(tmidr,av,atten,tgas,tdmean)
      real(kind=dp) att,av,atten,tmidr,tdmin,tdmax,tdmean
      real(kind=dp) xa,denk,tgas,denh,denc,f,a_eff
      real(kind=dp) term1,term2
      real(kind=dp), dimension (2001) :: t
      real(kind=dp), dimension (ndsmax,nbin) :: stint
      real(kind=dp), dimension (ndsmax,2000,nbin) :: dstint
      integer l,i,j,ilow,ihi
      integer it0,it1,itn,icrate
      integer :: jcall=0
      save stint,dstint,t

      if(jcall==0) then
!       Read in integrals for dust grains
            open(1,file="dustints.dat")
             do l=1,nds
              i=1
              read(1,*) t(i),(stint(l,j),j=1,nbin)
              do i=2,2001
                 read(1,*) t(i),(dstint(l,i-1,j),j=1,nbin)
              enddo
             enddo
             close(1) 
       end if
       jcall=1  ! After first call do not read dustints.dat again.

!  Initialize array tdust to midplane value
       tdust=tmidr
       ilow=tmidr 
       ilow=ilow-1
!  Set temperature equal to effective of direct (t(i) and diffuse (tmidr)
         att=atten*exp(-av)
         do l=1,nds
          do j=1,nbin
            term1=stint(l,j)*att
            i=tmidr-1
            if(dstint(l,i,j)<term1) then
              if(dstint(l,2000,j)>term1) then
                ihi=2000
                ilow=tmidr-1
tloop:          do 
                 i=(ilow+ihi)/2
                 if(dstint(l,i,j)>term1) then
                  ihi=i
                 else
                  ilow=i
                 endif
                 if((ihi-ilow)<2) exit tloop          
                end do tloop
                tdust(l,j)=(tmidr**4.0+t(i+1)**4.0)**0.25
              else 
                tdust(l,j)=(tmidr**4.0+t(2001)**4.0)**0.25
              endif
            endif
            tdust(l,j)=min(tsub(l),tdust(l,j))
          end do
         end do
              
! Dust grain size distribution constant,k, where n(a)= k a^-s da       
      xa=4.0-as
      denk=(xa/(amax**xa-amin**xa))/(1.0e-4**xa)
      denk=denk/(4.0*pi/3.0) 

! In this formulation, k = denk * eta * rho_gas; where k is the k above

!  Dust collisional heating/cooling calculation
!         denh=0.0
!         denc=0.0
!         do l=1,nds
!          f=denk*dfrac(l)/drho(l)
!           do i=1,nbin
!            f*a_eff=area*number density of grains of size a 
!             a_eff=(1.0d-4*adust(l,i))**(-as+3.0)*(10.0**da-1.0)
!             if(tdust(l,i)>tsub(l)) a_eff=0.0
!              if(tdust(l,i)<tgas) then
!               denc=denc+pi*a_eff*f*(tdust(l,i)-tgas)
!              else
!               denh=denh+pi*a_eff*f*(tdust(l,i)-tgas)         
!              endif
!              delh=denh+denc
!            end do
!         end do

       
! Finding the gas temperature that balances dust heating and cooling
       tdmin=minval(tdust)
       tdmax=maxval(tdust)
       tdmean=(tdmin+tdmax)*0.5
       if((tdmax-tdmin).lt.2.0) return
       do
         denh=0.0
         do i=1,nbin
          do l=1,nds
             f=denk*dfrac(l)/drho(l)
!!            f*a_eff=area*number density of grains of size a 
             a_eff=(1.0d-4*adust(l,i))**(-as+3.0)*(10.0**da-1.0)
             if(tdust(l,i)>tsub(l)) a_eff=0.0
             denh=denh+a_eff*f*(tdust(l,i)-tdmean)         
            end do
         end do
         if(denh>0) then
            tdmin=tdmean
         else
            tdmax=tdmean
         endif
         tdmean=(tdmin+tdmax)*0.5
         if((tdmax-tdmin)<2) exit
       enddo
    return
    END SUBROUTINE gettdust
     
!-----------------------------------------------------------------
      SUBROUTINE fitdust(tmidr,av,atten,tgas,tdmean)
      real(kind=dp) att,av,atten,tmidr,tdmin,tdmax,tdmean
      real(kind=dp) xa,denk,tgas,denh,denc,f,a_eff
      real(kind=dp) term1,term2,xjunk,xp,tp
      real(kind=dp) time1,time2,dummy
      real(kind=dp), dimension (2001) :: t
!     real(kind=dp), dimension (nds,nbin) :: stint
      real(kind=dp), dimension (10) :: stint,uplim,dnlim
      real(kind=dp), dimension (6,10) :: dstpoly
      integer l,i,j,ilow
      integer :: jcall=0
      save stint,dstpoly,t,uplim,dnlim

! Warning, this only works after fits have been made to the
! dust integrals. Only did this for silicates, and with the
! size range 50A-1cm dust size grains and 10 bins.

      l=1
      if(jcall==0) then
!       Read in integrals for dust grains
             open(1,file="dustints.dat")
             if(nds>1) STOP 'ERROR ' 
             read(1,*) t(1),(stint(j),j=1,nbin)
             close(1) 
!        Read in polynomial fits to dust integrals and temperature
            open(1,file='dustpoly')
             do j=1,nbin
               read(1,*) i,dummy,(dstpoly(l,j),l=1,6),dnlim(j),uplim(j)
             enddo
            close(1)
       end if
       jcall=1  ! After first call do not read dustints.dat again.

!  Initialize array tdust to midplane value
       tdust=tmidr
       att=atten*exp(-av)
         do j=1,nbin
           tp=0.0d0
           xp=log10(1.0d-50+stint(j)*att)
           if(xp.lt.dnlim(j)) then
            tp=0.0
           elseif(xp.gt.uplim(j)) then
            tp=2000.0
           else
            tp=dstpoly(1,j)*xp**5+dstpoly(2,j)*xp**4+dstpoly(3,j)*xp**3+         &
               dstpoly(4,j)*xp**2+dstpoly(5,j)*xp+dstpoly(6,j)
            tp=10**tp
           endif
!  Set temperature equal to effective of direct (tp) and diffuse (tmidr)
           tdust(1,j)=(tmidr**4.0+tp**4.0)**0.25
           tdust(1,j)=min(tsub(1),tdust(1,j))
         enddo
 
! Dust grain size distribution constant,k, where n(a)= k a^-s da       
      xa=4.0-as
      denk=(xa/(amax**xa-amin**xa))/(1.0e-4**xa)
      denk=denk/(4.0*pi/3.0) 
! In this formulation, k = denk * eta * rho_gas; where k is the k above

!  Dust collisional heating/cooling calculation
!         denh=0.0
!         denc=0.0
!         do l=1,nds
!          f=denk*dfrac(l)/drho(l)
!           do i=1,nbin
!            f*a_eff=area*number density of grains of size a 
!            a_eff=(1.0d-4*adust(l,i))**(-as+3.0)*(10.0**da-1.0)
!            if(tdust(l,i)<tsub(l)) then
!             if(tdust(l,i)<tgas) then
!              denc=denc+a_eff*f*(tdust(l,i)-tgas)
!             else
!              denh=denh+a_eff*f*(tdust(l,i)-tgas)         
!             endif
!             delh=denh+denc
!            endif
!           end do
!        end do


! Finding the gas temperature that balances dust heating and cooling
! This is being set as the tmidr that was input 
       tdmin=minval(tdust)
       tdmax=maxval(tdust)
       tdmean=(tdmin+tdmax)*0.5
       f=denk*dfrac(1)/drho(1)
       do
          denh=0.0
           do i=1,nbin
!            f*a_eff=area*number density of grains of size a 
             a_eff=(1.0d-4*adust(1,i))**(-as+3.0)*(10.0**da-1.0)
             if(tdust(1,i)>tsub(1)) a_eff=0.0
             denh=denh+a_eff*f*(tdust(1,i)-tdmean)         
           end do
         if(denh>0) then
            tdmin=tdmean
         else
            tdmax=tdmean
         endif
         tdmean=(tdmin+tdmax)*0.5
         if((tdmax-tdmin)<2) exit
       enddo
    return
    END SUBROUTINE fitdust

END MODULE DUST

