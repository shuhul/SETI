!=================================
!    routine to compute alpha
!================================
        
        subroutine alpha_func(nr,r,tprsigma,alpha,temp,time,sm)

        implicit none
        double precision dummy, dummy2, time,sm
        double precision r(nr),tprsigma(nr),alpha(nr),temp(nr)
        double precision alphamax,alphamin
        integer ir,j,k,nr

!  nr       =  size of spatial grid
!  r        =  radius in cm
!  tprsigma =  2 pi r Sigma, if you need sigma, you need to divide
!              this by 2 pi r, so sigma = tprsigma/(2 pi r)
!               [this has been defined as dummy2 below]
!  alpha    =  value of alpha as returned by this function
!  temp     =  Disk gas temperature 
!  time     =  Time in s at the instant this routine is being called 
!  sm       =  stellar mass


        alphamin=1.0d-6
        alphamax=0.1 

        do ir=1,nr

!        constant value of alpha
         alpha(ir)=0.0003

!        Example of power law with radius below
!        alpha = c r^(-a), c=0.0005, a=0.5
!        alpha(ir)=0.0005 * (r(ir)/1.5e13 )**0.5
         

!        dummy2=tprsigma(ir)/(2.0*3.14159*r(ir))       ! This is sigma_gas  

!-------------------------------------------------------------
!        Adding a new gravitational instability  criterion
!        First calculate Toomre Q = c_s omega /( pi * G * sigma) if < 1 unstable
!        Maximum angular mom transport possible => alpha=0.1
!        Add this mode to get the total effective alpha

!         dummy = 7.33e6*sqrt(sm/r(ir))*6.28*temp(ir)/tprsigma(ir)
!         alpha(ir)=alpha(ir)+ 0.1/(1.0+dummy**3.0)

!-------------------------------------------------------------
!        Leave this as is, ensures that your function never results in an alpha
!        that is unrealistic. 
         alpha(ir)=min(alphamax,alpha(ir))
         alpha(ir)=max(alphamin,alpha(ir))
!-------------------------------------------------------------

        enddo
        return
        end
       
!-------------------------------------------------------------
