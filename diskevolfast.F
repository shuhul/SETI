!     =============================================================


!        PROGRAM FOR TIME-DEPENDENT EVOLUTION OF PROTOSTELLAR / 
!               PROTOPLANETARY ACCRETION DISK WITH PHOTOEVAPORATION
!     ==============================================================

#include "configure.h"
#define THEEPS 1.d0
#define INPUT_MAX_LEN 100

!                             MAIN ROUTINE
      program diskevol
      implicit none
      integer ir,iz,iformat,idust,nr,itemp,idum,tsolve_ntemp,ia
      integer ngr
      logical ierror
      doubleprecision tdum,dummy,dum,teuv,rg,n0,rcr,vkep
      doubleprecision tsolve_temp0,tsolve_temp1
      doubleprecision au,cs2,diskmdotinit,pe_mdot_int
      doubleprecision dummy3,dummy4,sclht
      doubleprecision hp(FRSIZE_R)
      doubleprecision sigmadust(FRSIZE_R)
      doubleprecision sigdusta(NDUSTBIN,FRSIZE_R)
      logical fex
      parameter(au=1.496d13)
 
      doubleprecision pi,kk,mp,GG,mu,r,ss
      parameter(pi  = 3.14159265358979324)
      parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
      parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
      parameter(GG  = 6.672d-8)      ! Gravitational constant
      parameter(ss  = 5.6703d-5)     ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
 
#include "common_disk.h"
#include "common_dust.h"
#include "fuvG0.h"
#include "common_radialrt.h"

!     ADDED COMMMENT YESSSSSSS
 
!     -------------------------------------------------------------------------
 
!     Defaults for star, gas, dust 
 
      global_mstar      = 1.99d33    ! Default = Msun
      global_rstar      = 2.96d10    ! Default = Rsun
      global_tstar      = 2.00d3     ! Default = Tsun
      global_gtd        = 1.0d2       ! The default gas-to-dust ratio (not important!)
      global_flang      = 0.001       ! Fixed flaring angle
      global_adust      = 7.1d-4    ! Default grain size
      global_mugas      = 2.3        ! Fixed value for mugas
      global_prandtlinv = 1.d0       ! Default value for the inverse Prandtl number
      global_sigma_min  = 1d-60      ! Lower floor for surface density
      global_sigmad_min  =1d-62      ! Lower floor for dust surface density
      tempbg            = 3.d0      ! Background temperature
!
!     Defaults for time stepping and data dumping
!
      time_dtrel        = 0.01       ! Increase modeling time step by 1% each step
      time_dtmin        = 3.1536d7   ! 1 year is minimum time step
      time_end          = 3.1536d14  ! 10^7 years is finish
      time_tsave        = 3.1536d13  ! Save (by default) each 10^6 years
      time_inctsave     = -1d90      ! Increasing save time switched off by default
      time_maxstep      = 1000000000      ! Maximum number of time steps of simulation
      time_isave        = 11000      ! Max nr of steps before doing a dump
      idump_rt          = 20000     ! Default: no dump of structure
      idump_count       = 100000     ! If dump struct, then ensure always first save full 1+1D sol
      time_savref_mdot  = 0.d0       ! Critical mdot for switch to dt=const
      time_savref_dt    = 1.d5       ! If Mdot<Mdot_crit, then this is the dt
      time_savref_tmin  = 0.d0       ! But only do this if time > tmin
      time_stop_mdot    = 0.d0       ! Stop simulation if mdot drops below this
      time_stop_mdisk   = 1.0d-8     ! Stop simulation if mdisk drops below this
!
!     Default switches of processes in the disk
!
      isw_qvisc         = 0          ! Switch on the viscous heating
      isw_qirrad        = 1          ! Switch on the irradiative heating
      isw_growstar      = 1          ! Switch on the star mass growth
      isw_incmdisk      = 0          ! Switch on the inclusion of mass of disk in mstar
      isw_lumbl         = 0          ! Default: no irradiation by boundary layer
      isw_selfirr       = 0          ! Default: no self-irradiation
      isw_dustevolve    = 0          ! No Dust evolution 
      isw_backreacn     = 0          ! No back reaction of dust on gas
      idust_evap        = 1          ! Include the evaporation of dust in total opac
      dust_tevap        = 1500.      ! Dust evaporation temperature
      global_visheat_mode=0          ! How to compute T_c from T_eff for active accr.
      global_stlumgrow  = 0.d0       ! If >0, then grow Lstar with Mstar
      global_gravinst   = 0          ! Include enhanced viscosity by gravitational instab?
      global_ginst_plaw = 7.d0       ! The default p value for the grav-inst (1-Q)^p 
!
!     Defaults for dust processing stuff
!
      dpr_mode          = 0          ! Dust processing mode
      dpr_tempcryst     = 800.d0     ! Crystallization temperature of dust
      dpr_nu            = 0.d0       ! Rate constant for crystallization
!
!     Defaults for infalling envelope
!
      cloud_mass        = 0.d0       ! Mass of infalling cloud 
      cloud_asound      = 5991.*sqrt(tempbg)  ! Sound speed of infalling cloud (Shu model)
      cloud_omega       = 2.d-14     ! Rotation rate of cloud (Ulrich model)
      cloud_ismooth     = 0          ! =1 ---> Smooth the end of the infall phase
      cloud_smoothparam1= 5d-1       ! Smoothing parameter
      cloud_smoothparam2= 100.*1.4960000d+13 ! Smoothing parameter
      cloud_idistr      = 1          ! Default: Hueso & Guillot way of distrib 
                                     !   infalling matter onto the disk
!
!     Controls for the EUV and FUV photoevaporation
!
      isw_euv_evap      = 101          ! Switch for photoevaporation by EUV photons
      isw_fuv_evap      = 1          ! Switch for photoevaporation by FUV photons
                                     ! =1 use Gorti & Hollenbach gas temperatures
                                     ! =-1 use dust temperature as gas temperature
      fuvevap_nz        = 100         ! The default number of vertical grid points
      fuvevap_zrmin     = 0.02       ! Minimum z/r grid point
      fuvevap_zrmax     = 1.0        ! Maximum z/r grid point
      tsolve_ntemp      = 500        ! Number of pre-stored temperatures
      tsolve_temp0      = 1d-2       ! Lowest pre-stored temperature
      tsolve_temp1      = 12000.     ! Lowest pre-stored temperature
      euvp_phi          = 1.d40        ! EUV luminosity (in number of photons)
      vs_iter_nr        = 20         ! Maximum 40 vertical struct iterations for FUV evap
      vs_iter_crit      = 1.d-1      ! Convergence criterion for vertical struct iter
      pe_mdot           = 1.0d-08    ! PE Mass loss rate in msuns/year
      pe_mdot_int       = 1.0d-08    ! PE Mass loss rate in msuns/year
!
!     Dead zone stuff
!
      isw_deadzone      = 0          ! Default: no dead zone
      deadzone_siglay   = 2d2        ! Total active surface density (two sides of disk)
      deadzone_tempmri  = 2000.d0    ! If temperature is above this, always MRI
      deadzone_alpha    = 1d-6       ! Turbulence in dead zone
!
!     Technical switches (methods)
!
      ml_implicit       = 1          ! Switch for mass loading in implicit way
!
!     Default for infall stuff
!
      ml_mdot           = 0.d0       ! Default: no infall onto star+disk
      ml_mdotcap        = 0.d0       ! Default: no infall onto star directly



      open(unit=1,file='./input/diskevol.inp')
      call read_input_file()
      close(1)
      call parse_input_double('mstar@',6,global_mstar)
      call parse_input_double('rstar@',6,global_rstar)
      call parse_input_double('tstar@',6,global_tstar)
      call parse_input_double('stlumgrow@',10,global_stlumgrow)
      call parse_input_double('sigma_min@',10,global_sigma_min)
      call parse_input_double('a_dust@',7,global_adust)
      call parse_input_double('flang@',6,global_flang)
      call parse_input_double('dtrel@',6,time_dtrel)
      call parse_input_double('dtmin@',6,time_dtmin)
      call parse_input_double('tend@',5,time_end)
      call parse_input_double('tsave@',6,time_tsave)
      call parse_input_double('time_inctsave@',14,time_inctsave)
      call parse_input_integer('maxstep@',8,time_maxstep)
      call parse_input_integer('isave@',6,time_isave)
      call parse_input_double('tempbg@',7,tempbg)
      call parse_input_integer('incl_gravinst@',14,global_gravinst)
      call parse_input_double('ginst_plaw@',11,global_ginst_plaw)
      call parse_input_integer('incl_dustproc@',14,dpr_mode)
      call parse_input_double('dpr_tempcryst@',14,dpr_tempcryst)
      call parse_input_double('dpr_nu@',7,dpr_nu)
      call parse_input_integer('incl_dustevap@',14,idust_evap)
      call parse_input_double('dust_tevap@',11,dust_tevap)
      call parse_input_integer('visheat_mode@',13,global_visheat_mode)
      call parse_input_double('cloud_mass@',11,cloud_mass)
      call parse_input_double('cloud_asound@',13,cloud_asound)
      call parse_input_double('cloud_omega@',12,cloud_omega)
      call parse_input_integer('cloud_idistr@',13,cloud_idistr)
      call parse_input_integer('cloud_ismooth@',14,cloud_ismooth)
      call parse_input_double('cloud_smoothparam1@',19,
     %    cloud_smoothparam1)
      call parse_input_double('cloud_smoothparam2@',19,
     %    cloud_smoothparam2)
      call parse_input_integer('incl_lumbl@',11,isw_lumbl)
      call parse_input_integer('incl_selfirr@',13,isw_selfirr)
      call parse_input_integer('incl_qvisc@',11,isw_qvisc)
      call parse_input_integer('incl_qirrad@',12,isw_qirrad)
      call parse_input_double('prandtlinv@',11,global_prandtlinv)
      call parse_input_double('schmidtinv@',11,global_prandtlinv)
      call parse_input_integer('incl_growstar@',14,isw_growstar)
      call parse_input_integer('incl_incmdisk@',14,isw_incmdisk)
      call parse_input_integer('ml_implicit@',12,ml_implicit)
      call parse_input_integer('incl_euv_evap@',14,isw_euv_evap)
      call parse_input_double('euv_phi@',8,euvp_phi)
      call parse_input_integer('incl_fuv_evap@',14,isw_fuv_evap)
      call parse_input_integer('fuvevap_nz@',11,fuvevap_nz)
      call parse_input_double('fuvevap_zrmin@',14,fuvevap_zrmin)
      call parse_input_double('fuvevap_zrmax@',14,fuvevap_zrmax)
      call parse_input_integer('tsolve_ntemp@',13,tsolve_ntemp)
      call parse_input_double('tsolve_temp0@',13,tsolve_temp0)
      call parse_input_double('tsolve_temp1@',13,tsolve_temp1)
      call parse_input_integer('idump_rt@',9,idump_rt)
      call parse_input_integer('vs_iter_nr@',11,vs_iter_nr)
      call parse_input_double('vs_iter_crit@',13,vs_iter_crit)
      call parse_input_integer('incl_deadzone@',14,isw_deadzone)
      call parse_input_integer('incl_dustevolve@',16,isw_dustevolve)
      call parse_input_integer('incl_backreacn@',15,isw_backreacn)
      call parse_input_double('pe_mdot@',8,pe_mdot)
      call parse_input_double('deadzone_siglay@',16,deadzone_siglay)
      call parse_input_double('deadzone_tempmri@',17,deadzone_tempmri)
      call parse_input_double('deadzone_alpha@',15,deadzone_alpha)
      call parse_input_double('savref_mdot@',12,time_savref_mdot)
      call parse_input_double('savref_dt@',10,time_savref_dt)
      call parse_input_double('savref_tmin@',12,time_savref_tmin)
      call parse_input_double('stop_mdot@',10,time_stop_mdot)
      call parse_input_double('stop_mdisk@',11,time_stop_mdisk)
      call check_all_lines_ok(ierror)
      if(ierror) stop

!     Messages

      if(global_gravinst.ne.0) write(*,*) '--> Including ',
     %             'gravitational instability alpha-enhancement'
      if(dpr_mode.ne.0) write(*,*) '--> Including dead zone treatment'
      if(idust_evap.ne.0) write(*,*) '--> Including evaporation of ',
     %             'dust in computing midplane temperature'
      if(isw_lumbl.ne.0) write(*,*) '--> Including irradiation by ',
     %             'the stellar accretion column'
      if(isw_selfirr.ne.0) write(*,*) '--> Including self-irradiation ',
     %             'of outer disk by active accretion'
      if(isw_qvisc.ne.0) write(*,*) '--> Including heating by viscous ',
     %             'dissipation of accretion energy'
      if(isw_qirrad.ne.0) write(*,*) '--> Including heating by ',
     %             'irradiation by central star'
      if(isw_growstar.ne.0) write(*,*) '--> Including stellar mass ',
     %             'increase due to accretion'
      if(isw_incmdisk.ne.0) write(*,*) '--> Including Mdisk into ',
     %             'gravity'
      if(ml_implicit.eq.0) write(*,*) 'WARNING: dSigma/dt is not done ',
     %             'with implicit integration!!!!'
      if(isw_euv_evap.ne.0) write(*,*) '--> Including EUV ',
     %             'photoevaporation, using mode ',isw_euv_evap
      if(isw_fuv_evap.ne.0) write(*,*) '--> Including FUV ',
     %             'photoevaporation, using mode ',isw_fuv_evap
      if(isw_deadzone.ne.0) write(*,*) '--> Including dead zone ',
     %             'treatment'
!
!     Warnings
!
      if(global_visheat_mode.ne.0) then
          write(*,*) '--------------------------------------'
          write(*,*) 'WARNING: Non-standard viscous heating!'
          write(*,*) '--------------------------------------'
      endif
      if((isw_euv_evap.ne.0).and.(euvp_phi.eq.0.d0)) then
          write(*,*) 'EUV Photoevap is switched on, but Phi=0.d0...'
      endif
      write(*,*) '--> Photoevaporation rate ', pe_mdot
!
!     Read the radial grid
!
      open(unit=1,file='./input/rgrid.inp')
      read(1,*) grid_nr
      if(grid_nr.gt.FRSIZE_R) then
          write(*,*) 'Grid size too big'
          stop 13
      endif
      do ir=1,grid_nr
          read(1,*) grid_r(ir)
      enddo
      close(1)
!
! calculate stellar luminosity for gast.F
          starl=(4.0*3.14159*global_rstar**2.0*5.67d-5*
     %          global_tstar**4.0)


!     Read the surface density sigma, both the total one and the
!     dust species sigma
!
      ngr=NDUSTBIN
      open(unit=1,file='./input/sigma.inp')
      read(1,*) nr,global_ndust
      if(nr.ne.grid_nr) then
          write(*,*) 'Grid size sigma.inp not consistent with radius'
          stop 13
      endif
      do ir=1,grid_nr
         gastodust(ir)=global_gtd
         planetesimal(ir)=0.0d0
         read(1,*) dummy
         disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*dummy
         if(disk_tprsigma(ir).lt.2.d0*pi*grid_r(ir)*
     %        global_sigma_min) then
             disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*
     %            global_sigma_min
         endif
c        All water ice below with X(O)=3.0e-4
        do ia=1,ngr
         disk_tprsigmaice(ir,ia,1)=1.8e-4*18.0*disk_tprsigma(ir)/ngr
         disk_tprsigmaice(ir,ia,2)=1.4e-4*28.0*disk_tprsigma(ir)/ngr
         disk_tprsigmadust(ia,ir) = disk_tprsigma(ir)/gastodust(ir)/ngr
        enddo
      enddo
      close(1)
      do ir=1,grid_nr
      disk_tprsigmadtot(ir)=sum(disk_tprsigmadust(1:ngr,ir)) +
     %           sum(disk_tprsigmaice(ir,1:ngr,1))
     %         +sum(disk_tprsigmaice(ir,1:ngr,2))

      enddo

!
!     Read the viscous alpha
!
      open(unit=1,file='./input/alpha.inp')
      read(1,*) nr
      if(grid_nr.ne.nr) then
          write(*,*) 'Grid size of alpha.inp not equal'
          stop 13
      endif
      do ir=1,grid_nr
          read(1,*) global_alpha(ir)
      enddo
      close(1)
!
!     Read the Rosseland mean opacity table for dust 
!
      inquire(file='./input/rossmean_dust.inp',exist=fex)
      if(fex) then
!
!         New way
!
          inquire(file='./input/rossmean.inp',exist=fex)
          if(fex) then
              write(*,*) 'ERROR: Either rossmean.inp or ',
     %                   'rossmean_dust.inp'
              stop
          endif
          open(unit=1,file='./input/rossmean_dust.inp')
          read(1,*) ross_ntemp
          if(ross_ntemp.gt.FRSIZE_TEMP) then
              write(*,*) 'ERROR: Nr of temperatures in the rossmean.inp'
              write(*,*) '       greater than FRSIZE_TEMP...'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) ross_temp(itemp),ross_kappa(itemp,1)
          enddo
          close(1)
      else
!
!         Old way
!
          inquire(file='./input/rossmean_dust.inp',exist=fex)
          if(fex) then
              write(*,*) 'ERROR: Either rossmean.inp or ',
     %                   'rossmean_dust.inp'
              stop
          endif
          write(*,*) 'BACKWARD COMPATIBILITY MODE rossmean.inp...'
          open(unit=1,file='./input/rossmean.inp')
          read(1,*) ross_ntemp
          if(ross_ntemp.gt.FRSIZE_TEMP) then
              write(*,*) 'ERROR: Nr of temperatures in the rossmean.inp'
              write(*,*) '       greater than FRSIZE_TEMP...'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) ross_temp(itemp),ross_kappa(itemp,1)
          enddo
          close(1)
      endif
!
!     Check if rossmean_gas.inp exists. If not, and if dust evaporation
!     is included, then protest...
!
      if(idust_evap.ne.0) then 
          inquire(file='./input/rossmean_gas.inp',exist=fex)
          if(.not.fex) then 
              write(*,*) 'ERROR: Dust evaporation cannot be included'
              write(*,*) ' if you dont have a file rossmean_gas.inp.'
              stop
          endif
          open(unit=1,file='./input/rossmean_gas.inp')
          read(1,*) idum
          if(idum.ne.ross_ntemp) then
              write(*,*) 'ERROR: Nr of temperatures in the '
              write(*,*) '       rossmean_gas.inp not equal '
              write(*,*) '       to rosssmean.inp'
              stop
          endif
          do itemp=1,ross_ntemp
              read(1,*) dum,ross_kappa(itemp,2)
              if(abs(dum/ross_temp(itemp)-1.d0).gt.1d-6) then
                  write(*,*) 'ERROR: Temperature grid of '
                  write(*,*) '       rossmean_gas.inp not equal '
                  write(*,*) '       to that of rossmean.inp'
                  stop
              endif
          enddo
          close(1)
      endif

!     Read sigmadot profile
      open(unit=1,file='./input/sigmadot.inp')
      read(1,*) nr
      if(grid_nr.ne.nr) then
          write(*,*) 'Grid size of sigmadot.inp not equal'
          stop 
      endif
      do ir=1,grid_nr
          read(1,*) photoev_isigdot(ir)
          photoev_isigdot(ir)=photoev_isigdot(ir)*(pe_mdot/pe_mdot_int)
      enddo
      close(1)
!
!     Set grain size to that specified in input
      do ir=1,nr
       grsize(ir)= global_adust ! *sqrt(dexp(time_time/1.5d13) )
       kappadec(ir)=1.0
       gastodust(ir)=global_gtd
       pahscale(ir)=(100.0/global_gtd)*sqrt(3.16d-6/global_adust) ! = 1 for ISM
      enddo
      dust_agrain=global_adust


!     Initialuze mean particle mass to 2.3 will be computed later
       do ir=1,nr
         do iz=1,fuvevap_nz
           mup(iz,ir)=2.3d0
         enddo
       enddo

      if(isw_fuv_evap.ne.0) then
          call radialrt_init(tsolve_ntemp,tsolve_temp0,tsolve_temp1)
      endif
!
!     For the sake of the dumping of 1+1D radiative transfer solutions,
!     compute a first estimate of the temperature of the disk
!
      if(isw_fuv_evap.ne.0) then
          do ir=1,nr
              if(ir.gt.1) then
!
!                 First make a simple estimate of the midplane
!                 temperature without accounting for accretion
!                  
                  disk_temp(ir) = global_tstar*(0.5*global_flang*
     %                 (global_rstar/grid_r(ir))**2)**0.25
                  disk_temp(ir)=min(dust_tevap,disk_temp(ir))
                  hp(ir)  = sqrt(kk*disk_temp(ir)*grid_r(ir)**3/
     %                      (global_mugas*mp*GG*global_mstar))
              endif
              gastodust(ir) = global_gtd
              sigmadust(ir) = disk_tprsigma(ir)/(2*pi*grid_r(ir))/
     %                        gastodust(ir)
          enddo
          hp(1) = hp(2)
          disk_temp(1) = disk_temp(2)
      endif

!
!     Determine the fractions of water and CO ice from the partial
!     pressure criterion for both
!
      do ir=1,grid_nr
         vkep=sqrt(GG*global_mstar/grid_r(ir))
         sclht=(5994.0*sqrt(disk_temp(ir)))/(vkep/grid_r(ir))
         call geticefrac(disk_tprsigma(ir),sclht,disk_temp(ir),
     &     grid_r(ir),dummy3,dummy4)
         ! water
         dummy=disk_tprsigma(ir)*18.0*1.8d-4
         disk_tprsigmaice(ir,1:ngr,1) = dummy3*dummy/ngr ! ice fraction
         disk_tprsigmavapor(ir,1) = (1.0-dummy3)*dummy ! vapor fraction
         if(dummy3<0) snowrad = grid_r(ir) 
         ! co
         dummy=disk_tprsigma(ir)*28.0*1.4d-4
         disk_tprsigmaice(ir,1:ngr,2) = dummy4*dummy/ngr ! ice fraction
         disk_tprsigmavapor(ir,2) = (1.0-dummy4)*dummy ! vapor fraction
         disk_tprsigmadtot(ir) = sum(disk_tprsigmadust(1:ngr,ir)) + 
     %         sum(disk_tprsigmaice(ir,1:ngr,1)) +
     %         sum(disk_tprsigmaice(ir,1:ngr,2))
         gastodust(ir) = disk_tprsigma(ir)/disk_tprsigmadtot(ir)
      enddo



      if(time_inctsave.gt.0.d0) then
          varsdt = .true.
      else
          varsdt = .false.
      endif


      open(unit=1,file=outdir//'/sigma.dat',status='unknown')
      write(1,*) grid_nr,global_ndust
      close(1)
      open(unit=1,file=outdir//'/sigmadust.dat',status='unknown')
      write(1,*) grid_nr, NDUSTBIN
      close(1)
      open(unit=1,file=outdir//'/velo.dat',status='unknown')
      write(1,*) grid_nr,global_ndust
      close(1)
      open(unit=1,file=outdir//'/temperature.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/visc.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/alpha.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/mdot.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/mstar.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file=outdir//'/sigmadotevap.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/gastodust.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/ices.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/planetesimal.dat',status='unknown')
      write(1,*) grid_nr
      close(1)
      open(unit=1,file=outdir//'/j_outflow.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file=outdir//'/infall.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file=outdir//'/infall2.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file=outdir//'/time.dat',status='unknown')
      write(1,*)
      close(1)
      open(unit=1,file=outdir//'/rtdump.info',status='unknown')
      write(1,*) 0
      close(1)
      open(unit=1,file=outdir//'/rtdump_times.dat',status='unknown')
      write(1,*) 0
      close(1)
c
c     Do the run now
c
      call do_run()
c
      end

c     --------------------------------------------------------------
c                             DO THE RUN
c     --------------------------------------------------------------
      subroutine do_run
c
      use dust
c
      implicit none
      integer it,itsave,ir,idust,irs,ia,whichpf,nsi,iceformn
      doubleprecision timesave,year,dummy,xp,diskmass,tsig,epspf
      doubleprecision ncritco,ncrith2o,sclht,rh2o,dummy2
      doubleprecision dummy3, dummy4, zdiffusion,schmidtno
      doubleprecision diff_frac
      doubleprecision ckpr,ckqr,vkep,zactual,zcrit,sigdusttot
      doubleprecision sigdusta(NDUSTBIN,FRSIZE_R),sigda(FRSIZE_R)
      doubleprecision stokes(NDUSTBIN,FRSIZE_R),stkin(FRSIZE_R)
      doubleprecision photoev_ape(FRSIZE_R)
      doubleprecision fracd(FRSIZE_R) 
      doubleprecision fracw(FRSIZE_R) 
      doubleprecision fracco(FRSIZE_R) 
c      doubleprecision src(FRSIZE_R),crystrate(FRSIZE_R)
      logical simstop,dump,savref
      parameter(year= 3.1536d7)      ! Year                    [s]

#include "common_disk.h"
#include "common_radialrt.h"

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ! A. initialize, get sigmadot to find dt
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c
c     Reset time stuff
c
      itsave   = 1
      timesave = 0.d0
c         
c     Photoevaporation [ initial call to determine size of timestep]
c
      if((isw_euv_evap.ne.0).or.(isw_fuv_evap.ne.0)) then
          call photoevap(grid_nr,global_mstar,grid_r,disk_tprsigma,
     %         disk_temp,euvp_phi,
     %         isw_euv_evap,isw_fuv_evap,
     %         euvp_nbase,photoev_sigdot,photoev_izev)
      endif

c     Set water snowline location to 1au initially
       rh2o=snowrad
c
c     First save data
c
      call save_data()
c
c     Then do the time loop
c
      simstop = .false.
      dump    = .false.
      savref  = .false.

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      ! B. Start time loop                    
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      do it=1,time_maxstep

c         Reset flags
c
          global_flag_lowtemp = 0
c
c         Set the time step & Save current value of gas surface density 
c
          tsig=1.0d40
          do ir=1,grid_nr
            disk_tprsigma_prev(ir)=disk_tprsigma(ir)
            if(photoev_sigdot(ir).gt.1.0e-25.and.
     %        disk_tprsigma(ir).gt.1.0e-6) then
            tsig=min(tsig,disk_tprsigma(ir)/photoev_sigdot(ir))
            endif
          enddo
          tsig=0.5*tsig ! factor of 2 resolution at least
          time_dt = max(time_dtmin,time_time*time_dtrel)
          if(tsig.gt.1000*year.and.tsig.lt.time_dt) time_dt=tsig
c
c         Message
c
         if(mod(it,100).eq.0) then
          write(*,*) 'Time step ',it,', Time = ',time_time/year,' year'
         endif
c
c         Prepare mass loading sources
c
          if(cloud_mass.gt.0.d0) then
c
c             Check if cloud parameters are consistent 
c
              if((cloud_asound.eq.0.d0).or.(cloud_omega.eq.0.d0)) then
                  write(*,*) 'ERROR: Input parameters of cloud contain'
                  write(*,*) '       values of 0.'
                  stop
              endif
c
c             Compute the rotating Shu/Ulrich infall model, and compute the
c             mass accretion rate onto the surface of the disk (psi(R)) and
c             the mass accretion rate inward of the capture radius (the 
c             `effective radius' of the central star).
c
              call shuulrich_infall(grid_nr,grid_r,cloud_mass,
     %                      cloud_asound,cloud_omega,time_time,
     %                      ml_psi,ml_mdot,ml_mdotcap,ml_rcentr,
     %                      cloud_ismooth,cloud_smoothparam1,
     %                      cloud_smoothparam2,cloud_idistr)
          endif
c ------------------------------------------        
c     Photoevaporation
c-------------------------------------------
          if((isw_euv_evap.ne.0).or.(isw_fuv_evap.ne.0)) then

           call photoevap(grid_nr,global_mstar,grid_r,disk_tprsigma,
     %            disk_temp,euvp_phi,
     %            isw_euv_evap,isw_fuv_evap,
     %            euvp_nbase,photoev_sigdot,photoev_izev)

c   Gradually ramp up photoevaporation over 1e5 years
           if(time_time.lt.3.15d12) then
             do ir=1,grid_nr
              dummy=min(1.0,time_time/3.15d12)
              photoev_sigdot(ir)=photoev_sigdot(ir)*dummy
             enddo
           endif

          endif
c         Make new prescription for alpha

c         call alpha_func(grid_nr,grid_r,disk_tprsigma,global_alpha,
c     %         disk_temp,time_time,global_mstar)

c
c         First (re)compute the gravitational mass
c
          if(isw_incmdisk.ne.0) then
c
c             Include the mass of the disk inward of R to the total 
c             mass inducing the gravity.
c
              disk_mgrav(1)  = global_mstar
              do ir=2,grid_nr
                  disk_mgrav(ir) = disk_mgrav(ir-1) + 0.5 * 
     %                 ( disk_tprsigma(ir) + disk_tprsigma(ir-1) ) *
     %                 ( grid_r(ir) - grid_r(ir-1) )
              enddo
          else
c
c             Only use the star as the source of gravity
c     
              do ir=1,grid_nr
                  disk_mgrav(ir) = global_mstar
              enddo
          endif
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          ! C. Get the Sigmadust in each grain and set its ices
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c
c   some switches
          iceformn=1 ! 0 for no ices
        ! whichpf=0  Skip planetesimal formation completely
        ! whichpf=1  Streaming instability according to Carrera+2015
        ! whichpf=2  Estrada+2016 Eq 62 which says gastodust<sqrt(stmax/alpha)
        ! whichpf=3  Form pebbles instead
          whichpf=1 ! planetformation mode
          schmidtno = 1.0 ! schmidt no for vertical diffusion

c          
c    distribute the dust grains into size bins from collfrag
c
c     Fraction of dust is equal to gas surface density by gas to dust
c     times the fractional mass in this size bin


       do ir=1,grid_nr
         if(isw_dustevolve.gt.1) then ! =0 means one grain size
            dummy=disk_tprsigma(ir)/(2.0*pi*grid_r(ir))
            vkep=sqrt(GG*global_mstar/grid_r(ir))
            call dustcollision(global_alpha(ir),dummy,
     +         disk_tprsigmadtot(ir)/(2*pi*grid_r(ir)),
     +          disk_temp(ir),grsize(ir),grid_r(ir),rh2o)
         endif
!        from the gas sigma and gastodust, get the dust sigmas at each a
         dummy = sum(disk_tprsigmadust(1:nbin,ir))
         do ia=1,nbin
            sigdusta(ia,ir)= dummy*dmfrac(1,ia) 
            sigdusta(ia,ir)=max(sigdusta(ia,ir),global_sigmad_min/10.0)
         enddo

!        find the total ice that exists and distribute it with the
!        same mass fraction determined by collisions,
         dummy = sum(disk_tprsigmaice(ir,1:nbin,1))
         dummy2 = sum(disk_tprsigmaice(ir,1:nbin,2))

         ! normalize these by the dust since they evolved
         do ia=1,nbin
          disk_tprsigmaice(ir,ia,1)=dummy*dmfrac(1,ia)
          disk_tprsigmaice(ir,ia,2)=dummy2*dmfrac(1,ia)
         enddo

c           Here find the ratio of sigmadust for the smallest
c           size bin and normalize it for an equivalent ism
c           value to get a pah abundance scaling
c           (disk_tprsigma cancels out in scaling)
         pahscale(ir)=(100.0/gastodust(ir))*(dmfrac(1,1)/dmfracism)
       enddo



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! D. Advance gas, vapor  by one timestep
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

c         First compute radial velocity
c      
c         Now do the diffusion step for the gas+dust
          call compute_vr_v0_difcoef()
          call update_sigma_diff()

          if(iceformn.eq.1) then
c
c             Do this for water vapor, disk_vr is used again for dust
c
              call vapcompute_vr_v0_difcoef(1) ! water
              call update_sigma_vapor_diff(1)  ! water
              call vapcompute_vr_v0_difcoef(2) ! CO
              call update_sigma_vapor_diff(2)  ! CO
          endif

c         Check if we are about to cross the save time
c
          if(timesave+time_dt.ge.time_tsave) then
            time_dt = time_tsave - timesave
            dump = .true.
          endif
c
c         Check if we are about to cross the time-limit
c
          if(time_time+time_dt.gt.time_end) then
             time_dt = time_end - time_time
             simstop = .true.
          endif

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          ! E. Now advance dust by one timestep
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

      if(isw_dustevolve.gt.0) then

c    Dust radial evolution is added here in similar fashion to the gas,
c    with a call first to get sigmadot, and then compute radial velocity
c    and then do the diffusion step. 

c     now we need to see how much dust diffuses/advects
C     Need Stokes parameter
       do ir=1,grid_nr
         if(isw_backreacn > 0) decel(ir)=0.0d0
         do ia=1,nbin
           if(isw_dustevolve.gt.2) then
           stokes(ia,ir)=pi*1.0d-4*adust(1,ia)*drho(1)*0.5*
     %      (2.0*pi*grid_r(ir)/disk_tprsigma(ir))
           else
             stokes(ia,ir)=1.0d-15
           endif
         enddo
      enddo

c     Call for each dust grain size:
         do ia=1,nbin
          do ir=1,grid_nr
            stkin(ir)=stokes(ia,ir)
c           Check to see if this grain size photoevaporates or not
c           but first set default pe mass loss to zero. 
            photoev_ape(ir)=0.0
c             if(adust(1,ia)<fdustpe(ir)) photoev_ape(ir)=
c     %      photoev_sigdot(ir)*(sigdusta(ia,ir)/disk_tprsigma(ir)) 
          enddo

          call dust_compute_vr_v0_difcoef(stkin,ckpr,ckqr)

          sigda(1:grid_nr)=sigdusta(ia,1:grid_nr)+1.0d-60
          call dust_update_sigma_diff(sigda,ckpr,ckqr,photoev_ape)
          disk_tprsigmadust(ia,1:grid_nr)=sigda(1:grid_nr)
          
          sigda(1:grid_nr)=disk_tprsigmaice(1:grid_nr,ia,1)+1.0d-60
          call dust_update_sigma_diff(sigda,ckpr,ckqr,photoev_ape)
          disk_tprsigmaice(1:grid_nr,ia,1)=sigda(1:grid_nr)

          sigda(1:grid_nr)=disk_tprsigmaice(1:grid_nr,ia,2)+1.0d-60
          call dust_update_sigma_diff(sigda,ckpr,ckqr,photoev_ape)
          disk_tprsigmaice(1:grid_nr,ia,2)=sigda(1:grid_nr)
         enddo

c 10au=410,30au=480,100au=570,300au=643

         do ir=1,grid_nr
           disk_tprsigmadtot(ir)=sum(disk_tprsigmadust(1:nbin,ir))+
     %     sum(disk_tprsigmaice(ir,1:nbin,1))+
     %     sum(disk_tprsigmaice(ir,1:nbin,2))
           gastodust(ir) = disk_tprsigma(ir)/disk_tprsigmadtot(ir)
         enddo


!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         ! F. Find how much vapor diffused in z, add it to ices
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
         do ir=1,grid_nr
         vkep=sqrt(GG*global_mstar/grid_r(ir))
         sclht=(5994.0*sqrt(disk_temp(ir)))/(vkep/grid_r(ir))!scaleheight = cs/omega
          
         ! first release all the vapor from inside snowlines
         call geticefrac(disk_tprsigma(ir),sclht,disk_temp(ir),
     &     grid_r(ir),dummy3,dummy4)

         if(dummy3<0) then
             ! inside radial water snowline
             do ia=1,nbin
              disk_tprsigmavapor(ir,1)=disk_tprsigmavapor(ir,1)+
     &         disk_tprsigmaice(ir,ia,1)
              disk_tprsigmaice(ir,ia,1) = 0.0
             enddo
             rh2o = grid_r(ir) 
         endif
         if(dummy4<0) then
             ! inside radial co snowline
             do ia=1,nbin
              disk_tprsigmavapor(ir,2)=disk_tprsigmavapor(ir,2)+
     &         disk_tprsigmaice(ir,ia,2)
              disk_tprsigmaice(ir,ia,2) = 0.0
             enddo
         endif

         ! mixing length in this timestep, and fractional Sigma within
         zdiffusion = sqrt(4.0*time_dt*disk_nu(ir)/schmidtno)
         diff_frac = erf(0.7071*(1.0+zdiffusion/sclht))-erf(0.7071)
         diff_frac=0.0


         ! assume that all grains <10um un size (an estimate based on
         ! the scaleheight (hd ~ schlt/(1+st/alpha)^0.5 (estrada+)
         ! the st is small enough that these grains diffuse with gas
         ! and have the same diffusion fraction

         ! H2O
         ! total water ice and vapor at this radius and time 
         dummy=disk_tprsigmavapor(ir,1)
         dummy2=sum(disk_tprsigmaice(ir,1:nbin,1))
         if(dummy2 >0) then ! if ice then consider diffusion
              do ia=1,nbin
!               disk_tprsigmadust(ia,ir)=disk_tprsigmadust(ia,ir)+
!    %                   diff_frac*disk_tprsigmavapor(ir,1)/nbin 
                disk_tprsigmaice(ir,ia,1)=disk_tprsigmaice(ir,ia,1)+
     %                   diff_frac*disk_tprsigmavapor(ir,1)/nbin 
              enddo
              disk_tprsigmavapor(ir,1)=(1.0d0-diff_frac)*
     %                        disk_tprsigmavapor(ir,1)
              ! now release vapor from diffused grains
              do ia=1,nbin
              if(adust(1,ia)<10.0d-4) then
                ! below is the fraction (sigma) that diffuses to surface
                dummy = diff_frac*disk_tprsigmaice(ir,ia,1) 
                disk_tprsigmavapor(ir,1)= disk_tprsigmavapor(ir,1)+
     %             dummy
!               disk_tprsigmadust(ia,ir)=disk_tprsigmadust(ia,ir)-
!    %             dummy
                disk_tprsigmaice(ir,ia,1)=disk_tprsigmaice(ir,ia,1)-
     %             dummy
              endif
              enddo
          endif

         ! CO 
         dummy=disk_tprsigmavapor(ir,2)
         dummy2=sum(disk_tprsigmaice(ir,1:nbin,2))
         if(dummy2 <0) then ! if ice then consider diffusion
              do ia=1,nbin
!               disk_tprsigmadust(ia,ir)=disk_tprsigmadust(ia,ir)+
!    %                   diff_frac*disk_tprsigmavapor(ir,2)/nbin 
                disk_tprsigmaice(ir,ia,2)=disk_tprsigmaice(ir,ia,2)+
     %                   diff_frac*disk_tprsigmavapor(ir,2)/nbin 
              enddo
              disk_tprsigmavapor(ir,2)=(1.0d0-diff_frac)*
     %                        disk_tprsigmavapor(ir,2)
              ! now release vapor from diffused grains
              do ia=1,nbin
              if(adust(1,ia)<10.0d-4) then
                ! below is the fraction (sigma) that diffuses to surface
                dummy = diff_frac*disk_tprsigmaice(ir,ia,2) 
                disk_tprsigmavapor(ir,2)= disk_tprsigmavapor(ir,2)+
     %             dummy
!               disk_tprsigmadust(ia,ir)=disk_tprsigmadust(ia,ir)-
!    %             dummy
                disk_tprsigmaice(ir,ia,2)=disk_tprsigmaice(ir,ia,2)-
     %             dummy
              endif
              enddo
          endif
          disk_tprsigmadtot(ir)=sum(disk_tprsigmadust(1:nbin,ir))+
     %     sum(disk_tprsigmaice(ir,1:nbin,1)) +
     %     sum(disk_tprsigmaice(ir,1:nbin,2))

          gastodust(ir)=disk_tprsigma(ir)/disk_tprsigmadtot(ir)
       enddo ! close ir loop



!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
       ! Check if planetesimals have formed
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

c      Streaming instability
c    PLANET FORMATION BEGIN
      epspf=0.90
      if(whichpf.gt.0) then
       do ir=1,grid_nr
        zactual=log10(1.0d0/(1.0+gastodust(ir))) ! Total Z
c       mass weighted stokes, needs this sigdusttot below, reset later
        sigdusttot=0.0 ! Total Sigma_d
        dummy=0.0
        do ia=1,nbin
         dummy2=disk_tprsigmadust(ia,ir) + disk_tprsigmaice(ir,ia,1)+
     %        disk_tprsigmaice(ir,ia,2)
         dummy=dummy+sqrt(stokes(ia,ir)+alpha_turb)
     %       *dummy2/(disk_tprsigma(ir)+dummy2)
        enddo
        dummy=dummy/sqrt(alpha_turb)
c       check that midplane density is > 1
        if(dummy.ge.1.0) then
        sigdusttot=0.0
        do ia=nbin,1,-1
         dummy2=disk_tprsigmadust(ia,ir) + disk_tprsigmaice(ir,ia,1)+
     %        disk_tprsigmaice(ir,ia,2)
         sigdusttot=sigdusttot+dummy2
         zactual=log10(sigdusttot/(sigdusttot+disk_tprsigma(ir)))
         zcrit=-1.86d0+0.3*(log10(stokes(ia,ir))+0.98)**2.0
         if(whichpf.eq.2) zcrit= - 0.3 ! simple gastodust=1 criterion (log 1/2)
         if(zactual>zcrit) then
!         check if stokes bounds are met
          if(whichpf.eq.2.or.stokes(ia,ir).gt.3.0d-3) then
c          First put all this in planetesimals and then remove from sigmadust
            planetesimal(ir)=planetesimal(ir)+ dummy2*epspf
            disk_tprsigmaice(ir,ia,1)=disk_tprsigmaice(ir,ia,1)*
     %           (1.0d0-epspf)
            disk_tprsigmaice(ir,ia,2)=disk_tprsigmaice(ir,ia,2)*
     %           (1.0d0-epspf)
            disk_tprsigmadust(ia,ir)=
     %                 disk_tprsigmadust(ia,ir)*(1.0d0-epspf)
          endif !  stokes block
         endif !  zcrit block
        enddo ! ia loop
        endif ! end dummy>1 block
       enddo ! ir loop
c      Streaming instability end               

c      Pebble formation




c      Now update gastodust ratio after planetesimals have formed
       do ir=1,grid_nr
          disk_tprsigmadtot(ir)=sum(disk_tprsigmadust(1:nbin,ir))+
     %     sum(disk_tprsigmaice(ir,1:nbin,1)) +
     %     sum(disk_tprsigmaice(ir,1:nbin,2)) 
          gastodust(ir)=disk_tprsigma(ir)/disk_tprsigmadtot(ir)
       enddo

       endif ! whichpf if block

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         ! G. save if needed and loop back
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

         endif
c
c         Warn, when necessary
c
          if(global_flag_lowtemp.eq.1) then
              write(*,*) "WARNING: TOO LOW TEMP"
          endif
c
c         If matter infalling onto disk:
c         Update the mass of the central star by matter falling directly
c         from the envelope onto the star (i.e. into the capture radius
c         of the star).
c
          if(cloud_mass.gt.0.d0) then
              global_mstar = global_mstar + time_dt * ml_mdotcap
          endif
c
c
c         Update mass of the star (if active)
c
          if(isw_growstar.ne.0) then
              global_mstar = global_mstar + time_dt * abs(disk_mdot(1))
          endif
c
c         Update the total amount of outflowing angular momentum through
c         the inner and outer edge
c         NOTE 1: For the outer edge this is inaccurate, since there is no
c                 explicit no-friction BC imposed. But the purpose of this
c                 variable is anyway to make sure it is near 0. 
c         NOTE 2: I take the absolute value of this flux of angular 
c                 momentum, since both these variables are only used for
c                 testing.
c
          global_j_in  = global_j_in + abs(disk_mdot(1)) * 
     %                   sqrt(GG*disk_mgrav(1)*grid_r(1)) * 
     %                   time_dt
          global_j_out = global_j_out + abs(disk_mdot(grid_nr+1)) * 
     %                   sqrt(GG*disk_mgrav(grid_nr)*grid_r(grid_nr)) * 
     %                   time_dt
c
c         Update time
c
          time_time = time_time + time_dt
c
c         Check if we are about to drop below the mdot limit to stop the sim
c
          if((abs(disk_mdot(1))+ml_mdot).lt.time_stop_mdot) then
              write(*,*) 'Mdot dropped below critical value ',
     %                   'to stop simulation'
              simstop = .true.
          endif
c
c         Same for disk mass
c
          if(time_stop_mdisk.gt.0.d0) then
              diskmass = 0.d0
              do ir=2,grid_nr
                  diskmass = diskmass + 0.5 * 
     %                 ( disk_tprsigma(ir) + disk_tprsigma(ir-1) ) *
     %                 ( grid_r(ir) - grid_r(ir-1) )
              enddo
              if(diskmass.lt.time_stop_mdisk) then
                  write(*,*) 'Mdisk dropped below critical value ',
     %                   'to stop simulation'
                  simstop = .true.
              endif
          endif
c
c         If the accretion rate drops below a certain critical value,
c         then it may be useful to switch to constant save time stepping
c         again. This mode is by default switched off, i.e. the
c         time_savref_mdot=0.d0.
c
          if((abs(disk_mdot(1))+ml_mdot.lt.time_savref_mdot).and.
     %       (time_time.gt.time_savref_tmin).and.(.not.savref)) then
c
c             Make a message
c
              write(*,*) 'Mdot dropped below critical level: ',
     %                 'Now switching to constant time stepping'
c
c             Do this only once
c
              savref = .true.
c
c             Switch off the variable dumping time steps
c             
              varsdt = .false.
c
c             Set the new save time step and force a save dump
c
              time_tsave = time_savref_dt
              dump = .true.
          endif
c
c         Then save when necessary
c
          timesave = timesave + time_dt
          itsave   = itsave+1
          if(dump) then
c              write(*,*) 'Saving data'
              call save_data()
              timesave=0.d0
              itsave=0
          endif
c
c         If variable save time, then recompute the time_tsave
c
          if(varsdt.and.dump) then
              timesave   = timesave + time_tsave
              time_tsave = time_tsave * time_inctsave
c              write(*,*) '***********',time_tsave/year
          endif
c
c         Check if we're done
c
          if(simstop) goto 20
c
c         Reset things
c
          dump=.false.
c
      enddo
c
 20   continue
      if(.not.varsdt) then
          write(*,*) 'Saving last data'
          call save_data()
      endif
      write(*,*) 'Finished,     Time = ',time_time/year,
     %           ' year'
c
c     Done
c
      end



c     --------------------------------------------------------------
c                     UPDATE THE SIGMA FROM DIFFUSION
c
c     This is the main routine for the viscous evolution of the 
c     disk. 
c     --------------------------------------------------------------
      subroutine update_sigma_diff()
      implicit none
c
#include "common_disk.h"
c
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision null(FRSIZE_R),onee(FRSIZE_R)
      doubleprecision photev(FRSIZE_R)
      doubleprecision infall(FRSIZE_R)
      doubleprecision src(FRSIZE_R+1)
      doubleprecision backup(FRSIZE_R)
      doubleprecision r,tprsigma,signusr,temp,nu
      doubleprecision sigma_nu_sqrtr,deriv_signusr_sigr
      doubleprecision deriv_signusr_r
      doubleprecision pl,pr,ql,qr,rl,rr,dt
      doubleprecision pi
      parameter(pi=3.14159265358979324)
      integer ir,iflag,number_invalid
c
c     Set null to null
c
      do ir=1,grid_nr
          null(ir) = 0.d0
          onee(ir) = 1.d0
      enddo
c
c     Backup current tprsigma
c
      do ir=1,grid_nr
          backup(ir) = disk_tprsigma(ir)
      enddo
c
c     Set boundary conditions at inner edge
c
      pl      = 0.d0
      ql      = 2.d0
      rl      = 2.d0*global_sigma_min*2*pi*grid_r(1)
c
c     Set boundary conditions at outer edge
c
      r       = grid_r(grid_nr)
      tprsigma= disk_tprsigma(grid_nr)
      signusr = sigma_nu_sqrtr(grid_nr,r,tprsigma,temp,nu)
      pr      = deriv_signusr_sigr(grid_nr,r,tprsigma) * nu
      qr      = deriv_signusr_r(grid_nr,r,tprsigma) * nu / r
      rr      = 0.d0
c
c     First reset source term
c
      do ir=1,grid_nr
          infall(ir) = 0.d0
      enddo
c
c     Include mass loading elements 
c
      if((cloud_mass.gt.0.d0).and.(ml_implicit.eq.1)) then
          do ir=1,grid_nr
              infall(ir) = infall(ir) + 2*pi*grid_r(ir)*ml_psi(ir)
          enddo
      endif
c
c     Include photoevap sink
c
      if((isw_euv_evap.ne.0).or.(isw_fuv_evap.ne.0)) then
          do ir=1,grid_nr
c
              photev(ir) = -2*pi*grid_r(ir)*
     %                     photoev_sigdot(ir)/
     %                     disk_tprsigma(ir)
          enddo
      endif


          disk_vr(ir) = disk_mdot(ir) / ( 0.5*(disk_tprsigma(ir)+
     %                                    disk_tprsigma(ir-1)) )

      
      call cranknich_general(grid_nr,grid_r,disk_tprsigma,
     %           disk_iv0,disk_dcoef,onee,infall,photev,
     %           time_dt,pl,pr,ql,qr,rl,rr,a,b,c,q)
c
c     First check if the results are reasonable
c
      do ir=1,grid_nr
          if(number_invalid(disk_tprsigma(ir)).ne.0) then
              write(*,*) 'ERROR: Invalid disk_tprsigma()'
              write(*,*) '     ir    = ',ir
              write(*,*) '     nu    = ',nu
              write(*,*) '     alpha = ',global_alpha(ir)
              write(*,*) '     T     = ',disk_temp(ir)
              write(*,*) '     tprs  = ',disk_tprsigma(ir)
              write(*,*) pl,ql,rl,pr,qr,rr
              stop
          endif
      enddo
c
c     Check for positiveness
c
      iflag=0
      do ir=1,grid_nr
          if(disk_tprsigma(ir).lt.2.d0*pi*grid_r(ir)*
     %             global_sigma_min) then
              disk_tprsigma(ir) = 2.d0*pi*grid_r(ir)*
     %             global_sigma_min
              iflag=1
          endif
      enddo
      if(iflag.ne.0) then
c          write(*,*) 'Negative Sigma detected'
      endif
c
c     Now a-posteriori exact computation of mass flux at interfaces
c
      call aposteriori_flux_and_source(grid_nr,grid_r,
     %           backup,disk_tprsigma,disk_iv0,disk_dcoef,
     %           onee,null,null,time_dt,2,disk_mdot,src)
c
c     Compute the exact v_R from flux/(2*pi*r*sigma)
c
      do ir=2,grid_nr
          disk_vr(ir) = disk_mdot(ir) / ( 0.5*(disk_tprsigma(ir)+
     %                                    disk_tprsigma(ir-1)) )
      enddo
c      disk_vr(1)         = disk_mdot(1) / disk_tprsigma(1)
c      disk_vr(grid_nr+1) = disk_mdot(grid_nr+1) / 
c     %                     disk_tprsigma(grid_nr)
      disk_vr(1) = 0.d0
      disk_vr(grid_nr+1) = 0.d0
c
c     Done
c
      end





c     --------------------------------------------------------------
c             COMPUTE RADIAL VELOCITY AND DIFFUSION COEFFICIENT
c
c     This routine computes the radial velocity vr according to the
c     physics of Shakura-Sunyaev disks. It also computes the 
c     effective radial velocity v0 and the diffusion coefficient
c     D that you get when you write the continuity equation out
c     into a diffusion-advection type equation. Because the pure
c     advection form with vr is numerically unstable, we use the
c     diffusion-advection type equation for the model.
c     --------------------------------------------------------------
      subroutine compute_vr_v0_difcoef()
      implicit none
c
#include "common_disk.h"
c
      doubleprecision pi,sigma,mu,tillum,taues,temp_illum,r,tprsigma
      doubleprecision dumarray(FRSIZE_R),signusr,pd_r,pd_sigr,nu,temp
      doubleprecision deriv_signusr_r,deriv_signusr_sigr,sigma_nu_sqrtr
      integer ir,number_invalid
      parameter(pi=3.14159265358979324)
      parameter(tillum=3.481278d9)
c
      mu    = global_mugas
c
c     First check something
c
      if(global_sigma_min.eq.0.d0) then
         write(*,*) 'Error: Must define global_sigma_min'
         write(*,*) '       in order to keep log derivs intact'
         stop 13
      endif
c
c     For the numerical model we use v0 and D. These are found as
c     partial derivatives of Sigma * nu * sqrt(R) to (Sigma*R)_{const R}
c     and to (R)_{const Sigma*R}. This is simply using the chain rule.
c     
c                 / dlg(Sigma*nu*sqrt(R)) \           nu
c       v_0 = -3  | --------------------- |           --
c                 \        dlg(R)         /{Sigma*R}  R
c
c                 / dlg(Sigma*nu*sqrt(R)) \                        .
c       D   =  3  | --------------------- |     nu
c                 \      dlg(Sigma*R)     /{R}
c
      do ir=1,grid_nr
          r                 = grid_r(ir)
          tprsigma          = disk_tprsigma(ir) + 
     %                        2*pi*r*global_sigma_min
          signusr           = sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
          if(temp.le.0.d0) then 
              write(*,*) 'ERROR: Temperature <=0 found!'
              stop
          endif
          disk_temp(ir)     = temp
          disk_temp(ir)=min(dust_tevap,disk_temp(ir))
          disk_nu(ir)       = nu
          disk_teff(ir)     = st_tempeff

          global_alpnew(ir) = st_alpha
          pd_r              = deriv_signusr_r(ir,r,tprsigma)
          pd_sigr           = deriv_signusr_sigr(ir,r,tprsigma)
          disk_v0(ir)       = -3.d0*pd_r*nu/r
          disk_dcoef(ir)    = 3.d0*pd_sigr*nu
      enddo
c
      do ir=2,grid_nr
          disk_iv0(ir) = 0.5d0 * ( disk_v0(ir) + disk_v0(ir-1) )
      enddo
      disk_iv0(1) = 0.d0
c
c     Done
c
      end




c     --------------------------------------------------------------
c          THE CRANK-NICHOLSON SCHEME FOR DIFFUSION-ADVECTION
c
c     Perform one time step for the following PDE:
c
c        du    d /     \    d /           d  / u  \ \
c        -- + -- | u v | - -- | h(x) D(x) -- |----| | = K + L u
c        dt   dx \     /   dx \           dx \h(x)/ /
c
c     with boundary conditions
c
c         du/h |            |
c       p ---- |      + q u |       = r
c          dx  |x=xbc       |x=xbc
c
c     ARGUMENTS:
c       n             Nr of grid points
c       x             The grid
c       u             The current values of u(x)
c       vint          The values for v(x) [NOTE: vint defined at 
c                                                cell interface!]
c       d             The values for D(x)
c       h             The values for h(x)
c       k             The values for K(x)
c       l             The values for L(x)
c       dt            The time step
c       pl,pr         The p value for the lhs/rhs of the BC equation
c       ql,qr         The q value for the lhs/rhs of the BC equation
c       rl,rr         The r value for the lhs/rhs of the BC equation
c       a,b,c,q       Dummy arrays used by the cranknich() routine
c
c     NOTE: This is the new version, in which the bugs are fixed,
c           and in which a K and L term are added on the right-
c           hand-side (01-11-03).
c
c     NOTE: An even newer version (changed name to cranknich_general)
c           which also allows `relative diffusion': the h(x) factor
c           in the equation (27-01-04).
c
c     *** IMPORTANT NOTE: ***
c           This is an adapted version of the cranknich_general
c           routine in which the velocity v is given at the cell
c           boundaries!!! Specially designed for the radial mixing
c           problem... 
c           **** NOTE: **** 27-08-04
c                I changed the indexing of vint, so that it is 
c                consistent with the indexing of 'flux' in the
c                aposteriori_flux_and_source() routine!!!
c
c     --------------------------------------------------------------
      subroutine cranknich_general(n,x,u,vint,d,h,k,l,dt,
     %                     pl,pr,ql,qr,rl,rr,a,b,c,q)
      implicit none
      integer n
      doubleprecision pl,pr,ql,qr,rl,rr,dt
      doubleprecision x(n),u(n),vint(n+1),d(n),a(n),b(n),c(n),q(n)
      doubleprecision k(n),l(n),h(n)
c
      integer i
      doubleprecision eps,eps1,dum,dpl,dmn,vpl,vmn
      parameter(eps=THEEPS) ! eps=0 -> exp, eps=1 -> imp, eps=0.5 -> CN
c
      if(vint(1).ne.0.d0) then
          write(*,*) 'ERROR: This is a special version of '
          write(*,*) '   cranknich_general: v defined at interface!'
          write(*,*) '   Found vint(1).ne.0... Should not be!'
          stop
      endif
c
c     The main body of the equation
c
      eps1=1.d0-eps
      do i=2,n-1
          dum  = 2.0*dt/(x(i+1)-x(i-1))
          dpl  = dum*0.25*(h(i+1)+h(i))*(d(i+1)+d(i))/(x(i+1)-x(i))
          dmn  = dum*0.25*(h(i)+h(i-1))*(d(i)+d(i-1))/(x(i)-x(i-1))
          vpl  = dt*vint(i+1)/(x(i+1)-x(i-1))
          vmn  = dt*vint(i)/(x(i+1)-x(i-1))
          b(i) = 1.d0+eps*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))
          a(i) = -eps*dmn/h(i-1)-eps*vmn
          c(i) = -eps*dpl/h(i+1)+eps*vpl
          q(i) = u(i)+k(i)*dt
     %            -eps1*((dpl+dmn)/h(i)-vmn+vpl-dt*l(i))*u(i)
     %            +eps1*dmn*u(i-1)/h(i-1)+eps1*dpl*u(i+1)/h(i+1)
     %            +eps1*vmn*u(i-1)-eps1*vpl*u(i+1)
      enddo
c
c     The left (i=1) boundary condition
c
      b(1) = ql-pl/(h(1)*(x(2)-x(1)))
      c(1) = pl/(h(2)*(x(2)-x(1)))
      q(1) = rl
c
c     The right (i=n) boundary condition
c
      b(n) = qr+pr/(h(n)*(x(n)-x(n-1)))
      a(n) = -pr/(h(n-1)*(x(n)-x(n-1)))
      q(n) = rr
c
c     Solve the implicit differencing system
c      
      call tridag(a,b,c,q,u,n)
c
      end

c     --------------------------------------------------------------
c                    FUNCTION: SIGMA * NU * SQRT(R)
c
c     This function contains all the microphysics needed for the
c     Shakura-Sunyaev disk. It computes the Sigma * nu * sqrt(R)
c     value, given the radius R (r) and the Sigma*R (tprsigma). The
c     arguments temp and nu are meant for returning the values
c     of temperature and viscosity, in case these are needed. 
c        A crucial property of this function is that the resulting
c     value only depends on r and tprsigma! This is because the 
c     function will be called by the functions that compute the
c     double-log derivative coefficients for the linearized 
c     equations. The 'ir' argument is only there to make it
c     easier to locate where we are on the grid. 
c
c     ARGUMENTS:
c        ir                Index
c        r                 Radius in cm
c        tprsigma          Sigma*2*pi*R in g/cm
c        disk_temp         Initial guess of temperature (usually
c                          this is not necessary, but may become
c                          important for disk instabilities due to
c                          the multiple temperature solutions).
c
c     RETURNS:
c        sigma_nu_sqr      The value of Sigma*nu*sqrt(R)
c        temp              The midplane temperature of the disk
c        nu                The viscosity coefficient
c        global_alpnew     The new value of the alpha with grav inst
c     --------------------------------------------------------------
      function sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
      implicit none
      integer ir
      doubleprecision sigma_nu_sqrtr,r,tprsigma,temp,nu
c
#include "common_disk.h"
c
      doubleprecision pi,mu,zbrent,solvetemp
      parameter(pi=3.14159265358979324)
      integer number_invalid
c
      mu    = global_mugas
c
c     Compute the midplane temperature, and as a by-product compute
c     the modified alpha (modified by gravitational instability)
c
      temp  = solvetemp(ir,r,tprsigma,disk_temp(ir))
c
c     The viscosity coefficient, using the alpha value returned by 
c     the solvetemp() routine (=st_alpha)
c
      nu    = 3.1957972d11 * st_alpha * r**1.5 * temp /
     %              ( sqrt(disk_mgrav(ir)) * mu )


      if(number_invalid(nu).ne.0) then
          write(*,*) nu,ir,r,temp,disk_mgrav(ir),global_mstar
          stop
      endif
c
c     Do a check
c
      if(nu.le.0.d0) then 
          write(*,*) 'ERROR: nu zero or negative...'
          stop
      endif
c
c     Return the final value
c
      sigma_nu_sqrtr = (tprsigma/(2*pi*r)) * nu * sqrt(r)
c
      return
      end



c     --------------------------------------------------------------
c           DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(SIGMA*R)
c     --------------------------------------------------------------
      function deriv_signusr_sigr(ir,r,tprsigma)
      implicit none
      integer ir
      doubleprecision deriv_signusr_sigr,r,tprsigma
c
      doubleprecision eps,sigma_nu_sqrtr,y1,y2,dum1,dum2
      parameter(eps=1.d-4)
c
      y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2)
      y2 = sigma_nu_sqrtr(ir,r,tprsigma*(1.d0+eps),dum1,dum2)
      deriv_signusr_sigr = 2.d0*(y2-y1)/(eps*(y2+y1))
c
      return
      end


c     --------------------------------------------------------------
c             DERIVATIVE OF LOG(SIGMA NU SQRT(R)) TO LOG(R)
c     --------------------------------------------------------------
      function deriv_signusr_r(ir,r,tprsigma)
      implicit none
      integer ir
      doubleprecision deriv_signusr_r,r,tprsigma
c
      doubleprecision eps,sigma_nu_sqrtr,y1,y2,dum1,dum2
      parameter(eps=1.d-4)
c
      y1 = sigma_nu_sqrtr(ir,r,tprsigma,dum1,dum2)
      y2 = sigma_nu_sqrtr(ir,r*(1.d0+eps),tprsigma,dum1,dum2)
      deriv_signusr_r = 2.d0*(y2-y1)/(eps*(y2+y1))
c
      return
      end

c     --------------------------------------------------------------
c             DERIVATIVE OF LOG(Pressure) TO LOG(R)
c     --------------------------------------------------------------
      function deriv_pres_r(ir,r,tprsigma,cs2)
      implicit none
      integer ir
      doubleprecision deriv_pres_r,r,tprsigma,r2,cs2,hr
c
      doubleprecision eps,y1,y2,dum1,dum2,pi
      parameter(eps=1.d-4)
      parameter(pi=3.14159265358979324)
c
#include "common_disk.h"

      hr=sqrt(cs2*r**3/(6.67d-08*global_mstar))
      y1 = (tprsigma/(2.0*pi*r*2.5066283d0*hr))*cs2 ! Pressure
      r2=r*(1.d0+eps) 
      hr=sqrt(cs2*r2**3/(6.67d-08*global_mstar))
      y2 = (tprsigma/(2.0*pi*r2*2.5066283d0*hr))*cs2 ! Pressure
      deriv_pres_r = 2.d0*(y2-y1)/(eps*(y2+y1))
c
      return
      end

c     --------------------------------------------------------------
c                      FUNCTION: SOLVE TEMPERATURE 
c
c     This routine solves the local disk temperature including the
c     heating via viscous dissipation as well as the heating by
c     irradiation by the central star. As an important by-producy
c     this routine also computes the viscosity coefficient alpha
c     which may have been modified by gravitational instability.
c     --------------------------------------------------------------
      function solvetemp(ir,r,tprsigma,temp0)
      implicit none
      doubleprecision solvetemp,tprsigma,temp0,pi,kk,mp,GG,mu,r,ss
      parameter(pi  = 3.14159265358979324)
      parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
      parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
      parameter(GG  = 6.672d-8)      ! Gravitational constant
      parameter(ss  = 5.6703d-5)     ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
      integer ir
c
      doubleprecision sigma,zbrent,helptemp,temp_visc,lumbl,rdx
      doubleprecision temp_irrad,temp_selfirr,temp_blirr
      doubleprecision eps,dummy,dm1,dm2,find_temp_selfirr
      external helptemp
c
#include "common_disk.h"
c
      mu             = global_mugas
c
c     First check if the ir indeed belongs to the r
c     Just an internal consistency check.
c     
      if((ir.lt.1).or.(ir.gt.grid_nr)) stop 88324
      if(ir.lt.grid_nr) then
          if((r.lt.grid_r(ir)).or.(r.ge.grid_r(ir+1))) stop 88325
      endif
c
c     Compute the basic alpha at location r by linear interpolation
c     Note: the real alpha (including grav inst effects) is computed
c     on the fly during the temperature solving.
c
      if(ir.ne.grid_nr) then
          eps    = (r-grid_r(ir))/(grid_r(ir+1)-grid_r(ir))
          st_alpha_mri = (1.d0-eps)*global_alpha(ir)+
     %                   eps*global_alpha(ir+1)
      else
          st_alpha_mri = global_alpha(grid_nr)
      endif
c
c     Make sure to use the right Rosseland mean table for the 
c     temperature determination
c
      ist_icmp = 1
c
c     First the temperature for accretion only
c
      if(isw_qvisc.ne.0) then
c
c         ...prepare the helper function for zbrent
c
          sigma          = (tprsigma/(2*pi*r))
          st_fact_qpl    = (9./4.)*sigma*
     %                     (kk/(mu*mp))*sqrt(GG*disk_mgrav(ir)/r**3)
          if(global_gravinst.ne.0) then
              st_fact_toomre = sqrt(kk/(mu*mp))*
     %                         sqrt(GG*global_mstar/r**3)/
     %                         (pi*GG*sigma)
          else
              st_fact_toomre = 0.d0
          endif

          st_sigmakap     = sigma/global_gtd
c
c         ...call zbrent to find the temperature
c
          if(helptemp(1d0).gt.0.d0) then 
              dm1 = helptemp(1d0)
              dm2 = helptemp(1d10)
              if(dm2.gt.0.d0) then
                  temp_visc = 1d10
              else
                  temp_visc = zbrent(helptemp,1d0,1d10,1d-2)
              endif
          else 
              temp_visc=0.0d0
          endif
c
c         One more call to the temperature routine, to make sure that
c         we get the current values, both for the self-irradiation 
c         (the st_qplus) and for the gravitational instability alpha.
c
          if(temp_visc.gt.0.d0) then
              dummy      = helptemp(temp_visc)
          else
              st_qplus   = 0.d0
          endif
c
      else
          temp_visc = 0.d0
          st_qplus  = 0.d0
      endif
c
c     Now the temperature from the irradiation
c     
c     ...First the irradiation by the star
c
      if(isw_qirrad.ne.0) then 
c
c         Check if we wish the star luminosity to grow with star mass
c
          if((global_stlumgrow.gt.0.d0).and.
     %       (cloud_mass.ne.0.d0)) then
              rdx = (global_mstar/cloud_mass)**global_stlumgrow
              if(rdx.gt.1.d0) rdx=1.d0
          else
              rdx = 1.d0
          endif
c
c         Compute the irradiation temperature
c         
          temp_irrad = global_tstar*(0.5*global_flang*rdx*
     %                 (global_rstar/r)**2)**0.25
      else
          temp_irrad = 0.d0
      endif
c
c     ...Then the irradiation by the disk itself (self-irradiation)
c
      if(isw_selfirr.ne.0) then
          temp_selfirr = find_temp_selfirr(ir)
      else
          temp_selfirr = 0.d0
      endif
c
c     ...Then the irradiation by the boundary layer of the disk
c        (or in other words the `accretion shock' on the star)

c
      if(isw_lumbl.ne.0) then
          if(isw_lumbl.eq.1) then
c             
c             Use Calvet & Gullbring magnetospheric accretion model
c
              lumbl      = (1.d0-global_rstar/grid_r(1)) * 
     %                     ( GG * abs(disk_mdot(1)) * global_mstar / 
     %                       global_rstar )
              temp_blirr = (0.5*global_flang*lumbl/
     %                      (4*pi*r*r*ss))**0.25
          else
c
c             Simple boundary layer (not yet done)
c
              stop 8264
          endif
      else
          temp_blirr = 0.d0
      endif
c
c     Return the 'sum' of the two 
c
      solvetemp      = (temp_visc**4+temp_irrad**4+tempbg**4+
     %                  temp_blirr**4+temp_selfirr**4)**0.25

c
c     Compute (for use in self-irradiation) the effective temperature
c     assuming a blackbody disk.
c
      st_tempeff     = ( (st_qplus/(2*ss)) + temp_irrad**4 )**0.25
c
c     If gravitational instability mode active, then compute new alpha
c
      if(global_gravinst.ne.0) then
          st_alpha   = st_alpha_mri + st_alpha_grav 
      else
          st_alpha   = st_alpha_mri
      endif
c
      return
      end

c     --------------------------------------------------------------
c               FUNCTION: FIND ROSSELAND MEAN OPACITY
c     --------------------------------------------------------------
      function rossmean(temp,icmp)
      implicit none
      doubleprecision rossmean,temp
      integer icmp
c
#include "common_disk.h"
c
      doubleprecision eps
      integer jlo
c
c     Find the position of temp in the rosseland mean opac table
c
      call hunt(ross_temp,ross_ntemp,temp,jlo)
c
c     Check if out of bound
c     If so, then return the closest kappa
c
      if(jlo.le.0) then
          rossmean = ross_kappa(1,icmp)
          return
      endif
      if(jlo.ge.ross_ntemp) then
          rossmean = ross_kappa(ross_ntemp,icmp)
          return
      endif
c
c     Return the linear interpol of the kappa table
c
      eps    = (temp-ross_temp(jlo))/
     %         (ross_temp(jlo+1)-ross_temp(jlo))
      if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 'eps not between 0 and 1' 
      rossmean  = (1.d0-eps)*ross_kappa(jlo,icmp)+
     %            eps*ross_kappa(jlo+1,icmp)
c
      return
      end


c     --------------------------------------------------------------
c               HELPER FUNCTION FOR SOLVE TEMPERATURE
c                    FOR PURE ACCRETIONAL HEATING
c     --------------------------------------------------------------
      function helptemp(temp)
      implicit none
      doubleprecision helptemp,temp
c
#include "common_disk.h"
c
      doubleprecision rossmean,kappa,qmin,qtoomre,qcrit
      doubleprecision ss,kk,mp
      parameter(qcrit=1.d0)
      parameter(ss  = 5.6703d-5)     ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
      parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
      parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
c
c     Find the rosseland mean opacity for this temperature
c
      if(ist_icmp.le.0) stop 15543
      if(ist_icmp.gt.2) stop 15544
      kappa    = rossmean(temp,ist_icmp)
c
c     Find the Toomre parameter
c
      if(st_fact_toomre.ne.0.d0) then
          qtoomre       = st_fact_toomre * sqrt(temp)
          if(qtoomre.lt.qcrit) then
              if(qtoomre.lt.0.d0) stop 13
              st_alpha_grav = (1.d0-qtoomre)**global_ginst_plaw
          else
              st_alpha_grav = 0.d0
          endif
      else
          st_alpha_grav = 0.d0
      endif
c
c     Find the qplus for the accretion
c
      st_qplus = st_fact_qpl * temp * ( st_alpha_mri + st_alpha_grav )
c
c     Find the qmin 
c
      if(global_visheat_mode.eq.0) then

          qmin     = 8.*ss*temp**4 / (3.*st_sigmakap*kappa)

      elseif(global_visheat_mode.eq.-1) then
c
c         Absolute extreme case: T_c = T_eff. So ignore the blanket-effect
c
          qmin     = 2*ss*temp**4
      else
          stop 33543
      endif
c
c     Return the differences of the logs
c
      helptemp = log(st_qplus) - log(qmin)
      return
      end


c     --------------------------------------------------------------
c                   FUNCTION: COMPUTE SELF-RADIATION
c
c     A very simple self-irradiation recipe. We assume that the disk 
c     is optically thick. 
c     --------------------------------------------------------------
      function find_temp_selfirr(ir)
      implicit none
c
#include "common_disk.h"
c
      doubleprecision flux,ss,find_temp_selfirr
      integer ir,ir2
      parameter(ss  = 5.6703d-5)     ! Stefan-Boltzmann const  [erg/cm^2/K^4/s]
c
      flux = 0.d0
      if(ir.gt.2) then 
          do ir2=2,ir-1
              flux = flux + ss *
     %           (global_flang*(1.d0-(grid_r(ir2)/grid_r(ir))))**2 *
     %           0.5 * ( disk_teff(ir2-1)**4 + disk_teff(ir2)**4 ) *
     %           ( grid_r(ir2)**2 - grid_r(ir2-1)**2 ) / 
     %           grid_r(ir)**2
          enddo
      endif
      find_temp_selfirr = (flux/ss)**0.25
c
      end


c     --------------------------------------------------------------
c            A-POSTERIORI COMPUTATION OF FLUXES AND SOURCES
c
c     This routine computes to machine precision what the fluxes
c     at the cell interfaces were over the last time step and what
c     the source terms were. For the cell interfaces between two
c     cell centers this is straight-forward. For the cell interfaces
c     at the edges this is not straightforward, but the fluxes can
c     be computed indirectly. For this to work one has to specify
c     where the boundaries of the outermost cells are. This can be
c     directly at the location of the grid point (possibility 1) or
c     a linear extrapolation (possibility 2):
c
c          *---*-----*-----*
c          | |    |     |      (possibility 1)
c
c          *---*-----*-----*
c        |   |    |     |      (possibility 2)
c
c     The resulting flux and source are such that one can verify:
c
c                            2*dt     /             \
c       unew_i - uold_i = ----------- | F_i+1 - F_i | + dt * Src_i
c                         x_i+1-x_i-1 \             /
c
c     *** IMPORTANT NOTE: ***
c           This is an adapted version of the cranknich_general
c           routine in which the velocity v is given at the cell
c           boundaries!!! Specially designed for the radial mixing
c           problem... 
c           **** NOTE: **** 27-08-04
c                I changed the indexing of vint, so that it is 
c                consistent with the indexing of 'flux' in the
c                aposteriori_flux_and_source() routine!!!
c
c     --------------------------------------------------------------
      subroutine aposteriori_flux_and_source(n,x,uold,unew,vint,
     %                d,h,k,l,dt,ibnd,flux,src)
      implicit none
c
      integer n,ibnd
      doubleprecision dt,flux(n+1),src(n+1)
      doubleprecision x(n+1),uold(n+1),unew(n+1)
      doubleprecision vint(n+1),d(n+1),k(n+1),l(n+1),h(n+1)
c
      integer i
      doubleprecision eps,eps1,flold,flnew,v_av,d_av,uold_av,unew_av,dx
      doubleprecision dxl,dxr,h_av
      parameter(eps=THEEPS) ! eps=0 -> exp, eps=1 -> imp, eps=0.5 -> CN
c
      if(vint(1).ne.0.d0) then
          write(*,*) 'ERROR: This is a special version of '
          write(*,*) '   cranknich_general: v defined at interface!'
          write(*,*) '   Found vint(1).ne.0... Should not be!'
          stop
      endif
c
c     Compute the fluxes for i=2,n
c
      eps1=1.d0-eps
      do i=2,n
          v_av    = vint(i)
          d_av    = 0.5*(d(i-1)+d(i))
          h_av    = 0.5*(h(i-1)+h(i))
          uold_av = 0.5*(uold(i-1)+uold(i))
          unew_av = 0.5*(unew(i-1)+unew(i))
          dx      = x(i)-x(i-1)
          flold   = v_av*uold_av - h_av*d_av*
     %                   (uold(i)/h(i)-uold(i-1)/h(i-1))/dx 
          flnew   = v_av*unew_av - h_av*d_av*
     %                   (unew(i)/h(i)-unew(i-1)/h(i-1))/dx 
          flux(i) = eps1*flold + eps*flnew
          src(i)  = k(i) + l(i) * ( eps1*uold(i) + eps*unew(i) )
      enddo
c
c     For i=1 and i=n+1 (the boundaries) we know what the source
c     should be, and we know what the boundary conditions have 
c     produced (unew(1) and unew(n+1)). So we can compute implicitly
c     what the flux must have been.
c
c     First compute the source terms:
c
      src(1)   = k(1) + l(1) *  ( eps1*uold(1) + eps*unew(1) )
      src(n+1) = k(n+1) + l(n+1) *  ( eps1*uold(n+1) + eps*unew(n+1) )
c
c     Now specify the dx of the boundary cells according to possibility
c     1 or 2
c
      if(ibnd.eq.1) then
          dxl = 0.5*(x(2)-x(1))
          dxr = 0.5*(x(n)-x(n-1))
      elseif(ibnd.eq.2) then
          dxl = x(2)-x(1)
          dxr = x(n)-x(n-1)
      else
          write(*,*) 'ERROR in aposteriori_flux_and_source():'
          write(*,*) '   Dont know ibnd=',ibnd
          stop
      endif
c
c     Now find the flux at the left and right boundaries
c
      flux(1) = (dxl/dt)*(unew(1)-uold(1))-dxl*src(1)+flux(2)
      flux(n+1) = -(dxr/dt)*(unew(n)-uold(n))+dxr*src(n)+flux(n)
c
c     Done...
c
      end


c     --------------------------------------------------------------
c                   NUMERICAL RECIPES ROUTINE: TRIDAG
c     --------------------------------------------------------------
      SUBROUTINE tridag(a,b,c,r,u,n)
      INTEGER n
      DOUBLEPRECISION a(n),b(n),c(n),r(n),u(n)
      INTEGER j
      DOUBLEPRECISION bet,gam(FRSIZE_R+2)
      if(b(1).eq.0.)stop  'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.)stop  'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

c     --------------------------------------------------------------
c                        INTEGRATION ROUTINE 
c     --------------------------------------------------------------
      function integrate(n,x,f)
      implicit none
      integer n,i
      doubleprecision x(n),f(n),integrate,int
      int=0.d0
      do i=2,n
          int=int+0.5d0*(f(i)+f(i-1))*(x(i)-x(i-1))
      enddo
      if(x(n).gt.x(1)) then
          integrate = int
      else
          integrate = -int
      endif
      return
      end



c     --------------------------------------------------------------
c                COMPUTE THE PHOTOEVAPORATIVE SIGMADOT
c
c     Here a number of different photoevaporation recipes are 
c     implemented and more can be added.
c     --------------------------------------------------------------
      subroutine photoevap(nr,mstar,r,tprsigma,temp,phi,
     %                     ieuv_evap,ifuv_evap,
     %                     nbase,sigdot,izev)
      implicit none
c
#include "common_disk.h" 
#include "common_radialrt.h" 
#include "fuvG0.h"
c
      integer nr,ieuv_evap,ifuv_evap,ir
      doubleprecision mstar,phi,diskmdot
      doubleprecision r(FRSIZE_R),temp(FRSIZE_R),hp(FRSIZE_R)
      doubleprecision sigmadust(FRSIZE_R)
      doubleprecision sigdusta(NDUSTBIN,FRSIZE_R)
      doubleprecision psi(FRSIZE_R)
      doubleprecision nbase(FRSIZE_R),sigdot(FRSIZE_R)
      doubleprecision tprsigma(FRSIZE_R)
      doubleprecision teuv,cs2,rg,d0,rcr,kk,mp,GG,pi,au
      doubleprecision mug,mug_h2,sig_h,mug_fuv
      doubleprecision acst,cd,rin,nmid,rdx,tau,n0,hprin,reccoef
      doubleprecision ds,reff,rs,x
      double precision p,ri,ro,zh,nrim
      doubleprecision surfden
      integer izev(FRSIZE_R)
c
      parameter(pi=3.14159265358979324)
      parameter(au=1.496d13)
      parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
      parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
      parameter(GG  = 6.672d-8)      ! Gravitational constant
c


c     Mean molecular weight for ionized gas, and the cross section of
c     a HI atom to EUV
c
      sig_h = 6d-18     ! Taken from Osterbrock (lecture Dullemond Sept 05)
      mug_h2 = 2.3      ! The mean mol weight for bulk of the gas
      mug_fuv = 1.27
c
c     Clear the sigdot array
c
      do ir=1,nr
          sigdot(ir) = 0.d0
      enddo
c
c	ONLY DO EUV IF WIND IS THIN TO EUV
         if(abs(disk_mdot(1)).le.6.0d16) then
          if(ieuv_evap.eq.101) then
c         
c         Now do things according to our own version of Alexander 2006b
c         kind of evaporation. In contrast to Alexander we use a mugas
c         that is more consistent with ionized gas and we use the 
c         isothermal sound speed.
c         
c         First some basic numbers
c
          mug  = 0.64           ! Computed for ionized H and He mixture
          teuv = 1d4            ! Assumed temperature for EUV ionized layer
          cs2  = kk*teuv/(mug*mp)   ! We use isothermal sound speed
          rg   = GG*mstar/cs2
          rcr  = 0.1*rg         ! Critical radius (Liffmann 2003)
          rdx  = 1.d0           ! We do not include a reduction factor
c
c         Find the inner radius of the disk according to 
c         Alexander's criterion.
c         Note: sqrt(2*pi)~2.5

          tau   = 0.d0
          rin   = r(1)
          hprin = hp(1)/r(1)
          do ir=2,nr
              hp(ir) = sqrt(kk*temp(ir)*r(ir)**3/(mug_h2*mp*GG*mstar))
              nmid   = tprsigma(ir)/(2*pi*r(ir)*2.5*hp(ir)*mug_h2*mp)
              tau    = tau + nmid * sig_h * (r(ir)-r(ir-1))
              if(tau.lt.4.61d0) then 
                  rin   = r(ir)
                  hprin = hp(ir)/r(ir)
              endif
          enddo
c
c                --------- New scheme for rim --------------------------------

c                We find the Mdot_pe at the rim surface by calculating the
c                denisty from a Stromgren sphere estimate, with the radial
c               extent equal to the vertical scaleheight z_h at the hole r_h
c               We then spread this out between r_i and r_o with some steep
c               density dependence and calculate an effective sigmadot

                 p=2.5           ! Power law for density drop off from rim
                 ri=0.9*rin      ! Inner edge
                 ro=1.8* rin     ! Outer edge of photoevaporating rim
                 zh=0.5*rin      ! Scale height of neutral disk at rin
                 reccoef=2.53d-13 ! Case B recomb. coefficient
                 nrim=2.0d0*(2.0d0-p)*sqrt(zh*phi/(4.0*pi*reccoef))/
     &              (rin**p*(ro**(2.0d0-p) - rin**(2.0d0-p)))
                

c
c             Still shielding by the inner disk, so use the usual
c             photoevaporation rate
c
c             Clarke et al. 2001, Eqs. (7,8,9) but with Rg->Rcr,
c             and with 9.0d4 instead of 5.7d4 (email David).
c              n0   = 9.0d4 * sqrt(phi/1d41) * (rcr/1d14)**(-1.5)
              do ir=1,nr
                ds=0.0d0
                if(ri.le.rcr.or.r(ir).gt.ri) then
c From Font et al. 2004  nbase(ir)  = 0.754 * n0*(r(ir)/rcr)**(-2.5)
                  n0   = 9.0d4 * sqrt(phi/1d41) 
                  n0   = 0.14* sqrt(3.0*phi/ (4*pi*reccoef*rg**3.0)) 
                  nbase(ir)=n0*
     %            (2.0d0/((r(ir)/rg)**7.5+(r(ir)/rg)**12.5))**0.2
                 
c rdx is the reduction in the wind speed from alexander and armitage 2007, appendix A
c This is fit to Font et al. data
                  if(r(ir).ge.rcr) then
                    rdx=0.3423*exp(-0.3612*(r(ir)/rg -0.1))*
     %                  (r(ir)/rg-0.1)**0.2457
                      sigdot(ir) = 2*rdx*sqrt(cs2)*nbase(ir)*mp*1.35
                  else
                      sigdot(ir)=0.0
                  endif
                else
                  if(ri.gt.rcr.and.r(ir).le.ro) then
                      ds=min(nrim*(r(ir)/rin)**(-p),
     &                   phi/(sqrt(cs2)*4*3.1516*(r(ir)**2.0)))
                      sigdot(ir)=sigdot(ir)+ ds*sqrt(cs2)*mp*mug
c                    A factor of 2 already in nrim for both sides of disk
                  endif
                endif
             enddo
      endif
c
        endif

c	Skipped all above for EUV if mdot was too high !!

          do ir=1,nr
              sigdot(ir) = sigdot(ir) + photoev_isigdot(ir) ! Add profile rate
          enddo

c     ------------------------------------------------------


      do ir=1,nr
         surfden=tprsigma(ir)/(2.0*pi*r(ir))
         grsize(ir)= global_adust ! min(1.0,grsize(ir))
      enddo

      end


c     --------------------------------------------------------------
c               CALCULATE THE INFALL ONTO THE DISK SURFACE
c
c     Use Shu model coupled to Ulrich model to calculate the 
c     accretion rate onto the disk surface.
c     --------------------------------------------------------------
      subroutine shuulrich_infall(nr,r,mcloud,cs,omega,time,
     %                            psi,mdot,mdotcap,rcentr,ismooth,
     %                            smoothparam1,smoothparam2,
     %                            idistr)
      implicit none
      integer nr,nrcap,icen,ir,number_invalid,ismooth,idistr
      doubleprecision r(FRSIZE_R),psi(FRSIZE_R),rcap(FRSIZE_R)
      doubleprecision dum(FRSIZE_R),fact,integrate,dm
      doubleprecision smoothparam1,smoothparam2
      doubleprecision mcloud,cs,omega,time,mcentr,psicap
      doubleprecision m0,GG,pi
      parameter(m0     = 0.975)
      parameter(GG     = 6.672d-8)
      parameter(pi     = 3.14159265359d0)
      doubleprecision mdot,rcloud,tcloud,tcloudm,mfree,mcoll
      doubleprecision rcoll,rcoll0,mdotcap,rcentr,mmu0,rhocmueq,vz
      parameter(nrcap=30)
c
c     Calculate the global numbers for the Shu infall model with
c     rotation.
c
      mdot   = m0*cs**3/GG
      rcloud = GG * mcloud / ( 2*cs**2 )
      tcloud = rcloud/cs
      tcloudm= tcloud*(2./m0)
      if(time.lt.tcloudm) then
          mcentr = (m0/2.)*mcloud*(time/tcloud)
          if(time.lt.tcloud) then 
              mcoll  = mcloud*(time/tcloud)
          else
              mcoll  = mcloud
          endif
      else
          mcentr = mcloud
          mcoll  = mcloud
      endif
      mfree  = mcoll-mcentr
      rcoll  = cs*time 
      rcoll0 = rcoll * (m0/2.d0)
      if(time.gt.tcloud) then
          rcoll = rcloud
      endif
      if(time.gt.tcloudm) then
          rcoll0 = rcloud
      endif
      rcentr = omega**2*rcoll0**4/(GG*(mcentr+1d0))
c
c     Now include, if asked for, the smoothing-off of the infall
c
      if(ismooth.eq.0) then 
c
c         Original Shu type model: an abrupt end to the infall
c
          if(time.gt.tcloudm) then
              mdot   = 0.d0
          endif
      elseif(ismooth.eq.1) then
c
c         A linearly smoothed-off version
c
          if(smoothparam1.gt.1.d0) stop 7290
          if(smoothparam1.lt.0.d0) stop 7291
          if(time/tcloudm.gt.1.d0-smoothparam1) then
              dm = (time/tcloudm-1.d0+smoothparam1)/(2*smoothparam1)
              if(dm.lt.0.d0) stop 7292
              mdot = (1.d0-dm)*mdot
              if(mdot.lt.0.d0) mdot=0.d0
          endif
      elseif(ismooth.eq.11) then
c
c         A constant tail of infall
c
          if(time.gt.tcloudm) then
              mdot   = smoothparam1*mdot
              rcentr = smoothparam2
          endif
      elseif(ismooth.eq.21) then
c
c         A constant tail of infall
c
          if(time.gt.tcloudm) then
              mdot   = smoothparam1*mdot*(tcloudm/time)
              rcentr = smoothparam2
          endif
      endif
c
c     Compute the accretion rate onto the disk
c
      if(mdot.gt.0.d0) then
c
c         Yes, the mass-loading is active at present...
c
          icen = 0
          do ir=1,nr
              if(r(ir).lt.rcentr) icen=ir
          enddo
          psi(1) = 0.d0
          if(idistr.eq.0) then
c
c             Distribute the matter over the disk in the Ulrich way
c
              do ir=1,icen-1
                  mmu0      = sqrt(1-r(ir)/rcentr)
                  rhocmueq  = mdot*(r(ir)/rcentr)/
     %                 ((4*pi*sqrt(GG*mcentr*r(ir)**3))*(2*mmu0**2))
                  vz        = sqrt(GG*mcentr/r(ir))*mmu0
                  psi(ir)   = 2*rhocmueq*vz
              enddo
          elseif(idistr.eq.1) then
c
c             Distribute the matter over the disk in such a way that
c             there is no friction of the infalling matter with the disk
c             (see Hueso & Guillot 2005)
c
              do ir=1,icen-1
                  psi(ir)   = (mdot/(8*pi*r(ir)*rcentr)) *
     %                        ((r(ir)/rcentr)*
     %                        (1.d0-sqrt(r(ir)/rcentr)))**(-0.5)
              enddo
          else
              stop 8109
          endif
          do ir=max(icen,1),nr
              psi(ir) = 0.d0
          enddo
c
c         Compute rate of accretion inward of inner R coordinate
c         This mass is assumed to be directly loaded onto the star
c
          if(rcentr.le.r(1)) then
c
c             All matter falls within R_in, which is by definition the
c             capture radius of the central star.
c
              if(icen.ne.0) stop 19009
              mdotcap = mdot
          else
              if(nrcap.gt.FRSIZE_R) stop 71200
              if(idistr.eq.0) then
c
c                 Ulrich way
c
                  do ir=1,nrcap
                      rcap(ir)   = r(1)*ir/(1.d0*nrcap)
                      mmu0       = sqrt(1-rcap(ir)/rcentr)
                      rhocmueq   = mdot*(rcap(ir)/rcentr)/
     %                     ((4*pi*sqrt(GG*mcentr*rcap(ir)**3))*
     %                     (2*mmu0**2))
                      vz         = sqrt(GG*mcentr/rcap(ir))*mmu0
                      psicap     = 2*rhocmueq*vz
                      dum(ir)    = 2*pi*rcap(ir)*psicap
                  enddo
              elseif(idistr.eq.1) then
c
c                 Hueso & Guillot way
c
                  do ir=1,nrcap
                      rcap(ir)   = r(1)*ir/(1.d0*nrcap)
                      psicap     = (mdot/(8*pi*rcap(ir)*rcentr)) *
     %                        ((rcap(ir)/rcentr)*
     %                        (1.d0-sqrt(rcap(ir)/rcentr)))**(-0.5)
                      dum(ir)    = 2*pi*rcap(ir)*psicap
                  enddo
              else
                  stop 7204
              endif
              mdotcap = integrate(nrcap,rcap,dum)
              if(mdotcap.gt.2.*mdot) stop 90909
              if(mdotcap.lt.0.d0) then 
                  write(*,*) (rcap(ir),ir=1,nrcap)
                  write(*,*) (dum(ir),ir=1,nrcap)
                  stop 
              endif
          endif
c
c         Normalize this to make sure the mass loading is perfect
c
          do ir=1,nr
              dum(ir) = 2*pi*r(ir)*psi(ir)
          enddo
          fact = mdot / ( mdotcap + integrate(nr,r,dum) )
          if(number_invalid(fact).ne.0) then
              write(*,*) 'fact = ',fact
              stop
          endif
          if((fact.gt.1.5).or.(fact.lt.1.d0/1.5)) then
              write(*,*) 'ERROR: Something wrong with normalization'
              write(*,*) '       of the mass loading...'
              write(*,*) '       Factor wrong is: ',fact
              write(*,*) mdot,mdotcap,integrate(nr,r,dum)
              stop
          endif
          do ir=1,nr
              psi(ir) = psi(ir) * fact
          enddo
          mdotcap = mdotcap * fact
      else
c
c         Apparently there is no mass loading (anymore)
c
          mdotcap = 0.d0
          do ir=1,nr
              psi(ir) = 0.d0
          enddo
      endif
c
c     Done...
c
      end


c     --------------------------------------------------------------
c                 SAVE THE DATA BY APPENDING TO FILE
c     --------------------------------------------------------------
      subroutine save_data
      implicit none
c
#include "common_disk.h"
#include "fuvG0.h"
c
      integer ir,idust,k
      doubleprecision pi
      parameter(pi=3.14159265358979324)
c
      open(unit=1,file=outdir//'/sigma.dat',
     %  status='unknown',access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_tprsigma(ir)/
     %                (2*pi*grid_r(ir))+1.d-98
 10       format(E13.6)
 11       format(2(E13.6,1X))
      enddo
      if(global_ndust.gt.0) then
          do idust=1,global_ndust
              do ir=1,grid_nr
                  write(1,10) disk_tprsigmadust(idust,ir)/
     %                        (2*pi*grid_r(ir))+1.d-98
              enddo
          enddo
      endif
      close(1)
c
      open(unit=1,file=outdir//'/sigmadust.dat',
     %  status='unknown',access='append')
      write(1,*) 
      do ir=1,grid_nr
       write(1,'(100ES12.3)') (disk_tprsigmadust(k,ir)/
     %                (2*pi*grid_r(ir))+1.d-98,k=1,NDUSTBIN)
      enddo
c
	
      open(unit=1,file=outdir//'/fuvlum.dat',
     %  status='unknown',access='append')
	      write(1,'(E12.3)') lumbl
	close(1)

c
      open(unit=1,file=outdir//'/velo.dat',
     %  status='unknown',access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_vr(ir)
      enddo
      if(global_ndust.gt.0) then
          do idust=1,global_ndust
              do ir=1,grid_nr
                  write(1,10) disk_vrdust(ir,idust)
              enddo
          enddo
      endif
      close(1)
c
      open(unit=1,file=outdir//'/temperature.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_temp(ir)
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/visc.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_nu(ir)
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/alpha.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) global_alpnew(ir)
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/mdot.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) disk_mdot(ir)
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/sigmadotevap.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,12) photoev_sigdot(ir)+1d-90,photoev_izev(ir)
 12       format(1(E13.6,1X),1(I6,1X))
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/gastodust.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) gastodust(ir)
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/ices.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
      write(1,'(4ES16.2E4)') 
     %   sum(disk_tprsigmaice(ir,1:NDUSTBIN,1))/(2.0d0*pi*grid_r(ir)),
     %   disk_tprsigmavapor(ir,1)/(2.0d0*pi*grid_r(ir)),
     %   sum(disk_tprsigmaice(ir,1:NDUSTBIN,2))/(2.0d0*pi*grid_r(ir)),
     %   disk_tprsigmavapor(ir,2)/(2.0d0*pi*grid_r(ir))
      enddo
      close(1)
c
      open(unit=1,file=outdir//'/planetesimal.dat',status='unknown',
     %           access='append')
      write(1,*) 
      do ir=1,grid_nr
          write(1,10) planetesimal(ir)/(2.0d0*pi*grid_r(ir))
      enddo
      close(1)
c`
      open(unit=1,file=outdir//'/mstar.dat',status='unknown',
     %           access='append')
      write(1,10) global_mstar
      close(1)
c
      open(unit=1,file=outdir//'/j_outflow.dat',status='unknown',
     %           access='append')
      write(1,10) global_j_in,global_j_out
      close(1)
c
      open(unit=1,file=outdir//'/infall.dat',
     %  	status='unknown',access='append')
      write(1,11) ml_mdotcap,ml_rcentr
      close(1)
c
      open(unit=1,file=outdir//'/infall2.dat',
     %	status='unknown',access='append')
      write(1,11) ml_mdot
      close(1)
c
      open(unit=1,file=outdir//'/time.dat',
     %	status='unknown',access='append')
      write(1,10) time_time
      close(1)
c
      time_nr_save = time_nr_save + 1
      open(unit=1,file=outdir//'/time.info',status='unknown')
      write(1,*) time_nr_save
      close(1)
c
c     Now the Radiative Transfer dump
c
      if(idump_rt.gt.0) then
          if(idump_count+1.ge.idump_rt) then
              idump_isave = idump_isave + 1
              open(unit=1,file=outdir//'/rtdump.info')
              write(1,*) idump_isave
              close(1)
          open(unit=1,file=outdir//'/rtdump_times.dat',access='append')
              write(1,*) time_time
              close(1)
c Removed gastodust reset below
c              do ir=1,grid_nr
c                  gastodust(ir) = global_gtd
c              enddo
c             call write_diskstruct(idump_isave)
             idump_count=-1
          endif
          idump_count = idump_count + 1
      endif
c
c     Done
c
      end




c     --------------------------------------------------------------
c                   FUNCTION: TEST VALIDITY OF NUMBER
c
c     0 = Number is okay
c     1 = Number is INF
c     2 = Number is NAN
c     --------------------------------------------------------------
      function number_invalid(a)
      implicit none
      doubleprecision a,b,c
      integer number_invalid
      logical div,sub
c
      b=a*2.d0
      b=b/2.d0
      c=a-1.d100
c
      div = (b.eq.a)
      sub = (c.lt.a)
c
      if(div.and.sub) then
          number_invalid = 0
      elseif(div) then
          number_invalid = 1
      else
          number_invalid = 2
      endif
c
      return
      end


c     ==============================================================
c                   ROUTINES FOR PARAMETER LIST READING
c     ==============================================================


c     --------------------------------------------------------------
c            READ IN THE ENTIRE INPUT FILE INTO A STRING ARRAY
c     --------------------------------------------------------------
      subroutine read_input_file()
      implicit none
c
      common/parses/inpstring
      character*160 inpstring(INPUT_MAX_LEN)
      common/parsei/nrstring,lineok
      integer nrstring,lineok(INPUT_MAX_LEN)
c
      character*160 string
      integer iline,i
c
      do iline=1,INPUT_MAX_LEN
          read(1,'(A158)',end=10) string
          inpstring(iline) = string
          lineok(iline) = 0
      enddo
      write(*,*) 'Input file contains too many lines...'
      stop
 10   continue
      nrstring=iline-1
c      
      end



c     --------------------------------------------------------------
c         BASIC PARSING ROUTINE: SPLIT STRING INTO NAME AND VALUE
c
c     This routine interprets lines such as
c  
c       ; The following variable is useful:
c       var1   = 5.74d0  ; Here some comments
c
c     The first line is simply ignored (though marked as OK) and
c     the second line would be split into strname='var1' and
c     strvalue='5.74d0' with lenname being 4 and lenvalue being 6.
c     
c     --------------------------------------------------------------
      subroutine parse_name_value(string,strname,lenname,
     %                                   strvalue,lenvalue,icomm)
      implicit none
      character*160 string
      character*80 strname,strvalue
      integer lenname,lenvalue
      logical icomm
c
      integer ilen,i,istart,istop
c
      icomm    = .false.
      strname  = ''
      lenname  = 0
      strvalue = ''
      lenvalue = 0
      ilen     = len_trim(string)
      istart   = 1
      do i=1,ilen-1
          if(string(i:i).eq.' ') then
              istart=i+1
          else
              goto 20
          endif
      enddo
      icomm=.true.
      return
 20   continue
      if(string(istart:istart).eq.';') then
          icomm=.true.
          return
      endif
      if(istart.ge.ilen) return
      do i=istart+1,ilen-1
          if(string(i:i).eq.' ') goto 40
          if(string(i:i).eq.'=') goto 40
      enddo
      return
 40   continue
      istop=i-1
      lenname = istop-istart+1
      if(lenname.gt.80) then
          write(*,*) 'Variable name too long'
          stop
      endif
      strname = string(istart:istop)
      do i=istop+1,ilen-1
          if(string(i:i).eq.'=') goto 50
      enddo
      lenname=0
      strname=''
      lenvalue=0
      strvalue='' 
      return
 50   continue
      istart = i+1
      if(istart.gt.ilen) stop
      do i=istart,ilen
          if(string(i:i).eq.' ') then
              istart=i+1
          else
              goto 60
          endif
      enddo
 60   continue
      do i=istart+1,ilen
          if(string(i:i).eq.';') goto 80
          if(string(i:i).eq.' ') goto 80
      enddo
 80   continue
      istop = i-1
      lenvalue = istop-istart+1
      if(lenvalue.gt.80) then
          write(*,*) 'Value length too long'
          stop
      endif
      strvalue = string(istart:istop)
c
c      write(*,*) '##',strname(1:lenname),'##',strvalue(1:lenvalue),'##'
c     
      end


c     --------------------------------------------------------------
c                         STRING COMPARING
c     --------------------------------------------------------------
      function stringcompare(string1,string2,clen)
      implicit none
      integer clen
      character*80 string1
      character(len=clen):: string2
      integer i
      logical stringcompare
c
      if(clen.gt.80) stop
      if(clen.le.0) then
          stringcompare=.false.
          return
      endif
      stringcompare=.true.
      do i=1,clen
          if(string1(i:i).ne.string2(i:i)) stringcompare=.false.
      enddo
      return
      end


c     --------------------------------------------------------------
c                       PARSING ROUTINE FOR DOUBLE
c     --------------------------------------------------------------
      subroutine parse_input_double(name,iclen,value)
      implicit none
      integer iclen
      character(len=iclen) name
      doubleprecision value
c
      common/parses/inpstring
      character*160 inpstring(INPUT_MAX_LEN)
      common/parsei/nrstring,lineok
      integer nrstring,lineok(INPUT_MAX_LEN)
c
      integer iline,ichar
      character*80 strname,strvalue
      integer lenname,lenvalue,lenn
      logical stringcompare,found,icomm
c

      found=.false.
      do ichar=2,80
          if(name(ichar:ichar).eq.'@') goto 30
      enddo
      write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
      stop
 30   continue
      lenn=ichar-1
      lenn=iclen-1
      do iline=1,nrstring
          call parse_name_value(inpstring(iline),
     %             strname,lenname,strvalue,lenvalue,icomm)
          if(icomm) then
              lineok(iline) = 1
          elseif(stringcompare(strname,name,lenn).and.
     %           (lenn.eq.lenname)) then
              if(found) then
                  write(*,*) 'Found one parameter more than once'
                  write(*,*) 'in the input file: ',strname(1:lenname)
                  write(*,*) '*** ABORTING ***'
                  stop
              endif
              read(strvalue(1:lenvalue),*) value
              found=.true.
              lineok(iline) = 1
          endif
      enddo
c
      end

c     --------------------------------------------------------------
c                       PARSING ROUTINE FOR INTEGER
c     --------------------------------------------------------------
      subroutine parse_input_integer(name,iclen,value)
      implicit none
      integer iclen
      character(len=iclen) name
      integer value
c
      common/parses/inpstring
      character*160 inpstring(INPUT_MAX_LEN)
      common/parsei/nrstring,lineok
      integer nrstring,lineok(INPUT_MAX_LEN)
c
      integer iline,ichar
      character*80 strname,strvalue
      integer lenname,lenvalue,lenn
      logical stringcompare,found,icomm
c
      found=.false.
      do ichar=2,80
          if(name(ichar:ichar).eq.'@') goto 30
      enddo
      write(*,*) 'INTERNAL ERROR: Keywords must end with @!!!'
      stop
 30   continue
      lenn=ichar-1
      lenn=iclen-1
      do iline=1,nrstring
          call parse_name_value(inpstring(iline),
     %             strname,lenname,strvalue,lenvalue,icomm)
          if(icomm) then
              lineok(iline) = 1
          elseif(stringcompare(strname,name,lenn).and.
     %           (lenn.eq.lenname)) then
              if(found) then
                  write(*,*) 'Found one parameter more than once'
                  write(*,*) 'in the input file: ',strname(1:lenname)
                  write(*,*) '*** ABORTING ***'
                  stop
              endif
              read(strvalue(1:lenvalue),*) value
              found=.true.
              lineok(iline) = 1
          endif
      enddo
c
      end

c     --------------------------------------------------------------
c                   CHECK IF ALL LINES ARE OK
c     --------------------------------------------------------------
      subroutine check_all_lines_ok(ierror)
      implicit none
      logical ierror
c
      common/parses/inpstring
      character*160 inpstring(INPUT_MAX_LEN)
      common/parsei/nrstring,lineok
      integer nrstring,lineok(INPUT_MAX_LEN)
c
      character*160 string
      integer iline,ilen
c
      ierror = .false.
      do iline=1,nrstring
          if(lineok(iline).eq.0) then
              write(*,*) 'ERROR in input file: variable unknown:'
              ilen = len_trim(inpstring(iline))
              string = inpstring(iline)
              write(*,*) string(1:ilen)
              ierror = .true.
          endif
      enddo
c
      end


!-------------------------------------------------------------
!-------------------------------------------------------------
      ! REPEAT THE ADVECTION DIFFUSION EQUATION FOR DUST/VAPOR
!-------------------------------------------------------------
!-------------------------------------------------------------

c             COMPUTE RADIAL VELOCITY AND DIFFUSION COEFFICIENT
c
c     This routine computes the radial velocity and diffusion coefficient
c     for the dust particles, using an approach described in drift.pdf
c     notes, u_d and D_d are as in equation 11. 
c     --------------------------------------------------------------
      subroutine dust_compute_vr_v0_difcoef(st,ckpr,ckqr)
      implicit none
c
#include "common_disk.h"
c
      doubleprecision st(FRSIZE_R),signusr,sigma_nu_sqrtr
      doubleprecision pi,r,tprsigma,nu,temp,deriv_pres_r
      doubleprecision ckpr,ckqr,GG,cs2,mp,vkep
      doubleprecision pd_r,pd_sigr,deriv_signusr_r,deriv_signusr_sigr
      integer ir
      parameter(pi=3.14159265358979324,GG=6.672d-08)
      parameter(mp=1.6726d-24)


c     Diffusion coefficient
c       D=nu/(1.0+st*St)
c     Radial velocity
c       vr = (1.0/(St + (1.0/St)))*(cs*cs/vk)*dlnP/dlnr
c

      do ir=1,grid_nr
          r                 = grid_r(ir)
          tprsigma          = disk_tprsigma_prev(ir) + 
     %                        2*pi*r*global_sigma_min
          signusr           = sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
c  Above call returns nu
c          disk_dcoef(ir)=nu/(1.0+st(ir)**2.0)
          cs2=1.38d-16*temp/(global_mugas*mp)
          vkep=sqrt(GG*global_mstar/r)
          disk_v0(ir)=(1.0/(st(ir)+(1.0/st(ir))))*(cs2/vkep)*
     %        (1.0/(1.0+(1.0d0/gastodust(ir)))**2.0)*
     %        deriv_pres_r(ir,r,tprsigma,cs2) +
     %        disk_vr(ir)/(1.0+st(ir)**2) ! Drag component
c          dustvel(ir)=disk_v0(ir)
          dust_dcoef(ir)=disk_dcoef(ir)/(1.0+st(ir)**2.0)
      enddo
c
      do ir=2,grid_nr
          disk_iv0(ir) = 0.5d0 * ( disk_v0(ir) + disk_v0(ir-1) )
      enddo
      disk_iv0(1) = 0.d0
      
c
c     Save buoundary values for the crank-nicholson solutions
      ckpr=dust_dcoef(grid_nr)
      ckqr=disk_iv0(grid_nr)
       
c
      end

c     --------------------------------------------------------------
c                     UPDATE THE SIGMA FROM DIFFUSION
c
c     This is the main routine for the radial drift of dust 
c     in disk. 
c     --------------------------------------------------------------
      subroutine dust_update_sigma_diff(sigda,ckpr,ckqr,photoev_ape)
      implicit none
c
#include "common_disk.h"
c
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision null(FRSIZE_R),onee(FRSIZE_R)
      doubleprecision src(FRSIZE_R+1),ckh(FRSIZE_R)
      doubleprecision backup(FRSIZE_R),sigda(FRSIZE_R)
      doubleprecision photoev_ape(FRSIZE_R),photev(FRSIZE_R)
      doubleprecision r,tprsigma,signusr,temp,nu
      doubleprecision pl,pr,ql,qr,rl,rr,dt,ckpr,ckqr
      doubleprecision pi
      parameter(pi=3.14159265358979324)
      integer ir,iflag,number_invalid
c
c     Set null to null
c
      do ir=1,grid_nr
          null(ir) = 0.d0
          onee(ir) = 1.d0
          ckh(ir)= disk_tprsigma_prev(ir)
c          ckh(ir)= disk_tprsigma(ir)/(2.0d0*pi*grid_r(ir))
c          ckh(ir)=disk_tprsigma(ir)/(2.0d0*pi)
      enddo
c
c     Backup current tprsigma
c
      do ir=1,grid_nr
          backup(ir) = sigda(ir)
      enddo
c
c     Set boundary conditions at inner edge
c
      pl      = 0.d0
      ql      = 2.d0
      rl      = 2.d0*global_sigmad_min*1.0d-10*2*pi*grid_r(1)

c
c     Set boundary conditions at outer edge
c
      pr      = ckpr*ckh(grid_nr) 
      qr      = -ckqr
      rr      = 0.d0
c
c
c
c     Include photoevap sink
      if((isw_euv_evap.ne.0).or.(isw_fuv_evap.ne.0)) then
          do ir=1,grid_nr
              photev(ir) = -2.0*pi*grid_r(ir)*
     %                     photoev_ape(ir)/
     %                     sigda(ir)
          enddo
      endif
c
      
c
c     Perform Crank-Nicholson step and the sources
c

      call cranknich_general(grid_nr,grid_r,sigda,
     %           disk_iv0,dust_dcoef,ckh,null,photev,
     %           time_dt,pl,pr,ql,qr,rl,rr,a,b,c,q)
c
c     First check if the results are reasonable
c

      do ir=1,grid_nr
          if(number_invalid(sigda(ir)).ne.0) then
              write(*,*) 'ERROR: Invalid dust_tprsigma()'
              write(*,*) '     ir    = ',ir
              write(*,*) '     nu    = ',nu
              write(*,*) '     alpha = ',global_alpha(ir)
              write(*,*) '     T     = ',disk_temp(ir)
              write(*,*) '     tprs  = ',sigda(ir)
              write(*,*) pl,ql,rl,pr,qr,rr
              stop
          endif
      enddo
c
c     Check for positiveness

      iflag=0
      do ir=1,grid_nr
          if(sigda(ir).lt.2.d0*pi*grid_r(ir)*
     %             global_sigmad_min) then
              sigda(ir) = 2.d0*pi*grid_r(ir)*
     %             global_sigmad_min
              iflag=1
          endif
      enddo


       do ir=2,grid_nr
         dust_vr(ir)=disk_iv0(ir)
       enddo
      dust_vr(1) = 0.d0
      dust_vr(grid_nr+1) = 0.d0
c
c     Done
c
      end


c     --------------------------------------------------------------
c                     UPDATE THE SIGMA VAPOR FROM DIFFUSION
c
c      This is done exactly following how sigma is evolved in the 
c      disk. This does not account for the concentration gradient as in the
c      main code 
c     --------------------------------------------------------------
      subroutine update_sigma_vapor_diff(ivap)
      implicit none
c
#include "common_disk.h"
c
      doubleprecision a(FRSIZE_R),b(FRSIZE_R)
      doubleprecision c(FRSIZE_R),q(FRSIZE_R)
      doubleprecision null(FRSIZE_R),onee(FRSIZE_R)
      doubleprecision photev(FRSIZE_R),cvapor(FRSIZE_R)
      doubleprecision infall(FRSIZE_R)
      doubleprecision src(FRSIZE_R+1)
      doubleprecision backup(FRSIZE_R)
      doubleprecision r,tprsigmavapor,signusr,temp,nu
      doubleprecision sigma_nu_sqrtr,deriv_signusr_sigr
      doubleprecision deriv_signusr_r
      doubleprecision pl,pr,ql,qr,rl,rr,dt
      doubleprecision pi,dummy
      parameter(pi=3.14159265358979324)
      integer ir,iflag,number_invalid,ivap
c
c     Set null to null
c
      do ir=1,grid_nr
          null(ir) = 0.d0
          onee(ir) = 1.d0
      enddo
c
c     Backup current tprsigma and define cvapor
c
      do ir=1,grid_nr
          backup(ir) = disk_tprsigmavapor(ir,ivap)
          cvapor(ir) = backup(ir)
      enddo
c
c     Set boundary conditions at inner edge
c
      pl      = 0.d0
      ql      = 2.d0
      rl      = 2.d0*global_sigma_min*2*pi*grid_r(1)
c
c     Set boundary conditions at outer edge
c
      r       = grid_r(grid_nr)
      tprsigmavapor= disk_tprsigmavapor(grid_nr,ivap)
      signusr = sigma_nu_sqrtr(grid_nr,r,tprsigmavapor,temp,nu)
      pr      = deriv_signusr_sigr(grid_nr,r,tprsigmavapor) * nu
      qr      = deriv_signusr_r(grid_nr,r,tprsigmavapor) * nu / r
      rr      = 0.d0
c
c     First reset source term
c
      do ir=1,grid_nr
          infall(ir) = 0.d0
      enddo
c
c     Include photoevap sink
c
      if((isw_euv_evap.ne.0).or.(isw_fuv_evap.ne.0)) then
          do ir=1,grid_nr
c  This term is relative to sigma 
               photev(ir) = -2*pi*grid_r(ir)*
     %                     photoev_sigdot(ir)/
     %                     disk_tprsigma(ir)
c               photev(ir)=0.0

          enddo
      endif

       disk_vr(ir)=disk_mdot(ir)/( 0.5*(disk_tprsigma(ir)+
     %                             disk_tprsigma(ir-1)))

c
c     Perform Crank-Nicholson step and the sources
      
      call cranknich_general(grid_nr,grid_r,cvapor,
     %           disk_iv0,disk_dcoef,onee,infall,photev,
     %           time_dt,pl,pr,ql,qr,rl,rr,a,b,c,q)
c
c     First check if the results are reasonable, re-compute from cvapor
c
      do ir=1,grid_nr
          disk_tprsigmavapor(ir,ivap)=cvapor(ir)
          if(number_invalid(disk_tprsigmavapor(ir,ivap)).ne.0) then
              write(*,*) 'ERROR: Invalid disk_tprsigmavapor()'
              write(*,*) '     ir    = ',ir
              write(*,*) '     nu    = ',nu
              write(*,*) '     alpha = ',global_alpha(ir)
              write(*,*) '     T     = ',disk_temp(ir)
              write(*,*) '     tprs  = ',disk_tprsigmavapor(ir,ivap)
              write(*,*) pl,ql,rl,pr,qr,rr
              stop
          endif
      enddo
c
c     Check for positiveness
c
      iflag=0
      do ir=1,grid_nr
          if(disk_tprsigmavapor(ir,ivap).lt.2.d0*pi*grid_r(ir)*
     %             global_sigma_min) then
              disk_tprsigmavapor(ir,ivap) = 2.d0*pi*grid_r(ir)*
     %             global_sigma_min
              iflag=1
          endif
      enddo

      end

c---------------------------------------------------------------
c computing diffusion coefficients for vapor component

      subroutine vapcompute_vr_v0_difcoef(ivap)
      implicit none
c
#include "common_disk.h"
c
      doubleprecision pi,sigma,mu,tillum,taues,temp_illum,r
      doubleprecision tprsigmavapor,tprsigma
      doubleprecision dumarray(FRSIZE_R),signusr,pd_r,pd_sigr,nu,temp
      doubleprecision deriv_signusr_r,deriv_signusr_sigr,sigma_nu_sqrtr
      integer ir,number_invalid,ivap
      parameter(pi=3.14159265358979324)
      parameter(tillum=3.481278d9)
c
      mu    = global_mugas
c
      do ir=1,grid_nr
        r                 = grid_r(ir)
        tprsigma          = disk_tprsigmavapor(ir,ivap) + 
     %                     2*pi*r*global_sigma_min
        signusr           = sigma_nu_sqrtr(ir,r,tprsigma,temp,nu)
        if(temp.le.0.d0) then 
            write(*,*) 'ERROR: Temperature <=0 found!'
            stop
        endif
        disk_temp(ir)     = temp
        disk_temp(ir)=min(dust_tevap,disk_temp(ir))
        disk_nu(ir)       = nu
        disk_teff(ir)     = st_tempeff
        global_alpnew(ir) = st_alpha
        pd_r              = deriv_signusr_r(ir,r,tprsigma)
        pd_sigr           = deriv_signusr_sigr(ir,r,tprsigma)
        disk_v0(ir)       = -3.d0*pd_r*nu/r

c        Using exactly the same vr and D as for gas, but the vr is increased for vapor
c        disk_v0(ir)=disk_v0(ir)*disk_tprsigma(ir)/disk_tprsigmavapor(ir,ivap)

        disk_dcoef(ir)    = 3.d0*pd_sigr*nu
      enddo

      do ir=2,grid_nr
          disk_iv0(ir) = 0.5d0 * ( disk_v0(ir) + disk_v0(ir-1) )
      enddo
      disk_iv0(1) = 0.d0

      end


      subroutine geticefrac(sig,h,T,rad,ficew,ficeco)

          double precision sig,h,T,ficew,ficeco,rad
          double precision mu,xsnow,dummy,fac,ncrit,nmid
          double precision, parameter :: pi=acos(-1.0d0)

          ! sig =  2 pi r sigma; x = z/sqrt(2) h
          ! midplane number density at this radius
          nmid = sig/(3.8d-24*h*2.0*sqrt(2.0*pi))/(2.0*pi*rad)

          ! WATER
          !Vapor pressure for water = 2.53e13 * exp(-6070/T) dynes/cm2
          ficew=-1.0
          ncrit = (2.53d13*exp(-6070.0/T)) / (1.38d-16*T)
          ! assume that snowline is always at 0.1r
          if(nmid*1.8d-4>ncrit) then
                xsnow=0.1*rad/(sqrt(2.0)*h)
                dummy=xsnow
                ficew = erf(xsnow)
          else
                ficew=-1.0  ! inside radial snowline
          endif

          ! CO
          !Vapor pressure for CO = 7.714e11 * exp(-1030/T) dynes/cm2
          ficeco=-1.0
          ncrit = (7.714d11*exp(-1030.0/T)) / (1.38d-16*T)
          !  assume that snowline is always at 0.1r
          if(nmid*1.4d-4>ncrit) then
                xsnow=0.1*rad/(sqrt(2.0)*h)
                xsnow = dummy ! vertically CO freezes with water
                ficeco = erf(xsnow)
          else
                ficeco=-1.0  ! inside radial snowline
          endif
         
          return
      end subroutine geticefrac






