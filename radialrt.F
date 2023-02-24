c     ==============================================================
c        ROUTINES THAT PERFORM RADIAL RADIATIVE TRANSFER IN DISKS
c
c     This program has now been converted in a set of subroutines
c     for use with the disk evolution code with FUV-photoevaporation.  
c     Moreover, I have removed the flexible z-grid: I now assume 
c     that z propto R, but one can use log grid for z. 
c     ==============================================================

#include "main.h"

c#define GASTODUST (100.d0)
#define TSOLVE_FAST (1)

c     --------------------------------------------------------------
c                       INITIALIZATION ROUTINE
c     --------------------------------------------------------------
      subroutine radialrt_init(ntemp,temp0,temp1)
      implicit none
      integer ntemp
      doubleprecision temp0,temp1
c
#include "common_dust.h"
c
      call read_all_data()
      call make_star_spectrum()
#ifdef TSOLVE_FAST
      call init_solvetempfast(ntemp,temp0,temp1)
#endif
c
      dust_agrain =  0.3*1.0d-4
      dust_sw     = 3.6d0
      write(*,*) 'NOTE: I TAKE AGRAIN=',dust_agrain/1.0d-4,'um 
     +       AND SPECIFIC WEIGHT=',dust_sw
c
      end


c     --------------------------------------------------------------
c                   MAKE AND COMPUTE 1+1D MODEL
c
c     This routine computes the dust radiative transfer in a simple
c     way, in which the radial rays are simply extincted outward,
c     and a midplane temperature is assumed to be the minimum
c     dust temperature at each radius. It will calculate the dust
c     temperature at each r and z. To make the RT easier, the 
c     z-grid is always scaled with r, so that radial outward 
c     radiation (star is assumed to be a point-source) is always
c     going from one grid point to the next, without need of 
c     interpolation in z. 
c     
c     The routine will also, upon request, compute the gas temperature
c     using Uma Gorti's gas temperature routines. Otherwise it is
c     assumed that T_g=T_d. 
c
c     ARGUMENTS:
c
c       nz          Nr of z-grid points requested
c       sigmadust   The dust surface density Sigma_dust(R)
c       hp          The pressure scale height as a function of
c                   radius, based on the midplane temperature 
c                   computed with the flaring-angle recipe.
c       tmid        The midplane temperature that follows from the
c                   flaring-angle recipe. This number also acts
c                   as the lower boundary to the dust temperature
c                   that will be calculated here. This is a way
c                   to mimic the diffuse radiation field present
c                   in the disk.
c       zrmin       Minimum z/r of the z/r grid used.
c       zrmax       Maximum z/r of the z/r grid used.
c       ifuv_evap   =-1 means: T_g=T_d
c                   =1 means: use the routines of Uma Gorti
c       iter_nr     Max nr of iterations for vertical structure
c                   (=0 no iteration)
c       iter_crit   The convergence criterion for the vertical
c                   structure calculations.
c
c
c     NOTE: For iter>0 we assume that the density of
c           the previous iteration is still in memory! It will
c           check if a reasonable density is present, and if
c           not, it will create a new one. This is very important
c           in the first time step, but also if for some magic
c           reason a part of the disk that was originally empty
c           now suddenly gets refilled.
c
c     NOTE: In principle the midplane temperature could be calculated
c           using the shape of the disk that comes out of the 
c           calculations, but this is very dangerous since it could
c           easily give rise to numerical instabilities.
c     --------------------------------------------------------------
      subroutine radialrt_make_11d_model(nz,sigmadust,hp,
     %              tmid,zrmin,zrmax,ifuv_evap,iter_nr,
     %              iter_crit)
      implicit none
c
#include "common_disk.h"
#include "common_grid.h"
#include "common_radialrt.h"
#include "common_natconst.h"
c
      integer nr,nz,nzlog,iz,ir,ifuv_evap,iter,iter_nr,inu,inuv
      doubleprecision zrmin,zrmax,dumorig,dum,iter_crit
      doubleprecision sigmadust(FRSIZE_R),find_dust_kappa
      doubleprecision hp(FRSIZE_R),tmid(FRSIZE_R)
      doubleprecision const,dummy,cursig,error,errloc,epsv
      doubleprecision gravc,grv,dlgt,dz,sigdum,rhoold(FRSIZE_Z)      
      logical iredo_initdens,repeat
c
c     A constant
c   
      const=sqrt(2*pi)
c
c     Check
c
      if(nz.gt.FRSIZE_Z) stop 912
      nr = grid_nr
c
c     Get the dust opacity
c
      do inu=1,freq_nr
          dstkap_kappa(inu) = find_dust_kappa(inu,1,1,star_t,1,0)
      enddo
c
c     Find kappa_V
c
      inuv=-1
      epsv=-1.d0
      do inu=freq_nr-1,1,-1
          if((freq_nu(inu+1)-2.9979d14/0.55)*
     %       (freq_nu(inu)-2.9979d14/0.55).le.0.d0) then
              inuv=inu
              epsv=(2.9979d14/0.55-freq_nu(inuv))/
     %             (freq_nu(inuv+1)-freq_nu(inuv))
          endif
      enddo
      if(inuv.lt.1) then
          write(*,*) 'ERROR: Freq-grid does not include V-band'
          stop 8343
      endif
      if((epsv.lt.0.d0).or.(epsv.ge.1.d0)) stop 8344
      dstkap_kappav = (1.d0-epsv)*dstkap_kappa(inuv) + 
     %                epsv*dstkap_kappa(inuv+1)
c
c     First compute how to divide the nz point into the log-regime
c     and the lin regime. 
c
      dumorig=-1d6
      do iz=2,nz-2
          dum = 1.d0/((zrmax/zrmin)**(1.d0/(1.0*iz))-1.d0)-
     %          1.d0*(nz-iz)
          if(dum*dumorig.lt.0.d0) goto 10
          dumorig = dum
      enddo
c     stop 1122
 10   continue
      nzlog = iz
c
c     Now create the global zr array (z/H_p-array, which is 
c     independent of r)
c
c     ...the log part
c
      do iz=nz-nzlog+1,nz
          grid_zr(iz) = zrmin * (zrmax/zrmin)**
     %              ((iz-(nz-nzlog+1)+1)/(1.d0*nzlog))
      enddo
c
c     ...the lin part
c
      do iz=1,nz-nzlog
          grid_zr(iz) = iz*zrmin/(1.d0*(nz-nzlog))
      enddo
c
c     Now the cell interfaces
c
      do iz=2,nz
          grid_zri(iz) = 0.5 * ( grid_zr(iz) + grid_zr(iz-1) )
      enddo
      grid_zri(1) = 0.d0
      grid_zri(nz+1) = grid_zr(nz) + 0.5*(grid_zr(nz)-grid_zr(nz-1))
c
c     Now the integration weights
c
      do iz=1,nz
          grid_dzr(iz) = grid_zri(iz+1) - grid_zri(iz)
      enddo
c
c     Put the nz information in the commons.....
c
      grid_nz = nz
c
c     Loop over ir

c
      do ir=1,nr
c
c         Compute the surface density of the current density profile
c         
          if(iter_nr.gt.0) then
              cursig = 0.d0
              do iz=1,nz
c                 BUGFIX 28.09.07
c                  cursig = cursig + rrt_rhogas(iz,iz) * grid_dzr(iz)
                  cursig = cursig + rrt_rhogas(iz,ir) * grid_dzr(iz)
              enddo
              cursig = grid_r(ir) * cursig * 2 ! Two sides
          else
              cursig = -1d99   ! Safety number
          endif
c
c         Check if reasonable original density available
c
          if(iter_nr.gt.0) then 
              if(abs(cursig/(gastodust(ir)*sigmadust(ir))-1.d0)
     %                  .gt.0.2) then 
                  iredo_initdens = .true.
              else
                  iredo_initdens = .false.
              endif
          endif
c
c         Initial guess: put the rho_gas onto this grid
c
          if((iter_nr.eq.0).or.iredo_initdens) then
c
c             Make completely new (initial) guess of density
c             structure
c
              do iz=1,nz
                  rrt_rhogas(iz,ir) = gastodust(ir) * ( sigmadust(ir) /
     %                 (2.5066283d0*hp(ir)) ) *
     %                 exp(-0.5d0*(grid_zr(iz)*grid_r(ir)/hp(ir))**2)
              enddo
          else
c
c             Take earlier density structure, but renormalize to
c             the right surface density
c
              dummy = (gastodust(ir)*sigmadust(ir))/cursig
              do iz=1,nz
                  rrt_rhogas(iz,ir) = dummy * rrt_rhogas(iz,ir)
              enddo
          endif
c
c         Now do the radiative transfer
c
          repeat = .false.
          call integrate_radial_rt_onestep(ir,repeat,tmid(ir),
     %              ifuv_evap)
c
c         If iter_nr gt 0 then do iteration of the vertical structure
c         Otherwise, we have reached our solution
c
          do iter=1,iter_nr
c
c             Buffer old density
c
              do iz=1,nz
                  rhoold(iz) = rrt_rhogas(iz,ir)
              enddo
c
c             Compute the vertical density structure
c
              gravc    = (GG*global_mstar/grid_r(ir)**3)*2.3*mp/kk
              rrt_rhogas(1,ir) = 1.d0
              do iz=1,nz-1
                  dz   = grid_r(ir)*(grid_zr(iz+1)-grid_zr(iz))
                  dlgt = (log(rrt_tgas(iz+1,ir))-
     %                    log(rrt_tgas(iz,ir)))/dz
                  grv  = 0.5 * gravc * grid_r(ir) * 
     %                 ( grid_zr(iz+1)/rrt_tgas(iz+1,ir) + 
     %                   grid_zr(iz)/rrt_tgas(iz,ir) )
                  rrt_rhogas(iz+1,ir) = 
     %                 rrt_rhogas(iz,ir) * exp ( -dz*(grv+dlgt) )
              enddo
              sigdum = 0.d0
              do iz=1,nz
                  sigdum = sigdum + grid_r(ir) * grid_dzr(iz) * 
     %                 rrt_rhogas(iz,ir)
              enddo
              do iz=1,nz
                  rrt_rhogas(iz,ir) = rrt_rhogas(iz,ir) * 
     %                    gastodust(ir)*sigmadust(ir) / sigdum
              enddo
c
c             Recompute the dust temperature
c              
              repeat = .true.
              call integrate_radial_rt_onestep(ir,repeat,
     %                    tmid(ir),ifuv_evap)
c
c             Compute the error
c             Only do when the density is greater than 100, otherwise ignore (Uma)
c
              error = 0.d0
              do iz=1,nz
                 if(rrt_rhogas(iz,ir).ge.3.0d-22) then
                  errloc = abs(rrt_rhogas(iz,ir)-rhoold(iz)) / 
     %                     abs(rrt_rhogas(iz,ir)+rhoold(iz))
                  if(errloc.gt.error) error=errloc
                  endif
              enddo
c        write(*,*) '   Error = ',error,iter,iter_crit,iter_nr,grid_r(ir)
              if(error.lt.iter_crit) goto 15
          enddo
 15       continue
      enddo   ! Next ir
c
      end


c     --------------------------------------------------------------
c                 COMPUTE VERTICAL PHOTOEVAPORATION
c     --------------------------------------------------------------
      subroutine compute_photoevap(rin,mugas,psi,izev)
      implicit none
c     
      doubleprecision mugas,rin
      doubleprecision psi(FRSIZE_R)
      integer izev(FRSIZE_R)
c     
#include "common_disk.h"
#include "common_grid.h"
#include "common_dust.h"
#include "common_radialrt.h"
c
      doubleprecision GG,mp,kk
      parameter(GG  = 6.672d-8)      ! Gravitational constant
      parameter(mp  = 1.6726d-24)    ! Mass of proton          [g]
      parameter(kk  = 1.3807d-16)    ! Bolzmann's constant     [erg/K]
c
      doubleprecision evap,rg,rs,cs,x,ds,reff,taucol,pi
      doubleprecision den,dustden,zz,colstar,avstar,ri,ro,mgrain
      doubleprecision sigdotrim,fzh,taucol0,zh,mdotrim,vflow,ape
      doubleprecision ape0
      integer ir,iz,irim,nzone
c
      psi(1)  = 0.d0  ! Put to zero by hand at inner edge...
      izev(1) = 0
      do ir=2,grid_nr
          psi(ir)  = 0.d0
          izev(ir) = 0
          do iz=1,grid_nz
c
c             Calculate the `local' gravitational radius
c
              rg = GG * star_m * mugas * mp / ( kk*rrt_tgas(iz,ir))
ccccc
ccccc             Now the sonic radius
ccccc
cccc              rs = 0.5 * rg
c
c             Sound speed
c
              cs = sqrt(kk*rrt_tgas(iz,ir)/(mugas*mp))
c
c             Basic evaporation rate
c             (Factor 2 for two sides!)
c
              evap = rrt_rhogas(iz,ir) * cs * 2.d0 
c
c             Now compute the correction factor of Uma Gorti,
c             which replaces the original formulae.
c
              reff = grid_r(ir)
              x    = 1.0 - 8.0 * reff / rg
              if(x>0.0d0) then
                  rs = 0.25*rg*(1+sqrt(x))
                  ds=dexp((-0.5d0*rg/reff)*(1.0-reff/rs)**2.0)
     &                 *(rs/reff)**2.0
                  if(ds.gt.1.0d0) ds=1.0d0 ! Subsonic launch only
              else
                  if(reff.ge.rg) then
                      rs=reff
                      ds=1.0
c                  elseif(reff.ge.0.5*rg) then
                  else
                      rs=rg*(1.0+6.0*(reff/rg -1.0d0)/7.0)
                      ds=dexp((-0.5d0*rg/reff)*
     &                     (1.0-reff/rs)**2.0)
     &                     *(rs/reff)**2.0
                      if(ds.gt.1.0d0) ds=1.0d0 ! Subsonic launch only
                  endif
              endif
              evap = evap * ds

c   Calculation of the stopping time and dust size that is retained (a>ape)

c      Flow velocity from conservation of mass
             vflow=ds*cs ! ds already has the area factor (rs/reff)^2 above
c      If there is any flow then find ape, otherwise 
             vflow=vflow+cs
c      Epstein Drag force = (4 pi/3) rhogas * a^2 * cs * vflow 
c      Vertical grav force = G M* (z/reff^3) 4pi/3 rho_gr a^3 
c      Equating both gives ape
c      ape = rhogas *cs* vflow / (G M* (z/reff^3) rho_gr)
c      1/(G * rho_gr) = 5.55e6 

        ape = 5.55d6*rrt_rhogas(iz,ir)*cs*vflow
     &          *grid_r(ir)**2.0/(global_mstar*grid_zr(iz))

c WHAT IF YOU DO NOT ALLOW DUST TO CAUSE PHOTOEVAPORATION?
          if((rrt_tgas(iz,ir)-rrt_tdust(iz,ir)).le.1.d0)  evap=0.0d0
          if((rrt_tgas(iz,ir)-rrt_tdust(iz,ir)).le.1.d0)  ape=0.0d0

c             If evap is larger then psi, then put psi=evap
              if(evap.gt.psi(ir)) then
                  ape0=ape
                  psi(ir)  = evap
                  izev(ir) = iz
              endif
          enddo ! Z loop

c       default minimum grain size if 50A or 5.0d-7
c       fpe = fraction of dust mass carried with pe flow=(sqrt(ape)-sqrt(amin))/
c             (sqrt(amax)-sqrt(amin)) = (sqrt(ape*amin)-amin)/( <a>-amin)
c        if(ape0.le.5.0d-7) then
c          fdustpe(ir)=0.0
c        elseif(ape0.ge.sqrt(grsize(ir)**2.0/5.0d-7)) then
c          fdustpe(ir)=1.0
c        else
c          fdustpe(ir)=(sqrt(ape0*5.0d-7)-5.0d-7)/(grsize(ir)-5.0d-7)
c        endif
c        if(izev(ir).lt.1) fdustpe(ir)=1.0
c       ape0=0.0
        fdustpe(ir)=1.0 ! default: only grains >1cm are retained
        if(izev(ir).gt.1) fdustpe(ir) = ape0 !  

      enddo ! R loop
c     

c Rim photoevaporation. To maintain clarity in code, some of
c the relevant subroutines and sections appear below again.
 
c  we assume z_h is the scale height to a column of 1.0d22 (corresponding
c  Av=1 for our chosen dust opacity. We assume a radial scale lenght
c fiven by fz_h, where f<1, typically ~ 0.3

          fzh=0.3
          taucol0=1.0d22
          pi=3.14159

c         Check if hole has formed, analogous to EUV case
          taucol=0.0d0
          irim=0
          do ir=1,grid_nr
            taucol=taucol+(grid_r(ir)-grid_r(ir-1))
     %          *rrt_rhogas(1,ir-1)/(1.273*mp)
            if(ir.gt.2.and.taucol.le.taucol0) irim=ir
          enddo
          if(irim.gt.0)print*,irim,'irim',taucol,taucol0,'taucol@rim'
          if(irim.gt.2) then
              ir=irim
c             Find the height where the radial column is taucol0
              zh=0.0d0
              do iz=2,grid_nz
                 if(radf_colrad(iz).ge.taucol0) 
     %                    zh=grid_r(ir)*grid_zr(iz)
              enddo 

c             --------------------------------------------------------------
c             Find sigdot for each iz and sum over area to get total mdot off the rim. 
              mdotrim=0.0d0
c              mgrain   = (4*pi/3.)*dust_sw*dust_agrain**3     
              mgrain   = (4*pi/3.)*dust_sw*grsize(ir)**3     
              do iz=2,grid_nz
c             Convert to number density
              den      = rrt_rhogas(iz,ir)/(1.273*mp)
              den=min(den,taucol0/(fzh*zh))
              dustden  = ((den*1.273*mp)/gastodust(ir)) / mgrain
c             Call routine Gorti & Hollenbach
              zz = grid_r(ir) * grid_zr(iz)
c              call GasT2(grid_r(ir),zz,den,dustden,grsize(ir),
c     %          rrt_tdust(iz,ir),radf_colrad(iz),
c     %                radf_colvert(iz),radf_av(iz),rrt_tgas(iz,ir))
              colstar=min(taucol0,radf_colrad(iz)) 
              avstar=colstar/taucol0 
              call GasT2(grid_r(ir),zz,den,dustden,grsize(ir),
     %          rrt_tdust(iz,ir),colstar,colstar,avstar,
     %                rrt_tgas(iz,ir))
              
c             Repeat the photoevaporation calculation above at the beginning
c             of this subroutine, reinserted here without comments for clarity. 
c
              rg = GG * star_m * mugas * mp / ( kk*rrt_tgas(iz,ir))
              cs = sqrt(kk*rrt_tgas(iz,ir)/(mugas*mp))
              evap = (den*1.273*mp)    * cs * 2.d0  ! Basic rate
              reff = grid_r(ir)
              x    = 1.0 - 8.0 * reff / rg
              if(x>0.0d0) then
                  rs = 0.25*rg*(1+sqrt(x))
                  ds=dexp((-0.5d0*rg/reff)*(1.0-reff/rs)**2.0)
     &                 *(rs/reff)**2.0
                  if(ds.gt.1.0d0) ds=1.0d0 ! Subsonic launch only
              else
                  if(reff.ge.rg) then
                      rs=reff
                      ds=1.0
c                  elseif(reff.ge.0.5*rg) then
                  else
                      rs=rg*(1.0+6.0*(reff/rg -1.0d0)/7.0)
                      ds=dexp((-0.5d0*rg/reff)*
     &                     (1.0-reff/rs)**2.0)
     &                     *(rs/reff)**2.0
                      if(ds.gt.1.0d0) ds=1.0d0 ! Subsonic launch only
                  endif
              endif
              evap = evap * ds

              mdotrim=mdotrim+ evap*2.0*pi*grid_r(ir)*
     &                      ( grid_r(ir)*(grid_zr(iz)-grid_zr(iz-1))) ! This is dz for the cell
              enddo
c             --------------------------------------------------------------

c          Now distribute this total rim mass loss over a small radial zone similar to
c          that for the EUV, here we choose 10 radial zones
           nzone=10
           ri=grid_r(irim)
           ro=grid_r(irim+nzone)
c          mdotrim is for both surfaces, so is sigdotrim
           sigdotrim= mdotrim/(pi*(ro*ro -ri*ri))
     
           do ir=irim,(irim+nzone),1
                  psi(ir)  = psi(ir) + sigdotrim
           enddo
        endif
          
           print*,'--------- Rim at -----------',grid_r(irim)/1.5d13
      end



c     --------------------------------------------------------------
c                          READ ALL THE DATA
c     --------------------------------------------------------------
      subroutine read_all_data()
      implicit none
c
#include "common_disk.h"
#include "common_grid.h"
#include "common_radialrt.h"
#include "common_dust.h"
c
      integer iformat
      integer ir,iz,nf,inu,i
      doubleprecision a,b,c,d,e
      doubleprecision frq,dummy
      logical fexf,fexo,fexs
c
c     Read the star information
c
c      open(unit=1,file='starinfo.inp',status='old')
c      read(1,*) iformat
c      read(1,*) star_r
c      read(1,*) star_m
c      read(1,*) star_t
c      close(1)
      star_r = global_rstar
      star_m = global_mstar
      star_t = global_tstar
c c
c c     Read the current density structure
c c
c       open(unit=1,file='diskstruct_0.dat',status='old')
c       read(1,*) nr,nz
c       if(nr.gt.FRSIZE_R) then
c           write(*,*) 'ERROR in radialrt:'
c           write(*,*) 'R grid too large'
c           stop 13
c       endif
c       if(nz.gt.FRSIZE_Z) then
c           write(*,*) 'ERROR in radialrt:'
c           write(*,*) 'Z grid too large'
c           stop 13
c       endif
c       do ir=1,nr
c           do iz=1,nz
c               read(1,*) a,b,c,d,e
c               grid_r(ir) = a
c               grid_z(iz,ir) = b
c               rrt_rhodust(iz,ir) = c / GASTODUST
c               rrt_vz(iz,ir) = e
c           enddo
c       enddo
c       close(1)
c
c     Read the frequency array and the opacities
c
      inquire(file='./input/kappaplanck.inp',exist=fexf)
      if(fexf) then
          write(*,*) 'PROBLEM: kappaplanck.inp is obsolete...'
          write(*,*) 'TIP    : mv kappaplanck.inp planckopac.inp'
          write(*,*) '         and add two colums with 0.d0'
          write(*,*) '         (i.e. kap_scat and g_scat).'
          stop 13
      endif
      inquire(file='./input/planckopac.inp',exist=fexo)
      if(fexo) then
c
c         Planck Opacity file is there, so use it. If the frequency
c         dependent opacities are also there, then they are assumed
c         to be the opacities from which the planck opacities were
c         derived. The frequency-dependent opacitities can then be
c         used to make the SED.
c     
          freq_nr = 1
          call read_planckopac()
c
c         Check if also dust opacs present
c
          inquire(file='./input/frequency.inp',exist=fexf)
          inquire(file='./input/dustopac_1.inp',exist=fexo)
          if(fexf.and.fexo) then
              write(*,*) '   (Read planck opacity from planckopac.inp)'
          endif
      else
c
c         Find the frequency-dependent opacities, and do the radial
c         transfer with these opacitities (do it with frequency
c         dependence).
c
          inquire(file='./input/frequency.inp',exist=fexf)
          inquire(file='./input/dustopac_1.inp',exist=fexo)
          if((.not.fexf).or.(.not.fexo)) then
              write(*,*) 'Need both frequency.inp and dustopac_1.inp'
              stop 13
          endif
          open(unit=1,file='./input/frequency.inp',status='old')
          read(1,*) nf
          if(nf.gt.FRSIZE_FREQ) then
              write(*,*) 'ERROR in radialrt:'
              write(*,*) 'Frequency grid too large'
              stop 13
          endif
          do inu=1,nf
              read(1,*) a
              freq_nu(inu) = a
          enddo
          dust_setup_nrspecies=1
          dust_setup_nrsizes(1)=1
          freq_nr = nf
c
c         Write a reminder:
c
          write(*,*) '------------------------------------------'
          write(*,*) 'REMINDER: Make sure that the frequency-'
          write(*,*) '  dependent opacities are consistent with'
          write(*,*) '  the rosseland mean opacities...'
          write(*,*) '------------------------------------------'
c
c         Reading the dust data.
c
          write(*,*) 'Reading dust data...'
          call read_dustdata()
          write(*,*) 'done'
      endif
c
c     Read (if required) the input spectrum
c
      if(star_t.le.0.d0) then
          inquire(file='./input/starspectrum.inp',exist=fexs)
          if(.not.fexs) then
              write(*,*) 'ERROR: Must specify either star_t '
              write(*,*) '       or starspectrum.inp'
              stop
          endif
          write(*,*) 'Reading starspectrum.inp'
          open(unit=1,file='./input/starspectrum.inp')
          read(1,*) i
          if(i.ne.freq_nr) then
              write(*,*) 'Sorry, stellar spectrum in starspectrum.inp'
              write(*,*) 'does not have the same number of freqs as ',
     %             'other input files.'
              stop 13
          endif
          if(i.eq.1) then
              read(1,*) frq,dummy
              star_lum(1) = 1.1965280793d38 * dummy 
          else
              do i=1,freq_nr
                  read(1,*) frq,dummy 
                  if(abs(frq-freq_nu(i))/(frq+freq_nu(i)).gt.1.d-3) then
                      write(*,*) 'PROBLEM: Frequency grid of stellar ',
     %                     'spectrum unequal to frequency.inp'
                      write(*,*) frq,freq_nu(i)
                      stop 13
                  endif
                  star_lum(i) = 1.1965280793d38 * dummy 
              enddo
          endif
          close(1)
      endif
c
c      grid_nr = nr
c      grid_nz = nz
c
      end




c     --------------------------------------------------------------
c                       MAKE THE STARSPECTRUM
c     --------------------------------------------------------------
      subroutine make_star_spectrum()
      implicit none
c
#     include "common_natconst.h"
#     include "common_grid.h"
#     include "common_radialrt.h"
c
      doubleprecision bplanck
      integer inu
c
      if(freq_nr.eq.1) then
          star_lum(1) = 4 * pi * star_r**2 * ss * star_t**4
      elseif(star_t.gt.0.d0) then
          do inu=1,freq_nr
              star_lum(inu) = 4 * pi * star_r**2 * pi *
     %             bplanck(star_t,freq_nu(inu)) 
          enddo
      else
          star_t = 10.
      endif
c
      end





c     --------------------------------------------------------------
c               RADIATIVE TRANSFER ALONG RADIAL RAYS
c
c     This subroutine performs the continuum RT from one vertical
c     slice to the next. So call this routine first for ir=1, then
c     for ir=2, then ir=3 etc. This routine does all z and all
c     frequencies at once. 
c
c     ARGUMENTS:
c
c       ir          Which radius are we now?
c       repeat      = .false. --> normal
c                   = .true.  --> Redo for this ir, i.e. take
c                          all radiation stuff from the backup
c       tmid        The midplane temperature at this particular 
c                   radius, and it is meant to form a lower
c                   limit to the temperature of the dust, i.e. 
c                   to mimic the effect of the diffuse radiation 
c                   field. 
c       ifuv_evap   =-1  means: T_g=T_d
c                   =1   means: use the routines of Uma Gorti 
c     --------------------------------------------------------------
      subroutine integrate_radial_rt_onestep(ir,repeat,tmid,
     %           ifuv_evap)
      implicit none
c
      integer ir,ifuv_evap
      logical repeat
      doubleprecision tmid
c
      doubleprecision bb
c
#     include "common_disk.h"
#     include "common_dust.h"
#     include "common_grid.h"
#     include "common_radialrt.h"
#     include "common_natconst.h"
c
      doubleprecision ds,rho
      doubleprecision solve_dust_temperature
      doubleprecision solve_dust_temperature_fast
      doubleprecision t0,told,den,dustden,zz
      doubleprecision find_planckopac,mgrain
      doubleprecision meanint(1:FRSIZE_FREQ)
      doubleprecision alpha(1:FRSIZE_FREQ)
      integer iz,inu,ii
c
c     Some small calculations
c     
c      mgrain   = (4*pi/3.)*dust_sw*dust_agrain**3     
       dust_agrain=grsize(ir)
       mgrain   = (4*pi/3.)*dust_sw*grsize(ir)**3     

c
c     Make backup or restore backup
c     This is necessary for the vertical structure iteration, where
c     the radiative transfer integration has to be repeated multiple
c     times for the same ir.
c
      if(.not.repeat) then
c
c         Make backup (colvert not necessary)
c
          do iz=1,grid_nz
              do inu=1,freq_nr
                  radfbk_slum(inu,iz) = radf_slum(inu,iz)
              enddo
              radfbk_av(iz) = radf_av(iz)
              radfbk_colrad(iz) = radf_colrad(iz)
          enddo
      else
c
c         Restore backup
c
          do iz=1,grid_nz
              do inu=1,freq_nr
                  radf_slum(inu,iz) = radfbk_slum(inu,iz)
              enddo
              radf_av(iz) = radfbk_av(iz)
              radf_colrad(iz) = radfbk_colrad(iz)
          enddo
      endif
c
c     Compute the downward vertical column 
c     (necessary for Gorti & Hollenbach's gas temperature formula)
c
      radf_colvert(grid_nz) = 0.d0
      do iz=grid_nz-1,1,-1
          radf_colvert(iz) = radf_colvert(iz+1) + grid_r(ir) *
     %           abs(grid_zr(iz+1)-grid_zr(iz)) * 0.5 *
     %           (rrt_rhogas(iz+1,ir) + rrt_rhogas(iz,ir)) /
     %           ( 1.273 * mp )
c replacing below with mu for atomic gas 1.273
c     %           ( 2.3 * mp )
      enddo
c
c     Loop over z
c
      do iz=1,grid_nz
c
c         Only for ir.gt.1 do extinction;
c         For ir=1 start with the stellar luminosity
c 
          if(ir.gt.1) then
c
c             Find the opacity at the inner edge
c
              rho       = ( 0.5*rrt_rhogas(iz,ir) + 
     %                      0.5*rrt_rhogas(iz,ir-1) ) /
     %                    gastodust(ir)
              do inu=1,freq_nr
                  alpha(inu) = rho * dstkap_kappa(inu)
c    Opacity reductions is not done like below but instead when
c    dstkap_kappa is calculated through find_dust_kappa, the
c    grain size reduction is included in there.
C   Add opacity reduction
c               alpha(inu)=rho*dstkap_kappa(inu)*kappadec(ir)
              enddo
c
c             Find the ds between to successive radial shells
c
              ds = ( grid_r(ir) - grid_r(ir-1) ) *
     %             sqrt(1.d0+grid_zr(iz)**2)
c
c             Do extinction and compute also the A_V and the
c             radial column of hydrogen (last two necessary
c             for the gas temperature routines of Gorti)
c
              do inu=1,freq_nr
                  radf_slum(inu,iz) = radf_slum(inu,iz) * 
     %                     exp(-ds * alpha(inu))
              enddo
              radf_av(iz)     = radf_av(iz) + ds * rho * dstkap_kappav
C   Add opacity reduction
     %                         *kappadec(ir)
              radf_colrad(iz) = radf_colrad(iz) + ds * rho *
     %                            gastodust(ir) / (1.273 * mp )
c Replacing below with mu for atomic gas 1.273
c     %                            gastodust / (2.3 * mp )
c
          else
c
c             Copy the stellar luminosity into first bin
c
              do inu=1,freq_nr
                  radf_slum(inu,iz) = star_lum(inu)
              enddo
              radf_av(iz) = 0.d0
              radf_colrad(iz) = 0.d0
          endif
c
c         Compute the mean intensity
c
c                     F_nu         L_nu         
c            J_nu  =  ----  =   -----------   
c                     4 pi      16 pi^2 R^2 
c
          do inu=1,freq_nr
              meanint(inu) = 6.332574d-3 * radf_slum(inu,iz) / 
     %                       (grid_r(ir)**2)
          enddo
c
c         Compute the dust temperature
c
          if(freq_nr.gt.1) then
#ifndef TSOLVE_FAST
              rrt_tdust(iz,ir) = solve_dust_temperature(meanint,1,1)
#else
              rrt_tdust(iz,ir) = solve_dust_temperature_fast(
     %                                 meanint,1,1)
#endif
          else
              t0   = (meanint(1)*5.540d4)**0.25d0
              rrt_tdust(iz,ir) = t0
              do ii=1,60
                  told  = rrt_tdust(iz,ir)
c Add opacity reduction
                  rrt_tdust(iz,ir) = t0*(dstkap_kappa(1)*kappadec(ir)/
c                  rrt_tdust(iz,ir) = t0*(dstkap_kappa(1)/
     %               find_planckopac(rrt_tdust(iz,ir)))**0.25
                  if(abs(rrt_tdust(iz,ir)-told)/
     %                  (rrt_tdust(iz,ir)+told+1d-4)
     %                      .lt.TEMPACCUR) goto 201
              enddo
              write(*,*) 'ERROR: Temperature not converged'
              write(*,*) '   temp,told,ts0,kapstar,find_planckopac=',
     %             rrt_tdust(iz,ir),told,dstkap_kappa(1)*kappadec(ir),
     %               find_planckopac(rrt_tdust(iz,ir))
              stop 9323
 201          continue
          endif
c
c         Now ensure that the temperature does not drop below
c         the midplane temperature. This trick is a simple way to
c         mimic the midplane diffuse radiation field.
c         Uma: Only do this for high Av, at low Av we do not want this
c

          if(rrt_tdust(iz,ir).lt.tmid) then
             if(ir.eq.1.or.radf_av(iz).gt.3) then
              rrt_tdust(iz,ir) = tmid
             endif
          endif
c
c         Now calculate the gas temperature according to one of the
c         recipes.
c
          if(ifuv_evap.eq.-1) then
c
c             Simple use the dust temperature
c
              rrt_tgas(iz,ir) = rrt_tdust(iz,ir)
c
          elseif(ifuv_evap.eq.1) then
c
c             Use Uma Gorti & David Hollenbach's routines
c             for computing the gas temperature
c
c             Start with the initial guess
c
              rrt_tgas(iz,ir) = rrt_tdust(iz,ir)
c
c             Convert to input variables of Gorti & Hollenbach
c
              den      = rrt_rhogas(iz,ir)/(1.273*mp)
c replacing below with mu for atomic gas 
c              den      = rrt_rhogas(iz,ir)/(2.3*mp)
              dustden  = (rrt_rhogas(iz,ir)/gastodust(ir)) / mgrain
c
c             Call routine Gorti & Hollenbach
c
              zz = grid_r(ir) * grid_zr(iz)
	      if(ir.gt.1.and.radf_av(iz).le.1) then
              if(iz.gt.1) rrt_tgas(iz,ir)=rrt_tgas(iz-1,ir)
              call GasT2(grid_r(ir),zz,den,dustden,grsize(ir),
     %          rrt_tdust(iz,ir),radf_colrad(iz),
     %                radf_colvert(iz),radf_av(iz),rrt_tgas(iz,ir))
              endif
c
          else
              write(*,*) 'ERROR: Do not know ifuv_evap=',ifuv_evap
              stop
          endif
c
      enddo

c   Save midplane temperature for later use in grain evolution
        tmp(ir)=rrt_tgas(1,ir)
c
      end




c     --------------------------------------------------------------
c                THE BLACKBODY PLANCK FUNCTION B_nu(T)
c
c     This function computes the Blackbody function 
c
c                    2 h nu^3 / c^2
c        B_nu(T)  = -------------------    [ erg / cm^2 s ster Hz ]
c                   exp(-h nu / kT) - 1
c
c     ARGUMENTS:
c        nu    [Hz]            = Frequency
c        temp  [K]             = Electron temperature
c     --------------------------------------------------------------
      function bplanck(temp,nu)
      implicit none
      doubleprecision temp
      doubleprecision nu
      doubleprecision bplanck
c
      if(temp.eq.0.d0) then 
          bplanck = 0.d0
          return
      endif
c
      bplanck = 1.47455d-47 * nu * nu * nu / 
     %      (exp(4.7989d-11 * nu / temp)-1.d0) + 1.d-290
c
      return
      end


c     --------------------------------------------------------------
c                       WRITE THE TEMPERATURES
c     --------------------------------------------------------------
      subroutine write_diskstruct(idump)
      implicit none
c
      integer idump
c
#include "common_disk.h"
#include "common_grid.h"
#include "common_radialrt.h"
c
      character*80 base
      character*80 ext
      character*80 filename
      integer ir,iz
c
      base = outdir//'/diskstruct_'
      ext  = '.dat'
      call make_indexed_filename(base,idump,ext,filename)
      write(*,*) 'Writing disk structure' ! ' to file ',filename
      open(unit=1,file=filename,status='unknown')
      write(1,*) grid_nr,grid_nz
      write(1,*)
      do ir=1,grid_nr
          do iz=1,grid_nz
              write(1,10) grid_r(ir),grid_zr(iz)*grid_r(ir),
     %                   rrt_rhogas(iz,ir)+1d-90,
     %                   rrt_tgas(iz,ir)+1d-90,rrt_tdust(iz,ir)
     %                  ,grsize(ir),gastodust(ir)
 10           format(7(E13.6,1X))
          enddo
      enddo
      close(1)
      write(*,*) '    Done writing...'
c
      end


c     --------------------------------------------------------------
c                    READ THE PLANCK MEAN OPACITY
c     --------------------------------------------------------------
      subroutine read_planckopac
      implicit none
c
#     include "common_grid.h"
#     include "common_radialrt.h"
c
      doubleprecision dum1,dum2,dum3,dum4
      integer it,nt
      logical fex
c
      open(unit=1,file='./input/planckopac.inp')
      read(1,*) dopac_nt
      if(dopac_nt.gt.NR_TEMPERATURE) then
          write(*,*) 'Too many temperatures in kappaplanck.inp'
          write(*,*) ' Increase NR_TEMPERATURE and recompile'
          stop 13
      endif
      do it=1,dopac_nt
          read(1,*) dum1,dum2,dum3,dum4
          dopac_temp(it)=dum1
          dopac_kappa(it)=dum2
          dopac_kscat(it)=dum3
          dopac_gscat(it)=dum4
      enddo
      close(1)
c
      end



c     --------------------------------------------------------------
c                  INTERPOLATE THE PLANCK MEAN OPACITY
c     --------------------------------------------------------------
      function find_planckopac(temp)
      implicit none
      doubleprecision find_planckopac,temp
c
#     include "common_grid.h"
#     include "common_radialrt.h"
c
      doubleprecision dum1,dum2,dt
      integer it
c
      if(dopac_nt.le.1) then 
          find_planckopac = dopac_kappa(1)
          return
      endif
c
      if(temp.lt.dopac_temp(1)) then 
          find_planckopac = dopac_kappa(1)
      elseif(temp.gt.dopac_temp(dopac_nt)) then
          find_planckopac = dopac_kappa(dopac_nt)
      else
          call hunt(dopac_temp,dopac_nt,temp,it)
          dt = (temp-dopac_temp(it))/
     %            (dopac_temp(it+1)-dopac_temp(it))
          if((dt.lt.0.d0).or.(dt.gt.1.d0)) stop 552
          find_planckopac = (1.d0-dt)*dopac_kappa(it)+
     %       dt*dopac_kappa(it+1)
      endif
c
      return
      end


c     --------------------------------------------------------------
c                INTERPOLATE THE PLANCK MEAN SCAT OPAC
c     --------------------------------------------------------------
      subroutine find_planckall(temp,cabs,cscat,g)
      implicit none
      doubleprecision temp,cabs,cscat,g
c
#     include "common_grid.h"
#     include "common_radialrt.h"
c
      doubleprecision dum1,dum2,dt
      integer it
c
      if(dopac_nt.le.1) then 
          cabs  = dopac_kappa(1)
          cscat = dopac_kscat(1)
          g     = dopac_gscat(1)
          return
      endif
c
      if(temp.lt.dopac_temp(1)) then 
          cabs  = dopac_kappa(1)
          cscat = dopac_kscat(1)
          g     = dopac_gscat(1)
      elseif(temp.gt.dopac_temp(dopac_nt)) then
          cabs  = dopac_kappa(dopac_nt)
          cscat = dopac_kscat(dopac_nt)
          g     = dopac_gscat(dopac_nt)
      else
          call hunt(dopac_temp,dopac_nt,temp,it)
          dt = (temp-dopac_temp(it))/
     %            (dopac_temp(it+1)-dopac_temp(it))
          if((dt.lt.0.d0).or.(dt.gt.1.d0)) stop 552
          cabs  = (1.d0-dt)*dopac_kappa(it)+
     %       dt*dopac_kappa(it+1)
          cscat = (1.d0-dt)*dopac_kscat(it)+
     %       dt*dopac_kscat(it+1)
          g     = (1.d0-dt)*dopac_gscat(it)+
     %       dt*dopac_gscat(it+1)
      endif
c
      return
      end
