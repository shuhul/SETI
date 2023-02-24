c     ==============================================================
c                       PROTOSTELLAR DISK MODEL
c           MODULE: ADDITIONAL ROUTINES FOR DUST OPACITIES
c
c     ==============================================================
#include "main.h"



c     --------------------------------------------------------------
c          SOLVE TEMPERATURE OF DUST GRAIN FOR GIVEN MEAN INT
c     --------------------------------------------------------------
      function solve_dust_temperature(meanint,isize,ispec)
      implicit none
      doubleprecision solve_dust_temperature
      doubleprecision dust_integral_emisabs,zbrent
      integer isize,ispec
      doubleprecision meanint(FRSIZE_FREQ)
      external dust_integral_emisabs
c
#include "common_grid.h"
#include "common_dust.h"
c
      common/dustintegral/di_meanint
      doubleprecision di_meanint(FRSIZE_FREQ)
      common/idustintegral/di_isize,di_ispec
      integer di_isize,di_ispec
c
      doubleprecision dum1,dum2
      integer inu
c
c     Put the ispec, isize and meanint into the common
c
      di_isize   = isize
      di_ispec   = ispec
      do inu=1,freq_nr
          di_meanint(inu) = meanint(inu) 
      enddo
c
c     Now check if there is indeed sufficient radiation for a
c     solution.
c
      dum1 = dust_integral_emisabs(0.02d0)
      dum2 = dust_integral_emisabs(10000.d0)
      if(dum1.gt.0.d0) then
c          write(*,*) 'Error: dust temperature < 0.01 K'
c          write(*,*) 'The mean intensities are:'
c          write(*,*) (meanint(inu),inu=1,freq_nr)
c          stop 13
          solve_dust_temperature=0.d0
          return
      endif
      if(dum2.lt.0.d0) then
          write(*,*) 'Error: dust temperature > 10000 K'
          stop 13
      endif
c
c     Now solve for the dust temperature
c
      solve_dust_temperature = zbrent(dust_integral_emisabs,
     %                              2.d-2,1.d4,1.d-1) 
c
      return
      end



c     --------------------------------------------------------------
c                    INTEGRAL OF EMISSION-ABSORPTION
c     --------------------------------------------------------------
      function dust_integral_emisabs(temp)
      implicit none
      doubleprecision dust_integral_emisabs,temp
c
      common/dustintegral/di_meanint
      doubleprecision di_meanint(FRSIZE_FREQ)
      common/idustintegral/di_isize,di_ispec
      integer di_isize,di_ispec
c
#include "common_grid.h"
#include "common_dust.h"
c     
      integer inu
      doubleprecision dummy,kappa,bplanck,find_dust_kappa
c
      dummy = 0.d0
      do inu=1,freq_nr
          kappa = find_dust_kappa(inu,di_isize,di_ispec,temp,-1,0)
          dummy = dummy + kappa * dust_freq_wgt(inu)
     %      * ( bplanck(temp,freq_nu(inu)) - di_meanint(inu) )
      enddo
c
      dust_integral_emisabs = dummy
      return
      end





c     --------------------------------------------------------------
c          SOLVE TEMPERATURE OF DUST GRAIN FOR GIVEN MEAN INT
c     --------------------------------------------------------------
      function solve_dust_temperature_fast(meanint,isize,ispec)
      implicit none
      doubleprecision solve_dust_temperature_fast
      doubleprecision dust_integral_emisabs,zbrent
      integer isize,ispec
      doubleprecision meanint(FRSIZE_FREQ)
      external dust_integral_emisabs
c
#include "common_grid.h"
#include "common_dust.h"
c
      common/fdustintegral/di_qplus,di_temp,di_qmin
      doubleprecision di_qplus,di_temp(TSF_NTEMP)
      doubleprecision di_qmin(TSF_NTEMP,DUST_SIZE_MAX,DUST_SPECIES_MAX)
      common/ifdustintegral/di_ntemp
      integer di_ntemp
      common/idustintegral/di_isize,di_ispec
      integer di_isize,di_ispec
c
      doubleprecision dum1,dum2,temp,kappa,find_dust_kappa,eps
      integer inu,itemp
c
c     Check...
c
      if(di_qmin(1,isize,ispec).le.0.d0) then
          write(*,*) 'ERROR: Qmin array not initialized...'
          stop
      endif
c
c     Put the ispec, isize and meanint into the common
c
      di_isize   = isize
      di_ispec   = ispec
c
c     Compute Q+
c
      temp     = 10.  ! Dummy temperature
      di_qplus = 0.d0
      do inu=1,freq_nr
          kappa = find_dust_kappa(inu,di_isize,di_ispec,temp,-1,0)
          di_qplus = di_qplus + kappa * 
     %               dust_freq_wgt(inu) * meanint(inu)
      enddo
c
c     Find the index of the temperature
c
      call hunt(di_qmin(1,isize,ispec),di_ntemp,di_qplus,itemp)
      if(itemp.lt.1) then
          solve_dust_temperature_fast = di_temp(1)
      elseif(itemp.ge.di_ntemp) then
          solve_dust_temperature_fast = di_temp(di_ntemp)
      else
          eps  =  (di_qplus-di_qmin(itemp,isize,ispec)) / 
     %            (di_qmin(itemp+1,isize,ispec)-
     %             di_qmin(itemp,isize,ispec))
          if((eps.lt.0.d0).or.(eps.gt.1.d0)) stop 24751
          solve_dust_temperature_fast = 
     %         (1.d0-eps)*di_temp(itemp) + 
     %         eps*di_temp(itemp+1)
      endif
c
      return
      end



c     --------------------------------------------------------------
c               INITIALIZE THE DUST TEMPERATURE SOLUTION
c
c     NOTE: This routine should also be compatible with a future
c           implementation in RADMC. But the temperature-dependent
c           opacities are not allowed, so one must check for that.
c           This check is not yet included in this version.
c           09.07.06
c     --------------------------------------------------------------
      subroutine init_solvetempfast(ntemp,temp0,temp1)
      implicit none
      integer ntemp
      doubleprecision temp0,temp1
c
#include "common_grid.h"
#include "common_dust.h"
c
      common/fdustintegral/di_qplus,di_temp,di_qmin
      doubleprecision di_qplus,di_temp(TSF_NTEMP)
      doubleprecision di_qmin(TSF_NTEMP,DUST_SIZE_MAX,DUST_SPECIES_MAX)
      common/ifdustintegral/di_ntemp
      integer di_ntemp
      common/idustintegral/di_isize,di_ispec
      integer di_isize,di_ispec
c
      doubleprecision find_dust_kappa,tempdummy,kappa,dummy,bplanck
      integer itemp,isize,ispec,inu
c
c     Check...
c
      if(ntemp.gt.TSF_NTEMP) then
          write(*,*) 'ERROR: Nr of temperatures in the Q_min array'
          write(*,*) '  cannot exceed compiled array size. '
          write(*,*) '  Please recompile with larger TSF_NTEMP'
          stop
      endif
c
c     Now do the computation of the Q_min(T)
c
      tempdummy = 1.d0    ! Dummy; this routine does not work with T-dep kappa
      di_ntemp = ntemp
      do itemp=1,ntemp
          di_temp(itemp) = temp0 * (temp1/temp0)**
     %               ((itemp-1.d0)/(ntemp-1.d0))
          do ispec=1,dust_nr_species
              do isize=1,dust_nr_size(ispec)
                  dummy = 0.d0
                  do inu=1,freq_nr
                      kappa = find_dust_kappa(inu,isize,ispec,
     %                           tempdummy,-1,0)
                      dummy = dummy + kappa * dust_freq_wgt(inu)
     %                     * bplanck(di_temp(itemp),freq_nu(inu))
                  enddo
                  di_qmin(itemp,isize,ispec) = dummy
                  if(dummy.lt.1d-90) then
                      write(*,*) 'WARNING: lowest temperature of'
                      write(*,*) '   temp-array is lower than what'
                      write(*,*) '   can be consistent with'
                      write(*,*) '   the lowest freq-array point.'
                      write(*,*) '   TIP: chose larger temp0...'
                  endif
              enddo
          enddo
      enddo
cc###########################
c      write(*,*) temp0,temp1
c      open(unit=1,file='qmin_dummy.dat')
c      do itemp=1,ntemp
c          write(1,*) di_temp(itemp),di_qmin(itemp,1,1)
c      enddo
c      close(1)
c      stop
cc###########################
c
      end


c     --------------------------------------------------------------
c              FIND TOTAL OPACITY FROM ALL DUST SPECIES
c     --------------------------------------------------------------
      function find_total_dust_opacity(dustrho,dusttemp,inu,iabs,iscat)
      implicit none
      doubleprecision find_total_dust_opacity,find_dust_kappa
      doubleprecision dustrho(DUST_SPECIES_MAX)
      doubleprecision dusttemp(DUST_SIZE_MAX,DUST_SPECIES_MAX)
      integer inu,iabs,iscat
c
#     include "common_grid.h"
#     include "common_dust.h"
c
      doubleprecision rho,temp
      integer ispec,isize
c
      find_total_dust_opacity = 0.d0
c
      do ispec=1,dust_nr_species
          do isize=1,dust_nr_size(ispec)
              rho  = dustrho(ispec)
              temp = dusttemp(isize,ispec)
              find_total_dust_opacity = find_total_dust_opacity +
     %             rho * find_dust_kappa(inu,isize,ispec,temp,
     %             iabs,iscat)
          enddo
      enddo
c     
      return
      end



c     --------------------------------------------------------------
c              FIND TOTAL OPACITY FROM ALL DUST SPECIES
c     --------------------------------------------------------------
      function find_total_dust_emission(dustrho,dusttemp,inu,iabs,iscat)
      implicit none
      doubleprecision find_total_dust_emission,find_dust_kappa,bplanck 
      doubleprecision dustrho(DUST_SPECIES_MAX)
      doubleprecision dusttemp(DUST_SIZE_MAX,DUST_SPECIES_MAX)
      integer inu,iabs,iscat
c
#     include "common_grid.h"
#     include "common_dust.h"
c
      doubleprecision rho,temp
      integer ispec,isize
c
      find_total_dust_emission = 0.d0
c
      do ispec=1,dust_nr_species
          do isize=1,dust_nr_size(ispec)
              rho  = dustrho(ispec)
              temp = dusttemp(isize,ispec)
              find_total_dust_emission = find_total_dust_emission +
     %             rho * find_dust_kappa(inu,isize,ispec,temp,
     %             iabs,iscat) * bplanck(temp,freq_nu(inu))
          enddo
      enddo
c     
      return
      end


