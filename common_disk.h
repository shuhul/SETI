      common/igrid/grid_nr
      integer grid_nr
c
      common/grid/grid_r,grid_dr
      doubleprecision grid_r(FRSIZE_R)
      doubleprecision grid_dr(FRSIZE_R)
c
      common/disk/disk_tprsigma,disk_vr,disk_nu,
     %            disk_temp,disk_v0,disk_iv0,disk_dcoef,
     %            dust_dcoef,snowrad,disk_tprsigmadtot,
     %            disk_tprsigmadust,disk_vrdust,disk_mdot,
     %            disk_mgrav,disk_teff,gastodust,
     %            disk_tprsigma_prev,planetesimal,
     %            disk_tprsigmaice,disk_tprsigmavapor
      doubleprecision disk_tprsigma(FRSIZE_R)
      doubleprecision snowrad
      doubleprecision disk_tprsigmadust(NDUSTBIN,FRSIZE_R)
      doubleprecision disk_tprsigmadtot(FRSIZE_R)
      doubleprecision disk_vr(FRSIZE_R+1)
      doubleprecision disk_vrdust(FRSIZE_R+1,NRDUST)
      doubleprecision disk_temp(FRSIZE_R)
      doubleprecision disk_nu(FRSIZE_R)
      doubleprecision disk_v0(FRSIZE_R)
      doubleprecision disk_iv0(FRSIZE_R)
      doubleprecision dust_vr(FRSIZE_R)
      doubleprecision disk_dcoef(FRSIZE_R)
      doubleprecision dust_dcoef(FRSIZE_R)
      doubleprecision disk_mdot(FRSIZE_R+1)
      doubleprecision disk_mgrav(FRSIZE_R)
      doubleprecision disk_teff(FRSIZE_R)
      doubleprecision gastodust(FRSIZE_R)
      doubleprecision planetesimal(FRSIZE_R)
      doubleprecision disk_tprsigmaice(FRSIZE_R,NDUSTBIN,2)
      doubleprecision disk_tprsigmavapor(FRSIZE_R,2)
      doubleprecision disk_tprsigma_prev(FRSIZE_R)
c
      common/dustevap/dust_tevap
      doubleprecision dust_tevap
      common/idustevap/idust_evap
      integer idust_evap
c
      common/tmpmin/tempbg
      doubleprecision tempbg
c
      common/time/time_dt,time_dtrel,time_dtmin,
     %            time_time,time_end,time_tsave,time_inctsave,
     %            time_savref_mdot,time_savref_dt,time_stop_mdot,
     %            time_savref_tmin,time_stop_mdisk
      doubleprecision time_dt,time_dtrel,time_dtmin
      doubleprecision time_time,time_end,time_tsave,time_inctsave
      doubleprecision time_savref_mdot,time_savref_dt
      doubleprecision time_stop_mdot,time_savref_tmin
      doubleprecision time_stop_mdisk
c     
      common/itime/time_maxstep,time_isave,time_nr_save
      integer time_maxstep,time_isave,time_nr_save
c
      common/global/global_mstar,global_alpha,global_prandtlinv,
     %              global_sigma_min,global_gtd,global_adust,
     %              global_sigmad_min,
     %              global_rstar,global_tstar,global_flang,
     %              global_alpnew,global_mugas,global_stlumgrow,
     %              global_j_in,global_j_out,global_ginst_plaw
      doubleprecision global_mstar,global_mugas,global_prandtlinv
      doubleprecision global_alpha(FRSIZE_R),global_gtd
      doubleprecision global_ginst_plaw
      doubleprecision global_sigma_min,global_adust
      doubleprecision global_sigmad_min
      doubleprecision global_j_in,global_j_out
      doubleprecision global_rstar,global_tstar,global_flang
      doubleprecision global_alpnew(FRSIZE_R)
      doubleprecision global_stlumgrow
c
      common/iglobal/global_ndust,global_flag_lowtemp,global_gravinst,
     %               global_visheat_mode
      integer global_ndust,global_flag_lowtemp,global_gravinst,
     %               global_visheat_mode
c
      common/rossopac/ross_temp,ross_kappa
      doubleprecision ross_temp(FRSIZE_TEMP)
      doubleprecision ross_kappa(FRSIZE_TEMP,2)
      common/irossopac/ross_ntemp
      integer ross_ntemp
c
      common/vdt/varsdt
      logical varsdt
c
      common/slvtemp/st_fact_qpl,st_fact_toomre,st_sigmakap,
     %               st_alpha_grav,st_alpha_mri,st_alpha,
     %               st_qplus,st_tempeff
      doubleprecision st_fact_qpl,st_fact_toomre,st_sigmakap,
     %               st_alpha_grav,st_alpha_mri,st_alpha,
     %               st_qplus,st_tempeff
      common/islvtemp/ist_icmp
      integer ist_icmp
c     
      common/dpr/dpr_tempcryst,dpr_nu
      doubleprecision dpr_tempcryst,dpr_nu
      common/idpr/dpr_mode
      integer dpr_mode
c
      common/switch/isw_qvisc,isw_qirrad,isw_growstar,isw_incmdisk,
     %              isw_lumbl,isw_selfirr,isw_euv_evap,isw_fuv_evap
     %        ,isw_dustevolve,isw_backreacn
      integer isw_qvisc,isw_qirrad,isw_growstar,isw_incmdisk
      integer isw_lumbl,isw_selfirr,isw_euv_evap,isw_fuv_evap
      integer isw_dustevolve,isw_backreacn
c
      common/cloud/cloud_mass,cloud_asound,cloud_omega,
     %             cloud_smoothparam1,cloud_smoothparam2
      doubleprecision cloud_mass,cloud_asound,cloud_omega,
     %             cloud_smoothparam1,cloud_smoothparam2
c
      common/icloud/cloud_ismooth,cloud_idistr
      integer cloud_ismooth,cloud_idistr
c
      common/massload/ml_psi,ml_mdot,ml_mdotcap,ml_rcentr
      doubleprecision ml_psi(FRSIZE_R)
      doubleprecision ml_mdot,ml_mdotcap,ml_rcentr
      common/imassload/ml_implicit
      integer ml_implicit
c
      common/photevap/euvp_phi,euvp_nbase,photoev_sigdot,
     % photoev_isigdot,pe_mdot
      doubleprecision euvp_phi,pe_mdot,euvp_nbase(FRSIZE_R)
      doubleprecision photoev_sigdot(FRSIZE_R)
      doubleprecision photoev_isigdot(FRSIZE_R)
      common/iphotoevap/photoev_izev
      integer photoev_izev(FRSIZE_R)
c
      common/debug/debug_yes
      integer debug_yes
c
      common/ifuvevapdata/fuvevap_nz,idump_rt,idump_count,idump_isave
      integer fuvevap_nz,idump_rt,idump_count,idump_isave
      common/fuvevapdata/fuvevap_zrmin,fuvevap_zrmax
      doubleprecision fuvevap_zrmin,fuvevap_zrmax
c
      common/ideadzone/isw_deadzone
      integer isw_deadzone
      common/deadzone/deadzone_siglay,deadzone_tempmri,deadzone_alpha
      doubleprecision deadzone_siglay,deadzone_tempmri,deadzone_alpha
c
c      Added by Uma Gorti (09/08) to include gast.F, this header also included
c	   in the subroutine Thermbal, G0 is conputed in diskevol.F
c	   common /starhe/ xLx, fuvlum
	   doubleprecision xLx, fuvlum

c	An output directory to be specified each run
       character*6 outdir
       parameter(outdir='models')
c      Turbulent alpha
       doubleprecision alpha_turb 
       parameter(alpha_turb=1.0d-4)

!      array to hold mean dust size and opacity reduction for dust evolution
       doubleprecision taufac
       doubleprecision grsize(FRSIZE_R), kappadec(FRSIZE_R)
       doubleprecision tmp(FRSIZE_R),fdustpe(FRSIZE_R),pahscale(FRSIZE_R)
       common /devolve/ taufac,grsize, kappadec,tmp,fdustpe,pahscale


       common /particlemass/ mup
       doubleprecision mup(FRSIZE_Z,FRSIZE_R)

       common /ginst/ igi
       integer igi(FRSIZE_R)

       doubleprecision decel(FRSIZE_R),gasvel(FRSIZE_R),dustvel(FRSIZE_R)
       common /backreac/ decel,gasvel,dustvel




