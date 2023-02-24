c
c     Commons for the central star
c
      common/starrt/star_t,             ! Temperature of central star
     %            star_r,             ! Radius of central star
     %            star_m,             ! Mass of star
     %            star_lum            ! Spectrum of the star in lum
      doubleprecision star_r,star_t,star_m
      doubleprecision star_lum(FRSIZE_FREQ)
c
c     Commons for the disk
c
      common/diskrt/rrt_rhogas,rrt_tdust,rrt_tgas,rrt_vz
cc,rrt_hp,rrt_sigmadust
      doubleprecision rrt_rhogas(FRSIZE_Z,FRSIZE_R)
      doubleprecision rrt_tdust(FRSIZE_Z,FRSIZE_R)
      doubleprecision rrt_tgas(FRSIZE_Z,FRSIZE_R)
      doubleprecision rrt_vz(FRSIZE_Z,FRSIZE_R)
cc      doubleprecision rrt_hp(FRSIZE_R)
cc      doubleprecision rrt_sigmadust(FRSIZE_R)
c
c     The planck-mean dust opacity
c
      common/dopacrt/dopac_temp,dopac_kappa,dopac_plkap_star,
     %        dopac_kscat,dopac_gscat
      doubleprecision dopac_temp(NR_TEMPERATURE),dopac_plkap_star
      doubleprecision dopac_kappa(NR_TEMPERATURE)
      doubleprecision dopac_kscat(NR_TEMPERATURE)
      doubleprecision dopac_gscat(NR_TEMPERATURE)
      common/idopacrt/dopac_nt
      integer dopac_nt
c
      common/radfield/radf_slum,radf_av,radf_colrad,radf_colvert
      doubleprecision radf_slum(FRSIZE_FREQ,FRSIZE_Z)
      doubleprecision radf_av(FRSIZE_Z)
      doubleprecision radf_colrad(FRSIZE_Z)
      doubleprecision radf_colvert(FRSIZE_Z)
c
      common/radfieldbk/radfbk_slum,radfbk_av,radfbk_colrad,
     %                  radfbk_colvert
      doubleprecision radfbk_slum(FRSIZE_FREQ,FRSIZE_Z)
      doubleprecision radfbk_av(FRSIZE_Z)
      doubleprecision radfbk_colrad(FRSIZE_Z)
      doubleprecision radfbk_colvert(FRSIZE_Z)
c
      common/dstkap/dstkap_kappa,dstkap_kappav
      doubleprecision dstkap_kappa(FRSIZE_FREQ),dstkap_kappav
c
      common/ivertstruct/vs_iter_nr
      integer vs_iter_nr
      common/vertstruct/vs_iter_crit
      doubleprecision vs_iter_crit
c
