c
c     Commons for the grid
c
      common/gridrt/grid_zr,          ! Vertical grid in units of r
     %            grid_zri,           ! The cell-interfaces of the z-grid
     %            grid_dzr,           ! The cell-interfaces of the z-grid
     %            freq_nu             ! Frequency grid
      doubleprecision grid_zr(FRSIZE_Z)
      doubleprecision grid_zri(FRSIZE_Z+1)
      doubleprecision grid_dzr(FRSIZE_Z)
      doubleprecision freq_nu(FRSIZE_FREQ)
      common/igridrt/grid_nz,           ! Number of z grid points
     %             freq_nr            ! Number of frequency points
      integer grid_nz,freq_nr
c
