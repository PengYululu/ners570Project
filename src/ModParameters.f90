Module ModParameters

    ! MPI information
    integer             ::  iRank,nRanks

    ! Simulation box and grid
    integer             ::  Mx,My,Mz                !   Number of blocks in each direction.
    integer             ::  nBlocks_global          !   Total number of blocks.
    integer             ::  nBlocks_local           !   Local number of blocks.
    integer             ::  iBlock_range(2)         !   Indice of local blocks.
    integer             ::  nx,ny,nz                !   Number of grid points in one block
                                                    !   for each direction.
    integer             ::  ng                      !   Number of ghost cells layers.
    real                ::  x_range_global(2),&     !   The xyz range of the whole box.
                            y_range_global(2),&
                            z_range_global(2)
    
    ! RSST
    logical             ::  UseRSST=.True.
    real                ::  Xi=1.0                  !   Xi for reduced sound speed

    ! Diffusion
    character(len=100)  ::  TypeDiffusion           !   If use gravity
    real                ::  Diffusion_h             !   Used in artificial Diffusion
    
    ! Gravity
    logical             ::  UseGravity=.True.       !   If use gravity

    ! Gas
    character(len=100)  ::  TypeGas                 !   Ideal or nonIdeal
    real                ::  GasGamma=5./3.          !   Gamma of the ideal gas
    
    ! Superadiabaticity
    character(len=100)  ::  TypeSuperAdiabaticity   !   Type of superadiabaticity
    real                ::  SuperAdiabaticity=1.e-4
    
    ! Initiatino
    character(len=100)  ::  TypeInitiation
    real                ::  RandomScale=1.e-7

    ! TimeStepping
    integer             ::  rk_level=4
    real                ::  CFL=0.5
    
    ! RUN
    integer             ::  nsteps=100

    ! Save plot
    integer             ::  nStepSavePlot
end module ModParameters