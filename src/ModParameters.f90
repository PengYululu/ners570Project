Module ModParameters

    ! Simulation box and grid
    integer             ::  Mx,My,Mz                !   Number of blocks in each direction.
    integer             ::  nBlocks_global          !   Total number of blocks.
    integer             ::  nBlocks                 !   Local number of blocks.
    integer             ::  iBlock_range            !   Indice of local blocks.
    integer             ::  nx,ny,nz                !   Number of grid points in one block
                                                    !   for each direction.
    integer             ::  ng                      !   Number of ghost cells layers.
    real                ::  x_range_global(2),&     !   The xyz range of the whole box.
                            y_range_global(2),&
                            z_range_global(2)
    
    ! RSST
    logical             ::  UseRSST
    real                ::  Xi                      !   Xi for reduced sound speed

    ! Diffusion
    character(len=100)  ::  TypeDiffusion           !   If use gravity
    real                ::  Diffusion_h             !   Used in artificial Diffusion
    
    ! Gravity
    logical             ::  UseGravity              !   If use gravity

    ! Gas
    character(len=100)  ::  TypeGas                 !   Ideal or nonIdeal
    real                ::  GasGamma=5./3.          !   Gamma of the ideal gas
    
    ! Superadiabaticity
    character(len=100)  ::  TypeSuperAdiabaticity   !   Type of superadiabaticity
    real                ::  SuperAdiabaticity
    
end module ModParameters