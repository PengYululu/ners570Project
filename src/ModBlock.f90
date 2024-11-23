module ModBlock
    implicit none
    
    ! Block type definition
    type BlockType
        ! Local grid dimensions
        integer :: nx, ny, nz  ! number of grid points in each direction
        integer :: ng         ! number of ghost cells
        
        real, allocatable :: values(:,:,:,:)
        real, allocatable :: values_rk(:,:,:,:)
        real, allocatable :: rho0(:,:,:), P0(:,:,:)
        
        ! Grid spacing
        real :: dx, dy, dz
        
        ! Global index ranges for this block
        real :: x_range(2)  ! (start_index, end_index) in x direction
        real :: y_range(2)  
        real :: z_range(2)  

        logical :: if_top
        logical :: if_bottom
        
    end type BlockType

contains
    subroutine Init_OneBlock(block1, iBlock, Mx, My, Mz, nx, ny, nz, ng, nvars, x_range_global, y_range_global, z_range_global)
        ! Initialize a block based on MPI rank and block configuration
        !
        ! Arguments:
        ! block  - Block to initialize
        ! irank  - MPI rank number (0-based)
        ! Mx,My,Mz - Number of blocks in each direction
        ! nx,ny,nz - Number of grid points per block in each direction
        ! ng      - Number of ghost cells
        ! nvars   - Number of variables stored per grid point
        implicit none
        type(BlockType), target :: block1
        integer, intent(in) :: iBlock, Mx, My, Mz
        integer, intent(in) :: nx, ny, nz, ng, nvars
        real, intent(in) :: x_range_global(2), y_range_global(2), z_range_global(2)
        integer :: ix, iy, iz
        
        ! Store grid dimensions
        block1%nx = nx
        block1%ny = ny
        block1%nz = nz
        block1%ng = ng
        
        ! Calculate block indices (convert rank to i,j,k indices)
        ! rank = ix + Mx*(iy-1) + (Mx*My)*(iz-1)
        iz = iBlock / (Mx * My) +  1
        iy = (mod(iBlock - 1, Mx * My) + 1) / Mx + 1
        ix = mod(iBlock - 1, Mx) + 1
        
        block1%x_range(1) = x_range_global(1) * (Mx-ix+1)   + x_range_global(2) * (ix-1)
        block1%x_range(2) = x_range_global(1) * (Mx-ix)     + x_range_global(2) * ix
        block1%y_range(1) = y_range_global(1) * (My-iy+1)   + y_range_global(2) * (iy-1)
        block1%y_range(2) = y_range_global(1) * (My-iy)     + y_range_global(2) * iy
        block1%z_range(1) = z_range_global(1) * (Mz-iz+1)   + z_range_global(2) * (iz-1)
        block1%z_range(2) = z_range_global(1) * (Mz-iz)     + z_range_global(2) * iz

        block1%dx = (block1%x_range(2) - block1%x_range(1))/nx
        block1%dy = (block1%y_range(2) - block1%y_range(1))/ny
        block1%dz = (block1%z_range(2) - block1%z_range(1))/nz

        block1%if_bottom = (iz == 1)
        block1%if_top = (iz == Mz)
        
        ! Allocate memory for values array including ghost cells
        if (allocated(block1%values)) deallocate(block1%values)
        allocate(block1%values(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng, nvars))

        ! Allocate memory for values array including ghost cells
        if (allocated(block1%values_rk)) deallocate(block1%values_rk)
        allocate(block1%values_rk(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng, nvars))

        ! Set initial values to 0
        block1%values = 0.0
        block1%values_rk = 0.0
    end subroutine Init_OneBlock
    
    subroutine Clean_OneBlock(block)
        ! Deallocate memory for a block
        type(BlockType), intent(inout) :: block
        
        if (allocated(block%values)) deallocate(block%values)
    end subroutine Clean_OneBlock

    subroutine Init_AllBlocks(blocks, nBlocks, Mx, My, Mz, nx, ny, nz, ng, nvars, &
        x_range_global, y_range_global, z_range_global)
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer, intent(in) :: nBlocks, Mx, My, Mz, nx, ny, nz, ng, nvars
        real, intent(in) :: x_range_global(2), y_range_global(2), z_range_global(2)
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Init_OneBlock(blocks(iBlock+1), iBlock, Mx, My, Mz, nx, ny, nz, ng, nvars, &
                x_range_global, y_range_global, z_range_global)
        end do
        call All_Block_rho0P0(blocks, nBlocks)
    end subroutine Init_AllBlocks

    subroutine Clean_AllBlocks(blocks, nBlocks)
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer, intent(in) :: nBlocks
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Clean_OneBlock(blocks(iBlock+1))
        end do
    end subroutine Clean_AllBlocks

    subroutine Init_Block_rho0P0(block1, GasGamma)
        implicit none
        type(BlockType), intent(inout) :: block1
        real, intent(in) :: GasGamma
        integer ::  i, j, k
        real    ::  m,z

        m = 1.0 / (GasGamma - 1.0)

        do k = -block1%ng+1, block1%nz + block1%ng
            do j = -block1%ng+1, block1%ny + block1%ng
                do i = -block1%ng+1, block1%nx + block1%ng
                    z = block1%z_range(1) + (k-0.5)*block1%dz
                    block1%rho0(i,j,k) = (1 - z/(m+1))**m
                    block1%P0(i,j,k) = (1 - z/(m+1))**(m+1)
                end do
            end do
        end do
    
    end subroutine Init_Block_rho0P0

    subroutine All_Block_rho0P0(blocks, nBlocks)
        use ModParameters,only  :   GasGamma
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer, intent(in) :: nBlocks
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Init_Block_rho0P0(blocks(iBlock+1), GasGamma)
        end do
    end subroutine All_Block_rho0P0



end module ModBlock