module ModBlock
    implicit none
    
    ! Block type definition
    type BlockType
        ! Local grid dimensions
        integer :: nx, ny, nz  ! number of grid points in each direction
        integer :: ng         ! number of ghost cells
        
        real, allocatable :: values(:,:,:,:)
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
    subroutine Init_OneBlock(block, iBlock, Mx, My, Mz, nx, ny, nz, ng, nvars, x_range_global, y_range_global, z_range_global)
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
        type(BlockType), target :: block
        integer, intent(in) :: iBlock, Mx, My, Mz
        integer, intent(in) :: nx, ny, nz, ng, nvars
        real, intent(in) :: x_range_global, y_range_global, z_range_global
        integer :: ix, iy, iz
        
        ! Store grid dimensions
        block%nx = nx
        block%ny = ny
        block%nz = nz
        block%ng = ng
        
        ! Calculate block indices (convert rank to i,j,k indices)
        ! rank = ix + Mx*(iy-1) + (Mx*My)*(iz-1)
        iz = iBlock / (Mx * My) +  1
        iy = (mod(iBlock - 1, Mx * My) + 1) / Mx + 1
        ix = mod(iBlock - 1, Mx) + 1
        
        block%x_range(1) = (ix - 1) * nx * (x_range_global / real(Mx * nx)) + 1
        block%x_range(2) = ix * nx * (x_range_global / real(Mx * nx))
        block%y_range(1) = (iy - 1) * ny * (y_range_global / real(My * ny)) + 1
        block%y_range(2) = iy * ny * (y_range_global / real(My * ny))
        block%z_range(1) = (iz - 1) * nz * (z_range_global / real(Mz * nz)) + 1
        block%z_range(2) = iz * nz * (z_range_global / real(Mz * nz))

        block%dx = (x_range_global / real(Mx * nx))
        block%dy = (y_range_global / real(My * ny))
        block%dz = (z_range_global / real(Mz * nz)) 

        block%if_bottom = (iz == 1)
        block%if_top = (iz == Mz)
        
        ! Allocate memory for values array including ghost cells
        if (allocated(block%values)) deallocate(block%values)
        allocate(block%values(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng, nvars))
        block%values = 0.0
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
        real, intent(in) :: x_range_global, y_range_global, z_range_global
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Init_OneBlock(blocks(iBlock+1), iBlock, Mx, My, Mz, nx, ny, nz, ng, nvars, &
                x_range_global, y_range_global, z_range_global)
        end do
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

    subroutine Init_Block_rho0P0(block, iBlock, Mx, My, Mz, GasGamma, z_range_global)
        implicit none
        type(BlockType), intent(inout) :: block
        integer, intent(in) :: iBlock, Mx, My, Mz
        real, intent(in) :: GasGamma, z_range_global
        integer :: i, j, k, m, iz
        real :: z_start, z_end, z_local, z_global

        iz = iBlock / (Mx * My) + 1
        z_start = block%z_range(1)
        z_end = block%z_range(2)
        m = 1.0 / (GasGamma - 1.0)

        do k = -block%ng+1, block%nz + block%ng
            do j = -block%ng+1, block%ny + block%ng
                do i = -block%ng+1, block%nx + block%ng
                    z_local = z_start + (k-1)*block%dz
                    z_global = z_local / z_range_global
                    block%rho0(i,j,k) = (1 - z_global/(m+1))**m
                    block%P0(i,j,k) = (1 - z_global/(m+1))**(m+1)
                end do
            end do
        end do
    
    end subroutine Init_Block_rho0P0

    subroutine All_Block_rho0P0(blocks, nBlocks, Mx, My, Mz, GasGamma, z_range_global)
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer, intent(in) :: nBlocks, Mx, My, Mz
        real, intent(in) :: GasGamma, z_range_global
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Init_Block_rho0P0(blocks(iBlock+1), iBlock, Mx, My, Mz, GasGamma, z_range_global)
        end do
    end subroutine All_Block_rho0P0



end module ModBlock