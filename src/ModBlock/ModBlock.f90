module ModBlock
    implicit none
    
    ! Block type definition
    type BlockType
        ! Local grid dimensions
        integer :: nx, ny, nz  ! number of grid points in each direction
        integer :: ng         ! number of ghost cells
        
        real, allocatable :: values(:,:,:,:)
        
        ! Grid spacing
        real :: dx, dy, dz
        
        ! Global index ranges for this block
        integer :: x_range(2)  ! (start_index, end_index) in x direction
        integer :: y_range(2)  
        integer :: z_range(2)  
    end type BlockType

contains
    subroutine Init_OneBlock(block, irank, Mx, My, Mz, nx, ny, nz, ng, nvars)
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
        integer, intent(in) :: irank, Mx, My, Mz
        integer, intent(in) :: nx, ny, nz, ng, nvars
        integer :: ix, iy, iz
        
        ! Store grid dimensions
        block%nx = nx
        block%ny = ny
        block%nz = nz
        block%ng = ng
        
        ! Calculate block indices (convert rank to i,j,k indices)
        ! rank = ix + Mx*iy + (Mx*My)*iz
        iz = irank / (Mx * My)
        iy = mod(irank, Mx * My) / Mx
        ix = mod(irank, Mx)
        
        ! Calculate global index ranges for this block
        block%x_range(1) = ix * nx + 1
        block%x_range(2) = (ix + 1) * nx
        
        block%y_range(1) = iy * ny + 1
        block%y_range(2) = (iy + 1) * ny
        
        block%z_range(1) = iz * nz + 1
        block%z_range(2) = (iz + 1) * nz
        
        ! Allocate memory for values array including ghost cells
        if (allocated(block%values)) deallocate(block%values)
        allocate(block%values(-ng:nx+ng, -ng:ny+ng, -ng:nz+ng, nvars))
        block%values = 0.0
    end subroutine Init_OneBlock
    
    subroutine Clean_OneBlock(block)
        ! Deallocate memory for a block
        type(BlockType), intent(inout) :: block
        
        if (allocated(block%values)) deallocate(block%values)
    end subroutine Clean_OneBlock

    subroutine Init_AllBlocks(blocks, nBlocks, Mx, My, Mz, nx, ny, nz, ng, nvars)
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer, intent(in) :: nBlocks, Mx, My, Mz, nx, ny, nz, ng, nvars
        integer :: iBlock

        do iBlock = 0, nBlocks-1
            call Init_OneBlock(blocks(iBlock+1), iBlock, Mx, My, Mz, nx, ny, nz, ng, nvars)
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

end module ModBlock