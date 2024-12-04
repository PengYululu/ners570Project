module ModBlock

    use ModJobAssigner
    use ModParameters, only : nx,ny,nz,ng,Mx,My,Mz,&
        x_range_global,y_range_global,z_range_global,&
        nBlocks_local,iRank,nRanks,iBlock_range
    implicit none
    
    ! Block type definition
    type BlockType
        integer             ::  iBlock                  ! Global block index
        integer             ::  nx, ny, nz              ! number of grid points in each direction
        integer             ::  ng                      ! number of ghost cells
        
        real, allocatable   ::  values(:,:,:,:)
        real, allocatable   ::  values_rk(:,:,:,:)
        real, allocatable   ::  rho0(:,:,:),P0(:,:,:)   ! Background rho and p
        
        real                ::  dx, dy, dz              ! Grid spacing
        real                ::  x_range(2),&            ! Global index ranges for this block
                                y_range(2),&            ! (start_index, end_index) in x direction
                                z_range(2)
        
        real, allocatable   ::  x_list(:),y_list(:),&
                                z_list(:)

        logical             ::  if_top,if_bottom        ! If near top or bottom boundary
        integer             ::  iBlock_l,iBlock_r,&     ! Neighbours
                                iBlock_f,iBlock_b,&
                                iBlock_u,iBlock_d
        integer             ::  iRank_l,iRank_r,&
                                iRank_f,iRank_b,&
                                iRank_u,iRank_d

        integer             ::  n_sends                 ! Number of information to be sent each time

        real, allocatable   ::  Buffer_send_l(:,:,:,:),Buffer_send_r(:,:,:,:),&
                                Buffer_send_f(:,:,:,:),Buffer_send_b(:,:,:,:),&
                                Buffer_send_u(:,:,:,:),Buffer_send_d(:,:,:,:)
        
        real, allocatable   ::  Buffer_recv_l(:,:,:,:),Buffer_recv_r(:,:,:,:),&
                                Buffer_recv_f(:,:,:,:),Buffer_recv_b(:,:,:,:),&
                                Buffer_recv_u(:,:,:,:),Buffer_recv_d(:,:,:,:)
    end type BlockType

    contains

    subroutine Init_AllBlocks(blocks)
        implicit none
        type(BlockType)         ::  blocks(:)
        integer                 ::  iBlock_local,iBlock_global

        do iBlock_local = 1, nBlocks_local
            
            iBlock_global=iBlock_local+iBlock_range(1)-1
            call Init_OneBlock(blocks(iBlock_local),iBlock_global)
        end do
        call All_Block_rho0P0(blocks)
    end subroutine Init_AllBlocks

    subroutine Init_OneBlock(block1, iBlock)
        ! Initialize a block based on MPI rank and block configuration
        !
        ! Arguments:
        ! block  - Block to initialize
        ! irank  - MPI rank number (0-based)
        ! nvars   - Number of variables stored per grid point
        implicit none
        type(BlockType), target ::  block1
        integer, intent(in)     ::  iBlock
        integer                 ::  ix, iy, iz,i
        
        ! Label the block index
        Block1%iBlock=iBlock
        ! Store grid dimensions
        block1%nx = nx
        block1%ny = ny
        block1%nz = nz
        block1%ng = ng
        
        ! Calculate block indices (convert rank to i,j,k indices)
        ! rank = ix + Mx*(iy-1) + (Mx*My)*(iz-1)
        iz = (iBlock - 1) / (Mx * My) +  1
        iy = mod(iBlock - 1, Mx * My) / Mx + 1
        ix = mod(iBlock - 1, Mx) + 1
        
        block1%x_range(1) = (x_range_global(1) * (Mx-ix+1)   + x_range_global(2) * (ix-1)   )/Mx
        block1%x_range(2) = (x_range_global(1) * (Mx-ix)     + x_range_global(2) * ix       )/Mx
        block1%y_range(1) = (y_range_global(1) * (My-iy+1)   + y_range_global(2) * (iy-1)   )/My
        block1%y_range(2) = (y_range_global(1) * (My-iy)     + y_range_global(2) * iy       )/My
        block1%z_range(1) = (z_range_global(1) * (Mz-iz+1)   + z_range_global(2) * (iz-1)   )/Mz
        block1%z_range(2) = (z_range_global(1) * (Mz-iz)     + z_range_global(2) * iz       )/Mz

        ! Get the grid size
        block1%dx = (block1%x_range(2) - block1%x_range(1))/nx
        block1%dy = (block1%y_range(2) - block1%y_range(1))/ny
        block1%dz = (block1%z_range(2) - block1%z_range(1))/nz

        ! Get xyz_list

        allocate(block1%x_list(-ng+1:nx+ng))
        allocate(block1%y_list(-ng+1:ny+ng))
        allocate(block1%z_list(-ng+1:nz+ng))

        do i=-ng+1,nx+ng
            block1%x_list(i)=(block1%x_range(1)*(nx-i+0.5)+block1%x_range(2)*(i-0.5))/nx
        end do
        do i=-ng+1,ny+ng
            block1%y_list(i)=(block1%y_range(1)*(ny-i+0.5)+block1%y_range(2)*(i-0.5))/ny
        end do
        do i=-ng+1,nz+ng
            block1%z_list(i)=(block1%z_range(1)*(nz-i+0.5)+block1%z_range(2)*(i-0.5))/nz
        end do

        ! Determine if this block is near top or bottom boundary 
        block1%if_bottom = (iz == 1)
        block1%if_top = (iz == Mz)
        
        ! Allocate values array including ghost cells
        allocate(block1%values(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng, 5))
        allocate(block1%values_rk(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng, 5))

        ! Set initial values to 0
        block1%values = 0.0
        block1%values_rk = 0.0

        ! Allocate rho0 and p0
        allocate(Block1%p0(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng))
        allocate(Block1%rho0(-ng+1:nx+ng, -ng+1:ny+ng, -ng+1:nz+ng))

        ! Get the neighbours
        Block1%iBlock_l=iBlock-1;       Block1%iBlock_r=iBlock+1
        Block1%iBLock_b=iBlock-Mx;      Block1%iBlock_f=iBlock+Mx
        Block1%iBlock_u=iBlock+Mx*My;   Block1%iBlock_d=iBLock-Mx*My
        if (ix==1)  Block1%iBlock_l=Block1%iBlock_l+Mx
        if (ix==Mx) Block1%iBlock_r=Block1%iBlock_r-Mx
        if (iy==1)  Block1%iBlock_b=Block1%iBlock_b+Mx*My
        if (iy==My) Block1%iBlock_f=Block1%iBlock_f-Mx*My
        if (iz==1)  Block1%iBlock_d=-1
        if (iz==Mz) Block1%iBlock_u=-1

        ! Get the ranks containing these neighbours
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_l, Block1%iRank_l)
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_r, Block1%iRank_r)
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_b, Block1%iRank_b)
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_f, Block1%iRank_f)
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_u, Block1%iRank_u)
        call GetiRank(Mx*My*Mz, nRanks, Block1%iBlock_d, Block1%iRank_d)

        Block1%n_sends=0
        if (Block1%iRank_l/=iRank) Block1%n_sends=Block1%n_sends+1
        if (Block1%iRank_r/=iRank) Block1%n_sends=Block1%n_sends+1
        if (Block1%iRank_b/=iRank) Block1%n_sends=Block1%n_sends+1
        if (Block1%iRank_f/=iRank) Block1%n_sends=Block1%n_sends+1
        if (Block1%iRank_u/=iRank .and. .not. Block1%if_top) Block1%n_sends=Block1%n_sends+1
        if (Block1%iRank_d/=iRank .and. .not. Block1%if_bottom) Block1%n_sends=Block1%n_sends+1

        ! Allocate the buffer grids
        if (Block1%iRank_l/=iRank) &
        allocate(Block1%Buffer_send_l(1:ng,1:ny,1:nz,5),        Block1%Buffer_recv_l(1:ng,1:ny,1:nz,5)   )
        if (Block1%iRank_r/=iRank) &
        allocate(Block1%Buffer_send_r(nx-ng+1:nx,1:ny,1:nz,5),  Block1%Buffer_recv_r(nx+1:nx+ng,1:ny,1:nz,5))
        if (Block1%iRank_b/=iRank) &
        allocate(Block1%Buffer_send_b(1:nx,1:ng,1:nz,5),        Block1%Buffer_recv_b(1:nx,-ng+1:0,1:nz,5)   )
        if (Block1%iRank_f/=iRank) &
        allocate(Block1%Buffer_send_f(1:nx,ny-ng+1:ny,1:nz,5),  Block1%Buffer_recv_f(1:nx,ny+1:ny+ng,1:nz,5))
        if (Block1%iRank_d/=iRank .and. .not. Block1%if_bottom) &
        allocate(Block1%Buffer_send_d(1:nx,1:ny,1:ng,5),        Block1%Buffer_recv_d(1:nx,1:ny,-ng+1:0,5)   )
        if (Block1%iRank_u/=iRank .and. .not. Block1%if_top) &
        allocate(Block1%Buffer_send_u(1:nx,1:ny,nz-ng+1:nz,5),  Block1%Buffer_recv_u(1:nx,1:ny,nz+1:nz+ng,5))
    end subroutine Init_OneBlock
    
    subroutine Clean_OneBlock(block)
        ! Deallocate memory for a block
        type(BlockType), intent(inout) :: block
        
        if (allocated(block%values)) deallocate(block%values)
        if (allocated(block%values_rk)) deallocate(block%values_rk)
        if (allocated(block%rho0)) deallocate(block%rho0)
        if (allocated(block%p0)) deallocate(block%p0)
    end subroutine Clean_OneBlock

    subroutine Clean_AllBlocks(blocks)
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer :: iBlock

        do iBlock = 0, nBlocks_local-1
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

        do k = -ng+1, nz + ng
            do j = -ng+1, ny + ng
                do i = -ng+1, nx + ng
                    z = block1%z_list(k)
                    block1%rho0(i,j,k) = (1.0 - z/(m+1))**m
                    block1%P0(i,j,k) = (1.0 - z/(m+1))**(m+1)
                end do
            end do
        end do
    
    end subroutine Init_Block_rho0P0

    subroutine All_Block_rho0P0(blocks)
        use ModParameters,only  :   GasGamma
        implicit none
        type(BlockType), dimension(:), intent(inout) :: blocks
        integer :: iBlock_local

        do iBlock_local = 1, nBlocks_local
            call Init_Block_rho0P0(blocks(iBlock_local), GasGamma)
        end do
    end subroutine All_Block_rho0P0



end module ModBlock