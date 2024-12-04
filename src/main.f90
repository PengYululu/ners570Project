program main

    use ieee_arithmetic
    use ModAdvance
    use ModCommunication
    use ModParameters,only      :   iRank,nRanks,nBlocks_local,&
        nsteps,nx,ny,nz
    use ModParamReader,only     :   Read_paramters
    use ModInitiator
    use ModSavePlot
    use MPI

    implicit none

    type(BlockType),allocatable ::  Blocks(:)
    integer                     ::  ierr
    integer                     ::  istep
    integer                     ::  i,j,k
    character(len=8)            ::  iStep_char
    !integer                     ::  iBlock_local

    ! MPI initiation
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, iRank, ierr)
    
    ! Read parameters
    call Read_paramters('./PARAM.in.test')

    ! Initiation
    allocate(Blocks(nBlocks_local))
    call Init_AllBlocks(blocks)
    call Initiate_values(Blocks)

    ! First communication
    !do k=1,nz; do j=1,ny; do i=1,nx
    !    Blocks(1)%values(i,j,k,:)=Blocks(1)%x_list(i)
    !end do; end do; end do
    
    call ModCommunication_CommunicateAll(Blocks,.false.)
    call ModSave_2D(Blocks,'k',2.1,[100,100],filename='test.dat',logical_unit=2)
    !print *,iRank,Blocks(1)%x_list(0),Blocks(1)%values(0,1,1,:)
    !print *,"I'm",iRank,Blocks(1)%iBlock_r,Blocks(1)%iBlock_l,&
    !    Blocks(1)%iBlock_b,Blocks(1)%iBlock_f,&
    !    Blocks(1)%iBlock_u,Blocks(1)%iBlock_d
    ! Run
    
    do istep=1,nsteps
        call ModAdvance_All(Blocks)
        if(iRank==0)print *,'End Advancing at iStep=',iStep
        if(iRank==0)print *,maxval(Blocks(1)%values)

        write(iStep_char,'(I8)')iStep
        if (mod(iStep,500)==0) then
                
            do i=1,len_trim(iStep_char)
                if (iStep_char(i:i)==' ') iStep_char(i:i)='0'
            end do
            call ModSave_2D(Blocks,'k',2.1,[100,100],filename=iStep_char//'.dat',logical_unit=2)
        end if
        
    end do

    ! Finalize the program.
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_finalize(ierr)

end program main