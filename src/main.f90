program main

    use ieee_arithmetic
    use ModAdvance
    use ModCommunication
    use ModParameters,only      :   iRank,nRanks,nBlocks_local,&
        nsteps
    use ModParamReader,only     :   Read_paramters
    use ModInitiator
    use MPI

    implicit none

    type(BlockType),allocatable ::  Blocks(:)
    integer                     ::  ierr
    integer                     ::  istep
    !integer                     ::  iBlock_local

    ! MPI initiation
    call MPI_INIT(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, iRank, ierr)
    
    ! Read parameters
    call Read_paramters('/Users/qplazmfree/Documents/570/project/input/PARAM.in.test')

    ! Initiation
    allocate(Blocks(nBlocks_local))
    call Init_AllBlocks(blocks)
    call Initiate_values(Blocks)

    ! First communication
    call ModCommunication_CommunicateAll(Blocks,.false.)

    ! Run
    do istep=1,nsteps
        call ModAdvance_All(Blocks)
        if(iRank==0)print *,'End Advancing at iStep=',iStep
        print *,maxval(Blocks(1)%values)
    end do

    ! Finalize the program.
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call mpi_finalize(ierr)

end program main