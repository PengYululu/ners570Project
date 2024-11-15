program test_JobAssigner

    use mpi
    use ModJobAssigner

    implicit none

    integer :: nGridTotal, nRanks, irank, iBlock, iStart, iEnd, ierr

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)

    nGridTotal = 100

    call GetRange(irank, iBlock, nGridTotal, nRanks, iStart, iEnd)
    print *, "Ranks", irank, "for grid: ", iStart, iEnd

    call GetiBlock(iStart, iEnd, nGridTotal, nRanks, iBlock)
    print *, "Rank", irank, "for block", iBlock

    call MPI_Finalize(ierr)

end program test_JobAssigner