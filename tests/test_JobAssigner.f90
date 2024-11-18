program test_JobAssigner

    use mpi
    use ModJobAssigner

    implicit none

    integer         ::  nBlockTotal                 !   Toal N of BLocks
    integer         ::  nRanks                      !   Total N of ranks
    integer         ::  irank                       !   I of rank (starts from 0)
    integer         ::  iBlock_neighbour            !   I of neighbour block
    integer         ::  iRank_neighbour             !   Which rank the neighbour belongs to
    integer         ::  iBlock_range(2)             !   Local Block range
    integer         ::  ierr                        !   I of Error used for MPI

    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nRanks, ierr)

    nBlockTotal = 100

    call GetRange(irank, nBlockTotal, nRanks, iBlock_range)
    print *, "I'm ", irank, " and my range of Blocks is: ", iBlock_range

    iBlock_neighbour=55

    call GetiRank(nBlockTotal, nRanks, iBlock_neighbour, iRank_neighbour)
    print *, "The iBlock=", iBlock_neighbour, " belongs to ", iRank_neighbour

    call MPI_Finalize(ierr)

end program test_JobAssigner