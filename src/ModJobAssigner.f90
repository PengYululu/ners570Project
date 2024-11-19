module ModJobAssigner

    implicit none

contains

    ! Get the range of local iBlocks

    subroutine GetRange(irank, nBlockTotal, nRanks, iBlock_range)
        integer, intent(in)     ::  irank               !   I of rank
        integer, intent(in)     ::  nBlockTotal         !   Total N of blocks
        integer, intent(in)     ::  nRanks              !   Total N of ranks
        integer, intent(out)    ::  iBlock_range(2)     !   Local range of iBlocks

        integer :: nBlockPerRank

        nBlockPerRank = nBlockTotal/nRanks
        iBlock_range=[irank * nBlockPerRank + 1,(irank + 1 ) * nBlockPerRank]
    end subroutine GetRange

    ! Use iBlock to find which iRank it belongs to

    subroutine GetiRank(nBlockTotal, nRanks, iBlock, iRank)
        integer, intent(in)     ::  nBlockTotal     !   N of Blocks total
        integer, intent(in)     ::  nRanks          !   N of ranks total
        integer, intent(in)     ::  iBlock          !   i of Block to look at
        integer, intent(out)    ::  iRank           !   The result i of rank

        integer                 ::  nBlockPerRank   !   How many blocks per rank

        nBlockPerRank = nBlockTotal / nRanks
        iRank = (iBlock - 1) / nBlockPerRank
    end subroutine GetiRank

end module ModJobAssigner