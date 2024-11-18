module ModJobAssigner

    implicit none

contains

    ! Get the range of local iBlocks
    subroutine GetRange(irank, iBlock, nBlockTotal, nRanks, iStart, iEnd)
        integer, intent(in) :: irank, nBlockTotal, nRanks
        integer, intent(out) :: iBlock, iStart, iEnd

        integer :: nBlockPerRank

        nBlockPerRank = nBlockTotal/nRanks
        iStart = irank * nBlockPerRank + 1
        iEnd = (irank + 1 ) * nBlockPerRank

        if (irank == nRanks - 1) then
            iEnd = nBlockTotal
        endif 

        iBlock = irank + 1

    end subroutine GetRange

    ! Use iBlock to find which iRank
    ! it belongs to

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