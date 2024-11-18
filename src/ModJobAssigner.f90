module ModJobAssigner
    use mpi

    ! use ModParamReader
    ! use ModParameters

    implicit none

contains
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

    subroutine GetiRank(nBlockTotal, nRanks, iBlock)
        integer, intent(in) :: nBlockTotal, nRanks
        integer, intent(out) :: iRank

        integer :: nBlockPerRank
        nBlockPerRank = nBlockTotal / nRanks
        iRank = (iBlock - 1) / nBlockPerRank
    
    end subroutine GetiRank

end module ModJobAssigner