module ModJobAssigner
    use mpi

    ! use ModParamReader
    ! use ModParameters

    implicit none

contains
    subroutine GetRange(irank, iBlock, nGridTotal, nRanks, iStart, iEnd)
        integer, intent(in) :: irank, nGridTotal, nRanks
        integer, intent(out) :: iBlock, iStart, iEnd

        integer :: nGridPerRank

        nGridPerRank = nGridTotal/nRanks
        iStart = irank * nGridPerRank + 1
        iEnd = (irank + 1 ) * nGridPerRank

        if (irank == nRanks - 1) then
            iEnd = nGridTotal
        endif 

        iBlock = irank + 1

    end subroutine GetRange

    subroutine GetiBlock(iStart, iEnd, nGridTotal, nRanks, iBlock)
        integer, intent(in) :: iStart, iEnd, nGridTotal, nRanks
        integer, intent(out) :: iBlock

        integer :: nGridPerRank
        nGridPerRank = nGridTotal / nRanks
        iBlock = (iStart - 1) / nGridPerRank + 1
    
    end subroutine GetiBlock

end module ModJobAssigner