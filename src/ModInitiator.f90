module ModInitiator

    use ModParameters,only  :   TypeInitiation,RandomScale,nBlocks_local
    use ModBlock
    use ModJobAssigner

    contains

    ! Set the Initial Block%values for
    ! all the local blocks

    subroutine Initiate_values(Blocks)
        implicit none
        type(BlockType),intent(inout)   ::  Blocks(:)
        integer                         ::  iBlock  

        select case(TypeInitiation)
        case('random')
            do iBlock=1,nBlocks_local
                call random_number(Blocks(iBlock)%values)
                Blocks(iBlock)%values=(Blocks(iBlock)%values-0.5)*RandomScale
            end do
        end select
    end subroutine Initiate_values
end module ModInitiator