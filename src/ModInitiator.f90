module ModInitiator

    use ModParameters,only  :   TypeInitiation,RandomScale
    use ModBlock
    use ModJobAssigner

    contains

    ! Set the Initial Block%values for
    ! all the local blocks

    subroutine Initiate_all(Blocks)
        implicit none
        type(BlockType),intent(inout)   ::  Blocks(:)
        integer                         ::  iBlock  

        select case(TypeInitiation)
        case('random')
            do iBlock=1,size(Blocks)
                call random_number(Blocks(iBlock)%values)
                Blocks(iBlock)%values=(Blocks(iBlock)%values-0.5)*RandomScale
            end do
        end select
    end subroutine Initiate_all
end module ModInitiator