module ModAdvance

    use ModSolver
    use ModParameters,only  :   nx,ny,nz,ng,rk_level,&
        iRank,nRanks
    use ModCommunication
    use ModBlock
    
    contains

    subroutine ModAdvance_All(Blocks)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        type(BlockType),pointer ::  Block1
        integer                 ::  rk_index,iBlock_local
        real                    ::  EQN_update_R(1:nx,1:ny,1:nz,5)
        logical                 ::  if_rk_input,if_rk_output

        real                    ::  dt

        ! First get dt
        call ModCommunication_GetDtGlobal(Blocks,dt)
        print *,dt

        ! Runge Kutta 4 time stepping
        do rk_index=1,rk_level
            do iBlock_local=1,size(Blocks)
                Block1=>Blocks(iBlock_local)

                if_rk_input=(rk_index>1)
                if_rk_output=(rk_index<rk_level)
                
                call GetR(Block1, EQN_update_R, if_rk_input)

                ! Get the next RK
                if (if_rk_output) then
                    Block1%values_rk(1:nx,1:ny,1:nz,:)=&
                        Block1%values(1:nx,1:ny,1:nz,:)+dt*EQN_update_R/(rk_level+1.0-rk_index)
                else
                    Block1%values(1:nx,1:ny,1:nz,:)=&
                        Block1%values(1:nx,1:ny,1:nz,:)+dt*EQN_update_R/(rk_level+1.0-rk_index)
                end if
            end do

            ! Communicate ghost cells
            call ModCommunication_CommunicateAll(Blocks,if_rk_output)
        end do

    end subroutine ModAdvance_All

end module ModAdvance