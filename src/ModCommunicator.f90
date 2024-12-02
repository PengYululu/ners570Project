module ModCommunication

    use ModParameters, only :    nx,ny,nz,ng,CFL,Xi,GasGamma,&
        SuperAdiabaticity,iRank,nRanks,iBlock_range
    use ModBlock
    use mpi
    
    contains

    subroutine ModCommunication_GetDtGlobal(Blocks,dt_global)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        real,intent(out)        ::  dt_global
        real                    ::  dt_local
        integer                 ::  ierr

        call ModCommunication_GetDtLocal(Blocks,dt_local)
        call MPI_ALLReduce(dt_local,dt_global,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
    end subroutine ModCommunication_GetDtGlobal

    subroutine ModCommunication_GetDtLocal(Blocks,dt)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        type(BlockType),pointer ::  Block1
        real,intent(out)        ::  dt
        real                    ::  dt1
        integer                 ::  iBlock_local

        dt=1.e20
        do iBlock_local=1,size(Blocks)
            Block1=>Blocks(iBlock_local)
            dt1 = CFL*min(Block1%dx,Block1%dy,Block1%dz)/&
                (maxval(abs(Block1%values(:,:,:,2:4)))+&
                1./Xi*sqrt(GasGamma/(8.*SuperAdiabaticity))*maxval(Block1%p0/Block1%rho0))
            dt=min(dt,dt1)
        end do
    end subroutine ModCommunication_GetDtLocal

    subroutine ModCommunication_CommunicateAll(Blocks,if_rk)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        logical,intent(in)      ::  if_rk

        call ModCommunication_Local(Blocks,if_rk)
        call ModCommunication_SendRecv(Blocks,if_rk)
    end subroutine ModCommunication_CommunicateAll

    subroutine ModCommunication_Local(Blocks,if_rk)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        type(BlockType),pointer ::  Block1,Block2
        logical,intent(in)      ::  if_rk
        real,pointer            ::  values(:,:,:,:)
        integer                 ::  iBlock_Local

        do iBlock_local=1,size(BLocks)


            ! Block pointer
            Block1=>Blocks(iBlock_local)

            ! values pointer
            if (if_rk) then
                values=>Block1%values_rk
            else
                values=>Block1%values
            end if
            
            ! Left neighbour
            if (BLock1%iRank_l==iRank) then
                Block2=>Blocks(Block1%iBlock_l-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(nx+1:nx+ng,1:ny,1:nz,:)=&
                        Block1%values_rk(1:ng,1:ny,1:nz,:)
                else
                    Block2%values(nx+1:nx+ng,1:ny,1:nz,:)=&
                        Block1%values(1:ng,1:ny,1:nz,:)
                end if
            end if

            ! Right neighbour
            if (BLock1%iRank_r==iRank) then
                Block2=>Blocks(Block1%iBlock_r-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(-ng+1:0,1:ny,1:nz,:)=&
                        Block1%values_rk(nx-ng+1:nx,1:ny,1:nz,:)
                else
                    Block2%values(-ng+1:0,1:ny,1:nz,:)=&
                        Block1%values(nx-ng+1:nx,1:ny,1:nz,:)
                end if
            end if

            ! Back neighbour
            if (BLock1%iRank_b==iRank) then
                Block2=>Blocks(Block1%iBlock_b-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(1:nx,ny+1:ny+ng,1:nz,:)=&
                        Block1%values_rk(1:nx,1:ng,1:nz,:)
                else
                    Block2%values(1:nx,ny+1:ny+ng,1:nz,:)=&
                        Block1%values(1:nx,1:ng,1:nz,:)
                end if
            end if

            ! Front neighbour
            if (BLock1%iRank_f==iRank) then
                Block2=>Blocks(Block1%iBlock_f-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(1:nx,-ng+1:0,1:nz,:)=&
                        Block1%values_rk(1:nx,ny-ng+1:,1:nz,:)
                else
                    Block2%values(1:nx,-ng+1:0,1:nz,:)=&
                        Block1%values(1:nx,ny-ng+1:,1:nz,:)
                end if
            end if

            ! Up neighbour
            if (.not. Block1%if_top .and. BLock1%iRank_u==iRank) then
                Block2=>Blocks(Block1%iBlock_u-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(1:nx,1:ny,1:nz,-ng+1:0)=&
                        Block1%values_rk(1:nx,1:ny,nz-ng+1:nz,:)
                else
                    Block2%values(1:nx,1:ny,-ng+1:0,:)=&
                        Block1%values(1:nx,1:ny,nz-ng+1:nz,:)
                end if
            end if

            ! Down neighbour
            if (.not. Block1%if_bottom .and. BLock1%iRank_d==iRank) then
                Block2=>Blocks(Block1%iBlock_d-iBlock_range(1)+1)
                if (if_rk) then
                    Block2%values_rk(1:nx,1:ny,nz+1:nz+ng,:)=&
                        Block1%values_rk(1:nx,1:ny,1:ng,:)
                else
                    Block2%values(1:nx,1:ny,nz+1:nz+ng,:)=&
                        Block1%values(1:nx,1:ny,1:ng,:)
                end if
            end if
        end do

    end subroutine ModCommunication_Local

    subroutine ModCommunication_SendRecv(Blocks,if_rk)
        implicit none
        type(BlockType),target  ::  Blocks(:)
        type(BlockType),pointer ::  Block1
        logical,intent(in)      ::  if_rk
        real,pointer            ::  values(:,:,:,:)
        integer                 ::  iBlock_Local
        integer,allocatable     ::  requests(:)
        integer                 ::  irequest,ierr
        integer                 ::  status(MPI_STATUS_SIZE)

        ! First get the total number of sends
        allocate(requests(sum(Blocks(:)%n_sends)))
        irequest=0

        ! Send
        do iBlock_local=1,size(BLocks)

            ! Block pointer
            Block1=>Blocks(iBlock_local)

            ! values pointer
            if (if_rk) then
                values=>Block1%values_rk
            else
                values=>Block1%values
            end if
            
            ! Left neighbour
            if (BLock1%iRank_l/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(1:ng,1:ny,1:nz,:),&
                ny*nz*ng*5,mpi_real,BLock1%iRank_l,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if

            ! Right neighbour
            if (BLock1%iRank_r/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(nx-ng+1:nx,1:ny,1:nz,:),&
                ny*nz*ng*5,mpi_real,BLock1%iRank_r,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if

            ! Back neighbour
            if (BLock1%iRank_b/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(1:nx,1:ng,1:nz,:),&
                nx*nz*ng*5,mpi_real,BLock1%iRank_b,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if

            ! Front neighbour
            if (BLock1%iRank_f/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(1:nx,ny-ng+1:ny,1:nz,:),&
                nx*nz*ng*5,mpi_real,BLock1%iRank_f,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if

            ! Up neighbour
            if (.not. Block1%if_top .and. BLock1%iRank_u/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(1:nx,1:ny,nz-ng+1:nz,:),&
                ny*nx*ng*5,mpi_real,BLock1%iRank_u,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if

            ! Down neighbour
            if (.not. Block1%if_bottom .and. BLock1%iRank_d/=iRank) then
                irequest=irequest+1
                call MPI_ISEND(values(1:nx,1:ny,1:ng,:),&
                ny*nx*ng*5,mpi_real,BLock1%iRank_d,Block1%iBlock,MPI_COMM_WORLD,requests(irequest),ierr)
            end if
        end do

        ! Receive
        do iBlock_local=1,size(BLocks)

            ! Block pointer
            Block1=>Blocks(iBlock_local)

            ! values pointer
            if (if_rk) then
                values=>Block1%values_rk
            else
                values=>Block1%values
            end if

            ! Left neighbour
            if (BLock1%iRank_l/=iRank) then
                call MPI_RECV(values(-ng+1:0,1:ny,1:nz,:),&
                ny*nz*ng*5,mpi_real,Block1%iRank_l,Block1%iBlock_l,MPI_COMM_WORLD,status,ierr)
            end if

            ! Right neighbour
            if (BLock1%iRank_r/=iRank) then
                call MPI_RECV(values(nx+1:nx+ng,1:ny,1:nz,:),&
                ny*nz*ng*5,mpi_real,Block1%iRank_r,Block1%iBlock_r,MPI_COMM_WORLD,status,ierr)
            end if

            ! Back neighbour
            if (BLock1%iRank_b/=iRank) then
                call MPI_RECV(values(1:nx,-ng+1:0,1:nz,:),&
                nx*nz*ng*5,mpi_real,Block1%iRank_b,Block1%iBlock_b,MPI_COMM_WORLD,status,ierr)
            end if

            ! Front neighbour
            if (BLock1%iRank_f/=iRank) then
                call MPI_RECV(values(1:nx,ny+1:ny+ng,1:nz,:),&
                nx*nz*ng*5,mpi_real,Block1%iRank_f,Block1%iBlock_f,MPI_COMM_WORLD,status,ierr)
            end if

            ! Up neighbour
            if (.not. Block1%if_top .and. BLock1%iRank_u/=iRank) then
                call MPI_RECV(values(1:nx,1:ny,-ng+1:0,:),&
                ny*nx*ng*5,mpi_real,Block1%iRank_u,Block1%iBlock_u,MPI_COMM_WORLD,status,ierr)
            end if

            ! Down neighbour
            if (.not. Block1%if_bottom .and. BLock1%iRank_d/=iRank) then
                call MPI_RECV(values(1:nx,1:ny,nz+1:nz+ng,:),&
                ny*nx*ng*5,mpi_real,Block1%iRank_d,Block1%iBlock_d,MPI_COMM_WORLD,status,ierr)
            end if
        end do

        ! MPI_WAIT
        do irequest=1,size(requests)
            call MPI_WAIT(requests(irequest),status,ierr)
        end do
    end subroutine ModCommunication_SendRecv

end module ModCommunication