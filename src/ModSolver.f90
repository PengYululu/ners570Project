module ModSolver

    use ModBlock
    use ModDeviation

    contains

    subroutine GetR(iBlock, update_R)

        use ModParameter, only :: SuperAdiabacity, Xi, GasGamma
        implicit none
        type(BlockType), target :: iBlock
        integer :: direction
        real, intent(out) :: update_R(1:iBlock%nx, 1:iBlock%ny, 1:iBlock%nz, 1:5)
        real :: tmp(-iBlock%ng+1:iBlock%ng+iBlock%nx, -iBlock%ng+1:iBlock%ng+iBlock%ny, -iBlock%ng+1:iBlock%ng+iBlock%nz, 3)
        real :: p1(-iBlock%ng+1:iBlock%ng+iBlock%nx, -iBlock%ng+1:iBlock%ng+iBlock%ny, -iBlock%ng+1:iBlock%ng+iBlock%nz)

        p1 = GasGamma * iBlock%P0 * iBlock%values(:,:,:,1) / iBlock%rho0 + iBlock%values(:,:,:,5) 

        ! equation 1: partial(rho)/partial(t)
        do direction = 1, 3
            tmp(:,:,:,direction) = iBlock%rho0 * iBlock%values(:,:,:,direction)
            update_R(:,:,:,1) = update_R(:,:,:,1) - (1 / (8 * SuperAdiabacity * Xi**2) * &
                                ModDeviation_1st_O4_3D(tmp, iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng,&
                                iBlock%dx, iBlock%dy, iBlock%dz))
        end do

        ! equation 2: -(v dot nabla)v
        do direction = 1, 3
            update_R(:,:,:,direction+1) = iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,2) * ModDeviation_1st_O4_3D(iBlock%values(:,:,:,direction+1), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, direction) + &
                                        iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,3) * ModDeviation_1st_O4_3D(iBlock%values(:,:,:,direction+1), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, direction) + &
                                        iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,4) * ModDeviation_1st_O4_3D(iBlock%values(:,:,:,direction+1), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, direction)
        end do
        ! equation 2: -nabla(p1) /rho0 -rho1/rho0 * e_r
        do direction = 1, 3
            update_R(:,:,direction+1) = update_R(:,:,:,direction+1) - ModDeviation_1st_O4_3D(iBlock(:,:,:,1), , &
                                        iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, direction) / iBlock%rho0(1:iBlock%nx, 1:iBlock%ny, 1:iBlock%nz) &
                                        - iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,1) / iBlock%rho0(1:iBlock%nx, 1:iBlock%ny, 1:iBlock%nz) * iBlock%dz
        end do


        ! equation 3: vz/p0 - (v dot nabla)s1
        update_R(:,:,:,5) = iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,4) / iBlock%p0(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz) - &
                            iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,2) * ModDeviation_1st_O4_3D(iBlock(:,:,:,5), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, 1) - &
                            iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,3) * ModDeviation_1st_O4_3D(iBlock(:,:,:,5), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, 2) - &
                            iBlock%values(1:iBlock%nx,1:iBlock%ny,1:iBlock%nz,4) * ModDeviation_1st_O4_3D(iBlock(:,:,:,5), iBlock%nx, iBlock%ny, iBlock%nz, iBlock%ng, iBlock%dx, iBlock%dy, iBlock%dz, 3) 
        

    end subroutine GetR





end module ModSolver