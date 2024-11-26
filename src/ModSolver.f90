module ModSolver

    use ModParameters, only : SuperAdiabaticity, Xi, GasGamma, &
        nx,ny,nz,ng,Diffusion_h
    use ModBlock
    use ModDerivative
    use ModDiffusion

    contains

    subroutine GetR(Block1, update_R, if_rk)
        implicit none
        type(BlockType), target     ::  Block1
        logical,intent(in)          ::  if_rk
        real,intent(out)            ::  update_R(1:nx, 1:ny, 1:nz, 1:5)
        real,pointer                ::  values(:,:,:,:)
        integer                     ::  direction1,direction2
        real                        ::  tmp(-ng+1:ng+nx, -ng+1:ng+ny, -ng+1:ng+nz, 3)
        real                        ::  p1(-ng+1:ng+nx, -ng+1:ng+ny, -ng+1:ng+nz)

        ! Set if to por bottom boundary condition
        ! if necessary
        call SetBoundary(Block1,if_rk)

        ! Assign the values pointer
        if (if_rk) then
            values=>Block1%values_rk
        else
            values=>Block1%values
        end if

        ! Set p1
        p1 = GasGamma * Block1%P0 * values(:,:,:,1) / Block1%rho0 + values(:,:,:,5) 
        update_R=0.0

        ! equation 1: partial(rho)/partial(t)
        do direction1 = 1, 3
            tmp(:,:,:,direction1) = Block1%rho0 * values(:,:,:,direction1+1)
            update_R(:,:,:,1) = update_R(:,:,:,1) - 1.0 / (8.0 * SuperAdiabaticity * Xi**2) * &
                                ModDerivative_1st_O4_3D(tmp, nx, ny, nz, ng,&
                                Block1%dx, Block1%dy, Block1%dz,direction1)
        end do

        ! equation 2: -(v dot nabla)v
        do direction1 = 1, 3
            do direction2 = 1,3
                update_R(:,:,:,direction1+1) = &
                    values(1:nx,1:ny,1:nz,direction2+1) * &
                    ModDerivative_1st_O4_3D(values(:,:,:,direction1+1), &
                    nx, ny, nz, ng, Block1%dx, Block1%dy, Block1%dz, direction2)
            end do
        end do

        ! equation 2: -nabla(p1) /rho0
        do direction1 = 1, 3
            update_R(:,:,:,direction1+1) = &
                update_R(:,:,:,direction1+1) - &
                ModDerivative_1st_O4_3D(values(:,:,:,1), &
                nx, ny, nz, ng, Block1%dx, Block1%dy, Block1%dz, direction1) / &
                 Block1%rho0(1:nx, 1:ny, 1:nz)
        end do

        ! equation 2: -rho1/rho0 * e_z
        update_R(:,:,:,4)=update_R(:,:,:,4)-values(1:nx,1:ny,1:nz,1) / Block1%rho0(1:nx, 1:ny, 1:nz)

        ! equation 3: vz/p0
        update_R(:,:,:,5)=values(1:nx,1:ny,1:nz,4) / Block1%p0(1:nx,1:ny,1:nz)

        ! equation 3: - (v dot nabla)s1
        do direction2=1,3
            update_R(:,:,:,5) = update_R(:,:,:,5) - &
                values(1:nx,1:ny,1:nz,direction2+1) * &
                ModDerivative_1st_O4_3D(values(:,:,:,5), &
                nx, ny, nz, ng, Block1%dx, Block1%dy, Block1%dz, direction2)
        end do

        ! Artificial diffusion

        call ModDiffusion_Aritificial(values,nx,ny,nz,ng,&
            Block1%dx,Block1%dy,Block1%dz,update_R,&
            1./Xi*sqrt(GasGamma/(8.*SuperAdiabaticity)*Block1%p0/Block1%rho0),Diffusion_h)
    end subroutine GetR

    subroutine SetBoundary(Block1,if_rk)
        implicit none
        type(BlockType), target     ::  Block1
        real,pointer                ::  values(:,:,:,:)
        logical,intent(in)          ::  if_rk
        integer                     ::  i

        ! Assign the values pointer
        if (if_rk) then
            values=>Block1%values_rk
        else
            values=>Block1%values
        end if

        ! Top 
        if (Block1%if_top) then
            do i=nz+1,nz+ng
                values(:,:,i,[1,3,4])=values(:,:,2*nz+1-i,[1,3,4])
                values(:,:,i,[2,5])=-values(:,:,2*nz+1-i,[2,5])
            end do
        end if

        ! Bottom
        if (Block1%if_bottom) then
            do i=-ng+1,0
                values(:,:,i,[1,3,4])=values(:,:,1-i,[1,3,4])
                values(:,:,i,[2,5])=-values(:,:,1-i,[2,5])
            end do
        end if
    end subroutine SetBoundary
end module ModSolver