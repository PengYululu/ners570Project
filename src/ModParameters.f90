Module ModParameters

    integer             ::  Mx,My,Mz                !   Number of blocks in each direction.
    integer             ::  nx,ny,nz                !   Number of grid points in one block
                                                    !   for each direction.
    integer             ::  x_range_global(2),&     !   The xyz range of the whole blo
                            y_range_global(2),&
                            z_range_global(2)

end module ModParameters