module ModInterpolate
    contains

    ! Interpolate the 3D array into a 1D series of points
    !
    ! Input:    f(:,:,:)    : The array to be interpolated
    !           ni,nj,nk,ng : The grid numbers
    !           xi,xj,xk    : The coordinate of f
    !           nout        : N of output points
    !           xijk_out    : Xijk of output points
    !           fout        : Output
    !
    subroutine ModInterpolate_3D_to_1D(f,ni,nj,nk,ng,xi,xj,xk,nout,xijk_out,fout)
        implicit none
        integer,intent(in)  ::  ni,nj,nk,ng
        integer,intent(in)  ::  nout
        real,intent(in)     ::  xi(-ng+1:ni+ng),&
                                xj(-ng+1:nj+ng),&
                                xk(-ng+1:nk+ng),&
                                f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
        real,intent(in)     ::  xijk_out(nout,3)
        real,intent(out)    ::  fout(nout)

        integer             ::  i                       ! Index of the point
        real                ::  dxi,dxj,dxk             ! Grid size for three coordinates        
        real                ::  posi_i,posi_j,posi_k    ! The index positions for one point
        real                ::  w_i,w_j,w_k             ! The weight. Used for trilinear interpolation
        integer             ::  posi_i_integer,&        ! The integer part of index positions
                                posi_j_integer,&
                                posi_k_integer

        ! First get the grid sizes
        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        ! Loop all the points
        do i=1,nout

            ! Get the index positions of this point
            posi_i=-ng+1.+(xijk_out(i,1)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xijk_out(i,2)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xijk_out(i,3)-xk(-ng+1))/dxk

            ! Get their integer part
            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            ! The float part is what I called the weight
            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            ! Trilinear interpolation based on the surrounding 8 points
            fout(i)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  )*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  )*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  )*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  )*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1)*(   w_i)*(   w_j)*(   w_k)
            
        end do
    end subroutine ModInterpolate_3D_to_1D

    ! Interpolate the 3D(x,y,z)*1D(i of vars) array into 
    ! 1D(i of points)*1D(i of vars)
    !
    ! Input:    f(:,:,:,:)  : The array to be interpolated
    !           ni,nj,nk,ng : The grid numbers
    !           nvar        : N of variables. The fourth dimension of f
    !           xi,xj,xk    : The coordinate of f
    !           nout        : N of output points
    !           xijk_out    : Xijk of output points
    !           fout        : Output
    !
    subroutine ModInterpolate_3D1D_to_1D1D(f,nvar,ni,nj,nk,ng,xi,xj,xk,nout,xijk_out,fout)
        implicit none
        integer,intent(in)  ::  nvar,ni,nj,nk,ng,nout
        real,intent(in)     ::  xi(-ng+1:ni+ng),&
                                xj(-ng+1:nj+ng),&
                                xk(-ng+1:nk+ng)
        real,intent(in)     ::  f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvar)
        real,intent(in)     ::  xijk_out(nout,3)

        real                ::  fout(nout,nvar)

        integer             ::  i                       ! Index of the point
        real                ::  dxi,dxj,dxk             ! Grid size for three coordinates        
        real                ::  posi_i,posi_j,posi_k    ! The index positions for one point
        real                ::  w_i,w_j,w_k             ! The weight. Used for trilinear interpolation
        integer             ::  posi_i_integer,&        ! The integer part of index positions
                                posi_j_integer,&
                                posi_k_integer

        ! First get the grid sizes
        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        ! Loop all the points
        do i=1,nout

            ! Get the index positions of this point
            posi_i=-ng+1.+(xijk_out(i,1)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xijk_out(i,2)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xijk_out(i,3)-xk(-ng+1))/dxk

            ! Get their integer part
            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            ! The float part is what I called the weight
            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            ! Trilinear interpolation based on the surrounding 8 points
            fout(i,:)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  ,:)*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  ,:)*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  ,:)*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1,:)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  ,:)*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1,:)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1,:)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1,:)*(   w_i)*(   w_j)*(   w_k)
        end do
    end subroutine ModInterpolate_3D1D_to_1D1D

    ! Interpolate the 3D(x,y,z)*1D(i of vars) array into 
    ! 2D(i,j of map)*1D(i of vars)
    !
    ! Input:    f(:,:,:,:)  : The array to be interpolated
    !           ni,nj,nk,ng : The grid numbers
    !           nvar        : N of variables. The fourth dimension of f
    !           xi,xj,xk    : The coordinate of f
    !           nij_out     : ni & nj of the output map
    !           xi(j,k)_out : Coordinates of the map
    !           fout        : Output
    !
    subroutine ModInterpolate_3D1D_to_2D1D(f,nvar,ni,nj,nk,ng,xi,xj,xk,nij_out,xi_out,xj_out,xk_out,fout)
        implicit none
        integer,intent(in)  ::  nvar,ni,nj,nk,ng,nij_out(2)
        real,intent(in)     ::  xi(-ng+1:ni+ng),&
                                xj(-ng+1:nj+ng),&
                                xk(-ng+1:nk+ng)
        real,intent(in)     ::  f(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvar)
        real,intent(in)     ::  xi_out(nij_out(1),nij_out(2)),&
                                xj_out(nij_out(1),nij_out(2)),&
                                xk_out(nij_out(1),nij_out(2))
        real,intent(out)    ::  fout(nij_out(1),nij_out(2),nvar)

        integer             ::  i_out,j_out             ! Indice of point in map
        real                ::  dxi,dxj,dxk             ! Grid size for three coordinates        
        real                ::  posi_i,posi_j,posi_k    ! The index positions for one point
        real                ::  w_i,w_j,w_k             ! The weight. Used for trilinear interpolation
        integer             ::  posi_i_integer,&        ! The integer part of index positions
                                posi_j_integer,&
                                posi_k_integer
        
        ! First get the grid sizes
        dxi=xi(2)-xi(1); dxj=xj(2)-xj(1); dxk=xk(2)-xk(1)

        ! Loop the whole map
        do i_out=1,nij_out(1); do j_out=1,nij_out(2)

            ! Get the index positions of this point
            posi_i=-ng+1.+(xi_out(i_out,j_out)-xi(-ng+1))/dxi
            posi_j=-ng+1.+(xj_out(i_out,j_out)-xj(-ng+1))/dxj
            posi_k=-ng+1.+(xk_out(i_out,j_out)-xk(-ng+1))/dxk

            ! Get their integer part
            posi_i_integer=floor(posi_i)
            posi_j_integer=floor(posi_j)
            posi_k_integer=floor(posi_k)

            ! The float part is what I called the weight
            w_i=posi_i-posi_i_integer
            w_j=posi_j-posi_j_integer
            w_k=posi_k-posi_k_integer

            ! Trilinear interpolation based on the surrounding 8 points
            fout(i_out,j_out,:)=&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer  ,:)*(1.-w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer  ,:)*(   w_i)*(1.-w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer  ,:)*(1.-w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer  ,posi_k_integer+1,:)*(1.-w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer  ,:)*(   w_i)*(   w_j)*(1.-w_k)+&
                f(posi_i_integer  ,posi_j_integer+1,posi_k_integer+1,:)*(1.-w_i)*(   w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer  ,posi_k_integer+1,:)*(   w_i)*(1.-w_j)*(   w_k)+&
                f(posi_i_integer+1,posi_j_integer+1,posi_k_integer+1,:)*(   w_i)*(   w_j)*(   w_k)
        end do; end do
    end subroutine ModInterpolate_3D1D_to_2D1D
end module ModInterpolate