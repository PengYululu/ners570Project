module ModDiffusion

    use ieee_arithmetic

    contains

    subroutine ModDiffusion_Aritificial(primitive,ni,nj,nk,ng,dxi,dxj,dxk,EQN_update_R,&
        c_s,h)

        implicit none
        integer,intent(in)              ::  ni,nj,nk,ng
        real,intent(in)                 ::  h
        real,intent(in)                 ::  dxi,dxj,dxk
        real,intent(in)                 ::  primitive(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,5) 
        real,intent(in)                 ::  c_s(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)

        integer                         ::  direction1,direction2
        real                            ::  c(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng)
                                            
        real                            ::  d_primitive(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,5)
        real,allocatable                ::  flux(:,:,:,:),phi(:,:,:,:)
        real,intent(inout)              ::  EQN_update_R(1:ni,1:nj,1:nk,1:5)
        integer ::  ivar,i,j,k

        do direction1=1,3

            

            ! get total speed c
            c=abs(primitive(:,:,:,direction1+1))+c_s

            ! use minmod to find \Delta u
            call ModDiffusion_minmod(5,ni,nj,nk,ng,direction1,primitive,d_primitive)

            ! get phi & flux
            select case(direction1)
            case(1)
                allocate(flux(ni+1,nj,nk,5),phi(ni+1,nj,nk,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(1:ni+1,:,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)/&
                    (primitive(1:ni+1,1:nj,1:nk,:)-primitive(0:ni,1:nj,1:nk,:))-1))
                do i=1,ni+1; do j=1,nj; do k=1,nk; do ivar=1,5
                    phi(i,j,k,ivar)=merge(0.0,phi(i,j,k,ivar),ieee_is_nan(phi(i,j,k,ivar)))
                end do; end do; end do; end do 
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(1:ni+1,:,:,:)=&
                    -0.5*phi*(primitive(1:ni+1,1:nj,1:nk,:)-d_primitive(1:ni+1,1:nj,1:nk,:)*0.5-&
                    primitive(0:ni,1:nj,1:nk,:)-d_primitive(0:ni,1:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(1:ni+1,:,:)+c(0:ni,:,:))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(1:ni+1,:,:,direction2+1)+primitive(0:ni,:,:,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                ! update EQN_update_R
                do ivar=1,5
                    EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                        (flux(1:ni,:,:,ivar)-flux(2:ni+1,:,:,ivar))/dxi
                end do
            case(2)
                allocate(flux(ni,nj+1,nk,5),phi(ni,nj+1,nk,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,1:nj+1,:,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)/&
                    (primitive(1:ni,1:nj+1,1:nk,:)-primitive(1:ni,0:nj,1:nk,:))-1))
                do i=1,ni; do j=1,nj+1; do k=1,nk; do ivar=1,5
                    phi(i,j,k,ivar)=merge(0.0,phi(i,j,k,ivar),ieee_is_nan(phi(i,j,k,ivar)))
                end do; end do; end do; end do 
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,1:nj+1,:,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj+1,1:nk,:)-d_primitive(1:ni,1:nj+1,1:nk,:)*0.5-&
                    primitive(1:ni,0:nj,1:nk,:)-d_primitive(1:ni,0:nj,1:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,1:nj+1,:)+c(:,0:nj,:))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(:,1:nj+1,:,direction2+1)+primitive(:,0:nj,:,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                do ivar=1,5
                    EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                            (flux(:,1:nj,:,ivar)-flux(:,2:nj+1,:,ivar))/dxj
                end do
            case(3)
                allocate(flux(ni,nj,nk+1,5),phi(ni,nj,nk+1,5))

                ! first get \Phi_{h}
                ! phi = max[0, 1+((u_r-u_l)/(u_i+1-u_i)-1)],
                ! for the i+1/2 face.
                phi(:,:,1:nk+1,:)=&
                    max(0.0,1.0+h*((primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)/&
                    (primitive(1:ni,1:nj,1:nk+1,:)-primitive(1:ni,1:nj,0:nk,:))-1))
                
                do i=1,ni; do j=1,nj; do k=1,nk+1; do ivar=1,5
                    phi(i,j,k,ivar)=merge(0.0,phi(i,j,k,ivar),ieee_is_nan(phi(i,j,k,ivar)))
                end do; end do; end do; end do
                
                ! then get the flux.
                ! f_{i+1/2}=-0.5 * c_{i+1/2}*phi*(u_r-u_l)
                flux(:,:,1:nk+1,:)=&
                    -0.5*phi*(primitive(1:ni,1:nj,1:nk+1,:)-d_primitive(1:ni,1:nj,1:nk+1,:)*0.5-&
                    primitive(1:ni,1:nj,0:nk,:)-d_primitive(1:ni,1:nj,0:nk,:)*0.5)
                
                ! multiply flux by c_{i+1/2}
                do ivar=1,5
                    flux(:,:,:,ivar)=flux(:,:,:,ivar)*(c(:,:,1:nk+1)+c(:,:,0:nk))*0.5
                end do

                do direction2=1,3
                    flux(:,:,:,direction2+1)=flux(:,:,:,direction2+1)+&
                        (primitive(:,:,1:nk+1,direction2+1)+primitive(:,:,0:nk,direction2+1))*&
                        0.5*flux(:,:,:,1)
                end do

                

                do ivar=1,5
                    EQN_update_R(:,:,:,ivar)=EQN_update_R(:,:,:,ivar)+&
                            (flux(:,:,1:nk,ivar)-flux(:,:,2:nk+1,ivar))/dxk
                end do
            end select
            deallocate(flux,phi)
        end do
    end subroutine ModDiffusion_Aritificial

    subroutine ModDiffusion_minmod(nvars,ni,nj,nk,ng,direction,u,d_u)
        implicit none
        integer,intent(in)      ::  nvars,ni,nj,nk,ng,direction
        real,intent(in)         ::  u(-ng+1:ni+ng,-ng+1:nj+ng,-ng+1:nk+ng,nvars)
        real,intent(out)        ::  d_u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars)
        real                    ::  u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    a(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    b(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars),&
                                    c(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,nvars)

        ! first step: get the u at left and right
        ! for every grid we care 

        select case(direction)
        case(1)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+1:ni+ng-2,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+3:ni+ng,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)
        case(2)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+1:nj+ng-2,-ng+2:nk+ng-1,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+3:nj+ng,-ng+2:nk+ng-1,:)
        case(3)
            u_l(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+1:nk+ng-2,:)
            u_r(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)=&
                u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+3:nk+ng,:)
        end select

        ! then get the three elements in minmod:
        ! \Delta_u = minmod( (u_r-u_l)/2 , 2(u_r-u) , 2(u-u_l) )
            
        a=(u_r-u_l)*0.5
        b=2.0*(u_r-u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:))
        c=2.0*(u(-ng+2:ni+ng-1,-ng+2:nj+ng-1,-ng+2:nk+ng-1,:)-u_l)

        ! initialize d_u and get it using minmod

        d_u=0.5
        d_u=(sign(d_u,max(0.0,a*b)*max(0.0,b*c)*1.e10-1.e-30)+0.5)*&
            sign(min(abs(a),abs(b),abs(c)),a)
    end subroutine ModDiffusion_minmod


end module ModDiffusion
