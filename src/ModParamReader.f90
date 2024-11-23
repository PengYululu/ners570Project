module ModParamReader

    use ModParameters, only: Mx,My,Mz,nBlocks_global,&
        nx,ny,nz,ng,x_range_global,y_range_global,z_range_global,&
        Xi,TypeDiffusion,Diffusion_h,UseGravity,TypeGas,&
        TypeSuperAdiabaticity,SuperAdiabaticity,UseRSST,&
        TypeInitiation,RandomScale

    ! Generic interface to read a parameter from a line
    !
    ! Usage:
    !       call read_var(line,parameter)
    !
    ! Keywords:
    !       line        :   The character type of one line
    !       parameter   :   The target parameter to read into

    contains

    subroutine Read_paramters(filename)
        implicit none
        character(len=*),intent(in)     ::  filename            ! Parameter input filename
        character(len=256)              ::  line,var            ! Which line it is
        character(len=14)               ::  name_sub='Read_paramters'
        integer                         ::  unit, ios           ! For reading
    
        ! Open the input file
        unit = 10
        open(unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: ", filename
            stop
        end if
    
        ! Read the file line by line
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit  ! End of file
    
            ! If starts with # then enter this command block
            if (line(1:1) == "#") then
    
                ! See which command it is
                select case (trim(adjustl(line)))
                case ("#GRAVITY")

                    ! UseGravity
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("T","t","true","True","TRUE")
                        UseGravity = .true.
                    case("F","f","false","False","FALSE")
                        UseGravity = .false.
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown UseGravity=",trim(adjustl(var))
                        stop 1
                    end select
                    
                case ("#DOMAIN")
                    read(unit, *, iostat=ios) x_range_global(1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading xmin"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) y_range_global(1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ymin"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) z_range_global(1)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading zmin"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) x_range_global(2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading xmax"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) y_range_global(2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ymax"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) z_range_global(2)
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading zmax"
                        stop 1
                    end if
                
                case ("#BLOCKS")
                    read(unit, *, iostat=ios) Mx
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading Mx"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) My
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading My"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) Mz
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading Mz"
                        stop 1
                    end if

                case ("#GRID")
                    read(unit, *, iostat=ios) nx
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nx"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) ny
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ny"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) nz
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading nz"
                        stop 1
                    end if
                    read(unit, *, iostat=ios) ng
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading ng"
                        stop 1
                    end if

                case("#GAS")
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("Ideal",'ideal','IDEAL')
                        TypeGas="Ideal"
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown TypeGas=",trim(adjustl(var))
                        stop 1
                    end select

                case("#RSST")
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("T","t","true","True","TRUE")
                        UseRSST = .true.
                    case("F","f","false","False","FALSE")
                        UseRSST = .false.
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown UseRSST=",trim(adjustl(var))
                        stop 1
                    end select

                case("#DIFFUSION")
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("artificial","ARTIFICIAL","Aritificial")
                        TypeDiffusion="Aritificial"
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown TypeDiffusion=",trim(adjustl(var))
                        stop 1
                    end select

                    read(unit, *, iostat=ios) Diffusion_h
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading Diffusion_h"
                        stop 1
                    end if

                case("#SUPERADIABATICITY")
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("const","CONST","Const","Constant","CONSTANT","constant")
                        TypeSuperAdiabaticity="const"
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown TypeSuperAdiabaticity=",trim(adjustl(var))
                        stop 1
                    end select

                    read(unit, *, iostat=ios) SuperAdiabaticity
                    if (ios/=0) then
                        write(*,*) "Error from ",name_sub,": Error reading SuperAdiabaticity"
                        stop 1
                    end if
                
                case("#INITIATION")
                    read(unit, *, iostat=ios) var
                    select case (trim(adjustl(var)))
                    case("random","Random","RANDOM","r","R")
                        TypeInitiation="random"
                    case default
                        write(*,*) "Error from ",name_sub,": Unknown TypeInitiation=",trim(adjustl(var))
                        stop 1
                    end select

                    select case(TypeInitiation)
                    case("random")
                        read(unit, *, iostat=ios) RandomScale
                        if (ios/=0) then
                            write(*,*) "Error from ",name_sub,": Error reading RandomScale"
                            stop 1
                        end if
                    end select

                end select
            end if
        end do
    
        ! Close the file
        close(unit)
    end subroutine Read_paramters

end module ModParamReader