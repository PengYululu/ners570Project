! Unit test for ModParameters.f90
! and ModParamReader.f90

program test_Parameters

    use ModParameters
    use ModParamReader

    call Read_paramters('/Users/qplazmfree/Documents/570/project/input/PARAM.in.test')

    print *,Mx,My,Mz
    print *,nx,ny,nz,ng
    print *,superadiabaticity
    print *,Diffusion_h

end program test_Parameters