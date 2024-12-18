module proplemProperties
    use numberPrecisionMod
    use geometryMod, only : domainLength
    implicit none

    logical, parameter :: existsObsitacle = .true.

    ! sod shock-tube problem
    real(dp), parameter, dimension(4) :: shockDown = [1.0d0, 0.0d0, 0.0d0, 1000d0]   ! rho, u, v, p
    real(dp), parameter, dimension(4) :: shockUp   = [1.0d0, 0.0d0, 0.0d0, 0.01d0]
    real(dp), parameter :: endTime = 0.03d0
    real(dp), parameter :: x0 = domainLength*0.20

    real(dp), parameter :: shockFrontX = x0
    real(dp), parameter :: cfl = 0.9d0
    real(dp), parameter :: gammaFluid = 1.4d0

end module proplemProperties
