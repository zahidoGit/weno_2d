module wenoWeightsMod
    use numberPrecisionMod
    implicit none

    real(dp), parameter, dimension(-1:2,0:2) :: wenoCrj = reshape([11d0/6d0, -7d0/6d0,  1d0/3d0   &
                                                             & ,   1d0/3d0,   5d0/6d0, -1d0/6d0   &
                                                             & ,  -1d0/6d0,   5d0/6d0,  1d0/3d0   &
                                                             & ,   1d0/3d0,  -7d0/6d0,  11d0/6d0] &
                                                             & ,   [4,3], order=[2,1])

    real(dp), parameter, dimension(0:2) :: wenoD = [3d0/10d0, 3d0/5d0, 1d0/10d0]
    real(dp), parameter :: wenoEpsilon = 1d-6

    contains

    pure function wenoBeta(v)
        implicit none
        real(dp), dimension(-2:2), intent(in) :: v
        real(dp), dimension(0:2) :: wenoBeta
        real(dp) :: var1, var2
        real(dp) :: a, b, c, d, e

        var1 = ( v(0) - 2.d0*v(1) + v(2) )**2.d0
        var2 = ( 3.0d0*v(0) - 4.d0*v(1) + v(2) )**2.d0
        wenoBeta(0) = 13.0d0/12.0d0*var1 + 1.0d0/4.0d0*var2

        var1 = ( v(-1) - 2.d0*v(0) + v(1) )**2.d0
        var2 = ( v(-1) - v(1) )**2.d0
        wenoBeta(1) = 13.0d0/12.0d0*var1 + 1.0d0/4.0d0*var2

        var1 = ( v(-2) - 2.d0*v(-1) + v(0) )**2.d0
        var2 = ( v(-2) - 4.d0*v(-1) + 3.0d0*v(0) )**2.d0
        wenoBeta(2) = 13.0d0/12.0d0*var1 + 1.0d0/4.0d0*var2

    end function wenoBeta


    pure function wenoW(wBeta)
        implicit none
        real(dp), dimension(0:2), intent(in) :: wBeta
        real(dp), dimension(0:2) :: wenoAlpha
        real(dp), dimension(0:2) :: wenoW
        real(dp) :: temp
        integer :: r

        wenoAlpha = wenoD/(wenoEpsilon + wBeta)**2.0d0

        wenoW = wenoAlpha/sum(wenoAlpha)

    end function wenoW


    pure function wenoReconstruct(v) result(fHat)
        implicit none
        real(dp), dimension(-2:2), intent(in) :: v
        real(dp) :: fHat
        real(dp) :: vv
        real(dp), dimension(0:2) :: vr
        real(dp), dimension(0:2) :: beta
        real(dp), dimension(0:2) :: wenoWeights
        integer :: r, j

        do r = 0, 2
            vr(r) = sum(wenoCrj(r,:)*v(-r:-r+2))
        end do

        beta = wenoBeta(v)

        wenoWeights = wenoW(beta)

        vv = sum(wenoWeights * vr)

        fHat = vv

    end function wenoReconstruct


end module wenoWeightsMod
