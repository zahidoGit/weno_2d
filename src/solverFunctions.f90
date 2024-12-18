module solverFunctionMod
    use numberPrecisionMod
    use generalDataMod
    use wenoWeightsMod
    use proplemProperties, only : gammaFluid, cfl
    implicit none

    contains
    pure function getFaceReconstructedFlux(ff, uu, alp, neg) result(fFace)
        implicit none
        real(dp), dimension(-2:2), intent(in) :: ff, uu
        logical, intent(in) :: neg
        real(dp) :: fFace
        real(dp), dimension(-2:2) :: fp, fn
        real(dp), intent(in) :: alp

        if (neg) then
            fn = (ff - alp*uu)/2.0d0
            fn = fn(2:-2:-1)
            fFace = wenoReconstruct(fn)
        else
            fp = (ff + alp*uu)/2.0d0
            fFace = wenoReconstruct(fp)
        end if

    end function getFaceReconstructedFlux


    pure function calDomainAlpha() result(alp)
        implicit none
        real(dp) :: alp
        integer :: ii, jj

        alp = 0d0
        do jj = 1, nBlockY
            do ii = 1, nBlockX
                if (fluid(ii,jj)%zone .eq. "solid") cycle
                alp = max(alp, maxval(abs(fluid(ii,jj)%prim(:,:,idU)) + abs(fluid(ii,jj)%prim(:,:,idV))))
            end do
        end do

    end function calDomainAlpha


    pure function calDomainTimeStep() result(dt)
        implicit none
        real(dp) :: dt
        real(dp) :: sonicA, minGridSize, tempA, tempVel
        integer :: ii, jj, nX, nY

        sonicA = 0d0
        minGridSize = fluid(1,1)%dX
        do jj = 1, nBlockY
            do ii = 1, nBlockX
                if (fluid(ii,jj)%zone .eq. "solid") cycle
                nX = fluid(ii,jj)%nX
                nY = fluid(ii,jj)%nY
                ! tempA  = maxval(sqrt(gammaFluid * (fluid(ii,jj)%prim(1:nX,1:nY,idP)/fluid(ii,jj)%prim(1:nX,1:nY,idRho)) ))
                tempA  = maxval(fluid(ii,jj)%a(1:nX,1:nY))
                tempVel  = maxval(fluid(ii,jj)%prim(:,:,idU)**2.0d0 + fluid(ii,jj)%prim(:,:,idV)**2.0d0)
                tempVel = sqrt(tempVel)
                sonicA = max(sonicA, tempA + tempVel)
                minGridSize = min(minGridSize, min(fluid(ii,jj)%dX, fluid(ii,jj)%dY))
            end do
        end do

        dt = cfl*minGridSize/sonicA

    end function calDomainTimeStep




end module solverFunctionMod
