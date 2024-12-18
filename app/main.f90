program main
    use geometryMod
    use wenoWeightsMod
    use generalDataMod
    use primitiveMod
    use fluxVectorMod
    use solverFunctionMod
    use proplemProperties, only : endTime
    implicit none
    integer :: ii, jj, i, j, pp, tt
    integer :: lbY, ubY, lbX, ubX
    real(dp), dimension(-2:2) :: fp, fn, ff, uu, vel
    real(dp) :: time, dt, alp
    integer :: saveCounter, iterationCounter
    integer, parameter :: saveFrequency = 50
    real(dp), parameter :: zero = 0d0

    call makeMesh()

    call saveStructuredTecPlot(zero)

    call initializeDomain()

    call saveStructuredTecPlot(zero)

    dt = 1d-15
    time = 0d0
    saveCounter = 0
    iterationCounter = 0
    do while (.not. abs(time - endtime) .le. epsilon(endTime))

        dt = calDomainTimeStep()
        if (iterationCounter .le. 10) then
            dt = dt/10d0
        end if

        iterationCounter = iterationCounter + 1

        if (endTime-time .le. dt) then
            print*, "adjusting the time step"
            dt = abs(endTime-time)
        end if

        alp = 0d0
        alp = calDomainAlpha()
        do jj = 1, nBlockY
            do ii = 1, nBlockY

                if (fluid(ii,jj)%zone .eq. "solid") cycle

                fluid(ii,jj)%u0 = fluid(ii,jj)%compU

                call getWenoFluxes()
                call rhs()
                fluid(ii,jj)%compU = fluid(ii,jj)%u0 + fluid(ii,jj)%L

                call getWenoFluxes()
                call rhs()
                fluid(ii,jj)%compU = 0.75d0*fluid(ii,jj)%u0 + 0.25d0*(fluid(ii,jj)%compU + fluid(ii,jj)%L)

                call getWenoFluxes()
                call rhs()
                fluid(ii,jj)%compU = 1d0/3d0*fluid(ii,jj)%u0 + 2d0/3d0*(fluid(ii,jj)%compU + fluid(ii,jj)%L)
            end do
        end do
        time = time + dt
        print*, time, dt

        saveCounter = saveCounter + 1
        if (mod(saveCounter, saveFrequency) .eq. 0) then
            call saveStructuredTecPlot(time)
            saveCounter = 0
        end if
    end do

    call saveStructuredTecPlot(zero)
    call saveStructuredTecPlot(time)

    call deallocateAll()

    contains
    subroutine getWenoFluxes()
        implicit none
        integer :: nX, nY

        nX = fluid(ii,jj)%nX
        nY = fluid(ii,jj)%nY

        call blockBoundaries()


        do concurrent (j=0:fluid(ii,jj)%nY+1, i=0:fluid(ii,jj)%nX+1)

            fluid(ii,jj)%prim (i,j,:) = getPrimitive(fluid(ii,jj)%compU(i,j,:))

            fluid(ii,jj)%a    (i,j)   = getSonicSpeed(fluid(ii,jj)%prim(i,j,:))
            fluid(ii,jj)%compF(i,j,:) = getF(fluid(ii,jj)%prim(i,j,:))
            fluid(ii,jj)%compG(i,j,:) = getG(fluid(ii,jj)%prim(i,j,:))


            do concurrent (pp = 1 : nComp)

                uu = fluid(ii,jj)%compU(i-2:i+2, j, pp)
                ff = fluid(ii,jj)%compF(i-2:i+2, j, pp)

                fluid(ii,jj)%fUX(i,j,pp) = getFaceReconstructedFlux(ff, uu, alp, .false.)
                fluid(ii,jj)%fDX(i,j,pp) = getFaceReconstructedFlux(ff, uu, alp, .true. )

                uu = fluid(ii,jj)%compU(i, j-2:j+2, pp)
                ff = fluid(ii,jj)%compG(i, j-2:j+2, pp)

                fluid(ii,jj)%fUY(i,j,pp) = getFaceReconstructedFlux(ff, uu, alp, .false.)
                fluid(ii,jj)%fDY(i,j,pp) = getFaceReconstructedFlux(ff, uu, alp, .true. )

            end do
        end do

    end subroutine getWenoFluxes


    subroutine rhs()
        implicit none
        real(dp) :: dX, dY
        real(dp), DIMENSION(nComp) :: varX, varY

        dX = fluid(ii,jj)%dX
        dY = fluid(ii,jj)%dY

        varX = 0d0
        varY = 0d0

        ! do j = 1, fluid(ii,jj)%nY
        !     do i = 1, fluid(ii,jj)%nX

        do concurrent (j=0:fluid(ii,jj)%nY, i=0:fluid(ii,jj)%nX)

                varX = fluid(ii,jj)%fUX(i,j,:) - fluid(ii,jj)%fUX(i-1,j,:)
                varX = varX+ fluid(ii,jj)%fDX(i+1,j,:) - fluid(ii,jj)%fDX(i,j,:)

                varY = fluid(ii,jj)%fUY(i,j,:) - fluid(ii,jj)%fUY(i,j-1,:)
                varY = varY + fluid(ii,jj)%fDY(i,j+1,:) - fluid(ii,jj)%fDY(i,j,:)

                fluid(ii,jj)%L(i,j,:) = -dt*(varX/dX + varY/dY)
                ! fluid(ii,jj)%L(i,j,:) = -dt*(varX/approxDX + varY/approxDY)

            end do
        ! end do

    end subroutine rhs

end program main
