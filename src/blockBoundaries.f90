subroutine blockBoundaries()
    use numberPrecisionMod
    use meshMod, only : nHalo
    use generalDataMod, only : fluid
    use geometryMod
    use proplemProperties, only : existsObsitacle
    implicit none
    integer :: ii, jj, i, j, lbX, ubX, lbY, ubY
    integer :: nX, nY

    ! right face
    call zeroGradientBlockEast(3,1)
    call zeroGradientBlockEast(3,2)
    call zeroGradientBlockEast(3,3)

    ! left face
    call zeroGradientBlockWest(1,1)
    call zeroGradientBlockWest(1,2)
    call zeroGradientBlockWest(1,3)

    ! internal faces
    call interfaceBlockEast(1,1)
    call interfaceBlockEast(1,3)
    call interfaceBlockWest(2,1)
    call interfaceBlockWest(2,3)

    ! internal faces
    call interfaceBlockEast(2,1)
    call interfaceBlockEast(2,3)
    call interfaceBlockWest(3,1)
    call interfaceBlockWest(3,3)

    ! internal faces
    call interfaceBlockNorth(1,1)
    call interfaceBlockNorth(3,1)
    call interfaceBlockSouth(1,2)
    call interfaceBlockSouth(3,2)

    ! internal faces
    call interfaceBlockNorth(1,2)
    call interfaceBlockNorth(3,2)
    call interfaceBlockSouth(1,3)
    call interfaceBlockSouth(3,3)

    ! bottom face
    call zeroGradientBlockSouth(1,1)
    call zeroGradientBlockSouth(2,1)
    call zeroGradientBlockSouth(3,1)

    ! top face
    call zeroGradientBlockNorth(1,3)
    call zeroGradientBlockNorth(2,3)
    call zeroGradientBlockNorth(3,3)

    if (existsObsitacle) then
        ! around the solid
        call noSlipBlockEast (1,2)
        call noSlipBlockWest (3,2)
        call noSlipBlockNorth(2,1)
        call noSlipBlockSouth(2,3)
    else
        ! around the solid
        call interfaceBlockEast (1,2)
        call interfaceBlockWest (3,2)
        call interfaceBlockNorth(2,1)
        call interfaceBlockSouth(2,3)

        ! solid internal faces
        call interfaceBlockNorth(2,2)
        call interfaceBlockSouth(2,2)
        call interfaceBlockEast (2,2)
        call interfaceBlockWest (2,2)
    end if

end subroutine blockBoundaries
