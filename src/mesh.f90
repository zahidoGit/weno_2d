subroutine makeMesh()
    use numberPrecisionMod, only : dp
    use geometryMod
    use meshMod
    use generalDataMod, only : fluid, nComp, nPrim
    use proplemProperties, only : existsObsitacle
    implicit none
    real(dp) :: dx, dy
    real(dp) :: tempVar1, tempVar2
    integer :: ii, jj, i, j, lbX, ubX, lbY, ubY
    integer :: nCellBlockX, nCellBlockY
    integer :: iostat

    tempVar1 = (domainLength - obstacleLength)/2.0d0
    tempVar2 = (domainWidth - obstacleWidth)/2.0d0

    jj = 1
    fluid(1,jj)%startX = domainStartX
    fluid(1,jj)%startY = domainStartY
    fluid(1,jj)%endX = fluid(1,jj)%startX + tempVar1
    fluid(1,jj)%endY = fluid(1,jj)%startY + tempVar2

    fluid(2,jj)%startX = fluid(1,jj)%endX
    fluid(2,jj)%startY = fluid(1,jj)%startY
    fluid(2,jj)%endX   = fluid(2,jj)%startX + obstacleLength
    fluid(2,jj)%endY   = fluid(1,jj)%endY

    fluid(3,jj)%startX = fluid(2,jj)%endX
    fluid(3,jj)%startY = fluid(2,jj)%startY
    fluid(3,jj)%endX   = fluid(3,jj)%startX + tempVar1
    fluid(3,jj)%endY   = fluid(1,jj)%endY


    do jj = 2, nBlockY
        fluid(1,jj)%startX = domainStartX
        fluid(1,jj)%startY = fluid(1, jj-1)%endY
        fluid(1,jj)%endX = fluid(1,jj)%startX + tempVar1
        if (jj .eq. 2) then
            fluid(1,jj)%endY = fluid(1,jj)%startY + obstacleWidth
        else
            fluid(1,jj)%endY = fluid(1,jj)%startY + tempVar2
        end if

        fluid(2,jj)%startX = fluid(1,jj)%endX
        fluid(2,jj)%startY = fluid(1,jj)%startY
        fluid(2,jj)%endX   = fluid(2,jj)%startX + obstacleLength
        fluid(2,jj)%endY   = fluid(1,jj)%endY

        fluid(3,jj)%startX = fluid(2,jj)%endX
        fluid(3,jj)%startY = fluid(2,jj)%startY
        fluid(3,jj)%endX   = fluid(3,jj)%startX + tempVar1
        fluid(3,jj)%endY   = fluid(1,jj)%endY
    end do


    do jj = 1, nBlockY
        do ii = 1, nBlockX

            ! nCellBlockX = ceiling((fluid(ii,jj)%endX - fluid(ii,jj)%startX)/domainLength * nCellDomainX)
            ! nCellBlockY = ceiling((fluid(ii,jj)%endY - fluid(ii,jj)%startY)/domainWidth * nCellDomainY)

            dX = domainLength/nCellDomainX
            dY = domainWidth/nCellDomainY
            nCellBlockX = ceiling((fluid(ii,jj)%endX - fluid(ii,jj)%startX)/dX)
            nCellBlockY = ceiling((fluid(ii,jj)%endY - fluid(ii,jj)%startY)/dY)

            fluid(ii,jj)%nX = nCellBlockX
            fluid(ii,jj)%nY = nCellBlockY

            allocate(fluid(ii,jj)%vx(-nHalo : nCellBlockX+nHalo+1, -nHalo : nCellBlockY+nHalo+1))   ! node is added to the left
            allocate(fluid(ii,jj)%vy(-nHalo : nCellBlockX+nHalo+1, -nHalo : nCellBlockY+nHalo+1))

            allocate(fluid(ii,jj)%cx(-nHalo+1 : nCellBlockX+nHalo, -nHalo+1 : nCellBlockY+nHalo))
            allocate(fluid(ii,jj)%cy(-nHalo+1 : nCellBlockX+nHalo, -nHalo+1 : nCellBlockY+nHalo))

            dx = (fluid(ii,jj)%endX - fluid(ii,jj)%startX)/dfloat(nCellBlockX)
            dy = (fluid(ii,jj)%endY - fluid(ii,jj)%startY)/dfloat(nCellBlockY)

            lbX = -nHalo
            ubX = nCellBlockX+nHalo+1
            lbY = -nHalo
            ubY = nCellBlockY+nHalo+1

            fluid(ii,jj)%dx = dx
            fluid(ii,jj)%dy = dy

            tempVar1 = fluid(ii,jj)%startX
            tempVar2 = fluid(ii,jj)%startY
            do j = lbY, ubY
                fluid(ii,jj)%vx(:,j) = [(tempVar1 + dfloat(i)*dx, i = lbX, ubX)]
            end do

            do i = lbX, ubX
                fluid(ii,jj)%vy(i,:) = [(tempVar2 + dfloat(j)*dy, j = lbY, ubY)]
            end do

            fluid(ii,jj)%zone = "fluid"

        end do
    end do

    if(existsObsitacle) fluid(2,2)%zone = "solid"

    ! getting cell centers are vertices averages
    do jj = 1, nBlockY
        do ii = 1, nBlockX
            do j = -2, fluid(ii,jj)%nY + 3
                do i = -2, fluid(ii,jj)%nX + 3
                    fluid(ii,jj)%cx(i,j) = 0.25d0*(fluid(ii,jj)%vx(i-1,j-1)+fluid(ii,jj)%vx(i,j-1) &
                                          & + fluid(ii,jj)%vx(i,j)+fluid(ii,jj)%vx(i-1,j))

                    fluid(ii,jj)%cy(i,j) = 0.25d0*(fluid(ii,jj)%vy(i-1,j-1)+fluid(ii,jj)%vy(i,j-1) &
                                          & + fluid(ii,jj)%vy(i,j)+fluid(ii,jj)%vy(i-1,j))
                end do
            end do

        end do
    end do


    do jj = 1, nBlockY
        do ii = 1, nBlockX
            allocate(fluid(ii,jj)%prim (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nPrim))
            allocate(fluid(ii,jj)%compU(-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%compF(-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%compG(-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))

            allocate(fluid(ii,jj)%fHat (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%L    (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%fUX  (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%fDX  (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%fUY  (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%fDY  (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))
            allocate(fluid(ii,jj)%u0   (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo, nComp))

            allocate(fluid(ii,jj)%a    (-nHalo+1 : fluid(ii,jj)%nX+nHalo, -nHalo+1 : fluid(ii,jj)%nY+nHalo))
        end do
    end do

    call blockBoundaries()

end subroutine makeMesh
