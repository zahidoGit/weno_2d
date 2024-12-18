subroutine initializeDomain()
    use generalDataMod
    use blocksMod
    use proplemProperties
    use fluxVectorMod
    use meshMod, only : nHalo
    use primitiveMod, only : getSonicSpeed
    implicit none
    integer :: ii, jj, i, j
    real(dp) :: en, ke, energy

    do jj = 1, nBlockY
        do ii = 1, nBlockX
            if (fluid(ii,jj)%zone .eq. "solid") cycle
            do j = -nHalo+1, fluid(ii,jj)%nY + nHalo
                do i = -nHalo+1, fluid(ii,jj)%nX + nHalo
                    if (fluid(ii,jj)%cx(i,j) .le. shockFrontX) then
                        fluid(ii,jj)%prim(i,j,idRho) = shockDown(1)
                        fluid(ii,jj)%prim(i,j,idU)   = shockDown(2)
                        fluid(ii,jj)%prim(i,j,idV)   = shockDown(3)
                        fluid(ii,jj)%prim(i,j,idP)   = shockDown(4)
                    else
                        fluid(ii,jj)%prim(i,j,idRho) = shockUp(1)
                        fluid(ii,jj)%prim(i,j,idU)   = shockUp(2)
                        fluid(ii,jj)%prim(i,j,idV)   = shockUp(3)
                        fluid(ii,jj)%prim(i,j,idP)   = shockUp(4)
                    end if
                end do
            end do
        end do
    end do

    do jj = 1, nBlockY
        do ii = 1, nBlockX

            if (fluid(ii,jj)%zone .eq. "solid") cycle

            do j = -nHalo+1, fluid(ii,jj)%nY + nHalo
                do i = -nHalo+1, fluid(ii,jj)%nX + nHalo
                    fluid(ii,jj)%compU(i,j,idRho) = fluid(ii,jj)%prim(i,j,idRho)
                    fluid(ii,jj)%compU(i,j,idU)   = fluid(ii,jj)%prim(i,j,idRho) * fluid(ii,jj)%prim(i,j,idU)
                    fluid(ii,jj)%compU(i,j,idV)   = fluid(ii,jj)%prim(i,j,idRho) * fluid(ii,jj)%prim(i,j,idV)

                    en = fluid(ii,jj)%prim(i,j,idP)/(gammaFluid-1.0d0)
                    en = en/fluid(ii,jj)%prim(i,j,idRho)
                    ke = 0.5d0 * (fluid(ii,jj)%prim(i,j,idU)**2.0d0 + fluid(ii,jj)%prim(i,j,idV)**2.0d0)

                    energy = ke + en
                    energy = energy * fluid(ii,jj)%prim(i,j,idRho)

                    fluid(ii,jj)%compU(i,j,idE) = energy

                    fluid(ii,jj)%compF(i,j,:) = getF(fluid(ii,jj)%prim(i,j,:))
                    fluid(ii,jj)%compG(i,j,:) = getG(fluid(ii,jj)%prim(i,j,:))

                    fluid(ii,jj)%a(i,j) = getSonicSpeed(fluid(ii,jj)%prim(i,j,:))

                end do
            end do
        end do
    end do

    call blockBoundaries()

end subroutine initializeDomain
