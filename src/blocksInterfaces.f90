subroutine interfaceBlockEast(ii, jj)
    use numberPrecisionMod
    use meshMod, only : nHalo
    use generalDataMod, only : fluid
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubX

    ubX = fluid(ii,jj)%nX + nHalo

    fluid(ii,jj)%compU(ubX-2:ubX,:,:) = fluid(ii+1,jj)%compU(1:3,:,:)

end subroutine interfaceBlockEast


subroutine interfaceBlockWest(ii, jj)
    use numberPrecisionMod
    use meshMod, only : nHalo
    use generalDataMod, only : fluid
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbX
    integer :: nX

    lbX = -nHalo + 1
    nX = fluid(ii-1,jj)%nX

    fluid(ii,jj)%compU(lbX:0,:,:) = fluid(ii-1,jj)%compU(nX+lbX:nX,:,:)

end subroutine interfaceBlockWest


subroutine interfaceBlockNorth(ii, jj)
    use numberPrecisionMod
    use meshMod, only : nHalo
    use generalDataMod, only : fluid
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubY
    integer :: nY

    nY = fluid(ii,jj)%nY
    ubY = nY + nHalo

    fluid(ii,jj)%compU(:,nY+1:ubY,:) = fluid(ii,jj+1)%compU(:,1:3,:)

end subroutine interfaceBlockNorth


subroutine interfaceBlockSouth(ii, jj)
    use numberPrecisionMod
    use meshMod, only : nHalo
    use generalDataMod, only : fluid
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbY
    integer :: nY

    lbY = -nHalo + 1
    nY = fluid(ii,jj-1)%nY

    fluid(ii,jj)%compU(:,lbY:0,:) = fluid(ii,jj-1)%compU(:,nY+lbY:nY,:)

end subroutine interfaceBlockSouth


subroutine noSlipBlockEast(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubX
    integer :: nX

    nX = fluid(ii,jj)%nX
    ubX = nX + nHalo

    fluid(ii,jj)%compU(nX+1:ubX,:,idRho) =  fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idRho)
    fluid(ii,jj)%compU(nX+1:ubX,:,idU)   = -fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idU) * 1d0
    fluid(ii,jj)%compU(nX+1:ubX,:,idV)   = -fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idV)
    fluid(ii,jj)%compU(nX+1:ubX,:,idE)   =  fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idE)

end subroutine noSlipBlockEast


subroutine noSlipBlockWest(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbX
    integer :: nX

    nX = fluid(ii,jj)%nX
    lbX = -nHalo + 1

    fluid(ii,jj)%compU(0:lbX:-1,:,idRho) =  fluid(ii,jj)%compU(1:nHalo,:,idRho)
    fluid(ii,jj)%compU(0:lbX:-1,:,idU)   = -fluid(ii,jj)%compU(1:nHalo,:,idU) * 1d0
    fluid(ii,jj)%compU(0:lbX:-1,:,idV)   = -fluid(ii,jj)%compU(1:nHalo,:,idV)
    fluid(ii,jj)%compU(0:lbX:-1,:,idE)   =  fluid(ii,jj)%compU(1:nHalo,:,idE)

end subroutine noSlipBlockWest


subroutine noSlipBlockNorth(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubY
    integer :: nY

    nY = fluid(ii,jj)%nY
    ubY = nY + nHalo

    fluid(ii,jj)%compU(:,nY+1:ubY,idRho) =  fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idRho)
    fluid(ii,jj)%compU(:,nY+1:ubY,idU)   = -fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idU)
    fluid(ii,jj)%compU(:,nY+1:ubY,idV)   = -fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idV) * 1d0
    fluid(ii,jj)%compU(:,nY+1:ubY,idE)   =  fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idE)

end subroutine noSlipBlockNorth


subroutine noSlipBlockSouth(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbY
    integer :: nY

    nY = fluid(ii,jj)%nY
    lbY = -nHalo + 1

    fluid(ii,jj)%compU(:,0:lbY:-1,idRho) =  fluid(ii,jj)%compU(:,1:nHalo,idRho)
    fluid(ii,jj)%compU(:,0:lbY:-1,idU)   = -fluid(ii,jj)%compU(:,1:nHalo,idU)
    fluid(ii,jj)%compU(:,0:lbY:-1,idV)   = -fluid(ii,jj)%compU(:,1:nHalo,idV) * 1d0
    fluid(ii,jj)%compU(:,0:lbY:-1,idE)   =  fluid(ii,jj)%compU(:,1:nHalo,idE)

end subroutine noSlipBlockSouth


subroutine zeroGradientBlockEast(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubX
    integer :: nX

    nX = fluid(ii,jj)%nX
    ubX = nX + nHalo

    ! fluid(ii,jj)%compU(nX+1:ubX,:,:) =  spread(fluid(ii,jj)%compU(nX,:,:),1,nHalo)

    fluid(ii,jj)%compU(nX+1:ubX,:,idRho) = fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idRho)
    fluid(ii,jj)%compU(nX+1:ubX,:,idU)   = fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idU)
    fluid(ii,jj)%compU(nX+1:ubX,:,idV)   = fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idV)
    fluid(ii,jj)%compU(nX+1:ubX,:,idE)   = fluid(ii,jj)%compU(nX:nX-nHalo+1:-1,:,idE)

end subroutine zeroGradientBlockEast


subroutine zeroGradientBlockWest(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbX
    integer :: nX

    nX = fluid(ii,jj)%nX
    lbX = -nHalo + 1

    ! fluid(ii,jj)%compU(lbX:0,:,:) = spread(fluid(ii,jj)%compU(1,:,:), 1, nHalo)

    fluid(ii,jj)%compU(0:lbX:-1,:,idRho) = fluid(ii,jj)%compU(1:nHalo,:,idRho)
    fluid(ii,jj)%compU(0:lbX:-1,:,idU)   = fluid(ii,jj)%compU(1:nHalo,:,idU)
    fluid(ii,jj)%compU(0:lbX:-1,:,idV)   = fluid(ii,jj)%compU(1:nHalo,:,idV)
    fluid(ii,jj)%compU(0:lbX:-1,:,idE)   = fluid(ii,jj)%compU(1:nHalo,:,idE)


end subroutine zeroGradientBlockWest

subroutine zeroGradientBlockSouth(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbY
    integer :: nY

    nY = fluid(ii,jj)%nY
    lbY = -nHalo + 1

    fluid(ii,jj)%compU(:,0:lbY:-1,idRho) = fluid(ii,jj)%compU(:,1:nHalo,idRho)
    fluid(ii,jj)%compU(:,0:lbY:-1,idU)   = fluid(ii,jj)%compU(:,1:nHalo,idU)
    fluid(ii,jj)%compU(:,0:lbY:-1,idV)   = fluid(ii,jj)%compU(:,1:nHalo,idV)
    fluid(ii,jj)%compU(:,0:lbY:-1,idE)   = fluid(ii,jj)%compU(:,1:nHalo,idE)

end subroutine zeroGradientBlockSouth


subroutine zeroGradientBlockNorth(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    implicit none
    integer, intent(in) :: ii, jj
    integer :: ubY
    integer :: nY

    nY = fluid(ii,jj)%nY
    ubY = nY + nHalo

    fluid(ii,jj)%compU(:,nY+1:ubY,idRho) = fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idRho)
    fluid(ii,jj)%compU(:,nY+1:ubY,idU)   = fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idU)
    fluid(ii,jj)%compU(:,nY+1:ubY,idV)   = fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idV)
    fluid(ii,jj)%compU(:,nY+1:ubY,idE)   = fluid(ii,jj)%compU(:,nY:nY-nHalo+1:-1,idE)

end subroutine zeroGradientBlockNorth


subroutine dirichletBlockWest(ii, jj)
    use numberPrecisionMod
    use generalDataMod
    use meshMod, only : nHalo
    use proplemProperties
    implicit none
    integer, intent(in) :: ii, jj
    integer :: lbX, n = 0
    real(dp) :: energy, en, ke
    real(dp) :: u, v, rho, p

    lbX = -nHalo + 1

    ! rho = (shockDown(1)*2d0 - fluid(ii,jj)%prim(n+1,3,idRho))
    ! u   = (shockDown(2)*2d0 - fluid(ii,jj)%prim(n+1,3,idU  ))
    ! v   = (shockDown(3)*2d0 - fluid(ii,jj)%prim(n+1,3,idV  ))
    ! p   = (shockDown(4)*2d0 - fluid(ii,jj)%prim(n+1,3,idP  ))

    rho = shockDown(1)
    u   = shockDown(2)
    v   = shockDown(3)
    p   = shockDown(4)

    fluid(ii,jj)%compU(lbX:n,:,idRho) = rho
    fluid(ii,jj)%compU(lbX:n,:,idU)   = u
    fluid(ii,jj)%compU(lbX:n,:,idV)   = v

    en = p/(gammaFluid - 1.0d0)
    en = en/rho
    ke = 0.5d0 * sqrt(u**2.0d0 + v**2.0d0)

    energy =  ke + en
    energy = energy * rho

    fluid(ii,jj)%compU(lbX:n,:,idE) = energy


end subroutine dirichletBlockWest





! fluid(ii,jj)%compU(lbX:0,:,:)    = spread((fluid(ii,jj)%compU(1,:,:)), 1, 3)

