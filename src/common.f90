module numberPrecisionMod
    use iso_fortran_env, only : dp => real64
end module numberPrecisionMod


module blocksMod
    use numberPrecisionMod, only : dp
    implicit none
    private
    public :: blocks, nComp, nPrim, idRho, idU, idV, idE &
          & , idP

    integer, parameter :: nComp = 4
    integer, parameter :: nPrim = 5

    integer, parameter :: idRho = 1
    integer, parameter :: idU   = 2
    integer, parameter :: idV   = 3
    integer, parameter :: idE   = 4
    integer, parameter :: idP   = 5
    ! integer, parameter :: idA   = 6

    type blocks
        real(dp) :: startX, endX, startY, endY
        integer :: nX, nY   ! number of cells in a block
        real(dp) :: dx, dy   ! mesh size
        real(dp), allocatable, dimension(:,:) :: cx, cy   ! cell centers
        real(dp), allocatable, dimension(:,:) :: vx, vy   ! vertices

        real(dp), allocatable, dimension(:,:) :: a   ! sonic speed

        real(dp), allocatable, dimension(:,:,:) :: prim, compU, compF, compG
        real(dp), allocatable, dimension(:,:,:) :: L, fHat, u0, fUX, fDX, fUY, fDY

        character(5) :: zone    ! fluid or solid
    end type blocks

end module blocksMod


! this and geometryMod go together. Had to move domainAspectRatio here
module meshMod
    use numberPrecisionMod, only : dp
    implicit none
    private
    public :: nCellDomainX, nCellDomainY, nHalo, nDim, domainAspectRatio

    real(dp), parameter :: domainAspectRatio = 1.0d0

    integer, parameter :: nDim = 2      ! spatial dimenstions
    integer, parameter :: nCellDomainX = 400   ! number of cells, not nodes in full domain
    integer, parameter :: nCellDomainY = ceiling(nCellDomainX/domainAspectRatio)

    integer, parameter :: nHalo = 3   ! on each side of physical dimensions

end module meshMod


module geometryMod
    use numberPrecisionMod, only : dp
    use meshMod
    implicit none

    real(dp), parameter :: domainStartX = 0.0d0
    real(dp), parameter :: domainStartY = 0.0d0
    real(dp), parameter :: domainLength = 1d0
    real(dp), parameter :: domainWidth  = domainLength/domainAspectRatio

    real(dp), parameter :: approxDx = domainLength/nCellDomainX
    real(dp), parameter :: approxDy = domainWidth /nCellDomainY

    real(dp), parameter :: obstacleLength = domainLength/8d0
    real(dp), parameter :: obstacleWidth  = domainWidth/8d0

    integer, parameter :: nBlockX = 3
    integer, parameter :: nBlockY = 3

end module geometryMod


module generalDataMod
    use numberPrecisionMod, only : dp
    use geometryMod, only : nBlockX, nBlockY
    use blocksMod
    implicit none

    type(blocks), dimension(nBlockX, nBlockY) :: fluid

end module generalDataMod




