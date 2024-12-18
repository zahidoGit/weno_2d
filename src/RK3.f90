module rk3Mod
    use numberPrecisionMod
    use generalDataMod

    ! dt = 1d-4
    ! do while (time .le. 0.2d0)
    ! ! do tt = 1, 500
    !     do jj = 1, nBlockY
    !         do ii = 1, nBlockX
    !
    !             dx = fluid(ii,jj)%dx
    !             dy = fluid(ii,jj)%dy
    !
    !             do pp = 1, nComp
    !                 do j = 1, fluid(ii,jj)%nY
    !                     i = 1
    !                     fp = wenoRconstruct(fluid(ii,jj)%compU(-2+i:2+i, j, pp))
    !                     do i = 1, fluid(ii,jj)%nX
    !                         fn = fp
    !                         fp = wenoRconstruct(fluid(ii,jj)%compU(-2+i+1:2+i+1, j, pp))
    !                         a = fluid(ii,jj)%compU(i+1,j,pp)
    !                         b = fluid(ii,jj)%compU(i,j,pp)
    !                         alpha = max(abs(fluid(ii,jj)%compU(i+1,j,idU)), abs(fluid(ii,jj)%compU(i,j,idU)))
    !                         fluid(ii,jj)%compU(i,j,pp) = fluid(ii,jj)%compU(i,j,pp) + (fp - fn + alpha*(a-b))*dt/dx
    !                     end do
    !                 end do
    !             end do
    !
    !             lbY = lbound(fluid(ii,jj)%prim, 2)
    !             ubY = ubound(fluid(ii,jj)%prim, 2)
    !
    !             lbX = lbound(fluid(ii,jj)%prim, 1)
    !             ubX = ubound(fluid(ii,jj)%prim, 1)
    !
    !             do j = lbY, ubY
    !                 do i = lbx, ubX
    !                     fluid(ii,jj)%prim(i,j,:)  = getPrimitive(fluid(ii,jj)%compU(i,j,:))
    !                     fluid(ii,jj)%compF(i,j,:) = getF(fluid(ii,jj)%prim(i,j,:))
    !                     fluid(ii,jj)%compG(i,j,:) = getG(fluid(ii,jj)%prim(i,j,:))
    !                 end do
    !             end do
    !
    !         end do
    !     end do
    !     call blockBoundaries()
    !     print*, time
    !     time = time + dt
    ! end do
    !





end module rk3Mod
