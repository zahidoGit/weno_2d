subroutine deallocateAll()
    use generalDataMod
    implicit none
    integer :: ii, jj

    do jj = 1, nBlockY
        do ii = 1, nBlockX

            deallocate(fluid(ii,jj)%vx)
            deallocate(fluid(ii,jj)%vy)
            deallocate(fluid(ii,jj)%cx)
            deallocate(fluid(ii,jj)%cy)

            deallocate(fluid(ii,jj)%prim)
            deallocate(fluid(ii,jj)%compU)
            deallocate(fluid(ii,jj)%compF)
            deallocate(fluid(ii,jj)%compG)

            deallocate(fluid(ii,jj)%L)
            deallocate(fluid(ii,jj)%fHat)
            deallocate(fluid(ii,jj)%u0)
            deallocate(fluid(ii,jj)%fUX)
            deallocate(fluid(ii,jj)%fDX)
            deallocate(fluid(ii,jj)%fUY)
            deallocate(fluid(ii,jj)%fDY)

            deallocate(fluid(ii,jj)%a)

        end do
    end do

end subroutine
