subroutine saveStructuredTecPlot(t)
    use generalDataMod
    implicit none
    real(dp), intent(in) :: t
    character(len=128) :: fname
    character(len=128) :: stt
    integer :: ii, jj, i, j, lbX, ubX, lbY, ubY

    write(stt,*) t
    fname = trim("output/grid") // trim(adjustl(stt)) // trim(".txt")
    open(unit=10, file=trim(fname), status='replace')

    write(10, '(A)') "X Y Z u v rho p E"

    do jj = 1, nBlockY
        do ii = 1, nBlockX

            if( fluid(ii,jj)%zone .eq. "solid" ) cycle

            lbX = lbound(fluid(ii,jj)%cx, 1)
            ubX = ubound(fluid(ii,jj)%cx, 1)
            lbY = lbound(fluid(ii,jj)%cx, 2)
            ubY = ubound(fluid(ii,jj)%cx, 2)

            ! do j = lbY, ubY
            !     do i = lbX, ubX

            do j = 1, fluid(ii,jj)%nY
                do i = 1, fluid(ii,jj)%nX

                    ! write(10,'(*(ES16.8))') fluid(ii,jj)%cx(i,j), fluid(ii,jj)%cy(i,j), 0.0d0 &
                    write(10,*) fluid(ii,jj)%cx(i,j), fluid(ii,jj)%cy(i,j), 0.0d0 &
                            & , fluid(ii,jj)%prim(i,j,idU) &
                            & , fluid(ii,jj)%prim(i,j,idV) &
                            & , fluid(ii,jj)%prim(i,j,idRho) &
                            & , fluid(ii,jj)%prim(i,j,idP) &
                            & , fluid(ii,jj)%compU(i,j,idE)

                enddo
            enddo

        end do
    end do

    close(10)

end subroutine saveStructuredTecPlot
