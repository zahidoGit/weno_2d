module primitiveMod
    use numberPrecisionMod
    use generalDataMod
    use proplemProperties, only : gammaFluid

    contains
    pure function getPrimitive(u) result(p)
        implicit none
        real(dp), dimension(nPrim) :: p
        real(dp), dimension(:), intent(in) :: u
        real(dp) :: en, ke

        p(idRho) = u(idRho)
        p(idU)   = u(idU)/u(idRho)
        p(idV)   = u(idV)/u(idRho)

        ke       = 0.5d0 * (p(idU)**2 + p(idV)**2)
        en       = u(idE)/p(idRho) - ke

        p(idE)   = en
        p(idP)   = en*p(idRho)*(gammaFluid - 1.0d0)

    end function getPrimitive


    pure function getSonicSpeed(p) result(a)
    implicit none
    real(dp), dimension(:), intent(in) :: p
    real(dp) :: a

    a = sqrt( abs(gammaFluid * p(idP)/p(idRho)) )

    end function getSonicSpeed

end module primitiveMod


