module fluxVectorMod
    use numberPrecisionMod
    use generalDataMod
    implicit none

    contains
    pure function getF(p) result(f)
        implicit none
        real(dp), dimension(nComp) :: f
        real(dp), dimension(:), intent(in) :: p
        real(dp) :: ke, energy

        f(idRho) = p(idRho)*p(idU)
        f(idU)   = f(idRho)*p(idU) + p(idP)
        f(idV)   = f(idRho)*p(idV)
        ke       = 0.5d0 * (p(idU)**2 + p(idV)**2)
        energy   = (ke + p(idE)) * p(idrho)
        f(idE)   = p(idU)*(energy + p(idP))

    end function getF


    pure function getG(p) result(g)
        implicit none
        real(dp), dimension(nComp) :: g
        real(dp), dimension(:), intent(in) :: p
        real(dp) :: ke, energy

        g(idRho) = p(idRho)*p(idV)
        g(idU)   = g(idRho)*p(idU)
        g(idV)   = g(idRho)*p(idV) + p(idP)
        ke       = 0.5d0 * (p(idU)**2 + p(idV)**2)
        energy   = (ke + p(idE)) * p(idrho)
        g(idE)   = p(idV)*(energy + p(idP))

    end function getG


end module fluxVectorMod
