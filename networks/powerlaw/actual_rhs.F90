module actual_rhs_module

  use amrex_constants_module
  use network
  use burn_type_module
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine actual_rhs_init()

    implicit none

  end subroutine actual_rhs_init



  subroutine actual_rhs(state, ydot)

    use extern_probin_module, only: f_act, T_burn_ref, rho_burn_ref, rtilde, nu
    use temperature_integration_module, only: temperature_rhs, temperature_jac

    implicit none

    type (burn_t), intent(in)    :: state
    real(rt)        , intent(inout) :: ydot(neqs)
    real(rt)         :: xfueltmp
    real(rt)         :: dens, temp, rate, y(nspec)

    !$gpu

    ydot = ZERO

    xfueltmp = max(state % xn(ifuel_), ZERO)
    dens     = state % rho
    temp     = state % T

    ! Rate is expressed in mass fraction form

    if (temp < f_act * T_burn_ref) then
       rate = ZERO
    else
       rate = rtilde * (dens/rho_burn_ref) * xfueltmp**2 * (temp/T_burn_ref)**nu
    endif

    ydot(ifuel_)  = -rate
    ydot(iash_)   =  rate
    ydot(iinert_) = 0.0_rt

    ! Convert back to molar form

    ydot(1:nspec) = ydot(1:nspec) * aion_inv(1:nspec)

    call ener_gener_rate(ydot(1:nspec), ydot(net_ienuc))

    call temperature_rhs(state, ydot)

  end subroutine actual_rhs



  ! At present the analytical Jacobian is not implemented.

  subroutine actual_jac(state, jac)

    implicit none

    type (burn_t), intent(in) :: state
    real(rt)        , intent(inout) :: jac(njrows, njcols)

    !$gpu

    jac(:,:) = ZERO

  end subroutine actual_jac



  ! Computes the instantaneous energy generation rate

  subroutine ener_gener_rate(dydt, enuc)

    use network

    implicit none

    real(rt)         :: dydt(nspec), enuc

    !$gpu

    ! This is basically e = m c**2

    enuc = sum(dydt(:) * aion(1:nspec) * ebin(1:nspec))

  end subroutine ener_gener_rate

end module actual_rhs_module
