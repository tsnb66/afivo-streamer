#include "../afivo/src/cpp_macros.h"
!> Module for the coupling between the gas dynamics and the fluid model
module m_coupling
  use m_types
  use m_af_all
  use m_units_constants
  use m_gas
  use m_streamer

  implicit none
  private




  public :: coupling_add_fluid_source
  public :: coupling_update_gas_density


contains
 
  !> Add source terms form the fluid model to the Euler equations
  subroutine coupling_add_fluid_source(tree, dt)
    type(af_t), intent(inout) :: tree
    real(dp), intent(in)      :: dt

    call af_loop_box_arg(tree, add_heating_box, [dt], .true.)
  end subroutine coupling_add_fluid_source

  subroutine add_heating_box(box, dt_vec)
    use m_gas
    use m_lookup_table
    type(box_t), intent(inout) :: box
    real(dp), intent(in)       :: dt_vec(:)
    integer                    :: IJK, nc
    real(dp)                   :: J_dot_E
    !real(dp),parameter         :: gamma = 1.4_dp, gas_const = 8.314_dp
    real(dp)                   :: eta, eta_rt, eta_el, eta_vt
    eta_rt = 0.0_dp
    eta_el = 0.0_dp
    eta_vt = 0.0_dp
    eta = 1.0_dp
    nc = box%n_cell
    do KJI_DO(1, nc)
       ! Compute inner product flux * field over the cell faces
       J_dot_E = 0.5_dp * sum(box%fc(IJK, :, flux_elec) * box%fc(IJK, :, electric_fld))
#if NDIM == 1
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, 1, flux_elec) * box%fc(i+1, 1, electric_fld))
#elif NDIM == 2
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, 1, flux_elec) * box%fc(i+1, j, 1, electric_fld) + &
            box%fc(i, j+1, 2, flux_elec) * box%fc(i, j+1, 2, electric_fld))
#elif NDIM == 3
       J_dot_E = J_dot_E + 0.5_dp * (&
            box%fc(i+1, j, k, 1, flux_elec) * box%fc(i+1, j, k, 1, electric_fld) + &
            box%fc(i, j+1, k, 2, flux_elec) * box%fc(i, j+1, k, 2, electric_fld) + &
            box%fc(i, j, k+1, 3, flux_elec) * box%fc(i, j, k+1, 3, electric_fld))
#endif
      if (effic_table_use) then
        call LT_lin_interp_list(rt_efficiency_field,rt_efficiency_val, &
          box%cc(IJK, electric_fld), eta_rt)
        call LT_lin_interp_list(el_efficiency_field,el_efficiency_val, &
        box%cc(IJK, electric_fld), eta_el)
        call LT_lin_interp_list(vt_efficiency_field,vt_efficiency_val, &
        box%cc(IJK, electric_fld), eta_vt)
        eta = eta_rt + 0.3_dp*eta_el
        if (eta .ge. 1) error stop "Heating efficiency larger than 1"
      end if
      J_dot_E = J_dot_E * UC_elec_charge
      box%cc(IJK, i_vibration_energy) = box%cc(IJK, i_vibration_energy)+ &
      !   (0.5e9_dp*(box%cc(IJK, i_vibration_energy) - 1.0_dp))*dt_vec(1)
        (eta_vt*J_dot_E - box%cc(IJK, i_vibration_energy)/t_vt) * dt_vec(1)
      !box%cc(IJK, i_vibration_energy) = box%cc(IJK, i_vibration_energy)+ &
      !  (J_dot_E*(gamma-1)/(gamma*gas_const)) * dt_vec(1)
       box%cc(IJK, gas_vars(i_e)) = box%cc(IJK, gas_vars(i_e)) + &
           (eta *  J_dot_E+ box%cc(IJK, i_vibration_energy)/t_vt) &
            * dt_vec(1)
    end do; CLOSE_DO
  end subroutine add_heating_box

  !> Update gas number density in the fluid model
  subroutine coupling_update_gas_density(tree)
    type(af_t), intent(inout) :: tree

    call af_loop_box(tree, update_gas_density, .true.)
  end subroutine coupling_update_gas_density

  subroutine update_gas_density(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc
    real(dp)                   :: inv_weight

    nc         = box%n_cell
    inv_weight = 1/gas_molecular_weight

    box%cc(DTIMES(1:nc), i_gas_dens) = &
         box%cc(DTIMES(1:nc), gas_vars(i_rho)) * inv_weight
  end subroutine update_gas_density

end module m_coupling
