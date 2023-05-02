#include "../afivo/src/cpp_macros.h"
!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods
  use m_gas
  use m_streamer

  implicit none
  private

  ! Public methods
  public :: user_initialize
  !real(dp) :: density_ratio = 0.8_dp
  !real(dp) :: sphere_radius = 0.2_dp
  !real(dp) :: sphere_center(NDIM) = 0.5_dp

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree
    logical                   :: lower_density_region = .true.
    call CFG_add_get(cfg, "lower_density_region", lower_density_region, &
            "Whether to have a lower density region in the streamer path")

    !call CFG_add_get(cfg, "density_ratio", density_ratio, &
    !     "Density ratio (> 0)")
    !call CFG_add_get(cfg, "sphere_center", sphere_center, &
    !     "Center (relative to domain) of sphere")
    !call CFG_add_get(cfg, "sphere_radius", sphere_radius, &
    !     "Radius (relative to domain) of sphere")
    !call CFG_add_get(cfg, "density_ratio_outside_sphere", density_ratio_outside_sphere, &
    !     "Whether density ratio is outside sphere")

    user_log_variables => average_temperature
    !if (lower_density_region) then
    !    user_gas_density => gas_density_sphere
    !end if
  end subroutine user_initialize

  subroutine average_temperature(tree, n_vars, var_names, var_values)
       use m_streamer
       use m_gas
       type(af_t), intent(in)                 :: tree
       integer, intent(out)                   :: n_vars
       character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
       real(dp), intent(inout)                :: var_values(user_max_log_vars)
       real(dp)                               :: r(2)
       type(af_loc_t)                         :: loc_maxtemp
       
       
       n_vars = 5
       var_names(1) = 'average_temperature'
       var_names(2) = 'x'
       var_names(3) = 'y'
       var_names(4) = 'max_temperature'
       var_names(5) = 'J.E'
       
       call af_tree_sum_cc(tree, gas_prim_vars(i_e+1), var_values(1))
       call af_tree_sum_cc(tree, i_power_density, var_values(5))
       call af_tree_max_cc(tree, gas_prim_vars(i_e+1), var_values(4), loc_maxtemp)
       r= af_r_loc(tree, loc_maxtemp)
       var_values(2) = r(1)
       var_values(3) = r(2)
  end subroutine average_temperature
  
  
  real(dp) function gas_density_sphere(box, IJK)
    use m_gas

    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    integer                    :: nc
    real(dp)                   :: rr(NDIM), inv_weight


    nc = box%n_cell
    inv_weight = 1/gas_molecular_weight
    rr = af_r_cc(box, [IJK])
    if (rr(1) < 0.05*ST_DOMAIN_LEN(1)) then
            gas_density_sphere = gas_number_density 
    end if
  end function gas_density_sphere

end module m_user
