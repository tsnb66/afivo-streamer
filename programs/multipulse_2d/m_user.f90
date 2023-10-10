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

  real(dp) :: density_ratio       = 0.8_dp
  real(dp) :: shock_width         = 0.01_dp
  real(dp) :: line_coeff(NDIM+1)  = 0.0_dp
  real(dp) :: sphere_center(NDIM) = 0.5_dp
  real(dp) :: sphere_radius       = 0.1_dp

  ! Whether density ratio is inside sphere
  logical, protected, public :: density_ratio_inside_sphere = .false.

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout)  :: tree
    character(len=100)         :: gradient_type
    logical                    :: use_gradient = .false.

    user_log_variables => gas_heating_postProc
    call CFG_add_get(cfg, "use_gradient", use_gradient, &
      "Whether to use a gas density gradient")
    if (use_gradient) then
      gradient_type = "line"
      call CFG_add_get(cfg, "gradient_type", gradient_type, &
           "What type of gas gradient to use (line, sphere)")

      select case (gradient_type)
      case ("line")
         user_gas_density => gas_density_line
      case ("sphere")
         user_gas_density => gas_density_sphere
      case default
         error stop "Unknown gradient_type"
      end select

      call CFG_add_get(cfg, "density_ratio", density_ratio, &
           "Density ratio (<= 1)")
      call CFG_add_get(cfg, "shock_width", shock_width, &
           "Shock width (relative to domain size)")
      call CFG_add_get(cfg, "line_coeff", line_coeff, &
           "Coefficients a, b, c of a line a + bx + cy = 0")
      call CFG_add_get(cfg, "sphere_center", sphere_center, &
           "Center (relative to domain) of sphere")
      call CFG_add_get(cfg, "sphere_radius", sphere_radius, &
           "Radius (relative to domain) of sphere")
      call CFG_add_get(cfg, "density_ratio_inside_sphere", density_ratio_inside_sphere, &
           "Whether density ratio is inside sphere")
    end if
         
  end subroutine user_initialize

  !> Gas density is different on two sides of a line
  real(dp) function gas_density_line(box, IJK)
    use m_gas
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    real(dp)                :: r_rel(NDIM), q, tmp

    ! Get relative coordinates in the range [0, 1]
    r_rel = (af_r_cc(box, [IJK]) - ST_domain_origin)/ST_domain_len
    q = line_coeff(1) + sum(line_coeff(2:) * r_rel)

    if (q < -shock_width) then
       gas_density_line = gas_number_density
    else if (q > shock_width) then
       gas_density_line = gas_number_density * density_ratio
    else
       ! Linear interpolation
       tmp = (q + shock_width) / (2 * shock_width)
       gas_density_line = gas_number_density * (1 + (density_ratio-1) * tmp)
    end if
  end function gas_density_line

  !> Gas density is different inside and outside sphere
  real(dp) function gas_density_sphere(box, IJK)
    use m_gas
    type(box_t), intent(in) :: box
    integer, intent(in)     :: IJK
    real(dp)                :: q, tmp

    ! Get relative coordinates with regard to sphere center
    q = norm2((af_r_cc(box, [IJK]) - ST_domain_origin)/ST_domain_len &
         - sphere_center)

    if (q < sphere_radius - shock_width) then
       gas_density_sphere = gas_number_density
       if (density_ratio_inside_sphere) then
          gas_density_sphere = gas_number_density * density_ratio
       end if
    else if (q > sphere_radius + shock_width) then
       gas_density_sphere = gas_number_density * density_ratio
       if (density_ratio_inside_sphere) then
          gas_density_sphere = gas_number_density
       end if
    else
       ! Linear interpolation
       tmp = (q - sphere_radius + shock_width) / (2 * shock_width)
       if (density_ratio_inside_sphere) then
          tmp = (sphere_radius + shock_width - q) / (2 * shock_width)
       end if
       gas_density_sphere = gas_number_density * (1 + (density_ratio-1) * tmp)
    end if
  end function gas_density_sphere

  subroutine gas_heating_postProc(tree, n_vars, var_names, var_values)
       use m_streamer
       use m_gas
       use m_units_constants
       type(af_t), intent(in)                 :: tree
       integer, intent(out)                   :: n_vars
       character(len=name_len), intent(inout) :: var_names(user_max_log_vars)
       real(dp), intent(inout)                :: var_values(user_max_log_vars)
       real(dp)                               :: r(2), avg_temp, vol_area
       type(af_loc_t)                         :: loc_maxtemp
       type(af_loc_t)                         :: loc_maxpowerdens


       n_vars = 7
       var_names(1) = 'average_temperature'
       var_names(2) = 'x'
       var_names(3) = 'y'
       var_names(4) = 'max_temperature'
       var_names(5) = 'x'
       var_names(6) = 'y'
       var_names(7) = 'max_power_density'
       vol_area = (UC_pi*ST_domain_len(1)**2*ST_domain_len(NDIM))
       print *, "Volume: ", vol_area

       call af_tree_sum_cc(tree, gas_prim_vars(i_e+1), avg_temp)
       var_values(1) = avg_temp/vol_area
       call af_tree_max_cc(tree, gas_prim_vars(i_e+1), var_values(4), loc_maxtemp)
       r= af_r_loc(tree, loc_maxtemp)
       var_values(2) = r(1)
       var_values(3) = r(2)
       call af_tree_max_cc(tree, i_power_density, var_values(7), loc_maxpowerdens)
       r= af_r_loc(tree, loc_maxpowerdens)
       var_values(5) = r(1)
       var_values(6) = r(2)
  end subroutine gas_heating_postProc

end module m_user
