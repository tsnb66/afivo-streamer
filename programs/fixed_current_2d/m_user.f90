!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  real(dp) :: desired_current = -50e-3_dp
  real(dp) :: relaxation_time = 2e-9_dp

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_field_amplitude => my_field_amplitude

  end subroutine user_initialize

  real(dp) function my_field_amplitude(tree, time)
    use m_streamer
    use m_field
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: time

    real(dp)       :: resistance, goal_voltage, voltage_change, dt
    integer, save  :: counter   = 0
    real(dp), save :: prev_time = 0

    dt = time - prev_time

    if (time < 1e-9_dp) then
       my_field_amplitude = 2e6_dp
    else
       ! Estimate resistance
       resistance = current_voltage/ST_global_current

       goal_voltage = desired_current * resistance
       voltage_change = (goal_voltage - current_voltage) * dt / relaxation_time
       my_field_amplitude = (current_voltage + voltage_change) / (-ST_domain_len(NDIM))

       counter = counter + 1
       if (modulo(counter, 100) == 0) then
          print *, time, ST_global_current, my_field_amplitude, resistance
       end if
    end if

    prev_time = time

  end function my_field_amplitude


end module m_user
