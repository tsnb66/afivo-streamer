!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_log_variables => average_temperature
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
       
       
       n_vars = 4
       var_names(1) = 'average_temperature'
       var_names(2) = 'x'
       var_names(3) = 'y'
       var_names(4) = 'max_temperature'
       
       call af_tree_sum_cc(tree, gas_prim_vars(i_e+1), var_values(1))
       call af_tree_max_cc(tree, gas_prim_vars(i_e+1), var_values(4), loc_maxtemp)
       r= af_r_loc(tree, loc_maxtemp)
       var_values(2) = r(1)
       var_values(3) = r(2)
  end subroutine average_temperature
  
  

end module m_user
