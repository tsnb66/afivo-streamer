#include "../afivo/src/cpp_macros.h"
!> Template for user code, this one simply uses the default routines
module m_user
  use m_af_all
  use m_config
  use m_user_methods

  implicit none
  private

  ! Relative coordinates in the range [0, 1]
  real(dp) :: rmin_rel(NDIM), rmax_rel(NDIM)
  ! Density of the patch
  real(dp) :: patch_density

  ! Public methods
  public :: user_initialize

contains

  subroutine user_initialize(cfg, tree)
    type(CFG_t), intent(inout) :: cfg
    type(af_t), intent(inout) :: tree

    user_initial_conditions => cuboid_patch
    call CFG_add_get(cfg, "cuboid_patch_rmin", rmin_rel, &
         "Min coordinates for the patch")
    call CFG_add_get(cfg, "cuboid_patch_rmax", rmax_rel, &
         "Max coordinates for the patch")
    call CFG_add_get(cfg, "cuboid_patch_density", patch_density, &
         "Neutral seed density inside the patch")
    ! r_min = [12e-3_dp, 12e-3_dp]
    ! r_max = [20e-3_dp, 20e-3_dp]
    ! my_density = 1e15_dp

  end subroutine user_initialize

  subroutine cuboid_patch(box)
    use m_streamer
    type(box_t), intent(inout) :: box
    ! integer                    :: i, j
    integer                    :: IJK, nc
    real(dp)                   :: r(NDIM)
    real(dp)                   :: rmin_abs(NDIM), rmax_abs(NDIM)

    ! Converting from relative to absolute coordinates
    rmin_abs = rmin_rel*ST_domain_len + ST_domain_origin
    rmax_abs = rmax_rel*ST_domain_len + ST_domain_origin

    nc = box%n_cell
    do KJI_DO(1, nc)
    !do j = 1, box%n_cell
    !   do i = 1, box%n_cell
          !r = af_r_cc(box, [i, j])
          r = af_r_cc(box, [IJK])

          if (all(r >= rmin_abs .and. r <= rmax_abs)) then
             !box%cc(i, j, i_electron) = box%cc(i, j, i_electron) + my_density
             !box%cc(i, j, i_1pos_ion) = box%cc(i, j, i_electron) + my_density
             box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + patch_density
             box%cc(IJK, i_1pos_ion) = box%cc(IJK, i_1pos_ion) + patch_density
          end if
       !end do
    !end do
    end do; CLOSE_DO
  end subroutine cuboid_patch

end module m_user
