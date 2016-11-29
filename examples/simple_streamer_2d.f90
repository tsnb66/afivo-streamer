!> \example
!> A simplified model for ionization waves and/or streamers in 2D
program simple_streamer_2d

  use m_a2_types
  use m_a2_core
  use m_a2_ghostcell
  use m_a2_utils
  use m_a2_restrict
  use m_a2_multigrid
  use m_a2_output
  use m_write_silo
  use m_a2_prolong, onLy: a2_prolong_copy_from

  implicit none

  integer            :: i, n
  character(len=100) :: fname
  type(a2_t)         :: tree ! This contains the full grid information
  type(mg2_t)        :: mg   ! Multigrid option struct
  type(ref_info_t)   :: ref_info

  ! Indices of cell-centered variables
  integer, parameter :: n_var_cell = 7 ! Number of variables
  integer, parameter :: i_elec     = 1 ! Electron density
  integer, parameter :: i_pion     = 2 ! Positive ion density
  integer, parameter :: i_elec_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pion_old = 4 ! For time-stepping scheme
  integer, parameter :: i_phi      = 5 ! Electrical potential
  integer, parameter :: i_fld      = 6 ! Electric field norm
  integer, parameter :: i_rhs      = 7 ! Source term Poisson

  ! Indices of face-centered variables **
  integer, parameter :: n_var_face = 2 ! Number of variables
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  ! Names of the cell-centered variables
  character(len=10) :: ST_cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs"]

  ! Simulation parameters
  real(dp), parameter :: end_time      = 8e-9_dp
  real(dp), parameter :: dt_output     = 20e-11_dp
  real(dp), parameter :: dt_max        = 20e-11_dp
  integer, parameter  :: ref_per_steps = 1
  integer, parameter  :: box_size      = 8

  ! Physical parameters
  real(dp), parameter :: applied_field = -0.8e7_dp
  real(dp), parameter :: mobility      = 0.03_dp
  real(dp), parameter :: diffusion_c   = 0.2_dp

  ! Computational domain
  real(dp), parameter :: domain_length = 10e-3_dp
  real(dp), parameter :: refine_max_dx = 1e-3_dp
  real(dp), parameter :: refine_min_dx = 1e-9_dp

  ! Settings for the initial conditions
  real(dp), parameter :: init_density = 1e15_dp
  real(dp), parameter :: init_y_min   = 8.0e-3_dp
  real(dp), parameter :: init_y_max   = 9.0e-3_dp

  ! Simulation variables
  real(dp) :: dt
  real(dp) :: time
  integer  :: output_count

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs

  ! Routines to use for ...
  mg%sides_bc => sides_bc_pot ! Filling ghost cell on physical boundaries

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  output_count = 0 ! Number of output files written
  time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_init_cond)
     call compute_fld(tree, .false.)
     call a2_adjust_refinement(tree, refinement_routine, ref_info)
     if (ref_info%n_add == 0 .and. ref_info%n_rm == 0) exit
  end do

  call a2_print_info(tree)

  do
     ! Get a new time step, which is at most dt_max
     call a2_reduction(tree, max_dt, get_min, dt_max, dt)

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?", dt
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_count * dt_output <= time) then
        output_count = output_count + 1
        write(fname, "(A,I6.6)") "simple_streamer_2d_", output_count
        call a2_write_silo(tree, fname, output_count, time, dir="output")
     end if

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ref_per_steps
        time = time + dt

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! First calculate fluxes
           call a2_loop_boxes(tree, fluxes_koren, leaves_only=.true.)
           call a2_consistent_fluxes(tree, [f_elec])

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [dt], &
                leaves_only=.true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)

           ! Fill ghost cells
           call a2_gc_tree(tree, i_elec, a2_gc_interp_lim, a2_bc_dirichlet_zero)

           ! Compute new field on first iteration
           if (i == 1) call compute_fld(tree, .true.)
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens)

        ! Compute field with new density
        call compute_fld(tree, .true.)
     end do

     call a2_restrict_tree(tree, i_pion)
     call a2_gc_tree(tree, i_pion, a2_gc_interp_lim, a2_bc_dirichlet_zero)
     call a2_adjust_refinement(tree, refinement_routine, ref_info, 4)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_fld(tree, .true.)
     end if
  end do

  call a2_destroy(tree)

contains

  ! Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = domain_length / box_size

    ! Initialize tree
    call a2_init(tree, box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=2, n_boxes = n_boxes_init, cc_names=ST_cc_names)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...

    nb_list(a2_neighb_lowy, id)  = -1 ! physical boundary
    nb_list(a2_neighb_highy, id) = -1 ! idem
    nb_list(a2_neighb_lowx, id)  = id ! periodic boundary
    nb_list(a2_neighb_highx, id) = id ! idem

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! This routine sets the refinement flag for boxes(id)
  subroutine refinement_routine(box, cell_flags)
    type(box2_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: dx, dens, fld, adx, xy(2)

    nc      = box%n_cell
    dx      = box%dr

    do j = 1, nc
       do i = 1, nc
          xy   = a2_r_cc(box, [i,j])
          dens = box%cc(i, j, i_elec)
          fld = box%cc(i, j, i_fld)
          adx = get_alpha(fld) * dx

          if (dens > 1e0_dp .and. adx > 0.8_dp) then
             cell_flags(i, j) = af_do_ref
          else if (dx < 1.25e-5_dp .and. adx < 0.1_dp) then
             cell_flags(i, j) = af_rm_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    else if (dx > 0.5_dp * refine_max_dx) then
       where(cell_flags == af_rm_ref) cell_flags = af_keep_ref
    end if

  end subroutine refinement_routine

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2), normal_rands(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          if (xy(2) > init_y_min .and. xy(2) < init_y_max) then
             ! Approximate Poisson distribution with normal distribution
             normal_rands = two_normals(box%dr**3 * init_density, &
                   sqrt(box%dr**3 * init_density))
             ! Prevent negative numbers
             box%cc(i, j, i_elec) = abs(normal_rands(1)) / box%dr**3
          else
             box%cc(i, j, i_elec) = 0
          end if
       end do
    end do

    box%cc(:, :, i_pion) = box%cc(:, :, i_elec)
    box%cc(:, :, i_phi)  = 0 ! Inital potential set to zero

  end subroutine set_init_cond

  !> Return two normal random variates
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(mean, sigma) result(rands)
    real(dp), intent(in) :: mean, sigma
    real(dp) :: rands(2), sum_sq

    do
       call random_number(rands)
       rands = 2 * rands - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    rands = mean + rands * sigma
  end function two_normals

  real(dp) function get_min(a, b)
    real(dp), intent(in) :: a, b
    get_min = min(a, b)
  end function get_min

  ! Get maximum time step based on e.g. CFL criteria
  real(dp) function max_dt(box)
    type(box2_t), intent(in) :: box
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    integer :: i, j, nc
    real(dp)               :: fld(2)
    real(dp)               :: dt_cfl, dt_dif, dt_drt

    nc = box%n_cell
    dt_cfl = dt_max

    do j = 1, nc
       do i = 1, nc
          fld(1) = 0.5_dp * (box%fx(i, j, f_fld) + box%fx(i+1, j, f_fld))
          fld(2) = 0.5_dp * (box%fy(i, j, f_fld) + box%fy(i, j+1, f_fld))

          ! The 0.5 is here because of the explicit trapezoidal rule
          dt_cfl = min(dt_cfl, 0.5_dp / sum(abs(fld * mobility) / box%dr))
       end do
    end do

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * mobility * &
         maxval(box%cc(1:nc, 1:nc, i_elec)))

    ! Diffusion condition
    dt_dif = 0.25_dp * box%dr**2 / max(diffusion_c, epsilon(1.0_dp))

    max_dt = min(dt_cfl, dt_drt, dt_dif)
  end function max_dt

  ! Taken from: Spatially hybrid computations for streamer discharges: II. Fully
  ! 3D simulations, Chao Li, Ute Ebert, Willem Hundsdorfer, J. Comput. Phys.
  ! 231, 1020-1050 (2012), doi:10.1016/j.jcp.2011.07.023
  elemental function get_alpha(fld) result(alpha)
    real(dp), intent(in) :: fld
    real(dp)             :: alpha
    real(dp), parameter  :: c0 = 1.04e1_dp
    real(dp), parameter  :: c1 = 0.601_dp
    real(dp), parameter  :: c2 = 1.86e7_dp

    alpha = exp(c0) * (abs(fld) * 1e-5_dp)**c1 * exp(-c2 / abs(fld))
  end function get_alpha

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, have_guess)
    type(a2_t), intent(inout) :: tree
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = 1.6021766208e-19_dp / 8.8541878176e-12_dp
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_elec) - &
               tree%boxes(id)%cc(:, :, i_pion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Perform an FMG cycle (logicals: store residual, first call)
    call mg2_fas_fmg(tree, mg, .false., have_guess)

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call a2_gc_tree(tree, i_fld, a2_gc_interp, a2_bc_neumann_zero)
  end subroutine compute_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, f_fld) + box%fx(2:nc+1, 1:nc, f_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, f_fld) + box%fy(1:nc, 2:nc+1, f_fld))**2)
  end subroutine fld_from_pot

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp)                    :: inv_dr, tmp, gradp, gradc, gradn
    real(dp)                    :: v_drift
    real(dp)                    :: fld
    real(dp)                    :: gc_data(boxes(id)%n_cell, a2_num_neighbors)
    integer                     :: i, j, nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    call a2_gc2_box(boxes, id, i_elec, a2_gc2_prolong_linear, &
         a2_bc2_neumann_zero, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          fld        = boxes(id)%fx(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)

          if (v_drift < 0.0_dp) then
             if (i == nc+1) then
                tmp = gc_data(j, a2_neighb_highx)
             else
                tmp = boxes(id)%cc(i+1, j, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
          else                  ! v_drift > 0
             if (i == 1) then
                tmp = gc_data(j, a2_neighb_lowx)
             else
                tmp = boxes(id)%cc(i-2, j, i_elec)
             end if
             gradp = boxes(id)%cc(i-1, j, i_elec) - tmp
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i-1, j, i_elec) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_elec) - &
               diffusion_c * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          fld        = boxes(id)%fy(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)

          if (v_drift < 0.0_dp) then
             if (j == nc+1) then
                tmp = gc_data(i, a2_neighb_highy)
             else
                tmp = boxes(id)%cc(i, j+1, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
          else                  ! v_drift > 0
             if (j == 1) then
                tmp = gc_data(i, a2_neighb_lowy)
             else
                tmp = boxes(id)%cc(i, j-2, i_elec)
             end if
             gradp = boxes(id)%cc(i, j-1, i_elec) - tmp
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j-1, i_elec) + koren_mlim(gradc, gradp))
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_elec) - &
               diffusion_c * gradc * inv_dr
       end do
    end do

  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, sflux, fld
    real(dp)                    :: alpha
    integer                     :: i, j, nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    do j = 1, nc
       do i = 1, nc
          fld   = box%cc(i,j, i_fld)
          alpha = get_alpha(fld)
          src   = box%cc(i, j, i_elec) * mobility * abs(fld) * alpha

          sflux = (box%fx(i, j, f_elec) - box%fx(i+1, j, f_elec) + &
               box%fy(i, j, f_elec) - box%fy(i, j+1, f_elec)) * inv_dr

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (src + sflux) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src * dt(1)
       end do
    end do

  end subroutine update_solution

  ! For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_elec)
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pion)
          call a2_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box(tree%boxes, id, i_elec, &
               a2_gc_interp_lim, a2_bc_dirichlet_zero)
          call a2_gc_box(tree%boxes, id, i_pion, &
               a2_gc_interp_lim, a2_bc_dirichlet_zero)
          call a2_gc_box(tree%boxes, id, i_phi, &
               mg2_sides_rb, mg%sides_bc)
       end do
    end do
  end subroutine prolong_to_new_boxes

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(box, nb, iv, bc_type)
    type(box2_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
    case (a2_neighb_lowy)
       ! bc_type = af_bc_dirichlet
       ! box%cc(1:nc, 0, iv) = -applied_field * domain_length
       bc_type = af_bc_neumann
       box%cc(1:nc, 0, iv) = applied_field
    case (a2_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = 0
    case default
       stop "sides_bc_pot: unspecified boundary"
    end select
  end subroutine sides_bc_pot

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

end program
