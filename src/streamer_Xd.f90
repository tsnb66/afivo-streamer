#include "../afivo/src/cpp_macros_$Dd.h"
!> Program to perform $Dd streamer simulations in Cartesian and cylindrical coordinates
program streamer_$Dd

  use m_a$D_all
  use m_streamer
  use m_field_$Dd
  use m_init_cond_$Dd

  implicit none

  integer                :: i, it
  character(len=ST_slen) :: fname
  logical                :: write_out

  type(CFG_t)            :: cfg ! The configuration for the simulation
  type(a$D_t)            :: tree      ! This contains the full grid information
  type(mg$D_t)           :: mg        ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  real(dp) :: B0_hat(3)

  integer :: output_cnt = 0 ! Number of output files written

  call CFG_update_from_arguments(cfg)

  B0_hat = [1.0_dp, 0.0_dp, 0.0_dp]
  call CFG_add_get(cfg, "magnetic_field_unitvec", B0_hat, &
       "Unitvector of magnetic field (automatically normalized)")
  B0_hat = B0_hat / norm2(B0_hat)

  call ST_initialize(cfg, $D)
  call ST_load_transport_data(cfg)
  call field_initialize(cfg, mg)
  call init_cond_initialize(cfg, $D)

  fname = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_out.cfg"
  call CFG_write(cfg, trim(fname))

  ! Set the multigrid options. First define the variables to use
  mg%i_phi = i_phi
  mg%i_tmp = i_electric_fld
  mg%i_rhs = i_rhs

  ! This automatically handles cylindrical symmetry
  mg%box_op => mg$D_auto_op
  mg%box_gsrb => mg$D_auto_gsrb
  mg%box_corr => mg$D_auto_corr

  ! This routine always needs to be called when using multigrid
  call mg$D_init_mg(mg)

  ! Initialize the tree (which contains all the mesh information)
  if (ST_restart_file == "") then
     call init_tree(tree)

     ! Set up the initial conditions
     do
        call a$D_loop_box(tree, init_cond_set_box)
        call field_compute(tree, mg, .false.)
        call a$D_adjust_refinement(tree, refine_routine, ref_info, 4)
        if (ref_info%n_add + ref_info%n_rm == 0) exit
     end do
  else
     call a$D_read_tree(tree, trim(ST_restart_file))
  end if

  output_cnt = 0 ! Number of output files written
  ST_time    = 0 ! Simulation time (all times are in s)

  call a$D_print_info(tree)

  do it = 1, huge(1)
     ! Get a new time step, which is at most dt_amr
     call a$D_reduction_vec(tree, get_max_dt, get_min, &
          [ST_dt_max, ST_dt_max, ST_dt_max], ST_dt_vec, ST_dt_num_cond)
     ST_dt = minval(ST_dt_vec)

     if (ST_dt < ST_dt_min) then
        print *, "ST_dt getting too small, instability?"
        print *, ST_dt_vec
        exit
     end if

     if (ST_time > ST_end_time) exit

     ! Every ST_dt_output, write output
     if (output_cnt * ST_dt_output <= ST_time + ST_dt) then
        write_out  = .true.
        ST_dt      = output_cnt * ST_dt_output - ST_time
        output_cnt = output_cnt + 1
     else
        write_out = .false.
     end if

     ! Copy previous solution
     call a$D_tree_copy_cc(tree, i_electron, i_electron_old)
     call a$D_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)

     ! Two forward Euler steps over ST_dt
     do i = 1, 2
        ST_time = ST_time + ST_dt

        ! First calculate fluxes
        call a$D_loop_boxes(tree, fluxes_koren, .true.)
        call a$D_consistent_fluxes(tree, [flux_elec])

        ! Update the solution
        call a$D_loop_box_arg(tree, update_solution, [ST_dt], .true.)

        if (i == 1) then
           call field_compute(tree, mg, .true.)
           call a$D_restrict_tree(tree, i_electron)
        end if
     end do

     ST_time = ST_time - ST_dt        ! Go back one time step

     ! Take average (explicit trapezoidal rule)
     call a$D_loop_box(tree, average_density)

     ! Restrict electron density (so coarse block ghost cells are filled correctly)
     call a$D_restrict_tree(tree, i_electron)

     ! Compute field with new density
     call field_compute(tree, mg, .true.)

     if (write_out) then
        call a$D_loop_box(tree, box_set_src_rate)

        ! Fill ghost cells before writing output
        call a$D_gc_tree(tree, i_electron, a$D_gc_interp, a$D_bc_neumann_zero)
        call a$D_gc_tree(tree, i_pos_ion, a$D_gc_interp, a$D_bc_neumann_zero)
        call a$D_gc_tree(tree, i_src_rate, a$D_gc_interp, a$D_bc_neumann_zero)

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", output_cnt
        call a$D_write_silo(tree, fname, output_cnt, &
             ST_time, dir=ST_output_dir, add_vars=add_velocity, add_names=["v1", "v2"])

        if (ST_datfile_write) then
           call a$D_write_tree(tree, fname, ST_output_dir)
        end if

        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_log.txt"
        call write_log_file(tree, fname, output_cnt, ST_output_dir)

        if (ST_lineout_write) then
           write(fname, "(A,I6.6)") trim(ST_simulation_name) // &
                "_line_", output_cnt
           call a$D_write_line(tree, trim(fname), &
                [i_electron, i_pos_ion, i_phi, i_electric_fld], &
                r_min = ST_lineout_rmin(1:$D) * ST_domain_len, &
                r_max = ST_lineout_rmax(1:$D) * ST_domain_len, &
                n_points=ST_lineout_npoints, dir=ST_output_dir)
        end if
     end if

     if (mod(it, ST_refine_per_steps) == 1) then
        ! Restrict the electron and ion densities before refinement
        call a$D_restrict_tree(tree, i_electron)
        call a$D_restrict_tree(tree, i_pos_ion)

        ! Fill ghost cells before refinement
        call a$D_gc_tree(tree, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero)
        call a$D_gc_tree(tree, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero)

        call a$D_adjust_refinement(tree, refine_routine, ref_info, 4)

        if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
           ! For boxes which just have been refined, set data on their children
           call prolong_to_new_boxes(tree, ref_info)

           ! Compute the field on the new mesh
           call field_compute(tree, mg, .true.)
        end if

     end if
  end do

  call a$D_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a$D_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list($D, 1) ! Spatial indices of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a$D_init(tree, ST_box_size, n_var_cell_$Dd, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, coord=af_cyl, &
            cc_names=ST_cc_names(1:n_var_cell_$Dd))
    else
       call a$D_init(tree, ST_box_size, n_var_cell_$Dd, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, &
            cc_names=ST_cc_names(1:n_var_cell_$Dd))
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = 1          ! With index 1,1 ...

    ! Create the base mesh
    call a$D_set_base(tree, 1, ix_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    use m_init_cond_$Dd
    type(box$D_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld($D), vel($D)
    real(dp)                 :: dist, dns
    type(LT2_loc_t) :: loc

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_Ex:i_Ex+$D-1)

       call get_velocity(fld, vel, loc)
       alpha = LT2_get_col_at_loc(ST_td_tbl, i_alpha, loc)
       ! The refinement is based on the ionization length
       adx   = box%dr * alpha
       dns = box%cc(IJK, i_electron)

       ! The refinement is also based on the intensity of the source term.
       ! Here we estimate the curvature of phi (given by dx**2 *
       ! Laplacian(phi))
       cphi = dx2 * abs(box%cc(IJK, i_rhs))

       if (adx / ST_refine_adx + cphi / ST_refine_cphi > 1 .and. &
            dns > ST_refine_min_density) then
          cell_flags(IJK) = af_do_ref
       else if ((adx < 0.125_dp * ST_refine_adx .and. &
            cphi < 0.0625_dp * ST_refine_cphi .or. &
            dns < ST_refine_min_density) &
            .and. 2 * dx < ST_derefine_dx) then
          cell_flags(IJK) = af_rm_ref
       else
          cell_flags(IJK) = af_keep_ref
       end if

       ! Refine around the initial conditions
       if (ST_time < ST_refine_init_time) then
          do n = 1, init_conds%n_cond
             dist = GM_dist_line(a$D_r_cc(box, [IJK]), &
                  init_conds%seed_r0(:, n), &
                  init_conds%seed_r1(:, n), $D)
             if (dist - init_conds%seed_width(n) < 2 * dx &
                  .and. box%dr > ST_refine_init_fac * &
                  init_conds%seed_width(n)) then
                cell_flags(IJK) = af_do_ref
             end if
          end do
       end if

    end do; CLOSE_DO

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  !> Get maximum time step based on e.g. CFL criteria
  function get_max_dt(box, n_cond) result(dt_vec)
    use m_units_constants
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: n_cond
    integer                  :: IJK, nc
    real(dp)                 :: fld($D), vel($D), diffc, mu_par, mu_cross, tmp
    real(dp)                 :: dt_vec(n_cond)
    type(LT2_loc_t) :: loc

    nc = box%n_cell
    dt_vec = ST_dt_max

    do KJI_DO(1,nc)
       fld = box%cc(IJK, i_Ex:i_Ex+$D-1)

       call get_velocity(fld, vel, loc)
       diffc = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
       mu_par = LT2_get_col_at_loc(ST_td_tbl, i_mobility_B, loc)
       mu_cross = LT2_get_col_at_loc(ST_td_tbl, i_mobility_ExB, loc)

       ! The 0.5 is here because of the explicit trapezoidal rule
       tmp = 0.5_dp / max(sum(abs(vel) / box%dr), epsilon(1.0_dp))
       dt_vec(ST_ix_cfl) = min(dt_vec(ST_ix_cfl), tmp)

       ! Dielectric relaxation time (using highest mobility)
       dt_vec(ST_ix_drt) = min(dt_vec(ST_ix_drt), &
            UC_eps0 / (UC_elem_charge * max(mu_par, mu_cross) * &
            max(box%cc(IJK, i_electron), epsilon(1.0_dp))))

       ! Diffusion condition
       dt_vec(ST_ix_diff) = min(dt_vec(ST_ix_diff), &
            0.25_dp * box%dr**2 / max(diffc, epsilon(1.0_dp)))
    end do; CLOSE_DO

  end function get_max_dt

  real(dp) function angle_degrees(vec_par, vec_norm)
    real(dp), intent(in) :: vec_par($D), vec_norm
    real(dp), parameter :: fac = 180 / acos(-1.0_dp)

    angle_degrees = acos(norm2(vec_par)/vec_norm)
    angle_degrees = angle_degrees * fac

    ! Lies between 0-180 now, convert to 0-90
    if (angle_degrees > 90) then
       angle_degrees = 180 - angle_degrees
    end if
  end function angle_degrees

  subroutine split_field(E, E_par, E_perp, E_cross, E_norm)
    real(dp), intent(in) :: E($D)
    real(dp), intent(out) :: E_par($D), E_perp($D)
    real(dp), intent(out) :: E_cross($D), E_norm
    real(dp) :: E_hat($D)
#if $D == 2
    real(dp) :: tmp(3)
#endif

    E_norm = norm2(E)
    E_hat = E / E_norm

    E_par = sum(E * B0_hat(1:$D)) * B0_hat(1:$D)
    E_perp = E - E_par

#if $D == 2
    tmp(1:2) = E
    tmp(3) = 0
    tmp = cross_product(tmp, B0_hat)
    E_cross = tmp(1:2)
#elif $D == 3
    E_cross = cross_product(E, B0_hat)
#endif

  end subroutine split_field

  !> Return the cross product of vectors a and b
  pure function cross_product(a, b) result(vec)
    real(dp), intent(in) :: a(3)
    real(dp), intent(in) :: b(3)
    real(dp)             :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function cross_product

  function get_min(a, b, n) result(min_vec)
    integer, intent(in)  :: n
    real(dp), intent(in) :: a(n), b(n)
    real(dp)             :: min_vec(n)

    min_vec = min(a, b)
  end function get_min

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    use m_flux_schemes
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    type(LT2_loc_t) :: loc
    real(dp)                     :: inv_dr, fld($D), vel($D)
    real(dp), allocatable        :: v(DTIMES(:), :)
    real(dp), allocatable        :: dc(DTIMES(:), :)
    real(dp), allocatable        :: cc(DTIMES(:))
    integer                      :: nc, n, m
#if $D == 3
    integer                      :: l
#endif

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    allocate(v(DTIMES(1:nc+1), $D))
    allocate(dc(DTIMES(1:nc+1), $D))
    allocate(cc(DTIMES(-1:nc+2)))

    ! Fill ghost cells
    call a$D_gc_box(boxes, id, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero)
    call a$D_gc2_box(boxes, id, i_electron, a$D_gc2_prolong_linear, &
         a$D_bc2_neumann_zero, cc, nc)

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    do n = 1, nc+1
       do m = 1, nc
#if $D == 2
          fld = 0.5_dp * (boxes(id)%cc(n-1, m, i_Ex:i_Ey) + &
               boxes(id)%cc(n, m, i_Ex:i_Ey))
          fld(1) = inv_dr * (boxes(id)%cc(n-1, m, i_phi) - &
               boxes(id)%cc(n, m, i_phi))
          call get_velocity(fld, vel, loc)
          v(n, m, 1)  = vel(1)
          dc(n, m, 1) = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)

          fld       = 0.5_dp * (boxes(id)%cc(m, n-1, i_Ex:i_Ey) + &
               boxes(id)%cc(m, n, i_Ex:i_Ey))
          fld(2) = inv_dr * (boxes(id)%cc(m, n-1, i_phi) - &
               boxes(id)%cc(m, n, i_phi))
          call get_velocity(fld, vel, loc)
          v(m, n, 2)  = vel(2)
          dc(m, n, 2) = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
#elif $D == 3
          do l = 1, nc
             fld = 0.5_dp * (&
                  boxes(id)%cc(n-1, m, l, i_Ex:i_Ez) + &
                  boxes(id)%cc(n, m, l, i_Ex:i_Ez))
             fld(1) = inv_dr * (boxes(id)%cc(n-1, m, l, i_phi) - &
                  boxes(id)%cc(n, m, l, i_phi))
             call get_velocity(fld, vel, loc)
             v(n, m, l, 1)  = vel(1)
             dc(n, m, l, 1) = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, n-1, l, i_Ex:i_Ez) + &
                  boxes(id)%cc(m, n, l, i_Ex:i_Ez))
             fld(2) = inv_dr * (boxes(id)%cc(m, n-1, l, i_phi) - &
               boxes(id)%cc(m, n, l, i_phi))
             call get_velocity(fld, vel, loc)
             v(m, n, l, 2)  = vel(2)
             dc(m, n, l, 2) = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, l, n-1, i_Ex:i_Ez) + &
                  boxes(id)%cc(m, l, n, i_Ex:i_Ez))
             fld(3) = inv_dr * (boxes(id)%cc(m, l, n-1, i_phi) - &
               boxes(id)%cc(m, l, n, i_phi))
             call get_velocity(fld, vel, loc)
             v(m, l, n, 3)  = vel(3)
             dc(m, l, n, 3) = LT2_get_col_at_loc(ST_td_tbl, i_diffusion, loc)
          end do
#endif
       end do
    end do

    call flux_koren_$Dd(cc, v, nc, 2)
    call flux_diff_$Dd(cc, dc, inv_dr, nc, 2)

    boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc
  end subroutine fluxes_koren

  subroutine get_velocity(fld, vel, loc)
    real(dp), intent(in) :: fld($D)
    real(dp), intent(out) :: vel($D)
    type(LT2_loc_t), intent(out) :: loc

    real(dp)                 :: angle, E_par($D), E_perp($D), E_cross($D)
    real(dp)                 :: E_norm, mu_par, mu_perp, mu_cross

    call split_field(fld, E_par, E_perp, E_cross, E_norm)
    angle = angle_degrees(E_par, E_norm)

    loc = LT2_get_loc(ST_td_tbl, [angle, E_norm])
    mu_par = LT2_get_col_at_loc(ST_td_tbl, i_mobility_B, loc)
    mu_perp = LT2_get_col_at_loc(ST_td_tbl, i_mobility_xB, loc)
    mu_cross = LT2_get_col_at_loc(ST_td_tbl, i_mobility_ExB, loc)
    vel = -mu_par * E_par - mu_perp * E_perp + mu_cross * E_cross
  end subroutine get_velocity

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box$D_t), intent(inout) :: box
    box%cc(DTIMES(:), i_electron) = 0.5_dp *  &
         (box%cc(DTIMES(:), i_electron) + box%cc(DTIMES(:), i_electron_old))
    box%cc(DTIMES(:), i_pos_ion) = 0.5_dp * &
         (box%cc(DTIMES(:), i_pos_ion)  + box%cc(DTIMES(:), i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: dt(:)
    real(dp)                     :: inv_dr, src, fld($D),vel($D)
    real(dp)                     :: alpha, eta, sflux
#if $D == 2
    real(dp)                     :: rfac(2)
    integer                      :: ioff
#endif
    integer                      :: IJK, nc
    type(LT2_loc_t)               :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
#if $D == 2
    ioff   = (box%ix(1)-1) * nc
#endif

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_Ex:i_Ex+$D-1)
       call get_velocity(fld, vel, loc)

       alpha = LT2_get_col_at_loc(ST_td_tbl, i_alpha, loc)
       eta   = LT2_get_col_at_loc(ST_td_tbl, i_eta, loc)

       ! Contribution of flux
#if $D == 2
       if (ST_cylindrical) then
          ! Weighting of flux contribution for cylindrical coordinates
          rfac(:) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
       else
          rfac(:) = 1.0_dp
       end if

       sflux = (box%fc(i, j, 2, flux_elec) - box%fc(i, j+1, 2, flux_elec) + &
            rfac(1) * box%fc(i, j, 1, flux_elec) - &
            rfac(2) * box%fc(i+1, j, 1, flux_elec)) * inv_dr * dt(1)
#elif $D == 3
       sflux = (sum(box%fc(i, j, k, 1:3, flux_elec)) - &
            box%fc(i+1, j, k, 1, flux_elec) - &
            box%fc(i, j+1, k, 2, flux_elec) - &
            box%fc(i, j, k+1, 3, flux_elec)) * inv_dr * dt(1)
#endif

       ! Source term
       src = norm2(vel) * box%cc(IJK, i_electron) * (alpha - eta) * dt(1)

       ! Add flux and source term
       box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + sflux + src

       ! Add source term
       box%cc(IJK, i_pos_ion)  = box%cc(IJK, i_pos_ion) + src

    end do; CLOSE_DO
  end subroutine update_solution

  subroutine box_set_src_rate(box)
    type(box$D_t), intent(inout) :: box
    real(dp)                     :: fld($D), vel($D)
    integer                      :: IJK, nc
    type(LT2_loc_t)               :: loc

    nc     = box%n_cell

    do KJI_DO(1,nc)
       fld   = box%cc(IJK, i_Ex:i_Ex+$D-1)
       call get_velocity(fld, vel, loc)

       box%cc(IJK, i_src_rate) = &
            LT2_get_col_at_loc(ST_td_tbl, i_alpha, loc) * norm2(vel)
    end do; CLOSE_DO
  end subroutine box_set_src_rate

  subroutine add_velocity(box, v, n_var)
    type(box$D_t), intent(in) :: box
    integer, intent(in)       :: n_var
#if $D == 2
    real(dp)                  :: v(0:box%n_cell+1, 0:box%n_cell+1, n_var)
#elif $D == 3
    real(dp)                  :: v(0:box%n_cell+1, 0:box%n_cell+1, &
         0:box%n_cell+1, n_var)
#endif

    integer         :: IJK
    real(dp)        :: fld($D)
    type(LT2_loc_t) :: loc

    do KJI_DO(0,box%n_cell+1)
       fld   = box%cc(IJK, i_Ex:i_Ex+$D-1)
       call get_velocity(fld, v(IJK, 1:$D), loc)
    end do; CLOSE_DO
  end subroutine add_velocity

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a$D_prolong
    type(a$D_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       !$omp parallel do private(id, p_id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_electron)
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_pos_ion)
          call a$D_prolong_linear(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do
       !$omp end parallel do

       !$omp parallel do private(id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a$D_gc_box(tree%boxes, id, i_electron, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero)
          call a$D_gc_box(tree%boxes, id, i_pos_ion, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero)
          call a$D_gc_box(tree%boxes, id, i_phi, &
               mg%sides_rb, mg%sides_bc)
       end do
       !$omp end parallel do
    end do
  end subroutine prolong_to_new_boxes

  subroutine write_log_file(tree, filename, out_cnt, dir)
    type(a$D_t), intent(in)      :: tree
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: out_cnt
    character(len=*), intent(in) :: dir
    character(len=ST_slen)       :: fname
    character(len=20)            :: fmt
    integer, parameter           :: my_unit = 123
    real(dp)                     :: velocity
    real(dp), save               :: prev_pos($D) = 0
    real(dp)                     :: sum_elec, sum_pos_ion
    real(dp)                     :: max_elec, max_field
    type(a$D_loc_t)              :: loc_elec, loc_field

    call a$D_prepend_directory(filename, dir, fname)

    call a$D_tree_sum_cc(tree, i_electron, sum_elec)
    call a$D_tree_sum_cc(tree, i_pos_ion, sum_pos_ion)
    call a$D_tree_max_cc(tree, i_electron, max_elec, loc_elec)
    call a$D_tree_max_cc(tree, i_electric_fld, max_field, loc_field)

    if (out_cnt == 1) then
       open(my_unit, file=trim(fname), action="write")
#if $D == 2
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) "&
            &"max(E) x y max(n_e) x y"
#elif $D == 3
       write(my_unit, *) "# it time dt v sum(n_e) sum(n_i) "&
            &"max(E) x y z max(n_e) x y z"
#endif
       close(my_unit)

       ! Start with velocity zero
       prev_pos = a$D_r_loc(tree, loc_field)
    end if

#if $D == 2
    fmt = "(I6,11E16.8)"
#elif $D == 3
    fmt = "(I6,13E16.8)"
#endif

    velocity = norm2(a$D_r_loc(tree, loc_field) - prev_pos) / ST_dt_output
    prev_pos = a$D_r_loc(tree, loc_field)

    open(my_unit, file=trim(fname), action="write", &
         position="append")
    write(my_unit, fmt) out_cnt, ST_time, ST_dt, velocity, sum_elec, &
         sum_pos_ion, max_field, a$D_r_loc(tree, loc_field), max_elec, &
         a$D_r_loc(tree, loc_elec)
    close(my_unit)

  end subroutine write_log_file

end program streamer_$Dd
