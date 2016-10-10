!> Module  that contains subroutines to handle with configuration files
!!
!! * Creating configuration files
!! * Reading configuration files
!! * Loading configuration files
!! * Setting the initial conditions from the configuration
!!
!! Further it contains a subroutine for
!! * loading the transport data
!! * Setting  the voltage on fixed time
!!
!! Moreover, it contains several pre-defined variables like 
!!
!! * Indices of cell-centered variables
!! * Names of the cell-centered variables
!! * Indices of face-centered variables
!! * Indices of transport data
!! * How many multigrid FMG cycles we perform per time step
!!
!! The subroutines and variables are used by the examples:
!! * streamer_2d
!! * streamer_3d
!! * streamer_cyl
!!
module m_streamer

  use m_config
  use m_random
  use m_photons
  use m_lookup_table

  implicit none
  public

  integer, parameter, private :: dp = kind(0.0d0)

  ! Default length of strings
  integer, parameter :: ST_slen = 200

  ! ** Indices of cell-centered variables **
  integer, parameter :: n_var_cell     = 8 ! Number of variables
  integer, parameter :: i_electron     = 1 ! Electron density
  integer, parameter :: i_pos_ion      = 2 ! Positive ion density
  integer, parameter :: i_electron_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pos_ion_old  = 4 ! For time-stepping scheme
  integer, parameter :: i_phi          = 5 ! Electrical potential
  integer, parameter :: i_electric_fld = 6 ! Electric field norm
  integer, parameter :: i_rhs          = 7 ! Source term Poisson
  integer, parameter :: i_photo        = 8 ! Phototionization rate

  ! Names of the cell-centered variables
  character(len=12) :: ST_cc_names(n_var_cell) = &
       [character(len=12) :: "electron", "pos_ion", "electron_old", &
       "pos_ion_old", "phi", "electric_fld", "rhs", "pho"]

  ! ** Indices of face-centered variables **
  integer, parameter :: n_var_face   = 2 ! Number of variables
  integer, parameter :: flux_elec    = 1 ! Electron flux
  integer, parameter :: electric_fld = 2 ! Electric field vector

  ! ** Indices of transport data **
  integer, parameter :: n_var_td    = 4 ! Number of transport coefficients
  integer, parameter :: i_mobility  = 1 ! Electron mobility
  integer, parameter :: i_diffusion = 2 ! Electron diffusion constant
  integer, parameter :: i_alpha     = 3 ! Ionization coeff (1/m)
  integer, parameter :: i_eta       = 4 ! Attachment coeff (1/m)

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

  ! Type to store initial conditions in
  type initcnd_t
     real(dp)              :: background_density
     integer               :: n_cond
     real(dp), allocatable :: seed_r0(:, :)
     real(dp), allocatable :: seed_r1(:, :)
     real(dp), allocatable :: seed_density(:)
     integer, allocatable  :: seed_charge_type(:)
     real(dp), allocatable :: seed_width(:)
     integer, allocatable  :: seed_falloff(:)
  end type initcnd_t

  ! This will contain the initial conditions
  type(initcnd_t), protected :: ST_init_cond

  ! Table with transport data vs electric field
  type(lookup_table_t), protected :: ST_td_tbl

  ! The configuration for the simulation
  type(CFG_t) :: ST_config

  ! Random number generator
  type(RNG_t) :: ST_rng

  ! Whether we use phototionization
  logical, protected :: ST_photoi_enabled

  ! Oxygen fraction
  real(dp), protected :: ST_photoi_frac_O2

  ! Photoionization efficiency
  real(dp), protected :: ST_photoi_eta

  ! Number of photons to use
  integer, protected :: ST_photoi_num_photons

  ! Table for photoionization
  type(photoi_tbl_t), protected :: ST_photoi_tbl

  ! Start modifying the vertical background field after this time
  real(dp), protected :: ST_electric_fld_y_mod_t0

  ! Amplitude of sinusoidal modification
  real(dp), protected :: ST_electric_fld_y_sin_amplitude

  ! Frequency (Hz) of sinusoidal modification
  real(dp), protected :: ST_electric_fld_y_sin_freq

  ! Linear derivative of background field
  real(dp), protected :: ST_electric_fld_y_lin_deriv

  ! Decay time of background field
  real(dp), protected :: ST_electric_fld_y_decay

  ! Start modifying the horizontal background field after this time
  real(dp), protected :: ST_electric_fld_x_mod_t0

  ! Amplitude of sinusoidal modification
  real(dp), protected :: ST_electric_fld_x_sin_amplitude

  ! Frequency (Hz) of sinusoidal modification
  real(dp), protected :: ST_electric_fld_x_sin_freq

  ! Linear derivative of background field
  real(dp), protected :: ST_electric_fld_x_lin_deriv

  ! Decay time of background field
  real(dp), protected :: ST_electric_fld_x_decay

  ! Name of the simulations
  character(len=ST_slen), protected :: ST_simulation_name

  ! Output directory
  character(len=ST_slen), protected :: ST_output_dir

  ! The number of steps after which the mesh is updated
  integer, protected :: ST_refine_per_steps

  ! The grid spacing will always be larger than this value
  real(dp), protected :: ST_refine_min_dx

  ! The grid spacing will always be smaller than this value
  real(dp), protected :: ST_refine_max_dx

  ! Refine if alpha*dx is larger than this value
  real(dp), protected :: ST_refine_adx

  ! Refine if the curvature in phi is larger than this value
  real(dp), protected :: ST_refine_cphi

  ! Only derefine if grid spacing if smaller than this value
  real(dp), protected :: ST_derefine_dx

  ! Refine around initial conditions up to this time
  real(dp), protected :: ST_refine_init_time

  ! Refine until dx is smaller than this factor times the seed width
  real(dp), protected :: ST_refine_init_fac

  ! Number of output files written
  integer :: ST_out_cnt

  ! Current time step
  real(dp) :: ST_dt

  ! Maximum allowed time step
  real(dp), protected :: ST_dt_max

  ! Time between writing output
  real(dp), protected :: ST_dt_out

  ! Current time
  real(dp)  :: ST_time

  ! End time of the simulation
  real(dp), protected :: ST_end_time

  ! The size of the boxes that we use to construct our mesh
  integer, protected :: ST_box_size

  ! The length of the (square) domain
  real(dp), protected :: ST_domain_len

  ! The applied electric field (vertical direction)
  real(dp), protected :: ST_applied_electric_fld_y

  ! The applied electric field (horizontal direction)
  real(dp), protected :: ST_applied_electric_fld_x

  ! The applied voltage (vertical direction)
  real(dp), protected :: ST_applied_voltage

  ! Pressure of the gas in bar
  real(dp), protected :: ST_gas_pressure

  ! Dielectric constant
  real(dp), protected :: ST_epsilon_diel

contains

  !> Create the configuration file with default values
  subroutine ST_create_config()

    call CFG_add(ST_config, "end_time", 10.0d-9, &
         "The desired endtime (s) of the simulation")
    call CFG_add(ST_config, "simulation_name", "sim", &
         "The name of the simulation")
    call CFG_add(ST_config, "output_dir", "", &
         "Directory where the output should be written")
    call CFG_add(ST_config, "box_size", 8, &
         "The number of grid cells per coordinate in a box")
    call CFG_add(ST_config, "domain_len", 32e-3_dp, &
         "The length of the domain (m)")
    call CFG_add(ST_config, "gas_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(ST_config, "gas_pressure", 1.0_dp, &
         "The gas pressure (bar), used for photoionization")
    call CFG_add(ST_config, "applied_electric_fld_y", 1.0d7, &
         "The applied electric field (V/m) (vertical)")
    call CFG_add(ST_config, "applied_electric_fld_x", 0.0d0, &
         "The applied electric field (V/m) (horizontal)")
    call CFG_add(ST_config, "epsilon_diel", 1.5_dp, &
         "The relative dielectric constant (if a dielectric is specified)")

    call CFG_add(ST_config, "electric_fld_y_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time (s)")
    call CFG_add(ST_config, "electric_fld_y_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add(ST_config, "electric_fld_y_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add(ST_config, "electric_fld_y_lin_deriv", 0.0_dp, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add(ST_config, "electric_fld_y_decay", huge(1.0_dp), &
         "Decay time of field (s)")

    call CFG_add(ST_config, "electric_fld_x_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time (s)")
    call CFG_add(ST_config, "electric_fld_x_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification (V/m)")
    call CFG_add(ST_config, "electric_fld_x_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification (Hz)")
    call CFG_add(ST_config, "electric_fld_x_lin_deriv", 0.0_dp, &
         "Linear derivative of field [V/(ms)]")
    call CFG_add(ST_config, "electric_fld_x_decay", huge(1.0_dp), &
         "Decay time of field (s)")

    call CFG_add(ST_config, "background_density", 1.0d12, &
         "The background ion and electron density (1/m3)")
    call CFG_add(ST_config, "seed_density", [5.0d19], &
         "Initial density of the seed (1/m3)", .true.)
    call CFG_add(ST_config, "seed_charge_type", [0], &
         "Type of seed: neutral (0), ions (1) or electrons (-1)", .true.)
    call CFG_add(ST_config, "seed_rel_r0", [0.5d0, 0.4d0], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(ST_config, "seed_rel_r1", [0.5d0, 0.6d0], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(ST_config, "seed_width", [0.5d-3], &
         "Seed width (m)", .true.)
    call CFG_add(ST_config, "seed_falloff", [1], &
         "Fallof type for seed, see m_geometry.f90", .true.)

    call CFG_add(ST_config, "dt_output", 1.0d-10, &
         "The timestep for writing output (s)")
    call CFG_add(ST_config, "dt_max", 1.0d-11, &
         "The maximum timestep (s)")

    call CFG_add(ST_config, "refine_per_steps", 2, &
         "The number of steps after which the mesh is updated")
    call CFG_add(ST_config, "refine_min_dx", 1.0e-6_dp, &
         "The grid spacing will always be larger than this value")
    call CFG_add(ST_config, "refine_max_dx", 1.0e-3_dp, &
         "The grid spacing will always be smaller than this value")

    call CFG_add(ST_config, "refine_adx", 0.8_dp, &
         "Refine if alpha*dx is larger than this value")
    call CFG_add(ST_config, "refine_cphi", 1e99_dp, &
         "Refine if the curvature in phi is larger than this value")
    call CFG_add(ST_config, "derefine_dx", 1e-4_dp, &
         "Only derefine if grid spacing if smaller than this value")
    call CFG_add(ST_config, "refine_init_time", 10.0e-9_dp, &
         "Refine around initial conditions up to this time")
    call CFG_add(ST_config, "refine_init_fac", 0.25_dp, &
         "Refine until dx is smaller than this factor times the seed width")

    call CFG_add(ST_config, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(ST_config, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen (0-1)")
    call CFG_add(ST_config, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor, typically around 0.05")
    call CFG_add(ST_config, "photoi_num_photons", 50*1000, &
         "Number of discrete photons to use for photoionization")
    call CFG_add(ST_config, "photoi_rng_seed", [8123, 91234, 12399, 293434], &
         "Seed for the photoionization random number generator")

    call CFG_add(ST_config, "input_file", "transport_data_file.txt", &
         "Input file with transport data")
    call CFG_add(ST_config, "lookup_table_size", 1000, &
         "The transport data table size in the fluid model")
    call CFG_add(ST_config, "lookup_table_max_electric_fld", 3.0d7, &
         "The maximum electric field in the fluid model coefficients")
    call CFG_add(ST_config, "td_mobility_name", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(ST_config, "td_diffusion_name", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(ST_config, "td_alpha_name", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(ST_config, "td_eta_name", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")

    call CFG_add(ST_config, "td_alpha_fac", 1.0_dp, &
         "Modify alpha by this factor")
    call CFG_add(ST_config, "td_eta_fac", 1.0_dp, &
         "Modify eta by this factor")
    call CFG_add(ST_config, "td_mobility_fac", 1.0_dp, &
         "Modify mobility by this factor")
    call CFG_add(ST_config, "td_diffusion_fac", 1.0_dp, &
         "Modify diffusion by this factor")
  end subroutine ST_create_config

  ! Set the initial conditions from the configuration
  subroutine ST_get_init_cond(n_dim)
    integer, intent(in)            :: n_dim
    integer                        :: n_cond, varsize
    real(dp)                       :: dlen
    real(dp), allocatable          :: tmp_vec(:)
    type(initcnd_t) :: ic

    call CFG_get(ST_config, "background_density", ic%background_density)
    call CFG_get(ST_config, "domain_len", dlen)

    call CFG_get_size(ST_config, "seed_density", n_cond)

    call CFG_get_size(ST_config, "seed_rel_r0", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r0 variable has incompatible size"

    call CFG_get_size(ST_config, "seed_rel_r1", varsize)
    if (varsize /= n_dim * n_cond) &
         stop "seed_rel_r1 variable has incompatible size"

    call CFG_get_size(ST_config, "seed_charge_type", varsize)
    if (varsize /= n_cond) &
         stop "seed_charge_type variable has incompatible size"

    call CFG_get_size(ST_config, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_width variable has incompatible size"

    ic%n_cond = n_cond
    allocate(ic%seed_density(n_cond))
    allocate(ic%seed_charge_type(n_cond))
    allocate(ic%seed_r0(n_dim, n_cond))
    allocate(ic%seed_r1(n_dim, n_cond))
    allocate(ic%seed_width(n_cond))
    allocate(ic%seed_falloff(n_cond))

    allocate(tmp_vec(n_dim * n_cond))
    call CFG_get(ST_config, "seed_rel_r0", tmp_vec)
    ic%seed_r0 = dlen * reshape(tmp_vec, [n_dim, n_cond])
    call CFG_get(ST_config, "seed_rel_r1", tmp_vec)
    ic%seed_r1 = dlen * reshape(tmp_vec, [n_dim, n_cond])

    call CFG_get(ST_config, "seed_density", ic%seed_density)
    call CFG_get(ST_config, "seed_charge_type", ic%seed_charge_type)
    call CFG_get(ST_config, "seed_width", ic%seed_width)
    call CFG_get(ST_config, "seed_falloff", ic%seed_falloff)
    ST_init_cond = ic
  end subroutine ST_get_init_cond

  !> Initialize the transport coefficients
  subroutine ST_load_transport_data()
    use iso_fortran_env, only: int64
    use m_transport_data
    use m_config

    character(len=ST_slen) :: input_file, gas_name
    integer                :: table_size, rng_int4_seed(4)
    integer(int64)         :: rng_int8_seed(2)
    real(dp)               :: max_electric_fld, alpha_fac, eta_fac
    real(dp)               :: mobility_fac, diffusion_fac
    real(dp), allocatable  :: x_data(:), y_data(:)
    character(len=ST_slen) :: data_name

    call CFG_get(ST_config, "input_file", input_file)
    call CFG_get(ST_config, "gas_name", gas_name)

    call CFG_get(ST_config, "lookup_table_size", table_size)
    call CFG_get(ST_config, "lookup_table_max_electric_fld", max_electric_fld)

    call CFG_get(ST_config, "td_alpha_fac", alpha_fac)
    call CFG_get(ST_config, "td_eta_fac", eta_fac)
    call CFG_get(ST_config, "td_mobility_fac", mobility_fac)
    call CFG_get(ST_config, "td_diffusion_fac", diffusion_fac)

    ! Create a lookup table for the model coefficients
    ST_td_tbl = LT_create(0.0_dp, max_electric_fld, table_size, n_var_td)

    ! Fill table with data
    call CFG_get(ST_config, "td_mobility_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(ST_td_tbl, i_mobility, x_data, y_data)

    call CFG_get(ST_config, "td_diffusion_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(ST_td_tbl, i_diffusion, x_data, y_data)

    call CFG_get(ST_config, "td_alpha_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(ST_td_tbl, i_alpha, x_data, y_data)

    call CFG_get(ST_config, "td_eta_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(ST_td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    if (ST_photoi_enabled) then
       call CFG_get(ST_config, "photoi_rng_seed", rng_int4_seed)
       rng_int8_seed = transfer(rng_int4_seed, rng_int8_seed)
       call ST_rng%set_seed(rng_int8_seed)

       call photoi_get_table_air(ST_photoi_tbl, ST_photoi_frac_O2 * &
            ST_gas_pressure, 2 * ST_domain_len)
    end if

  end subroutine ST_load_transport_data

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

  !> Read the values of the configuration file.
  ! Only values which differ from the default values must be listed
  ! in the configuration file
  subroutine ST_read_config_files()
    character(len=ST_slen)     :: tmp_name, prev_name, config_name
    integer                    :: n

    ST_simulation_name = ""
    prev_name = ""
    do n = 1, command_argument_count()
       call get_command_argument(n, config_name)
       call CFG_read_file(ST_config, trim(config_name))

       call CFG_get(ST_config, "simulation_name", tmp_name)
       if (ST_simulation_name == "") then
          ST_simulation_name = tmp_name
       else if (tmp_name /= "" .and. tmp_name /= prev_name) then
          ST_simulation_name = trim(ST_simulation_name) // "_" // trim(tmp_name)
       end if
       prev_name = tmp_name
    end do

  end subroutine ST_read_config_files

  !> Update the configuration file with the values read in 
  ! previous call of ST_read_config_files()
  subroutine ST_load_config()
    character(len=ST_slen)    :: tmp_name

    call CFG_get(ST_config, "end_time", ST_end_time)
    call CFG_get(ST_config, "box_size", ST_box_size)
    call CFG_get(ST_config, "output_dir", ST_output_dir)
    call CFG_get(ST_config, "domain_len", ST_domain_len)
    call CFG_get(ST_config, "applied_electric_fld_y", ST_applied_electric_fld_y)
    call CFG_get(ST_config, "applied_electric_fld_x", ST_applied_electric_fld_x)
    call CFG_get(ST_config, "dt_output", ST_dt_out)

    call CFG_get(ST_config, "refine_per_steps", ST_refine_per_steps)
    call CFG_get(ST_config, "refine_min_dx", ST_refine_min_dx)
    call CFG_get(ST_config, "refine_max_dx", ST_refine_max_dx)
    call CFG_get(ST_config, "refine_adx", ST_refine_adx)
    call CFG_get(ST_config, "refine_cphi", ST_refine_cphi)
    call CFG_get(ST_config, "derefine_dx", ST_derefine_dx)

    call CFG_get(ST_config, "refine_init_time", ST_refine_init_time)
    call CFG_get(ST_config, "refine_init_fac", ST_refine_init_fac)

    call CFG_get(ST_config, "dt_max", ST_dt_max)
    call CFG_get(ST_config, "epsilon_diel", ST_epsilon_diel)

    call CFG_get(ST_config, "gas_pressure", ST_gas_pressure)
    call CFG_get(ST_config, "photoi_enabled", ST_photoi_enabled)
    call CFG_get(ST_config, "photoi_frac_O2", ST_photoi_frac_O2)
    call CFG_get(ST_config, "photoi_eta", ST_photoi_eta)
    call CFG_get(ST_config, "photoi_num_photons", ST_photoi_num_photons)

    call CFG_get(ST_config, "electric_fld_y_mod_t0", ST_electric_fld_y_mod_t0)
    call CFG_get(ST_config, "electric_fld_y_sin_amplitude", ST_electric_fld_y_sin_amplitude)
    call CFG_get(ST_config, "electric_fld_y_sin_freq", ST_electric_fld_y_sin_freq)
    call CFG_get(ST_config, "electric_fld_y_lin_deriv", ST_electric_fld_y_lin_deriv)
    call CFG_get(ST_config, "electric_fld_y_decay", ST_electric_fld_y_decay)

    call CFG_get(ST_config, "electric_fld_x_mod_t0", ST_electric_fld_x_mod_t0)
    call CFG_get(ST_config, "electric_fld_x_sin_amplitude", ST_electric_fld_x_sin_amplitude)
    call CFG_get(ST_config, "electric_fld_x_sin_freq", ST_electric_fld_x_sin_freq)
    call CFG_get(ST_config, "electric_fld_x_lin_deriv", ST_electric_fld_x_lin_deriv)
    call CFG_get(ST_config, "electric_fld_x_decay", ST_electric_fld_x_decay)

    ST_applied_voltage = -ST_domain_len * ST_applied_electric_fld_y

    tmp_name = trim(ST_output_dir) // "/" // trim(ST_simulation_name) // "_output.cfg"
    print *, "Settings written to ", trim(tmp_name)
    call CFG_write(ST_config, trim(tmp_name))

  end subroutine ST_load_config

  !> Compute the electric field at a given time
  function ST_get_electric_field(time) result(electric_fld)
    use m_units_constants

    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - ST_electric_fld_y_mod_t0
    if (t_rel > 0) then
       electric_fld = ST_applied_electric_fld_y * exp(-t_rel/ST_electric_fld_y_decay) + &
                      t_rel * ST_electric_fld_y_lin_deriv + &
                      ST_electric_fld_y_sin_amplitude * &
                      sin(t_rel * ST_electric_fld_y_sin_freq * 2 * UC_pi)
    else
       electric_fld = ST_applied_electric_fld_y
    end if
  end function ST_get_electric_field

  !> Compute the voltage at a given time
  subroutine ST_set_voltage(time)
    use m_units_constants

    real(dp), intent(in) :: time
    real(dp)             :: electric_fld, t_rel

    t_rel = time - ST_electric_fld_y_mod_t0
    if (t_rel > 0) then
       electric_fld = ST_applied_electric_fld_y * exp(-t_rel/ST_electric_fld_y_decay) + &
                      t_rel * ST_electric_fld_y_lin_deriv + &
                      ST_electric_fld_y_sin_amplitude * &
                      sin(t_rel * ST_electric_fld_y_sin_freq * 2 * UC_pi)
    else
       electric_fld = ST_applied_electric_fld_y
    end if
    ST_applied_voltage = -ST_domain_len * electric_fld

  end subroutine ST_set_voltage

end module m_streamer
