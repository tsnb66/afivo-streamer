#include "../afivo/src/cpp_macros.h"
module m_chemistry
  use m_types
  use m_af_all
  use m_lookup_table
  use m_table_data

  implicit none
  private

  type reaction_t
     integer, allocatable  :: ix_in(:)
     integer, allocatable  :: ix_out(:)
     integer, allocatable  :: multiplicity_out(:)
     integer               :: rate_type
     real(dp)              :: rate_factor
     real(dp), allocatable :: x_data(:)
     real(dp), allocatable :: y_data(:)
     character(len=50)     :: description
  end type reaction_t

  type fast_react_t
     integer, allocatable  :: ix_in(:)
     integer, allocatable  :: ix_out(:)
     integer, allocatable  :: multiplicity_out(:)
  end type fast_react_t

  integer, parameter :: constant_rate        = 1
  integer, parameter :: field_dependent_rate = 2
  integer, parameter :: max_num_species      = 100
  integer, parameter :: max_num_reactions    = 100

  integer, public, protected :: n_species = 0
  integer, public, protected :: n_reactions = 0

  character(len=comp_len), public, protected :: species_list(max_num_species)
  integer, public, protected                 :: species_charge(max_num_species) = 0
  integer, public, protected                 :: species_ix(max_num_species)
  type(reaction_t), public, protected        :: reactions(max_num_reactions)
  type(fast_react_t)                         :: fast_react(max_num_reactions)
  type(LT_t)                                 :: chemtbl

  integer, allocatable, protected :: charged_species_ix(:)
  integer, allocatable, protected :: charged_species_charge(:)

  public :: charged_species_ix
  public :: charged_species_charge

  public :: chemistry_initialize
  public :: get_rates
  public :: get_derivatives

  public :: species_index

contains

  subroutine chemistry_initialize(tree, cfg)
    use m_config
    use m_units_constants
    use m_table_data
    use m_transport_data
    use m_gas
    type(af_t), intent(inout)  :: tree
    type(CFG_t), intent(inout) :: cfg
    integer                    :: n, i
    character(len=string_len)  :: reaction_file
    character(len=comp_len)    :: tmp_name

    reaction_file = "UNDEFINED"
    call CFG_add_get(cfg, "chemistry%reaction_file", reaction_file, &
         "File with a list of reactions")

    if (reaction_file == "UNDEFINED") then
       print *, "m_chemistry: no reactions defined, using standard model"

       species_list(1) = "e"
       species_list(2) = "M+"
       species_list(3) = "M-"
       n_species       = 3
       n_reactions     = 2

       ! Ionization reaction
       if (gas_constant_density) then
          reactions(1)%ix_in = [1]
          reactions(1)%ix_out = [1, 2]
          reactions(1)%multiplicity_out = [2, 1]
          reactions(1)%rate_type = field_dependent_rate
          reactions(1)%rate_factor = 1.0_dp
          reactions(1)%x_data = &
               LT_get_xdata(td_tbl%x_min, td_tbl%dx, td_tbl%n_points)
          reactions(1)%y_data = td_tbl%rows_cols(:, td_alpha) * &
               td_tbl%rows_cols(:, td_mobility) * reactions(1)%x_data * &
               Townsend_to_SI
          reactions(1)%description = "e + M > e + e + M+"

          ! Attachment reaction
          reactions(2)%ix_in = [1]
          reactions(2)%ix_out = [3]
          reactions(2)%multiplicity_out = [1]
          reactions(2)%rate_type = field_dependent_rate
          reactions(2)%rate_factor = 1.0_dp
          reactions(2)%x_data = &
               LT_get_xdata(td_tbl%x_min, td_tbl%dx, td_tbl%n_points)
          reactions(2)%y_data = td_tbl%rows_cols(:, td_eta) * &
               td_tbl%rows_cols(:, td_mobility) * reactions(2)%x_data * &
               Townsend_to_SI
          reactions(2)%description = "e + M > M-"
       else
          error stop "Varying gas density not yet supported"
       end if
    else
       call read_reactions(reaction_file)
    end if

    ! Convert names to simple ascii
    do n = 1, n_species
       tmp_name = species_list(n)
       call to_simple_ascii(trim(tmp_name), species_list(n), &
            species_charge(n))
    end do

    chemtbl = LT_create(table_min_townsend, table_max_townsend, &
         table_size, n_reactions)

    do n = 1, n_reactions
       select case (reactions(n)%rate_type)
       case (constant_rate, field_dependent_rate)
          call LT_set_col(chemtbl, n, reactions(n)%x_data, &
               reactions(n)%rate_factor * reactions(n)%y_data)
       case default
          error stop "Unknown type of reaction rate"
       end select
    end do

    ! Also store in more memory-efficient structure
    do n = 1, n_reactions
       print *, reactions(n)%description
       print *, reactions(n)%ix_in
       print *, reactions(n)%ix_out
       print *, reactions(n)%multiplicity_out
       fast_react(n)%ix_in            = reactions(n)%ix_in
       fast_react(n)%ix_out           = reactions(n)%ix_out
       fast_react(n)%multiplicity_out = reactions(n)%multiplicity_out
    end do

    do n = 1, n_species
       call af_add_cc_variable(tree, trim(species_list(n)), &
            n_copies=2, ix=species_ix(n))
    end do

    ! Store list with only charged species
    n = count(species_charge(1:n_species) /= 0)
    allocate(charged_species_ix(n))
    allocate(charged_species_charge(n))

    i = 0
    do n = 1, n_species
       if (species_charge(n) /= 0) then
          i = i + 1
          charged_species_ix(i) = species_ix(n)
          charged_species_charge(i) = species_charge(n)
       end if
    end do

  end subroutine chemistry_initialize

  !> Compute reaction rates
  subroutine get_rates(fields, rates, n_cells)
    integer, intent(in)   :: n_cells
    real(dp), intent(in)  :: fields(n_cells)
    real(dp), intent(out) :: rates(n_cells, n_reactions)
    integer               :: n

    ! This order looks inefficient, but it is faster to look up multiple
    ! reactions at the same field
    do n = 1, n_cells
       rates(n, :) = LT_get_mcol(chemtbl, fields(n))
    end do
  end subroutine get_rates

  !> Compute derivatives due to chemical reactions
  subroutine get_derivatives(dens, rates, derivs, n_cells)
    integer, intent(in)   :: n_cells
    real(dp), intent(in)  :: dens(n_cells, n_species)
    real(dp), intent(in)  :: rates(n_cells, n_reactions)
    real(dp), intent(out) :: derivs(n_cells, n_species)
    real(dp)              :: prate(n_cells, n_reactions)
    integer               :: n, i, ix

    derivs(:, :) = 0.0_dp

    ! Loop over reactions and add to derivative
    do n = 1, n_reactions
       ! Determine production rate of 'full' reaction
       prate(:, n) = rates(:, n) * &
            product(dens(:, fast_react(n)%ix_in), dim=2)

       ! Input species are removed
       do i = 1, size(fast_react(n)%ix_in)
          ix = fast_react(n)%ix_in(i)
          derivs(:, ix) = derivs(:, ix) - prate(:, n)
       end do

       ! Output species are created
       do i = 1, size(fast_react(n)%ix_out)
          ix = fast_react(n)%ix_out(i)
          derivs(:, ix) = derivs(:, ix) + prate(:, n) * &
               fast_react(n)%multiplicity_out(i)
       end do
    end do
  end subroutine get_derivatives

  subroutine read_reactions(filename)
    character(len=*), intent(in) :: filename
    character(len=string_len)    :: line, data_value
    character(len=50)            :: reaction, how_to_get
    type(reaction_t)             :: new_reaction
    integer                      :: my_unit
    integer, parameter           :: n_fields = 3
    integer                      :: i0(n_fields), i1(n_fields)
    integer                      :: n_found, i

    open(newunit=my_unit, file=filename, action="read")

    n_reactions = 0

    do
       read(my_unit, "(A)", end=999) line
       line = adjustl(line)
       if (line(1:1) == "#") cycle
       call get_fields_string(line, ",", n_fields, n_found, i0, i1)

       if (n_found /= n_fields) then
          print *, trim(line)
          error stop "Invalid chemistry syntax"
       end if

       reaction = line(i0(1):i1(1))
       how_to_get = line(i0(2):i1(2))
       data_value = line(i0(3):i1(3))

       call get_reaction(trim(reaction), new_reaction)
       new_reaction%description = trim(reaction)

       select case (how_to_get)
       case ("table")
          ! Format is filename:reaction
          i = index(data_value, ":")
          call read_reaction_table(data_value(1:i-1), &
               trim(data_value(i+1:)), new_reaction)
       case ("constant")
          call read_reaction_constant(trim(data_value), new_reaction)
       end select

       n_reactions            = n_reactions + 1
       reactions(n_reactions) = new_reaction
    end do
    close(my_unit)

999 continue
  end subroutine read_reactions

  subroutine read_reaction_constant(text, rdata)
    character(len=*), intent(in)    :: text
    type(reaction_t), intent(inout) :: rdata
    real(dp)                        :: tmp

    rdata%rate_type = constant_rate
    read(text, *) tmp

    rdata%x_data = [table_min_townsend, table_max_townsend]
    rdata%y_data = [tmp, tmp]
  end subroutine read_reaction_constant

  subroutine read_reaction_table(filename, dataname, rdata)
    use m_transport_data
    character(len=*), intent(in)    :: filename
    character(len=*), intent(in)    :: dataname
    type(reaction_t), intent(inout) :: rdata

    rdata%rate_type = field_dependent_rate
    call table_from_file(filename, dataname, rdata%x_data, rdata%y_data)
  end subroutine read_reaction_table

  subroutine get_reaction(reaction_text, reaction)
    use m_gas
    character(len=*), intent(in)  :: reaction_text
    type(reaction_t), intent(out) :: reaction
    integer, parameter            :: max_components = 100
    character(len=comp_len)       :: component
    integer                       :: i, ix, n, n_found, multiplicity
    integer                       :: n_in, n_out
    integer                       :: i0(max_components)
    integer                       :: i1(max_components)
    integer                       :: ix_in(max_components)
    integer                       :: ix_out(max_components)
    integer                       :: multiplicity_out(max_components)
    logical                       :: left_side
    real(dp)                      :: rfactor

    call get_fields_string(reaction_text, " ", max_components, n_found, i0, i1)

    left_side = .true.
    n_in      = 0
    n_out     = 0
    rfactor   = 1.0_dp

    do n = 1, n_found
       component = reaction_text(i0(n):i1(n))

       if (component == "+") cycle

       if (component == "->") then
          left_side = .false.
          cycle
       end if

       ! Assume we have a multiplicity less than 10
       if (lge(component(1:1), '1') .and. lle(component(1:1), '9')) then
          read(component(1:1), *) multiplicity
          component = component(2:)
       else
          multiplicity = 1
       end if

       ix = gas_index(component)
       if (ix /= -1 .and. gas_constant_density) then
          ! Multiply reaction by constant density
          if (left_side) rfactor = rfactor * gas_densities(ix)
          cycle
       end if

       ix = species_index(component)
       if (ix == -1) then
          n_species        = n_species + 1
          ix               = n_species
          species_list(ix) = trim(component)
       end if

       if (left_side) then
          do i = 1, multiplicity
             n_in = n_in + 1
             ix_in(n_in) = ix
          end do
       else
          ! Check if species is already present in right-hand side
          do i = 1, n_out
             if (ix == ix_out(i)) exit
          end do

          if (i <= n_out) then
             multiplicity_out(i) = multiplicity_out(i) + multiplicity
          else
             ! If not already present, add the species
             n_out = n_out + 1
             ix_out(n_out) = ix
             multiplicity_out(n_out) = multiplicity
          end if
       end if
    end do

    reaction%ix_in            = ix_in(1:n_in)
    reaction%ix_out           = ix_out(1:n_out)
    reaction%multiplicity_out = multiplicity_out(1:n_out)
    reaction%rate_factor      = rfactor
  end subroutine get_reaction

  !> Find index of a species, return -1 if not found
  elemental integer function species_index(name)
    character(len=*), intent(in) :: name
    do species_index = 1, n_species
       if (species_list(species_index) == name) exit
    end do
    if (species_index == n_species+1) species_index = -1
  end function species_index

  !> Routine to find the indices of entries in a string
  subroutine get_fields_string(line, delims, n_max, n_found, ixs_start, ixs_end)
    !> The line from which we want to read
    character(len=*), intent(in)  :: line
    !> A string with delimiters. For example delims = " ,'"""//char(9)
    character(len=*), intent(in)  :: delims
    !> Maximum number of entries to read in
    integer, intent(in)           :: n_max
    !> Number of entries found
    integer, intent(inout)        :: n_found
    !> On return, ix_start(i) holds the starting point of entry i
    integer, intent(inout)        :: ixs_start(n_max)
    !> On return, ix_end(i) holds the end point of entry i
    integer, intent(inout)        :: ixs_end(n_max)

    integer                       :: ix, ix_prev

    ix_prev = 0
    n_found = 0

    do while (n_found < n_max)

       ! Find the starting point of the next entry (a non-delimiter value)
       ix = verify(line(ix_prev+1:), delims)
       if (ix == 0) exit

       n_found            = n_found + 1
       ixs_start(n_found) = ix_prev + ix ! This is the absolute position in 'line'

       ! Get the end point of the current entry (next delimiter index minus one)
       ix = scan(line(ixs_start(n_found)+1:), delims) - 1

       if (ix == -1) then              ! If there is no last delimiter,
          ixs_end(n_found) = len(line) ! the end of the line is the endpoint
       else
          ixs_end(n_found) = ixs_start(n_found) + ix
       end if

       ix_prev = ixs_end(n_found) ! We continue to search from here
    end do

  end subroutine get_fields_string

  !> An inefficient routine to replace *^+- characters in a string
  subroutine to_simple_ascii(text, simple, charge)
    character(len=*), intent(in)    :: text
    character(len=*), intent(inout) :: simple
    integer, intent(out)            :: charge
    integer                         :: n

    charge = 0
    simple = ""

    do n = 1, len_trim(text)
       select case (text(n:n))
       case ('*')
          simple = trim(simple) // "_star"
       case ('+')
          charge = charge + 1
          simple = trim(simple) // "_plus"
       case ('-')
          charge = charge - 1
          simple = trim(simple) // "_min"
       case ('^')
          simple = trim(simple) // "_hat"
       case default
          simple = trim(simple) // text(n:n)
       end select
    end do

    ! Handle some species separately
    if (simple == "e") charge = -1
  end subroutine to_simple_ascii

end module m_chemistry
