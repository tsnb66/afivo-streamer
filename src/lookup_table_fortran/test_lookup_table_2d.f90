program test
  use m_lookup_table_2d
  implicit none

  integer, parameter   :: dp             = kind(0.0d0)
  integer, parameter   :: in_size(2)   = [3, 11]
  integer, parameter   :: tbl_size(2)   = [7, 15]

  integer              :: i, j, cntr, table_size
  type(LT2_t)          :: lkp_tbl
  real(dp)             :: x1(in_size(1))
  real(dp) :: x2(in_size(2))
  real(dp)             :: y_values(in_size(1), in_size(2))

  print *, 'start test_lookup_table_2d'

  do i = 1, in_size(1)
     x1(i) = (i-1) * acos(-1.0_dp) / (in_size(1) - 1)
  end do

  do i = 1, in_size(2)
     x2(i) = (i-1) * acos(-1.0_dp) / (in_size(2) - 1)
  end do


  ! Create some testing data
  do j = 1, in_size(2)
     do i = 1, in_size(1)
        ! y_values(i, j) = cos(x1(i)) * cos(x2(j))
        y_values(i, j) = 2 * x1(i) + 3 * x2(j)
     end do
  end do

  lkp_tbl = LT2_create([x1(1), x2(1)], &
       [x1(in_size(1)), x2(in_size(2))], tbl_size, 1)
  call LT2_set_col(lkp_tbl, 1, x1, x2, y_values)

  do j = 1, in_size(2)
     do i = 1, in_size(1)
        print *, x1(i), x2(j), y_values(i, j), &
             LT2_get_col(lkp_tbl, 1, [x1(i), x2(j)])
     end do
  end do

  print *, "Testing m_lookup_table_2d.f90 implementation..."

end program test
