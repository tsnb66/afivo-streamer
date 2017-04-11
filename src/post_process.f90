program post_process
  use m_a2_all

  type(a2_t) :: tree
  type(a2_t) :: tree_in
  integer, parameter :: box_size = 128
  double precision :: dr = 8.0e-3 / box_size

  call a2_read_tree(tree_in, "output/streamer_cyl_000001.dat")

  call a2_init(tree, box_size, 2, 0, dr, r_min=[0.0d0, 4.0d-3])
  call a2_set_base(tree, 1, reshape([1,1], [1,2]))

  call a2_tree_copy_variable(tree_in, [1,6], tree, [1,2])
  call a2_write_silo(tree, "test", 1, 0.0d0)

end program post_process
