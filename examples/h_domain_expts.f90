#include "../src/cpp_macros.h"
program computational_domain
  use m_af_all

  implicit none

  integer    :: n_cell = 4 ! Boxes contain n_cell^dim cells
  integer, parameter :: n_boxes = 2
  type(af_t) :: tree       ! Will contain the quad/octree grid
  integer    :: grid_size(NDIM)
  real(dp)   :: domain_size(NDIM)
  logical    :: periodic(NDIM)

  ! Create mesh
  domain_size(1) = 1.0e-2_dp
  domain_size(2) = 5.0e-3_dp
  grid_size(:) = n_cell !One box in the y direction
  grid_size(1) = n_boxes * n_cell !n_boxes in the x direction

  call af_add_cc_variable(tree, "phi")
  call af_init(tree, n_cell, domain_size, grid_size)
  call af_write_vtk(tree, "comd_x" // DIMNAME // "_1", dir="output")
  call af_destroy(tree)

  ! Create mesh 2: two boxes along y-direction
  !grid_size(:) = n_cell
  !grid_size(2) = 2 * n_cell
  !domain_size = 1.0_dp * grid_size

  !call af_add_cc_variable(tree, "phi")
  !call af_init(tree, n_cell, domain_size, grid_size)
  !call af_write_vtk(tree, "compd_y" // DIMNAME // "_2", dir="output")
  !call af_destroy(tree)

  ! Create mesh 3: Two boxes along x-direction that are fully periodic
  !grid_size(:) = n_cell
  !grid_size(1) = 2 * n_cell
  !periodic(:)  = .true.
  !domain_size  = 1.0_dp * grid_size

  !call af_add_cc_variable(tree, "phi")
  !call af_init(tree, n_cell, domain_size, grid_size, periodic=periodic)
  !call af_write_vtk(tree, "computational_domain_" // DIMNAME // "_3", dir="output")
  !call af_destroy(tree)

end program computational_domain
