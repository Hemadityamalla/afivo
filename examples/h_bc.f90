#include "../src/cpp_macros.h"

program h_bcs


      use m_af_all
      implicit none

      integer, parameter :: n_cells = 4

      integer :: grid_size(NDIM)
      real(dp) :: domain_size(NDIM)

      integer :: i_phi

      type(af_t) :: tree

      domain_size(:) = 1.0_dp
      domain_size(1) = 2.0_dp

      grid_size(:) = n_cells
      grid_size(1) = 2*n_cells

      print *, "Running boundary conditions_" // DIMNAME // ""
      print *, "Number of threads", af_get_max_threads()

      !i_phi is the index of the variable
      call af_add_cc_variable(tree, "phi", ix=i_phi)

      !Setting the boundary condition for the variable with index i_phi
      call af_set_cc_methods(tree, i_phi, boundary_method)


      !Initialize tree- domain unit cube/square with 2boxes per dim.
      call af_init(tree, n_cells, domain_size, grid_size)

      call af_print_info(tree)

      !initialize the variables to zero
      call af_loop_box(tree, set_phi_zero)

      !Calling the ghost cell routine for i_phi
      call af_gc_tree(tree, [i_phi])

      call af_write_vtk(tree, 'bctest', dir="output")

      

     contains

      subroutine set_phi_zero(box)
              type(box_t), intent(inout) :: box
              box%cc(DTIMES(:), i_phi) = 0.0_dp
      end subroutine set_phi_zero

      subroutine boundary_method(box, nb, iv, coords, bc_val, bc_type)
              type(box_t), intent(in) :: box
              integer, intent(in) :: nb !neighbor indexing
              integer, intent(in) :: iv
              real(dp), intent(in) :: coords(NDIM, box%n_cell**(NDIM-1))
              real(dp), intent(out) :: bc_val(box%n_cell**(NDIM-1))
              integer, intent(out) :: bc_type
              integer :: nc


              nc = box%n_cell
              !coords may be used to set space dependent boundary conditions
              !print *, 'Coords of box ', box%ix, ": ", coords 

              select case (nb)

                case(af_neighb_highx)!af_neighb_highx == 2
                        bc_type = af_bc_neumann
                        bc_val = 0.0_dp

                case(af_neighb_lowx)!af_neighb_lowx == 1
                        bc_type = af_bc_dirichlet
                        bc_val = 1.0_dp
                
                case(af_neighb_lowy)!af_neighb_lowy == 3
                        bc_type = af_bc_dirichlet
                        bc_val = 1.0_dp
              
                case(af_neighb_highy)!af_neighb_highy == 4
                        bc_type = af_bc_neumann
                        bc_val = 0.0_dp
              end select

      end subroutine boundary_method



end program h_bcs
