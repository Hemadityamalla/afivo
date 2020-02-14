#include "../src/cpp_macros.h"

program h_advection
        
        use m_af_all
        implicit none

        integer, parameter :: ncells = 4
        integer :: iphi
        integer :: iphi_old
        integer :: ierr
        integer :: iflux
        integer, parameter :: coord_type = af_xyz
        character(len=100) :: fname

        real(dp) :: domain_len(NDIM)
        integer :: grid_size(NDIM)
        logical :: periodicBCs(NDIM)

        type(af_t) :: tree
        type(ref_info_t) :: refine_info
        integer :: time_steps, output_cnt
        integer :: i, n, n_steps
        real(dp) :: dt, time, end_time, err, sum_err2
        real(dp) :: sum_phi, sum_phi_t0
        real(dp) :: velocity(NDIM), dr_min(NDIM)

        domain_len(:) = 1.0_dp
        domain_len(1) = 2.0_dp

        grid_size(:) = 4*ncells
        grid_size(1) = 8*ncells

        periodicBCs(:) = .true.
        periodicBCs(2) = .false.

        print *, "Running advection_" // DIMNAME // "."
        print *, "Number of threads", af_get_max_threads()

        !Initializing the variables- based on the time stepping
        call af_add_cc_variable(tree, "phi", ix=iphi, n_copies=2)

        iphi_old = iphi+1
        call af_add_cc_variable(tree, "err", ix=ierr)
        call af_add_fc_variable(tree, "flux", ix=iflux)

        !Specifying the boundary conditions
        call af_set_cc_methods(tree, iphi, boundaries)

        call af_init(tree, ncells, domain_len, grid_size, &
             periodic=periodicBCs, coord=coord_type)



        output_cnt = 0
        time = 0.0_dp
        dt = 0.01_dp
        end_time = 2.0_dp
        velocity(:) = 0.0_dp
        velocity(1) = 1.0_dp
        
        !Setup initial conditions
        call af_loop_box(tree, set_ics)
        !Fill ghost cells
        call af_gc_tree(tree, [iphi])

        call af_loop_box(tree, computeCurvature)


        call af_write_silo(tree, "h_advection_2d_0", dir="output")
        !One level of refinement
        call af_adjust_refinement(tree, refine_routine, refine_info, 1)

        call af_restrict_tree(tree, iphi)

        call af_gc_tree(tree, [iphi])
        
        !Add the initial stuff
        call af_tree_sum_cc(tree, iphi, sum_phi_t0)

        call af_write_silo(tree, "h_advection_2d_1", dir="output")
        !Start integration
        n_steps = 1!int(end_time/dt)
        do i=1,n_steps
                call af_tree_copy_cc(tree, iphi, iphi_old)
                call  af_loop_boxes(tree, korenFlux)
                call af_loop_box_arg(tree, updateSolution, [dt])
                call af_gc_tree(tree, [iphi])
                if (mod(i,10) .eq. 0) then
                       write(fname, "(A,I0)") "h_advection_" // DIMNAME // "_", i
                        call af_write_silo(tree, trim(fname), dir="output")
                end if
                time = time+dt
        end do


        call af_destroy(tree)


        contains

         subroutine set_ics(box)
                 type(box_t), intent(inout) :: box
                 integer :: IJK, nc
                 real(dp) :: rr(NDIM)

                 nc = box%n_cell
                 do KJI_DO(1,nc)
                        !Obtaining coods of a cell center
                        rr = af_r_cc(box, [IJK])
                        box%cc(IJK, iphi) = solution(rr, 0.0_dp)
                  end do; CLOSE_DO
         end subroutine set_ics

         function solution(rr, t) result(sol)
                real(dp), intent(in) :: rr(NDIM), t
                real(dp) :: sol, rr_t(NDIM)

                rr_t = rr - velocity*t
                sol = 1.0_dp*exp(-100.0_dp*(rr(1) - 0.25_dp)**2 &
                                 -100.0_dp*(rr(2) - 0.5_dp)**2)
         end function solution
        
        subroutine korenFlux(boxes, id)
                use m_af_flux_schemes
                type(box_t), intent(inout) :: boxes(:)
                integer, intent(in) :: id
                integer :: nc
                real(dp), allocatable :: cc(DTIMES(:),:)
                real(dp), allocatable :: v(DTIMES(:), :)

                nc = boxes(id)%n_cell
                allocate(cc(DTIMES(-1:nc+2),1))
                allocate(v(DTIMES(1:nc+1),NDIM))
                call af_gc2_box(tree, id, [iphi], cc)

                v(:,:,1) = velocity(1)
                v(:,:,2) = velocity(2)

                call flux_koren_2d(cc(DTIMES(:), 1), v, nc, 2)
                tree%boxes(id)%fc(:,:,:, iphi) = v

        end subroutine korenFlux

        subroutine computeCurvature(box)
                type(box_t), intent(inout) :: box
                integer :: IJK, nc
               
                nc = box%n_cell 
                do KJI_DO(1,nc)
                        box%cc(IJK, ierr) = box%dr(1)**2*abs(box%cc(i+1,j,iphi)+box%cc(i-1,j,iphi)-2*box%cc(i,j,iphi)) + &
                               box%dr(2)**2*abs(box%cc(i,j+1,iphi)+box%cc(i,j-1,iphi)-2*box%cc(i,j,iphi))
                end do; CLOSE_DO
        end subroutine computeCurvature


        subroutine updateSolution(box, dt)
                type(box_t), intent(inout) :: box
                real(dp), intent(in) :: dt(:)
                real(dp) :: inv_dr(NDIM)
                integer :: IJK, nc

                nc = box%n_cell
                inv_dr = 1.0/box%dr

                do j=1,nc
                        do i=1,nc
                                box%cc(i,j,iphi) = box%cc(i,j,iphi)- dt(1)*( &
                                inv_dr(1)*( &
                                box%fc(i+1,j,1,iphi)-box%fc(i,j,1,iphi)) &
                                - inv_dr(2)*( &
                                box%fc(i,j+1,2,iphi)-box%fc(i,j,2,iphi)))
                        end do
                end do
        end subroutine updateSolution

        subroutine boundaries(box, nb, iv, coords, bc_val, bc_type)
                type(box_t), intent(in) :: box
                integer, intent(in) :: nb
                integer, intent(in) :: iv
                real(dp), intent(in) :: coords(NDIM, box%n_cell**(NDIM-1))
                real(dp), intent(out) :: bc_val(box%n_cell**(NDIM-1))
                integer, intent(out) :: bc_type
                integer :: nc

                nc = box%n_cell
                select case(nb)
                        case(af_neighb_lowy)
                                bc_type = af_bc_dirichlet
                                bc_val = 0.0_dp
                        case(af_neighb_highy)
                                bc_type = af_bc_dirichlet
                                bc_val = 0.0_dp
                end select

        end subroutine boundaries


        subroutine refine_routine(box, cell_flags)
                type(box_t), intent(in) :: box
                integer, intent(out) :: cell_flags(DTIMES(box%n_cell))
                real(dp) :: diff
                integer :: IJK, nc

                nc = box%n_cell

                cell_flags(:,:) = af_keep_ref
                !do KJI_DO(1,nc)
                !        if (box%lvl < 5 .and. box%ix(1) < 2) then
                !                cell_flags(IJK) = af_do_ref
                !        end if
                !end do; CLOSE_DO
                do KJI_DO(1,nc)
                       diff = box%dr(1)**2*abs(box%cc(i+1,j,iphi)+box%cc(i-1,j,iphi)-2*box%cc(i,j,iphi)) + &
                              box%dr(2)**2*abs(box%cc(i,j+1,iphi)+box%cc(i,j-1,iphi)-2*box%cc(i,j,iphi))
                       if (box%cc(i,j,ierr) > 1.0e-4_dp .and. box%lvl < 5) then
                                cell_flags(IJK) = af_do_ref
                               !print *, "Marked for refining"
                        else if (diff < 0.1_dp*1.0e-1_dp) then
                                cell_flags(IJK) = af_rm_ref
                               !print *, "Removed refining"
                        else
                                cell_flags(IJK) = af_keep_ref
                        end if
                end do; CLOSE_DO

        end subroutine refine_routine

end program h_advection
