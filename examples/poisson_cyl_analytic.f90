!> \example poisson_cyl_analytic_analytic.f90
!>
!> Example showing how to use multigrid and compare with an analytic solution. A
!> standard 5-point Laplacian is used in cylindrical coordinates.
program poisson_cyl_analytic
  use m_af_all

  implicit none

  integer, parameter :: box_size     = 8
  integer, parameter :: n_iterations = 10
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_err
  integer            :: i_sol
  integer            :: i_tmp

  real(dp), parameter :: domain_len = 1.25e-2_dp
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: sigma = 4e-4_dp * sqrt(0.5_dp)
  real(dp), parameter :: rz_source(2) = [0.0_dp, 0.5_dp] * domain_len
  real(dp), parameter :: epsilon_source = 8.85e-12_dp
  real(dp), parameter :: Q_source = 3e18_dp * 1.6022e-19_dp * &
       sigma**3 * sqrt(2 * pi)**3

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: mg_iter
  real(dp)           :: dr, residu(2), anal_err(2), max_val
  character(len=100) :: fname
  type(mg_t)        :: mg
  integer            :: count_rate,t_start, t_end

  print *, "Running poisson_cyl_analytic"
  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 1.25e-2_dp / box_size

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "err", ix=i_err)
  call af_add_cc_variable(tree, "sol", ix=i_sol)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       dr, &           ! Distance between cells on base level
       coord=af_cyl)   ! Cylindrical coordinates

  call af_set_coarse_grid(tree, [box_size, box_size])
  call af_print_info(tree)

  call system_clock(t_start, count_rate)
  do
     ! For each box, set the initial conditions
     call af_loop_box(tree, set_init_cond)

     ! This updates the refinement of the tree, by at most one level per call.
     call af_adjust_refinement(tree, ref_routine, ref_info)

     ! If no new boxes have been added, exit the loop
     if (ref_info%n_add == 0) exit
  end do
  call system_clock(t_end, count_rate)

  write(*,"(A,Es10.3,A)") " Wall-clock time generating AMR grid: ", &
       (t_end-t_start) / real(count_rate,dp), " seconds"

  call af_print_info(tree)

  ! Set the multigrid options.
  mg%i_phi        = i_phi       ! Solution variable
  mg%i_rhs        = i_rhs       ! Right-hand side variable
  mg%i_tmp        = i_tmp       ! Variable for temporary space
  mg%sides_bc     => sides_bc   ! Method for boundary conditions Because we use

  ! Automatically detect the right methods
  mg%box_op       => mg_auto_op
  mg%box_gsrb     => mg_auto_gsrb
  mg%box_corr     => mg_auto_corr

  ! Initialize the multigrid options. This performs some basics checks and sets
  ! default values where necessary.
  ! This routine does not initialize the multigrid variables i_phi, i_rhs
  ! and i_tmp. These variables will be initialized at the first call of mg_fas_fmg
  call mg_init(tree, mg)

  print *, "Multigrid iteration | max residual | max rel. error"
  call system_clock(t_start, count_rate)

  do mg_iter = 1, n_iterations
     ! Perform a FAS-FMG cycle (full approximation scheme, full multigrid). The
     ! third argument controls whether the residual is stored in i_tmp. The
     ! fourth argument controls whether to improve the current solution.
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)

     ! Compute the error compared to the analytic solution
     call af_loop_box(tree, set_err)

     ! Determine the minimum and maximum residual and error
     call af_tree_min_cc(tree, i_tmp, residu(1))
     call af_tree_max_cc(tree, i_tmp, residu(2))
     call af_tree_min_cc(tree, i_err, anal_err(1))
     call af_tree_max_cc(tree, i_err, anal_err(2))
     call af_tree_max_cc(tree, i_phi, max_val)
     anal_err = anal_err / max_val

     write(*,"(I8,2Es14.5)") mg_iter, maxval(abs(residu)), &
          maxval(abs(anal_err))

     write(fname, "(A,I0)") "poisson_cyl_analytic_", mg_iter
     call af_write_vtk(tree, trim(fname), dir="output")
  end do
  call system_clock(t_end, count_rate)

  write(*, "(A,I0,A,E10.3,A)") &
       " Wall-clock time after ", n_iterations, &
       " iterations: ", (t_end-t_start) / real(count_rate, dp), &
       " seconds"

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call af_destroy(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: crv

    nc = box%n_cell

    ! Compute the "curvature" in phi
    do j = 1, nc
       do i = 1, nc
          crv = box%dr**2 * abs(box%cc(i, j, i_rhs))

          ! And refine if it exceeds a threshold
          if (crv > 1e-1_dp .and. box%lvl < 10) then
             cell_flags(i, j) = af_do_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          rz = af_r_cc(box, [i,j]) - rz_source

          ! Gaussian source term
          box%cc(i, j, i_rhs) = - Q_source * &
               exp(-sum(rz**2) / (2 * sigma**2)) / &
               (sigma**3 * sqrt(2 * pi)**3 * epsilon_source)
       end do
    end do
  end subroutine set_init_cond

  real(dp) function solution(rz_in)
    real(dp), intent(in) :: rz_in(2)
    real(dp) :: d, erf_d

    d = norm2(rz_in - rz_source)

    ! Take limit when d -> 0
    if (d < sqrt(epsilon(1.0d0))) then
       erf_d = sqrt(2.0_dp/pi) / sigma
    else
       erf_d = erf(d * sqrt(0.5_dp) / sigma) / d
    end if

    solution = erf_d * Q_source / (4 * pi * epsilon_source)
  end function solution

  ! Compute error compared to the analytic solution
  subroutine set_err(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          rz = af_r_cc(box, [i,j])
          box%cc(i, j, i_sol) = solution(rz)
          box%cc(i, j, i_err) = box%cc(i, j, i_phi) - box%cc(i, j, i_sol)
       end do
    end do
  end subroutine set_err

  ! This routine sets boundary conditions for a box, by filling its ghost cells
  ! with approriate values. Note that on the axis (a boundary in the lower-x
  ! direction) we should use a Neumann zero condition in cylindrical
  ! coordinates.
  subroutine sides_bc(box, nb, iv, bc_type)
    type(box_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    real(dp)                    :: rz(2)
    integer                     :: n, nc

    nc = box%n_cell

    select case (nb)
    case (af_neighb_lowx)             ! Neumann zero on axis
       bc_type = af_bc_neumann
       box%cc(0, 1:nc, iv) = 0
    case (af_neighb_highx)             ! Use solution on other boundaries
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = af_rr_cc(box, [nc+0.5_dp, real(n, dp)])
          box%cc(nc+1, n, iv) = solution(rz)
       end do
    case (af_neighb_lowy)
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = af_rr_cc(box, [real(n, dp), 0.5_dp])
          box%cc(n, 0, iv) = solution(rz)
       end do
    case (af_neighb_highy)
       bc_type = af_bc_dirichlet
       do n = 1, nc
          rz = af_rr_cc(box, [real(n, dp), nc+0.5_dp])
          box%cc(n, nc+1, iv) = solution(rz)
       end do
    end select
  end subroutine sides_bc

end program poisson_cyl_analytic
