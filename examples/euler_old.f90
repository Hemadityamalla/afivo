#include "../src/cpp_macros.h"

program euler
  use m_af_all
  implicit none
  
  integer, parameter :: ncells = 8
  integer :: ic1, ic2, ic3, ic4
  integer :: if1x, if1y, if2x, if2y, if3x, if3y, if4x, if4y
  integer, parameter :: coord_type = af_xyz
  real(dp) :: l_max(NDIM), l_min(NDIM)
  integer :: grid(NDIM)
  logical :: periodicBC(NDIM)
  real(dp), parameter :: Y = 1.4_dp
  real(dp):: p(4), rho(4), u(4), v(4)
  real(dp) :: wSp, dr_min, error, test
  
  type(af_t) :: tree
  real(dp) :: dt, time, end_time
  integer :: t_iter
  character(len=100) :: fname
  print *, "Running Euler 2D"
  print *, "Number of threads", af_get_max_threads()
  
  !p = (/1.0_dp, 0.4_dp, 0.0439_dp, 0.15_dp/)
  !rho = (/1.0_dp, 0.5197_dp, 0.1072_dp, 0.2579_dp/)
  !u = (/0.0_dp, -0.7259_dp, -0.7259_dp, 0.0_dp/)
  !v = (/0.0_dp, 0.0_dp, -1.4045_dp, -1.4045_dp/)

  rho = (/0.125_dp, 1.0_dp, 1.0_dp, 0.125_dp/)
  p = (/0.1_dp, 1.0_dp, 1.0_dp, 0.1_dp/)
  u = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
  v = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

  
  wSp = calc_speed(rho, u, v, p)
  
  
  
  grid(:) = 50*ncells
  l_max(:) = 1.0_dp
  l_min(:) = 0.0_dp
  periodicBC(:) = .false.
  
  call af_add_cc_variable(tree, "rho", ix=ic1)
  call af_add_cc_variable(tree, "rhoU", ix=ic2)
  call af_add_cc_variable(tree, "rhoV", ix=ic3)
  call af_add_cc_variable(tree, "E", ix=ic4)

  call af_add_fc_variable(tree, "flux1x", ix=if1x)
  call af_add_fc_variable(tree, "flux1y", ix=if1y)
  call af_add_fc_variable(tree, "flux2x", ix=if2x)
  call af_add_fc_variable(tree, "flux2y", ix=if2y)
  call af_add_fc_variable(tree, "flux3x", ix=if3x)
  call af_add_fc_variable(tree, "flux3y", ix=if3y)
  call af_add_fc_variable(tree, "flux4x", ix=if4x)
  call af_add_fc_variable(tree, "flux4y", ix=if4y)
  
  call af_set_cc_methods(tree, ic1, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic2, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic3, af_bc_neumann_zero)
  call af_set_cc_methods(tree, ic4, af_bc_neumann_zero)
  
  call af_init(tree, ncells, l_max, grid, &
               periodic = periodicBC, r_min= l_min, &
               coord = coord_type)
               
  call af_loop_box(tree, setInitConds)
  
  call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
  
  !call af_write_silo(tree, 'eulerInit', dir='output')
  
  !Setting the timestep data
  time = 0.0_dp
  end_time = 0.2_dp
  t_iter = 0
  !dt = 0.01
  do 
   dr_min = af_min_dr(tree)
   dt = 0.1_dp/(wSp/dr_min) + epsilon(1.0_dp)
    if (mod(t_iter, 10) == 0) then
      write(fname, "(A,I0)") "euler_" // DIMNAME // "_", t_iter
      call af_write_silo(tree, trim(fname), t_iter, time, dir="output")
    end if
  
   call af_loop_tree(tree, korenFlux)
   call af_loop_box_arg(tree, updateSoln, [dt])
   call af_gc_tree(tree, [ic1, ic2, ic3, ic4])
    
   t_iter = t_iter + 1
   time = time + dt
   
   
   call af_tree_maxabs_cc(tree, ic1, error)
   call af_tree_max_fc(tree, 1, if4x, test)
   print *, "Max value of fc: ", test
   if (error > 2.0_dp) then
    print *, "solution diverging!"
    exit
   end if
   
   if (time > end_time) exit
  end do
  
  
  
  
  call af_destroy(tree)
  !=====================================================================
  contains
  
  !=====================================================================
  
  function calc_speed(rho, u, v, p) result(waveSpeed)
    real(dp), intent(in) :: rho(4), u(4), v(4), p(4)
    real(dp):: waveSpeed
    integer :: i
    
    waveSpeed = 0.0_dp
    do i=1,4
      waveSpeed = max(waveSpeed, abs(sqrt((Y*p(i))/(rho(i))) + &
                                        sqrt(u(i)**2 + v(i)**2)))
    end do 

  end function calc_speed
  !=====================================================================
  function convert_to_conservatives( rho, u, v, p ) result(c)
    real(dp), intent(in) :: rho, u, v, p
    real(dp):: c(4)
    
    c(1) = rho
    c(2) = rho*u
    c(3) = rho*v
    c(4) = (p/(Y-1.0_dp)) + 0.5_dp*rho*(u**2 + v**2)
    
  end function convert_to_conservatives
  !=====================================================================
  
  subroutine setInitConds( box )
    type(box_t), intent(inout) :: box
    integer :: IJK, nc
    real(dp) :: rr(NDIM)
    !real(dp) :: p(4), rho(4), u(4), v(4)
    real(dp) :: conservatives(4)
    
    
    nc = box%n_cell
    do KJI_DO(0, nc+1)
      rr = af_r_cc(box, [IJK])
      if (rr(1) > 0.5_dp .and. rr(2) > 0.5_dp) then
        conservatives = convert_to_conservatives(rho(1), u(1), v(1), p(1))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .ge. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(2), u(2), v(2), p(2))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
      elseif (rr(1) .le. 0.5_dp .and. rr(2) .le. 0.5_dp) then
        conservatives = convert_to_conservatives(rho(3), u(3), v(3), p(3))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
      else
        conservatives = convert_to_conservatives(rho(4), u(4), v(4), p(4))
        box%cc(IJK, ic1) = conservatives(1)
        box%cc(IJK, ic2) = conservatives(2)
        box%cc(IJK, ic3) = conservatives(3)
        box%cc(IJK, ic4) = conservatives(4)
      end if
    end do; CLOSE_DO
  end subroutine setInitConds
  !=====================================================================
  subroutine korenFlux( tree, id )
    use m_af_flux_schemes
    type(af_t), intent(inout) :: tree
    integer, intent(in) :: id
    integer :: nc
    real(dp), allocatable :: cc(DTIMES(:), :)
    real(dp), allocatable :: fc1(DTIMES(:), :), &
                             fc2(DTIMES(:), :), &
                             fc3(DTIMES(:), :), &
                             fc4(DTIMES(:), :)
    
    nc = tree%boxes(id)%n_cell
    allocate(cc(DTIMES(-1:nc+2), 4))
    allocate(fc1(DTIMES(1:nc+1), NDIM))
    allocate(fc2(DTIMES(1:nc+1), NDIM))
    allocate(fc3(DTIMES(1:nc+1), NDIM))
    allocate(fc4(DTIMES(1:nc+1), NDIM))
    
    call af_gc2_box(tree, id, [ic1, ic2, ic3, ic4], cc)
    fc1 = 1.0_dp
    fc2 = 1.0_dp
    fc3 = 1.0_dp
    fc4 = 1.0_dp
    
    !There are some redundant calculations being done here.
    !For instance, I do not need the values of fc(:,:,3,ifx)....
    
    call flux_koren_2d(cc(DTIMES(:), ic1), fc1, nc, 2)
    call flux_koren_2d(cc(DTIMES(:), ic2), fc2, nc, 2)
    call flux_koren_2d(cc(DTIMES(:), ic3), fc3, nc, 2)
    call flux_koren_2d(cc(DTIMES(:), ic4), fc4, nc, 2)
    
    tree%boxes(id)%fc(:,:,:,if1x) = fc2
    tree%boxes(id)%fc(:,:,:,if1y) = fc3
    
    
    tree%boxes(id)%fc(:,:,:,if2x) = (fc2**2/fc1) + &
                                    (Y-1.0_dp)*( fc4 - &
                                    (fc2**2 + fc3**2)/ &
                                    (2.0_dp*fc1))
    tree%boxes(id)%fc(:,:,:,if2y) = (fc2*fc3)/fc1
    
    
    tree%boxes(id)%fc(:,:,:,if3x) = (fc2*fc3)/fc1
    tree%boxes(id)%fc(:,:,:,if3y) = (fc3**2/fc1) + &
                                    (Y-1.0_dp)*( fc4 - &
                                    (fc2**2 + fc3**2)/ &
                                    (2.0_dp*fc1))
                                    
                                    
    tree%boxes(id)%fc(:,:,:,if4x) = (fc2/fc1)*(fc4 + &
                                     (Y - 1.0_dp)*( fc4 - &
                                     (fc2**2 + fc3**2)/ &
                                     (2.0_dp*fc1)))
                                     
    tree%boxes(id)%fc(:,:,:,if4y) = (fc3/fc1)*(fc4 + &
                                     (Y - 1.0_dp)*( fc4 - &
                                     (fc2**2 + fc3**2)/ &
                                     (2.0_dp*fc1)))

    
  end subroutine korenFlux
  !=====================================================================
  
  subroutine updateSoln( box, dt )
    type(box_t), intent(inout) :: box
    real(dp), intent(in) :: dt(:)
    real(dp) :: inv_dr(NDIM)
    integer ::IJK, nc
    nc = box%n_cell
    inv_dr = 1.0_dp/box%dr
  
    do j=1,nc
      do i=1,nc
        box%cc(i,j,ic1) = box%cc(i,j,ic1) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if1x) - box%fc(i,j,1,if1x)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if1y) - box%fc(i,j,2,if1y)))
       box%cc(i,j,ic2) = box%cc(i,j,ic2) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if2x) - box%fc(i,j,1,if2x)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if2y) - box%fc(i,j,2,if2y)))
       
       box%cc(i,j,ic3) = box%cc(i,j,ic3) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if3x) - box%fc(i,j,1,if3x)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if3y) - box%fc(i,j,2,if3y)))
       
       box%cc(i,j,ic4) = box%cc(i,j,ic4) - dt(1) *( &
                          inv_dr(1)*( &
                          box%fc(i+1,j,1,if4x) - box%fc(i,j,1,if4x)) + &
                          inv_dr(2)*( &
                          box%fc(i,j+1,2,if4y) - box%fc(i,j,2,if4y)))
                          
        
      end do
    end do
  end subroutine updateSoln
  

end program euler
