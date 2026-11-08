! Transport core with processes that are not enabled (the shared stubs, or a real process module
! that was not initialized): assembly and tendency addition must leave every tendency bit unchanged.
program test_tracer_transport_stubs

   use mpas_kind_types, only: RKIND, I8KIND
   use unit_assert, only: check, finish
   use mpas_atmphys_tracer_transport
   use mpas_atmphys_pbl_gas, only: pbl_gas_enabled
   use mpas_atmphys_cu_gas, only: cu_gas_enabled

   implicit none

   integer, parameter :: ns = 6, nz = 3, nc = 4
   real(kind=RKIND) :: scalars(ns,nz,nc), tend(ns,nz,nc), tend0(ns,nz,nc)
   real(kind=RKIND) :: zz(nz,nc), rho_zz(nz,nc), zgrid(nz+1,nc), area(nc), mass(nz,nc)
   type(transport_budget_t) :: b
   integer :: s, k, i

   do i = 1, nc
      do k = 1, nz
         do s = 1, ns
            scalars(s,k,i) = 1.0e-9_RKIND*(s + 10*k + 100*i)
            tend0(s,k,i) = 1.0e-6_RKIND*(s - 3*k + 7*i)
         end do
      end do
   end do
   tend0(1,1,1) = -0.0_RKIND
   zz = 1.0_RKIND
   rho_zz = 1.0_RKIND
   do i = 1, nc
      zgrid(:,i) = [0.0_RKIND, 100.0_RKIND, 300.0_RKIND, 600.0_RKIND]
   end do
   area = 1.0_RKIND
   mass = 2.0_RKIND

   call check(.not. pbl_gas_enabled() .and. .not. cu_gas_enabled(), 'processes are not enabled before initialization')

   call tracer_transport_setup([5, 3], nz, nc)
   call check(tracer_transport_n() == 2, 'map is stored')
   call check(.not. tracer_transport_enabled(), 'stub processes leave transport disabled')
   call check(tracer_transport_buffer_size() == 0, 'disabled transport allocates no tendency buffers')

   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, 3.0_RKIND)
   call tracer_transport_add_tend(mass, tend)
   call check(bitwise_equal(tend, tend0), 'disabled transport leaves every tendency bit unchanged')
   b = tracer_transport_budget(TRANSPORT_PBL, 1)
   call check(b%active_columns == 0, 'disabled transport records no active columns')

   call tracer_transport_teardown()
   call tracer_transport_setup([integer ::], nz, nc)
   call check(tracer_transport_n() == 0 .and. .not. tracer_transport_enabled(), 'empty map is a no-op')
   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, 3.0_RKIND)
   call tracer_transport_add_tend(mass, tend)
   call check(bitwise_equal(tend, tend0), 'empty map leaves tendencies unchanged')

   call tracer_transport_teardown()
   call tracer_transport_setup([4], nz, nc)
   call check(tracer_transport_n() == 1, 'setup after teardown starts clean')
   call tracer_transport_teardown()
   call check(tracer_transport_n() == 0, 'teardown resets the map')

   call finish('test_tracer_transport_stubs')

contains

   logical function bitwise_equal(a, b)
      real(kind=RKIND), intent(in) :: a(:,:,:), b(:,:,:)

      bitwise_equal = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
   end function bitwise_equal

end program test_tracer_transport_stubs
