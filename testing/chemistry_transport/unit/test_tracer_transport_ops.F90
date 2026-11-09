! Transport core driven by test-double processes (P1). The doubles in
! doubles/test_tracer_transport_ops/ replace the production stubs:
!   PBL: active in odd cells; moves half of layer-1 mass into layer 2.
!   CU : active in cells 1-2; species 1 moves layer 2 into layer 3,
!        species 2 doubles its top layer.
! Expected values are written out by hand for the column below.
program test_tracer_transport_ops

   use mpas_kind_types, only: RKIND, I8KIND
   use unit_assert, only: check, finish
   use mpas_atmphys_tracer_transport
   use mpas_atmphys_pbl_gas, only: pbl_double_active
   use mpas_atmphys_cu_gas, only: cu_double_active

   implicit none

   integer, parameter :: ns = 6, nz = 3, nc = 4
   integer, parameter :: SP1 = 5, SP2 = 3          ! map(1) = 5, map(2) = 3: noncontiguous, reversed
   real(kind=RKIND), parameter :: dt = 3.0_RKIND
   real(kind=RKIND), parameter :: tol = 1.0e-12_RKIND
   real(kind=RKIND) :: scalars(ns,nz,nc), scalars0(ns,nz,nc)
   real(kind=RKIND) :: tend(ns,nz,nc), tend0(ns,nz,nc)
   real(kind=RKIND) :: zz(nz,nc), rho_zz(nz,nc), zgrid(nz+1,nc), area(nc), mass(nz,nc)
   real(kind=RKIND) :: q2(2,nz,nc), expect(nz)
   type(transport_budget_t) :: b
   integer :: s, k, i

   ! m_d = zz*rho_zz*dz = [100, 160, 150] kg m-2 in every cell.
   zz = 1.0_RKIND
   do i = 1, nc
      rho_zz(:,i) = [1.0_RKIND, 0.8_RKIND, 0.5_RKIND]
      zgrid(:,i) = [0.0_RKIND, 100.0_RKIND, 300.0_RKIND, 600.0_RKIND]
   end do
   area = [1.0_RKIND, 2.0_RKIND, 3.0_RKIND, 4.0_RKIND]
   mass = 2.0_RKIND

   do i = 1, nc
      do k = 1, nz
         do s = 1, ns
            scalars(s,k,i) = 7.0e-6_RKIND*(s + 10*k + 100*i)
            tend0(s,k,i) = 1.0e-6_RKIND*(s - 3*k + 7*i)
         end do
      end do
      scalars(SP1,:,i) = [4.0e-9_RKIND, 2.0e-9_RKIND, 1.0e-9_RKIND]
      scalars(SP2,:,i) = [3.0e-9_RKIND, 3.0e-9_RKIND, 3.0e-9_RKIND]
   end do
   scalars0 = scalars

   ! Hand-computed q2 = CU(PBL(q0)) for each cell.
   q2(1,:,1) = [2.0e-9_RKIND, 0.0_RKIND, 1.0e-9_RKIND + 3.25e-9_RKIND*160.0_RKIND/150.0_RKIND]
   q2(2,:,1) = [1.5e-9_RKIND, 3.9375e-9_RKIND, 6.0e-9_RKIND]
   q2(1,:,2) = [4.0e-9_RKIND, 0.0_RKIND, 1.0e-9_RKIND + 2.0e-9_RKIND*160.0_RKIND/150.0_RKIND]
   q2(2,:,2) = [3.0e-9_RKIND, 3.0e-9_RKIND, 6.0e-9_RKIND]
   q2(1,:,3) = [2.0e-9_RKIND, 3.25e-9_RKIND, 1.0e-9_RKIND]
   q2(2,:,3) = [1.5e-9_RKIND, 3.9375e-9_RKIND, 3.0e-9_RKIND]
   q2(1,:,4) = scalars0(SP1,:,4)
   q2(2,:,4) = scalars0(SP2,:,4)

   call tracer_transport_setup([SP1, SP2], nz, nc)
   call check(tracer_transport_enabled(), 'bound processes enable transport')
   call check(tracer_transport_buffer_size() == 2*2*nz*nc, 'enabled transport allocates both tendency buffers')

   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, dt)
   call check(bitwise_equal(scalars, scalars0), 'assembly never writes the prognostic scalars')

   call tracer_transport_add_tend(mass, tend)

   do i = 1, nc
      expect = tend0(SP1,:,i) + (q2(1,:,i) - scalars0(SP1,:,i))/dt*mass(:,i)
      call check(all(abs(tend(SP1,:,i) - expect) <= tol*maxval(abs(expect))), 'species 1 tendency, cell '//char(48+i))
      expect = tend0(SP2,:,i) + (q2(2,:,i) - scalars0(SP2,:,i))/dt*mass(:,i)
      call check(all(abs(tend(SP2,:,i) - expect) <= tol*maxval(abs(expect))), 'species 2 tendency, cell '//char(48+i))
   end do
   call check(bitwise_equal(tend(:,:,4:4), tend0(:,:,4:4)), 'inactive cell keeps every tendency bit')
   call check(bitwise_equal(tend(1:2,:,:), tend0(1:2,:,:)) .and. bitwise_equal(tend(4:4,:,:), tend0(4:4,:,:)) &
              .and. bitwise_equal(tend(6:6,:,:), tend0(6:6,:,:)), 'unmapped scalars keep every tendency bit')

   b = tracer_transport_budget(TRANSPORT_PBL, 1)
   call check(b%active_columns == 2 .and. b%residual <= 1.0e-12_RKIND, 'PBL budget: two conservative columns')
   call check(abs(b%increment) <= 1.0e-12_RKIND*b%mass_before .and. b%dt_applied == dt, 'PBL budget: no increment, dt recorded')
   b = tracer_transport_budget(TRANSPORT_CU, 1)
   call check(b%active_columns == 2 .and. b%residual <= 1.0e-12_RKIND, 'CU budget species 1: conservative')
   b = tracer_transport_budget(TRANSPORT_CU, 2)
   ! doubling the top layer adds area*m_d(3)*q1(3) in cells 1 (area 1) and 2 (area 2): 150*3e-9*(1+2)
   call check(abs(b%increment - 1.35e-6_RKIND) <= 1.0e-12_RKIND*1.35e-6_RKIND, 'CU budget species 2: source is measured')
   ! per column: added 150*3e-9 over an initial 100*q(1) + 160*q(2) + 150*q(3) = 1.23e-6 in both cells
   call check(abs(b%residual - 4.5e-7_RKIND/1.23e-6_RKIND) <= 1.0e-12_RKIND, 'CU budget species 2: residual flags the source')

   ! A second assembly replaces, rather than accumulates, the previous tendencies and budgets.
   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, dt)
   call tracer_transport_add_tend(mass, tend)
   expect = tend0(SP1,:,1) + (q2(1,:,1) - scalars0(SP1,:,1))/dt*mass(:,1)
   call check(all(abs(tend(SP1,:,1) - expect) <= tol*maxval(abs(expect))), 'repeated assembly is not cumulative')
   b = tracer_transport_budget(TRANSPORT_PBL, 1)
   call check(b%active_columns == 2, 'budgets restart each assembly')

   ! Only owned cells (1:nCellsSolve) are assembled and updated, even where a process would act: with
   ! nCellsSolve = 2, cell 3 (odd, so the PBL double is active there) and cell 4 keep every tendency bit.
   call tracer_transport_teardown()
   call tracer_transport_setup([SP1, SP2], nz, 2)
   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, dt)
   call tracer_transport_add_tend(mass, tend)
   call check(bitwise_equal(tend(:,:,3:4), tend0(:,:,3:4)), 'cells beyond nCellsSolve keep every tendency bit')
   expect = tend0(SP1,:,1) + (q2(1,:,1) - scalars0(SP1,:,1))/dt*mass(:,1)
   call check(all(abs(tend(SP1,:,1) - expect) <= tol*maxval(abs(expect))), 'owned cells are still updated')

   ! A rank whose columns are all inactive (P5 MPI decomposition): no tendency is touched and the budgets are
   ! empty, with nothing divided by its zero mass.
   call tracer_transport_teardown()
   call tracer_transport_setup([SP1, SP2], nz, nc)
   pbl_double_active = .false.
   cu_double_active = .false.
   tend = tend0
   call tracer_transport_assemble_columns(scalars, zz, rho_zz, zgrid, area, dt)
   call tracer_transport_add_tend(mass, tend)
   call check(bitwise_equal(tend, tend0), 'an all-inactive rank keeps every tendency bit')
   b = tracer_transport_budget(TRANSPORT_PBL, 1)
   call check(b%active_columns == 0 .and. b%mass_before == 0.0_RKIND .and. b%mass_after == 0.0_RKIND .and. &
              b%increment == 0.0_RKIND .and. b%residual == 0.0_RKIND, 'an all-inactive rank has an empty PBL budget')
   b = tracer_transport_budget(TRANSPORT_CU, 2)
   call check(b%active_columns == 0 .and. b%increment == 0.0_RKIND, 'an all-inactive rank has an empty CU budget')
   pbl_double_active = .true.
   cu_double_active = .true.

   call tracer_transport_teardown()
   call finish('test_tracer_transport_ops')

contains

   logical function bitwise_equal(a, b)
      real(kind=RKIND), intent(in) :: a(:,:,:), b(:,:,:)

      bitwise_equal = all(transfer(a, [0_I8KIND]) == transfer(b, [0_I8KIND]))
   end function bitwise_equal

end program test_tracer_transport_ops
