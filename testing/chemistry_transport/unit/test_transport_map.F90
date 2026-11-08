! Unit tests for the MICM gas transport map (P1): which runtime chemistry
! scalars are transported independently, in what order, and which inputs fail.
program test_transport_map

   use unit_assert, only: check, finish
   use mpas_chem_transport_map

   implicit none

   integer, allocatable :: map(:)
   integer :: n, err, bad, i
   integer :: idx64(64)

   ! Chemistry group occupies scalars 7..9 unless a case says otherwise.

   call run([integer ::], [logical ::], [logical ::], 9, 7, 9)
   call check(err == MAP_OK .and. n == 0, 'empty mechanism is a valid no-op')

   call run([7], [.false.], [.true.], 7, 7, 7)
   call check(err == MAP_OK .and. n == 1, 'one species: count')
   call check(map_is([7]), 'one species: index')

   call run([7, 8, 9], [.false., .false., .false.], [.true., .true., .true.], 9, 7, 9)
   call check(err == MAP_OK .and. map_is([7, 8, 9]), 'three species in input order')

   idx64 = [(7 + i, i = 1, 64)]
   call run(idx64, spread(.false., 1, 64), spread(.true., 1, 64), 71, 8, 71)
   call check(err == MAP_OK .and. n == 64, '64 species: count')
   call check(map_is(idx64), '64 species: every index once, in order')

   call run([9, 7, 8], [.false., .false., .false.], [.true., .true., .true.], 9, 7, 9)
   call check(err == MAP_OK .and. map_is([9, 7, 8]), 'reordered mechanism keeps input order')

   ! Group 7..14; scalars 8, 10, 12, 13 are chemistry scalars not in this list.
   call run([14, 9, 11], [.false., .false., .false.], [.true., .true., .true.], 14, 7, 14)
   call check(err == MAP_OK .and. map_is([14, 9, 11]), 'noncontiguous indices interleaved with unmapped ones')

   ! H2O, O2 and N2 bound to the native qv index (1); NO and NO2 transported.
   call run([1, 7, 1, 8, 1], [.true., .false., .true., .false., .true.], &
            [.false., .true., .false., .true., .false.], 9, 7, 9)
   call check(err == MAP_OK .and. map_is([7, 8]), 'host-bound H2O, O2, N2 are excluded')

   call run([7, 8], [.false., .false.], [.true., .false.], 9, 7, 9)
   call check(err == MAP_OK .and. map_is([7]), 'species without write-back is excluded')

   call run([7, 0], [.false., .false.], [.true., .true.], 9, 7, 9)
   call check(err == MAP_OUT_OF_RANGE .and. bad == 2 .and. n == 0, 'unresolved index 0 is out of range')

   call run([7, 10], [.false., .false.], [.true., .true.], 9, 7, 9)
   call check(err == MAP_OUT_OF_RANGE .and. bad == 2, 'index above num_scalars is out of range')

   call run([7, 1], [.false., .false.], [.true., .true.], 9, 7, 9)
   call check(err == MAP_NOT_CHEMISTRY .and. bad == 2 .and. n == 0, 'native scalar index is rejected')

   call run([7], [.false.], [.true.], 9, 0, -1)
   call check(err == MAP_NOT_CHEMISTRY .and. bad == 1, 'missing chemistry group is rejected')

   call run([7, 8, 7], [.false., .false., .false.], [.true., .true., .true.], 9, 7, 9)
   call check(err == MAP_DUPLICATE .and. bad == 3 .and. n == 0, 'duplicate index is rejected at its second use')

   call check(len_trim(transport_map_error_text(MAP_OUT_OF_RANGE)) > 0 .and. &
              transport_map_error_text(MAP_OUT_OF_RANGE) /= transport_map_error_text(MAP_NOT_CHEMISTRY) .and. &
              transport_map_error_text(MAP_NOT_CHEMISTRY) /= transport_map_error_text(MAP_DUPLICATE), &
              'error texts are distinct')

   call finish('test_transport_map')

contains

   subroutine run(mpas_index, host_bound, write_back, num_scalars, chem_first, chem_last)
      integer, intent(in) :: mpas_index(:)
      logical, intent(in) :: host_bound(:), write_back(:)
      integer, intent(in) :: num_scalars, chem_first, chem_last

      if (allocated(map)) deallocate(map)
      allocate(map(size(mpas_index)))
      map(:) = -999
      call build_transport_map(mpas_index, host_bound, write_back, num_scalars, &
                               chem_first, chem_last, map, n, err, bad)
   end subroutine run

   logical function map_is(expected)
      integer, intent(in) :: expected(:)

      map_is = (n == size(expected))
      if (map_is) map_is = all(map(1:n) == expected)
   end function map_is

end program test_transport_map
