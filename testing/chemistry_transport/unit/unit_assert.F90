! Minimal assertion helpers for the standalone transport unit tests.
module unit_assert

   implicit none
   private
   public :: check, finish

   integer, save :: n_checks = 0
   integer, save :: n_failed = 0

contains

   subroutine check(condition, name)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name

      n_checks = n_checks + 1
      if (condition) then
         write(*, '(a,a)') 'PASS ', name
      else
         n_failed = n_failed + 1
         write(*, '(a,a)') 'FAIL ', name
      end if
   end subroutine check

   subroutine finish(suite)
      character(len=*), intent(in) :: suite

      write(*, '(a,a,i0,a,i0,a)') trim(suite), ': ', n_checks - n_failed, ' of ', n_checks, ' checks passed'
      if (n_failed > 0 .or. n_checks == 0) error stop 1
   end subroutine finish

end module unit_assert
