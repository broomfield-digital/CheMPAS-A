! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!-----------------------------------------------------------------------
!  mpas_chem_transport_map
!
!> \brief Select the MICM gases that physics transports independently
!> \details
!>  Builds the compact list of MPAS scalar indices for every prognostic MICM
!>  species that owns its own runtime chemistry tracer. Host-bound species
!>  (H2O read from qv, O2/N2 anchored on qv) and species without write-back
!>  are excluded. Nothing is filtered by name, ordering, contiguity, initial
!>  value or chemical lifetime. An empty result is a valid no-op.
!>
!>  The kernel has no MPAS or MICM dependencies so it can be unit tested;
!>  mpas_musica supplies the species table and logs the result.
!-----------------------------------------------------------------------
module mpas_chem_transport_map

    implicit none
    private

    public :: build_transport_map, transport_map_error_text
    public :: MAP_OK, MAP_OUT_OF_RANGE, MAP_NOT_CHEMISTRY, MAP_DUPLICATE

    integer, parameter :: MAP_OK            = 0
    integer, parameter :: MAP_OUT_OF_RANGE  = 1   ! index < 1 or > num_scalars
    integer, parameter :: MAP_NOT_CHEMISTRY = 2   ! index outside the runtime chemistry group
    integer, parameter :: MAP_DUPLICATE     = 3   ! index already claimed by an earlier species

contains

    !-----------------------------------------------------------------------
    !  routine build_transport_map
    !
    !> \brief Compact map from transported gas n to its MPAS scalar index
    !> \details
    !>  Species are taken in input order. Every transported index must lie in
    !>  the runtime chemistry constituent group [chem_first, chem_last]; an
    !>  index outside it would alias a native or other scalar. On error,
    !>  n_transport is 0 and bad_species is the offending input position.
    !-----------------------------------------------------------------------
    pure subroutine build_transport_map(mpas_index, host_bound, write_back, num_scalars, &
                                        chem_first, chem_last, transport_to_scalar, n_transport, &
                                        error_code, bad_species)

        integer, intent(in) :: mpas_index(:)          ! MPAS scalar index of each MICM species
        logical, intent(in) :: host_bound(:)          ! species reads an existing host field
        logical, intent(in) :: write_back(:)          ! MICM result is copied back to MPAS
        integer, intent(in) :: num_scalars
        integer, intent(in) :: chem_first, chem_last  ! runtime chemistry group bounds
        integer, intent(out) :: transport_to_scalar(:) ! size >= size(mpas_index)
        integer, intent(out) :: n_transport
        integer, intent(out) :: error_code
        integer, intent(out) :: bad_species

        integer :: i, idx
        logical :: claimed(max(num_scalars, 0))

        n_transport = 0
        error_code = MAP_OK
        bad_species = 0
        claimed(:) = .false.

        do i = 1, size(mpas_index)
            if (host_bound(i) .or. .not. write_back(i)) cycle

            idx = mpas_index(i)
            if (idx < 1 .or. idx > num_scalars) then
                error_code = MAP_OUT_OF_RANGE
            else if (chem_first < 1 .or. idx < chem_first .or. idx > chem_last) then
                error_code = MAP_NOT_CHEMISTRY
            else if (claimed(idx)) then
                error_code = MAP_DUPLICATE
            end if
            if (error_code /= MAP_OK) then
                n_transport = 0
                bad_species = i
                return
            end if

            claimed(idx) = .true.
            n_transport = n_transport + 1
            transport_to_scalar(n_transport) = idx
        end do

    end subroutine build_transport_map

    !-----------------------------------------------------------------------
    !  function transport_map_error_text
    !
    !> \brief Short description of a build_transport_map error code
    !-----------------------------------------------------------------------
    pure function transport_map_error_text(code) result(text)

        integer, intent(in) :: code
        character(len=64) :: text

        select case (code)
        case (MAP_OK)
            text = 'ok'
        case (MAP_OUT_OF_RANGE)
            text = 'scalar index out of range'
        case (MAP_NOT_CHEMISTRY)
            text = 'scalar index outside the chemistry tracer group'
        case (MAP_DUPLICATE)
            text = 'scalar index already used by another species'
        case default
            text = 'unknown transport map error'
        end select

    end function transport_map_error_text

end module mpas_chem_transport_map
