! Copyright (c) 2026 The University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at https://mpas-dev.github.io/license.html .
!
!=================================================================================================================
 module mpas_atmphys_cu_gas
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: block_type

 implicit none
 private
 public:: cu_gas_init,      &
          cu_gas_enabled,   &
          cu_gas_active,    &
          cu_gas_apply,     &
          cu_gas_finalize


!Convective transport of MICM gases: shared stub (P1 of docs/plans/pbl-transport.md). Never enabled; the
!scheme-specific body replaces this module on chempas-a-convect. Interface used by
!mpas_atmphys_tracer_transport:
!   cu_gas_init(block,nTransport)   resolve the scheme and check capability, once after chemistry_init
!   cu_gas_enabled()                true when convective gas transport can act in any column
!   cu_gas_active(iCell)            true when the held operator acts in column iCell this step
!   cu_gas_apply(iCell,dt,m_d,q)    q(nTransport,nVertLevels) -> q after dt [kg/kg dry air], closed column
!   cu_gas_finalize()


 contains


!=================================================================================================================
 subroutine cu_gas_init(block,nTransport)
!=================================================================================================================
!block gives configs, structs and dimensions; nTransport is the number of independently transported gases.
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport

 end subroutine cu_gas_init

!=================================================================================================================
 logical function cu_gas_enabled()
!=================================================================================================================

 cu_gas_enabled = .false.

 end function cu_gas_enabled

!=================================================================================================================
 logical function cu_gas_active(iCell)
!=================================================================================================================
 integer,intent(in):: iCell

 cu_gas_active = .false.

 end function cu_gas_active

!=================================================================================================================
 subroutine cu_gas_apply(iCell,dt,m_d,q)
!=================================================================================================================
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(:):: m_d
 real(kind=RKIND),intent(inout),dimension(:,:):: q

!identity.

 end subroutine cu_gas_apply

!=================================================================================================================
 subroutine cu_gas_finalize
!=================================================================================================================

 end subroutine cu_gas_finalize

!=================================================================================================================
 end module mpas_atmphys_cu_gas
!=================================================================================================================
