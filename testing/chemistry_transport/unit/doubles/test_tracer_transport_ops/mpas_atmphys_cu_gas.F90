! Test double for the convective gas process: active in cells 1 and 2. Species 1
! moves all of layer 2 into layer 3 (conservative); species 2 doubles its top
! layer (a deliberate source, so the budget increment is nonzero).
 module mpas_atmphys_cu_gas
 use mpas_kind_types,only: RKIND
 use mpas_derived_types,only: block_type

 implicit none
 private
 public:: cu_gas_init,cu_gas_enabled,cu_gas_active,cu_gas_apply,cu_gas_finalize,cu_double_active
 logical,save:: cu_double_active = .true.   ! .false.: no column is active (an all-inactive rank)

 contains

 subroutine cu_gas_init(block,nTransport)
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport
 end subroutine cu_gas_init

 logical function cu_gas_enabled()
 cu_gas_enabled = .true.
 end function cu_gas_enabled

 logical function cu_gas_active(iCell)
 integer,intent(in):: iCell
 cu_gas_active = cu_double_active .and. iCell <= 2
 end function cu_gas_active

 subroutine cu_gas_apply(iCell,dt,m_d,q)
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in):: dt,m_d(:)
 real(kind=RKIND),intent(inout):: q(:,:)
 q(1,3) = q(1,3) + q(1,2)*m_d(2)/m_d(3)
 q(1,2) = 0._RKIND
 q(2,3) = 2._RKIND*q(2,3)
 end subroutine cu_gas_apply

 subroutine cu_gas_finalize
 end subroutine cu_gas_finalize

 end module mpas_atmphys_cu_gas
