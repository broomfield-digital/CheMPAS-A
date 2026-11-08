! Model test double for mpas_atmphys_pbl_gas (P2 realized-increment proof; never part of a production build).
! check_realized_increment.sh copies it over the production module in a scratch worktree. The mode comes from
! the environment variable CHEMPAS_TEST_PBL_GAS at init:
!   source     q(n,k) += dt*S(n,k), with S(n,k) = S0*n*(1 + 0.5*sin(k)) [kg/kg/s]
!   diffusion  closed-column implicit diffusion with G = G0 at every interior interface
!   otherwise  disabled (the production stub behaviour)
!check_realized_increment.py uses the same formulas, independently, to predict the realized increments.
 module mpas_atmphys_pbl_gas
 use mpas_kind_types,only: RKIND,StrKIND
 use mpas_derived_types,only: block_type
 use mpas_log,only: mpas_log_write
 use mpas_atmphys_pbl_diffusion,only: pbl_diffuse_column

 implicit none
 private
 public:: pbl_gas_init,pbl_gas_enabled,pbl_gas_active,pbl_gas_apply,pbl_gas_finalize

 integer,parameter:: MODE_OFF = 0, MODE_SOURCE = 1, MODE_DIFFUSION = 2
 real(kind=RKIND),parameter:: S0 = 1.e-13_RKIND           ! [kg/kg/s]
 real(kind=RKIND),parameter:: G0 = 0.5_RKIND              ! [kg m-2 s-1]
 integer,save:: mode = MODE_OFF

 contains

 subroutine pbl_gas_init(block,nTransport)
 type(block_type),intent(in):: block
 integer,intent(in):: nTransport
 character(len=StrKIND):: value
 integer:: status

 value = ''
 call get_environment_variable('CHEMPAS_TEST_PBL_GAS',value,status=status)
 select case(trim(value))
    case('source')
       mode = MODE_SOURCE
    case('diffusion')
       mode = MODE_DIFFUSION
    case default
       mode = MODE_OFF
 end select
 call mpas_log_write('[transport-test] PBL gas test double: mode = '//trim(value)//', nTransport = $i', &
                     intArgs=(/nTransport/))
 end subroutine pbl_gas_init

 logical function pbl_gas_enabled()
 pbl_gas_enabled = mode /= MODE_OFF
 end function pbl_gas_enabled

 logical function pbl_gas_active(iCell)
 integer,intent(in):: iCell
 pbl_gas_active = mode /= MODE_OFF
 end function pbl_gas_active

 subroutine pbl_gas_apply(iCell,dt,m_d,q)
 integer,intent(in):: iCell
 real(kind=RKIND),intent(in):: dt
 real(kind=RKIND),intent(in),dimension(:):: m_d
 real(kind=RKIND),intent(inout),dimension(:,:):: q
 real(kind=RKIND),dimension(size(m_d)+1):: G
 integer:: k,n

 select case(mode)
    case(MODE_SOURCE)
       do k = 1, size(q,2)
       do n = 1, size(q,1)
          q(n,k) = q(n,k) + dt*S0*n*(1._RKIND + 0.5_RKIND*sin(real(k,RKIND)))
       enddo
       enddo
    case(MODE_DIFFUSION)
       G(:) = G0
       call pbl_diffuse_column(dt,m_d,G,q)
 end select
 end subroutine pbl_gas_apply

 subroutine pbl_gas_finalize
 end subroutine pbl_gas_finalize

 end module mpas_atmphys_pbl_gas
