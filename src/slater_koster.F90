

module slater_koster
use modcom
implicit none
private
real(dp) :: qpz_pi
real(dp) :: qpz_sig
real(dp) :: tpz_pi0
real(dp) :: tpz_sig0

type, public :: SK
    contains
    procedure :: init=>init_sk_pars
    procedure, nopass :: tij
endtype 


contains


subroutine init_sk_pars(THIS,option)
class(SK), intent(out) :: THIS
character(len=*), intent(in) :: option
if (trim(adjustl(option)).eq."tbgsk1") then
  ! in this case in-planle hoppins should be read from the table
  qpz_pi=2.484913272_dp
  qpz_sig=3.031741524_dp
  tpz_pi0=-32.872337923_dp
  tpz_sig0=0.306297655_dp
else if (trim(adjustl(option)).eq."tbgsk") then
  ! pure SK model with parameters from the literature (like Koshino)
  qpz_pi=3.14_dp
  qpz_sig=7.43_dp
  tpz_pi0=-2.7_dp
  tpz_sig0=0.48_dp
else
  call throw("SK%init_sk_pars()","unknown option")
end if
end subroutine

!real(dp) function tij_sk(THIS,dvec)

real(dp) function tij(option,lmr1,lmr2,dvec)
!class(SK), intent(in) :: THIS
character(len=*), intent(in) :: option
integer, intent(in)  :: lmr1(2),lmr2(2)
real(dp), intent(in) :: dvec(NDIM)
real(dp) rr,zz
if (lmr1(1).ne.1.or.lmr2(2).ne.1) call throw("SK%tij_sk","this subroutine is currently for pz-pz hoppings only")
if (lmr2(1).ne.1.or.lmr2(2).ne.1) call throw("SK%tij_sk","this subroutine is currently for pz-pz hoppings only")
if (NDIM.ne.3) call throw("SK%tij_sk","this subroutine is for 3D case only")
rr=sqrt(sum(dvec(:)**2))
if (abs(rr).gt.epslat) then
  zz=dvec(ZAXIS)/rr
  if (trim(adjustl(option)).eq.'tbgsk') then
    ! full SK from the literature
    tij=(  tpz_pi(rr)*(1._dp-zz**2)+tpz_sig(rr)*zz**2  )*fcut(rr)
  else if (trim(adjustl(option)).eq.'tbgsk1') then
    ! full mixed SK (out-of-plane) with ab-initio (in-plane)
    if (abs(dvec(ZAXIS)).gt.0.5_dp*tbg_ab_distance) then
       ! out-of-plane
       tij=(  tpz_pi(rr)*(1._dp-zz**2)+tpz_sig(rr)*zz**2  )*fcut(rr)
    else
       ! in-plane
       tij=tbg_inplane_table(rr)
    end if
  else
    call throw("slater_koster%tij()","unknown input otion")
  end if
  !if (rr.lt.7.and.abs(dvec(ZAXIS)).gt.0.5_dp*tbg_ab_distance) then
  !  write(*,'(10F10.4)') dvec,rr,tij
  !end if
else
  tij=0._dp
end if
end function

real(dp) function tbg_inplane_table(rr)
real(dp), intent(in) :: rr
real(dp) units_of_lvec
units_of_lvec=rr/graphene_lvec_length
!    hop_inplt(1)=-2.88421920d0
!    hop_inplt(2)= 0.21383630d0
!    hop_inplt(3)=-0.19803719d0
!    hop_inplt(4)= 0.01992761d0
!    hop_inplt(5)= 0.02505946d0
!    hop_inplt(6)=-0.01083505d0
!    hop_inplt(7)=-0.00384209d0
!    hop_inplt(8)=-0.00513276d0
!    hop_inplt(9)= 0.00226530d0
!    hop_inplt(10)=0.00030354d0
if (units_of_lvec.lt.epslat) then
  tbg_inplane_table=0._dp
else if (units_of_lvec.lt.0.90_dp) then
  ! first nearest neighbor (at graphene cc distanse ~ 0.5)
  tbg_inplane_table=-2.88421920_dp
else if (units_of_lvec.lt.1.15_dp) then
  ! secpnd nn (at cc dist ~1.0)
  tbg_inplane_table= 0.21383630_dp
else if (units_of_lvec.lt.1.52_dp) then
  ! and so on
  tbg_inplane_table=-0.19803719_dp
else if (units_of_lvec.lt.1.73_dp) then
  tbg_inplane_table= 0.01992761_dp
else if (units_of_lvec.lt.1.99_dp) then
  tbg_inplane_table= 0.02505946_dp
else if (units_of_lvec.lt.2.07_dp) then
  tbg_inplane_table=-0.01083505_dp
else 
  tbg_inplane_table= 0._dp
end if
  
end function

real(dp) function tpz_pi(dd)
real(dp), intent(in) :: dd
if (abs(dd).lt.0.5_dp*graphene_cc_distance) then
  tpz_pi=0._dp
else
  tpz_pi=tpz_pi0*exp(qpz_pi*(1._dp-dd/graphene_cc_distance))
end if
end function

real(dp) function tpz_sig(dd)
real(dp), intent(in) :: dd
if (abs(dd).lt.0.5_dp*tbg_ab_distance) then
   tpz_sig=0._dp
else
   tpz_sig=tpz_sig0*exp(qpz_sig*(1._dp-dd/tbg_ab_distance))
end if
end function

real(dp) function fcut(dd)
  real(dp), intent(in) :: dd
  real(dp), parameter :: rc=6.14_dp
  real(dp), parameter :: lc=0.265_dp
  real(dp) t1,t2
  t1=dd-rc
  t2=lc
  fcut=1._dp/(1._dp+exp(t1/t2))
end function fcut

end module
