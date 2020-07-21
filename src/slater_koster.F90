

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
    procedure, nopass :: init=>init_sk_pars
    procedure, nopass :: tij
endtype 


contains


subroutine init_sk_pars(option)
character(len=*), intent(in) :: option
if (trim(adjustl(option)).eq."tbgsk") then
  ! pure SK model with parameters from the literature (like Koshino)
  qpz_pi=3.14_dp
  qpz_sig=7.43_dp
  tpz_pi0=-2.7_dp
  tpz_sig0=0.48_dp
else if (trim(adjustl(option)).eq."tbgsk1") then
  ! SK parameters obtained by fitting Siesta's single zeta result
  qpz_pi=2.5603_dp
  qpz_sig=3.2911_dp
  tpz_pi0=-35.6714_dp
  tpz_sig0=0.3073_dp
else if (trim(adjustl(option)).eq."tbgsk2") then
  ! SK parameters obtained by fitting Siesta's double zeta-p result
  qpz_pi=2.6123_dp
  qpz_sig=2.8313_dp
  tpz_pi0=-57.7144_dp
  tpz_sig0=0.3372_dp
else if (trim(adjustl(option)).eq."tbgsk3") then
  ! SK parameters obtained by fitting Siesta's signle-zeta result, but with PBE XC
  qpz_pi=2.6218_dp
  qpz_sig=3.4203_dp
  tpz_pi0=-36.7532_dp
  tpz_sig0=0.2897_dp
else if (trim(adjustl(option)).eq."tbgsk4") then
  ! SK parameters obtained by fitting Siesta's double zeta-p result, but with PBE XC
  qpz_pi=2.6690_dp
  qpz_sig=3.1294_dp
  tpz_pi0=-57.0311_dp
  tpz_sig0=0.3317_dp
else if (trim(adjustl(option)).eq."tbgsk_old") then
  ! old out of-plane hopping parameterisation
  qpz_pi=2.484913272_dp
  qpz_sig=3.031741524_dp
  tpz_pi0=-32.872337923_dp
  tpz_sig0=0.306297655_dp
else if (trim(adjustl(option)).eq."tbgsk1piold") then
  ! new tbgsk1, but with tpi from tbgsk_old
  qpz_pi=2.5603_dp
  qpz_sig=3.2911_dp
  tpz_pi0=-32.0000_dp
  tpz_sig0=0.3073_dp
else if (trim(adjustl(option)).eq."tbgsk1sigmaold") then
  qpz_pi=2.5603_dp
  qpz_sig=3.0300_dp
  tpz_pi0=-35.6714_dp
  tpz_sig0=0.3073_dp
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
  select case(trim(adjustl(option))) 
  case('tbgsk')
    ! full SK from the literature
    tij=(  tpz_pi(rr)*(1._dp-zz**2)+tpz_sig(rr)*zz**2  )!*fcut(rr)
  case('tbgsk1','tbgsk2','tbgsk3','tbgsk4','tbgsk_old','tbgsk1piold','tbgsk1sigmaold')
    ! full mixed SK (out-of-plane) with ab-initio (in-plane)
    if (abs(dvec(ZAXIS)).gt.0.5_dp*tbg_ab_distance) then
       ! out-of-plane
       tij=(  tpz_pi(rr)*(1._dp-zz**2)+tpz_sig(rr)*zz**2  )!*fcut(rr)
    else
       ! in-plane
       tij=tbg_inplane_table(option,rr)
    end if
  case default
    call throw("slater_koster%tij()","unknown input otion")
  end select
else
  tij=0._dp
end if
end function

real(dp) function tbg_inplane_table(option,rr)
character(len=*), intent(in) :: option
real(dp), intent(in) :: rr
real(dp) units_of_lvec
units_of_lvec=rr/graphene_lvec_length
select case(trim(adjustl(option))) 
case('tbgsk1','tbgsk_old','tbgsk1piold','tbgsk1sigmaold')
  ! single-zeta bais
  if (units_of_lvec.lt.epslat) then
    tbg_inplane_table=0._dp ! 0nn
  else if (units_of_lvec.lt.0.90_dp) then
    tbg_inplane_table=-2.7783_dp ! 1nn
  else if (units_of_lvec.lt.1.15_dp) then
    tbg_inplane_table= 0.2292_dp ! 2nn
  else if (units_of_lvec.lt.1.52_dp) then
    tbg_inplane_table=-0.1719_dp ! 3nn
  else if (units_of_lvec.lt.1.73_dp) then
    tbg_inplane_table= 0.0075_dp ! 4nn
  else if (units_of_lvec.lt.1.99_dp) then
    tbg_inplane_table= 0.0218_dp ! 5nn
  else if (units_of_lvec.lt.2.07_dp) then
    tbg_inplane_table=-0.0078_dp ! 6nn
  else if (units_of_lvec.lt.2.29_dp) then
    tbg_inplane_table=-0.0035_dp ! 7nn
  else if (units_of_lvec.lt.2.50_dp) then
    tbg_inplane_table=-0.0042_dp ! 8nn
  else if (units_of_lvec.lt.2.63_dp) then
    tbg_inplane_table= 0.0018_dp ! 9nn
  else 
    tbg_inplane_table= 0._dp
  end if
case ('tbgsk2')
  ! double-zeta bais
  if (units_of_lvec.lt.epslat) then
    tbg_inplane_table=0._dp ! 0nn
  else if (units_of_lvec.lt.0.90_dp) then
    tbg_inplane_table=-2.8464_dp ! 1nn
  else if (units_of_lvec.lt.1.15_dp) then
    tbg_inplane_table= 0.2176_dp ! 2nn
  else if (units_of_lvec.lt.1.52_dp) then
    tbg_inplane_table=-0.2375_dp ! 3nn
  else if (units_of_lvec.lt.1.73_dp) then
    tbg_inplane_table= 0.0184_dp ! 4nn
  else if (units_of_lvec.lt.1.99_dp) then
    tbg_inplane_table= 0.0400_dp ! 5nn
  else if (units_of_lvec.lt.2.07_dp) then
    tbg_inplane_table=-0.0141_dp ! 6nn
  else if (units_of_lvec.lt.2.29_dp) then
    tbg_inplane_table=-0.0179_dp ! 7nn
  else if (units_of_lvec.lt.2.50_dp) then
    tbg_inplane_table=-0.0031_dp ! 8nn
  else if (units_of_lvec.lt.2.63_dp) then
    tbg_inplane_table= 0.0070_dp ! 9nn
  else if (units_of_lvec.lt.2.87_dp) then
    tbg_inplane_table= 0.0018_dp ! 10nn
  else if (units_of_lvec.lt.2.98_dp) then
    tbg_inplane_table=-0.0046_dp ! 11nn
  else if (units_of_lvec.lt.3.045_dp) then
    tbg_inplane_table=-0.0018_dp ! 12nn
  else 
    tbg_inplane_table= 0._dp
  end if
case('tbgsk3')
  ! single-zeta PBE
  if (units_of_lvec.lt.epslat) then
    tbg_inplane_table=0._dp ! 0nn
  else if (units_of_lvec.lt.0.90_dp) then
    tbg_inplane_table=-2.7903_dp ! 1nn
  else if (units_of_lvec.lt.1.15_dp) then
    tbg_inplane_table= 0.2282_dp ! 2nn
  else if (units_of_lvec.lt.1.52_dp) then
    tbg_inplane_table=-0.1662_dp ! 3nn
  else if (units_of_lvec.lt.1.73_dp) then
    tbg_inplane_table= 0.0073_dp ! 4nn
  else if (units_of_lvec.lt.1.99_dp) then
    tbg_inplane_table= 0.0207_dp ! 5nn
  else if (units_of_lvec.lt.2.07_dp) then
    tbg_inplane_table=-0.0075_dp ! 6nn
  else if (units_of_lvec.lt.2.29_dp) then
    tbg_inplane_table=-0.0033_dp ! 7nn
  else if (units_of_lvec.lt.2.50_dp) then
    tbg_inplane_table=-0.0040_dp ! 8nn
  else if (units_of_lvec.lt.2.63_dp) then
    tbg_inplane_table= 0.0017_dp ! 9nn
  else 
    tbg_inplane_table= 0._dp
  end if
case('tbgsk4')
  ! double-zeta PBE
  if (units_of_lvec.lt.epslat) then
    tbg_inplane_table=0._dp ! 0nn
  else if (units_of_lvec.lt.0.90_dp) then
    tbg_inplane_table=-2.8400_dp ! 1nn
  else if (units_of_lvec.lt.1.15_dp) then
    tbg_inplane_table= 0.2188_dp ! 2nn
  else if (units_of_lvec.lt.1.52_dp) then
    tbg_inplane_table=-0.2293_dp ! 3nn
  else if (units_of_lvec.lt.1.73_dp) then
    tbg_inplane_table= 0.0146_dp ! 4nn
  else if (units_of_lvec.lt.1.99_dp) then
    tbg_inplane_table= 0.0390_dp ! 5nn
  else if (units_of_lvec.lt.2.07_dp) then
    tbg_inplane_table=-0.0126_dp ! 6nn
  else if (units_of_lvec.lt.2.29_dp) then
    tbg_inplane_table=-0.0173_dp ! 7nn
  else if (units_of_lvec.lt.2.50_dp) then
    tbg_inplane_table=-0.0032_dp ! 8nn
  else if (units_of_lvec.lt.2.63_dp) then
    tbg_inplane_table= 0.0066_dp ! 9nn
  else if (units_of_lvec.lt.2.87_dp) then
    tbg_inplane_table= 0.0017_dp ! 10nn
  else if (units_of_lvec.lt.2.98_dp) then
    tbg_inplane_table=-0.0042_dp ! 11nn
  else if (units_of_lvec.lt.3.045_dp) then
    tbg_inplane_table=-0.0019_dp ! 12nn
  else 
    tbg_inplane_table= 0._dp
  end if
case default
  call throw("slater_koster%tij()","unknown input otion")
end select
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
