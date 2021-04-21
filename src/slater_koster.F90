

module slater_koster
use modcom
implicit none
private

! SK parameters correspond to the line [q_pi,q_sig,t_pi0,t_sig0], currently for TBG only
!
! the SK model is 
! t(r) = t_pi(r)*(1-n^2) + t_sig(r)*n^2
! where n = z/r  (r=[x,y,z])
! t_pi(r) = t_pi0 *exp[q_pi*(1-r/r_cc)]
! t_sig(r) = t_sig0 *exp[q_sig*(1-r/d_ab)]
! r_cc, d_ab carbon-carbon and AB graphen interplain distance

real(dp) :: skpar_ipl(4)
real(dp) :: skpar_opl(4)
real(dp) :: sktab_ipl(15)
real(dp) :: sktab_opl(15)

type, public :: SK
    character(len=100) :: sktype
    contains
    procedure :: init=>init_sk_pars
    procedure, nopass :: tij
endtype 


contains


subroutine init_sk_pars(self, option, inpl, oupl)
class(SK), intent(inout) :: self
character(len=*), intent(in) :: option, inpl, oupl
skpar_ipl=0._dp
skpar_opl=0._dp
sktab_ipl=0._dp
sktab_opl=0._dp
self%sktype = option
select case(trim(adjustl(option))) 
case('tbgsk')

  ! init in-plane parameters for TBG
  select case(trim(adjustl(inpl))) 
  case('sz_table')
    sktab_ipl(1:9)  = (/-2.7783_dp, 0.2292_dp,-0.1719_dp,0.0075_dp,0.0218_dp,&
                        -0.0078_dp,-0.0035_dp,-0.0042_dp,0.0018_dp/)
  case('dzp_table')
    sktab_ipl(1:12) = (/-2.8464_dp, 0.2176_dp,-0.2375_dp,0.0184_dp,0.0400_dp,&
                        -0.0141_dp,-0.0179_dp,-0.0031_dp,0.0070_dp,0.0018_dp,-0.0046_dp,-0.0018_dp/)
  case('original')
    ! pure SK model with parameters from the literature (like Koshino)
    skpar_ipl = (/3.14_dp,7.43_dp,-2.7_dp,.48_dp/)
  case default
    call throw("slater_koster%init_sk_pars()","unknown inpl parameter")
  end select

  ! init out-of-plane parameters for TBG
  select case(trim(adjustl(oupl))) 
  case('original')
    ! pure SK model with parameters from the literature (like Koshino)
    skpar_opl = (/3.14_dp,7.43_dp,-2.7_dp,.48_dp/)
  case('sz_fit')
    skpar_opl = (/2.5603_dp,3.2911_dp,-35.6714_dp,0.3073_dp/)
  case('dzp_fit')
    skpar_opl = (/2.6123_dp,2.8313_dp,57.7144_dp,0.3372_dp/)
  case default
    call throw("slater_koster%init_sk_pars()","unknown inpl parameter")
  end select

case default

  call throw("slater_koster%init_sk_pars()","unknown option parameter")

end select
end subroutine

!real(dp) function tij_sk(THIS,dvec)

real(dp) function tij(option,lmr1,lmr2,dvec,inpl,oupl)
!class(SK), intent(in) :: THIS
character(len=*), intent(in) :: option
character(len=*), intent(in) :: inpl,oupl
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

    ! full mixed SK (out-of-plane) with ab-initio (in-plane)
    if (abs(dvec(ZAXIS)).gt.0.5_dp*tbg_ab_distance) then
       ! out-of-planep
       select case(trim(adjustl(oupl))) 
       case('original','sz_fit','dzp_fit')
         tij = t_tbg(rr,zz,skpar_opl)
       case default
         call throw("slater_koster%tij()","unknown oupl value")
       end select
    else
       ! in-plane
       select case(trim(adjustl(inpl))) 
       case('original','sz_fit')
         tij = t_tbg(rr,zz,skpar_ipl)
       case('sz_table','dzp_table')
         tij = tbg_inplane_table(rr,zz)
       case default
         call throw("slater_koster%tij()","unknown inpl value")
       end select
    end if

  case default
    call throw("slater_koster%tij()","unknown input otion")
  end select

else
  tij=0._dp
end if
end function

real(dp) function tbg_inplane_table(rr,zz)
real(dp), intent(in) :: rr,zz
real(dp) units_of_lvec
units_of_lvec=rr/graphene_lvec_length
if (units_of_lvec.lt.epslat)   then ; tbg_inplane_table = 0._dp         ; return ; endif ! 0nn
if (units_of_lvec.lt.0.90_dp)  then ; tbg_inplane_table = sktab_ipl(1)  ; return ; endif ! 1nn
if (units_of_lvec.lt.1.15_dp)  then ; tbg_inplane_table = sktab_ipl(2)  ; return ; endif ! 2nn
if (units_of_lvec.lt.1.52_dp)  then ; tbg_inplane_table = sktab_ipl(3)  ; return ; endif ! 3nn
if (units_of_lvec.lt.1.73_dp)  then ; tbg_inplane_table = sktab_ipl(4)  ; return ; endif ! 4nn
if (units_of_lvec.lt.1.99_dp)  then ; tbg_inplane_table = sktab_ipl(5)  ; return ; endif ! 5nn
if (units_of_lvec.lt.2.07_dp)  then ; tbg_inplane_table = sktab_ipl(6)  ; return ; endif ! 6nn
if (units_of_lvec.lt.2.29_dp)  then ; tbg_inplane_table = sktab_ipl(7)  ; return ; endif ! 7nn
if (units_of_lvec.lt.2.50_dp)  then ; tbg_inplane_table = sktab_ipl(8)  ; return ; endif ! 8nn
if (units_of_lvec.lt.2.63_dp)  then ; tbg_inplane_table = sktab_ipl(9)  ; return ; endif ! 9nn
if (units_of_lvec.lt.2.87_dp)  then ; tbg_inplane_table = sktab_ipl(10) ; return ; endif ! 10nn
if (units_of_lvec.lt.2.98_dp)  then ; tbg_inplane_table = sktab_ipl(11) ; return ; endif ! 11nn
if (units_of_lvec.lt.3.045_dp) then ; tbg_inplane_table = sktab_ipl(12) ; return ; endif ! 12nn
                                      tbg_inplane_table = 0._dp
end function

real(dp) function tpz_pi(tt,qq,dd)
real(dp), intent(in) :: tt,qq,dd
if (abs(dd).lt.0.5_dp*graphene_cc_distance) then
  tpz_pi=0._dp
else
  tpz_pi=tt*exp(qq*(1._dp-dd/graphene_cc_distance))
end if
end function

real(dp) function tpz_sig(tt,qq,dd)
real(dp), intent(in) :: tt,qq,dd
if (abs(dd).lt.0.5_dp*tbg_ab_distance) then
   tpz_sig=0._dp
else
   tpz_sig=tt*exp(qq*(1._dp-dd/tbg_ab_distance))
end if
end function

real(dp) function t_tbg(rr,zz,pars)
real(dp), intent(in) :: rr,zz
real(dp), intent(in) :: pars(*)
real(dp) qpz_pi,qpz_sig,tpz_pi0,tpz_sig0
qpz_pi=pars(1)
qpz_sig=pars(2)
tpz_pi0=pars(3)
tpz_sig0=pars(4)
t_tbg = (  tpz_pi(tpz_pi0,qpz_pi,rr) *(1._dp-zz**2) + tpz_sig(tpz_sig0,qpz_sig,rr) *zz**2  )!*fcut(rr)
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
