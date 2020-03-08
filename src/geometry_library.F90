

module geometry_library
#ifdef MPI
  use mpi
#endif
use modcom
implicit none
private
type, public :: geomlib
  integer nspec
  integer nmaxatm_pspec
  integer, allocatable ::  nat_per_spec(:)
  real(dp) avec(NDIM,NDIM)
  real(dp) bvec(NDIM,NDIM)
  real(dp), allocatable :: atml(:,:,:)
  contains 
  procedure :: init=>generate_structure
endtype

contains

subroutine generate_structure(THIS,istruct,mode)
class(geomlib), intent(inout) :: THIS
character(len=*), intent(in) :: mode
integer, intent(in) :: istruct
if (trim(adjustl(mode)).eq."tbg") then
   call generate_structure_tbg3D(THIS,istruct)
else if (trim(adjustl(mode)).eq."slg") then
   call generate_structure_slg3D(THIS)
else
   call throw("geometry_library%generate_structure()","unknown mode parameter")
end if
end subroutine



subroutine generate_structure_slg3D(THIS)
class(geomlib), intent(inout) :: THIS
integer, parameter :: number_of_atoms=2
real(dp) tvec(3,3)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
THIS%avec(1,:)=graphene_lvec_length*(/ 0.5_dp,0.5_dp*sqrt(3._dp),0._dp/)
THIS%avec(2,:)=graphene_lvec_length*(/-0.5_dp,0.5_dp*sqrt(3._dp),0._dp/)
THIS%avec(3,:)=     tbg_ab_distance*(/  0._dp,            0._dp,1._dp/)
! reciprocal lattice vectors
tvec=THIS%avec
call dmatrix_inverse(tvec,THIS%bvec,NDIM)
THIS%bvec=transpose(THIS%bvec)*twopi
THIS%nspec=1
THIS%nmaxatm_pspec=2
allocate(THIS%nat_per_spec(THIS%nspec))
allocate(THIS%atml(3,THIS%nmaxatm_pspec,THIS%nspec))
THIS%nat_per_spec(:)=THIS%nmaxatm_pspec
THIS%atml(:,1,1)=(/ 0._dp  , 0._dp ,0._dp/)
THIS%atml(:,2,1)=(/ twothrd , twothrd,0._dp/)
if (mp_mpi) then
  call info("generate_structure_slg()","")
  write(*,*) "avec:"
  write(*,'(10F16.6)')THIS%avec(1,:)
  write(*,'(10F16.6)')THIS%avec(2,:)
  write(*,'(10F16.6)')THIS%avec(3,:)
  write(*,*) "bvec:"
  write(*,'(10F16.6)')THIS%bvec(1,:)
  write(*,'(10F16.6)')THIS%bvec(2,:)
  write(*,'(10F16.6)')THIS%bvec(3,:)
  write(*,'("Number of atoms: ",I6)') THIS%nspec*THIS%nmaxatm_pspec
  write(*,*)
end if
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
return
end subroutine



subroutine generate_structure_tbg3D(THIS,jstruct)
class(geomlib), intent(inout) :: THIS
integer, intent(in) :: jstruct
integer, parameter :: ntrans=140
integer i,j,iat,istruct
integer counter1,counter2,counter
integer number_of_atoms,spec0(4)
real(dp) t1,t2,t3,z,zmid,vac
real(dp) costheta,theta,thetadeg
real(dp) a33,tbg_interplane_dist,zlat,dvar
real(dp) atl(3),atc(3),atl_super(3)
real(dp) rot(3,3),atxy(4,3)
real(dp) ave(3,3),ave2(3,3),tvec(3,3),bvect(3,3)
integer, allocatable :: spec(:)
real(dp), allocatable :: atmlt(:,:)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (jstruct.lt.0) call throw("geomlib%generate_structure_tbg()","integer structure value of TBG should be positive or zero")
if (jstruct.ge.100) then
  istruct=jstruct-100
else
  istruct=jstruct
end if
if (istruct.gt.99) call throw("geomlib%generate_structure_tbg()","integer structure value of TBG is too large")
if (NDIM.ne.3) call throw("geomlib%generate_structure_tbg()","this code works only for three dimensions")
ave(1,:)=graphene_lvec_length*(/ 0.5_dp,0.5_dp*sqrt(3._dp), 0._dp/)
ave(2,:)=graphene_lvec_length*(/-0.5_dp,0.5_dp*sqrt(3._dp), 0._dp/)
ave(3,:)=   (/  0._dp,            0._dp, 1._dp/)
costheta=dble(3*istruct**2+3*istruct+0.5)/dble(3*istruct**2+3*istruct+1)
theta=acos(costheta)
thetadeg=theta*57.29577951307854999853_dp
! superlattice
THIS%avec=ave
THIS%avec(1,1)=istruct*ave(1,1)+(istruct+1)*ave(2,1)
THIS%avec(1,2)=istruct*ave(1,2)+(istruct+1)*ave(2,2)
THIS%avec(2,1)=-(istruct+1)*ave(1,1)+(2*istruct+1)*ave(2,1)
THIS%avec(2,2)=-(istruct+1)*ave(1,2)+(2*istruct+1)*ave(2,2)
! renormalize z
t1=sqrt(THIS%avec(1,1)**2+THIS%avec(1,2)**2)
t2=sqrt(THIS%avec(2,1)**2+THIS%avec(2,2)**2)
vac=max(t1,t2)
tbg_interplane_dist=tbg_ab_distance
a33=tbg_interplane_dist+vac
zlat=tbg_interplane_dist/a33
THIS%avec(3,3)=a33
! rescale z of ave, generate ave2 rotated with respect to ave
ave(3,3)=a33
call getrot(.true.,theta,rot)
ave2=matmul(ave,rot)
! reciprocal lattice vectors
tvec=THIS%avec
call dmatrix_inverse(tvec,THIS%bvec,NDIM)
THIS%bvec=transpose(THIS%bvec)*twopi
! additional vectors for supercell construction
bvect=transpose(THIS%bvec)/twopi
! compute number of atoms
t1=abs(THIS%avec(1,1)*THIS%avec(2,2)-THIS%avec(1,2)*THIS%avec(2,1))
t2=abs(ave(1,1)*ave(2,2)-ave(1,2)*ave(2,1))
number_of_atoms=nint(4._dp*t1/t2)
if (mp_mpi) then
  call info("generate_structure_tbg()","")
  write(*,'("istruct: ",I6)') istruct
  write(*,'("costheta: ",F10.6)') costheta
  write(*,'("theta (radians): ",F10.6)') theta
  write(*,'("theta (degrees): ",F10.6)') thetadeg
  write(*,*) "avec:"
  write(*,'(10F16.6)')THIS%avec(1,:)
  write(*,'(10F16.6)')THIS%avec(2,:)
  write(*,'(10F16.6)')THIS%avec(3,:)
  write(*,*) "bvec:"
  write(*,'(10F16.6)')THIS%bvec(1,:)
  write(*,'(10F16.6)')THIS%bvec(2,:)
  write(*,'(10F16.6)')THIS%bvec(3,:)
  write(*,'("Number of atoms: ",I6)') number_of_atoms
  write(*,*)
end if
! AB stacking
atxy(1,:)=(/ 0._dp  , 0._dp ,0._dp/)
atxy(2,:)=(/ twothrd , twothrd,0._dp/)
atxy(3,:)=(/-twothrd ,-twothrd,zlat/)
atxy(4,:)=(/ 0._dp  , 0._dp ,zlat/)
! AA stacking
atxy(1,:)=(/ 0._dp  , 0._dp ,0._dp/)  
atxy(2,:)=(/ twothrd , twothrd,0._dp/)
atxy(3,:)=(/ 0._dp  , 0._dp ,zlat/)
atxy(4,:)=(/ twothrd , twothrd,zlat/)
spec0(1)=1
spec0(2)=2
spec0(3)=1
spec0(4)=2
counter=number_of_atoms*(2*ntrans+1)*(2*ntrans+1)
allocate(atmlt(3,counter))
allocate(spec(counter))
counter=0
do i=-ntrans-1,ntrans
  do j=-ntrans-1,ntrans
    do iat=1,4
      atl=(/atxy(iat,1)+dble(i),atxy(iat,2)+dble(j),atxy(iat,3)/)
      if (abs(atxy(iat,3)).lt.epslat)then
        ! first layer
        atc=matmul(atl,ave)
      else
        ! second layer
        atc=matmul(atl,ave2)
      end if
      ! lattice corrdinates with respect to Moire lattice vecors
      atl_super=matmul(atc,bvect)
      ! main array
      if (( (atl_super(1).ge.-epslat    ).and.(atl_super(2).ge.-epslat    ) ) .and.&
          ( (atl_super(1).lt.1._dp-epslat).and.(atl_super(2).lt.1._dp-epslat) ) ) then
        ! lattice corrdinates with respect to Moire lattice vecors
        counter=counter+1
        atmlt(:,counter)=atl_super(:)
        spec(counter)=spec0(iat)
      end if
    end do
  end do
end do
if (number_of_atoms.ne.counter) then
  call throw("geomlib%generate_structure_tbg()","number of atoms counted is not correct")
end if
! allocate geometry library arrays
THIS%nspec=2
THIS%nmaxatm_pspec=number_of_atoms/2
allocate(THIS%nat_per_spec(THIS%nspec))
allocate(THIS%atml(3,THIS%nmaxatm_pspec,THIS%nspec))
THIS%nat_per_spec(:)=THIS%nmaxatm_pspec
zmid=0.5_dp*tbg_ab_distance
counter1=0
counter2=0
do iat=1,number_of_atoms
  if (spec(iat).eq.1) then
    counter1=counter1+1
    THIS%atml(:,counter1,spec(iat))=atmlt(:,iat)
  else if (spec(iat).eq.2) then
    counter2=counter2+1
    THIS%atml(:,counter2,spec(iat))=atmlt(:,iat)
  else
    call throw("geomlib%generate_structure_tbg()","unexpected species index")
  end if
  atc(:)=matmul(atmlt(:,iat),THIS%avec)
  if (jstruct.gt.100) then
    ! add distortion in Z
    t1=dot_product(THIS%bvec(1,:),atc)
    t2=dot_product(THIS%bvec(2,:),atc)
    t3=dot_product(THIS%bvec(1,:)+THIS%bvec(2,:),atc)
    tbg_interplane_dist=(tbg_aa_distance+2._dp*tbg_ab_distance)/3._dp
    dvar=(tbg_aa_distance-tbg_ab_distance)/9._dp
    z=tbg_interplane_dist+2._dp*dvar*(cos(t1)+cos(t2)+cos(t3))
  else
    z=tbg_interplane_dist
  end if
  if (atc(3).lt.zmid) then
    t1=-0.5_dp*z
  else
    t1=0.5_dp*z
  end if
  t2=t1/a33
  if (spec(iat).eq.1) then
    THIS%atml(3,counter1,spec(iat))=t2
  else if (spec(iat).eq.2) then
    THIS%atml(3,counter2,spec(iat))=t2
  end if
end do
deallocate(atmlt,spec)
#ifdef MPI
  call MPI_barrier(mpi_com,mpi_err)
#endif
return
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! side subroutines 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getrot(trans,theta,rotmat)
implicit none
logical, intent(in)  :: trans
real(dp), intent(in)  :: theta
real(dp), intent(out) :: rotmat(3,3)
if (trans) then
  rotmat(1,:)=(/ cos(theta), sin(theta),0._dp/)
  rotmat(2,:)=(/-sin(theta), cos(theta),0._dp/)
else
  rotmat(1,:)=(/ cos(theta),-sin(theta),0._dp/)
  rotmat(2,:)=(/ sin(theta), cos(theta),0._dp/)
end if
rotmat(3,:)=(/0._dp,0._dp,1._dp/)
return
end subroutine


end module
