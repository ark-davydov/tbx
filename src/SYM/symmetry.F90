
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmetry(symtype_,nspec_,nat_per_spec_,nmaxatm_pspec_,natmtot_,avec_,atml_,&
    nsymcrys_,lsplsymc_,lspnsymc_,invmap_,symlat_,ieqatom_,ieqatom_inv_,vtlsymc_)
use modelk
implicit none
integer, intent(in)    :: symtype_,nspec_,nat_per_spec_(nspec_),nmaxatm_pspec_,natmtot_
real(8), intent(in)    :: avec_(3,3)
real(8), intent(inout) :: atml_(3,nmaxatm_pspec_,nspec_)
integer, intent(out)   :: nsymcrys_,lsplsymc_(48),lspnsymc_(48),invmap_(48),symlat_(3,3,48)
integer, intent(out)   :: ieqatom_(nmaxatm_pspec_,nspec_,48),ieqatom_inv_(nmaxatm_pspec_,nspec_,48)
real(8), intent(out)   :: vtlsymc_(3,48)
integer lspl,ispec
integer isym,jsym,iat,jat
integer ispl,jspl,ispl_inv
real(8) t1
real(8) vl(3),sl(3,3),si(3,3)
symtype=symtype_
nspecies=nspec_
maxatoms=nmaxatm_pspec_
natmmax=natmtot_
natmtot=natmtot_
do ispec=1,nspecies
  natoms_arr(ispec)=nat_per_spec_(ispec)
  spsymb(ispec)='XX'
end do
allocate(ieqatom(natmmax,nspecies,maxsymcrys))
allocate(ieqatom_inv(natmmax,nspecies,maxsymcrys))
! atomic positions in lattice coordinates
allocate(atposl(3,maxatoms,maxspecies))
! atomic positions in Cartesian coordinates
allocate(atposc(3,maxatoms,maxspecies))
atposl(:,:,:)=atml_(:,:,:)
! ELK works with columns of lattice vectors, thus->transpose
avec=transpose(avec_)
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! find Bravais lattice symmetries
call findsymlat
! find the crystal symmetries and shift atomic positions if required
call findsymcrys
! find equivalent atoms for direct and inverse crystal symmetry operation
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  sl(:,:)=dble(symlat(:,:,lspl))
  si(:,:)=dble(symlat(:,:,isymlat(lspl)))
  ! fix translation vectors such that the textbook {S,\tau}x=Sx+\tau would work
  !   (currently {S,\tau}x=S(x+\tau)
  vtlsymc(:,isym)=matmul(si,vtlsymc(:,isym))
  do jat=1,natoms_arr(1)
    vl(:)=atposl(:,jat,1)
    vl(:)=matmul(sl,vl)+vtlsymc(:,isym)
    call r3fracz05(epslat,vl)
    do iat=1,natoms_arr(1)
      t1=abs(atposl(1,iat,1)-vl(1))+abs(atposl(2,iat,1)-vl(2))+abs(atposl(3,iat,1)-vl(3))
      if (t1.lt.epslat) then
        ieqatom(jat,1,isym)=iat
        go to 9
      end if
    end do
    write(*,*) vtlsymc(:,isym)
    write(*,*) atposl(:,jat,1)
    write(*,*) sl(1,:)
    write(*,*) sl(2,:)
    write(*,*) sl(3,:)
    stop
    9 continue
  end do
  do jat=1,natoms_arr(1)
    vl(:)=atposl(:,jat,1)
    vl(:)=matmul(si,vl)+vtlsymc(:,isym)
    call r3fracz05(epslat,vl)
    do iat=1,natoms_arr(1)
      t1=abs(atposl(1,iat,1)-vl(1))+abs(atposl(2,iat,1)-vl(2))+abs(atposl(3,iat,1)-vl(3))
      if (t1.lt.epslat) then
        ieqatom_inv(jat,1,isym)=iat
        go to 10
      end if
    end do
    write(*,*) "Error(symmetry): atom not found"
    write(*,*) atposl(:,jat,1)
    write(*,*) si(1,:)
    write(*,*) si(2,:)
    write(*,*) si(3,:)
    stop
    10 continue
  end do
end do
do isym=1,nsymcrys
  ispl=lsplsymc(isym)
  ispl_inv=isymlat(ispl)
  do jsym=1,nsymcrys
    jspl=lsplsymc(jsym)
    if (jspl.eq.ispl_inv) then
      invmap(isym)=jsym
      exit
    end if
  end do
end do
if (sum(abs(atml_-atposl)).gt.epslat) then
  write(*,*) "Error(elk/symmetry.f90): atoms positions were changed"
  stop
end if
atml_        =atposl
ieqatom_     =ieqatom
ieqatom_inv_ =ieqatom_inv
nsymcrys_    =nsymcrys
lsplsymc_    =lsplsymc
lspnsymc_    =lspnsymc
invmap_      =invmap
symlat_      =symlat
ieqatom_     =ieqatom
ieqatom_inv_ =ieqatom_inv
vtlsymc_     =vtlsymc
call writesym
deallocate(atposl,atposc,ieqatom,ieqatom_inv)
return
end subroutine


