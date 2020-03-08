
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmetry(symtype_,nspec_,nat_per_spec_,nmaxatm_pspec_,natmtot_,avec_,atml_,&
    nsymcrys_,lsplsymc_,lspnsymc_,invmap_,symlat_,ieqatom_,vtlsymc_,lattice_shift_,lwrite)
use modelk
implicit none
integer, intent(in)    :: symtype_,nspec_,nat_per_spec_(nspec_),nmaxatm_pspec_,natmtot_
real(8), intent(in)    :: avec_(3,3)
real(8), intent(inout) :: atml_(3,nmaxatm_pspec_,nspec_)
integer, intent(out)   :: nsymcrys_,lsplsymc_(48),lspnsymc_(48),invmap_(48),symlat_(3,3,48)
integer, intent(out)   :: ieqatom_(nmaxatm_pspec_,nspec_,48)
real(8), intent(out)   :: vtlsymc_(3,48)
real(8), intent(out)   :: lattice_shift_(3)
logical, intent(in)    :: lwrite
integer lspl,ispec,jspec
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
allocate(ieqatom(maxatoms,nspecies,48))
! atomic positions in lattice coordinates
allocate(atposl(3,maxatoms,nspecies))
! atomic positions in Cartesian coordinates
allocate(atposc(3,maxatoms,nspecies))
atposl(:,:,:)=atml_(:,:,:)
! ELK works with columns of lattice vectors, thus->transpose
avec=transpose(avec_)
! inverse of the lattice vector matrix
call r3minv(avec,ainv)
! find Bravais lattice symmetries
call findsymlat
! find the crystal symmetries and shift atomic positions if required
call findsymcrys(lwrite)
! find equivalent atoms for direct and inverse crystal symmetry operation
do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  sl(:,:)=dble(symlat(:,:,lspl))
  si(:,:)=dble(symlat(:,:,isymlat(lspl)))
  ! fix translation vectors such that the textbook {S,\tau}x=Sx+\tau would work
  !   (currently {S,\tau}x=S(x+\tau)
 ! vtlsymc(:,isym)=matmul(si,vtlsymc(:,isym))
  do jspec=1,nspecies
    do jat=1,natoms_arr(jspec)
      vl(:)=atposl(:,jat,jspec)
      !vl(:)=matmul(sl,vl)+vtlsymc(:,isym)
      vl(:)=matmul(sl,vl+vtlsymc(:,isym))
      call r3fracz05(epslat,vl)
      do ispec=1,nspecies
        do iat=1,natoms_arr(ispec)
          t1=sum(abs(atposl(:,iat,ispec)-vl(:)))
          if (t1.lt.epslat) then
            ieqatom(jat,jspec,isym)=iat
            go to 9
          end if
        end do
      end do
      write(*,*) "Error(symmetry): atom not found"
      write(*,*) vtlsymc(:,isym)
      write(*,*) atposl(:,jat,jspec)
      write(*,*) sl(1,:)
      write(*,*) sl(2,:)
      write(*,*) sl(3,:)
      stop
      9 continue
    end do
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
ieqatom_       =ieqatom
nsymcrys_      =nsymcrys
lsplsymc_      =lsplsymc
lspnsymc_      =lspnsymc
invmap_        =invmap
symlat_        =symlat
vtlsymc_       =vtlsymc
lattice_shift_ =atposl(:,1,1)-atml_(:,1,1)
if (lwrite) call writesym
deallocate(atposl,atposc,ieqatom)
return
end subroutine


