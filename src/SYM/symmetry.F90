
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symmetry(symtype_,nspec_,nat_per_spec_,nmaxatm_pspec_,natmtot_,avec_,atml_,&
    nsymcrys_,lsplsymc_,lspnsymc_,invmap_,symlat_,vtlsymc_,lattice_shift_,shifted_,tshift_,lwrite)
use modelk
implicit none
integer, intent(in)    :: symtype_,nspec_,nat_per_spec_(nspec_),nmaxatm_pspec_,natmtot_
real(8), intent(in)    :: avec_(3,3)
real(8), intent(in) :: atml_(3,nmaxatm_pspec_,nspec_)
integer, intent(out)   :: nsymcrys_,lsplsymc_(48),lspnsymc_(48),invmap_(48),symlat_(3,3,48)
real(8), intent(out)   :: vtlsymc_(3,48)
real(8), intent(inout) :: lattice_shift_(3,nmaxatm_pspec_,nspec_)
logical, intent(out)   :: shifted_
logical, intent(in)    :: tshift_
logical, intent(in)    :: lwrite
integer ispec
integer isym,jsym
integer ispl,jspl,ispl_inv
tshift=tshift_
symtype=symtype_
nspecies=nspec_
maxatoms=nmaxatm_pspec_
natmmax=natmtot_
natmtot=natmtot_
do ispec=1,nspecies
  natoms_arr(ispec)=nat_per_spec_(ispec)
  spsymb(ispec)='XX'
end do
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
call findsymcrys(lwrite,shifted_)
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
lattice_shift_ =atposl-atml_
nsymcrys_      =nsymcrys
lsplsymc_      =lsplsymc
lspnsymc_      =lspnsymc
invmap_        =invmap
symlat_        =symlat
vtlsymc_       =vtlsymc

if (lwrite) call writesym
deallocate(atposl,atposc)
return
end subroutine


