
! Copyright (C) 2007-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: findsym
! !INTERFACE:
subroutine findsym(apl1,apl2,nsym,lspl,lspn,iea)
! !USES:
use modelk
! !INPUT/OUTPUT PARAMETERS:
!   apl1 : first set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   apl2 : second set of atomic positions in lattice coordinates
!          (in,real(3,maxatoms,maxspecies))
!   nsym : number of symmetries (out,integer)
!   lspl : spatial rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   lspn : spin rotation element in lattice point group for each symmetry
!          (out,integer(48))
!   iea  : equivalent atom index for each symmetry
!          (out,integer(iea(natmmax,nspecies,48))
! !DESCRIPTION:
!   Finds the symmetries which rotate one set of atomic positions into another.
!   Both sets of positions differ only by a translation vector and have the same
!   muffin-tin magnetic fields (stored in the global array {\tt bfcmt}). Any
!   symmetry element consists of a spatial rotation of the atomic position
!   vectors followed by a global magnetic rotation: $\{\alpha_S|\alpha_R\}$. In
!   the case of spin-orbit coupling $\alpha_S=\alpha_R$. The symmetries are
!   returned as indices of elements in the Bravais lattice point group. An
!   index to equivalent atoms is stored in the array {\tt iea}.
!
! !REVISION HISTORY:
!   Created April 2007 (JKD)
!   Fixed use of proper rotations for spin, February 2008 (L. Nordstrom)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: apl1(3,maxatoms,maxspecies)
real(8), intent(in) :: apl2(3,maxatoms,maxspecies)
integer, intent(out) :: nsym
integer, intent(out) :: lspl(48)
integer, intent(out) :: lspn(48)
integer, intent(out) :: iea(natmmax,nspecies,48)
! local variables
integer isym,jsym
integer is,ia,ja
real(8) sl(3,3),v(3),t1
! allocatable arrays
integer, allocatable :: jea(:,:)
real(8), allocatable :: apl3(:,:)
allocate(jea(natmmax,nspecies))
allocate(apl3(3,natmmax))
nsym=0
! loop over lattice symmetries (spatial rotations)
do isym=1,nsymlat
! make real copy of lattice rotation symmetry
  sl(:,:)=dble(symlat(:,:,isym))
! loop over species
  do is=1,nspecies
! map apl1 coordinates to [0,1) and store in apl3
    do ia=1,natoms_arr(is)
      apl3(:,ia)=apl1(:,ia,is)
      call r3fracz05(epslat,apl3(:,ia))
    end do
    do ja=1,natoms_arr(is)
! apply lattice symmetry to atomic positions
      v(:)=sl(:,1)*apl2(1,ja,is)+sl(:,2)*apl2(2,ja,is)+sl(:,3)*apl2(3,ja,is)
! map coordinates to [0,1)
      call r3fracz05(epslat,v)
! check if atomic positions are invariant
      do ia=1,natoms_arr(is)
        t1=abs(apl3(1,ia)-v(1))+abs(apl3(2,ia)-v(2))+abs(apl3(3,ia)-v(3))
        if (t1.lt.epslat) then
! equivalent atom index
          jea(ia,is)=ja
          goto 10
        end if
      end do
! not invariant so try new spatial rotation
      goto 40
10 continue
    end do
  end do
! all atomic positions invariant at this point
  jsym=1
! everything invariant so add symmetry to set
  nsym=nsym+1
  lspl(nsym)=isym
  lspn(nsym)=jsym
  do is=1,nspecies
    do ia=1,natoms_arr(is)
      iea(ia,is,nsym)=jea(ia,is)
    end do
  end do
40 continue
! end loop over spatial rotations
end do
deallocate(jea,apl3)
return
end subroutine
!EOC

