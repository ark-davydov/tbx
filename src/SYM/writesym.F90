
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writesym
! !INTERFACE:
subroutine writesym
! !USES:
use modelk
! !DESCRIPTION:
!   Outputs the Bravais, crystal and site symmetry matrices to files
!   {\tt SYMLAT.OUT}, {\tt SYMCRYS.OUT} and {\tt SYMSITE.OUT}, respectively.
!   Also writes out equivalent atoms and related crystal symmetries to
!   {\tt EQATOMS.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
integer i
integer isym,lspl,lspn
! output the Bravais lattice symmetries
open(50,file='SYMLAT.OUT',form='FORMATTED')
write(50,'(I4," : nsymlat")') nsymlat
do isym=1,nsymlat
  write(50,*)
  write(50,'(I4)') isym
  do i=1,3
    write(50,'(3I4)') symlat(i,:,isym)
  end do
end do
close(50)
! output the crystal symmetries
open(50,file='SYMCRYS.OUT',form='FORMATTED')
write(50,*)
write(50,'("(translation vectors and rotation matrices are in lattice &
 &coordinates)")')
write(50,*)
write(50,'(I4," : nsymcrys")') nsymcrys
do isym=1,nsymcrys
  write(50,*)
  write(50,'("Crystal symmetry : ",I4)') isym
  write(50,'(" spatial translation :")')
  write(50,'(3G18.10)') vtlsymc(:,isym)
  write(50,'(" spatial rotation :")')
  lspl=lsplsymc(isym)
  do i=1,3
    write(50,'(3I4)') symlat(i,:,lspl)
  end do
  write(50,'(" global spin rotation :")')
  lspn=lspnsymc(isym)
  do i=1,3
    write(50,'(3I4)') symlat(i,:,lspn)
  end do
end do
close(50)
return
end subroutine
!EOC

