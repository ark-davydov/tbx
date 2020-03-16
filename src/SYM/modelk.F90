
module modelk

real(8), parameter :: epslat=1.d-6

real(8) avec(3,3)
real(8) ainv(3,3)
!----------------------------!
!     symmetry variables     !
!----------------------------!
! type of symmetry allowed for the crystal
!  0 : only the identity element is used
!  1 : full symmetry group is used
!  2 : only symmorphic symmetries are allowed
integer symtype
! number of Bravais lattice point group symmetries
integer nsymlat
! Bravais lattice point group symmetries
integer symlat(3,3,48)
! determinants of lattice symmetry matrices (1 or -1)
integer symlatd(48)
! index to inverses of the lattice symmetries
integer isymlat(48)
! lattice point group symmetries in Cartesian coordinates
real(8) symlatc(3,3,48)
integer invmap(48)
! tshift is .true. if atomic basis is allowed to be shifted
logical :: tshift
! tsyminv is .true. if the crystal has inversion symmetry
logical tsyminv
! maximum of symmetries allowed
integer, parameter :: maxsymcrys=48
! number of crystal symmetries
integer nsymcrys
! crystal symmetry translation vector in lattice coordinates
real(8) vtlsymc(3,maxsymcrys)
! tv0symc is .true. if the translation vector is zero
logical tv0symc(maxsymcrys)
! spatial rotation element in lattice point group for each crystal symmetry
integer lsplsymc(maxsymcrys)
! global spin rotation element in lattice point group for each crystal symmetry
integer lspnsymc(maxsymcrys)
! equivalent atom index for each crystal symmetry
integer, allocatable :: ieqatom(:,:,:)
! eqatoms(ia,ja,is) is .true. if atoms ia and ja are equivalent
logical, allocatable :: eqatoms(:,:,:)
! number of site symmetries
integer, allocatable :: nsymsite(:)
! site symmetry spatial rotation element in lattice point group
integer, allocatable :: lsplsyms(:,:)
! site symmetry global spin rotation element in lattice point group
integer, allocatable :: lspnsyms(:,:)

!--------------------------!
!     atomic variables     !
!--------------------------!
! maximum allowed species
integer, parameter :: maxspecies=100
! species symbol
character(256) spsymb(maxspecies)
! maximum allowed atoms per species
integer maxatoms
! number of species
integer nspecies
! number of atoms for each species
integer natoms_arr(maxspecies)
! maximum number of atoms over all the species
integer natmmax
! total number of atoms
integer natmtot
! molecule is .true. is the system is an isolated molecule
logical molecule
! primcell is .true. if primitive unit cell is to be found automatically
logical primcell
! index to atoms and species
integer, allocatable :: idxas(:,:)
! inverse atoms and species indices
integer, allocatable :: idxis(:)
integer, allocatable :: idxia(:)
! atomic positions in lattice coordinates
real(8), allocatable :: atposl(:,:,:)
! atomic positions in Cartesian coordinates
real(8), allocatable ::  atposc(:,:,:)
! magnitude of random displacements added to the atomic positions
real(8) rndatposc


CONTAINS 

subroutine r3frac(eps,v)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed iv, September 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
! local variables
integer i
do i=1,3
  v(i)=v(i)-int(v(i))
  if (v(i).lt.0.d0) v(i)=v(i)+1.d0
  if ((1.d0-v(i)).lt.eps) v(i)=0.d0
 ! if (v(i).ge.0.5d0) v(i)=v(i)-1.d0
  if (v(i).lt.eps) v(i)=0.d0
end do
return
end subroutine

subroutine r3frac1(eps,v)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed iv, September 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
! local variables
integer i
do i=1,3
  v(i)=v(i)-int(v(i))
  if (v(i).lt.0.d0) v(i)=v(i)+1.d0
  if ((1.d0-v(i)).lt.eps) v(i)=0.d0
  if (v(i).lt.eps) v(i)=0.d0
end do
return
end subroutine

subroutine r3fracz05(eps,v)
! !INPUT/OUTPUT PARAMETERS:
!   eps : zero component tolerance (in,real)
!   v   : input vector (inout,real(3))
! !DESCRIPTION:
!   Finds the fractional part of each component of a real 3-vector using the
!   function ${\rm frac}\,(x)=x-\lfloor x\rfloor$. A component is taken to be
!   zero if it lies within the intervals $[0,\epsilon)$ or $(1-\epsilon,1]$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Removed iv, September 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: eps
real(8), intent(inout) :: v(3)
! local variables
integer i
do i=1,3
  v(i)=v(i)-int(v(i))
  if (v(i).lt.0.d0) v(i)=v(i)+1.d0
  if ((1.d0-v(i)).lt.eps) v(i)=0.d0
  if (v(i).lt.eps) v(i)=0.d0
end do
if (v(3).gt.0.5d0) v(3)=v(3)-1.d0
return
end subroutine

subroutine r3cross(x,y,z)
! !INPUT/OUTPUT PARAMETERS:
!   x : input vector 1 (in,real(3))
!   y : input vector 2 (in,real(3))
!   z : output cross-product (out,real(3))
! !DESCRIPTION:
!   Returns the cross product of two real 3-vectors.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: x(3),y(3)
real(8), intent(out) :: z(3)
z(1)=x(2)*y(3)-x(3)*y(2)
z(2)=x(3)*y(1)-x(1)*y(3)
z(3)=x(1)*y(2)-x(2)*y(1)
return
end subroutine

subroutine r3minv(a,b)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
!   b : output matrix (out,real(3,3))
! !DESCRIPTION:
!   Computes the inverse of a real $3\times 3$ matrix.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3)
real(8), intent(out) :: b(3,3)
! local variables
real(8) t1
t1=a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
  -a(1,1)*a(2,3)*a(3,2)+a(1,1)*a(2,2)*a(3,3)-a(1,2)*a(2,1)*a(3,3)
if (abs(t1).lt.1.d-40) then
  write(*,*)
  write(*,'("Error(r3minv): singular matrix")')
  write(*,*)
  stop
end if
t1=1.d0/t1
b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*t1
b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))*t1
b(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))*t1
b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))*t1
b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))*t1
b(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))*t1
b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))*t1
b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))*t1
b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))*t1
return
end subroutine

subroutine r3mv(a,x,y)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,real(3,3))
!   x : input vector (in,real(3))
!   y : output vector (out,real(3))
! !DESCRIPTION:
!   Multiplies a real $3\times 3$ matrix with a vector.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3)
real(8), intent(in) :: x(3)
real(8), intent(out) :: y(3)
y(1)=a(1,1)*x(1)+a(1,2)*x(2)+a(1,3)*x(3)
y(2)=a(2,1)*x(1)+a(2,2)*x(2)+a(2,3)*x(3)
y(3)=a(3,1)*x(1)+a(3,2)*x(2)+a(3,3)*x(3)
return
end subroutine

integer function i3mdet(a)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,integer(3,3))
! !DESCRIPTION:
!   Returns the determinant of an integer $3\times 3$ matrix $A$.
!
! !REVISION HISTORY:
!   Created October 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: a(3,3)
i3mdet=a(1,1)*(a(2,2)*a(3,3)-a(3,2)*a(2,3)) &
      +a(2,1)*(a(3,2)*a(1,3)-a(1,2)*a(3,3)) &
      +a(3,1)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
return
end function

subroutine i3minv(a,b)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix (in,integer(3,3))
!   b : output matrix (in,integer(3,3))
! !DESCRIPTION:
!   Computes the inverse of a integer $3\times 3$ matrix: $B=A^{-1}$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: a(3,3)
integer, intent(out) :: b(3,3)
! local variables
integer m
m=a(1,2)*a(2,3)*a(3,1)-a(1,3)*a(2,2)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
 -a(1,1)*a(2,3)*a(3,2)+a(1,1)*a(2,2)*a(3,3)-a(1,2)*a(2,1)*a(3,3)
if ((m.ne.1).and.(m.ne.-1)) then
  write(*,*)
  write(*,'("Error(i3minv): cannot invert matrix")')
  write(*,'(" Determinant : ",I8)') m
  write(*,*)
  stop
end if
b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*m
b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))*m
b(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))*m
b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))*m
b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))*m
b(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))*m
b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))*m
b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))*m
b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))*m
return
end subroutine

subroutine r3mtm(a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix 1 (in,real(3,3))
!   b : input matrix 2 (in,real(3,3))
!   c : output matrix (out,real(3,3))
! !DESCRIPTION:
!   Multiplies the transpose of one real $3\times 3$ matrix with another.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3),b(3,3)
real(8), intent(out) :: c(3,3)
c(1,1)=a(1,1)*b(1,1)+a(2,1)*b(2,1)+a(3,1)*b(3,1)
c(2,1)=a(1,2)*b(1,1)+a(2,2)*b(2,1)+a(3,2)*b(3,1)
c(3,1)=a(1,3)*b(1,1)+a(2,3)*b(2,1)+a(3,3)*b(3,1)
c(1,2)=a(1,1)*b(1,2)+a(2,1)*b(2,2)+a(3,1)*b(3,2)
c(2,2)=a(1,2)*b(1,2)+a(2,2)*b(2,2)+a(3,2)*b(3,2)
c(3,2)=a(1,3)*b(1,2)+a(2,3)*b(2,2)+a(3,3)*b(3,2)
c(1,3)=a(1,1)*b(1,3)+a(2,1)*b(2,3)+a(3,1)*b(3,3)
c(2,3)=a(1,2)*b(1,3)+a(2,2)*b(2,3)+a(3,2)*b(3,3)
c(3,3)=a(1,3)*b(1,3)+a(2,3)*b(2,3)+a(3,3)*b(3,3)
return
end subroutine

subroutine r3mm(a,b,c)
! !INPUT/OUTPUT PARAMETERS:
!   a : input matrix 1 (in,real(3,3))
!   b : input matrix 2 (in,real(3,3))
!   c : output matrix (out,real(3,3))
! !DESCRIPTION:
!   Multiplies two real $3\times 3$ matrices.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: a(3,3),b(3,3)
real(8), intent(out) :: c(3,3)
c(1,1)=a(1,1)*b(1,1)+a(1,2)*b(2,1)+a(1,3)*b(3,1)
c(2,1)=a(2,1)*b(1,1)+a(2,2)*b(2,1)+a(2,3)*b(3,1)
c(3,1)=a(3,1)*b(1,1)+a(3,2)*b(2,1)+a(3,3)*b(3,1)
c(1,2)=a(1,1)*b(1,2)+a(1,2)*b(2,2)+a(1,3)*b(3,2)
c(2,2)=a(2,1)*b(1,2)+a(2,2)*b(2,2)+a(2,3)*b(3,2)
c(3,2)=a(3,1)*b(1,2)+a(3,2)*b(2,2)+a(3,3)*b(3,2)
c(1,3)=a(1,1)*b(1,3)+a(1,2)*b(2,3)+a(1,3)*b(3,3)
c(2,3)=a(2,1)*b(1,3)+a(2,2)*b(2,3)+a(2,3)*b(3,3)
c(3,3)=a(3,1)*b(1,3)+a(3,2)*b(2,3)+a(3,3)*b(3,3)
return
end subroutine

subroutine sphcrd(v,r,tp)
! !INPUT/OUTPUT PARAMETERS:
!   v  : input vector (in,real(3))
!   r  : length of v (out,real)
!   tp : (theta, phi) coordinates (out,real(2))
! !DESCRIPTION:
!   Returns the spherical coordinates $(r,\theta,\phi)$ of a vector
!   $$ {\bf v}=(r\sin(\theta)\cos(\phi), r\sin(\theta)\sin(\phi),
!    r\cos(\theta)). $$
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: v(3)
real(8), intent(out) :: r,tp(2)
! local variables
real(8), parameter :: pi=3.1415926535897932385d0
real(8), parameter :: eps=1.d-14
real(8) t1
r=sqrt(v(1)**2+v(2)**2+v(3)**2)
if (r.gt.eps) then
  t1=v(3)/r
  if (t1.ge.1.d0) then
    tp(1)=0.d0
  else if (t1.le.-1.d0) then
    tp(1)=pi
  else
    tp(1)=acos(t1)
  end if
  if ((abs(v(1)).gt.eps).or.(abs(v(2)).gt.eps)) then
    tp(2)=atan2(v(2),v(1))
  else
    tp(2)=0.d0
  end if
else
  tp(1)=0.d0
  tp(2)=0.d0
end if
return
end subroutine


end module
