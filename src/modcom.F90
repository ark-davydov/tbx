
module modcom
#ifdef MPI
  use mpi
#endif
implicit none


logical, parameter :: centered_kgrid=.false.
logical, parameter :: centered_qgrid=.false.
integer, parameter :: dp=kind(0.d0)
integer, parameter :: NDIM=3
! careful, this is not the constant, it will be re-set in input file read
integer :: NDIM_COUL=3
integer, parameter :: ZAXIS=3
integer, parameter :: nmaxspec=10
integer, parameter :: nmaxvert=100
real(dp), parameter :: abohr=0.52917721067_dp
real(dp), parameter :: twothrd=0.666666666666666666666666_dp
real(dp), parameter :: onethrd=0.333333333333333333333333_dp
real(dp), parameter :: pi=3.141592653589793115997963468544185161590576171875_dp
real(dp), parameter :: twopi=6.283185307179586231995926937088370323181152343750_dp
real(dp), parameter :: fourpi=12.5663706143591724639918538741767406463623046875_dp
real(dp), parameter :: pihalf=1.570796326794896557998981734272092580795288085937_dp
real(dp), parameter :: sqrtpi=1.77245385091_dp
real(dp), parameter :: sqrt2=1.41421356237
real(dp), parameter :: sqrt12=0.707106781187
real(dp), parameter :: epslat=1.e-6_dp
real(dp), parameter :: epsengy=1.e-7_dp
real(dp), parameter :: Hartree_to_ev=27.211386245988_dp
real(dp), parameter :: CoulombForceConstant=abohr*Hartree_to_ev
complex(dp), parameter :: cmplx_0=cmplx(0._dp,0._dp,kind=dp)
complex(dp), parameter :: cmplx_1=cmplx(1._dp,0._dp,kind=dp)
real(dp) :: graphene_lvec_length=2.46_dp
real(dp) :: graphene_cc_distance=1.42_dp
real(dp) :: tbg_aa_distance=3.60_dp
real(dp) :: tbg_ab_distance=3.35_dp

! Error message from MPI
integer mpi_err
! MPI communicator for main code
integer mpi_com
! number of MPI processes
integer np_mpi
! local MPI process number
integer lp_mpi
! mp_mpi is .true. if the local MPI process is the master (0)
logical mp_mpi




contains 

  ! split a string into 2 either side of a delimiter token
  SUBROUTINE split_string(instring, string1, string2, crit)
    CHARACTER :: crit
    CHARACTER(256) :: instring
    CHARACTER(256),INTENT(OUT):: string1,string2
    INTEGER :: index

    instring = TRIM(adjustl(instring))

    index = SCAN(instring,crit)
!    if (index.eq.0) then
!      string1=instring
!      string2=''
!    else
      string1 = instring(1:index-1)
      string2 = instring(index+1:)
!    end if

  END SUBROUTINE split_string

subroutine  throw(place,what)
   character(len=*), intent(in) :: place
   character(len=*), intent(in) :: what
#ifdef MPI
   integer error_msg_mpi
#endif
   write(*,*) "ERROR[ ",place," ] ","MESSAGE: ",what
#ifdef MPI
   call mpi_abort(mpi_com,error_msg_mpi,mpi_err)
#endif
   stop
end subroutine

subroutine  info(place,what)
   character(len=*), intent(in) :: place
   character(len=*), intent(in) :: what
   if (mp_mpi) write(*,*) "INFO[ ",place," ] ",": ",what
end subroutine

subroutine  message(what)
   character(len=*), intent(in) :: what
   if (mp_mpi) write(*,*) trim(adjustl(what))
end subroutine


subroutine dmatrix_inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine dmatrix_inverse

function rotation3D(phi,ax) result (rot)
  real(dp), intent(in) :: phi,ax(3)
  real(dp) rot(3,3)
  real(dp) norm,ux,uy,uz,co1,si,co
  norm=sqrt(sum(ax**2))
  ux=ax(1)/norm
  uy=ax(2)/norm
  uz=ax(3)/norm
  co=cos(phi)
  co1=1_dp-co
  si=sin(phi)
  ! wiki
  rot(1,:) = (/ co+ux**2 *co1   , ux*uy*co1-uz*si , ux*uz*co1+uy*si /)
  rot(2,:) = (/ uy*ux*co1+uz*si , co+uy**2 *co1   , uy*uz*co1-ux*si /)
  rot(3,:) = (/ uz*ux*co1-uy*si , uz*uy*co1+ux*si , co+uz**2 *co1   /)
end function

subroutine eigenv_problem(ndim,matrix,eval)
integer, intent(in) :: ndim
complex(dp), intent(inout) :: matrix(:,:)
real(dp), intent(out) :: eval(:)
integer info,ld,lwork
real(dp), allocatable :: rwork(:)
complex(dp), allocatable :: work(:)
ld=ndim
lwork=2*ndim-1
allocate(rwork(3*ndim-2),work(lwork))
! diagonalise the matrix in double precision
call zheev('V','U',ndim,matrix,ld,eval,work,lwork,rwork,info)
if (info.ne.0) call throw("modcom%eigenv_problem()","zheev returned nonzero")
deallocate(rwork,work)
end subroutine


  !> the IL-th through IU-th eigenvalues will be found.
  subroutine zheevx_pack(JOBZ, UPLO, N, il, iu, A, eigval, eigvec)
     implicit none
     !> inout variables
     character(1), intent(in) :: JOBZ ! 'V' compute eigenvectors and eigenvalues; 'N' only eigenvalues
     character(1), intent(in) :: UPLO
     integer, intent(in) :: N  !> dimension of matrix A
     integer, intent(in) :: il !> lowest band to be calculated
     integer, intent(in) :: iu !> highest band to be calculated
     real(dp), intent(inout) :: eigval(iu-il+1) !> eigenvalues
     complex(dp), intent(inout) :: A(N, N) !> the input hermite complex matrix
     complex(dp), optional, intent(inout) :: eigvec(N, iu-il+1) !> eigenvectors

     !> local variables
     real(dp), allocatable :: eigenvalues(:)
     integer , allocatable :: iwork(:)
     integer , allocatable :: ifail(:)
     real(dp), allocatable :: rwork(:)
     complex(dp), allocatable :: work(:)

     integer :: mdim
     integer :: lwork
     integer :: info
     real(dp) :: vl  !> not referenced in this subroutine
     real(dp) :: vu  !> not referenced in this subroutine

     real(dp), parameter :: abstol= 1D-10

     !real(dp), external :: DLAMCH
     !abstol= 2*DLAMCH('S')

     mdim= iu-il+1
     lwork= 64*N

     allocate(ifail(N))
     allocate(iwork(5*N))
     allocate(rwork(7*N))
     allocate(work(lwork))
     allocate(eigenvalues(N))
     ifail = 0
     iwork = 0
     rwork = 0_dp
     work = 0_dp

     eigenvalues= 0_dp
     eigvec= 0_dp
     vl=0_dp; vu=0_dp

     call zheevx(JOBZ,'I',UPLO,N,A,N,vl,vu,il,iu,abstol,&
         mdim,eigenvalues,eigvec,N,work,lwork,rwork,iwork,ifail,info)

     if (info/=0) then
        print *, ' Error info in zheevx: ', info
        if (info.gt.0) then
          print *, ' zheevx failed to converge ',info, ' eigenvalues' 
          print *, ' increase the number of states to compute to a full spectrum, in order to use zheev '
        end if
        stop ' Error happens in zheev_pack'
     endif

     eigval(1:mdim)= eigenvalues(1:mdim)

     deallocate(eigenvalues, work, iwork, ifail, rwork)

     return
  end subroutine  zheevx_pack


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sortidx
! !INTERFACE:
subroutine sortidx(n,a,idx)
! !INPUT/OUTPUT PARAMETERS:
!   n   : number of elements in array (in,integer)
!   a   : real array (in,real(n))
!   idx : permutation index (out,integer(n))
! !DESCRIPTION:
!   Finds the permutation index {\tt idx} which sorts the real array {\tt a}
!   into ascending order. No sorting of the array {\tt a} itself is performed.
!   Uses the heapsort algorthim.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Included tolerance eps, April 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: a(n)
integer, intent(out) :: idx(n)
! local variables
integer i,j,k,l,m
! tolerance for deciding if one number is smaller than another
real(8), parameter :: eps=1.e-14_dp
if (n.le.0) then
  write(*,*)
  write(*,'("Error(sortidx): n <= 0 : ",I8)') n
  write(*,*)
  stop
end if
do i=1,n
  idx(i)=i
end do
if (n.eq.1) return
l=n/2+1
k=n
10 continue
if (l.gt.1) then
  l=l-1
  m=idx(l)
else
  m=idx(k)
  idx(k)=idx(1)
  k=k-1
  if (k.eq.1) then
    idx(1)=m
    return
  end if
end if
i=l
j=l+l
20 continue
if (j.le.k) then
  if (j.lt.k) then
    if (a(idx(j)).lt.a(idx(j+1))+eps) j=j+1
  end if
  if (a(m).lt.a(idx(j))+eps) then
    idx(i)=idx(j)
    i=j
    j=j+j
  else
    j=k+1
  end if
  goto 20
end if
idx(i)=m
goto 10
end subroutine
!EOC

subroutine io_evec(ik,action,fname,norb,nstates,evec)
integer, intent(in) :: ik,norb,nstates
character(len=*), intent(in) :: action,fname
complex(dp), intent(inout) :: evec(norb,nstates)
! local
logical exs
character(len=4) num
integer iostat,recl,unt
unt=1000+ik
inquire(iolength=recl) evec(:,:)
write(num,'(I4.4)') ik
if (action.eq."write") then
  call system("mkdir -p _evec")
  open(unt,file="_evec/"//trim(adjustl(fname))//num,form='unformatted',access='direct',&
    action=trim(adjustl(action)),recl=recl,iostat=iostat)
  if (iostat.ne.0) call throw("modcom%io_evec","problem reading evec file")
  ! write a record
  write(unt,rec=1,iostat=iostat) evec(:,:)
  if (iostat.ne.0) call throw("modcom%io_evec","problem writing evec file")
  close(unt)
else if (action.eq."read") then
  inquire(file="_evec/",exist=exs)
  if (.not.exs) call throw("modcom%io_evec","directory _evec does not exist")
  open(unt,file="_evec/"//trim(adjustl(fname))//num,form='unformatted',access='direct',&
    action=trim(adjustl(action)),recl=recl,iostat=iostat)
  if (iostat.ne.0) call throw("modcom%io_evec","problem reading evec file")
  ! read a record
  read(unt,rec=1,iostat=iostat) evec(:,:)
  if (iostat.ne.0) call throw("modcom%io_evec","problem reading evec file")
  close(unt)
else
  call throw("modcom%io_evec()","unknown read/write action")
end if
end subroutine

subroutine io_eval(unt,action,fname,wannier,nstates,nkpt,efermi,vkl,eval)
integer, intent(in) :: unt,nstates,nkpt
character(len=*), intent(in) :: action,fname
logical, intent(in) :: wannier
real(dp), intent(in) :: efermi
real(dp), intent(in) :: vkl(NDIM,nkpt)
real(dp), intent(inout) :: eval(nstates,nkpt)
! local
integer ik,ip,ist,jst,nk,nst
real(dp) vpl(NDIM),ef,t1
if (trim(adjustl(action)).ne."write".and.trim(adjustl(action)).ne."read") &
  call throw("modcom%io_eval()","unknown read/write action")
open(unt,file=trim(adjustl(fname)),action=trim(adjustl(action)))
if (.not.wannier) then
  if (trim(adjustl(action)).eq."write") then
    write(unt,'(2I6)') nkpt,nstates
  else
    read(unt,*) nk,nst
    if (nk.ne.nkpt) call throw("modcom%io_eval()","number of k-points in eigenvalue file is not what expected")
    if (nst.ne.nstates) call throw("modcom%io_eval()","number of states in eigenvalue file is not what expected")
  end if
end if
do ik=1,nkpt
  if (.not.wannier) then
     if (trim(adjustl(action)).eq."write") then
       write(unt,'(I6,10G18.10)') ik,vkl(:,ik)
     else
       read(unt,*) ip,vpl(:)
       if (ip.ne.ik) call throw("modcom%io_eval()","k-point sequence in eigenvalue file is not what expected")
       if (sum(abs(vpl(:)-vkl(:,ik))).gt.epslat) &
         call throw("modcom%io_eval()","k-point sequence in eigenvalue file is not what expected")
     end if
  end if
  do ist=1,nstates
     if (.not.wannier) then
       if (trim(adjustl(action)).eq."write") then
          write(unt,'(I6,2G18.10)') ist,eval(ist,ik),eval(ist,ik)-efermi
       else
          read(unt,*) jst,eval(ist,ik),t1
          if (ist.ne.jst) call throw("modcom%io_eval()","state index sequence in eigenvalue file is not what expected")
          ef=eval(ist,ik)-t1
          if (abs(ef-efermi).gt.epslat) call throw("modcom%io_eval()","Fermi energy derived from the input is different")
       end if
     else
       if (trim(adjustl(action)).ne."write") &
         call throw("modcom%io_eval()"," this subroutine is not supposed to read wannier files")
       write(unt,'(2I6,G18.10)') ist,ik,eval(ist,ik)-efermi
     end if
   end do
end do
close(unt)
end subroutine

real(dp) function gauss(x,sig)
implicit none
real(dp), intent(in) :: x,sig
real(dp) t1
t1=(x/sig)**2
gauss=exp(-t1)/sig/sqrtpi
end function

subroutine clmr(ll,mm,l,m1,m2,sgn)
integer, intent(in) :: ll,mm
integer, intent(out) :: l,m1,m2,sgn
if (ll.lt.10) call throw("modcom%clmr","ll should be > 10 here")
l=ll-10
sgn=sign(1,mm)
m1=abs(mm)/10
m2=mod(abs(mm),10)
return
end subroutine

function string_to_lmr(str_orbital) result(lmr)
character(len=*), intent(in) :: str_orbital
integer lmr(2)
select case(trim(adjustl(str_orbital))) 
case('s')          ;  lmr=(/0,1/)
case('pz')         ;  lmr=(/1,1/)
case('px')         ;  lmr=(/1,2/)
case('py')         ;  lmr=(/1,3/)
case('dz2')        ;  lmr=(/2,1/)
case('dxz')        ;  lmr=(/2,2/)
case('dyz')        ;  lmr=(/2,3/)
case('dx2-y2')     ;  lmr=(/2,4/)
case('dxy')        ;  lmr=(/2,5/)
case('fz3')        ;  lmr=(/3,1/)
case('fxz2')       ;  lmr=(/3,2/)
case('fyz2')       ;  lmr=(/3,3/)
case('fz(x2-y2)')  ;  lmr=(/3,4/)
case('fxyz')       ;  lmr=(/3,5/)
case('fx(x2-3y2)') ;  lmr=(/3,6/)
case('fy(3x2-y2)') ;  lmr=(/3,7/)
case('pp')         ;  lmr=(/11,23/)
case('pm')         ;  lmr=(/11,-23/)
! hybridized
case('sp-1')       ;  lmr=(/-1,1/)
case('sp-2')       ;  lmr=(/-1,2/)
case('sp2-1')      ;  lmr=(/-2,1/)
case('sp2-2')      ;  lmr=(/-2,2/)
case('sp2-3')      ;  lmr=(/-2,3/)
case('sp3-1')      ;  lmr=(/-3,1/)
case('sp3-2')      ;  lmr=(/-3,2/)
case('sp3-3')      ;  lmr=(/-3,3/)
case('sp3-4')      ;  lmr=(/-3,4/)
case('sp3d-1')     ;  lmr=(/-4,1/)
case('sp3d-2')     ;  lmr=(/-4,2/)
case('sp3d-3')     ;  lmr=(/-4,3/)
case('sp3d-4')     ;  lmr=(/-4,4/)
case('sp3d-5')     ;  lmr=(/-4,5/)
case('sp3d2-1')    ;  lmr=(/-5,1/)
case('sp3d2-2')    ;  lmr=(/-5,2/)
case('sp3d2-3')    ;  lmr=(/-5,3/)
case('sp3d2-4')    ;  lmr=(/-5,4/)
case('sp3d2-5')    ;  lmr=(/-5,5/)
case('sp3d2-6')    ;  lmr=(/-5,6/)
case default       ;  call throw("modcom%string_to_lmr","unknown orbital string")
end select
end function

function lmr_to_string(lmr) result(str_orbital)
integer, intent(in) :: lmr(2)
character(len=10) str_orbital
integer ll,mr
ll=lmr(1)
mr=lmr(2)
select case(ll)
case(11)
  select case(mr)
  case(23)          ; str_orbital='pp'
  case(-23)          ; str_orbital='pm'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(0)
  select case(mr)
  case(1)          ; str_orbital='s'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(1)
  select case(mr)
  case(1)          ; str_orbital='pz'
  case(2)          ; str_orbital='px'
  case(3)          ; str_orbital='py'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(2)
  select case(mr)
  case(1)          ; str_orbital='dz2'
  case(2)          ; str_orbital='dxz'
  case(3)          ; str_orbital='dyz'
  case(4)          ; str_orbital='dx2-y2'
  case(5)          ; str_orbital='dxy'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(3)
  select case(mr)
  case(1)          ; str_orbital='fz3'
  case(2)          ; str_orbital='fxz2'
  case(3)          ; str_orbital='fyz2'
  case(4)          ; str_orbital='fz(x2-y2)'
  case(5)          ; str_orbital='fxyz'
  case(6)          ; str_orbital='fx(x2-3y2)'
  case(7)          ; str_orbital='fy(3x2-y2)'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(-1)
  select case(mr)
  case(1)          ; str_orbital='sp-1'
  case(2)          ; str_orbital='sp-2'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(-2)
  select case(mr)
  case(1)          ; str_orbital='sp2-1'
  case(2)          ; str_orbital='sp2-2'
  case(3)          ; str_orbital='sp2-3'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(-3)
  select case(mr)
  case(1)          ; str_orbital='sp3-1'
  case(2)          ; str_orbital='sp3-2'
  case(3)          ; str_orbital='sp3-3'
  case(4)          ; str_orbital='sp3-4'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(-4)
  select case(mr)
  case(1)          ; str_orbital='sp3d-1'
  case(2)          ; str_orbital='sp3d-2'
  case(3)          ; str_orbital='sp3d-3'
  case(4)          ; str_orbital='sp3d-4'
  case(5)          ; str_orbital='sp3d-5'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case(-5)
  select case(mr)
  case(1)          ; str_orbital='sp3d2-1'
  case(2)          ; str_orbital='sp3d2-2'
  case(3)          ; str_orbital='sp3d2-3'
  case(4)          ; str_orbital='sp3d2-4'
  case(5)          ; str_orbital='sp3d2-5'
  case(6)          ; str_orbital='sp3d2-6'
  case default     ; call throw("modcom%lmr_to_string","unknown mr")
  end select
case default 
  call throw("modcom%lmr_to_string","unknown l ")
end select
end function

real(dp) function pwave_ovlp(dc)
real(dp), intent(in) :: dc
real(dp), parameter :: charge_pz=3.18_dp
pwave_ovlp=( 1._dp/( 1._dp+(dc*abohr/charge_pz)**2 ) )**3
end function

subroutine find_degroups(nx,xx,ngr,idx,eps)
integer :: nx,ngr
real(dp) :: xx(nx),eps
integer, intent(out) :: idx(nx,2)
integer nlast,ix
ngr=1
nlast=0
idx(ngr,1)=1
idx(ngr,2)=1
do ix=2,nx
  if (abs(xx(ix)-xx(ix-1)).lt.eps) then
     nlast=nlast+1
     idx(ngr,2)=idx(ngr,1)+nlast
  else
    ngr=ngr+1
    nlast=0
    idx(ngr,1)=ix
    idx(ngr,2)=ix
  end if
end do
end subroutine
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: r3cross
! !INTERFACE:
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
!EOC


! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gengclq
! !INTERFACE:
subroutine gengclq(ngridq,bvec,gclq)
! !DESCRIPTION:

!   The Fock matrix elements
!   $$ V_{ij{\bf k}}\equiv\sum_{l{\bf k'}}\int
!    \frac{\Psi^{\dag}_{i{\bf k}}({\bf r})\cdot\Psi_{l{\bf k}'}({\bf r})
!    \Psi^{\dag}_{l{\bf k}'}({\bf r}')\cdot\Psi_{j{\bf k}}({\bf r}')}
!    {|{\bf r}-{\bf r'}|}\,d^3r\,d^3r' $$
!   contain a divergent term in the sum over ${\bf k}'$ which behaves as
!   $1/q^2$, where ${\bf q}\equiv{\bf k}-{\bf k}'$ is in the first Brillouin
!   zone. The resulting convergence with respect to the number of discrete
!   $q$-points, $N_q$, is very slow. This routine computes the regularised
!   Coulomb Green's function
!   \begin{align}
!    g({\bf q}_i)=\frac{4\pi}{V}\int_{V_i}\frac{1}{q^2}\,d^3q,
!   \end{align}
!   where the integral is over the small parallelepiped with volume
!   $V=\Omega_{\rm BZ}/N_q$ and centered on the discrete point ${\bf q}_i$.
!   This dramatically increases the rate of convergence of methods which involve
!   a summation over the $1/q^2$ part of the Coulomb interaction. The above
!   integral is evaluated numerically on increasingly finer grids and then
!   extrapolated to the continuum.
!
! !REVISION HISTORY:
!   Created August 2004 (JKD,SS)
!   Changed from genwiq2, July 2017 (JKD)
!   Adopted to TBX code June 2020 (Arkadiy Davydov)
!EOP
!BOC
implicit none
integer, intent(in) :: ngridq(NDIM)
real(dp), intent(in) :: bvec(NDIM,NDIM)
real(dp), intent(out) :: gclq
! local variables
integer, parameter :: np=5
integer, parameter :: ns0=10,nss=20
integer ns,i1,i2,i3,ip
real(8) d(NDIM),sum,t1,t2
real(8) v1(NDIM),v2(NDIM),v3(NDIM)
real(8) xa(np),ya(np),c(np)
if (NDIM_COUL.eq.1) then
   call throw("modcom%gengclq()","in 1D Fuorier transform of 1/r is a step function, dont know how to deal with that yet")
else if (NDIM_COUL.gt.3) then
   call throw("modcom%gengclq()","NDIM_COUL is greater than 3")
end if
! loop over different subdivisions
ns=ns0
do ip=1,np
! subdivision vectors in lattice coordinates
  d(:)=1._dp/dble(ngridq(:)*2*ns)
! compute the integral of 1/q^2
  sum=0._dp
  if (NDIM_COUL.eq.2) then
    do i1=-ns,ns-1
      t1=dble(i1)*d(1)
      v1(:)=t1*bvec(1,:)
      do i2=-ns,ns-1
        t1=dble(i2)*d(2)
        v2(:)=v1(:)+t1*bvec(2,:)
        t2=v2(1)**2+v2(2)**2
        ! V~1/|q| in 2D
        if (t2.gt.1.e-14_dp) sum=sum+0.5_dp/sqrt(t2)
      end do
    end do
  else if (NDIM_COUL.eq.3) then
    do i1=-ns,ns-1
      t1=dble(i1)*d(1)
      v1(:)=t1*bvec(1,:)
      do i2=-ns,ns-1
        t1=dble(i2)*d(2)
        v2(:)=v1(:)+t1*bvec(2,:)
        do i3=-ns,ns-1
          t1=dble(i3)*d(3)
          v3(:)=v2(:)+t1*bvec(3,:)
          t2=v3(1)**2+v3(2)**2+v3(3)**2
          ! V~1/|q|**2 in 3D
          if (t2.gt.1.e-14_dp) sum=sum+1._dp/t2
        end do
      end do
    end do
  end if
  t1=1._dp/dble(2*ns)
  xa(ip)=t1
  ya(ip)=fourpi*CoulombForceConstant*sum*t1**NDIM_COUL
  ! increment number of subdivisions
  ns=ns+nss
end do
! extrapolate the volume element to zero with a polynomial
gclq=polynm(0,np,xa,ya,0._dp,c)
return
end subroutine
!EOC


! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: polynm
! !INTERFACE:
real(8) function polynm(m,np,xa,ya,x,c)
! !INPUT/OUTPUT PARAMETERS:
!   m  : order of derivative (in,integer)
!   np : number of points to fit (in,integer)
!   xa : abscissa array (in,real(np))
!   ya : ordinate array (in,real(np))
!   x  : evaluation abscissa (in,real)
! !DESCRIPTION:
!   Fits a polynomial of order $n_p-1$ to a set of $n_p$ points. If $m\ge 0$ the
!   function returns the $m$th derviative of the polynomial at $x$, while for
!   $m<0$ the integral of the polynomial from the first point in the array to
!   $x$ is returned.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! argmuments
integer, intent(in) :: m,np
real(8), intent(in) :: xa(np),ya(np)
real(8), intent(in) :: x
real(8), intent(inout) :: c(np)
! local variables
integer i,j,k
real(8) x0,x1,x2,x3,y0,y1,y2,y3
real(8) t0,t1,t2,t3,t4,t5,t6
real(8) c1,c2,c3,sum
! fast evaluations for small np
select case(np)
case(1)
  select case(m)
  case(:-1)
    polynm=ya(1)*(x-xa(1))
  case(0)
    polynm=ya(1)
  case default
    polynm=0.d0
  end select
  return
case(2)
  c1=(ya(2)-ya(1))/(xa(2)-xa(1))
  t1=x-xa(1)
  select case(m)
  case(:-1)
    polynm=t1*(ya(1)+0.5d0*c1*t1)
  case(0)
    polynm=c1*t1+ya(1)
  case(1)
    polynm=c1
  case default
    polynm=0.d0
  end select
  return
case(3)
  x0=xa(1)
  x1=xa(2)-x0; x2=xa(3)-x0
  y0=ya(1)
  y1=ya(2)-y0; y2=ya(3)-y0
  t0=1.d0/(x1*x2*(x2-x1))
  t1=x1*y2; t2=x2*y1
  c1=x2*t2-x1*t1
  c2=t1-t2
  t1=x-x0
  select case(m)
  case(:-1)
    polynm=t1*(y0+t0*t1*(0.5d0*c1+0.3333333333333333333d0*c2*t1))
  case(0)
    polynm=y0+t0*t1*(c1+c2*t1)
  case(1)
    polynm=t0*(2.d0*c2*t1+c1)
  case(2)
    polynm=t0*2.d0*c2
  case default
    polynm=0.d0
  end select
  return
case(4)
  x0=xa(1)
  x1=xa(2)-x0; x2=xa(3)-x0; x3=xa(4)-x0
  y0=ya(1)
  y1=ya(2)-y0; y2=ya(3)-y0; y3=ya(4)-y0
  t4=x1-x2; t5=x1-x3; t6=x2-x3
  t1=x1*x2*y3; t2=x2*x3*y1; t3=x1*x3
  t0=1.d0/(x2*t3*t4*t5*t6)
  t3=t3*y2
  c3=t1*t4+t2*t6-t3*t5
  t4=x1**2; t5=x2**2; t6=x3**2
  y1=t3*t6-t1*t5; y2=t1*t4-t2*t6; y3=t2*t5-t3*t4
  c2=y1+y2+y3
  c1=x1*y1+x2*y2+x3*y3
  t1=x-x0
  select case(m)
  case(:-1)
    polynm=t1*(y0+t0*t1*(0.5d0*c1+t1*(0.25d0*c3*t1-0.3333333333333333333d0*c2)))
  case(0)
    polynm=y0+t0*t1*(c1+t1*(c3*t1-c2))
  case(1)
    polynm=t0*(c1+t1*(3.d0*c3*t1-2.d0*c2))
  case(2)
    polynm=t0*(6.d0*c3*t1-2.d0*c2)
  case(3)
    polynm=t0*6.d0*c3
  case default
    polynm=0.d0
  end select
  return
end select
if (np.le.0) then
  write(*,*)
  write(*,'("Error(polynm): np <= 0 : ",I8)') np
  write(*,*)
  stop
end if
if (m.ge.np) then
  polynm=0.d0
  return
end if
! find the polynomial coefficients in divided differences form
c(:)=ya(:)
do i=2,np
  do j=np,i,-1
    c(j)=(c(j)-c(j-1))/(xa(j)-xa(j+1-i))
  end do
end do
! special case m=0
if (m.eq.0) then
  sum=c(1)
  t1=1.d0
  do i=2,np
    t1=t1*(x-xa(i-1))
    sum=sum+c(i)*t1
  end do
  polynm=sum
  return
end if
x0=xa(1)
! convert to standard form
do j=1,np-1
  do i=1,np-j
    k=np-i
    c(k)=c(k)+(x0-xa(k-j+1))*c(k+1)
  end do
end do
if (m.gt.0) then
! take the m th derivative
  do j=1,m
    do i=m+1,np
      c(i)=c(i)*dble(i-j)
    end do
  end do
  t1=c(np)
  t2=x-x0
  do i=np-1,m+1,-1
    t1=t1*t2+c(i)
  end do
  polynm=t1
else
! find the integral
  t1=c(np)/dble(np)
  t2=x-x0
  do i=np-1,1,-1
    t1=t1*t2+c(i)/dble(i)
  end do
  polynm=t1*t2
end if
return
end function
!EOC

  !=============================================================!
  subroutine utility_zgemm(a, b, c, transa_opt, transb_opt)
    !=============================================================!
    !                                                             !
    ! Return matrix product of complex matrices a and b:          !
    !                                                             !
    !                       C = Op(A) Op(B)                       !
    !                                                             !
    ! transa = 'N'  ==> Op(A) = A                                 !
    ! transa = 'T'  ==> Op(A) = transpose(A)                      !
    ! transa = 'C'  ==> Op(A) = congj(transpose(A))               !
    !                                                             !
    ! similarly for B                                             !
    !                                                             !
    ! Due to the use of assumed shape arrays, this routine is a   !
    ! safer and more general replacement for the above routine    !
    ! utility_zgemm. Consider removing utility_zgemm and using    !
    ! utility_zgemm_new throughout.                               !
    !                                                             !
    !=============================================================!


    implicit none

    complex(kind=dp), intent(in)            :: a(:, :)
    complex(kind=dp), intent(in)            :: b(:, :)
    complex(kind=dp), intent(out)           :: c(:, :)
    character(len=1), intent(in), optional  :: transa_opt
    character(len=1), intent(in), optional  :: transb_opt

    integer          :: m, n, k
    character(len=1) :: transa, transb

    transa = 'N'
    transb = 'N'
    if (present(transa_opt)) transa = transa_opt
    if (present(transb_opt)) transb = transb_opt

    ! m ... number of rows in Op(A) and C
    ! n ... number of columns in Op(B) and C
    ! k ... number of columns in Op(A) resp. rows in Op(B)
    m = size(c, 1)
    n = size(c, 2)

    if (transa /= 'N') then
      k = size(a, 1)
    else
      k = size(a, 2)
    end if

    call zgemm(transa, transb, m, n, k, cmplx_1, a, size(a, 1), b, size(b, 1), cmplx_0, c, m)

  end subroutine utility_zgemm

  !=============================================================!
  subroutine utility_zgemmm(a, transa, b, transb, c, transc, &
                            prod1, eigval, prod2)
    !===============================================================!
    ! Returns the complex matrix-matrix-matrix product              !
    ! --> prod1 = op(a).op(b).op(c),                                !
    ! where op(a/b/c) are defined according to transa/transb/transc !
    ! (see also documentation of utility_zgemm above)               !
    !                                                               !
    ! If eigval and prod2 are present, also                         !
    ! --> prod2 = op(a).diag(eigval).op(b).op(c)                    !
    ! is returned.                                                  !
    !===============================================================!

    complex(kind=dp), dimension(:, :), intent(in)  :: a, b, c
    character(len=1), intent(in)                  :: transa, transb, transc
    real(kind=dp), dimension(:), optional, &
      intent(in)       :: eigval
    complex(kind=dp), dimension(:, :), optional, &
      intent(out) :: prod1, prod2

    complex(kind=dp), dimension(:, :), allocatable :: tmp
    integer                                       :: nb, mc, i, j

    ! query matrix sizes
    ! naming convention:
    ! matrix op(a) [resp. op(b) and op(c)] is of size na x ma [resp. nb x mb and nc x mc]
    ! only nb (=ma) and mc are explicitly needed
    if (transb /= 'N') then
      nb = size(b, 2)
    else
      nb = size(b, 1)
    end if
    if (transc /= 'N') then
      mc = size(c, 1)
    else
      mc = size(c, 2)
    end if

    ! tmp = op(b).op(c)
    allocate (tmp(nb, mc))
    call utility_zgemm(b, c, tmp, transb, transc)

    ! prod1 = op(a).tmp
    if (present(prod1)) then
      call utility_zgemm(a, tmp, prod1, transa, 'N')
    end if

    if (present(prod2) .and. present(eigval)) then
      ! tmp = diag(eigval).tmp
      forall (i=1:nb, j=1:mc)
      tmp(i, j) = eigval(i)*tmp(i, j)
      end forall
      ! prod2 = op(a).tmp
      call utility_zgemm(a, tmp, prod2, transa, 'N')
    end if
  end subroutine

subroutine utility_zgesvd(arr,u,s,v)
      complex(kind=dp), intent(in)  :: arr(:,:)
      complex(kind=dp), intent(out) :: u(:,:)
      real(kind=dp), intent(out) :: s(:)
      complex(kind=dp), intent(out) :: v(:,:)

      real(kind=dp), allocatable :: rwork(:)
      complex(kind=dp), allocatable :: cwork(:)
     
      integer nsize,info
      nsize=size(arr,1)

      allocate (rwork(5*nsize))
      allocate (cwork(4*nsize))

      call ZGESVD('A', 'A', nsize, nsize, arr, &
                  nsize, s, u, nsize, v, nsize, cwork, &
                  4*nsize, rwork, info)

      deallocate(rwork,cwork)

end subroutine

subroutine utility_zgetri(arr)
      complex(kind=dp), intent(inout)  :: arr(:,:)

      complex(kind=dp), allocatable :: cwork(:)
      integer, allocatable :: ipiv(:)
      integer nsize,info
      nsize=size(arr,1)
      
      allocate (ipiv(nsize))
      allocate (cwork(5*nsize))

      call zgetrf(nsize,nsize,arr,nsize,ipiv,info)
      if (info.eq.0) call zgetri(nsize,arr,nsize,ipiv,cwork,nsize,info)
      if (info.ne.0) call throw("modcom%utility_zgetri","unable to invert a matrix")

      deallocate(cwork,ipiv)

end subroutine
  subroutine unitarize(leftright,m,n,cmat)
 
     integer, intent(in) :: leftright,m,n
     complex(dp), intent(inout) :: cmat(m,n)
     complex(dp) :: U(m,m),V(n,n)
     complex(dp) :: tmp(max(m,n),max(m,n))
     real(dp) :: S(min(n,m))
     integer i,j,k
     ! fix unitarity by SVD. Note that first argument of zgesvd is destroyed
     tmp(1:m,1:n)=cmat
     call utility_zgesvd(tmp(1:m,1:n),U,S,V)
     k=min(m,n)
     do i=1,k
       if (abs(S(i))<epslat) then
         S(i)=0._dp
       else
         S(i)=1._dp/abs(S(i))
       end if
     end do
     ! 1/|A|^1/2
     do i=1,n
       do j=1,n
         if (leftright>0) then
           tmp(i,j)=sum(S(1:k)*conjg(V(1:k,i)*V(1:k,j)))
         else
           tmp(i,j)=sum(S(1:k)*U(i,1:k)*conjg(U(j,1:k)))
         end if
       end do
     end do
     if (leftright>0) then
        cmat(:,:)=matmul(cmat,tmp(1:n,1:n))
     else
        cmat(:,:)=matmul(tmp(1:m,1:m),cmat)
     end if
  end subroutine



end module
