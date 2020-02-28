
module modcom
#ifdef MPI
  use mpi
#endif
implicit none


logical, parameter :: centered_kgrid=.false.
logical, parameter :: centered_qgrid=.false.
integer, parameter :: dp=kind(0.d0)
integer, parameter :: NDIM=3
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
real(dp), parameter :: epslat=1.e-6_dp
real(dp), parameter :: epsengy=1.e-6_dp
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
  SUBROUTINE split_string(instring, string1, string2)
    CHARACTER(256) :: instring
    CHARACTER(256),INTENT(OUT):: string1,string2
    INTEGER :: index

    instring = TRIM(adjustl(instring))

    index = SCAN(instring,' ')
    string1 = instring(1:index-1)
    string2 = instring(index+1:)

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



end module
