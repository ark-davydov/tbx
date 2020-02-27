subroutine arpack_interface(rvec,bmat,which,nev,nx,npack,ia,ja,matrix,eval,evec)
use modcom, only : throw,info
implicit none
logical, intent(in)           :: rvec
character(len=1), intent(in)  :: bmat
character(len=2), intent(in)  :: which
integer, intent(in)           :: nev,nx,npack
integer, intent(in)           :: ia(nx+1)
integer, intent(in)           :: ja(npack)
complex*16, intent(in)        :: matrix(npack)
double precision, intent(out) :: eval(nev)
complex*16, intent(out)       :: evec(nx,nev)
! local
character(len=100) strmsg
Complex*16, parameter  :: sigma=cmplx(0.d0,0.d0,kind=8)
Double precision, parameter     :: eps=1.d-6
integer           maxn, maxnev, maxncv, ldv
integer           iparam(11), ipntr(14)
integer           ido, n, ncv, lworkl, arinfo, ierr, nconv, maxitr, ishfts, mode, counter
integer           iev
Double precision  tol
logical, allocatable          :: select(:)
integer, allocatable          :: ipiv(:)
Complex*16, allocatable       :: D(:), V(:,:), resid(:)
Complex*16, allocatable       :: workd(:), workev(:), workl(:)
Double precision, allocatable :: rwork(:)
#ifdef DEBUG
  integer j
  Complex*16, allocatable  :: ax(:)
  Double precision, allocatable ::  rd(:)
#endif
!!!! common PARDISO and LAPACK vars  !!!
COMPLEX*16, ALLOCATABLE :: BB(:)
INTEGER, PARAMETER :: NRHS=1
!!!! PARDISO VARIABLES (CAPS LOCK) !!!
#ifdef PARDI
integer(8) PT(64)
INTEGER IPARM(64)
INTEGER ERROR
INTEGER, PARAMETER :: MTYPE=-4 ! HERMITIAN INDEFINITE
DOUBLE PRECISION DPARM(64)
INTEGER, PARAMETER :: MAXFCT=1
INTEGER, PARAMETER :: MNUM=1
INTEGER :: PHASE=13
INTEGER, ALLOCATABLE :: PERM(:)
INTEGER, PARAMETER :: MSGLVL=0 ! output messages, 0 otherwise
COMPLEX*16, ALLOCATABLE :: XX(:)
#else
!!!! conventional LA variables !!!
integer lda,ldb,ii,jj,lainfo
integer, allocatable :: laipiv(:)
complex*16, allocatable :: aa(:,:)
#endif

! set dimensions
n      = nx
maxn   = n
ncv    = min(5*nev,maxn)
maxnev = nev
maxncv = ncv
ldv    = maxn
lworkl = 3*maxncv*maxncv+5*maxncv

if (nev+2.gt.ncv) then
   call throw("arpack_interface",&
  "number of requested eigen states should be less or equal the size of the matrix minus 2, which is not the case now")
end if

#ifdef DEBUG
allocate(ax(maxn))
allocate(rd(maxncv))
#endif

allocate(select(maxncv))
allocate(workd(3*maxn))
allocate(resid(maxn))
allocate(V(ldv, maxncv))
allocate(workev(2*maxncv))
allocate(workl(lworkl))
allocate(D(maxncv))
allocate(ipiv(maxn))
allocate(rwork(maxn))
allocate(BB(NX))

tol    = 1.d-8
ido    = 0
arinfo = 0
ishfts = 1
maxitr = 800
mode   = 3
iparam(1) = ishfts 
iparam(3) = maxitr 
iparam(7) = mode 

#ifdef PARDI
  !!!!!!!!! PARDISO !!!!!!!!
  IPARM(3)=1
  IPARM(6)=1
  CALL PARDISOINIT(PT, MTYPE, IPARM)
  allocate(PERM(NX),XX(NX))
#else
  !!!!!!!!! conventional LAPACK  !!!!!!!!!!
  lda=nx
  ldb=nx
  allocate(laipiv(maxn))
  allocate(aa(nx,nx))
  ! unpack the full matrix
  aa=0.d0
  do jj=1,npack
    do ii=1,nx
      if (jj.ge.ia(ii).and.jj.lt.ia(ii+1)) then
         aa(ii,ja(jj))=matrix(jj)
         if (ii.ne.ja(jj)) aa(ja(jj),ii)=conjg(matrix(jj))
      end if
    end do
  end do
  CALL zgetrf( nx, nx, aa, lda, laipiv, lainfo )
#endif

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) | 
!     %-------------------------------------------%
 counter=0
 20   continue
         counter=counter+1

         call znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,&
              v, ldv, iparam, ipntr, workd, workl, lworkl,&
              rwork,arinfo )
         if (ido .eq. -1 .or. ido .eq. 1 ) then
            call zcopy ( NX, workd(ipntr(1)), 1, BB, 1)

#ifdef PARDI
            ! PARDISO
            ! Solve the system A*X = B with PARDISO
            if (counter.eq.1) then
              PT=0
              PHASE=13
            else
              PHASE=33
            end if
            call PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, NX, MATRIX, IA, JA,&
                   PERM, NRHS, IPARM, MSGLVL, BB, XX, ERROR, DPARM)
            ! check for errors
            if ( ERROR .ne. 0 ) then
               write(strmsg,*) ERROR
               call throw("arpack_interface%PARDISO()",strmsg)
            end if
            call zcopy ( NX, XX(1), 1, workd(ipntr(2)), 1)
#else
            ! CONVENTIONAL LAPACK
            call zcopy ( NX, workd(ipntr(1)), 1, BB, 1)
            ! Solve the system A*X = B with conventional LAPACK
            CALL zgetrs( 'N', nx, NRHS, aa, lda, laipiv, bb, ldb, lainfo )
            ! check for errors
            if ( lainfo .ne. 0 ) then
               write(strmsg,*) lainfo
               call throw("arpack_interface%zgetrs()",strmsg)
            end if
            call zcopy ( NX, BB, 1, workd(ipntr(2)), 1)
#endif
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
            go to 20
         end if
!     %-----------------------------------------%
!     | Either we have convergence, or there is |
!     | an error.                               |
!     %-----------------------------------------%
      if ( arinfo .lt. 0 ) then
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
         write(strmsg,*) arinfo
         call throw("arpack_interface%znaupd()",strmsg)
      else 
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD .                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |  
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    | 
!        %-------------------------------------------%
 !        select=.true.
         call zneupd  (rvec, 'A', select, D, V, ldv, sigma, &
                      workev, bmat, n, which, nev, tol,&
                      resid, ncv, v, ldv, iparam, ipntr, workd,&
                      workl, lworkl, rwork, ierr)

!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
         if ( ierr .ne. 0) then
!         %------------------------------------%
!         | Error condition:                   | 
!         | Check the documentation of ZNEUPD . |
!         %------------------------------------%
            write(*,*) "ZNEUPD gave an error: ", ierr
         else
             nconv = iparam(5) 
#ifdef DEBUG
             do 60 j=1, nconv
!             %---------------------------%
!             | Compute the residual norm |
!             |                           |
!             |   ||  A*x - lambda*x ||   |
!             |                           |
!             | for the NCONV accurately  |
!             | computed eigenvalues and  |
!             | eigenvectors.  (iparam(5) |
!             | indicates how many are    |
!             | accurate to the requested |
!             | tolerance)                |
!             %---------------------------%
              call av(nx, v(:,j), ax)
              call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
              rd(j) = dble (d(j))
 60         continue
!         %-----------------------------%
!         | Display computed residuals. |
!         %-----------------------------%
             call dmout (6, nconv, 1, rd, maxncv, -6, &
                  'Ritz values (Real, Imag) and relative residuals')
#endif
         end if
!     %-------------------------------------------%
!     | Print additional convergence information. |
!     %-------------------------------------------%
         if ( arinfo .eq. 1) then
             call info("arpack_interface","Maximum number of iterations reached")
         else if ( arinfo .eq. 3) then
             call throw("arpack_interface",&
             "No shifts could be applied during implicit Arnoldi update, try increasing NCV")
         end if      

#ifdef DEBUG
         print *, ' '
         print *, '_NDRV2 '
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated', &
                  ' (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', &
                    nconv 
         print *, ' The number of Implicit Arnoldi update', &
                  ' iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '
#endif
      end if
      if (nev.ne.nconv) then
         call info("arpack_interface","Warning, NEV not equal to NCONV")
      end if
      do iev=1,nev
        eval(iev)=dble(D(iev))
        if (abs(eval(iev)-D(iev)).gt.1.d0) then
           write(strmsg,*) "eigenvalues have large imaginary part ",abs(eval(iev)-D(iev))
           call throw("arpack_interface",strmsg)
        end if
        if (rvec) evec(1:nx,iev)=V(1:nx,iev)
      end do

! ARPACK arrays
#ifdef DEBUG
deallocate(ax)
deallocate(rd)
#endif
deallocate(select)
deallocate(workd)
deallocate(resid)
deallocate(V)
deallocate(workev)
deallocate(workl)
deallocate(D)
deallocate(ipiv)
deallocate(rwork)
#ifdef PARDI
  ! deallocate all memory with phase -1
  PHASE=-1
  call PARDISO(PT, MAXFCT, MNUM, MTYPE, PHASE, NX, MATRIX, IA, JA,&
         PERM, NRHS, IPARM, MSGLVL, BB, XX, ERROR, DPARM)
  deallocate(PERM,XX)
#else
  deallocate(laipiv)
  deallocate(aa)
#endif
deallocate(BB)

contains

#ifdef DEBUG
      subroutine av (nx, vv, ww)
      integer, intent(in) :: nx
      Complex*16, intent(inout)    ::    vv(:), ww(:)
      integer ii,jj
      ww(1:nx)=0.d0
      do jj=1,npack
        do ii=1,nx
          if (jj.ge.ia(ii).and.jj.lt.ia(ii+1)) then
             ! ij matrix compontent
             ww(ii)=ww(ii)+matrix(jj)*vv(ja(jj))
             if (ii.ne.ja(jj)) then
               ww(ja(jj))=ww(ja(jj))+conjg(matrix(jj))*vv(ii)
             end if
          end if
        end do
      end do

      return
      end subroutine
#endif

end subroutine
