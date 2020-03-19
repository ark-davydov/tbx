module wannier_supplementary
use modcom, only : throw,NDIM,twopi
use parameters, only : CLpars
use symmetryclass
implicit none
private
integer, parameter :: maxallowd_orbs=46340
integer, parameter :: stdout=6
integer, parameter :: dp=kind(0.d0)
real(dp), parameter :: eps6=1.e-6_dp
real(dp), parameter :: eps8=1.e-8_dp
real(dp), parameter :: eps14=1.e-14_dp
real(dp), parameter :: pi=3.141592653589793115997963468544185161590576171875_dp
real(dp), parameter :: tpi=6.283185307179586231995926937088370323181152343750_dp
real(dp), parameter :: fpi=12.5663706143591724639918538741767406463623046875_dp


public :: raddec,ylm_wannier,set_u_matrix

type, public :: wbase
   integer :: norb
   integer :: ncenters
   integer :: nmaxorb_pcent
   integer, allocatable :: lmr(:,:)
   integer, allocatable :: norb_ic(:)
   integer, allocatable :: orb_icio(:,:)
   integer, allocatable :: icio_orb(:,:)
   integer, allocatable :: ics2c(:,:)
   real(dp), allocatable :: vcs2t(:,:,:)
   real(dp), allocatable :: waxis(:,:,:)
   real(dp), allocatable :: centers(:,:)
   real(dp), allocatable :: centers_cart(:,:)
   contains 
   procedure :: init
   procedure :: init_smap
   procedure :: wws
   procedure, nopass :: wws_full
endtype

contains

subroutine init(THIS,pars,ncenters,norb,norb_ic,lmr,waxis,centers)
! allocate a basis for TB model
class(wbase), intent(inout) :: THIS
class(CLpars), intent(in) :: pars
integer, intent(in) :: ncenters,norb,norb_ic(ncenters)
integer, optional, intent(in) :: lmr(2,norb)
real(dp), optional, intent(in) :: waxis(NDIM,2,norb)
real(dp), optional, intent(in) :: centers(NDIM,ncenters)
! local
integer iorb,ic,ios
THIS%ncenters=ncenters
allocate(THIS%norb_ic(ncenters))
THIS%norb_ic=norb_ic
THIS%nmaxorb_pcent=maxval(THIS%norb_ic(:))
allocate(THIS%icio_orb(THIS%ncenters,THIS%nmaxorb_pcent))
! compute the total number of orbitals
THIS%norb=0
do ic=1,THIS%ncenters
  do ios=1,THIS%norb_ic(ic)
     THIS%norb=THIS%norb+1
     if(THIS%norb.gt.maxallowd_orbs) then
       call throw("wbase%init()","maximum allowed basis orbitals exceeded (integer overflow in NxN value)")
     end if
     ! map: center,orbialindex to the final basis orbital
     THIS%icio_orb(ic,ios)=THIS%norb
  end do
end do
allocate(THIS%orb_icio(THIS%norb,2))
do ic=1,THIS%ncenters
  do ios=1,THIS%norb_ic(ic)
     iorb=THIS%icio_orb(ic,ios)
     ! map:  final basis orbital to center,orbialindex
     THIS%orb_icio(iorb,1)=ic
     THIS%orb_icio(iorb,2)=ios
  end do
end do
allocate(THIS%lmr(2,THIS%norb))
allocate(THIS%waxis(NDIM,2,THIS%norb))
allocate(THIS%centers(NDIM,THIS%ncenters))
allocate(THIS%centers_cart(NDIM,THIS%ncenters))
if (THIS%norb.ne.norb) call throw("wbase%init","input orbital list has unexpected size")
do iorb=1,THIS%norb
  THIS%lmr(:,iorb)=lmr(:,iorb)
  THIS%waxis(:,:,iorb)=waxis(:,:,iorb)
end do
do ic=1,THIS%ncenters
  THIS%centers(:,ic)=centers(:,ic)
  THIS%centers_cart(:,ic)=matmul(centers(:,ic),pars%avec)
end do
end subroutine

subroutine init_smap(THIS,sym,pars)
class(wbase), intent(inout) :: THIS
class(CLsym), intent(in) :: sym
class(CLpars), intent(in) :: pars
integer isym,ic,jc
logical lfound(THIS%ncenters)
real(dp) v1(NDIM),v2(NDIM),v3(NDIM),v4(NDIM)
allocate(THIS%ics2c(THIS%ncenters,sym%nsym))
THIS%ics2c=-999
!write(*,"(a,i5)") "  Number of symmetry operators = ", sym%nsym
do isym=1,sym%nsym
!  write(*,"(2x,i5,a)") isym, "-th symmetry operators is"
  !Writing rotation matrix and translation vector in Cartesian (QE units) coordinates
!  write(*,"(3f15.7)") &
!     sym%car(:,:,isym), sym%vtc(:,isym)/sqrt(sum(pars%avec(1,:)**2))
  lfound=.false.
  do ic=1,THIS%ncenters
    v1=THIS%centers_cart(:,ic)
    v2=matmul(sym%car(:,:,isym),v1+sym%vtc(:,isym))
 !   write(*,*) matmul(sym%lat(:,:,isym),THIS%centers(:,ic)+sym%vtl(:,isym))
    do jc=1,THIS%ncenters
      if(lfound(jc)) cycle
      v3=THIS%centers_cart(:,jc)
      v4=matmul(pars%bvec,v3-v2)/twopi
      if(sum(abs(dble(nint(v4))-v4)).lt.1d-2) then
         lfound(jc)=.true.
         THIS%ics2c(ic,isym)=jc
         exit !Sym.op.(isym) moves position(ips2p(ip,isym)) to position(ip) + T, where
      end if                                       !T is given by vps2t(:,ip,isym).
    end do
    if(THIS%ics2c(ic,isym).le.0) then
       write(*,"(a,3f18.10)")"  Could not find ",matmul(pars%bvec,v2)/twopi
       write(*,"(a,3f18.10)")"  coming from    ",matmul(pars%bvec,v1)/twopi
       write(*,"(a,i5,a               )")"  of Wannier site",ic,"."
       write(*,"(a,i5,a               )")"  symmetry ",isym,"."
       write(*,"(3F10.6)") sym%car(1,:,isym)
       write(*,"(3F10.6)") sym%car(2,:,isym)
       write(*,"(3F10.6)") sym%car(3,:,isym)
       call throw("wannier_supplemetart%init_smap", "missing Wannier sites")
    end if
  end do
end do
allocate(THIS%vcs2t(NDIM,THIS%ncenters,sym%nsym))
THIS%vcs2t=0._dp
do isym=1,sym%nsym
  do ic=1,THIS%ncenters
    v1=THIS%centers_cart(:,ic)
    jc=THIS%ics2c(ic,isym)
    v2=THIS%centers_cart(:,jc)
    v3=matmul(v2,sym%car(:,:,isym))-sym%vtc(:,isym)
    THIS%vcs2t(:,ic,isym)=v3-v1
  end do
end do
end subroutine

real(dp) function wws(THIS,sr,iorb,jorb)
class(wbase), intent(in) :: THIS
real(dp), intent(in) :: sr(3,3)
integer, intent(in) :: iorb,jorb
integer l_w1,mr_w1,l_w2,mr_w2
real(dp) xaxis1(3),zaxis1(3),xaxis2(3),zaxis2(3)
! assign l,mr,xaxis,zaxis from orbital indices
l_w1=THIS%lmr(1,iorb)
l_w2=THIS%lmr(1,jorb)
mr_w1=THIS%lmr(2,iorb)
mr_w2=THIS%lmr(2,jorb)
xaxis1=THIS%waxis(:,1,iorb)
xaxis2=THIS%waxis(:,1,jorb)
zaxis1=THIS%waxis(:,2,iorb)
zaxis2=THIS%waxis(:,2,jorb)
wws=wws_full(sr,l_w1,mr_w1,l_w2,mr_w2,xaxis1,zaxis1,xaxis2,zaxis2)
end function


real(dp) function wws_full(sr,l_w1,mr_w1,l_w2,mr_w2,xaxis1,zaxis1,xaxis2,zaxis2)
real(dp), intent(in) :: sr(3,3)
integer, intent(in) ::  l_w1,mr_w1,l_w2,mr_w2
real(dp), intent(in) :: xaxis1(3),zaxis1(3),xaxis2(3),zaxis2(3)
! local
real(dp), parameter :: pwg(2)=(/2.976190476190479d-2,3.214285714285711d-2/)
integer ip,ir
real(DP) dvec(3,32),dwgt(32),dvec2(3,32),dylm(32,5),vaxis1(3,3),vaxis2(3,3)
real(dp), parameter :: p12(3,12)=reshape(                            &
   (/0d0, 0d0, 1.00000000000000d0,                                   &
     0.894427190999916d0, 0d0, 0.447213595499958d0,                  &
     0.276393202250021d0, 0.850650808352040d0, 0.447213595499958d0,  &
    -0.723606797749979d0, 0.525731112119134d0, 0.447213595499958d0,  &
    -0.723606797749979d0, -0.525731112119134d0, 0.447213595499958d0, &
     0.276393202250021d0, -0.850650808352040d0, 0.447213595499958d0, &
     0.723606797749979d0, 0.525731112119134d0, -0.447213595499958d0, &
    -0.276393202250021d0, 0.850650808352040d0, -0.447213595499958d0, &
    -0.894427190999916d0, 0d0, -0.447213595499958d0,                 &
    -0.276393202250021d0, -0.850650808352040d0, -0.447213595499958d0,&
     0.723606797749979d0, -0.525731112119134d0, -0.447213595499958d0,&
     0d0, 0d0, -1.00000000000000d0/),(/3,12/))
real(dp), parameter :: p20(3,20)=reshape(                            &
   (/0.525731112119134d0, 0.381966011250105d0, 0.850650808352040d0,  &
    -0.200811415886227d0, 0.618033988749895d0, 0.850650808352040d0,  &
    -0.649839392465813d0, 0d0, 0.850650808352040d0,                  &
    -0.200811415886227d0, -0.618033988749895d0, 0.850650808352040d0, &
     0.525731112119134d0, -0.381966011250105d0, 0.850650808352040d0, &
     0.850650808352040d0, 0.618033988749895d0, 0.200811415886227d0,  &
    -0.324919696232906d0, 1.00000000000000d0, 0.200811415886227d0,   &
    -1.05146222423827d0, 0d0, 0.200811415886227d0,                   &
   -0.324919696232906d0, -1.00000000000000d0, 0.200811415886227d0,   &
    0.850650808352040d0, -0.618033988749895d0, 0.200811415886227d0,  &
    0.324919696232906d0, 1.00000000000000d0, -0.200811415886227d0,   &
   -0.850650808352040d0, 0.618033988749895d0, -0.200811415886227d0,  &
   -0.850650808352040d0, -0.618033988749895d0, -0.200811415886227d0, &
    0.324919696232906d0, -1.00000000000000d0, -0.200811415886227d0,  &
    1.05146222423827d0, 0d0, -0.200811415886227d0,                   &
    0.200811415886227d0, 0.618033988749895d0, -0.850650808352040d0,  &
   -0.525731112119134d0, 0.381966011250105d0, -0.850650808352040d0,  &
   -0.525731112119134d0, -0.381966011250105d0, -0.850650808352040d0, &
    0.200811415886227d0, -0.618033988749895d0, -0.850650808352040d0, &
   0.649839392465813d0, 0d0, -0.850650808352040d0/),(/3,20/))


! code
dvec(:,1:12)=p12
dvec(:,13:32)=p20
do ip=1,32
   dvec(:,ip)=dvec(:,ip)/sqrt(sum(dvec(:,ip)**2))
end do
dwgt(1:12)=pwg(1)
dwgt(13:32)=pwg(2)
!write(stdout,*) sum(dwgt) !Checking the weight sum to be 1.
dylm=0d0
vaxis1=0d0
vaxis2=0d0
do ip=1,5
   CALL ylm_wannier(dylm(1,ip),2,ip,dvec,32)
end do
!do ip=1,5
!   write(stdout,"(5f25.15)") (sum(dylm(:,ip)*dylm(:,jp)*dwgt)*2d0*tpi,jp=1,5)
!end do !Checking spherical integral.
call set_u_matrix(xaxis1,zaxis1,vaxis1)
call set_u_matrix(xaxis2,zaxis2,vaxis2)
CALL ylm_wannier(dylm(1,1),l_w1,mr_w1,matmul(vaxis1,dvec),32)
do ir=1,32
   dvec2(:,ir)=matmul(sr,dvec(:,ir))
end do
CALL ylm_wannier(dylm(1,2),l_w2,mr_w2,matmul(vaxis2,dvec2),32)
wws_full=sum(dylm(:,1)*dylm(:,2)*dwgt)*2d0*tpi !<Rotated Y(jw)|Not rotated Y(iw)> for sym.op.(isym).
end function




SUBROUTINE set_u_matrix(x,z,u)
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   IF (xx < eps6) CALL throw ('set_u_matrix',' |xaxis| < eps ')
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   IF (zz < eps6) CALL throw ('set_u_matrix',' |zaxis| < eps ')
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   IF (abs(coseno) > eps6) CALL throw('set_u_matrix',' xaxis and zaxis are not orthogonal !')

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)

   RETURN

END SUBROUTINE set_u_matrix

SUBROUTINE ylm_wannier(ylm,l,mr,r,nr)
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
! of the spherical harmonic identified  by indices (l,mr)
! in table 3.1 of the wannierf90 specification.
!
! No reference to the particular ylm ordering internal to Quantum ESPRESSO
! is assumed.
!
! If ordering in wannier90 code is changed or extended this should be the
! only place to be modified accordingly
!

! I/O variables
!
   INTEGER :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
 !  real(DP) :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy
 !  real(DP) :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   INTEGER :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   IF (l > 3 .or. l < -5 ) CALL throw('ylm_wannier',' l out of range ')
   IF (l>=0) THEN
      IF (mr < 1 .or. mr > 2*l+1) CALL throw('ylm_wannier','mr out of range' )
   ELSE
      IF (mr < 1 .or. mr > abs(l)+1 ) CALL throw('ylm_wannier','mr out of range')
   ENDIF

   DO ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      IF (rr < eps8) CALL throw('ylm_wannier',' rr too small ')

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1,ir) > eps8) THEN
         phi = atan( r(2,ir)/r(1,ir) )
      ELSEIF (r(1,ir) < -eps8 ) THEN
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      ELSE
         phi = sign( pi/2.d0,r(2,ir) )
      ENDIF


      IF (l==0) THEN   ! s orbital
                    ylm(ir) = s()
      ENDIF
      IF (l==1) THEN   ! p orbitals
         IF (mr==1) ylm(ir) = p_z(cost)
         IF (mr==2) ylm(ir) = px(cost,phi)
         IF (mr==3) ylm(ir) = py(cost,phi)
      ENDIF
      IF (l==2) THEN   ! d orbitals
         IF (mr==1) ylm(ir) = dz2(cost)
         IF (mr==2) ylm(ir) = dxz(cost,phi)
         IF (mr==3) ylm(ir) = dyz(cost,phi)
         IF (mr==4) ylm(ir) = dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = dxy(cost,phi)
      ENDIF
      IF (l==3) THEN   ! f orbitals
         IF (mr==1) ylm(ir) = fz3(cost)
         IF (mr==2) ylm(ir) = fxz2(cost,phi)
         IF (mr==3) ylm(ir) = fyz2(cost,phi)
         IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
         IF (mr==5) ylm(ir) = fxyz(cost,phi)
         IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      ENDIF
      IF (l==-1) THEN  !  sp hybrids
         IF (mr==1) ylm(ir) = bs2 * ( s() + px(cost,phi) )
         IF (mr==2) ylm(ir) = bs2 * ( s() - px(cost,phi) )
      ENDIF
      IF (l==-2) THEN  !  sp2 hybrids
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi)
      ENDIF
      IF (l==-3) THEN  !  sp3 hybrids
         IF (mr==1) ylm(ir) = 0.5d0*(s()+px(cost,phi)+py(cost,phi)+p_z(cost))
         IF (mr==2) ylm(ir) = 0.5d0*(s()+px(cost,phi)-py(cost,phi)-p_z(cost))
         IF (mr==3) ylm(ir) = 0.5d0*(s()-px(cost,phi)+py(cost,phi)-p_z(cost))
         IF (mr==4) ylm(ir) = 0.5d0*(s()-px(cost,phi)-py(cost,phi)+p_z(cost))
      ENDIF
      IF (l==-4) THEN  !  sp3d hybrids
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi)
         IF (mr==4) ylm(ir) = bs2*p_z(cost)+bs2*dz2(cost)
         IF (mr==5) ylm(ir) =-bs2*p_z(cost)+bs2*dz2(cost)
      ENDIF
      IF (l==-5) THEN  ! sp3d2 hybrids
         IF (mr==1) ylm(ir) = bs6*s()-bs2*px(cost,phi)-bs12*dz2(cost)+.5d0*dx2my2(cost,phi)
         IF (mr==2) ylm(ir) = bs6*s()+bs2*px(cost,phi)-bs12*dz2(cost)+.5d0*dx2my2(cost,phi)
         IF (mr==3) ylm(ir) = bs6*s()-bs2*py(cost,phi)-bs12*dz2(cost)-.5d0*dx2my2(cost,phi)
         IF (mr==4) ylm(ir) = bs6*s()+bs2*py(cost,phi)-bs12*dz2(cost)-.5d0*dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = bs6*s()-bs2*p_z(cost)+bs3*dz2(cost)
         IF (mr==6) ylm(ir) = bs6*s()+bs2*p_z(cost)+bs3*dz2(cost)
      ENDIF

   ENDDO

   RETURN

END SUBROUTINE ylm_wannier

!======== l = 0 =====================================================================
real(DP) FUNCTION s()
   s = 1.d0/ sqrt(fpi)
END FUNCTION s
!======== l = 1 =====================================================================
FUNCTION p_z(cost)
   real(DP) ::p_z, cost
   p_z =  sqrt(3.d0/fpi) * cost
END FUNCTION p_z
FUNCTION px(cost,phi)
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)
END FUNCTION px
FUNCTION py(cost,phi)
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)
END FUNCTION py
!======== l = 2 =====================================================================
FUNCTION dz2(cost)
   real(DP) ::dz2, cost
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)
END FUNCTION dz2
FUNCTION dxz(cost,phi)
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)
END FUNCTION dxz
FUNCTION dyz(cost,phi)
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)
END FUNCTION dyz
FUNCTION dx2my2(cost,phi)
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)
END FUNCTION dx2my2
FUNCTION dxy(cost,phi)
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)
END FUNCTION dxy
!======== l = 3 =====================================================================
FUNCTION fz3(cost)
   real(DP) ::fz3, cost
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost
END FUNCTION fz3
FUNCTION fxz2(cost,phi)
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
END FUNCTION fxz2
FUNCTION fyz2(cost,phi)
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
END FUNCTION fyz2
FUNCTION fzx2my2(cost,phi)
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)
END FUNCTION fzx2my2
FUNCTION fxyz(cost,phi)
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)
END FUNCTION fxyz
FUNCTION fxx2m3y2(cost,phi)
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)
END FUNCTION fxx2m3y2
FUNCTION fy3x2my2(cost,phi)
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)
END FUNCTION fy3x2my2

!
!
!-----------------------------------------------------------------------
!SUBROUTINE radialpart(omega,ng, q, alfa, rvalue, lmax, radial)
!  !-----------------------------------------------------------------------
!  !
!  ! This routine computes a table with the radial Fourier transform
!  ! of the radial functions.
!  !
!  !
!  ! I/O
!  INTEGER :: ng, rvalue, lmax
!  real(DP) :: omega, q(ng), alfa, radial(ng,0:lmax)
!  ! local variables
!  real(DP), PARAMETER :: xmin=-6.d0, dx=0.025d0, rmax=10.d0
!
!  real(DP) :: rad_int, pref, x
!  INTEGER :: l, ir, ig, mesh_r
!  real(DP), ALLOCATABLE :: bes(:), func_r(:), r(:), rij(:), aux(:)
!
!  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
!  ALLOCATE ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
!  ALLOCATE ( aux(mesh_r))
!  !
!  !    compute the radial mesh
!  !
!  DO ir = 1, mesh_r
!     x = xmin  + dble (ir - 1) * dx
!     r (ir) = exp (x) / alfa
!     rij (ir) = dx  * r (ir)
!  ENDDO
!  !
!  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
!  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * &
!                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
!  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(3.0d0/2.0d0) * &
!                     (1.d0 - 2.0d0/3.0d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
!                                           exp(-alfa*r(:)/3.0d0)
!  pref = fpi/sqrt(omega)
!  !
!  DO l = 0, lmax
!     DO ig=1,ng
!       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
!       aux(:) = bes(:) * func_r(:) * r(:) * r(:)
!       ! second r factor added upo suggestion by YY Liang
!       CALL simpson (mesh_r, aux, rij, rad_int)
!       radial(ig,l) = rad_int * pref
!     ENDDO
!  ENDDO
!
!  DEALLOCATE (bes, func_r, r, rij, aux )
!  RETURN
!END SUBROUTINE radialpart

!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
!subroutine sph_bes (msh, r, q, l, jl)
!  !--------------------------------------------------------------------
!  !
!  ! ... input:
!  ! ...   msh     = number of grid points points
!  ! ...   r(1:msh)= radial grid
!  ! ...   q       = q
!  ! ...   l       = angular momentum (-1 <= l <= 6)
!  ! ... output:
!  ! ...   jl(1:msh) = j_l(q*r(i))  (j_l = spherical bessel function)
!  !
!  !
!  !
!  integer :: msh, l
!  real(DP) :: r (msh), q, jl (msh)
!  !
!  ! xseries = convergence radius of the series for small x of j_l(x)
!  real(DP) :: x, xl, xseries = 0.05_dp
!  integer :: ir, ir0
!  !
!#if defined (__MASS)
!  real(DP) :: qr(msh), sin_qr(msh), cos_qr(msh)
!#endif
! 
!  !  case q=0
!
!  if (abs (q) < eps14) then
!     if (l == -1) then
!        call throw ('sph_bes', 'j_{-1}(0) ?!?')
!     elseif (l == 0) then
!        jl(:) = 1.d0
!     else
!        jl(:) = 0.d0
!     endif
!     return
!  end if 
!
!  !  case l=-1
!
!  if (l == - 1) then
!     if (abs (q * r (1) ) < eps14) call throw ('sph_bes', 'j_{-1}(0) ?!?')
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     jl = cos_qr / qr
!
!#else
!
!     jl (:) = cos (q * r (:) ) / (q * r (:) )
!
!#endif
!
!     return
!
!  end if
!
!  ! series expansion for small values of the argument
!  ! ir0 is the first grid point for which q*r(ir0) > xseries
!  ! notice that for small q it may happen that q*r(msh) < xseries !
!
!  ir0 = msh+1
!  do ir = 1, msh
!     if ( abs (q * r (ir) ) > xseries ) then
!        ir0 = ir
!        exit
!     end if
!  end do
!
!  do ir = 1, ir0 - 1
!     x = q * r (ir)
!     if ( l == 0 ) then
!        xl = 1.0_dp
!     else
!        xl = x**l
!     end if
!     jl (ir) = xl/semifact(2*l+1) * &
!                ( 1.0_dp - x**2/1.0_dp/2.0_dp/(2.0_dp*l+3) * &
!                ( 1.0_dp - x**2/2.0_dp/2.0_dp/(2.0_dp*l+5) * &
!                ( 1.0_dp - x**2/3.0_dp/2.0_dp/(2.0_dp*l+7) * &
!                ( 1.0_dp - x**2/4.0_dp/2.0_dp/(2.0_dp*l+9) ) ) ) )
!  end do
!
!  ! the following shouldn't be needed but do you trust compilers
!  ! to do the right thing in this special case ? I don't - PG
!
!  if ( ir0 > msh ) return
!
!  if (l == 0) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = sin_qr(ir0:) / (q * r (ir0:) )
!
!#else
!
!     jl (ir0:) = sin (q * r (ir0:) ) / (q * r (ir0:) )
!
!#endif
!
!  elseif (l == 1) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = ( sin_qr(ir0:) / (q * r (ir0:) ) - &
!                   cos_qr(ir0:) ) / (q * r (ir0:) )
!
!#else
!
!     jl (ir0:) = (sin (q * r (ir0:) ) / (q * r (ir0:) ) - &
!                  cos (q * r (ir0:) ) ) / (q * r (ir0:) )
!
!#endif
!
!  elseif (l == 2) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = ( (3.d0 / (q*r(ir0:)) - (q*r(ir0:)) ) * sin_qr(ir0: ) - &
!                    3.d0 * cos_qr(ir0:) ) / (q*r(ir0:))**2
!
!#else
!
!     jl (ir0:) = ( (3.d0 / (q*r(ir0:)) - (q*r(ir0:)) ) * sin (q*r(ir0:)) - &
!                    3.d0 * cos (q*r(ir0:)) ) / (q*r(ir0:))**2
!
!#endif
!
!  elseif (l == 3) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = (sin_qr (ir0:) * &
!                  (15.d0 / (q*r(ir0:)) - 6.d0 * (q*r(ir0:)) ) + &
!                  cos_qr (ir0:) * ( (q*r(ir0:))**2 - 15.d0) ) / &
!                  (q*r(ir0:))**3
!
!#else
!
!     jl (ir0:) = (sin (q*r(ir0:)) * &
!                  (15.d0 / (q*r(ir0:)) - 6.d0 * (q*r(ir0:)) ) + &
!                  cos (q*r(ir0:)) * ( (q*r(ir0:))**2 - 15.d0) ) / &
!                  (q*r(ir0:)) **3
!
!#endif
!
!  elseif (l == 4) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = (sin_qr (ir0:) * &
!                  (105.d0 - 45.d0 * (q*r(ir0:))**2 + (q*r(ir0:))**4) + &
!                  cos_qr (ir0:) * &
!                  (10.d0 * (q*r(ir0:))**3 - 105.d0 * (q*r(ir0:))) ) / &
!                    (q*r(ir0:))**5
!
!#else
!
!     jl (ir0:) = (sin (q*r(ir0:)) * &
!                  (105.d0 - 45.d0 * (q*r(ir0:))**2 + (q*r(ir0:))**4) + &
!                  cos (q*r(ir0:)) * &
!                  (10.d0 * (q*r(ir0:))**3 - 105.d0 * (q*r(ir0:))) ) / &
!                     (q*r(ir0:))**5
!#endif
!
!  elseif (l == 5) then
!
!#if defined (__MASS)
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = (-cos_qr(ir0:) - &
!                  (945.d0*cos_qr(ir0:)) / (q*r(ir0:)) ** 4 + &
!                  (105.d0*cos_qr(ir0:)) / (q*r(ir0:)) ** 2 + &
!                  (945.d0*sin_qr(ir0:)) / (q*r(ir0:)) ** 5 - &
!                  (420.d0*sin_qr(ir0:)) / (q*r(ir0:)) ** 3 + &
!                  ( 15.d0*sin_qr(ir0:)) / (q*r(ir0:)) ) / (q*r(ir0:))
!#else
!     jl (ir0:) = (-cos(q*r(ir0:)) - &
!                  (945.d0*cos(q*r(ir0:))) / (q*r(ir0:)) ** 4 + &
!                  (105.d0*cos(q*r(ir0:))) / (q*r(ir0:)) ** 2 + &
!                  (945.d0*sin(q*r(ir0:))) / (q*r(ir0:)) ** 5 - &
!                  (420.d0*sin(q*r(ir0:))) / (q*r(ir0:)) ** 3 + &
!                  ( 15.d0*sin(q*r(ir0:))) / (q*r(ir0:)) ) / (q*r(ir0:))
!#endif
!
!  elseif (l == 6) then
!
!#if defined (__MASS)
!
!     qr = q * r
!     call vcos( cos_qr, qr, msh)
!     call vsin( sin_qr, qr, msh)
!     jl (ir0:) = ((-10395.d0*cos_qr(ir0:)) / (q*r(ir0:))**5 + &
!                  (  1260.d0*cos_qr(ir0:)) / (q*r(ir0:))**3 - &
!                  (    21.d0*cos_qr(ir0:)) / (q*r(ir0:))    - &
!                             sin_qr(ir0:)                   + &
!                  ( 10395.d0*sin_qr(ir0:)) / (q*r(ir0:))**6 - &
!                  (  4725.d0*sin_qr(ir0:)) / (q*r(ir0:))**4 + &
!                  (   210.d0*sin_qr(ir0:)) / (q*r(ir0:))**2 ) / (q*r(ir0:))
!#else
!
!     jl (ir0:) = ((-10395.d0*cos(q*r(ir0:))) / (q*r(ir0:))**5 + &
!                  (  1260.d0*cos(q*r(ir0:))) / (q*r(ir0:))**3 - &
!                  (    21.d0*cos(q*r(ir0:))) / (q*r(ir0:))    - &
!                             sin(q*r(ir0:))                   + &
!                  ( 10395.d0*sin(q*r(ir0:))) / (q*r(ir0:))**6 - &
!                  (  4725.d0*sin(q*r(ir0:))) / (q*r(ir0:))**4 + &
!                  (   210.d0*sin(q*r(ir0:))) / (q*r(ir0:))**2 ) / (q*r(ir0:))
!#endif
!
!  else
!
!     call throw ('sph_bes', 'not implemented')
!
!  endif
!  !
!  return
!end subroutine sph_bes
!
integer function semifact(n)
  ! semifact(n) = n!!
  integer :: n, i

  semifact = 1
  do i = n, 1, -2
     semifact = i*semifact
  end do
  return
end function semifact

!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!SUBROUTINE simpson(mesh, func, rab, asum)
!  !-----------------------------------------------------------------------
!  !
!  !     simpson's rule integration. On input:
!  !       mesh = the number of grid points (should be odd)
!  !       func(i)= function to be integrated
!  !       rab(i) = r(i) * dr(i)/di * di
!  !     For the logarithmic grid not including r=0 :
!  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
!  !     For the logarithmic grid including r=0 :
!  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
!  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
!  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
!  !
!  INTEGER, INTENT(in) :: mesh
!  real(DP), INTENT(in) :: rab (mesh), func (mesh)
!  real(DP), INTENT(out):: asum
!  !
!  real(DP) :: f1, f2, f3, r12
!  INTEGER :: i
!  !
!  asum = 0.0d0
!  r12 = 1.0d0 / 3.0d0
!  f3 = func (1) * rab (1) * r12
!
!  DO i = 2, mesh - 1, 2
!     f1 = f3
!     f2 = func (i) * rab (i) * r12
!     f3 = func (i + 1) * rab (i + 1) * r12
!     asum = asum + f1 + 4.0d0 * f2 + f3
!  ENDDO
!  !
!  ! if mesh is not odd, use open formula instead:
!  ! ... 2/3*f(n-5) + 4/3*f(n-4) + 13/12*f(n-3) + 0*f(n-2) + 27/12*f(n-1)
!  !!! Under testing
!  !
!  !IF ( MOD(mesh,2) == 0 ) THEN
!  !   print *, 'mesh even: correction:', f1*5.d0/4.d0-4.d0*f2+23.d0*f3/4.d0, &
!  !                                      func(mesh)*rab(mesh), asum
!  !   asum = asum + f1*5.d0/4.d0 - 4.d0*f2 + 23.d0*f3/4.d0
!  !END IF
!
!  RETURN
!END SUBROUTINE simpson

!=-----------------------------------------------------------------------
!SUBROUTINE simpson_cp90( mesh, func, rab, asum )
!  !-----------------------------------------------------------------------
!  !
!  !    This routine computes the integral of a function defined on a
!  !    logaritmic mesh, by using the open simpson formula given on
!  !    pag. 109 of Numerical Recipes. In principle it is used to
!  !    perform integrals from zero to infinity. The first point of
!  !    the function should be the closest to zero but not the value
!  !    in zero. The formula used here automatically includes the
!  !    contribution from the zero point and no correction is required.
!  !
!  !    Input as "simpson". At least 8 integrating points are required.
!  !
!  !    last revised 12 May 1995 by Andrea Dal Corso
!  !
!  INTEGER, INTENT(in) :: mesh
!  real(DP), INTENT(in) :: rab (mesh), func (mesh)
!  real(DP), INTENT(out):: asum
!  !
!  real(DP) :: c(4)
!  INTEGER ::i
!  !
!  IF ( mesh < 8 ) CALL throw ('simpson_cp90','few mesh points')
!
!  c(1) = 109.0d0 / 48.d0
!  c(2) = -5.d0 / 48.d0
!  c(3) = 63.d0 / 48.d0
!  c(4) = 49.d0 / 48.d0
!
!  asum = ( func(1)*rab(1) + func(mesh  )*rab(mesh  ) )*c(1) &
!       + ( func(2)*rab(2) + func(mesh-1)*rab(mesh-1) )*c(2) &
!       + ( func(3)*rab(3) + func(mesh-2)*rab(mesh-2) )*c(3) &
!       + ( func(4)*rab(4) + func(mesh-3)*rab(mesh-3) )*c(4)
!  DO i=5,mesh-4
!     asum = asum + func(i)*rab(i)
!  ENDDO
!
!  RETURN
!END SUBROUTINE simpson_cp90
!
!-----------------------------------------------------------------------
!SUBROUTINE herman_skillman_int(mesh,func,rab,asum)
!!-----------------------------------------------------------------------
!  !     simpson rule integration for herman skillman mesh (obsolescent)
!  !     Input as in "simpson". BEWARE: "func" is overwritten!!!
!  !
!  INTEGER, INTENT(in) :: mesh
!  real(DP), INTENT(in) :: rab (mesh)
!  real(DP), INTENT(inout) :: func (mesh)
!  real(DP), INTENT(out):: asum
!  !
!  INTEGER :: i, j, k, i1, nblock
!  REAL(DP) :: a1, a2e, a2o, a2es
!  !
!  a1=0.0d0
!  a2e=0.0d0
!  asum=0.0d0
!  nblock=mesh/40
!  i=1
!  func(1)=0.0d0
!  DO j=1,nblock
!     DO k=1,20
!        i=i+2
!        i1=i-1
!        a2es=a2e
!        a2o=func(i1)/12.0d0
!        a2e=func(i)/12.0d0
!        a1=a1+5.0d0*a2es+8.0d0*a2o-a2e
!        func(i1)=asum+a1*rab(i1)
!        a1=a1-a2es+8.0d0*a2o+5.0d0*a2e
!        func(i)=asum+a1*rab(i)
!     ENDDO
!     asum=func(i)
!     a1=0.0d0
!  ENDDO
!  !
!  RETURN
!END SUBROUTINE herman_skillman_int

real(8) function raddec(alpha,r)
implicit none
real(8), intent(in) :: alpha,r
raddec=2.d0*alpha**(3.d0/2.d0)*exp(-alpha*r)
end function raddec

end module

  
  
