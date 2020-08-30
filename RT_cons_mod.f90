module cons
use, intrinsic :: iso_fortran_env, only: real32, real64, real128, int32, int64
implicit none
public

integer, parameter :: rkd = real64
real(kind=rkd), parameter :: pi = 4.*atan(1.d0)
real(kind=rkd), parameter :: gamma = 6.265d8
real(kind=rkd), parameter :: k = 1.38064852d-16 ! Boltzmann Constant : cgs
real(kind=rkd), parameter :: f12 = 0.416400     ! Oscillator Strengths
real(kind=rkd), parameter :: m_H = 1.6737236d-24 ! Hydrogen Mass : g 
real(kind=rkd), parameter :: m_e = 9.10938356d-28! Electron Mass : g 
real(kind=rkd), parameter :: e = 1.602d-19*3e9  ! Electric charge : emu
real(kind=rkd), parameter :: c = 2.99792458d10  ! Light Speed : cm/s
real(kind=rkd), parameter :: c_km = 2.99792458d5  ! Light Speed : km/s
real(kind=rkd), parameter :: Lylim = 911.7531d0 
real(kind=rkd), parameter :: Lya = Lylim/(1.d0 - 1.d0/2.d0**2) 
real(kind=rkd), parameter :: Lyb = Lylim/(1.d0 - 1.d0/3.d0**2) 
real(kind=rkd), parameter :: Lyg = Lylim/(1.d0 - 1.d0/4.d0**2) 
real(kind=rkd), parameter :: Lyd = Lylim/(1.d0 - 1.d0/5.d0**2) 
real(kind=rkd), parameter :: Lye = Lylim/(1.d0 - 1.d0/6.d0**2) 
real(kind=rkd), parameter :: Ha = Lylim/(1.d0/4.d0 - 1.d0/3.d0**2) 
real(kind=rkd), parameter :: Hb = Lylim/(1.d0/4.d0 - 1.d0/4.d0**2) 
real(kind=rkd), parameter :: Hg = Lylim/(1.d0/4.d0 - 1.d0/5.d0**2) 
real(kind=rkd), parameter :: He_II_Lyg = 972.134338 
real(kind=rkd), parameter :: He_II_Lyb = 1025.27893
integer,parameter :: master = 0                                                 
integer :: ierr 

type photon_type

real(kind=rkd) :: kx,ky,kz
real(kind=rkd) :: weight                ! weight
real(kind=rkd) :: x,y,z
real(kind=rkd) :: nu
real(kind=rkd) :: vel1, vel2
logical :: esc
integer :: ix,iy,iz
integer :: NS

real(kind=rkd) :: lambda
real(kind=rkd) :: den11, den12, den22 
real(kind=rkd) :: mu,phi 
integer :: level

integer :: ip

real(kind=rkd) :: lambda_i
real(kind=rkd) :: k_i(3) 
real(kind=rkd) :: r_i(3) 

end type photon_type

type grid_type
integer :: N_XYZ
integer :: N_X, N_Y, N_Z

!real(kind=rkd), allocatable :: den(:,:,:)
!real(kind=rkd), allocatable :: vx_grid(:,:,:), vy_grid(:,:,:), vz_grid(:,:,:)
!real(kind=rkd), allocatable :: X_grid(:), Y_grid(:), Z_grid(:)
real(kind=rkd), pointer :: den(:,:,:)
integer :: w_den
real(kind=rkd), pointer :: vx(:,:,:), vy(:,:,:), vz(:,:,:)
integer :: w_vx, w_vy, w_vz
real(kind=rkd), pointer :: X(:), Y(:), Z(:)
integer :: w_x, w_y, w_z
real(kind=rkd) :: dx, dy, dz 

real(kind=rkd) :: Ri, Rxy, Rz, Ro, H
real(kind=rkd) :: NH
end type grid_type



type obs_type
real(kind=rkd) :: ratio
integer :: ndec
integer :: nx, ny 
integer :: nx_spec, ny_spec
integer :: nmu
integer :: nlevel

integer :: npar
integer :: np
real(kind=rkd), pointer :: r_i(:,:),r_f(:,:)
real(kind=rkd), pointer :: k_i(:,:),k_f(:,:)
integer :: w_r_i, w_r_f, w_k_i, w_k_f
real(kind=rkd), pointer :: den11(:), den12(:), den22(:)
real(kind=rkd), pointer :: lambda_i(:), lambda_f(:) 
integer :: w_den11, w_den12, w_den22, w_lambda_i, w_lambda_f
real(kind=rkd), pointer :: par(:,:) 
integer :: w_par
real(kind=rkd), pointer :: photon(:,:)
integer :: w_photon
real(kind=rkd), allocatable :: photon_sum(:,:)



real(kind=rkd), allocatable :: spec(:,:,:), spec_sum(:,:,:)
real(kind=rkd), allocatable :: Qspec(:,:,:), Qspec_sum(:,:,:)
real(kind=rkd), allocatable :: Uspec(:,:,:), Uspec_sum(:,:,:)
real(kind=rkd), allocatable :: spec_total(:,:)
real(kind=rkd), allocatable :: Qspec_total(:,:)
real(kind=rkd), allocatable :: Uspec_total(:,:)

real(kind=rkd), allocatable :: wlmin(:), wlmax(:)
!real(kind=rkd) :: wlmin, wlmax
end type obs_type

type cross_type
integer :: ncross
real(kind=rkd), pointer :: lambda(:) 
integer :: w_lambda
real(kind=rkd), pointer :: sigtot(:),BR2(:),BR3(:),BR4(:)
integer :: w_sigtot, w_BR2, w_BR3, w_BR4
end type cross_type

type mpi_type
integer :: rank, nproc
integer :: h_rank, h_nproc, h_comm 
integer :: SUB_COMM, sub_nproc
end type mpi_type

end module cons
