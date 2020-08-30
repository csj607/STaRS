program radiative_transfer
use random
use mpi
use cons
use RT_grid
use RT_photon
use RT_obs
use RT_memory
implicit none

integer, parameter :: nphoton = 9e7
integer, parameter :: nNH = 7
integer, parameter :: nN_XYZ = 7
integer, parameter :: nlevel = 4
integer :: iphoton
integer :: iNH
integer :: iN_XYZ
integer :: idx,ilevel
integer :: NS
logical :: temp_phi
!	MPI
integer :: nproc, rank, irank, ip
integer :: err, num_photon
integer :: ans , tag
integer :: anstype, sender 
integer :: status(MPI_STATUS_SIZE) 
real(kind=rkd) :: time1, time2
real(kind=rkd) :: tau_s

real(kind=rkd), dimension(nNH) :: NH
integer, dimension(nN_XYZ) :: N_XYZ
!real(kind=rkd), dimension(nlevel,nNH,nv) :: flux
!real(kind=rkd), dimension(nlevel,nNH,nv) :: peak_max, peak_wl
real(kind=rkd) :: wlmin, wlmax
real(kind=rkd) :: wlc

real(kind=rkd), parameter :: dNH = 2.d0	
real(kind=rkd),parameter :: dv = 50.d5		! [cm/s]
real(kind=rkd) :: NH_min
real(kind=rkd) :: v_min

type(grid_type) :: grid
type(photon_type) :: photon
type(obs_type) :: obs
type(cross_type) :: cross
type(mpi_type) :: mpar



CALL MPI_INIT(err)

call MPI_COMM_SIZE(MPI_COMM_WORLD,mpar%nproc,err)
call MPI_COMM_RANK(MPI_COMM_WORLD,mpar%rank,err)
call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, mpar%h_comm, ierr)

! slited by each node
call MPI_COMM_RANK(mpar%h_comm, mpar%h_rank,ierr)
call MPI_COMM_SIZE(mpar%h_comm, mpar%h_nproc, ierr)


! slited by same rank in each node
call MPI_COMM_SPLIT(MPI_COMM_WORLD, mpar%h_rank, mpar%rank, mpar%SUB_COMM,ierr)
call MPI_COMM_SIZE(mpar%SUB_COMM, mpar%sub_nproc,ierr)


!	H alpha wing
wlmin = 984.d0
wlmax = 1060.d0

!	H Beta wing
wlmin = 954.d0
wlmax = 984.d0


call set_cross(cross,mpar,wlmin,wlmax)


NH(1) = 1.d20		! [/cm^2]
NH(2) = 1.d20*sqrt(10.d0)
NH(3) = 1.d21
NH(4) = 1.d21*sqrt(10.d0)
NH(5) = 1.d22
NH(6) = 1.d22*sqrt(10.d0)
NH(7) = 1.d23


N_XYZ(1) = 3
N_XYZ(2) = 5
N_XYZ(3) = 11
N_XYZ(4) = 21
N_XYZ(5) = 51
N_XYZ(6) = 101
N_XYZ(7) = 201


!do iNH = 1, nNH
do iNH = 3, 7, 2
do iN_XYZ = 1, nN_XYZ

num_photon = 0

time1 = MPI_WTIME()

!if(rank .ne. master) then
call set_grid(grid,mpar,NH(iNH),N_XYZ(iN_XYZ))
!endif
call set_obs(obs,mpar,nphoton)

if(mpar%rank .eq. master) print*,'start calculation'

call MPI_BARRIER(MPI_COMM_WORLD,err)

if(mpar%rank .eq. master) then	! Master

	do irank = 1, min(nphoton, mpar%nproc-1) 
	ip = irank
	tag = 1
	call MPI_SEND(ip,1,MPI_INTEGER8, irank, tag, MPI_COMM_WORLD,err)
	num_photon = num_photon + 1
	enddo

	do irank = 1, nphoton
	call MPI_RECV(ans,1,MPI_INTEGER8,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,err)
	sender  = status(MPI_SOURCE)
	anstype = status(MPI_TAG)

	if(num_photon .lt. nphoton) then

	if(mod(num_photon,nphoton/10) .eq. 0) then
	time2 = MPI_WTIME()
	print*,num_photon*100/nphoton,sender,time2-time1!,'aa'
	endif


	tag = 1
	num_photon = num_photon + 1
        ip = num_photon
	call MPI_SEND(ip,1,MPI_INTEGER8, sender, tag, MPI_COMM_WORLD,err)


	else

	tag = 0
	call MPI_SEND(MPI_BOTTOM,0,MPI_INTEGER8,sender,tag,MPI_COMM_WORLD,err)

	endif

	enddo

else
call init_random_seed()
	do 607
        call MPI_RECV(ip,1,MPI_INTEGER8,master,MPI_ANY_TAG,MPI_COMM_WORLD,status,err)
	photon%ip = ip
!	print*,rank,'slave'
        if (status(MPI_TAG) .eq.  0) exit
	ans = 1

!	num_photon = num_photon + 1
	call gen_photon(photon, wlmin, wlmax ,grid)

		do 608
		tau_s = -dlog(rand_number())
		call tracing_tau(photon,grid,cross,tau_s)
			if(photon%esc .eqv. .true.) then
			exit
			else
			call scattering(photon, cross, grid)
			endif
		if(photon%esc .eqv. .true.) exit

608		continue

		call collecting_photon(photon,obs)

tag = 1
call MPI_SEND(ans,1,MPI_INTEGER8,master,tag,MPI_COMM_WORLD,err)
607	continue

endif	!	Slave

call MPI_BARRIER(MPI_COMM_WORLD,err)
if(mpar%rank .eq. master) print*,'end of calculation'
call write_file(obs,mpar,iNH,iN_XYZ)



time2 = MPI_WTIME()

if(mpar%rank .eq. master) then
print*,time2-time1,iNH,iN_XYZ
endif


end do	! velocity
end do	! column density



CALL MPI_FINALIZE(err)



end program radiative_transfer
