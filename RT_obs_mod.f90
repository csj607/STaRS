module RT_obs
use cons
use RT_memory
implicit none
public

public collecting_photon
public write_file 
public set_obs 

contains

subroutine set_obs(obs,mpar,np)
type(obs_type) :: obs
type(mpi_type) :: mpar
integer, intent(in) :: np 

obs%np = np

obs%npar = 1                    ! 1     lambda_i : initial wavelength
obs%npar = obs%npar +   1       ! 2     lambda_f : final wavelength
obs%npar = obs%npar +   1       ! 3     final energy level 
obs%npar = obs%npar +   3       ! 4-6   r_i  : initial position
obs%npar = obs%npar +   3       ! 7-9   r_f  : final position
obs%npar = obs%npar +   3       ! 10-12  k_i  : inital wavevector 
obs%npar = obs%npar +   3       ! 13-15 k_f  : final wavevector 
obs%npar = obs%npar +   3       ! 16-18 density matrix : den11, den12, den22 
obs%npar = obs%npar +   4       ! 19-22 parameters that you get
call shared_mem(obs%photon, [obs%npar,obs%np], obs%w_photon, mpar%h_rank, mpar%h_comm,ierr)

obs%photon = 0.d0

call MPI_BARRIER(mpar%h_comm,ierr)



end subroutine set_obs

subroutine collecting_photon(photon,obs)
type(obs_type) :: obs
type(photon_type), intent(inout) :: photon
integer :: iph

iph = photon%ip
photon%lambda = photon%lambda*(1.d0 - photon%vel1/c)
!print*,iph,'bbbbb'

        obs%photon(1,iph)  = photon%lambda_i
        obs%photon(2,iph)  = photon%lambda
        obs%photon(3,iph)  = photon%level

        obs%photon(4,iph)  = photon%r_i(1)
        obs%photon(5,iph)  = photon%r_i(2)
        obs%photon(6,iph)  = photon%r_i(3)

        obs%photon(7,iph)  = photon%x
        obs%photon(8,iph)  = photon%y
        obs%photon(9,iph)  = photon%z

        obs%photon(10,iph)  = photon%k_i(1)
        obs%photon(11,iph)  = photon%k_i(2)
        obs%photon(12,iph)  = photon%k_i(3)

        obs%photon(13,iph) = photon%kx
        obs%photon(14,iph) = photon%ky
        obs%photon(15,iph) = photon%kz

        obs%photon(16,iph)  = photon%den11
        obs%photon(17,iph)  = photon%den12
        obs%photon(18,iph)  = photon%den22

        obs%photon(18,iph)  = photon%NS
!        obs%photon(19,iph)  = photon%nclump
!        obs%photon(20,iph)  = photon%path
!        obs%photon(21,iph)  = photon%path




end subroutine collecting_photon

subroutine write_file(obs,mpar,iNH,iv)
use mpi
type(obs_type), intent(inout) :: obs
type(mpi_type), intent(in) :: mpar
integer, intent(in) :: iNH, iv
character :: fn*30
integer :: iph


if(mpar%rank .eq. master) then
if(.not. allocated(obs%photon_sum)) allocate(obs%photon_sum(obs%npar,obs%np))
endif

call MPI_BARRIER(MPI_COMM_WORLD,ierr)

if(mpar%h_rank .eq. master) then
call MPI_REDUCE(obs%photon, obs%photon_sum, obs%npar*obs%np, MPI_DOUBLE_PRECISION,MPI_SUM, master, mpar%sub_comm, ierr)
call MPI_BARRIER(mpar%sub_comm,ierr)
obs%photon = 0.d0
endif


if(mpar%rank .eq. master) then
write(fn,112) 10+iNH,10+iv,'photon.dat'
112	format(I2,I2,A10)
open(19,file=fn)
do iph = 1, obs%np
write(19,*) obs%photon_sum(:,iph) 
enddo
close(19)
endif
call MPI_BARRIER(MPI_COMM_WORLD,ierr)



end subroutine
end module RT_obs
