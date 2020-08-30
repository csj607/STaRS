module RT_grid
use RT_memory
use cons
implicit none
public

public set_grid

contains

subroutine set_grid(grid,mpar,NH,N_XYZ)
type(grid_type) :: grid
type(mpi_type) :: mpar
real(kind=rkd), intent(in) :: NH
integer, intent(in) :: N_XYZ
real(kind=rkd) :: vx,vy,vz
real(kind=rkd) :: Xmin,Xmax,Ymax,Ymin,Zmin,Zmax
real(kind=rkd) :: v_exp
real(kind=rkd) :: temp1 
real(kind=rkd) :: temp2 
real(kind=rkd) :: temp3 
real(kind=rkd) :: x,y,z 
real(kind=rkd) :: pc
integer :: ix,iy,iz

!pc = 3.086d+18		! pc - cm

grid%N_XYZ = N_XYZ
grid%N_X = grid%N_XYZ 
grid%N_Y = grid%N_XYZ 
grid%N_Z = grid%N_XYZ 

call shared_mem(grid%den, [grid%N_X,grid%N_Y,grid%N_Z], grid%w_den, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%vx, [grid%N_X,grid%N_Y,grid%N_Z], grid%w_vx, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%vy, [grid%N_X,grid%N_Y,grid%N_Z], grid%w_vy, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%vz, [grid%N_X,grid%N_Y,grid%N_Z], grid%w_vz, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%x, [grid%N_X+1], grid%w_x, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%y, [grid%N_Y+1], grid%w_y, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(grid%z, [grid%N_Z+1], grid%w_z, mpar%h_rank, mpar%h_comm,ierr)


	grid%NH = NH	! column density [/cm^2]
	grid%Ri = 1.d0
	grid%Ro = 1.d0
	grid%H = 4.d0
!	v_exp = v  ! cm/s
	v_exp = 0.d5 

if(mpar%h_rank .eq. master) then

grid%den = 0.d0
grid%X = 0.d0
grid%Y = 0.d0
grid%Z = 0.d0
grid%vx = 0.d0
grid%vy = 0.d0
grid%vz = 0.d0

	Xmax = grid%Ro
	Xmin = -grid%Ro
	do ix = 1,grid%N_X+1
	grid%X(ix) = Xmin + (Xmax - Xmin)*(ix-1)/(grid%N_X)
	enddo
	Ymax = grid%Ro
	Ymin = -grid%Ro
	do iy = 1,grid%N_Y+1
	grid%Y(iy) = Ymin + (Ymax - Ymin)*(iy-1)/(grid%N_Y)
	enddo
	Zmax = grid%Ro
	Zmin = -grid%Ro
	do iz = 1,grid%N_Z+1
	grid%Z(iz) = Zmin + (Zmax - Zmin)*(iz-1)/(grid%N_Z)
	enddo


	do ix = 1, grid%N_X
	x = (grid%X(ix) + grid%X(ix+1))/2.
	do iy = 1, grid%N_Y
	y = (grid%Y(iy) + grid%Y(iy+1))/2.
	do iz = 1, grid%N_Z
	z = (grid%Z(iz) + grid%Z(iz+1))/2.
	temp1 = sqrt(x**2 + y**2 + z**2)
!	temp2 = abs(z)


!	if(temp1 .gt. grid%Ri .and. temp1 .le. grid%Ro) then
	if(temp1 .le. grid%Ro) then
!	if(temp1 .gt. grid%Ri .and. temp1 .le. grid%Ro .and. temp2 .le. grid%H) then
	grid%den(ix,iy,iz) = grid%NH/grid%Ro
	grid%vx(ix,iy,iz) = (x/temp1)*v_exp
	grid%vy(ix,iy,iz) = (y/temp1)*v_exp
	grid%vz(ix,iy,iz) = (z/temp1)*v_exp 
	endif

	enddo
	enddo
	enddo

endif

call MPI_BARRIER(mpar%h_comm,ierr)


end subroutine set_grid

end module RT_grid
