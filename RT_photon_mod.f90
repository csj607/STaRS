module RT_photon
use cons
use RT_grid
use random
use RT_cross
implicit none
public

public gen_photon
public tracing_tau 
public scattering 

contains

subroutine gen_photon(photon,wlmin,wlmax,grid)
!use random
implicit none
type(photon_type) :: photon
type(grid_type) :: grid
real(kind=rkd) vx,vy,vz
real(kind=rkd) v_cir
real(kind=rkd) dnuH, dnuK
real(kind=rkd) vel 
real(kind=rkd) radius, mu, phi
real(kind=rkd) temp1, temp2 
real(kind=rkd), intent(in) :: wlmin,wlmax

v_cir = 0.d5  ! cm/s (1km = 1e5 cm)

!        photon%ix = int((photon%x - grid%X_grid(1))/(grid%X_grid(grid%N_X+1) - grid%X_grid(1))*grid%N_X) + 1
!        photon%iy = int((photon%y - grid%Y_grid(1))/(grid%Y_grid(grid%N_Y+1) - grid%Y_grid(1))*grid%N_Y) + 1
!        photon%iz = int((photon%z - grid%Z_grid(1))/(grid%Z_grid(grid%N_Z+1) - grid%Z_grid(1))*grid%N_Z) + 1
        photon%x = 0.d0
        photon%y = 0.d0
        photon%z = 0.d0
        photon%ix = grid%N_X/2 + 1
        photon%iy = grid%N_Y/2 + 1
        photon%iz = grid%N_Z/2 + 1
	photon%r_i(1) = photon%x
	photon%r_i(2) = photon%y
	photon%r_i(3) = photon%z
!       initial wavevector
!        photon%mu = (2.d0*rand_number() - 1.d0)*grid%H/sqrt(grid%Ri**2+grid%H**2)
        photon%mu = 2.d0*rand_number() - 1.d0
        photon%phi = 2.d0*pi*rand_number()
        photon%kx = sqrt(1.d0-photon%mu**2)*cos(photon%phi)
        photon%ky = sqrt(1.d0-photon%mu**2)*sin(photon%phi)
        photon%kz = photon%mu
	photon%k_i(1) = photon%kx
	photon%k_i(2) = photon%ky
	photon%k_i(3) = photon%kz

        photon%esc = .false.

	photon%den11 = 0.5d0
	photon%den22 = 0.5d0
	photon%den12 = 0.0d0


!photon%lambda = 1025.27637
!        vx = rand_gauss()
!        vy = rand_gauss()
!        vz = rand_gauss()
!        vel = (photon%kx*vx + photon%ky*vy + photon%kz*vz)*v_cir
!	photon%lambda = rand_number()*(Lyb-Lye) + Lye
	photon%lambda = wlmin + (wlmax - wlmin)*rand_number()
	photon%lambda_i = photon%lambda


        photon%vel1 =    photon%kx*grid%vx(photon%ix,photon%iy,photon%iz) &
                        +photon%ky*grid%vy(photon%ix,photon%iy,photon%iz) &
                        +photon%kz*grid%vz(photon%ix,photon%iy,photon%iz)
	photon%vel2 = 0.d0


	photon%weight = 1.d0
	photon%level = 1
	photon%NS = 0



!!!			photon%lambda = photon%lambda*(1.d0 + photon%vel1/c)		!!!

end subroutine gen_photon




subroutine tracing_tau(photon,grid,cross,tau)
type(photon_type) :: photon
type(cross_type), intent(in) :: cross
type(grid_type), intent(in) :: grid
real(kind=rkd), intent(in) :: tau 
real(kind=rkd) dx,dy,dz,dl
real(kind=rkd) xp,yp,zp 
real(kind=rkd) sigma
real(kind=rkd) tauD,tauDf 


        tauD = 0.d0                                                                                  
	xp = photon%x
	yp = photon%y
	zp = photon%z

                do 202

                if(photon%kx .gt. 0.d0) then
                dx = (grid%X(photon%ix+1)-xp)/photon%kx
                else
                dx = (grid%X(photon%ix)-xp)/photon%kx
                endif
        
                if(photon%ky .gt. 0.d0) then
                dy = (grid%Y(photon%iy+1)-yp)/photon%ky
                else
                dy = (grid%Y(photon%iy)-yp)/photon%ky
                endif

                if(photon%kz .gt. 0.d0) then
                dz = (grid%Z(photon%iz+1)-zp)/photon%kz
                else
                dz = (grid%Z(photon%iz)-zp)/photon%kz
                endif


                dl = min(dx,dy,dz)

if(grid%den(photon%ix,photon%iy,photon%iz) .eq. 0.d0) then

                if(dx .eq. dl) then
                if(photon%kx .gt. 0.d0) then
                photon%ix = photon%ix + 1
                else if(photon%kx .lt. 0.d0) then
                photon%ix = photon%ix - 1
                endif
                endif

                if(dy .eq. dl) then
                if(photon%ky .gt. 0.d0) then
                photon%iy = photon%iy + 1
                else if(photon%ky .lt. 0.d0) then
                photon%iy = photon%iy - 1
                endif
                endif

                if(dz .eq. dl) then
                if(photon%kz .gt. 0.d0) then
                photon%iz = photon%iz + 1
                else if(photon%kz .lt. 0.d0) then
                photon%iz = photon%iz - 1
                endif
                endif

                if(photon%ix .gt. grid%N_x .or. photon%ix .lt. 1) then
                photon%esc = .true.
                exit
                endif

                if(photon%iy .gt. grid%N_y .or. photon%iy .lt. 1) then
                photon%esc = .true.
                exit
                endif

                if(photon%iz .gt. grid%N_z .or. photon%iz .lt. 1) then
                photon%esc = .true.
                exit
                endif


                xp = xp + photon%kx*dl
                yp = yp + photon%ky*dl
                zp = zp + photon%kz*dl

else    ! den .ne. 0

        photon%vel2 =    photon%kx*grid%vx(photon%ix,photon%iy,photon%iz) &
			+photon%ky*grid%vy(photon%ix,photon%iy,photon%iz) &
			+photon%kz*grid%vz(photon%ix,photon%iy,photon%iz)
        photon%lambda = photon%lambda*(1.d0 - photon%vel1/c + photon%vel2/c)
!        sigma = sigma_0/grid%dnu_th*(photon%sigmaK/3.d0 + 2.d0*photon%sigmaH/3.d0)
        photon%vel1 = photon%vel2

	call sigma_tot(photon%lambda,cross,sigma)

                tauDf = tauD + dl*grid%den(photon%ix,photon%iy,photon%iz)*sigma

        if(tauDf .gt. tau) then

             xp = xp + photon%kx*(tau-tauD)/grid%den(photon%ix,photon%iy,photon%iz)/sigma
             yp = yp + photon%ky*(tau-tauD)/grid%den(photon%ix,photon%iy,photon%iz)/sigma
             zp = zp + photon%kz*(tau-tauD)/grid%den(photon%ix,photon%iy,photon%iz)/sigma

        exit

        else


                if(dx .eq. dl) then
                if(photon%kx .gt. 0.d0) then
                photon%ix = photon%ix + 1
                else if(photon%kx .lt. 0.d0) then
                photon%ix = photon%ix - 1
                endif
                endif

                if(dy .eq. dl) then
                if(photon%ky .gt. 0.d0) then
                photon%iy = photon%iy + 1
                else if(photon%ky .lt. 0.d0) then
                photon%iy = photon%iy - 1
                endif
                endif

                if(dz .eq. dl) then
                if(photon%kz .gt. 0.d0) then
                photon%iz = photon%iz + 1
                else if(photon%kz .lt. 0.d0) then
                photon%iz = photon%iz - 1
                endif
                endif


                if(photon%ix .gt. grid%N_x .or. photon%ix .lt. 1) then
                photon%esc = .true.
                exit
                endif

                if(photon%iy .gt. grid%N_y .or. photon%iy .lt. 1) then
                photon%esc = .true.
                exit
                endif

                if(photon%iz .gt. grid%N_z .or. photon%iz .lt. 1) then
                photon%esc = .true.
                exit
                endif


                xp = xp + photon%kx*dl
                yp = yp + photon%ky*dl
                zp = zp + photon%kz*dl
                tauD = tauDf
endif
        endif

202     continue

		photon%x = xp
		photon%y = yp
		photon%z = zp


end subroutine tracing_tau


subroutine scattering(photon, cross, grid)
type(photon_type) :: photon
type(cross_type) :: cross
type(grid_type) :: grid
real(kind=rkd) :: mu, phi
real(kind=rkd) :: mup, phip
real(kind=rkd) :: kxp,kyp,kzp 
real(kind=rkd) :: den11, den12, den22
real(kind=rkd) :: den11p, den12p, den22p
real(kind=rkd) :: trace 
real(kind=rkd) :: temp
real(kind=rkd) :: dphi, sthe, sthep, cdphi, sdphi, c2dphi, s2dphi 
real(kind=rkd) :: BR2, BR3, BR4 


	mu = photon%mu
	phi = photon%phi
	den11 = photon%den11
	den22 = photon%den22
	den12 = photon%den12

	photon%NS = photon%NS + 1


	do 20
 
mup = 2.d0*rand_number() - 1.d0
phip = 2.d0*pi*rand_number()

	dphi  = phip-phi
	sthe  = sqrt(1.-mu**2)      !  sin(theta)
	sthep = sqrt(1.-mup**2)     !  sin(thetap)
	cdphi = cos(dphi)
	sdphi = sin(dphi)
	c2dphi = cos(2.*dphi)
	s2dphi = sin(2.*dphi)

	den11p = den11*cdphi**2 - den12*mu*s2dphi + den22*mu**2*sdphi**2
	den12p = 0.5*den11*mup*s2dphi &
	+ den12*mu*mup*c2dphi + den12*sthe*sthep*cdphi &
	- den22*mu*sthe*sthep*sdphi  &
	- 0.5*den22*mu*mu*mup*s2dphi
	den22p = den11*mup**2*sdphi**2 &
	+ den12*mup*(2.*sthe*sthep*sdphi + mu*mup*s2dphi) &
	+ den22*(mu*mup*cdphi+sthe*sthep)**2

trace =  den11p + den22p
temp =  rand_number()
if(temp .le. trace) exit



20 	continue

	photon%den11 = den11p/trace
	photon%den22 = den22p/trace
	photon%den12 = den12p/trace

	photon%mu = mup
	photon%phi = phip

        photon%kx = sthep*cos(phip)
        photon%ky = sthep*sin(phip)
        photon%kz = mup 

        photon%vel1 =    photon%kx*grid%vx(photon%ix,photon%iy,photon%iz) &
			+photon%ky*grid%vy(photon%ix,photon%iy,photon%iz) &
			+photon%kz*grid%vz(photon%ix,photon%iy,photon%iz)


	call branching(photon%lambda, cross, BR2, BR3, BR4)

	temp = rand_number()
	if(temp .lt. BR2) then
	photon%level = 2
	photon%lambda = 1.d0/(1.d0/photon%lambda - 1.d0/Lya)			!Raman-scattered
        photon%esc = .true.
	else if(temp .ge. BR2 .and. temp .lt. BR2 + BR3) then
	photon%level = 3
	photon%lambda = 1.d0/(1.d0/photon%lambda - 1.d0/Lyb)			!Raman-scattered
        photon%esc = .true.
	else if(temp .ge. BR2 + BR3 .and. temp .lt. BR2 + BR3 + BR4) then
	photon%level = 4
	photon%lambda = 1.d0/(1.d0/photon%lambda - 1.d0/Lyg)			!Raman-scattered
        photon%esc = .true.
	endif




end subroutine scattering



end module RT_photon
