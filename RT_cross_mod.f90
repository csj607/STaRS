module RT_cross
use cons
use RT_memory
implicit none
integer, parameter :: imax =  1000
integer, parameter :: jmax =  100
public

public sigma_tot
public branching
public set_cross
public sigma_all
public BDSUM
public CONTINT 
public BDMSUM
public CONTMINT
public BD3S
public CONT3S
public BD3D
public CONT3D
public BD4S
public CONT4S
public BD4D
public CONT4D

contains

subroutine sigma_tot(lambda,cross,sigma)
real(kind=rkd), intent(in) :: lambda
type(cross_type), intent(in) :: cross
real(kind=rkd) :: sigma
integer :: idx

idx = int((lambda - cross%lambda(1))/(cross%lambda(cross%ncross) - cross%lambda(1))*(cross%ncross - 1.d0) + 1.d0)
sigma = (cross%sigtot(idx+1) - cross%sigtot(idx))/(cross%lambda(idx+1) - cross%lambda(idx)) &
	*(lambda - cross%lambda(idx)) + cross%sigtot(idx)

end subroutine sigma_tot

subroutine branching(lambda,cross,BR2,BR3,BR4)
real(kind=rkd), intent(in) :: lambda
type(cross_type), intent(in) :: cross
real(kind=rkd) BR2,BR3,BR4
integer :: idx

idx = int((lambda - cross%lambda(1))/(cross%lambda(cross%ncross) - cross%lambda(1))*(cross%ncross - 1.d0) + 1.d0)
BR2 = (cross%BR2(idx+1) - cross%BR2(idx))/(cross%lambda(idx+1) - cross%lambda(idx)) &
	*(lambda - cross%lambda(idx)) + cross%BR2(idx)
BR3 = (cross%BR3(idx+1) - cross%BR3(idx))/(cross%lambda(idx+1) - cross%lambda(idx)) &
	*(lambda - cross%lambda(idx)) + cross%BR3(idx)
BR4 = (cross%BR4(idx+1) - cross%BR4(idx))/(cross%lambda(idx+1) - cross%lambda(idx)) &
	*(lambda - cross%lambda(idx)) + cross%BR4(idx)

end subroutine branching

subroutine set_cross(cross,mpar,wlmin,wlmax)
type(cross_type) :: cross
type(mpi_type), intent(in) :: mpar
real(kind=rkd), intent(in) :: wlmin, wlmax
real(kind=rkd) :: wl
integer :: icross

!wlmin = wlc*(1.d0 - 200.d0/c_km)
!wlmax = wlc*(1.d0 + 300.d0/c_km)

cross%ncross = int(wlmax - wlmin)*10 

call shared_mem(cross%lambda, [cross%ncross], cross%w_lambda, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(cross%sigtot, [cross%ncross], cross%w_sigtot, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(cross%BR2, [cross%ncross], cross%w_BR2, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(cross%BR3, [cross%ncross], cross%w_BR3, mpar%h_rank, mpar%h_comm,ierr)
call shared_mem(cross%BR4, [cross%ncross], cross%w_BR4, mpar%h_rank, mpar%h_comm,ierr)

!allocate(cross%lambda(cross%ncross))
!allocate(cross%sigtot(cross%ncross))
!allocate(cross%BR2(cross%ncross))
!allocate(cross%BR3(cross%ncross))
!allocate(cross%BR4(cross%ncross))

if(mpar%h_rank .eq. master) then

	do icross = 1,cross%ncross
	wl =  wlmin + (wlmax - wlmin)*(icross - 1.d0)/(cross%ncross - 1.d0)
	cross%lambda(icross) = wl
	call sigma_all(wl, cross%sigtot(icross), cross%BR2(icross), cross%BR3(icross), cross%BR4(icross))
	cross%BR2(icross) = cross%BR2(icross)/cross%sigtot(icross)
	cross%BR3(icross) = cross%BR3(icross)/cross%sigtot(icross)
	cross%BR4(icross) = cross%BR4(icross)/cross%sigtot(icross)
	enddo

	if(mpar%rank .eq. master) then
	open(19,file='cross.dat')
	do icross = 1,cross%ncross
	wl =  wlmin + (wlmax - wlmin)*(icross - 1.d0)/(cross%ncross - 1.d0)
	write(19,114) wl, cross%sigtot(icross), cross%BR2(icross), cross%BR3(icross), cross%BR4(icross)
114     format(5ES16.7)
	enddo
	close(19)
	endif

endif

call MPI_BARRIER(mpar%h_comm,ierr)


end subroutine set_cross


SUBROUTINE SIGMA_ALL(WLTH, SIGTOT, SIGRAM, RAM3, RAM4)
!     Cross section for Rayleigh scattering
real(kind=rkd) :: WLTH, SIGTOT, SIGRAM, RAM3, RAM4 
real(kind=rkd) :: WLLIM, VEL_LITE, WI, THOM
real(kind=rkd) :: WRTOWI, WR3TOWI, WR4TOWI 
real(kind=rkd) :: SIG1, SIGB, CINT 
real(kind=rkd) :: SIGMA3D, SIGMA3S 
real(kind=rkd) :: SIGMA4D, SIGMA4S 
real(kind=rkd) :: SIG3D, SIG3S 
real(kind=rkd) :: SIG4D, SIG4S 
real(kind=rkd) :: CINT3D, CINT3S 
real(kind=rkd) :: CINT4D, CINT4S 
      WLLIM = Lylim
         THOM = .66525d-24
         VEL_LITE = 2.998d5
         WI=WLLIM/WLTH
         WRTOWI=1.d0-0.75d0/WI
         WR3TOWI=1.d0-8.d0/9.d0/WI
         WR4TOWI=1.d0-15.d0/16.d0/WI
         CALL BDSUM(WI,SIGB)
         CALL CONTINT(WI,CINT)
         SIG1=THOM * (SIGB+CINT)**2.
!  Branching into n=2s states
         CALL BDMSUM(WI,SIGB)
         CALL CONTMINT(WI,CINT)
         SIGRAM = WRTOWI* THOM *(SIGB+CINT)**2.
       IF (WI .LT. .88889d0) THEN
         SIGTOT = SIG1 + SIGRAM
         RAM3 = 0.
         RAM4 = 0.
         RETURN
       END IF
!  Branching into n=3s, 3d states
         CALL BD3S(WI,SIG3S)
         CALL CONT3S(WI,CINT3S)
         SIGMA3S = WR3TOWI* THOM *(SIG3S+CINT3S)**2.
         CALL BD3D(WI,SIG3D)
         CALL CONT3D(WI,CINT3D)
         SIGMA3D = WR3TOWI* THOM *(SIG3D+CINT3D)**2.
!  Ram3 is the scattering into 3s+3d states.
         RAM3 = SIGMA3S + SIGMA3D
       IF (WI .LT. .93750d0) THEN
         SIGTOT = SIG1 + SIGRAM +RAM3
         RAM4 = 0.d0
         RETURN
       END IF
!  Branching into n=4s, 4d states
         CALL BD4S(WI,SIG4S)
         CALL CONT4S(WI,CINT4S)
         SIGMA4S = WR4TOWI* THOM *(SIG4S+CINT4S)**2.
         CALL BD4D(WI,SIG4D)
         CALL CONT4D(WI,CINT4D)
         SIGMA4D = WR4TOWI* THOM *(SIG4D+CINT4D)**2.
!  Ram4 is the scattering into 4s+4d states.
         RAM4 = SIGMA4S + SIGMA4D
         SIGTOT = SIG1 + SIGRAM + RAM3 + RAM4
!        VEL = (WLLIM/WI - 1025.) / 1025. * VEL_LITE
	RETURN
END SUBROUTINE SIGMA_ALL


SUBROUTINE BDMSUM(WI,SIGB)
REAL(kind=rkd) :: WI, SIGB
REAL(kind=rkd) :: RI, POWER, TEMP 
INTEGER :: J
REAL(kind=rkd) :: A1(jmax),A2(jmax)

      DO J=3,Jmax
         RI=REAL(J)
         POWER=((RI-1.d0)/(RI+1.d0)*(RI-2.d0)/(RI+2.d0))**(RI)
         TEMP=RI**3.d0/(RI*RI-1.d0)/(RI*RI-4.d0)**2
         A1(J)=TEMP*POWER/(1.d0-1.d0/RI/RI-WI)
         A2(J)=TEMP*POWER/(.25d0-1.d0/RI/RI+WI)
      END DO
      SIGB=0.d0
      DO J=3,Jmax
         SIGB=SIGB+A1(J)+A2(J)
      END DO
      SIGB=512.d0*SQRT(2.d0)/3.d0*SIGB
	RETURN

END SUBROUTINE BDMSUM

SUBROUTINE CONTMINT(WI,CINT)
REAL(kind=rkd) :: WI,CINT 
REAL(kind=rkd) :: DI, RI, RIS, POWER, TEMP, TA
REAL(kind=rkd) ::  C1(imax), C2(imax)
INTEGER :: I
      DI=real(Imax/10)
      DO I=1,Imax
         RI=REAL(I)/DI
         RIS=RI*RI
         POWER=ATAN(1.d0/RI)+ATAN(2.d0/RI)
         TEMP=EXP(-2.d0*RI*POWER)
         TA=RI**3.d0*TEMP/(RIS+1.d0)/(RIS+4.d0)**2/(1.-EXP(-2.d0*PI*RI))
         C1(I)=TA/(1.d0+1.d0/RIS-WI)
         C2(I)=TA/(.25+1.d0/RIS+WI)
      END DO
      CINT=0.d0
      DO I=1,Imax
         CINT=CINT+(C1(I)+C2(I))/DI
      END DO
      CINT=512.d0*SQRT(2.d0)/3.d0*CINT
	RETURN

END SUBROUTINE CONTMINT


SUBROUTINE BD3S(WI,SIG3S)
REAL(kind=rkd) :: WI, SIG3S
REAL(kind=rkd) :: RI, RIS, POWER, TEMP 
INTEGER :: J
REAL(kind=rkd) :: A3S1(jmax),A3S2(jmax)

      RI=2.d0
      RIS = RI*RI
      POWER=((RI-1.d0)/(RI+1.d0)*(RI-3.d0)/(RI+3.d0))**(RI)
      TEMP=RI**3.d0*(7.d0*RIS-27.d0)/(RIS-1.d0)/(RIS-9.d0)**3.
      A3S1(2)=TEMP*POWER/(1.d0-1.d0/RI/RI-WI)
      A3S2(2)=TEMP*POWER/(.25d0-1.d0/RI/RI+WI)

      SIG3S = A3S1(2) + A3S2(2)
      DO J=4,Jmax
         RI=REAL(J)
         RIS = RI*RI
         POWER=((RI-1.d0)/(RI+1.d0)*(RI-3.d0)/(RI+3.d0))**(RI)
         TEMP=RI**3*(7.d0*RIS-27.d0)/(RIS-1.d0)/(RIS-9.d0)**3.
         A3S1(J)=TEMP*POWER/(1.d0-1.d0/RI/RI-WI)
         A3S2(J)=TEMP*POWER/(.25d0-1.d0/RI/RI+WI)
      END DO

      DO J=4,Jmax
         SIG3S=SIG3S+A3S1(J)+A3S2(J)
      END DO
!      SIG3S=256.*3.*SQRT(3.)/2./3.*SIG3S
      SIG3S = 128.d0*SQRT(3.d0)*SIG3S
!     angular integration sqrt(1./9.)
	RETURN

END SUBROUTINE BD3S

SUBROUTINE CONT3S(WI,CINT3S)
REAL(kind=rkd) :: WI,CINT3S 
REAL(kind=rkd) :: DI, RI, RIS, POWER, TEMP, TA
REAL(kind=rkd) ::  C3S1(imax), C3S2(imax)
INTEGER :: I
      DI=real(Imax/10)
      DO I=1,Imax
         RI=REAL(I)/DI
         RIS=RI*RI
         POWER=ATAN(1.d0/RI)+ATAN(3.d0/RI)
         TEMP=EXP(-2.*RI*POWER)
         TA = RI**3*(7.d0*RIS+27.d0)*TEMP/(RIS+1.d0)/(RIS+9.d0)**3
         TA = TA/(1.d0-EXP(-2.d0*PI*RI))
         C3S1(I)=TA/(1.d0+1.d0/RIS-WI)
         C3S2(I)=TA/(.25d0+1.d0/RIS+WI)
      END DO
      CINT3S = 0.d0
      DO I=1,Imax
         CINT3S = CINT3S+(C3S1(I)+C3S2(I))/DI
      END DO
!     CINT3S=256.*3.*SQRT(3.)/2./3.*CINT3S
      CINT3S = 128.d0*SQRT(3.d0)*CINT3S
!     angular integration sqrt(1./9.)
	RETURN

END SUBROUTINE CONT3S


SUBROUTINE BD3D(WI,SIG3D)
REAL(kind=rkd) :: WI, SIG3D
REAL(kind=rkd) :: RI, RIS, POWER, TEMP 
INTEGER :: J
REAL(kind=rkd) :: A3D1(jmax),A3D2(jmax)
      RI=2.
      RIS = RI*RI
      POWER=((RI-1.d0)/(RI+1.d0)*(RI-3.d0)/(RI+3.d0))**(RI)
      TEMP=RI**5/(RIS-1.d0)/(RIS-9.d0)**3
      A3D1(2)=TEMP*POWER/(1.d0-1.d0/RI/RI-WI)
      A3D2(2)=TEMP*POWER/(.25d0-1.d0/RI/RI+WI)
      SIG3D = A3D1(2) + A3D2(2)
      DO J=4,Jmax
         RI=REAL(J)
         RIS = RI*RI
         POWER=((RI-1.d0)/(RI+1.d0)*(RI-3.d0)/(RI+3.d0))**(RI)
         TEMP=RI**5/(RIS-1.d0)/(RIS-9.d0)**3
         A3D1(J)=TEMP*POWER/(1.d0-1.d0/RI/RI-WI)
         A3D2(J)=TEMP*POWER/(.25d0-1.d0/RI/RI+WI)
      END DO

      DO J=4,Jmax
         SIG3D = SIG3D + A3D1(J) + A3D2(J)
      END DO
!     SIG3D=512.*SQRT(3.)*/SQRT(10.)*SQRT(2.)*SIG3D
      SIG3D=512.d0*SQRT(3.d0)/SQRT(5.d0)*SIG3D
	RETURN

END SUBROUTINE BD3D

SUBROUTINE CONT3D(WI,CINT3D)
REAL(kind=rkd) :: WI,CINT3D 
REAL(kind=rkd) :: DI, RI, RIS, POWER, TEMP, TA
REAL(kind=rkd) ::  C3D1(imax), C3D2(imax)
INTEGER :: I
      DI=REAL(Imax/10)
      DO I=1,Imax
         RI=REAL(I)/DI
         RIS=RI*RI
         POWER=ATAN(1.d0/RI)+ATAN(3.d0/RI)
         TEMP=EXP(-2.d0*RI*POWER)
         TA=RI**5.d0*TEMP/(RIS+1.d0)/(RIS+9.d0)**3/(1.d0-EXP(-2.d0*PI*RI))
         C3D1(I)=TA/(1.d0+1.d0/RIS-WI)
         C3D2(I)=TA/(.25d0+1.d0/RIS+WI)
      END DO
      CINT3D=0.d0
      DO I=1,Imax
         CINT3D = CINT3D + (C3D1(I)+C3D2(I))/DI
      END DO
!     CINT3D=512.*SQRT(3)*/SQRT(10.)*SQRT(2.)*CINT3D
      CINT3D=512.d0*SQRT(3.d0)/SQRT(5.d0)*CINT3D
	RETURN

END SUBROUTINE



SUBROUTINE BDSUM(WI,SIGB)
REAL(kind=rkd) :: WI, SIGB
REAL(kind=rkd) :: RI, RIS, POWER
INTEGER :: J
REAL(kind=rkd) :: A1(jmax),A2(jmax)
      DO J=2,jmax
         RI=REAL(J)
         RIS=RI*RI
         POWER=RI**3*((RI-1.d0)/(RI+1.d0))**(2.d0*RI)/(RIS-1.d0)**3.
         A1(J)=WI*POWER/(1.d0-1.d0/RIS-WI)*RIS/(RIS-1.d0)
         A2(J)=WI*POWER/(1.d0-1.d0/RIS+WI)*RIS/(RIS-1.d0)
      END DO
      SIGB=0.d0
      DO J=2,100
         SIGB=SIGB+A1(J)-A2(J)
      END DO
      SIGB=128.d0/3.d0*SIGB
	RETURN

END SUBROUTINE BDSUM
      

SUBROUTINE CONTINT(WI,CINT)
REAL(kind=rkd) :: WI,CINT
REAL(kind=rkd) :: DI, RI, RIS, TEMP, TA
REAL(kind=rkd) ::  C1(imax), C2(imax)
INTEGER :: I
      DI=real(Imax/10)
      DO I=1,Imax
         RI=REAL(I)/DI
         RIS=RI*RI
         TEMP=EXP(-4.d0*RI*ATAN(1.d0/RI))
         TA=RI*TEMP*(RI/(RIS+1.d0))**4/(1.d0-EXP(-2.d0*PI*RI))
         C1(I)=WI*TA/((RIS+1.d0)/RIS-WI)
         C2(I)=WI*TA/((RIS+1.d0)/RIS+WI)
      END DO
      CINT=0.d0
      DO I=1,1000
         CINT=CINT+(C1(I)-C2(I))/DI
      END DO
      CINT=128.d0/3.d0*CINT
	RETURN

END SUBROUTINE CONTINT


SUBROUTINE BD4S(WI,SIG4S)
REAL(kind=rkd) :: WI, SIG4S
REAL(kind=rkd) :: RI, RIS, POWER, TEMP 
INTEGER :: J
REAL(kind=rkd) :: A4S1(jmax),A4S2(jmax)

      RI=2.d0
      RIS = RI*RI
      POWER=((RI-1.d0)/(RI+1.d0)*(4.d0-RI)/(RI+4.d0))**(RI)
      TEMP=RI**3.d0/(RIS-1.d0)/(RIS-16.d0)**4.d0*(768.d0-288.d0*RIS+23.d0*RIS*RIS)
      A4S1(2)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
      A4S2(2)=TEMP*POWER/(.25d0-1.d0/RIS+WI)

      RI=3.d0
      RIS = RI*RI
      POWER=-1.d0*((RI-1.d0)/(RI+1.d0)*(4.d0-RI)/(RI+4.d0))**(RI)
      TEMP=RI**3/(RIS-1.d0)/(RIS-16.d0)**4*(768.d0-288.d0*RIS+23.d0*RIS*RIS)
      A4S1(3)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
      A4S2(3)=TEMP*POWER/(.25d0-1.d0/RIS+WI)

      SIG4S = A4S1(2)+A4S2(2) + A4S1(3)+A4S2(3)
      DO J=5,jmax
        RI=REAL(J)
        RIS = RI*RI
        POWER=((RI-1.d0)/(RI+1.d0)*(RI-4.d0)/(RI+4.d0))**(RI)
        TEMP=RI**3/(RIS-1.d0)/(RIS-16.d0)**4*(768.d0-288.d0*RIS+23.d0*RIS*RIS)
        A4S1(J)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
        A4S2(J)=TEMP*POWER/(.25d0-1.d0/RIS+WI)
      END DO

      DO J=5,Jmax
         SIG4S=SIG4S+A4S1(J)+A4S2(J)
      END DO
      SIG4S = 4096.d0/9.d0*SIG4S
!     angular integration sqrt(1./9.)
	RETURN

END SUBROUTINE BD4S

SUBROUTINE CONT4S(WI,CINT4S)
REAL(kind=rkd) :: WI,CINT4S 
REAL(kind=rkd) :: DI, RI, RIS, POWER, TEMP, TA, RISS
REAL(kind=rkd) ::  C4S1(imax), C4S2(imax)
INTEGER :: I
      DI=REAL(Imax/10)
      DO I=1,Imax
        RI=REAL(I)/DI
        RIS=RI*RI
        RISS=RIS*RIS
        POWER=ATAN(1.d0/RI)+ATAN(4.d0/RI)
        TEMP=EXP(-2.d0*RI*POWER)
        TA = RI**3*TEMP/(RIS+1.d0)/(RIS+16.d0)**4*(768.d0+288.d0*RIS+23.d0*RISS)
        TA = TA/(1.d0-EXP(-2.d0*PI*RI))
        C4S1(I)=TA/(1.d0+1.d0/RIS-WI)
        C4S2(I)=TA/(.25d0+1.d0/RIS+WI)
      END DO
      CINT4S = 0.d0
      DO I=1,Imax
        CINT4S = CINT4S+(C4S1(I)+C4S2(I))/DI
      END DO
      CINT4S = 4096.d0/9.d0*CINT4S
      RETURN
END SUBROUTINE CONT4S

SUBROUTINE BD4D(WI,SIG4D)
REAL(kind=rkd) :: WI, SIG4D
REAL(kind=rkd) :: RI, RIS, POWER, TEMP 
INTEGER :: J
REAL(kind=rkd) :: A4D1(jmax),A4D2(jmax)
      RI=2.d0
      RIS = RI*RI
      POWER=((RI-1.d0)/(RI+1.d0)*(4.d0-RI)/(RI+4.d0))**(RI)
      TEMP=RI**5/(RIS-1.d0)/(RIS-16.d0)**4*(7.d0*RIS-48.d0)
      A4D1(2)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
      A4D2(2)=TEMP*POWER/(.25d0-1.d0/RIS+WI)
      SIG4D = A4D1(2) + A4D2(2)

      RI=3.d0
      RIS = RI*RI
      POWER=-1.d0*((RI-1.d0)/(RI+1.d0)*(4.d0-RI)/(RI+4.d0))**(RI)
      TEMP=RI**5/(RIS-1.d0)/(RIS-16.d0)**4.*(7.d0*RIS-48.d0)
      A4D1(2)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
      A4D2(2)=TEMP*POWER/(.25d0-1.d0/RIS+WI)
      SIG4D = SIG4D+ A4D1(3) + A4D2(3)
      DO J=5,jmax
         RI=REAL(J)
         RIS = RI*RI
         POWER=((RI-1.d0)/(RI+1.d0)*(RI-4.d0)/(RI+4.d0))**(RI)
         TEMP=RI**5/(RIS-1.d0)/(RIS-16.d0)**4*(7.d0*RIS-48.d0)
         A4D1(J)=TEMP*POWER/(1.d0-1.d0/RIS-WI)
         A4D2(J)=TEMP*POWER/(.25d0-1.d0/RIS+WI)
      END DO
      DO J=5,100
         SIG4D = SIG4D + A4D1(J) + A4D2(J)
      END DO
      SIG4D=8192.d0/9.d0/SQRT(5.d0)*SQRT(2.d0)*SIG4D
	RETURN

END SUBROUTINE BD4D


SUBROUTINE CONT4D(WI,CINT4D)
REAL(kind=rkd) :: WI,CINT4D 
REAL(kind=rkd) :: DI, RI, RIS, POWER, TEMP, TA
REAL(kind=rkd) ::  C4D1(imax), C4D2(imax)
INTEGER :: I
      DI=real(Imax/10)
      DO I=1,Imax
         RI=REAL(I)/DI
         RIS=RI*RI
         POWER=ATAN(1.d0/RI)+ATAN(4.d0/RI)
         TEMP=EXP(-2.d0*RI*POWER)*(7.d0*RIS+48.d0)
         TA=RI**5/(RIS+1.d0)/(RIS+16.d0)**4*TEMP/(1.d0-EXP(-2.d0*PI*RI))
         C4D1(I)=TA/(1.d0+1.d0/RIS-WI)
         C4D2(I)=TA/(.25d0+1.d0/RIS+WI)
      END DO
      CINT4D=0.d0
      DO I=1,IMAX
         CINT4D = CINT4D + (C4D1(I)+C4D2(I))/DI
      END DO
      CINT4D=8192.d0/9.d0/SQRT(5.d0)*SQRT(2.d0)*CINT4D
	RETURN

END SUBROUTINE CONT4D


end module RT_cross
