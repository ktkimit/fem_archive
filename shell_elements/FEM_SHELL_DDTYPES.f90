!==============================================================================
! Derived Data TYPES for shell finite elements
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
MODULE FEM_SHELL_DDTYPES
USE PARAMETERS

! *****************
! * Types for FEM *
! *****************
!------------------------!
! Material property type !
!------------------------!
TYPE :: MPROP_TYPE
	REAL(KP) :: YOUNG
	REAL(KP) :: POISSON
	REAL(KP) :: DENS
END TYPE

!-----------!
! Node type !
!-----------!
TYPE :: NODE_TYPE
	INTEGER(IP) :: EQN(6)   ! equation numbers,
                            ! initially indicate BC's of
                            ! /Ux, Uy, Uz, Tx (Ta), Ty (Tb), Tz/
                            ! where 1 = fixed, 0 = free
	REAL(KP) :: X(3)
    REAL(KP) :: V(3,3)    ! V(1,:) = Shell director vectors
                          ! V(2,:) = base vector for rotation Ta
                          ! V(3,:) = base vector for rotation Tb
    REAL(KP) :: THICK
END TYPE

!--------------!
! Element type !
!--------------!
TYPE :: SHELL4N_TYPE
	INTEGER(IP) :: CONN(4)
	INTEGER(IP) :: MPROP
END TYPE

! ****************************
! * Types for storage format *
! ****************************
!-----------------!
! SSK format type !
!-----------------!
TYPE SSK_TYPE
	INTEGER(IP) :: N    ! matrix dimension
	INTEGER(IP), ALLOCATABLE :: DA(:)
	REAL(KP), ALLOCATABLE :: AA(:)
END TYPE SSK_TYPE

END MODULE FEM_SHELL_DDTYPES

