!==============================================================================
! Derived Data TYPES for 3D Solid elements
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
MODULE FEM_3DSOLID_DDTYPES
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
	INTEGER(IP) :: EQN(3)   ! equation numbers,
                            ! initially indicate BC's of /Ux, Uy, Uz/
                            ! where 1 = fixed, 0 = free
	REAL(KP) :: X(3)
END TYPE

!--------------!
! Element type !
!--------------!
TYPE :: HEX8N_TYPE
	INTEGER(IP) :: CONN(8)
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

END MODULE FEM_3DSOLID_DDTYPES

