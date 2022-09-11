!==============================================================================
! Derived Data TYPES for plane stress & plane strain FEM
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
MODULE FEM_PLANE_DDTYPES
USE PARAMETERS

! *****************
! * Types for FEM *
! *****************
!-----------------------
! Material property type 
!-----------------------
TYPE :: MPROP_TYPE
	REAL(KP) :: YOUNG
	REAL(KP) :: POISSON
	REAL(KP) :: THICK
	REAL(KP) :: DENS
END TYPE

!----------
! Node type
!----------
TYPE :: NODE_TYPE
	INTEGER(ip) :: BC(2)	 ! 1 : fixed, 0 : free
							 ! /Ux, Uy/ 
	INTEGER(ip) :: EQN(2)   ! equation numbers
	REAL(KP) :: X(2)
END TYPE

!----------------------
! Neumann boundary type
!----------------------
TYPE :: NBOUNDARY3N_TYPE
	INTEGER(ip) :: EN    ! element #
	INTEGER(ip) :: CONN(3)
	REAL(KP) :: UF(2)    ! Uniform surface force; [fx, fy]
END TYPE

! ------------
! Element type
! ------------
TYPE :: QUAD9N_TYPE
	INTEGER(ip) :: CONN(9)
	INTEGER(ip) :: MPROP
END TYPE

TYPE :: QUAD4N_TYPE
	INTEGER(ip) :: CONN(4)
	INTEGER(ip) :: MPROP
END TYPE

! ****************************
! * Types for storage format *
! ****************************
! ---------------
! SSK format type
! ---------------
TYPE SSK_TYPE
	INTEGER(ip) :: N    ! matrix dimension
	INTEGER(ip), ALLOCATABLE :: DA(:)
	REAL(KP), ALLOCATABLE :: AA(:)
END TYPE SSK_TYPE

END MODULE FEM_PLANE_DDTYPES

