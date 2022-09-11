!==============================================================================
! Derived Data TYPES for 1D Beam FEM
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
MODULE FEM_1DBEAM_DDTYPES
USE PARAMETERS

! *****************
! * Types for FEM *
! *****************
!-----------------------
! Material property type 
!-----------------------
TYPE :: MPROP_TYPE
	REAL(KP) :: YOUNG    ! Young's modulus
	REAL(KP) :: I		 ! second moment of area
	REAL(KP) :: SHEAR	 ! shaer modulus
	REAL(KP) :: AREA	 ! sectional area
	REAL(KP) :: K		 ! Timoshenko shear coefficient
	REAL(KP) :: RHO		 ! mass density
END TYPE

!----------
! Node type
!----------
TYPE :: NODE_TYPE
	INTEGER(I4B) :: BC(2)	 ! 1 : fixed, 0 : free
						 	 ! (/ Dw, Dt /)
	INTEGER(I4B) :: EQN(2)   ! equation numbers
						 	 ! (/ Dw, Dt /)
	REAL(KP) :: X
END TYPE

! ------------
! Element type
! ------------
TYPE :: LIN2N_TYPE
	INTEGER(I4B) :: CONN(2)
	INTEGER(I4B) :: MPROP
END TYPE

! ****************************
! * Types for storage format *
! ****************************
! ---------------
! SSK format type
! ---------------
TYPE SSK_TYPE
	INTEGER(I4B) :: N    ! matrix dimension
	INTEGER(I4B), ALLOCATABLE :: DA(:)
	REAL(KP), ALLOCATABLE :: AA(:)
END TYPE SSK_TYPE

END MODULE FEM_1DBEAM_DDTYPES
