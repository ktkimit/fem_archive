!==============================================================================
! FEA of plane stress & plane strain problems
! - 2D 4-nodes elements
!
! CONTENTS:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!==============================================================================
PROGRAM FEM_PLANE_4N_MAIN
USE PARAMETERS
USE FEM_PLANE_DDTYPES
USE FEM_PLANE_4N, ONLY : FEM_PLANE_STIFF4N, FEM_PLANE_MASS4N
IMPLICIT NONE

! File name
CHARACTER(20) :: FNAME

! Analysis starting time
CHARACTER(10) :: DATE(3)    ! date, time, zone
INTEGER(ip)  :: DATEVAL(8) ! yy, mm, dd, time diff wrt UTC,
							! h, m, sec, millisec

! Problem type
CHARACTER(12) :: PTYPE

TYPE(MPROP_TYPE), ALLOCATABLE :: MATERIAL(:)
TYPE(NODE_TYPE), ALLOCATABLE :: NODE(:)
TYPE(QUAD4N_TYPE), ALLOCATABLE :: QUAD4N(:)

! Total number of -
INTEGER(ip) :: TN_NODE    ! nodes
INTEGER(ip) :: TN_ELEMENT    ! elements
INTEGER(ip) :: TN_MSET    ! material sets

! Degrees of freedom
INTEGER(ip) :: TN_DOF    ! # of total dof
INTEGER(ip) :: TN_FDOF   ! # of free dof
INTEGER(ip) :: TN_PDOF   ! # of prescribed (Dirichlet) dof
INTEGER(ip) :: NDOF   ! # of dof per node
PARAMETER( NDOF = 2 )

! Array for stiffness and mass matrices
TYPE(SSK_TYPE) :: KT, MT

! Timer
REAL(KP) :: ITIME, ETIME

! Compute analysis starting time
call DATE_AND_TIME( date(1), date(2), date(3), dateval )

! Execute main routine
call MAIN_FEM()

CONTAINS
!==========================================================
!----------------------------------------------------------
! Main routine
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE MAIN_FEM()
	!--------------
	! Print heading
	!--------------
	call WRITE_HEADING(0)

	!----------------
	! Read data files
	!----------------
	write(*, '("A. Reading date files")')
	call READ_INPUTDATA()    ! allocate: NODE, MATERIAL, QUAD4N

	!--------------------------------------------
	! Construct linear system of equations by MFR
	!--------------------------------------------
	write(*, '("B. Contructing linear system of equations by FEM")')

	! Set equation number
	write(*, '(4X, "1. Setting equation numbers")')
	call SET_EQN()

	! Allocate arrays
	write(*, '(4X, "2. Allocating arrays")')
	call ALL_KMSSK()    ! allocate: KT.DA, KT.AA, MT.DA, MT.AA
	KT.AA = 0.0_kp
	MT.AA = 0.0_kp

	call WRITE_MEMORY(0)

	! Calculate and assemble total stiffness matrix
	write(*, '(4X, "3. Calculating and assembling stiffness matrix")')
	call cpu_time(itime)

	call ASSEMBLE_KMT_FEM_PLANE()

	call cpu_time(etime)
	print *, "Elapsed time: ", etime - itime

	!-----------------------
	! Export data into files
	!-----------------------
	call WRITE_KM()

	deallocate( NODE )
	deallocate( MATERIAL )
	deallocate( QUAD4N )
	deallocate( KT.DA, KT.AA )
	deallocate( MT.DA, MT.AA )
END SUBROUTINE MAIN_FEM

! ****************************************
! * Construct linear system of equations *
! ****************************************
!----------------------------------------------------------
! Allocate DA and AA arrays (SSK format) of K and M
!
! ALLOCATED VARIABLES:
!	KT.DA, MT.DA
!	KT.AA, MT.AA
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 10, February, 2015
!----------------------------------------------------------
SUBROUTINE ALL_KMSSK()
	integer(ip) :: el, i, j
	integer(ip) :: eqn, meqn
	integer(ip), allocatable :: first_nzeroi(:)    ! first nonzero eqn in each column

	KT.N = TN_DOF; MT.N = TN_DOF
	allocate( KT.DA(TN_DOF + 1), MT.DA(TN_DOF + 1) )
	allocate( first_nzeroi(TN_DOF) )

	first_nzeroi(:) = TN_DOF
	do el = 1, TN_ELEMENT
		! Find the minimum eqn
		meqn = TN_DOF
		do i = 1, 4
			do j = 1, NDOF
				eqn = NODE( QUAD4N(el).CONN(i) ).EQN(j)
				if ( eqn < meqn ) then
					meqn = eqn
				end if
			end do
		end do

		do i = 1, 4
			do j = 1, NDOF
				eqn = node( QUAD4N(el).CONN(i) ).EQN(j)
				if ( first_nzeroi(eqn) > meqn ) then
					first_nzeroi(eqn) = meqn
				end if
			end do
		end do
	end do

	KT.DA(1) = 1; MT.DA(1) = 1
	do i = 1, TN_DOF
		KT.DA(i+1) = KT.DA(i) + i - first_nzeroi(i) + 1
		MT.DA(i+1) = MT.DA(i) + i - first_nzeroi(i) + 1
	end do

	allocate( KT.AA( KT.DA(TN_DOF+1) - 1 ) )
	allocate( MT.AA( MT.DA(TN_DOF+1) - 1 ) )

        deallocate( first_nzeroi )
END SUBROUTINE ALL_KMSSK
!----------------------------------------------------------
! Calculate & assemble KT & MT arrays
!
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!----------------------------------------------------------
SUBROUTINE ASSEMBLE_KMT_FEM_PLANE()
	integer(ip) :: n, i, j, k
	integer(ip) :: id, jd, er, ec, r, c
	integer(ip) :: address
	real(kp) :: enode(4,2)
	real(kp) :: Ke(8,8), Me(8,8)

	do n = 1, TN_ELEMENT
		do k = 1, 4
			enode(k,:) = NODE( QUAD4N(n).CONN(k) ).X(:)
		end do

		call FEM_PLANE_STIFF4N(enode,PTYPE,&
			MATERIAL( QUAD4N(n).MPROP ).YOUNG,&
		 	MATERIAL( QUAD4N(n).MPROP ).POISSON,&
		 	MATERIAL( QUAD4N(n).MPROP ).THICK,&
			Ke)
		call FEM_PLANE_MASS4N(enode,&
			MATERIAL( QUAD4N(n).MPROP ).DENS,&
			MATERIAL( QUAD4N(n).MPROP ).THICK,&
			Me)

		do i = 1, 4
			do id = 1, NDOF
				er = (id-1)*4 + i
				r  = NODE( QUAD4N(n).CONN(i) ).EQN(id)
				do j = 1, 4
					do jd = 1, NDOF
						ec = (jd-1)*4 + j
						c  = NODE( QUAD4N(n).CONN(j) ).EQN(jd)

						if ( c >= r ) then
							address = KT.DA(c) + c - r
							KT.AA(address) = KT.AA(address) + Ke(er,ec)
							address = MT.DA(c) + c - r
							MT.AA(address) = MT.AA(address) + Me(er,ec)
						end if
					end do
				end do
			end do
		end do
	end do
END SUBROUTINE ASSEMBLE_KMT_FEM_PLANE

! **************************************
! * Miscellaneous computation routines *
! **************************************
!----------------------------------------------------------
! Set eqn #
!
! Author:
!	Ki-Tae Kim (qlsn55@gmail.com), 6, March, 2014
!----------------------------------------------------------
SUBROUTINE SET_EQN()
	integer(ip) :: i, j

	TN_DOF  = 0
	TN_FDOF = 0
	TN_PDOF = 0

	do i = 1, TN_NODE
		do j = 1, NDOF
			TN_DOF = TN_DOF + 1
		end do
	end do

	do i = 1, TN_NODE
		do j = 1, NDOF
			if (NODE(i).BC(j) == 0) then
				TN_FDOF = TN_FDOF + 1
				NODE(i).EQN(j) = TN_FDOF
			else if (NODE(i).BC(j) == 1) then
				TN_PDOF = TN_PDOF + 1
				NODE(i).EQN(j) = TN_DOF - TN_PDOF + 1
			end if
		end do
	end do
END SUBROUTINE SET_EQN

! ********************
! * Reading routines *
! ********************
!----------------------------------------------------------
! Read input data from formatted files
!
! ALLOCATED VARIABLES:
!	NODE
!	MATERIAL
!	QUAD4N
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 20, January, 2015
!----------------------------------------------------------
SUBROUTINE READ_INPUTDATA()
	character(25) :: dums
	integer(ip) :: i, n

	!-----------------
	! Input files name
	!-----------------
	write(*, '(4X, "Input the file name", &
		1X, "(without extension)")')
	write(*, '(8X, ">")', advance='no')
	read(*,*) FNAME

	!-----------
	! Open files
	!-----------
	open(unit=1, file=trim(FNAME)//".node", status='old', &
		action='read', form='formatted')
	open(unit=2, file=trim(FNAME)//".element", status='old', &
		action='read', form='formatted')

	!---------------
	! Read node data
	!---------------
	do i = 1, 6
		read(1,*)
	end do

	! Total # of nodes
	read(1,*) TN_NODE

	do i = 1, 10
		read(1,*)
	end do

	allocate( NODE(TN_NODE) )

	! Read
	do i = 1, TN_NODE
		read(1,*) n, NODE(n).X(:), NODE(n).BC(:)
	end do

	!------------------------------
	! Read material properties data
	!------------------------------
	do i = 1, 7
		read(2,*)
	end do

	! Problem type
	read(2,'(A)') dums
	PTYPE = ADJUSTL(TRIM(dums))

	do i = 1, 3
		read(2,*)
	end do

	! Total # of material sets
	read(2,*) TN_MSET

	allocate( MATERIAL(TN_MSET) )

	do i = 1, 5
		read(2,*)
	end do

	! Read
	do i = 1, TN_MSET
		read(2,*) n, MATERIAL(n).YOUNG, MATERIAL(n).POISSON, &
			MATERIAL(n).THICK, MATERIAL(n).DENS
	end do

	!------------------
	! Read element data
	!------------------
	do i = 1, 7
		read(2,*)
	end do

	! Total # of elements
	read(2,*) TN_ELEMENT

	allocate( QUAD4N(TN_ELEMENT) )

	do i = 1, 5
		read(2,*)
	end do

	! Read
	do i = 1, TN_ELEMENT
		read(2,*) n, QUAD4N(n).CONN(:), QUAD4N(n).MPROP
	end do
	
	!------------
	! Close files
	!------------
	close(unit=1)
	close(unit=2)
END SUBROUTINE READ_INPUTDATA

! ********************
! * Writing routines *
! ********************
!----------------------------------------------------------
! Write heading
!
! INPUT:
!	n : unit specifier in write statement
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 19, January, 2015
!----------------------------------------------------------
SUBROUTINE WRITE_HEADING(n)
	integer(ip), intent(in) :: n

	write(n, '( 8X, 61("*") )')
	write(n, '( 8X, "*", 59X, "*" )')
	write(n, '( 8X, "* Plane Stress & Strain Analysis by 4-nodes ", &
		"Finite Elements", X, "*")')
	write(n, '( 8X, "*", 59X, "*" )')
	write(n, '( 8X, "* Program author : Ki-Tae Kim,", 30X, "*" )')
	write(n, '( 8X, "*", 18X, "Department of Mechanical Engineering,", 4X, &
		"*" )')
	write(n, '( 8X, "*", 18X, "Massachusetts Institute of Technology", 4X, &
		"*" )')
	write(n, '( 8X, "*", 59X, "*" )')
	write(n, '( 8X, "* E-mail         : ktkim@mit.edu or qlsn55@gmail.com", &
		        8X, "*" )')
	write(n, '( 8X, "*", 59X, "*" )')
	Write(n, '( 8X, "* Analysis conducting time : ", I2, ":", I2, ", ", I2, & 
		"/", I2, "/", I4, 13X, " *" )' ) dateval(5), dateval(6), &
		dateval(2), dateval(3), dateval(1)
	write(n, '( 8X, 61("*") )')
END SUBROUTINE WRITE_HEADING
!----------------------------------------------------------
! Write the size of memory allocated to major arrays
!
! INPUT:
!	n : unit specifier in write statement
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 10, February, 2015
!----------------------------------------------------------
SUBROUTINE WRITE_MEMORY(n)
	integer(ip), intent(in) :: n

	write(*, '(8X, "---------------------")')
	write(*, '(8X, "TN_DOF  : ", I7)') TN_DOF
	write(*, '(8X, "TN_FDOF : ", I7)') TN_FDOF
	write(*, '(8X, "TN_PDOF : ", I7)') TN_PDOF
	write(*, '()')
	write(*, '(8X, "Allocated memory size")')
	write(*, 10) int(size(KT.AA) * 8.0E-3_kp +&
				 size(KT.DA) * 4.0E-3_kp)
	write(*, 11) int(size(MT.AA) * 8.0E-3_kp +&
				 size(MT.DA) * 4.0E-3_kp)
!	write(*, 20) int(size(RT) * 8.0E-3_kp)
!	write(*, 30) int(size(UT) * 8.0E-3_kp)
	10 format (8X, "KT :", 1X, I8, "KB")
	11 format (8X, "MT :", 1X, I8, "KB")
!	20 format (8X, "RT :", 1X, I7, "KB")
!	30 format (8X, "UT :", 1X, I7, "KB")
	write(*, '(8X, "---------------------")')
END SUBROUTINE WRITE_MEMORY
!----------------------------------------------------------
! Export stiffness and mass arrays to binary files
! - save as:
!	TN_DOF = *.N
!	nnzero = *.DA( *.N + 1 ) - 1
!	*.DA(1:TN_DOF)
!	*.AA(1:nnzero)
!
! OUTPUT:
!	*.stiff
!	*.mass
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 17, April, 2016
!----------------------------------------------------------
SUBROUTINE WRITE_KM()
	integer(ip) :: nwk, nwm

	open(unit=1, file=trim(FNAME)//".stiff", status='replace', &
		form='unformatted', action='write', access='stream')
	open(unit=2, file=trim(FNAME)//".mass", status='replace', &
		form='unformatted', action='write', access='stream')

	nwk = KT.DA(TN_FDOF+1) - 1
	nwm = nwk

	write(1) TN_FDOF
	write(1) nwk
	write(1) KT.DA(1:TN_FDOF+1)
	write(1) KT.AA(1:nwk)

	write(2) TN_FDOF
	write(2) nwm
	write(2) MT.DA(1:TN_FDOF+1)
	write(2) MT.AA(1:nwm)

	close(unit=1); close(unit=2)
END SUBROUTINE WRITE_KM

!==========================================================
END PROGRAM FEM_PLANE_4N_MAIN

