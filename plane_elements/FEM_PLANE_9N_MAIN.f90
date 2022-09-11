!==============================================================================
! FEA of plane stress & plane strain problems
! - 2D 9-nodes elements
!
! CONTENTS:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
PROGRAM FEM_PLANE_9N_MAIN
USE PARAMETERS
USE FEM_PLANE_DDTYPES
USE FEM_PLANE_9N, ONLY : FEM_PLANE_STIFF9N, FEM_PLANE_USLOAD3N, &
	MAT_H
USE FACSPAR, ONLY : LDLTSSK
IMPLICIT NONE

! File name
CHARACTER(20) :: FNAME

! Analysis starting time
CHARACTER(10) :: DATE(3)    ! date, time, zone
INTEGER(I4B)  :: DATEVAL(8) ! yy, mm, dd, time diff wrt UTC,
							! h, m, sec, millisec

! Problem type
CHARACTER(12) :: PTYPE

TYPE(MPROP_TYPE), ALLOCATABLE :: MATERIAL(:)
TYPE(NODE_TYPE), ALLOCATABLE :: NODE(:)
TYPE(QUAD9N_TYPE), ALLOCATABLE :: QUAD9N(:)
TYPE(NBOUNDARY3N_TYPE), ALLOCATABLE :: LINE3N(:)

! Total number of -
INTEGER(I4B) :: TN_NODE    ! nodes
INTEGER(I4B) :: TN_ELEMENT    ! elements
INTEGER(I4B) :: TN_MSET    ! material sets
INTEGER(I4B) :: TN_NLINE    ! Neumann boundary (uniform surface load) lines

! Degrees of freedom
INTEGER(I4B) :: TN_DOF    ! # of total dof
INTEGER(I4B) :: TN_FDOF   ! # of free dof
INTEGER(I4B) :: TN_PDOF   ! # of prescribed (Dirichlet) dof
INTEGER(I4B) :: NDOF   ! # of dof per node
PARAMETER( NDOF = 2 )

! Array for stiffness matrix
TYPE(SSK_TYPE) :: KT

! Array for load vector
REAL(KP), ALLOCATABLE :: RT(:)

! Array for unkown diaplacement vector
REAL(KP), ALLOCATABLE :: UT(:)

! Strain energy
REAL(KP) :: STRAIN_ENERGY

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
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE MAIN_FEM()
	integer(i4b) :: i

	!--------------
	! Print heading
	!--------------
	call WRITE_HEADING(0)

	!----------------
	! Read data files
	!----------------
	write(*, '("A. Reading date files")')
	call READ_INPUTDATA()    ! allocate: NODE, MATERIAL, QUAD9N
	call READ_NEUMANN()    ! allocate: LINE3N

	!--------------------------------------------
	! Construct linear system of equations by MFR
	!--------------------------------------------
	write(*, '("B. Contructing linear system of equations by FEM")')

	! Set equation number
	write(*, '(4X, "1. Setting equation numbers")')
	call SET_EQN()

	! Allocate arrays
	write(*, '(4X, "2. Allocating arrays")')
	call ALL_KSSK()    ! allocate: KT.DA, KT.AA
	KT.AA = 0.0_kp

	allocate( RT(TN_DOF) )    ! allocate: RT
	RT = 0.0_kp

	allocate( UT(TN_DOF) )    ! allocate: UT
	UT = 0.0_kp

	call WRITE_MEMORY(0)

	! Calculate and assemble total stiffness matrix
	write(*, '(4X, "3. Calculating and assembling stiffness matrix")')

	call cpu_time(itime)

	call ASSEMBLE_KT_FEM_PLANE()
	
	! Calculate and assemble unifrom surface load
	write(*, '(4X, "4. Calculating and assembling load vector")')
	call ASSEMBLE_RT_UNIFROM_FEM_PLANE()

	call cpu_time(etime)
	print *, "Elapsed time: ", etime - itime

	!--------------------------------
	! Solve linear system of equation
	!--------------------------------
	write(*, '("C. Sonvilng linear system of equation")')

	call cpu_time(itime)

	CALL SOLVE_LSE_DIRECT()

	call cpu_time(etime)
	print *, "Elapsed time: ", etime - itime

	write(*, '(8X, "STRAIN ENERGY = ", ES13.6)') STRAIN_ENERGY

	!---------------------------------
	! Export data into formatted files
	!---------------------------------
	! Equation number
!	call WRITE_EQN()

	! Displacements computed
!	call WRITE_U(11,11)

	deallocate( NODE )
	deallocate( MATERIAL )
	deallocate( QUAD9N )
	deallocate( LINE3N )
	deallocate( KT.DA, KT.AA )
	deallocate( RT )
	deallocate( UT )
END SUBROUTINE MAIN_FEM

! ****************************************
! * Construct linear system of equations *
! ****************************************
!----------------------------------------------------------
! Allocate DA and AA arrays (SSK format) of stiffness matrix K
!
! ALLOCATED VARIABLES:
!	KT.DA
!	KT.AA
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 10, February, 2015
!----------------------------------------------------------
SUBROUTINE ALL_KSSK()
	integer(i4b) :: el, i, j
	integer(i4b) :: eqn, meqn
	integer(i4b) :: first_nzeroi(TN_DOF)    ! first nonzero eqn in each column

	KT.N = TN_DOF
	allocate( KT.DA(TN_DOF + 1) )

	first_nzeroi(:) = TN_DOF
	do el = 1, TN_ELEMENT
		! Find the minimum eqn
		meqn = TN_DOF
		do i = 1, 9
			do j = 1, NDOF
				eqn = NODE( QUAD9N(el).CONN(i) ).EQN(j)
				if ( eqn < meqn ) then
					meqn = eqn
				end if
			end do
		end do

		do i = 1, 9
			do j = 1, NDOF
				eqn = node( QUAD9N(el).CONN(i) ).EQN(j)
				if ( first_nzeroi(eqn) > meqn ) then
					first_nzeroi(eqn) = meqn
				end if
			end do
		end do
	end do

	KT.DA(1) = 1
	do i = 1, TN_DOF
		KT.DA(i+1) = KT.DA(i) + i - first_nzeroi(i) + 1
	end do

	allocate( KT.AA( KT.DA(TN_DOF+1) - 1 ) )
END SUBROUTINE ALL_KSSK
!----------------------------------------------------------
! Calculate & assemble KT matrix 
!
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!----------------------------------------------------------
SUBROUTINE ASSEMBLE_KT_FEM_PLANE()
	integer(i4b) :: n, i, j, k
	integer(i4b) :: id, jd, er, ec, r, c
	integer(i4b) :: address
	real(kp) :: enode(9,2)
	real(kp) :: Ke(18,18)

	do n = 1, TN_ELEMENT
		do k = 1, 9
			enode(k,:) = NODE( QUAD9N(n).CONN(k) ).X(:)
		end do

		call FEM_PLANE_STIFF9N(enode,PTYPE,&
			MATERIAL( QUAD9N(n).MPROP ).YOUNG,&
		 	MATERIAL( QUAD9N(n).MPROP ).POISSON,&
		 	MATERIAL( QUAD9N(n).MPROP ).THICK,&
			Ke)

		do i = 1, 9
			do id = 1, NDOF
				er = (id-1)*9 + i
				r  = NODE( QUAD9N(n).CONN(i) ).EQN(id)
				do j = 1, 9
					do jd = 1, NDOF
						ec = (jd-1)*9 + j
						c  = NODE( QUAD9N(n).CONN(j) ).EQN(jd)

						if ( c >= r ) then
							address = KT.DA(c) + c - r
							KT.AA(address) = KT.AA(address) + Ke(er,ec)
						end if
					end do
				end do
			end do
		end do
	end do
END SUBROUTINE ASSEMBLE_KT_FEM_PLANE
!----------------------------------------------------------
! Calculate & assemble uniform surface load vector
!
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
!----------------------------------------------------------
SUBROUTINE ASSEMBLE_RT_UNIFROM_FEM_PLANE()
	integer(i4b) :: n, k
	integer(i4b) :: j, jd, ec, c
	real(kp) :: bnode(2,2)
	real(kp) :: fs(2)
	real(kp) :: Re(6)

	do n = 1, TN_NLINE
		do k = 1, 2
			bnode(k,:) = NODE( LINE3N(n).CONN(k) ).X(:)
		end do
		fs = LINE3N(n).UF

		call FEM_PLANE_USLOAD3N(bnode,fs,&
			MATERIAL( QUAD9N(LINE3N(n).EN).MPROP ).THICK,&
			Re)

		do j = 1, 3
			do jd = 1, NDOF
				ec = (jd-1)*3 + j
				c  = NODE( LINE3N(n).CONN(j) ).EQN(jd)

				RT(c) = RT(c) + Re(ec)
			end do
		end do
	end do
END SUBROUTINE ASSEMBLE_RT_UNIFROM_FEM_PLANE

! ********************
! * Solving routines *
! ********************
!----------------------------------------------------------
! Solve the linear system of equations by direct method
! and calculate strain energy
! - "KT*UT = RT"
! - LDL^T method
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
!----------------------------------------------------------
SUBROUTINE SOLVE_LSE_DIRECT()
	integer(i4b) :: na

	na = KT.DA(TN_FDOF+1) - 1

	! Copy
	UT(1:TN_FDOF) = RT(1:TN_FDOF)

	! LDL^T decomposition of KT
	call LDLTSSK(0, TN_FDOF, KT.AA(1:na), KT.DA(1:TN_FDOF+1), UT(1:TN_FDOF))

	! Solve
	call LDLTSSK(1, TN_FDOF, KT.AA(1:na), KT.DA(1:TN_FDOF+1), UT(1:TN_FDOF))

	! Calculate strain energy
	STRAIN_ENERGY = 0.5_kp * dot_product(UT, RT)
END SUBROUTINE SOLVE_LSE_DIRECT

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
	integer(i4b) :: i, j

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
	integer(i4b) :: i, n

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
	read(2,'(A12)') PTYPE

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
		read(2,*) n, MATERIAL(n).YOUNG, MATERIAL(n).POISSON, MATERIAL(n).THICK
	end do

	!------------------
	! Read element data
	!------------------
	do i = 1, 7
		read(2,*)
	end do

	! Total # of elements
	read(2,*) TN_ELEMENT

	allocate( QUAD9N(TN_ELEMENT) )

	do i = 1, 5
		read(2,*)
	end do

	! Read
	do i = 1, TN_ELEMENT
		read(2,*) n, QUAD9N(n).CONN(:), QUAD9N(n).MPROP
	end do
	
	!------------
	! Close files
	!------------
	close(unit=1)
	close(unit=2)
END SUBROUTINE READ_INPUTDATA

!----------------------------------------------------------
! Read Neumann boundary condition data from formatted files 
!
! ALLOCATED VARIABLES:
!	LINE3N	
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
!----------------------------------------------------------
SUBROUTINE READ_NEUMANN()
	integer(i4b) :: i, n

	open(unit=1, file=trim(FNAME)//".neumann", status='old', &
		action='read', form='formatted')

	do i = 1, 3
		read(1,*)
	end do

	! Total # of Neumann boundary lines (uniform load)
	read(1,*) TN_NLINE

	allocate( LINE3N(TN_NLINE) )

	do i = 1, 5
		read(1,*)
	end do

	! Read
	do i = 1, TN_NLINE
		read(1,*) n, LINE3N(n).EN, LINE3N(n).CONN(:), LINE3N(n).UF(:)
	end do

	close(unit=1)
END SUBROUTINE READ_NEUMANN

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
	integer(i4b), intent(in) :: n

	write(n, '( 8X, 61("*") )')
	write(n, '( 8X, "*", 59X, "*" )')
	write(n, '( 8X, "* Plane Stress & Strain Analysis by 9-nodes ", &
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
! Write equation number
!	File extension : eqn
!	Format : [(node #), eqns for Ux, eqns for Uy]
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, March, 2015
!----------------------------------------------------------
SUBROUTINE WRITE_EQN()
	integer(i4b) :: i, j, k

	open(unit=1, file=trim(FNAME)//".eqn", status='replace', &
		form='formatted', action='write', position='rewind')

	do i = 1, TN_NODE
		do j = 1, 2
			write(1,10,advance='no') NODE(i).EQN(j)
		end do
		write(1,11)
	end do

	close(unit=1)

	10 format(I6, 1X)
	11 format()
END SUBROUTINE WRITE_EQN
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
	integer(i4b), intent(in) :: n

	write(*, '(8X, "---------------------")')
	write(*, '(8X, "TN_DOF  : ", I7)') TN_DOF
	write(*, '(8X, "TN_FDOF : ", I7)') TN_FDOF
	write(*, '(8X, "TN_PDOF : ", I7)') TN_PDOF
	write(*, '()')
	write(*, '(8X, "Allocated memory size")')
	write(*, 10) int(size(KT.AA) * 8.0E-3_kp +&
				 size(KT.DA) * 4.0E-3_kp)
	write(*, 20) int(size(RT) * 8.0E-3_kp)
	write(*, 30) int(size(UT) * 8.0E-3_kp)
	10 format (8X, "KT :", 1X, I7, "KB")
	20 format (8X, "RT :", 1X, I7, "KB")
	30 format (8X, "UT :", 1X, I7, "KB")
	write(*, '(8X, "---------------------")')
END SUBROUTINE WRITE_MEMORY
!----------------------------------------------------------
! Write displacement computed
!	File extension : ugridx, ugridy, ux, uy
!
! INPUT:
!	m : # of points in r direction per element
!	n : # of points in s direction per element
!
! OUTPUT:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!----------------------------------------------------------
SUBROUTINE WRITE_U(m,n)
	integer(i4b), intent(in) :: m, n

	integer(i4b) :: en, i, j
	real(kp) :: dr, ds, r, s
	real(kp) :: x, y, u, v
	real(kp) :: MH(2,18)
	real(kp) :: Xe(9), Ye(9), Ue(9), Ve(9)

	open(unit=1, file=trim(FNAME)//".ugridx", status='replace', &
		form='formatted', action='write', position='rewind')
	open(unit=2, file=trim(FNAME)//".ugridy", status='replace', &
		form='formatted', action='write', position='rewind')
	open(unit=3, file=trim(FNAME)//".ux", status='replace', &
		form='formatted', action='write', position='rewind')
	open(unit=4, file=trim(FNAME)//".uy", status='replace', &
		form='formatted', action='write', position='rewind')

	write(1,'(I7, 1X, I7, 1X, I7)') TN_ELEMENT, m, n
	write(2,'(I7, 1X, I7, 1X, I7)') TN_ELEMENT, m, n

	dr = 2.0_kp / (m - 1)
	ds = 2.0_kp / (n - 1)
	do en = 1, TN_ELEMENT
		do i = 1, 9
			Xe(i) = NODE( QUAD9N(en).CONN(i) ).X(1)
			Ye(i) = NODE( QUAD9N(en).CONN(i) ).X(2)
			Ue(i) = UT( NODE( QUAD9N(en).CONN(i) ).EQN(1) )
			Ve(i) = UT( NODE( QUAD9N(en).CONN(i) ).EQN(2) )
		end do

		do j = 1, n 
			s = -1.0_kp + ds*(j - 1)
			do i = 1, m
				s = -1.0_kp + ds*(j - 1)
				r = -1.0_kp + dr*(i - 1)

				call MAT_H(r,s,MH)
				x = dot_product( MH(1,1:9), Xe )
				y = dot_product( MH(2,10:18), Ye )
				u = dot_product( MH(1,1:9), Ue )
				v = dot_product( MH(2,10:18), Ve )
				
				write(1,100,advance='no') x
				write(2,100,advance='no') y
				write(3,100,advance='no') u
				write(4,100,advance='no') v
			end do

			write(1,200)
			write(2,200)
			write(3,200)
			write(4,200)
		end do
	end do

	close(unit=1)
	close(unit=2)
	close(unit=3)
	close(unit=4)

	100 format(ES13.6, 1X)
	200 format()
END SUBROUTINE WRITE_U
!==========================================================
END PROGRAM FEM_PLANE_9N_MAIN

