!==============================================================================
! ANALYSIS OF SHELL STRUCTURES BY 4-NODES SHELL ELEMENTS
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 26, June, 2016
!==============================================================================
PROGRAM FEM_SHELL_4N_MAIN
USE PARAMETERS
USE FEM_SHELL_DDTYPES
USE FEM_SHELL_MITC4N, ONLY : FEM_SHELLMITC4N_MASS, FEM_SHELLMITC4N_STIFF, CROSS3
IMPLICIT NONE

! File name
CHARACTER(20) :: FNAME

! Analysis starting time
CHARACTER(10) :: DATE(3)    ! date, time, zone
INTEGER(IP)   :: DATEVAL(8) ! yy, mm, dd, time diff wrt UTC,
							! h, m, sec, millisec

TYPE(MPROP_TYPE), ALLOCATABLE :: MATERIAL(:)
TYPE(NODE_TYPE), ALLOCATABLE  :: NODE(:)
TYPE(SHELL4N_TYPE), ALLOCATABLE :: SHELL4N(:)

! Total number of -
INTEGER(IP) :: TN_NODE       ! nodes
INTEGER(IP) :: TN_ELEMENT    ! elements
INTEGER(IP) :: TN_MSET       ! material sets

! Degrees of freedom
INTEGER(IP) :: TN_DOF    ! # of total dof
INTEGER(IP) :: TN_FDOF   ! # of free dof
INTEGER(IP) :: TN_PDOF   ! # of prescribed (Dirichlet) dof
INTEGER(IP), PARAMETER :: NDOF = 5   ! # of dof per node

! Array for stiffness and mass matrices
TYPE(SSK_TYPE) :: KT, MT

! Timer
REAL(KP) :: ITIME, ETIME

! Units for formatted writing
INTEGER(IP), parameter :: ILOG = 401

! Compute analysis starting time
call DATE_AND_TIME( date(1), date(2), date(3), dateval )

! Drive program
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
	!---------------!
	! Print heading !
	!---------------!
	call WRITE_HEADING(0)

	!-----------------!
	! Read data files !
	!-----------------!
	write(*, '("A. Reading date files")')
	call READ_DAT()    ! allocate: NODE, SHELL4N, MATERIAL

    ! Compute DOF
    call COMP_DOF()

    ! Open files for formatted writing
    call OPEN_WFILE()
    call WRITE_HEADING(ILOG)

    ! Miscellaneous computations
    call SHELL_V12()

	!---------------------------------------------!
	! Construct linear system of equations by FEM !
	!---------------------------------------------!
	write(*, '("B. Contructing linear system of equations by FEM")')

	! Allocate arrays
	write(*, '(4X, "1. Allocating arrays")')
	call ALL_KMSSK()    ! allocate: KT.DA, KT.AA, MT.DA, MT.AA
	KT.AA = 0.0_kp
	MT.AA = 0.0_kp

	call WRITE_MEMORY(0)

	! Calculate and assemble total stiffness matrix
	write(*, '(4X, "3. Calculating and assembling stiffness matrix")')
	call cpu_time(itime)

    call ASSEMBLE_KMT_FEM_SHELL4N()

	call cpu_time(etime)
	print *, "Elapsed time: ", etime - itime

	!------------------------!
	! Export data into files !
	!------------------------!
	call WRITE_KM()

    ! close files for formatted writing
    call CLOSE_WFILE()

	deallocate( NODE )
	deallocate( MATERIAL )
	deallocate( SHELL4N )
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
	integer(ip) :: FIRST_NZEROI(TN_FDOF)    ! first nonzero eqn in each column

	KT.N = TN_FDOF; MT.N = TN_FDOF
	allocate( KT.DA(TN_FDOF + 1), MT.DA(TN_FDOF + 1) )

	FIRST_NZEROI(:) = TN_FDOF
	do el = 1, TN_ELEMENT
		! Find the minimum eqn
		meqn = TN_FDOF
		do i = 1, 4
			do j = 1, NDOF
				eqn = NODE( SHELL4N(el).CONN(i) ).EQN(j)
				if ((eqn > 0) .and. (eqn < meqn)) then
					meqn = eqn
				end if
			end do
		end do

		do i = 1, 4
			do j = 1, NDOF
				eqn = node( SHELL4N(el).CONN(i) ).EQN(j)
				if ((eqn > 0) .and. (FIRST_NZEROI(eqn) > meqn)) then
					FIRST_NZEROI(eqn) = meqn
				end if
			end do
		end do
	end do

	KT.DA(1) = 1; MT.DA(1) = 1
	do i = 1, TN_FDOF
		KT.DA(i+1) = KT.DA(i) + i - FIRST_NZEROI(i) + 1
		MT.DA(i+1) = MT.DA(i) + i - FIRST_NZEROI(i) + 1
	end do

	allocate( KT.AA( KT.DA(TN_FDOF+1) - 1 ) )
	allocate( MT.AA( MT.DA(TN_FDOF+1) - 1 ) )
END SUBROUTINE ALL_KMSSK
!----------------------------------------------------------
! Calculate & assemble KT & MT arrays
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 05, July, 2016
!----------------------------------------------------------
SUBROUTINE ASSEMBLE_KMT_FEM_SHELL4N()
	integer(ip) :: n, i, j, k
	integer(ip) :: id, jd, er, ec, r, c
	integer(ip) :: address
	real(kp) :: ENODE(4,3), ETHICK(4), EV(4,3,3)
	real(kp) :: KE(20,20), ME(20,20)

	do n = 1, TN_ELEMENT
		do k = 1, 4
			ENODE(k,:) = NODE( SHELL4N(n).CONN(k) ).X(:)
            ETHICK(k)  = NODE( SHELL4N(n).CONN(k) ).THICK
            EV(k,:,:)  = NODE( SHELL4N(n).CONN(k) ).V(:,:)
		end do

        call FEM_SHELLMITC4N_STIFF(ENODE,ETHICK,EV,&
                                   MATERIAL( SHELL4N(n).MPROP ).YOUNG,&
                                   MATERIAL( SHELL4N(n).MPROP ).POISSON,KE)
        call FEM_SHELLMITC4N_MASS(ENODE,ETHICK,EV,&
                                  MATERIAL( SHELL4N(n).MPROP ).DENS,ME)

        call WRITE_KME(ILOG,n,KE,ME)

		do i = 1, 4
			do id = 1, NDOF
				er = (id-1)*4 + i
				r  = NODE( SHELL4N(n).CONN(i) ).EQN(id)
                if (r > 0) then
				    do j = 1, 4
				    	do jd = 1, NDOF
				    		ec = (jd-1)*4 + j
				    		c  = NODE( SHELL4N(n).CONN(j) ).EQN(jd)

				    		if ( c >= r ) then
				    			address = KT.DA(c) + c - r
				    			KT.AA(address) = KT.AA(address) + KE(er,ec)
				    			address = MT.DA(c) + c - r
				    			MT.AA(address) = MT.AA(address) + ME(er,ec)
				    		end if
				    	end do
				    end do
                end if
			end do
		end do
	end do
END SUBROUTINE ASSEMBLE_KMT_FEM_SHELL4N

! ******************************
! * Miscellaneous computations *
! ******************************
!----------------------------------------------------------
! Compute degrees of freedom
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!----------------------------------------------------------
SUBROUTINE COMP_DOF()
    integer(ip) :: i, j

    TN_DOF  = NDOF*TN_NODE
    TN_FDOF = 0
    do i = 1, TN_NODE
        do j = 1, 6
            if (NODE(i).EQN(j) > 0) then
                TN_FDOF = TN_FDOF + 1

                if (j == 6) then
                    TN_DOF = TN_DOF + 1
                end if
            end if
        end do
    end do

    TN_PDOF = TN_DOF - TN_FDOF
END SUBROUTINE COMP_DOF

!----------------------------------------------------------
! Compute V1 and V2 in each node
!
! NEED:
!   CROSS3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 05, July, 2016
!----------------------------------------------------------
SUBROUTINE SHELL_V12()
    integer(ip) :: i, k
    real(kp), parameter :: tol = 1.0E-10_kp
    real(kp), parameter :: EY(3) = (/0.0_kp, 1.0_kp, 0.0_kp/)
    real(kp) :: ERR(3)


    do i = 1, TN_NODE
        if (NODE(i).V(1,1) == 0.0_kp .and.&
            NODE(i).V(1,2) == 0.0_kp .and.&
            NODE(i).V(1,3) == 0.0_kp) then
            cycle
        end if

        do k = 1, 3
            ERR(k) = abs(abs(NODE(i).V(1,k)) - EY(k))
        end do

        if ( ERR(1) < tol .and. ERR(2) < tol .and. ERR(3) < tol ) then
            NODE(i).V(2,:) = (/0.0_kp, 0.0_kp, 1.0_kp/)
        else
            call CROSS3(EY,NODE(i).V(1,:),NODE(i).V(2,:))
            NODE(i).V(2,:) = NODE(i).V(2,:) / NORM2(NODE(i).V(2,:))
        end if
        call CROSS3(NODE(i).V(1,:),NODE(i).V(2,:),NODE(i).V(3,:))
    end do
END SUBROUTINE SHELL_V12

! ********************
! * Reading routines *
! ********************
!----------------------------------------------------------
! Read input discretization data
!
! ALLOCATE:
!   NODE, SHELL4N, MATERIAL
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 26, June, 2016
!----------------------------------------------------------
SUBROUTINE READ_DAT()
    integer(ip) :: i, dumi
    character(40) :: dumc

	!-----------------!
	! Input file name !
	!-----------------!
	write(*, '(4X, "Input the file name", &
		1X, "(without extension)")')
	write(*, '(8X, ">")', advance='no')
	read(*,*) FNAME

	!------------!
	! Open files !
	!------------!
	open(unit=1, file=trim(FNAME)//".dis", status='old', &
		action='read', form='formatted')

    !-----------!
    ! Read data !
    !-----------!
    do
        read(1,'(A14)') dumc(1:14)
        if (dumc(1:14) == "* # OF NODES *") then
            exit
        end if
    end do
    read(1,*) TN_NODE
    allocate( NODE(TN_NODE) )

    do
        read(1,'(A34)') dumc(1:34)
        if (dumc(1:34) == "* NODE# X Y Z Vx Vy Vz EQN THICK *") then
            exit
        end if
    end do

    do i = 1, TN_NODE
        read(1,*) dumi, NODE(i).X
        read(1,*) NODE(i).V(1,:)
        read(1,*) NODE(i).EQN
        read(1,*) NODE(i).THICK
    end do

    do
        read(1,'(A30)') dumc(1:30)
        if (dumc(1:30) == "* # OF MATERIAL PROPERTY SET *") then
            exit
        end if
    end do
    read(1,*) TN_MSET
    allocate( MATERIAL(TN_MSET) )

    do
        read(1,'(A30)') dumc(1:30)
        if (dumc(1:30) == "* SET# YOUNG POISSON DENSITY *") then
            exit
        end if
    end do

    do i = 1, TN_MSET
        read(1,*) dumi, MATERIAL(i).YOUNG, MATERIAL(i).POISSON, MATERIAL(i).DENS
    end do

    do
        read(1,'(A17)') dumc(1:17)
        if (dumc(1:17) == "* # OF ELEMENTS *") then
            exit
        end if
    end do
    read(1,*) TN_ELEMENT
    allocate( SHELL4N(TN_ELEMENT) )

    do
        read(1, '(A35)') dumc(1:35)
        if (dumc(1:35) == "* ELEMENT# CONNECTIVITY MPROPSET# *") then
            exit
        end if
    end do

    do i = 1, TN_ELEMENT
        read(1,*) dumi, SHELL4N(i).CONN, SHELL4N(i).MPROP
    end do
   
    !-------------!
    ! Close files !
    !-------------!
    close(unit=1)
END SUBROUTINE READ_DAT

! ********************
! * Writing routines *
! ********************
!----------------------------------------------------------
! Open files to write results
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 07, June, 2016
!----------------------------------------------------------
SUBROUTINE OPEN_WFILE()
    open(unit=ILOG, file=trim(FNAME)//".log", status='replace', &
		 form='formatted', action='write', position='rewind')
END SUBROUTINE OPEN_WFILE
!----------------------------------------------------------
! Close files
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 07, June, 2016
!----------------------------------------------------------
SUBROUTINE CLOSE_WFILE()
    close(unit=ILOG)
END SUBROUTINE CLOSE_WFILE
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

	write(n, '( 8X, 57("*") )')
	write(n, '( 8X, "*", 55X, "*" )')
	write(n, '( 8X, "* ANALYSIS OF SHELL STRUCTURES BY MITC4 SHELL ELEMENTS",&
		2X, "*")')
	write(n, '( 8X, "*", 55X, "*" )')
	write(n, '( 8X, "* PROGRAM AUTHOR: KI-TAE KIM", 28X, "*" )')
	write(n, '( 8X, "* AFFILIATION   :", X, &
        "DEPARTMENT OF MECHANICAL ENGINEERING,", X, "*" )')
	write(n, '( 8X, "*", 17X, "MASSACHUSETTS INSTITUTE OF TECHNOLOGY", X,&
		"*" )')
	write(n, '( 8X, "* E-MAIL        : ktkim@mit.edu or qlsn55@gmail.com",&
		        5X, "*" )')
	write(n, '( 8X, "*", 55X, "*" )')
	Write(n, '( 8X, "* TIME OF ANALYSIS: ", I2, ":", I2, ", ", I2, & 
		"/", I2, "/", I4, 18X, " *" )' ) dateval(5), dateval(6), &
		dateval(2), dateval(3), dateval(1)
	write(n, '( 8X, "*", 55X, "*" )')
	write(n, '( 8X, 57("*") )')
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

	write(n, '(8X, "---------------------")')
	write(n, '(8X, "TN_DOF  : ", I7)') TN_DOF
	write(n, '(8X, "TN_FDOF : ", I7)') TN_FDOF
	write(n, '(8X, "TN_PDOF : ", I7)') TN_PDOF
	write(n, '()')
	write(n, '(8X, "Allocated memory size")')
	write(n, 10) int(size(KT.AA) * 8.0E-3_kp +&
				 size(KT.DA) * 4.0E-3_kp)
	write(n, 11) int(size(MT.AA) * 8.0E-3_kp +&
				 size(MT.DA) * 4.0E-3_kp)
!	write(n, 20) int(size(RT) * 8.0E-3_kp)
!	write(n, 30) int(size(UT) * 8.0E-3_kp)
	10 format (8X, "KT :", 1X, I8, "KB")
	11 format (8X, "MT :", 1X, I8, "KB")
!	20 format (8X, "RT :", 1X, I8, "KB")
!	30 format (8X, "UT :", 1X, I8, "KB")
	write(n, '(8X, "---------------------")')
END SUBROUTINE WRITE_MEMORY
!----------------------------------------------------------
! Export stiffness and mass arrays to binary files
! - save as:
!	TN_FDOF = *.N
!	nnzero = *.DA( *.N + 1 ) - 1
!	*.DA(1:TN_FDOF)
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
!----------------------------------------------------------
! Write element stiffness and mass matrices
!
! INPUT:
!   n         = unit specifier in write statement
!   el        = element #
!   KE(20,20) = element stiffness matrix
!   ME(20,20) = element mass matrix
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!----------------------------------------------------------
SUBROUTINE WRITE_KME(n,el,KE,ME)
    integer(ip), intent(in) :: n, el
    real(kp), intent(in)    :: KE(20,20), ME(20,20)

    integer(ip) :: i

    write(n,100) el
    write(n,101)
    do i = 1, 20
        write(n,110) KE(i,:)
    end do

    write(n,102)
    do i = 1, 20
        write(n,110) ME(i,:)
    end do

    100 format(/," ELEMENT ", I5)
    101 format(/," STIFFNESS MATRIX:")
    102 format(/," MASS MATRIX:")
    110 format(20(X, ES22.15))
    111 format(20(X, ES11.4))
END SUBROUTINE WRITE_KME

!==========================================================
END PROGRAM FEM_SHELL_4N_MAIN

