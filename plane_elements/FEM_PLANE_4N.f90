!==============================================================================
! Four-nodes elements for plane stress & plane strain FEA
!
! 	    s
!       *
!	2 * * * 1
!	*   *   *
!	*   * * * *r
!	*       *
!	3 * * * 4
!
!	U^T = [u1,u2,...,u4, v1,v2,...,v4]
!
! CONTENTS:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!==============================================================================
MODULE FEM_PLANE_4N
USE PARAMETERS
USE FEM_PLANE_BASICFS, ONLY : INVERSE2, DET2, GAUSSLEGENDRE_QUAD
IMPLICIT NONE

PUBLIC FEM_PLANE_STIFF4N
PUBLIC MAT_H, MAT_BJ

PRIVATE STRESS_STRAIN_LAW_PLANE
PRIVATE SHAPE4N, DSHAPE4N
PRIVATE DETJACOB

CONTAINS
!==========================================================
!----------------------------------------------------------
! Calculate element stiffness matrix
!
! INPUT:
!	node(4,2) : node position
!	ptype : problem type
!			= "Plane Stress" or "Plane strain"
!	young : Young's modulus
!	poisson : Poisson's ratio
!	thick : thickness
!
! OUTPUT:
!	Ke(8,8) : element stiffness matrix
!			  [u1,u2,...,u4, v1,v2,...,v4]
!
! NEED:
!	STRESS_STRAIN_LAW_PLANE
!	GAUSSLEGENDRE_QUAD
!	MAT_BJ
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE FEM_PLANE_STIFF4N(node,ptype,young,poisson,thick,Ke)
	real(kp), intent(in) :: node(4,2)
	character(len=12), intent(in) :: ptype
	real(kp), intent(in) :: young, poisson, thick
	real(kp), intent(out) :: Ke(8,8)

	integer(ip) :: n
	parameter (n = 2)

	integer(ip) :: i, j
	real(kp) :: C(3,3)
	real(kp) :: MB(3,8)
	real(kp) :: r(n), wr(n)
	real(kp) :: detJ

	call STRESS_STRAIN_LAW_PLANE(ptype,young,poisson,C)
	call GAUSSLEGENDRE_QUAD(n,r,wr) 

	Ke = 0.0_kp
	do i = 1, n
		do j = 1, n
			call MAT_BJ(node,r(i),r(j),MB,detJ)

			Ke = Ke + matmul(transpose(MB), matmul(C,MB))*thick*&
				detJ*wr(i)*wr(j)
		end do
	end do
END SUBROUTINE FEM_PLANE_STIFF4N
!----------------------------------------------------------
! Calculate element mass matrix
!
! INPUT:
!	node(4,2) : node position
!	density : density
!	thick : thickness
!
! OUTPUT:
!	Me(8,8) : element mass matrix
!			  [u1,u2,...,u4, v1,v2,...,v4]
!
! NEED:
!	GAUSSLEGENDRE_QUAD
!	MAT_H
!	DETJACOB
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE FEM_PLANE_MASS4N(node,dens,thick,Me)
	real(kp), intent(in) :: node(4,2)
	real(kp), intent(in) :: dens
	real(kp), intent(in) :: thick
	real(kp), intent(out) :: Me(8,8)

	integer(ip) :: n
	parameter (n = 2)

	integer(ip) :: i, j
	real(kp) :: MH(2,8)
	real(kp) :: r(n), wr(n)
	real(kp) :: detJ

	call GAUSSLEGENDRE_QUAD(n,r,wr) 

	Me = 0.0_kp
	do i = 1, n
		do j = 1, n
			call MAT_H(r(i),r(j),MH)
			call DETJACOB(node,r(i),r(j),detJ)

			Me = Me + matmul(transpose(MH), MH)*thick*&
				detJ*wr(i)*wr(j)
		end do
	end do
END SUBROUTINE FEM_PLANE_MASS4N

! **********************
! * Stress-strain laws *
! **********************
!----------------------------------------------------------
! Construct stress-strain laws for plane stress & strain problems
!
! INPUT:
!	ptype : problem type 
!			= "Plane Stress" or "Plane Strain"
!	young : Young's modulus
!	poisson : Poisson's ratio
!	
! OUTPUT:
!	C(3,3) : stress-strain matrix
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, February, 2015
!----------------------------------------------------------
SUBROUTINE STRESS_STRAIN_LAW_PLANE(ptype,young,poisson,C)
	character(len=12), intent(in) :: ptype
	real(kp), intent(in) :: young, poisson
	real(kp), intent(out) :: C(3,3)

	real(kp) :: d1, d2

	select case(ptype)
	case("Plane Stress")
		d1 = young / (1.0_kp - poisson*poisson)
		C(1,:) = (/ d1, d1*poisson, 0.0_kp /)
		C(2,:) = (/ C(1,2), C(1,1), 0.0_kp /)
		C(3,:) = (/ 0.0_kp, 0.0_kp, d1*0.5_kp*(1.0_kp - poisson) /)
	case("Plane Strain")
		d1 = young / (1 + poisson)
		d2 = d1 / (1 - 2*poisson)
		C(1,:) = (/ d2*(1 - poisson), d2*poisson, 0.0_kp /)
		C(2,:) = (/ C(1,2), C(1,1), 0.0_kp /)
		C(3,:) = (/ 0.0_kp, 0.0_kp, 0.5_kp*d1 /)
	end select
END SUBROUTINE STRESS_STRAIN_LAW_PLANE

! **************************
! * Interpolation matrices *
! **************************
!----------------------------------------------------------
! Construct matrix H
!
! INPUT:
!	r, s : parametric coordinates
!
! OUTPUT:
!	MH(2,8) : interpolation matrix for displacement
!			  MH(1,:) : wrt u
!			  MH(2,:) : wrt v
!
! NEED:
!	SHAPE4N
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE MAT_H(r,s,MH)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: MH(2,8)

	real(kp) :: h(4)

	call SHAPE4N(r,s,h)

	MH = 0.0_kp
	MH(1,1:4) = h
	MH(2,5:8) = h
END SUBROUTINE MAT_H
!----------------------------------------------------------
! Calculate the determinant of Jacobian
!
! INPUT:
!	node(4,2) : nodal positions
!	r, s : parametric coordinates
!
! OUTPUT:
!	detJ : determinant of Jacobian
!
! NEED:
!	DSHAPE4N
!	DET2
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE DETJACOB(node,r,s,detJ)
	real(kp), intent(in) :: node(4,2)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: detJ

	real(kp) :: dh(4,2), J(2,2)

	call DSHAPE4N(r,s,dh)
	J = matmul(transpose(dh), node)
	call DET2(J,detJ)
END SUBROUTINE DETJACOB
!----------------------------------------------------------
! Construct matrix B and calculate determinant of Jacobian
!
! INPUT:
!	node(4,2) : nodal positions
!	r, s : parametric coordinates
!
! OUTPUT:
!	MB(3,8) : interpolation matrix for strain
!			  MB(1,:) : wrt exx
!			  MB(2,:) : wrt eyy
!			  MB(3,:) : wrt rxy
!	detJ : determinant of Jacobian
!
! NEED:
!	DSHAPE4N
!	INVERSE2
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE MAT_BJ(node,r,s,MB,detJ)
	real(kp), intent(in) :: node(4,2)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: MB(3,8)
	real(kp), intent(out) :: detJ

	real(kp) :: dh(4,2)
	real(kp) :: F(2,2), invF(2,2)

	call DSHAPE4N(r,s,dh)
	F = matmul(transpose(dh), node)
	call INVERSE2(F,invF,detJ)

	MB = 0.0_kp
	MB(1,1:4) = invF(1,1)*dh(:,1) + invF(1,2)*dh(:,2)
	MB(2,5:8) = invF(2,1)*dh(:,1) + invF(2,2)*dh(:,2)
	MB(3,1:4) = MB(2,5:8)
	MB(3,5:8) = MB(1,1:4)
END SUBROUTINE MAT_BJ

! *******************
! * Shape functions *
! *******************
!----------------------------------------------------------
! Shape function
!
! INPUT:
!	r, s : parametric coordinates
!
! OUTPUT:
!	h(4) : shape functions
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE SHAPE4N(r,s,h)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: h(4)

	real(kp) :: cr(2), cs(2)

	cr(1) = 1.0_kp - r
	cr(2) = 1.0_kp + r

	cs(1) = 1.0_kp - s
	cs(2) = 1.0_kp + s

	h(1) = 0.25_kp * cr(2) * cs(2)
	h(2) = 0.25_kp * cr(1) * cs(2)
	h(3) = 0.25_kp * cr(1) * cs(1)
	h(4) = 0.25_kp * cr(2) * cs(1)
END SUBROUTINE SHAPE4N
!----------------------------------------------------------
! Derivative of shape functions
!
! INPUT:
!	r, s : parametric coordinates
!
! OUTPUT:
!	dh(4,2) : derivative of shape functions
!			  dh(:,1) : wrt r
!			  dh(:,2) : wrt s
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE DSHAPE4N(r,s,dh)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: dh(4,2)

	real(kp) :: cr(2), cs(2)

	cr(1) = 1.0_kp - r
	cr(2) = 1.0_kp + r

	cs(1) = 1.0_kp - s
	cs(2) = 1.0_kp + s

	! wrt r
	dh(1,1) =  0.25_kp * cs(2)
	dh(2,1) = -0.25_kp * cs(2)
	dh(3,1) = -0.25_kp * cs(1)
	dh(4,1) =  0.25_kp * cs(1)

	! wrt s
	dh(1,2) =  0.25_kp * cr(2)
	dh(2,2) =  0.25_kp * cr(1)
	dh(3,2) = -0.25_kp * cr(1)
	dh(4,2) = -0.25_kp * cr(2)
END SUBROUTINE DSHAPE4N
!==========================================================
END MODULE FEM_PLANE_4N
