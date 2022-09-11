!==============================================================================
! Nine-nodes elements for plane stress & plane strain FEA
!
! 	    s
!       *
!	2 * 5 * 1
!	*   *   *
!	6   9 * 8 *r
!	*       *
!	3 * 7 * 4
!
!	U^T = [u1,u2,...,u9, v1,v2,...,v9]
!
! CONTENTS:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 1982
!==============================================================================
MODULE FEM_PLANE_9N
USE PARAMETERS
USE FEM_PLANE_BASICFS, ONLY : INVERSE2, GAUSSLEGENDRE_QUAD
IMPLICIT NONE

PUBLIC FEM_PLANE_STIFF9N
PUBLIC FEM_PLANE_USLOAD3N
PUBLIC MAT_H, MAT_BJ

PRIVATE STRESS_STRAIN_LAW_PLANE
PRIVATE SHAPE9N, DSHAPE9N
PRIVATE SHAPE3N

CONTAINS
!==========================================================
!----------------------------------------------------------
! Calculate element stiffness matrix
!
! INPUT:
!	node(9,2) : node position
!	ptype : problem type
!			= "Plane Stress" or "Plane strain"
!	young : Young's modulus
!	poisson : Poisson's ratio
!	thick : thickness
!
! OUTPUT:
!	Ke(18,18) : element stiffness matrix
!				[u1,u2,...,u9, v1,v2,...,v9]
!
! NEED:
!	STRESS_STRAIN_LAW_PLANE
!	GAUSSLEGENDRE_QUAD
!	MAT_BJ
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE FEM_PLANE_STIFF9N(node,ptype,young,poisson,thick,Ke)
	real(kp), intent(in) :: node(9,2)
	character(len=12), intent(in) :: ptype
	real(kp), intent(in) :: young, poisson, thick
	real(kp), intent(out) :: Ke(18,18)

	integer(i4b) :: n
	parameter (n = 3)

	integer(i4b) :: i, j
	real(kp) :: C(3,3)
	real(kp) :: MB(3,18)
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
END SUBROUTINE FEM_PLANE_STIFF9N

!----------------------------------------------------------
! Calculate element uniform surface load
!
! INPUT:
!	bnode(2,2) : boundary nodal coordinates
!	fs(2) : uniform surface load
!			[px, py]
!	thick : thickness
!
! OUTPUT:
!	Re(6)
!	
! NEED:
!	GAUSSLEGENDRE_QUAD
!	SHAPE3N
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 25, March, 2015
!----------------------------------------------------------
SUBROUTINE FEM_PLANE_USLOAD3N(bnode,fs,thick,Re)
	real(kp), intent(in) :: bnode(2,2)
	real(kp), intent(in) :: fs(2)
	real(kp), intent(in) :: thick
	real(kp), intent(out) :: Re(6)

	integer(i4b) :: n
	parameter (n = 3)

	integer(i4b) :: i, j, k, m
	real(kp) :: r(n), wr(n)
	real(kp) :: h(3)
	real(kp) :: l(2), detJ

	l = bnode(2,:) - bnode(1,:)
	detJ = 0.5_kp * sqrt( l(1)*l(1) + l(2)*l(2) )

	call GAUSSLEGENDRE_QUAD(n,r,wr) 

	Re = 0.0_kp
	do i = 1, n
		call SHAPE3N(r(i),h)

		do j = 1, 2
			do k = 1, 3
				m = (j-1)*3 + k
				Re(m) = Re(m) + h(k)*fs(j)*thick*detJ*wr(i)
			end do
		end do
	end do
END SUBROUTINE FEM_PLANE_USLOAD3N

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
!	MH(2,18) : interpolation matrix for displacement
!			   MH(1,:) : wrt u
!			   MH(2,:) : wrt v
!
! NEED:
!	SHAPE9N
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE MAT_H(r,s,MH)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: MH(2,18)

	real(kp) :: h(9)

	call SHAPE9N(r,s,h)

	MH = 0.0_kp
	MH(1,1:9)   = h
	MH(2,10:18) = h
END SUBROUTINE MAT_H
!----------------------------------------------------------
! Construct matrix B and calculate determinant of Jacobian
!
! INPUT:
!	node(9,2) : nodal positions
!	r, s : parametric coordinates
!
! OUTPUT:
!	MB(3,18) : interpolation matrix for strain
!			   MB(1,:) : wrt exx
!			   MB(2,:) : wrt eyy
!			   MB(3,:) : wrt rxy
!	detJ : determinant of Jacobian
!
! NEED:
!	DSHAPE9N
!	INVERSE2
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE MAT_BJ(node,r,s,MB,detJ)
	real(kp), intent(in) :: node(9,2)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: MB(3,18)
	real(kp), intent(out) :: detJ

	integer(i4b) :: i, j, k
	real(kp) :: dh(9,2)
	real(kp) :: F(2,2), invF(2,2)

	call DSHAPE9N(r,s,dh)

	F = 0.0_kp
	do i = 1, 2
		do j = 1, 2
			do k = 1, 9
				F(i,j) = F(i,j) + dh(k,i) * node(k,j)
			end do
		end do
	end do

	call INVERSE2(F,invF,detJ)

	MB = 0.0_kp
	MB(1,1:9)   = invF(1,1)*dh(:,1) + invF(1,2)*dh(:,2)
	MB(2,10:18) = invF(2,1)*dh(:,1) + invF(2,2)*dh(:,2)
	MB(3,1:9)   = MB(2,10:18)
	MB(3,10:18) = MB(1,1:9)
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
!	h(9) : shape functions
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE SHAPE9N(r,s,h)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: h(9)

	real(kp) :: cr(3), cs(3), ch

	cr(1) = 1.0_kp - r
	cr(2) = 1.0_kp + r
	cr(3) = 1.0_kp - r*r

	cs(1) = 1.0_kp - s
	cs(2) = 1.0_kp + s
	cs(3) = 1.0_kp - s*s

	h(9) = cr(3) * cs(3)
	ch   = 0.5_kp*h(9)

	h(8) = 0.5_kp*cs(3)*cr(2) - ch
	h(7) = 0.5_kp*cr(3)*cs(1) - ch
	h(6) = 0.5_kp*cs(3)*cr(1) - ch
	h(5) = 0.5_kp*cr(3)*cs(2) - ch

	ch   = 0.25_kp*h(9)

	h(4) = 0.25_kp*cr(2)*cs(1) - 0.5_kp*(h(7) + h(8)) - ch
	h(3) = 0.25_kp*cr(1)*cs(1) - 0.5_kp*(h(6) + h(7)) - ch
	h(2) = 0.25_kp*cr(1)*cs(2) - 0.5_kp*(h(5) + h(6)) - ch
	h(1) = 0.25_kp*cr(2)*cs(2) - 0.5_kp*(h(5) + h(8)) - ch
END SUBROUTINE SHAPE9N
!----------------------------------------------------------
! Derivative of shape functions
!
! INPUT:
!	r, s : parametric coordinates
!
! OUTPUT:
!	dh(9,2) : derivative of shape functions
!			  dh(:,1) : wrt r
!			  dh(:,2) : wrt s
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE DSHAPE9N(r,s,dh)
	real(kp), intent(in) :: r, s
	real(kp), intent(out) :: dh(9,2)

	real(kp) :: cr(3), cs(3), ch

	cr(1) = 1.0_kp - r
	cr(2) = 1.0_kp + r
	cr(3) = 1.0_kp - r*r

	cs(1) = 1.0_kp - s
	cs(2) = 1.0_kp + s
	cs(3) = 1.0_kp - s*s

	! wrt r
	dh(9,1) = -2.0_kp*r*cs(3)
	ch = 0.5_kp*dh(9,1)

	dh(8,1) =  0.5_kp*cs(3) - ch
	dh(7,1) = -r*cs(1) - ch
	dh(6,1) = -0.5_kp*cs(3) - ch
	dh(5,1) = -r*cs(2) - ch

	ch = 0.25_kp*dh(9,1)

	dh(4,1) =  0.25_kp*cs(1) - 0.5_kp*(dh(7,1) + dh(8,1)) - ch
	dh(3,1) = -0.25_kp*cs(1) - 0.5_kp*(dh(6,1) + dh(7,1)) - ch
	dh(2,1) = -0.25_kp*cs(2) - 0.5_kp*(dh(5,1) + dh(6,1)) - ch
	dh(1,1) =  0.25_kp*cs(2) - 0.5_kp*(dh(5,1) + dh(8,1)) - ch

	! wrt s
	dh(9,2) = -2.0_kp*s*cr(3)
	ch = 0.5_kp*dh(9,2)

	dh(8,2) = -s*cr(2) - ch
	dh(7,2) = -0.5_kp*cr(3) - ch
	dh(6,2) = -s*cr(1) - ch
	dh(5,2) =  0.5_kp*cr(3) - ch

	ch = 0.25_kp*dh(9,2)

	dh(4,2) = -0.25_kp*cr(2) - 0.5_kp*(dh(7,2) + dh(8,2)) - ch
	dh(3,2) = -0.25_kp*cr(1) - 0.5_kp*(dh(6,2) + dh(7,2)) - ch
	dh(2,2) =  0.25_kp*cr(1) - 0.5_kp*(dh(5,2) + dh(6,2)) - ch
	dh(1,2) =  0.25_kp*cr(2) - 0.5_kp*(dh(5,2) + dh(8,2)) - ch
END SUBROUTINE DSHAPE9N

!----------------------------------------------------------
! Shape function
!
!       * * r
!       *
!	2 * 3 * 1
!
! INPUT:
!	r : parametric coordinates
!
! OUTPUT:
!	h(3) : shape functions
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!----------------------------------------------------------
SUBROUTINE SHAPE3N(r,h)
	real(kp), intent(in) :: r
	real(kp), intent(out) :: h(3)

	real(kp) :: ch

	h(3) = 1.0_kp - r*r
	ch = 0.5_kp*h(3)

	h(2) = 0.5_kp*(1.0_kp + r) - ch
	h(1) = 0.5_kp*(1.0_kp - r) - ch
END SUBROUTINE SHAPE3N

!==========================================================
END MODULE FEM_PLANE_9N

