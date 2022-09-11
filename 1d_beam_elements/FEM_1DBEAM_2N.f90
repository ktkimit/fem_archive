!==============================================================================
! 1D two-nodes beam elements
!
!		*
!		*
!	1 * * * 2 *r
!
!	U^T = [w1, w2, t1, t2]
!
! CONTENTS:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 11, May, 2015
!==============================================================================
MODULE FEM_1DBEAM_2N
USE PARAMETERS
IMPLICIT NONE

PUBLIC FEM_1DBEAM_STIFF2N, FEM_1DBEAM_MASS2N

PRIVATE SHEAR_BMH, SHAPE2N, GAUSSLEGENDRE_QUAD

CONTAINS
!==========================================================
!----------------------------------------------------------
! Calculate element stiffness matrix
!
! INPUT:
!	node(2) : node position
!	etype : element type 
!		   "DISP" = displacement based
!		   "MITC" = mixed interpolation
!	EI  : E*I where E = Young's modulus and I = second moment of 
!		  sectional area
!	GAk : G*A*k where G = shear modulus, A = sectional area and 
!		  k = Timoshenko shear coefficient
!	
! OUTPUT:
!	Ke(4,4) : element stiffness matrix
!
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 13, April, 2016
!----------------------------------------------------------
SUBROUTINE FEM_1DBEAM_STIFF2N(node, etype, EI, GAk, Ke)
	real(kp), intent(in) :: node(2)
	character(len=4), intent(in) :: etype
	real(kp), intent(in) :: EI, GAk
	real(kp), intent(out) :: Ke(4,4)

	integer(i4b) :: i
	real(kp) :: dx
	real(kp) :: c1, c2, c3
	real(kp) :: r(2), w(2)
	real(kp) :: BMH(4), KB(2,2), KS(4,4)

	dx = abs(node(2) - node(1))
	
	!------------------
	! bending stiffness
	!------------------
	c1 = EI / dx
	KB(1,:) = (/ c1, -c1 /)
	KB(2,:) = (/-c1,  c1 /)

	!-------------------
	! shearing stiffness
	!-------------------
	select case(etype)
	case("DISP")
		c1 = dx * 0.5_kp
		c2 = 2.0_kp / dx
		call GAUSSLEGENDRE_QUAD(2, r, w) 

		KS = 0.0_kp
		do i = 1, 2
			call SHEAR_BMH(r(i), BMH)

			BMH(1:2) = c2 * BMH(1:2)
			KS = KS + spread(BMH, 2, 4) * spread(BMH, 1, 4) * c1 * w(i)
		end do
	case("MITC")
		c1 = GAk / dx
		c2 = dx * 0.25_kp
		c3 = GAk * 0.5_kp
		KS(1,:) = (/ c1, -c1,  c3,  c3 /)
		KS(2,:) = (/-c1,  c1, -c3, -c3 /)
		KS(3,:) = (/ c3, -c3,  c2,  c2 /)
		KS(4,:) = (/ c3, -c3,  c2,  c2 /)
	end select

	Ke = KS
	Ke(3:4,3:4) = Ke(3:4,3:4) + KB
END SUBROUTINE FEM_1DBEAM_STIFF2N
!----------------------------------------------------------
! Calculate element mass matrix
!
! INPUT:
!	node(2) : node position
!	rhoA : rho * A
!	rhoI : rho * I
!	
! OUTPUT:
!	Me(4,4) : element mass matrix

! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 14, April, 2015
!----------------------------------------------------------
SUBROUTINE FEM_1DBEAM_MASS2N(node, rhoA, rhoI, Me)
	real(kp), intent(in) :: node(2)
	real(kp), intent(in) :: rhoA, rhoI
	real(kp), intent(out) :: Me(4,4)

	integer(i4b) :: i
	real(kp) :: detJ
	real(kp) :: r(2), w(2)
	real(kp) :: WM(2,2)
	real(kp) :: h(2), HH(2,2)

	detJ = abs(node(2) - node(1)) * 0.5_kp

	call GAUSSLEGENDRE_QUAD(2, r, w) 

	WM = 0.0_kp
	do i = 1, 2
		call SHAPE2N(r(i), h)
		HH(1,:) = (/ h(1)*h(1), h(1)*h(2) /)
		HH(2,:) = (/ HH(1,2), h(2)*h(2) /)

		WM = WM + HH * detJ * w(i)
	end do

	Me = 0.0_kp
	Me(1:2,1:2) = WM * rhoA
	Me(3:4,3:4) = WM * rhoI
END SUBROUTINE FEM_1DBEAM_MASS2N

! **************************
! * Interpolation matrices *
! **************************
!----------------------------------------------------------
! Construct shear strain interpolation matrix
!
! INPUT:
!	r : parametric coordinate
!
! OUTPUT:
!	BMH : [dh1/dr, dh2/dr, 0, 0] - [0, 0, h1, h2]
!	
! NEED:
!	SHAPE2N
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 13, April, 2016
!----------------------------------------------------------
SUBROUTINE SHEAR_BMH(r, BMH)
	real(kp), intent(in) :: r
	real(kp), intent(out) :: BMH(4)

	real(kp) :: h(2)

	call SHAPE2N(r, h)
	BMH = (/ -0.5_kp, 0.5_kp, h(1), h(2) /)
END SUBROUTINE SHEAR_BMH
!----------------------------------------------------------
! Shape function for two-nodes
!
! INPUT:
!	r : parametric coordinate
!
! OUTPUT:
!	h(2) : shape functions
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 13, April, 2016
!----------------------------------------------------------
SUBROUTINE SHAPE2N(r, h)
	real(kp), intent(in) :: r
	real(kp), intent(out) :: h(2)

	h(1) = 0.5_kp * (1.0_kp - r)
	h(2) = 0.5_kp * (1.0_kp + r)
END SUBROUTINE SHAPE2N

! *******************
! * Basic functions *
! *******************
!----------------------------------------------------------
! Assign Gauss points and their weights for Gauss-Legendre quadrature
!
! INPUT:
!	n : # of points
!
! OUTPUT:
!	Gp(n) : abscissa
!	w(n)  : weight
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 08, February, 2015
!----------------------------------------------------------
SUBROUTINE GAUSSLEGENDRE_QUAD(n,Gp,w)
	integer(i4b), intent(in) :: n
	real(kp), intent(out) :: Gp(n)
	real(kp), intent(out) :: w(n)

	select case(n)
        case(1) ; Gp = 0.0_kp
				   w = 2.0_kp
        case(2) ; Gp = (/-0.5773502691896257_kp, &
				          0.5773502691896257_kp/)
                   w = (/ 1.0_kp, &
				   		  1.0_kp/)
	end select
END SUBROUTINE GAUSSLEGENDRE_QUAD

!==========================================================
END MODULE FEM_1DBEAM_2N

