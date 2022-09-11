!==============================================================================
! Basic functions for plane stress & plane strain FEA
!
! CONTENTS:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 24, March, 2015
!==============================================================================
MODULE FEM_PLANE_BASICFS
USE PARAMETERS
IMPLICIT NONE

PUBLIC INVERSE2, DET2
PUBLIC GAUSSLEGENDRE_QUAD

CONTAINS
!==========================================================
! ******************
! * Linear algebra *
! ******************
!----------------------------------------------------------
! Calculate 2x2 inverse matrix
!
! INPUT:
!	A(2,2)
!
! OUTPUT:
!	B(2,2) : inverse of A
!	detA : determinant of A
!
! NEED:
!	DET2
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 30, March, 2015
!----------------------------------------------------------
SUBROUTINE INVERSE2(A,B,detA)
	real(kp), intent(in) :: A(2,2)
	real(kp), intent(out) :: B(2,2)
	real(kp), intent(out) :: detA

	call DET2(A, detA)

	B(1,1) = A(2,2)/detA
	B(2,2) = A(1,1)/detA
	B(1,2) = -A(1,2)/detA
	B(2,1) = -A(2,1)/detA
END SUBROUTINE INVERSE2
!----------------------------------------------------------
! Calculate the determinant of 2x2 matrix
!
! INPUT:
!	A(2,2)
!
! OUTPUT:
!	detA
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 28, April, 2016
!----------------------------------------------------------
SUBROUTINE DET2(A,detA)
	real(kp), intent(in) :: A(2,2)
	real(kp), intent(out) :: detA

	detA = A(1,1)*A(2,2) - A(1,2)*A(2,1)
END SUBROUTINE DET2

! ***************************
! * Gauss quadrature tables *
! ***************************
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
	integer(ip), intent(in) :: n
	real(kp), intent(out) :: Gp(n)
	real(kp), intent(out) :: w(n)

	select case(n)
        case(1) ; Gp = 0.0_kp
				   w = 2.0_kp
        case(2) ; Gp = (/-0.5773502691896257_kp, &
				          0.5773502691896257_kp/)
                   w = (/ 1.0_kp, &
				   		  1.0_kp/)
        case(3) ; Gp = (/-0.7745966692414834_kp, &
						  0.0_kp, &
					      0.7745966692414834_kp/)
                   w = (/ 0.5555555555555556_kp, &
				   		  0.8888888888888889_kp, &
						  0.5555555555555556_kp/)
        case(4) ; Gp = (/-0.8611363115940526_kp, &
						 -0.3399810435848563_kp, &
						  0.3399810435848563_kp, &
						  0.8611363115940526_kp/)
                   w = (/ 0.3478548451374538_kp, &
					      0.6521451548625461_kp, &
						  0.6521451548625461_kp, &
						  0.3478548451374538_kp/)
	end select
END SUBROUTINE GAUSSLEGENDRE_QUAD

!==========================================================
END MODULE FEM_PLANE_BASICFS

