!==============================================================================
! 8-nodes element for 3D solids FEA
! 
!          t
!          #
!       3**#******2
!     *    #    * *
!   4*********1   *
!   *         * # # # s
!   *   7#    *   6
!   *  #      * *
!   8#********5
!  r
!
! - Isoparametric FE
! - U = [u1, ..., u8, v1, ..., v8, w1, ..., w8]^T
! - X = [x1, ..., x8, y1, ..., y8, z1, ..., z8]^T
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!==============================================================================
MODULE FEM_3DSOLID_8N
USE PARAMETERS
IMPLICIT NONE

PUBLIC FEM_3DS8N_STIFF, FEM_3DS8N_MASS
PRIVATE INTEGRANT_KE, CONSTITUTIVE_3DS, INTEGRANT_ME
PRIVATE PULIN_8N, DPULIN_8N, JACOBIAN3, DET3, ADJDET3

CONTAINS
!==========================================================
! *******************************
! * Construct element ME and KE *
! *******************************
!----------------------------------------------------------
! Compute element stiffness matrix
!
! INPUT:
!   young     = Young's modulus
!   poisson   = Poisson's ratio
!   NODE(8,3) = nodal points
!
! OUTPUT:
!   KE(24,24) = element stiffness matrix
!
! NEED:
!   CONSTITUTIVE_3DS
!   INTEGRANT_KE
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE FEM_3DS8N_STIFF(young,poisson,NODE,KE)
    real(kp), intent(in) :: young, poisson
    real(kp), intent(in) :: NODE(8,3)
    real(kp), intent(out) :: KE(24,24)

    integer(ip) :: i, j, k
    real(kp) :: C(6,6), KINT(24,24)
    real(kp), parameter :: GAUSS2(2) = (/-0.5773502691896257645091_kp,&
                                          0.5773502691896257645091_kp /)

    call CONSTITUTIVE_3DS(young,poisson,C)

    KE = 0.0_kp
    do i = 1, 2
        do j = 1, 2
            do k = 1, 2
                call INTEGRANT_KE(C,NODE,&
                     (/GAUSS2(i),GAUSS2(j),GAUSS2(k)/),KINT)
                KE = KE + KINT
            end do
        end do
    end do
END SUBROUTINE FEM_3DS8N_STIFF
!----------------------------------------------------------
! Compute element mass matrix
!
! INPUT:
!   density   = mass density
!   NODE(8,3) = nodal points
!
! OUTPUT:
!   ME(24,24) = element mass matrix
!
! NEED:
!   INTEGRANT_ME
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE FEM_3DS8N_MASS(density,NODE,ME)
    real(kp), intent(in) :: density
    real(kp), intent(in) :: NODE(8,3)
    real(kp), intent(out) :: ME(24,24)

    integer(ip) :: i, j, k
    real(kp) :: MINT(8,8)
    real(kp), parameter :: GAUSS2(2) = (/-0.5773502691896257645091_kp,&
                                          0.5773502691896257645091_kp /)

    ME = 0.0_kp
    do i = 1, 2
        do j = 1, 2
            do k = 1, 2
                call INTEGRANT_ME(density,NODE,&
                     (/GAUSS2(i),GAUSS2(j),GAUSS2(k)/),MINT)
                ME(1:8,1:8) = ME(1:8,1:8) + MINT
            end do
        end do
    end do

    ME(9:16,9:16)   = ME(1:8,1:8)
    ME(17:24,17:24) = ME(1:8,1:8)
END SUBROUTINE FEM_3DS8N_MASS

! **************************
! * Integrant of ME and KE *
! **************************
!----------------------------------------------------------
! Compute the integrant of element stiffness matrix
!
! INPUT:
!   C(3,3)    = Constitutive matrix
!   NODE(8,3) = nodal points
!   R(3)      = parametric coordinates
!
! OUTPUT:
!   A(24,24) = the integrant
!
! NEED:
!   DPULIN_8N
!   JACOBIAN3
!   ADJDET3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE INTEGRANT_KE(C,NODE,R,A)
    real(kp), intent(in) :: C(6,6), NODE(8,3), R(3)
    real(kp), intent(out) :: A(24,24)

    real(kp) :: BX(3,8), B(6,24)
    real(kp) :: J3(3,3), ADJ(3,3), detJ

    call DPULIN_8N(R,BX)
    call JACOBIAN3(NODE,R,J3)
    call ADJDET3(J3,ADJ,detJ)
    BX = MATMUL(ADJ,BX)

    B          = 0.0_kp
    B(1,1:8)   = BX(1,:)
    B(2,9:16)  = BX(2,:)
    B(3,17:24) = BX(3,:)
    B(4,1:8)   = BX(2,:); B(4,9:16)  = BX(1,:)
    B(5,9:16)  = BX(3,:); B(5,17:24) = BX(2,:)
    B(6,1:8)   = BX(3,:); B(6,17:24) = BX(1,:)

    A = MATMUL( MATMUL(TRANSPOSE(B),C), B )
    A = A / detJ
END SUBROUTINE INTEGRANT_KE
!----------------------------------------------------------
! Compute constitutive laws for 3D solids element
!
! INPUT:
!   young   = Young's modulus
!   poisson = Poisson's ratio
!
! OUTPUT:
!   C(6,6) = element constitutive matrix
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE CONSTITUTIVE_3DS(young,poisson,C)
    real(kp), intent(in) :: young, poisson
    real(kp), intent(out) :: C(6,6)

    real(kp) :: tempv

    C = 0.0_kp

    tempv  = young / (1.0_kp + poisson)
    C(4,4) = tempv*0.5_kp

    tempv  = tempv / (1.0_kp - 2.0_kp*poisson)
    C(1,2) = tempv*poisson
    C(1,1) = tempv*(1.0_kp - poisson)

    C(2,2) = C(1,1); C(3,3) = C(1,1)
    C(1,3) = C(1,2)
    C(2,1) = C(1,2); C(2,3) = C(1,2)
    C(3,1) = C(1,2); C(3,2) = C(1,2)
    C(5,5) = C(4,4); C(6,6) = C(4,4)
END SUBROUTINE CONSTITUTIVE_3DS
!----------------------------------------------------------
! Compute the integrant of element mass matrix
!
! INPUT:
!   density   = mass density
!   NODE(8,3) = nodal points
!   R(3)      = parametric coordinates
!
! OUTPUT:
!   A(8,8) = the integrant
!
! NEED:
!   PULIN_8N
!   JACOBIAN3
!   DET3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE INTEGRANT_ME(density,NODE,R,A)
    real(kp), intent(in) :: density, NODE(8,3), R(3)
    real(kp), intent(out) :: A(8,8)

    integer(ip) :: i, j
    real(kp) :: H8(8)
    real(kp) :: J3(3,3), detJ3

    call PULIN_8N(R,H8)
    call JACOBIAN3(NODE,R,J3)
    call DET3(J3,detJ3)
    detJ3 = detJ3*density

    do i = 1, 8
        do j = i, 8
            A(i,j) = H8(i)*H8(j)*detJ3
            if (j /= i) then
                A(j,i) = A(i,j)
            end if
        end do
    end do
END SUBROUTINE INTEGRANT_ME

! **********************
! * Auxiliary routines *
! **********************
!----------------------------------------------------------
! Partition of unity linear (hat) functions for 8-nodes
!
! INPUT:
!   R(3) = parametric coordinates defined in [-1,1]
!
! OUTPUT:
!   H8(8) = PU linear functions
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE PULIN_8N(R,H8)
    real(kp), intent(in) :: R(3)
    real(kp), intent(out) :: H8(8)

    integer(ip) :: i
    real(kp) :: HI(3,2)

    do i = 1, 3
        HI(i,1) = 1.0_kp - R(i)
        HI(i,2) = 1.0_kp + R(i)
    end do
    HI(1,:) = HI(1,:)*0.125_kp

    H8(1) = HI(1,2)*HI(2,2)*HI(3,2)
    H8(2) = HI(1,1)*HI(2,2)*HI(3,2)
    H8(3) = HI(1,1)*HI(2,1)*HI(3,2)
    H8(4) = HI(1,2)*HI(2,1)*HI(3,2)

    H8(5) = HI(1,2)*HI(2,2)*HI(3,1)
    H8(6) = HI(1,1)*HI(2,2)*HI(3,1)
    H8(7) = HI(1,1)*HI(2,1)*HI(3,1)
    H8(8) = HI(1,2)*HI(2,1)*HI(3,1)
END SUBROUTINE PULIN_8N
!----------------------------------------------------------
! Derivatives of Partition of unity linear (hat) functions for 8-nodes
!
! INPUT:
!   R(3) = parametric coordinates defined in [-1,1]
!
! OUTPUT:
!   DH8(3,8) = derivatives of PU linear functions wrt
!              r(1), r(2) and r(3), respectively
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE DPULIN_8N(R,DH8)
    real(kp), intent(in) :: R(3)
    real(kp), intent(out) :: DH8(3,8)

    integer(ip) :: i
    real(kp) :: HI(3,2)
    real(kp) :: WV(2)

    do i = 1, 3
        HI(i,1) = 1.0_kp - R(i)
        HI(i,2) = 1.0_kp + R(i)
    end do

    ! d/dr
    WV(:)    =  0.125_kp*HI(2,:)
    DH8(1,1) =  WV(2)*HI(3,2)
    DH8(1,2) = -WV(2)*HI(3,2)
    DH8(1,3) = -WV(1)*HI(3,2)
    DH8(1,4) =  WV(1)*HI(3,2)
    DH8(1,5) =  WV(2)*HI(3,1)
    DH8(1,6) = -WV(2)*HI(3,1)
    DH8(1,7) = -WV(1)*HI(3,1)
    DH8(1,8) =  WV(1)*HI(3,1)

    ! d/ds
    WV(:)    =  0.125_kp*HI(1,:)
    DH8(2,1) =  WV(2)*HI(3,2)
    DH8(2,2) =  WV(1)*HI(3,2)
    DH8(2,3) = -WV(1)*HI(3,2)
    DH8(2,4) = -WV(2)*HI(3,2)
    DH8(2,5) =  WV(2)*HI(3,1)
    DH8(2,6) =  WV(1)*HI(3,1)
    DH8(2,7) = -WV(1)*HI(3,1)
    DH8(2,8) = -WV(2)*HI(3,1)

    ! d/dt
    DH8(3,1) =  WV(2)*HI(2,2)
    DH8(3,2) =  WV(1)*HI(2,2)
    DH8(3,3) =  WV(1)*HI(2,1)
    DH8(3,4) =  WV(2)*HI(2,1)
    DH8(3,5) = -WV(2)*HI(2,2)
    DH8(3,6) = -WV(1)*HI(2,2)
    DH8(3,7) = -WV(1)*HI(2,1)
    DH8(3,8) = -WV(2)*HI(2,1)
END SUBROUTINE DPULIN_8N
!----------------------------------------------------------
! Compute 3x3 Jacobian matrix
!
! INPUT:
!   NODE(8,3) = nodal points
!   R(3)      = parametric coordinates
!
! OUTPUT:
!   J3(3,3) = Jacobian matrix
!
! NEED:
!   DPULIN_8N
!   
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE JACOBIAN3(NODE,R,J3)
    real(kp), intent(in) :: NODE(8,3), R(3)
    real(kp), intent(out) :: J3(3,3)

    real(kp) :: DH8(3,8)

    call DPULIN_8N(R,DH8)
    J3 = MATMUL(DH8,NODE)
END SUBROUTINE JACOBIAN3
!----------------------------------------------------------
! Compute determinant of 3x3 matrix A
!
! INPUT:
!   A(3,3)
!
! OUTPUT:
!   detA
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE DET3(A,detA)
    real(kp), intent(in) :: A(3,3)
    real(kp), intent(out) :: detA

    detA = A(1,1)*(A(2,2)*A(3,3) - A(2,3)*A(3,2)) &
         - A(1,2)*(A(2,1)*A(3,3) - A(2,3)*A(3,1)) &
         + A(1,3)*(A(2,1)*A(3,2) - A(2,2)*A(3,1))
END SUBROUTINE DET3
!----------------------------------------------------------
! Compute adjugate and determinant of 3x3 matrix A
!
! INPUT:
!   A(3,3)
!
! OUTPUT:
!   ADJ(3,3) = adjugate of A
!   detA     = determinant of A
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 27, June, 2016
!----------------------------------------------------------
SUBROUTINE ADJDET3(A,ADJ,detA)
    real(kp), intent(in) :: A(3,3)
    real(kp), intent(out) :: ADJ(3,3), detA

    ADJ(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    ADJ(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    ADJ(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
    ADJ(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    ADJ(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    ADJ(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
    ADJ(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    ADJ(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    ADJ(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    detA = A(1,1)*ADJ(1,1) + A(1,2)*ADJ(2,1) + A(1,3)*ADJ(3,1)
END SUBROUTINE ADJDET3
!==========================================================
END MODULE FEM_3DSOLID_8N

