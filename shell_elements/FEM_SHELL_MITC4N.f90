!==============================================================================
! MITC4N shell element
!
! 	    s
!       *
!	2 * * * 1
!	*   *   *
!	*   * * * *r
!	*       *
!	3 * * * 4
!
!	U = [u1,...,u4,v1,...,v4,w1,...,w4,a1,...,a4,b1,...,b4]^T
!
! CONTENTS:
! NEED:
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!==============================================================================
MODULE FEM_SHELL_MITC4N
USE PARAMETERS
IMPLICIT NONE

PUBLIC CROSS3
PUBLIC FEM_SHELLMITC4N_MASS, FEM_SHELLMITC4N_STIFF

PRIVATE INTEGRANT_ME, INTEGRANT_KE, CONSTITUTIVE_SHELL, TRANS_QSH
PRIVATE STRAINS_XYZ, COVSTRAINS, DUDR_MATRIX
PRIVATE PUDPULIN_4N, JACOBIAN3, DET3, ADJDET3

CONTAINS
!==========================================================
! *******************************
! * Construct element ME and KE *
! *******************************
!----------------------------------------------------------
! Compute element stiffness matrix
!
! INPUT:
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!               V(nodal #, basis #, coordinate #)
!   young     = Young's modulus
!   poisson   = Poisson's ratio
!
! OUTPUT:
!   KE(20,20) = element stiffness matrix
!
! NEED:
!   CONSTITUTIVE_SHELL
!   INTEGRANT_KE
!   
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE FEM_SHELLMITC4N_STIFF(NODE,THICK,V,young,poisson,KE)
    real(kp), intent(in) :: young, poisson
    real(kp), intent(in) :: NODE(4,3), THICK(4), V(4,3,3)
    real(kp), intent(out) :: KE(20,20)

    integer(ip) :: i, j, k
    real(kp) :: KINT(20,20), C(6,6)
    real(kp), parameter :: GAUSS2(2) = (/-0.5773502691896257645091_kp,&
                                          0.5773502691896257645091_kp /)

    call CONSTITUTIVE_SHELL(young,poisson,C)

    KE = 0.0_kp
    do i = 1, 2
        do j = 1, 2
            do k = 1, 2
                call INTEGRANT_KE(NODE,THICK,V,C,&
                     (/GAUSS2(i),GAUSS2(j),GAUSS2(k)/),KINT)
                KE = KE + KINT
            end do
        end do
    end do
END SUBROUTINE FEM_SHELLMITC4N_STIFF
!----------------------------------------------------------
! Compute element mass matrix
!
! INPUT:
!   density   = element mass density
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!               V(nodal #, basis #, coordinate #)
!
! OUTPUT:
!   ME(20,20) = element mass matrix
!
! NEED:
!   INTEGRANT_ME
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE FEM_SHELLMITC4N_MASS(NODE,THICK,V,density,ME)
    real(kp), intent(in) :: density
    real(kp), intent(in) :: NODE(4,3), THICK(4), V(4,3,3)
    real(kp), intent(out) :: ME(20,20)

    integer(ip) :: i, j, k
    real(kp) :: MINT(20,20)
    real(kp), parameter :: GAUSS2(2) = (/-0.5773502691896257645091_kp,&
                                          0.5773502691896257645091_kp /)

    ME = 0.0_kp
    do i = 1, 2
        do j = 1, 2
            do k = 1, 2
                call INTEGRANT_ME(NODE,THICK,V,&
                     (/GAUSS2(i),GAUSS2(j),GAUSS2(k)/),density,MINT)
                ME = ME + MINT
            end do
        end do
    end do
END SUBROUTINE FEM_SHELLMITC4N_MASS

! **************************
! * Integrant of ME and KE *
! **************************
!----------------------------------------------------------
! Compute the integrant of element mass matrix
!
! INPUT:
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!   R(3)      = parametric coordinates
!   C(6,6)    = stress-strain laws matrix (shell-aligned Cartesian coordinates)
!
! OUTPUT:
!   A(20,20) = the integrant
!
! NEED:
!   CONSTITUTIVE_SHELL
!   STRAINS_XYZ
!   TRANS_QSH
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 2, July, 2016
!----------------------------------------------------------
SUBROUTINE INTEGRANT_KE(NODE,THICK,V,C,R,A)
    real(kp), intent(in) :: R(3), C(6,6)
    real(kp), intent(in) :: NODE(4,3), THICK(4), V(4,3,3)
    real(kp), intent(out) :: A(20,20)

    real(kp) :: QSH(6,6), EXYZ(6,20), J3(3,3), detJ

    call STRAINS_XYZ(NODE,THICK,V,R,EXYZ,J3,detJ)
    call TRANS_QSH(J3,QSH)

    EXYZ = MATMUL(QSH,EXYZ)
    A = MATMUL(TRANSPOSE(EXYZ), MATMUL(C,EXYZ))
    A = A*detJ
END SUBROUTINE INTEGRANT_KE
!----------------------------------------------------------
! Compute the constitutive laws for shell elements
! - defined in shell-aligned Cartesian coordinates
!
! INPUT:
!   young   = Young's modulus
!   poisson = Poisson's ratio
!
! OUTPUT:
!   C(6,6) = constitutive laws matrix
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 2, July, 2016
!----------------------------------------------------------
SUBROUTINE CONSTITUTIVE_SHELL(young,poisson,C)
    real(kp), intent(in) :: young, poisson
    real(kp), intent(out) :: C(6,6)

    real(kp), parameter :: k = 1.0_kp    ! shear correction factor

    C      = 0.0_kp
    C(1,1) = young / (1.0_kp - poisson**2)
    C(2,2) = C(1,1)
    C(1,2) = poisson * C(1,1)
    C(2,1) = C(1,2)
    C(4,4) = 0.5_kp*(1.0_kp - poisson) * C(1,1)
    C(5,5) = k * C(4,4)
    C(6,6) = C(5,5)
END SUBROUTINE CONSTITUTIVE_SHELL
!----------------------------------------------------------
! Compute transformation matrix QSH
! - shell-aligned Cartesian coordinates to global Cartesian coordinates
!
! INPUT:
!   COVV(3,3) = base vector of shell-aligned Cartesian coordinats
!               [q1, q2, q3]^T
!
! OUTPUT:
!   QSH(6,6) = the transformation matrix
!
! NEED:
!   CROSS3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE TRANS_QSH(COVV,QSH)
    real(kp), intent(in) :: COVV(3,3)
    real(kp), intent(out) :: QSH(6,6)

    real(kp) :: ER(3), ES(3), ET(3)

    call CROSS3(COVV(2,:),COVV(3,:),ER)
    ER = ER / NORM2(ER)

    call CROSS3(COVV(3,:),ER,ES)
    ES = ES / NORM2(ES)

    ET = COVV(3,:) / NORM2(COVV(3,:))

    QSH(1,1) = ER(1)*ER(1); QSH(1,2) = ER(2)*ER(2); QSH(1,3) = ER(3)*ER(3)
    QSH(1,4) = ER(1)*ER(2); QSH(1,5) = ER(2)*ER(3); QSH(1,6) = ER(3)*ER(1)

    QSH(2,1) = ES(1)*ES(1); QSH(2,2) = ES(2)*ES(2); QSH(2,3) = ES(3)*ES(3)
    QSH(2,4) = ES(1)*ES(2); QSH(2,5) = ES(2)*ES(3); QSH(2,6) = ES(3)*ES(1)

    QSH(3,1) = ET(1)*ET(1); QSH(3,2) = ET(2)*ET(2); QSH(3,3) = ET(3)*ET(3)
    QSH(3,4) = ET(1)*ET(2); QSH(3,5) = ET(2)*ET(3); QSH(3,6) = ET(3)*ET(1)

    QSH(4,1) = 2.0_kp*ER(1)*ES(1)
    QSH(4,2) = 2.0_kp*ER(2)*ES(2)
    QSH(4,3) = 2.0_kp*ER(3)*ES(3)
    QSH(4,4) = ER(1)*ES(2) + ES(1)*ER(2)
    QSH(4,5) = ER(2)*ES(3) + ES(2)*ER(3)
    QSH(4,6) = ER(3)*ES(1) + ES(3)*ER(1)

    QSH(5,1) = 2.0_kp*ES(1)*ET(1)
    QSH(5,2) = 2.0_kp*ES(2)*ET(2)
    QSH(5,3) = 2.0_kp*ES(3)*ET(3)
    QSH(5,4) = ES(1)*ET(2) + ET(1)*ES(2)
    QSH(5,5) = ES(2)*ET(3) + ET(2)*ES(3)
    QSH(5,6) = ES(3)*ET(1) + ET(3)*ES(1)

    QSH(6,1) = 2.0_kp*ET(1)*ER(1)
    QSH(6,2) = 2.0_kp*ET(2)*ER(2)
    QSH(6,3) = 2.0_kp*ET(3)*ER(3)
    QSH(6,4) = ET(1)*ER(2) + ER(1)*ET(2)
    QSH(6,5) = ET(2)*ER(3) + ER(2)*ET(3)
    QSH(6,6) = ET(3)*ER(1) + ER(3)*ET(1)
END SUBROUTINE TRANS_QSH

!----------------------------------------------------------
! Compute the integrant of element mass matrix
!
! INPUT:
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!   R(3)      = parametric coordinates
!   density   = element mass density
!
! OUTPUT:
!   A(20,20) = the integrant
!
! NEED:
!   PUDPULIN_4N
!   JACOBIAN3
!   DET3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE INTEGRANT_ME(NODE,THICK,V,R,density,A)
    real(kp), intent(in) :: R(3), density
    real(kp), intent(in) :: NODE(4,3), THICK(4), V(4,3,3)
    real(kp), intent(out) :: A(20,20)

    integer(ip) :: i, j, k, l, rs, re, cs, ce
    real(kp) :: hr3, WM(4,4)
    real(kp) :: J3(3,3), detJ
    real(kp) :: H4(4), DH4(2,4)
    real(kp) :: HT(3,8)

    call PUDPULIN_4N(R(1:2),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call DET3(J3,detJ)

    hr3 = 0.5_kp*R(3)
    do i = 1, 4
        HT(1,i) = hr3*THICK(i)*H4(i)
    end do
    HT(2,1:4) = HT(1,1:4)
    HT(3,1:4) = HT(1,1:4)
    HT(:,5:8) = HT(:,1:4)

    do i = 1, 3
        do j = 1, 4
            HT(i,j)   = -HT(i,j)*V(j,3,i)
            HT(i,4+j) =  HT(i,4+j)*V(j,2,i)
        end do
    end do

    ! Construct A
    A = 0.0_kp
    do i = 1, 4
        do j = i, 4
            A(i,j) = H4(i)*H4(j)*density*detJ
            if (j /= i) then
                A(j,i) = A(i,j)
            end if
        end do
    end do
    A(5:8,5:8)   = A(1:4,1:4)
    A(9:12,9:12) = A(1:4,1:4)

    do i = 1, 3
        re = 4*i
        rs = re - 3
        do j = 1, 2
            ce = 12 + 4*j
            cs = ce - 3
            do k = 1, 4
                do l = 1, 4
                    WM(k,l) = H4(k)*HT(i,l+4*(j-1))*density*detJ
                end do
            end do
            A(rs:re,cs:ce) = WM
            A(cs:ce,rs:re) = TRANSPOSE(WM)
        end do
    end do

    do i = 1, 8
        rs = 12 + i
        do j = i, 8
            cs = 12 + j
            A(rs,cs) = DOT_PRODUCT(HT(:,i),HT(:,j))*density*detJ
            if (cs /= rs) then
                A(cs,rs) = A(rs,cs)
            end if
        end do
    end do
END SUBROUTINE INTEGRANT_ME

! *************************
! * Compute strain matrix *
! *************************
!----------------------------------------------------------
! Compute strain matrix in global Cartesian coordinates
!
! INPUT:
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!   R(3)      = parametric coordinates
!
! OUTPUT:
!   EXYZ(6,20) = strain matrix
!                [Exx, Eyy, Ezz, 2*Exy, 2*Eyz, 2*Ezx]^T
!   J3(3,3)    = Jacobian
!   detJ       = determinant of Jacobian
!
! NEED:
!   PUDPULIN_4N
!   JACOBIAN3
!   COVSTRAINS
!   ADJDET3
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE STRAINS_XYZ(NODE,THICK,V,R,EXYZ,J3,detJ)
    real(kp), intent(in) :: NODE(4,3), THICK(4), V(4,3,3)
    real(kp), intent(in) :: R(3)
    real(kp), intent(out) :: EXYZ(6,20), J3(3,3), detJ

    integer(ip) :: i
    integer(ip), parameter :: trid(6,2) = (/1,2,3,1,2,3, 1,2,3,2,3,1/)
    real(kp) :: H4(4), DH4(2,4), CONTV(3,3)
    real(kp) :: ECOV(5,20)

    !----------------------------------------------------!
    ! compute assumed transverse shear covariant strains !
    !----------------------------------------------------!
    call PUDPULIN_4N((/1.0_kp,0.0_kp/),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call COVSTRAINS(THICK,V,H4,DH4,R(3),J3,ECOV(1:3,:),4)
    ECOV(4,:) = 0.5_kp*(1.0_kp + R(1))*ECOV(1,:)

    call PUDPULIN_4N((/-1.0_kp,0.0_kp/),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call COVSTRAINS(THICK,V,H4,DH4,R(3),J3,ECOV(1:3,:),4)
    ECOV(4,:) = ECOV(4,:) + 0.5_kp*(1.0_kp - R(1))*ECOV(1,:)

    call PUDPULIN_4N((/0.0_kp,1.0_kp/),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call COVSTRAINS(THICK,V,H4,DH4,R(3),J3,ECOV(1:3,:),5)
    ECOV(5,:) = 0.5_kp*(1.0_kp + R(2))*ECOV(1,:)

    call PUDPULIN_4N((/0.0_kp,-1.0_kp/),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call COVSTRAINS(THICK,V,H4,DH4,R(3),J3,ECOV(1:3,:),5)
    ECOV(5,:) = ECOV(5,:) + 0.5_kp*(1.0_kp - R(2))*ECOV(1,:)

    !------------------------------------!
    ! compute in-layer covariant strains !
    !------------------------------------!
    call PUDPULIN_4N(R(1:2),H4,DH4)
    call JACOBIAN3(NODE,THICK,V(:,1,:),H4,DH4,R(3),J3)
    call COVSTRAINS(THICK,V,H4,DH4,R(3),J3,ECOV(1:3,:),1)   
    
    !---------------------------------------------!
    ! compute detJ and contravariant base vectors !
    !---------------------------------------------!
    call ADJDET3(J3,CONTV,detJ)
    CONTV = CONTV / detJ

    !--------------!
    ! compute EXYZ !
    !--------------!
    do i = 1, 6
        EXYZ(i,:) = ECOV(1,:) * CONTV(trid(i,1),1)*CONTV(trid(i,2),1) +&
                    ECOV(2,:) * CONTV(trid(i,1),2)*CONTV(trid(i,2),2) +&
                    ECOV(3,:) *(CONTV(trid(i,1),1)*CONTV(trid(i,2),2) +&
                                CONTV(trid(i,1),2)*CONTV(trid(i,2),1))+&
                    ECOV(4,:) *(CONTV(trid(i,1),2)*CONTV(trid(i,2),3) +&
                                CONTV(trid(i,1),3)*CONTV(trid(i,2),2))+&
                    ECOV(5,:) *(CONTV(trid(i,1),3)*CONTV(trid(i,2),1) +&
                                CONTV(trid(i,1),1)*CONTV(trid(i,2),3))
    end do
    EXYZ(4,:) = 2.0_kp*EXYZ(4,:)
    EXYZ(5,:) = 2.0_kp*EXYZ(5,:)
    EXYZ(6,:) = 2.0_kp*EXYZ(6,:)
END SUBROUTINE STRAINS_XYZ
!----------------------------------------------------------
! Compute covariant strain components
!
! INPUT:
!   THICK(4)   = nodal thickness
!   V(4,3,3)   = nodal directional basis vectors
!   H4(4)      = PU
!   DH4(2,4)   = DPU
!   t          = the parametric coordniate in thickness direction
!   COVV(3,3)  = covariant base vectors stored as
!                [q1, q2, q3]^T
!   flag       = 1: in-layer strains
!                2: transverse shear strains
!
! OUTPUT:
!   ECOV(3,20) = covariant strain components matrix
!                [Err, Ess, Ers]^T when flag = 1 or 2 or 3
!                [Est, 0, 0]^T when flag = 4
!                [Ert, 0, 0]^T when flag = 5
! NEED:
!   DUDR_MATRIX
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE COVSTRAINS(THICK,V,H4,DH4,t,COVV,ECOV,flag)
    integer(ip), intent(in) :: flag
    real(kp), intent(in)    :: t
    real(kp), intent(in)    :: THICK(4), V(4,3,3)
    real(kp), intent(in)    :: H4(4), DH4(2,4)
    real(kp), intent(in)    :: COVV(3,3)
    real(kp), intent(out) :: ECOV(3,20)

    real(kp) :: DUDR(3,3,20)

    call DUDR_MATRIX(THICK,V,H4,DH4,t,DUDR)

    ECOV = 0.0_kp
    select case (flag)
    case(1:3)
        ECOV(1,:) = MATMUL(COVV(1,:),DUDR(:,1,:))
        ECOV(2,:) = MATMUL(COVV(2,:),DUDR(:,2,:))
        ECOV(3,:) = 0.5_kp*( MATMUL(COVV(1,:),DUDR(:,2,:)) +&
                             MATMUL(COVV(2,:),DUDR(:,1,:)) )
    case(4)
        ECOV(1,:) = 0.5_kp*( MATMUL(COVV(2,:),DUDR(:,3,:)) +&
                             MATMUL(COVV(3,:),DUDR(:,2,:)) )
    case(5)
        ECOV(1,:) = 0.5_kp*( MATMUL(COVV(1,:),DUDR(:,3,:)) +&
                             MATMUL(COVV(3,:),DUDR(:,1,:)) )
    end select
END SUBROUTINE COVSTRAINS
!----------------------------------------------------------
! Compute dU/dR matrix
!
! INPUT:
!   THICK(4)  = nodal thickness
!   V(4,3,3)  = nodal directional basis vectors
!   H4(4)     = PU
!   DH4(2,4)  = DPU
!   t         = the parametric coordniate in thickness direction
!
! OUTPUT:
!   DUDR(3,3,20) = dU/DR matrix
!                  DUDR(displacement #, coordinate #, :)
!   
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 2, July, 2016
!----------------------------------------------------------
SUBROUTINE DUDR_MATRIX(THICK,V,H4,DH4,t,DUDR)
    real(kp), intent(in) :: t
    real(kp), intent(in) :: THICK(4), V(4,3,3)
    real(kp), intent(in) :: H4(4), DH4(2,4)
    real(kp), intent(out) :: DUDR(3,3,20)

    integer(ip) :: i, j, k, st, en
    integer(ip), parameter :: idv(2) = (/3,2/)
    real(kp) :: ht

    DUDR = 0.0_kp
    do i = 1, 3
        en = 4*i
        st = en - 3
        DUDR(i,1,st:en) = DH4(1,:)
        DUDR(i,2,st:en) = DH4(2,:)
    end do

    ht = 0.5_kp*t
    do k = 1, 4
        st = 12 + k
        DUDR(1,1,st) = ht*THICK(k)
        DUDR(1,2,st) = DUDR(1,1,st)
        DUDR(1,3,st) = 0.5_kp*THICK(k)

        en = st + 4
        do j = 1, 3
            DUDR(1,j,en) = DUDR(1,j,st)
        end do
    end do
    DUDR(2,:,13:20) = DUDR(1,:,13:20)
    DUDR(3,:,13:20) = DUDR(1,:,13:20)

    do i = 1, 3
        do j = 1, 2
            en = 8 + 4*j
            do k = 1, 4
                st = en + k
                DUDR(i,1,st) = DUDR(i,1,st)*DH4(1,k)*V(k,idv(j),i)*(-1)**j
                DUDR(i,2,st) = DUDR(i,2,st)*DH4(2,k)*V(k,idv(j),i)*(-1)**j
                DUDR(i,3,st) = DUDR(i,3,st)*H4(k)*V(k,idv(j),i)*(-1)**j
            end do
        end do
    end do
END SUBROUTINE DUDR_MATRIX

! **********************
! * Auxiliary routines *
! **********************
!----------------------------------------------------------
! Partition of unity linear (hat) functions and derivatives of them for 4-nodes
!
! INPUT:
!   R(2) = parametric coordinates defined in [-1,1]
!
! OUTPUT:
!   H4(4)    = PU linear functions
!   DH4(2,4) = derivatives of PU linear functions wrt
!              r(1) and r(2), respectively
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE PUDPULIN_4N(R,H4,DH4)
    real(kp), intent(in) :: R(2)
    real(kp), intent(out) :: H4(4), DH4(2,4)

    integer(ip) :: i
    real(kp) :: HI(2,2)
    real(kp) :: WV(2)

    do i = 1, 2
        HI(i,1) = 1.0_kp - R(i)
        HI(i,2) = 1.0_kp + R(i)
    end do

    WV = HI(1,:)*0.25_kp

    H4(1) = WV(2)*HI(2,2)
    H4(2) = WV(1)*HI(2,2)
    H4(3) = WV(1)*HI(2,1)
    H4(4) = WV(2)*HI(2,1)

    DH4(2,1) =  WV(2)
    DH4(2,2) =  WV(1)
    DH4(2,3) = -WV(1)
    DH4(2,4) = -WV(2)

    WV = HI(2,:)*0.25_kp

    DH4(1,1) =  WV(2)
    DH4(1,2) = -WV(2)
    DH4(1,3) = -WV(1)
    DH4(1,4) =  WV(1)
END SUBROUTINE PUDPULIN_4N
!----------------------------------------------------------
! Compute 3x3 Jacobian matrix
!
! INPUT:
!   NODE(4,3) = nodal points
!   THICK(4)  = nodal thickness
!   V0(4,3)   = nodal direct vectors
!   H4(4)     = PU
!   DH4(2,4)  = DPU
!   t         = the parametric coordniate in thickness direction
!
! OUTPUT:
!   J3(3,3) = Jacobian matrix
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE JACOBIAN3(NODE,THICK,V0,H4,DH4,t,J3)
    real(kp), intent(in) :: t
    real(kp), intent(in) :: NODE(4,3)
    real(kp), intent(in) :: THICK(4), V0(4,3)
    real(kp), intent(in) :: H4(4), DH4(2,4)
    real(kp), intent(out) :: J3(3,3)

    integer(ip) :: i
    real(kp) :: ht
    real(kp) :: WM(3,4)

    ht = 0.5_kp*t
    do i = 1, 4
        WM(1,i) = ht*THICK(i)*DH4(1,i)
        WM(2,i) = ht*THICK(i)*DH4(2,i)
        WM(3,i) = 0.5*THICK(i)*H4(i)
    end do

    J3 = MATMUL(WM,V0)
    J3(1:2,:) = J3(1:2,:) + MATMUL(DH4,NODE)
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
!----------------------------------------------------------
! Compute a cross product
!
! INPUT:
!   A(3), B(3) = input vectors
!
! OUTPUT:
!   C(3) = output vector
!
! AUTHOR:
!   Ki-Tae Kim (qlsn55@gmail.com), 02, July, 2016
!----------------------------------------------------------
SUBROUTINE CROSS3(A,B,C)
    real(kp), intent(in) :: A(3), B(3)
    real(kp), intent(out) :: C(3)

    C(1) = A(2)*B(3) - A(3)*B(2)
    C(2) = A(3)*B(1) - A(1)*B(3)
    C(3) = A(1)*B(2) - A(2)*B(1)
END SUBROUTINE CROSS3
!==========================================================
END MODULE FEM_SHELL_MITC4N

