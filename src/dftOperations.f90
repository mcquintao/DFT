!----------------------------------------------------------------------!
! MÓDULO COM AS PRINCIPAIS FUNÇÕES E SUBROTINAS DOS MÉTODOS DFT E HF   !
!                                                                      !
! SUBROTINAS PRESENTES:                                                !
!                                                                      !
!----------------------------------------------------------------------!

MODULE DftOperations

USE MathLib
USE DiskOperations
CONTAINS

SUBROUTINE overlapMatrix(BASIS,NBASIS,SMAT)

      IMPLICIT NONE
      
      integer, intent(in) :: NBASIS
      integer :: i, j
      real*8, dimension(NBASIS), intent(in) :: BASIS
      real*8, dimension(NBASIS,NBASIS), intent(out) :: SMAT

      SMAT=0.d0
      do i=1,NBASIS
        do j=1,NBASIS
            SMAT(i,j) = overlap(BASIS(i),BASIS(j))
        end do
      end do


END SUBROUTINE overlapMatrix

!=======================================================================
SUBROUTINE TransferenceMatrix(NBASIS,SMAT,XMAT)
! Subrotina para o calculo de matriz de transferencia (X)
! X = U*S^(-1/2)
! XSX' = 1
!=======================================================================
IMPLICIT NONE

REAL*8, PARAMETER :: Tol = 1.0E-12
INTEGER :: i, j
INTEGER, INTENT(IN) :: NBASIS
REAL*8, DIMENSION(NBASIS,NBASIS), INTENT(IN) :: SMAT
REAL*8, DIMENSION(NBASIS,NBASIS), INTENT(OUT) :: XMAT
REAL*8, DIMENSION(NBASIS,NBASIS) :: LOCAL_SMAT, SMAT_VEC
REAL*8 WALBER

LOCAL_SMAT = SMAT
XMAT = 0.0D0
SMAT_VEC = 0.0D0
WALBER = 0.0D0

CALL JACOBI(LOCAL_SMAT,SMAT_VEC,Tol,NBASIS)

DO j = 1,NBASIS
    WALBER = SQRT(LOCAL_SMAT(j,j))  ! AUTO VALORES
    WALBER = 1.0d0/WALBER
    DO i = 1, NBASIS                ! AUTO VETORES
        XMAT(i,j) = SMAT_VEC(i,j)*WALBER
    END DO
END DO

END SUBROUTINE TransferenceMatrix

SUBROUTINE kineticMatrix(BASIS,NBASIS,TMAT)

      IMPLICIT NONE
      
      integer, intent(in) :: NBASIS
      integer :: i, j
      real*8, dimension(NBASIS), intent(in) :: BASIS
      real*8, dimension(NBASIS,NBASIS), intent(out) :: TMAT

      TMAT=0.d0
      do i=1,NBASIS
        do j=1,NBASIS
            TMAT(i,j) = kinetic(BASIS(i),BASIS(j))
        end do
      end do

END SUBROUTINE kineticMatrix


SUBROUTINE potentialMatrix(ZATOM,BASIS,NBASIS,Potmat) 
    IMPLICIT NONE 

    integer, intent(in) :: NBASIS,ZATOM 
    integer :: i, j 
    real*8, dimension(NBASIS), intent(in) :: BASIS 
    real*8, dimension(NBASIS,NBASIS), intent(out) :: Potmat 

    do i=1,NBASIS 
    do j=1,NBASIS 
    Potmat(i,j) = intgPot(ZATOM,BASIS(i),BASIS(j)) 
    end do 
    end do 

END SUBROUTINE potentialMatrix

SUBROUTINE HCor(T,U,H,NBASIS)

    integer, intent(in) :: NBASIS
    integer :: i, j
    real*8, dimension(NBASIS,NBASIS), intent(in) :: T, U
    real*8, dimension(NBASIS,NBASIS), intent(out) :: H
    
    H = 0.d0
    DO i=1,NBASIS
        DO j=1,NBASIS
            H(i,j) = T(i,j) + U(i,j)
        END DO
    END DO

END SUBROUTINE HCor

SUBROUTINE PGuess(PMAT,NBASIS)
    integer :: NBASIS
    real*8, dimension(NBASIS,NBASIS), intent(out) :: PMAT

    PMAT = 0.d0 ! Grande chute inicial  :) !

END SUBROUTINE PGuess



! FUNÇÕES ---------------------------------------------

REAL*8 FUNCTION JIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    JIntegral = JKIntegrals(a,b,c,d)
END FUNCTION JIntegral

REAL*8 FUNCTION KIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    KIntegral = JKIntegrals(a,d,c,b)
END FUNCTION KIntegral

END MODULE DftOperations

