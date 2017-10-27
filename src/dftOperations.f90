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



SUBROUTINE TransfFMAT(FMAT,XMAT,NBASIS)
    INTEGER, intent(in) :: NBASIS
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: XMAT
    REAL*8, dimension(NBASIS,NBASIS), intent(inout) :: FMAT
    REAL*8, dimension(NBASIS,NBASIS) :: TXMAT

    TXMAT = TRANSPOSE(XMAT)
    FMAT = MATMUL(TXMAT,MATMUL(FMAT,XMAT))

END SUBROUTINE TransfFMAT

SUBROUTINE newDensity(CMAT,P,NBASIS,NE)
    INTEGER, intent(in) :: NBASIS, NE
    INTEGER :: i, j, a
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: CMAT
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: P

    P = 0.d0
    DO i=1,NBASIS
        DO j=1,NBASIS
            DO a=1,NE/2
                P(i,j) = P(i,j) + 2.d0 * CMAT(i,a)*CMAT(j,a)  !! PAG. 139 EQ:3.145; SZABO (CONFERIR TRANSPOSTA)
            END DO
        END DO
    END DO

END SUBROUTINE newDensity



REAL*8 FUNCTION ENERGY(H,F,P,NB)
    
    INTEGER, intent(in) :: NB
    INTEGER :: i, j
    REAL*8, dimension(NB,NB), intent(in) :: H, F, P
    
    ENERGY = 0.d0
    DO i = 1, NB
        DO j = 1, NB
            ENERGY = 0.5D0 * P(j,i)*(H(i,j) + F(i,j)) !! PAG. 150 EQ: 3.184; SZABO. P(j,i) => Eh isso msm ??
        END DO
    END DO

END FUNCTION ENERGY


! SUBROTINA SCF!! -----------------------------------

SUBROUTINE SCF(XMAT,HCORE,BASIS,NBASIS,NE,PMAT,CMAT)
    REAL*8, PARAMETER :: Tol = 1.0E-12
    INTEGER, intent(in) :: NBASIS, NE
    INTEGER :: i,j,k,l,loop
    REAL*8, dimension(NBASIS), intent(in) :: BASIS
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: HCORE
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: PMAT, CMAT
    REAL*8, dimension(NBASIS,NBASIS) :: XMAT, FMAT, GMAT, XFMAT, NewPMAT
    REAL*8 :: a,b,c,d, TotEnergy

    TotEnergy = 0.d0
    NewPMAT = 0.d0

    DO loop=1,100
    
    PRINT *, "SCF Cycle ", loop

    XFMAT   = 0.d0
    GMAT    = 0.d0
    FMAT    = 0.d0
!    PMAT    = 1.0d0
    DO i=1,NBASIS
        a = BASIS(i)
        DO j=1,NBASIS
            b = BASIS(j)

            DO k=1,NBASIS
                c = BASIS(k)
                DO l=1,NBASIS
                    d = BASIS(l)
                    GMAT(i,j) = GMAT(i,j) + PMAT(k,l)*(JIntegral(a,b,c,d) - 0.5d0*KIntegral(a,b,c,d))
                END DO
            END DO

        END DO
    END DO

    call dbgMatrix(GMAT,NBASIS,"G MATRIX",8)
    
    DO i=1,NBASIS
        DO j=1,NBASIS
            FMAT(i,j) = HCORE(i,j) + GMAT(i,j)
        END DO
    END DO

call dbgMatrix(FMAT,NBASIS,"F MATRIX",8)

    XFMAT = FMAT
    call TransfFMAT(XFMAT,XMAT,NBASIS)  ! XFMAT = XMAT' * FMAT * XMAT
    
call dbgMatrix(XFMAT,NBASIS,"XF MATRIX",9)

    call FLUSH(99)
    
    call JACOBI(XFMAT,CMAT,Tol,NBASIS) ! XFMAT -> ENERGY(A. VAL); CMAT -> XCOEFF(A. VEC)
    
call dbgMatrix(XFMAT,NBASIS,"ENERGY MATRIX",13)
call dbgMatrix(CMAT,NBASIS,"XCOEFF MATRIX",13)

    CMAT = MATMUL(XMAT,CMAT)

call dbgMatrix(CMAT,NBASIS,"COEFF MATRIX",12)

    call newDensity(CMAT,NewPMAT,NBASIS,NE)

call dbgMatrix(NewPMAT,NBASIS,"NEW DENSITY MATRIX",18)

    PMAT = NewPMAT


    TotEnergy = ENERGY(HCORE,FMAT,PMAT,NBASIS)
    

    PRINT *, "--------------------------------"
    PRINT *, "TOTAL ENERGY: ", TotEnergy
    PRINT *, "--------------------------------"

    IF(ABS(TotEnergy - OldEnergy).lt.10E-4) THEN
        GOTO 150
    END IF

    OldEnergy = TotEnergy

    END DO

100 PRINT *, "O CICLO SCF NÃO CONVERGIU"
GOTO 200

150 PRINT *, "O CICLO SCF CONVERGIU!"
GOTO 200

200 PRINT *, "FIM DO SCF"
END SUBROUTINE SCF



! FUNÇÕES ---------------------------------------------

REAL*8 FUNCTION JIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    JIntegral = JKIntegrals(a,b,c,d)
END FUNCTION JIntegral

REAL*8 FUNCTION KIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    KIntegral = JKIntegrals(a,d,c,b)
END FUNCTION KIntegral

REAL*8 FUNCTION G(a,b,P,NB)

    INTEGER :: a, b, NB, i, j
    REAL*8, dimension(NB,NB) :: P

END FUNCTION G

END MODULE DftOperations

