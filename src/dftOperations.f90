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
INTEGER :: i, j, Nrot
INTEGER, INTENT(IN) :: NBASIS
REAL*8, DIMENSION(NBASIS,NBASIS), INTENT(IN) :: SMAT
REAL*8, DIMENSION(NBASIS,NBASIS), INTENT(OUT) :: XMAT
REAL*8, DIMENSION(NBASIS,NBASIS) :: LOCAL_SMAT, SMAT_VEC
REAL*8, DIMENSION(NBASIS) :: EVALUES
REAL*8 WALBER

LOCAL_SMAT = SMAT
XMAT = 0.0D0
SMAT_VEC = 0.0D0
WALBER = 0.0D0

! CALL JACOBI(LOCAL_SMAT,SMAT_VEC,Tol,NBASIS)
call Jacobi(LOCAL_SMAT,NBASIS,NBASIS,EVALUES,SMAT_VEC,Nrot)

DO j = 1,NBASIS
!    WALBER = SQRT(LOCAL_SMAT(j,j))  ! AUTO VALORES
    WALBER = SQRT(EVALUES(j))  ! AUTO VALORES
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
    H = T + U

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

SUBROUTINE newDensityClose(CMAT,P,NBASIS,NE)
    INTEGER, intent(in) :: NBASIS, NE
    INTEGER :: i, j, a
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: CMAT
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: P

    P = 0.d0
    DO i=1,NBASIS
        DO j=1,NBASIS
            IF(NE.eq.1) THEN
                P(i,j) = P(i,j) + 1.d0 * CMAT(i,1)*CMAT(j,1)  !! PAG. 139 EQ:3.145; SZABO (CONFERIR TRANSPOSTA)
            ELSE
                DO a=1,NE/2
                    P(i,j) = P(i,j) + 2.d0 * CMAT(i,a)*CMAT(j,a)  !! PAG. 139 EQ:3.145; SZABO (CONFERIR TRANSPOSTA)
                END DO
            END IF
        END DO
    END DO

END SUBROUTINE newDensityClose

SUBROUTINE newDensityOpen(CMAT,P,NBASIS,NE)
    INTEGER, intent(in) :: NBASIS, NE
    INTEGER :: i, j, a
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: CMAT
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: P
    
    P = 0.d0
    DO i=1,NBASIS
        DO j=1,NBASIS
                DO a=1,NE
                    P(i,j) = P(i,j) + 1.d0 * CMAT(i,a)*CMAT(j,a)  !! PAG. 139 EQ:3.145; SZABO (CONFERIR TRANSPOSTA)
                END DO
        END DO
    END DO

END SUBROUTINE newDensityOpen

REAL*8 FUNCTION ENERGY(H,F,P,NB)
    
    INTEGER, intent(in) :: NB
    INTEGER :: i, j
    REAL*8, dimension(NB,NB), intent(in) :: H, F, P
    REAL*8, dimension(NB,NB) :: HP, FP
   
    HP = MATMUL(H,P)
    FP = MATMUL(F,P)
    ENERGY = 0.d0
    DO i = 1, NB
            ENERGY = ENERGY + 0.5D0 * (HP(i,i) + FP(i,i))
    END DO

END FUNCTION ENERGY

REAL*8 FUNCTION ENERGYOPEN(H,FA,FB,PA,PB,NB)
    INTEGER, intent(in) :: NB
    INTEGER :: i, j
    REAL*8, dimension(NB,NB), intent(in) :: H, FA, FB, PA, PB
    REAL*8, dimension(NB,NB) :: P
    REAL*8, dimension(NB,NB) :: HP, FAPA, FBPB
    
    P = PA + PB
    HP = MATMUL(H,P)
    FAPA = MATMUL(FA,PA)
    FBPB = MATMUL(FB,PB)

    ENERGYOPEN = 0.d0
    DO i = 1, NB
            ENERGYOPEN = ENERGYOPEN + 0.5D0 * (HP(i,i) + FAPA(i,i) + FBPB(i,i)) 
    END DO

END FUNCTION ENERGYOPEN

subroutine newGMAT(G,P,BASIS,NB,NE,XHF)
CHARACTER(3), intent(in) :: XHF
INTEGER, intent(in) :: NB, NE
REAL*8, dimension(NB,NB), intent(in) :: P
REAL*8, dimension(NB), intent(in) :: BASIS
REAL*8, dimension(NB,NB), intent(out) :: G
REAL*8, dimension(NB,NB) :: JMAT, KMAT

INTEGER :: i,j,k,l
REAL*8 :: a,b,c,d,m
JMAT = 0.d0
KMAT = 0.d0
IF(NE.eq.1) THEN
    m = 1.0d0
ELSE
    m = 0.5d0
END IF

DO i=1,NB
    a = BASIS(i)
    DO j=1,NB
        b = BASIS(j)
        DO k=1,NB
            c = BASIS(k)
            DO l=1,NB
                d = BASIS(l)
                JMAT(i,j) = JMAT(i,j) + P(k,l)*JIntegral(a,b,c,d)
                KMAT(i,j) = KMAT(i,j) - m*P(k,l)*KIntegral(a,b,c,d)
            END DO
        END DO
        G(i,j) =JMAT(i,j) + KMAT(i,j)
    END DO
END DO
call dbgMatrix(JMAT,NB,"JMAT",4)
call dbgMatrix(KMAT,NB,"KMAT",4)
call dbgMatrix(JMAT + KMAT,NB,"J+KMAT",6)
PRINT *, "J e K: ", expectValue(P,JMAT,NB), expectValue(P,KMAT,NB)
PRINT *, "J + K: ", expectValue(P,JMAT,NB) + expectValue(P,KMAT,NB)

END SUBROUTINE newGMAT


! SUBROTINA SCF!! -----------------------------------

SUBROUTINE SCFCLOSE(XMAT,HCORE,SMAT,BASIS,NBASIS,NE,mix,PMAT,CMAT,GMAT,TotEnergy)
    REAL*8, PARAMETER :: Tol = 1.0E-12
    INTEGER, intent(in) :: NBASIS, NE
    INTEGER :: i,j,k,l,loop, Nrot
    REAL*8, dimension(NBASIS), intent(in) :: BASIS
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: HCORE, SMAT
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: PMAT, CMAT, GMAT
    REAL*8, dimension(NBASIS,NBASIS) :: XMAT, FMAT, XFMAT, NewPMAT
    REAL*8, DIMENSION(NBASIS) :: EVALUES
    REAL*8 :: a,b,c,d, oldenergy
    REAL*8, intent(in) :: mix
    REAL*8, intent(out) :: TotEnergy

    TotEnergy = 0.d0
    NewPMAT = 1.d0
    OldEnergy = 1000.00
    trPG = 0.d0

    DO loop=1,100
    
    PRINT *, "SCF Cycle ", loop

    XFMAT   = 0.d0
    GMAT    = 0.d0
    FMAT    = 0.d0
    CMAT    = 0.d0

    call newGMAT(GMAT,PMAT,BASIS,NBASIS,NE,"RHF")

    call dbgMatrix(GMAT,NBASIS,"G MATRIX",8)
    FMAT = HCORE + GMAT

    TotEnergy = ENERGY(HCORE,FMAT,PMAT,NBASIS)

    PRINT*, "##", loop, TotEnergy, expectValue(PMAT,FMAT,NBASIS)
    WRITE(99,*) "Energia Total:", TotEnergy    
    PRINT *, "--------------------------------"
    PRINT *, "TOTAL ENERGY: ", TotEnergy
    PRINT *, "--------------------------------"

!    IF(ABS(densityConv(PMAT,newPMAT,NBASIS)).lt.10E-12) THEN
!        GOTO 150
!    END IF

    IF(ABS(TotEnergy - OldEnergy).lt.10E-6) THEN
        GOTO 150
    END IF

    OldEnergy = TotEnergy



call dbgMatrix(FMAT,NBASIS,"F MATRIX",8)

    XFMAT = FMAT
    call TransfFMAT(XFMAT,XMAT,NBASIS)  ! XFMAT = XMAT' * FMAT * XMAT
    
call dbgMatrix(XFMAT,NBASIS,"XF MATRIX",9)

    call FLUSH(99)
    
!    call JACOBI(XFMAT,CMAT,Tol,NBASIS) ! XFMAT -> ENERGY(A. VAL); CMAT -> XCOEFF(A. VEC)
     call Jacobi(XFMAT,NBASIS,NBASIS,EVALUES,CMAT,Nrot) 
    
     write(99,*) "auto valores" 
     write(99,*) (Evalues(i),i=1,NBASIS) 
     XFMAT = 0.0
     DO i = 1, NBASIS
        XFMAT(i,i) = EVALUES(i)
     ENDDO


!call dbgMatrix(XFMAT,NBASIS,"ENERGY MATRIX",13)
call dbgMatrix(CMAT,NBASIS,"XCOEFF MATRIX",13)

    CMAT = MATMUL(XMAT,CMAT)

call dbgMatrix(CMAT,NBASIS,"COEFF MATRIX",12)

    call newDensityClose(CMAT,NewPMAT,NBASIS,NE)

call dbgMatrix(NewPMAT,NBASIS,"NEW DENSITY MATRIX",18)

    if (loop < 2) then
      PMAT = NewPMAT
    else
      PMAT = mix*PMAT + (1 - mix)*NEWPMAT
    ENDIF
    

    END DO

100 PRINT *, "O CICLO SCF NÃO CONVERGIU"
GOTO 200

150 PRINT *, "O CICLO SCF CONVERGIU!"
GOTO 200

200 PRINT *, "FIM DO SCF"

END SUBROUTINE SCFCLOSE

SUBROUTINE SCFOPEN(XMAT,HCORE,BASIS,NBASIS,NE,mix,PMAT,CMAT_ALFA,CMAT_BETA,GMAT,Totenergy)
    REAL*8, PARAMETER :: Tol = 1.0E-12
    INTEGER, intent(in) :: NBASIS 
    INTEGER :: i,j,k,l,loop, Nrot, NALFA, NBETA
    REAL*8, dimension(NBASIS), intent(in) :: BASIS
    REAL*8, dimension(NBASIS,NBASIS), intent(in) :: HCORE
    REAL*8, dimension(NBASIS,NBASIS), intent(out) :: PMAT, CMAT_ALFA, CMAT_BETA, GMAT
    REAL*8, dimension(NBASIS,NBASIS) :: XMAT, FMAT_ALFA, FMAT_BETA, GALFA, GBETA, XFMAT_ALFA, XFMAT_BETA, PALFA, PBETA
    REAL*8, dimension(NBASIS,NBASIS) :: NewPMAT_ALFA, NewPMAT_BETA
    REAL*8, DIMENSION(NBASIS) :: EVALUES_ALFA, EVALUES_BETA
    REAL*8 :: a,b,c,d, oldenergy
    REAL*8, intent(in) :: mix
    REAL*8, intent(out) :: TotEnergy

    TotEnergy = 0.d0
    NewPMAT_ALFA = 0.d0
    NewPMAT_BETA = 0.d0
    OldEnergy = 1000.00
    GMAT = 0.d0

    NALFA = NE/2
    NBETA = NALFA + MOD(NE,2)
    PALFA = PMAT
    PBETA = PMAT

    DO loop=1,100
    
    PRINT *, "SCF Cycle ", loop

    XFMAT_ALFA   = 0.d0
    XFMAT_BETA   = 0.d0
    GALFA        = 0.d0
    GBETA        = 0.d0
    FMAT_ALFA    = 0.d0
    FMAT_BETA    = 0.d0

    DO i=1,NBASIS
        a = BASIS(i)
        DO j=1,NBASIS
            b = BASIS(j)
            DO k=1,NBASIS
                c = BASIS(k)
                DO l=1,NBASIS
                    d = BASIS(l)
                    GALFA(i,j) = GALFA(i,j) + PMAT(k,l)*(JIntegral(a,b,c,d))  - PALFA(k,l)*(KIntegral(a,b,c,d))
                    GBETA(i,j) = GBETA(i,j) + PMAT(k,l)*(JIntegral(a,b,c,d))  - PBETA(k,l)*(KIntegral(a,b,c,d))
                END DO
            END DO
        END DO
    END DO
    
    GMAT = 0.5d0*(GALFA + GBETA)
call dbgMatrix(GMAT,NBASIS,"G MATRIX",8)

    FMAT_ALFA = HCORE + GALFA
    FMAT_BETA = HCORE + GBETA

    TotEnergy = ENERGYOPEN(HCORE,FMAT_ALFA,FMAT_BETA,NewPMAT_ALFA,NewPMAT_BETA,NBASIS)

    WRITE(99,*) "Energia Total:", TotEnergy    

    PRINT *, "--------------------------------"
    PRINT *, "TOTAL ENERGY: ", TotEnergy
    PRINT *, "--------------------------------"

    IF(ABS(TotEnergy - OldEnergy).lt.10E-6) THEN
        GOTO 160
    END IF
    OldEnergy = TotEnergy

    

!call dbgMatrix(FMAT,NBASIS,"F MATRIX",8)

    XFMAT_ALFA = FMAT_ALFA
    XFMAT_BETA = FMAT_BETA
    call TransfFMAT(XFMAT_ALFA,XMAT,NBASIS)  ! XFMAT = XMAT' * FMAT * XMAT
    call TransfFMAT(XFMAT_BETA,XMAT,NBASIS)  ! XFMAT = XMAT' * FMAT * XMAT
    
!call dbgMatrix(XFMAT,NBASIS,"XF MATRIX",9)

    call FLUSH(99)
    
!    call JACOBI(XFMAT,CMAT,Tol,NBASIS) ! XFMAT -> ENERGY(A. VAL); CMAT -> XCOEFF(A. VEC)

     call Jacobi(XFMAT_ALFA,NBASIS,NBASIS,EVALUES_ALFA,CMAT_ALFA,Nrot) 
     call Jacobi(XFMAT_BETA,NBASIS,NBASIS,EVALUES_BETA,CMAT_BETA,Nrot) 
    
!     write(99,*) "auto valores" 
!     write(99,*) (Evalues(i),i=1,NBASIS) 
     XFMAT_ALFA = 0.0
     XFMAT_BETA = 0.0
     do i = 1, NBASIS
        XFMAT_ALFA(i,i) = EVALUES_ALFA(i)
        XFMAT_BETA(i,i) = EVALUES_BETA(i)
     ENDDO


!call dbgMatrix(XFMAT,NBASIS,"ENERGY MATRIX",13)
!call dbgMatrix(CMAT,NBASIS,"XCOEFF MATRIX",13)

    CMAT_ALFA = MATMUL(XMAT,CMAT_ALFA)
    CMAT_BETA = MATMUL(XMAT,CMAT_BETA)


!call dbgMatrix(CMAT,NBASIS,"COEFF MATRIX",12)

    call newDensityOpen(CMAT_ALFA,NewPMAT_ALFA,NBASIS,NALFA)
    call newDensityOpen(CMAT_BETA,NewPMAT_BETA,NBASIS,NBETA)

!call dbgMatrix(NewPMAT,NBASIS,"NEW DENSITY MATRIX",18)

    if (loop < 2) then
      PALFA = NewPMAT_ALFA
      PBETA = NewPMAT_BETA
      PMAT = PALFA + PBETA
    else
      PALFA = mix*PALFA + (1 - mix)*NewPMAT_ALFA
      PBETA = mix*PBETA + (1 - mix)*NewPMAT_BETA
      PMAT = PALFA + PBETA
    ENDIF

    END DO

    PRINT *, "O CICLO SCF NÃO CONVERGIU"
GOTO 210

160 PRINT *, "O CICLO SCF CONVERGIU!"
GOTO 210

210 PRINT *, "FIM DO SCF"
END SUBROUTINE SCFOPEN

! FUNÇÕES ---------------------------------------------

REAL*8 FUNCTION JIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    JIntegral = JKIntegrals(a,b,c,d)
END FUNCTION JIntegral

REAL*8 FUNCTION KIntegral(a,b,c,d)
    REAL*8, intent(in) :: a,b,c,d
    KIntegral = JKIntegrals(a,d,c,b)
END FUNCTION KIntegral

REAL*8 FUNCTION expectValue(PMAT,MAT,NB)
    INTEGER, INTENT(in) :: NB
    INTEGER :: i
    REAL*8, DIMENSION(NB,NB) :: PM
    REAL*8, DIMENSION(NB,NB), intent(in) :: PMAT,MAT

    expectValue = 0.d0
    PM = MATMUL(MAT,PMAT)
    DO i=1,NB
        expectValue = expectValue + PM(i,i)
    END DO
END FUNCTION

REAL*8 FUNCTION mulliken(PMAT,SMAT,ZATOM,NB)
    INTEGER, INTENT(in) :: NB, ZATOM
    INTEGER :: i
    REAL*8, DIMENSION(NB,NB) :: PS
    REAL*8, DIMENSION(NB,NB), intent(in) :: PMAT,SMAT

    mulliken = 0.d0
    PS = MATMUL(PMAT,SMAT)
    DO i=1,NB
        mulliken = mulliken + PS(i,i)
    END DO
    mulliken = ZATOM - mulliken

END FUNCTION


REAL*8 FUNCTION checkOrto(MAT,SMAT,NB)
    INTEGER, INTENT(in) :: NB
    INTEGER :: i
    REAL*8, DIMENSION(NB,NB) :: CSC
    REAL*8, DIMENSION(NB,NB), intent(in) :: MAT,SMAT

    checkOrto = 0.d0
    CSC = MATMUL(TRANSPOSE(MAT),SMAT)
    CSC = MATMUL(CSC,MAT)
    DO i=1,NB
        checkOrto = checkOrto + CSC(i,i)
    END DO
    checkOrto = checkOrto / NB
END FUNCTION

REAL*8 FUNCTION HCORENERGY(HCOR,NE,NB)
    INTEGER, INTENT(in) :: NB, NE
    INTEGER :: i
    REAL*8, DIMENSION(NB,NB), intent(in) :: HCOR
    HCORENERGY = 0.d0
    DO i=1,NE
        HCORENERGY = HCORENERGY +  2.d0 * HCOR(i,i)
    END DO

END FUNCTION

END MODULE DftOperations

