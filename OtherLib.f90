!=======================================================================
!ESTA SUBROTINA FAZ A MULTIPLICACAO DE DUAS MATRIZES QUADRADAS C = A*B
!-----------------------------------------------------------------------
! PARAMETROS DE ENTRDA:
! A Matrize de entrada
! B Matrize de entrada
! N Dimensão da matriz (apenas matrizes quadradas)
! PARAMETROS DE SAIDA:
! C Matriz que armazenará a multiplicacao de A com B
!=======================================================================

SUBROUTINE MATMULT(SMAT_VAL,NEW_SMAT,XMAT,NBASIS)
IMPLICIT NONE
INTEGER :: NBASIS , i , j, k
REAL*8 : : SMAT_VAL(NBASIS,NBASIS),NEW_SMAT(NBASIS,NBASIS),XMAT(NBASIS,NBASIS)

XMAT = 0.0
DO i = 1 , N
  DO j = 1 , N
    DO k = 1 , N
      XMAT(i,j) = XMAT(i,j) + SMAT_VAL(i,k)∗NEW_SMAT(k,j)
    ENDDO
  ENDDO
ENDDO
END SUBROUTINE MATMULT

!=======================================================================
!ESTA SUBROTINA CALCULA A TRANSPOSTA DE UMA MATRIZ A
!-----------------------------------------------------------------------
! PARAMETROS DE ENTRDA:
! A Matrize de entrada
! N Dimensão da matriz (apenas matrizes quadradas)
! PARAMETROS DE SAIDA:
! C Matriz que armazenará a multiplicacao de A com B
!=======================================================================

SUBROUTINE MATRANS(A,N)
IMPLICIT NONE
INTEGER :: N, i, j
REAL*8 :: A(N,N), temp

DO i = 1,N
  DO j = 1,N
    IF(i<j) THEN
      temp = A(i,j)
      A(i,j) = A(j,i)
      A(j,i) = temp
    ENDIF
   ENDDO
ENDDO
END SUBROUTINE MATRANS

!=======================================================================
!ESTA SUBROTINA FAZ A SOMA DE DUAS MATRIZES QUADRADAS C = A + B
!-----------------------------------------------------------------------
! PARAMETROS DE ENTRDA:
! A Matrize de entrada
! B Matrize de entrada
! N Dimensão da matriz (apenas matrizes quadradas)
! PARAMETROS DE SAIDA:
! C Matriz que armazenará a soma de A com B
!=======================================================================

SUBROUTINE MATSUM(A,B,C,N)
IMPLICIT NONE
INTEGER :: N , i , j
REAL*8 : : A(N,N),B(N,N),C(N,N)

DO i = 1 , N
  DO j = 1 , N
    C(i,j) = A(i,j) + B(i,j)
  ENDDO
ENDDO
END SUBROUTINE MATSUM
