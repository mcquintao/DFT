!----------------------------------------------------------------------!
! MÓDULO COM AS PRINCIPAIS FUNÇÕES E SUBROTINAS DOS MÉTODOS DFT E HF   !
!                                                                      !
! SUBROTINAS PRESENTES:                                                !
!                                                                      !
!----------------------------------------------------------------------!

MODULE DftOperations

USE MathLib
CONTAINS

SUBROUTINE overlapMatrix(BASIS,NBASIS,SMAT)

      IMPLICIT NONE
      
      integer, intent(in) :: NBASIS
      integer :: i, j
      real*8, dimension(NBASIS), intent(in) :: BASIS
      real*8, dimension(NBASIS,NBASIS), intent(out) :: SMAT

      do i=1,NBASIS
        do j=1,NBASIS
            SMAT(i,j) = overlap(BASIS(i),BASIS(j))
        end do
      end do


END SUBROUTINE overlapMatrix

SUBROUTINE xOverlapMatrix(SMAT,NBASIS,SMAT_VEC,SMAT_VAL)

    real*8 :: error = 0.000001d0

    call jacobi(SMAT,SMAT_VEC,error,NBASIS)

END SUBROUTINE xOverlapMatrix


SUBROUTINE potentialMatrix(BASIS,NBASIS,UMAT)




END SUBROUTINE potentialMatrix


END MODULE DftOperations

