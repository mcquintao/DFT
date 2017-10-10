!----------------------------------------------------------------------!
! MÓDULO VOLTADO PARA OPERAÇÕES MATEMÁTICAS BÁSICAS                    !
!                                                                      !
! SUBROTINAS PRESENTES:                                                !
!                                                                      !
!----------------------------------------------------------------------!

MODULE MathLib
CONTAINS

FUNCTION overlap(alpha,beta)
    
    IMPLICIT NONE
    real*8, intent(in) :: alpha, beta
    real*8 :: C1, C2, PI, overlap
    
    PI = 4.D0*DATAN(1.D0)

    C1 = ((2.d0*alpha)/PI)**(3.d0/4.d0)
    C2 = ((2.d0*beta)/PI)**(3.d0/4.d0)

    overlap = C1*C2*(PI/(alpha+beta))**(3.d0/2.d0)


END FUNCTION overlap


END MODULE MathLib
