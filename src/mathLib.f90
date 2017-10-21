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

!=======================================================================
! ESTA SUBROTINA CALCULA OS AUTOVALORES E AUTOVETORES PELO
! METODO DE JACOBI
!-----------------------------------------------------------
! PARAMETROS DE ENTRADA:
! SMAT - Matriz de Overlap
! NBASIS      - Numero de elementos da matriz (Numeor de bases)
! abserr - tolerancia absoluta [Soma de elementos fora da diagonal ao quadrado]
! PARAMETROS DE SAIDA:
! SMAT_VAL - autovalor
! SMAT_VEC - autovetor
!=======================================================================

SUBROUTINE JACOBI(a,x,abserr,n)


!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================


IMPLICIT NONE
INTEGER :: i, j, k, n
DOUBLE PRECISION :: a(n,n),x(n,n)
DOUBLE PRECISION :: abserr, b2, bar
DOUBLE PRECISION :: beta, coeff, c, s, cs, sc

! Inicialização com x(i,j)=0, x(i,i)=1
x = 0.0
DO i=1,n
  x(i,i) = 1.0
ENDDO

! Encontrando os elemetos fora da diagonal
b2 = 0.0
DO i=1,n
  DO j=1,n
    IF (i.ne.j) b2 = b2 + a(i,j)**2
  ENDDO
ENDDO

IF (b2 .LE. abserr) RETURN

! média para os elementos da diagonal
bar = 0.5*b2/float(n*n)

DO WHILE (b2.gt.abserr)
  DO i=1,n-1
    DO j=i+1,n
      IF (a(j,i)**2 <= bar) CYCLE 
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! Calcula os coeficientes c e s para a Matriz de  Givens
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! Recalcula i e j
      DO k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      ENDDO
! Nova matriz a_{k+1} de a_{k}, e autovetores 
      DO k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      ENDDO
    ENDDO
  ENDDO
ENDDO

END SUBROUTINE JACOBI

FUNCTION kinetic(alpha,beta)

    
    IMPLICIT NONE
    real*8, intent(in) :: alpha, beta
    real*8 :: C1, C2, PI, kinetic
    
    PI = 4.D0*DATAN(1.D0)

    kinetic = 6.d0*sqrt(2.d0)*((alpha**3.d0)*(beta**7.d0))**(1.d0/4.d0)
    kinetic = kinetic/((alpha+beta)**(3.d0/2.d0))
    kinetic = kinetic*(1 - (beta/(alpha + beta)))



END FUNCTION kinetic



FUNCTION intgPot(ZATOM,alpha,beta) 

    IMPLICIT NONE 
    real*8, intent(in) :: alpha, beta 
    real*8 :: PI, C1, C2, intgPot 
    integer, intent(in) :: ZATOM 
    PI = 4.D0*DATAN(1.D0) 
    C1 = ((2.d0*alpha)/PI)**(3.d0/4.d0) 
    C2 = ((2.d0*beta)/PI)**(3.d0/4.d0) 

    intgPot =-2.d0*(ZATOM*PI/(alpha+beta))*C1*C2 
 
END FUNCTION intgPot

real*8 function MULTMAT(A,B,n)


    implicit none
    real*8, dimension(n,n), intent(in) :: A, B
    dimension  MULTMAT(n,n)
    integer :: i, j, k,n
    
    MULTMAT = 0.d0
    do i=1,n
        do j=1,n
            do k=1,n
                MULTMAT(i,j) = MULTMAT(i,j) + A(i,k)*B(k,j)
            end do
        end do
    end do
    
end function MULTMAT

END MODULE MathLib
