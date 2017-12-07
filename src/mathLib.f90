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



! Subrotina Jacobi
! Retirada do Numerical Recipes - Pag-460
!  
! Diagonaliza a matriz de Overlap (S) que nesta subrotina e o (A)
!
! Computa todos os autovalores e autovetores de uma matriz simetrica real A,
! que de tamanho NxN. No output, os elementos de A acima da diagonal sao desprezados.
! A matriz D retorna os autovalores de A em seus primeiros N elementos. V e uma matriz com 
! a mesma dimensao logica e fisica de A, na qual as colunas contem , no output, os autovetores normalizados de 
! A. 
! Nrot retorna o numero de rotacoes de Jacobi que foram requeridos.
!
! A => Matriz a ser diagonalizada.(NxN)
! D => Autovalores de A.
! V => Autovetores de A.
!

      Subroutine Jacobi(A,N,NP,D,V,Nrot)
      INTEGER n,np,nrot,NMAX
      REAL*8 A(np,np),d(np),V(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

!     Inicicializa  a matriz identidade

        DO ip=1,n
           DO iq=1,n
              V(ip,iq)=0
           ENDDO
           V(ip,ip)=1
        ENDDO
        
!     Inicializa B e D para a diagonal de A
 
        DO ip=1,n
           B(ip)=A(ip,ip)
           D(ip)=B(ip)

!     Este vetor ira acumular termos da forma tapq como na equacao da Pag 458
!     do Numerical Recipes - eq (11.1.14) 
 
           Z(ip)=0
        ENDDO
        nrot=0   
        DO i=1,50     
           sm=0.

!     Soma dos elementos fora da diagonal

           DO ip=1,n-1 
              DO   iq=ip+1,n 
                   sm=sm+abs(a(ip,iq)) 
              ENDDO  
           ENDDO
           if(sm.eq.0.) goto 102
             if(i.lt.4)then
             tresh=0.2*sm/n**2
             else
             tresh=0
             ENDif
           DO   ip=1,n-1
                DO  iq=ip+1,n
                    g=100.*abs(a(ip,iq))

!      Apos quatro varreduras,pular a rotacao 
!      se o elemento fora da diagonal for pequeno.

                if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then                                       
                 a(ip,iq)=0.
                else if(abs(a(ip,iq)).gt.tresh)then 
                h=d(iq)-d(ip)
                   if(abs(h)+g.eq.abs(h))then
                   t=a(ip,iq)/h   
                   else

!      Equacao (11.1.10) da Pag-457 do Numerical Recipes

                   theta=0.5*h/a(ip,iq)
                   t=1./(abs(theta)+sqrt(1.+theta**2))
                   if(theta.lt.0.)t=-t
                   ENDif 
                   c=1./sqrt(1+t**2)
                   s=t*c
                   tau=s/(1.+c)
                   h=t*a(ip,iq)
                   z(ip)=z(ip)-h
                   z(iq)=z(iq)+h
                   d(ip)=d(ip)-h
                   d(iq)=d(iq)+h
                   a(ip,iq)=0

!       Caso de rotacoes 1<=j< p

                   DO   j=1,ip-1
                   g=a(j,ip)
                   h=a(j,iq)
                   a(j,ip)=g-s*(h+g*tau)        
                   a(j,iq)=h+s*(g-h*tau)
                   ENDDO

!       Caso de rotacoes p < j < q

                   DO   j=ip+1,iq-1    
                        g=a(ip,j)
                        h=a(j,iq)
                        a(ip,j)=g-s*(h+g*tau)
                        a(j,iq)=h+s*(g-h*tau)
                   ENDDO

!       Caso de rotacoes q < j <= n

                   DO   j=iq+1,n
                        g=a(ip,j)
                        h=a(iq,j)
                        a(ip,j)=g-s*(h+g*tau) 
                        a(iq,j)=h+s*(g-h*tau)
                   ENDDO 
                   DO   j=1,n
                        g=v(j,ip)
                        h=v(j,iq)
                        v(j,ip)=g-s*(h+g*tau)
                        v(j,iq)=h+s*(g-h*tau)
                   ENDDO    
                        nrot=nrot+1
           ENDif    
           ENDDO
           ENDDO
           DO  ip=1,n
               b(ip)=b(ip)+z(ip)

!       Atualizar d com a soma de tapq e reinicializar z.
 
               d(ip)=b(ip)
               z(ip)=0
           ENDDO
       ENDDO 
       Stop 'muitas interacoes no JACOBI'

! A subrotina abaixo ordena os autovalores encontrados na Jacobi. 

102    CONTINUE
       Call EIGSRT(d,v,n,np)

       Return
       END Subroutine  

! Subrotina para ordenar os auovalores de Jacobi.(Ordem decrescente)

       Subroutine EIGSRT(d,v,n,np)
       Integer n,np
       Real * 8 d(np),v(np,np)
       Integer I,J,K
       Real * 8 P
       DO I=1,n-1
         K=I
         P=d(I)
         DO J=I+1,n
           If (d(J).le.P) then
               K=j        
               P=d(J)
           ENDIF
         ENDDO
         IF (K.ne.i) then
             d(K)=d(I)
             d(I)=P
             DO J=1,n
               P=v(J,I)
               v(J,I)=v(J,K)
               v(J,K)=P
             ENDDO
         ENDIF
       ENDDO
       RETURN
       END Subroutine
 


 



!SUBROUTINE JACOBI(a,x,abserr,n)
!
!
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
!
!
!IMPLICIT NONE
!INTEGER :: i, j, k, n
!DOUBLE PRECISION :: a(n,n),x(n,n)
!DOUBLE PRECISION :: abserr, b2, bar
!DOUBLE PRECISION :: beta, coeff, c, s, cs, sc
!
!! Inicialização com x(i,j)=0, x(i,i)=1
!x = 0.0
!DO i=1,n
!  x(i,i) = 1.0
!ENDDO
!
!! Encontrando os elemetos fora da diagonal
!b2 = 0.0
!DO i=1,n
!  DO j=1,n
!    IF (i.ne.j) b2 = b2 + a(i,j)**2
!  ENDDO
!ENDDO
!
!IF (b2 .LE. abserr) RETURN
!
!! média para os elementos da diagonal
!bar = 0.5*b2/float(n*n)
!
!DO WHILE (b2.gt.abserr)
!  DO i=1,n-1
!    DO j=i+1,n
!      IF (a(j,i)**2 <= bar) CYCLE 
!      b2 = b2 - 2.0*a(j,i)**2
!      bar = 0.5*b2/float(n*n)
!! Calcula os coeficientes c e s para a Matriz de  Givens
!      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
!      coeff = 0.5*beta/sqrt(1.0+beta**2)
!      s = sqrt(max(0.5+coeff,0.0))
!      c = sqrt(max(0.5-coeff,0.0))
!! Recalcula i e j
!      DO k=1,n
!        cs =  c*a(i,k)+s*a(j,k)
!        sc = -s*a(i,k)+c*a(j,k)
!        a(i,k) = cs
!        a(j,k) = sc
!      ENDDO
!! Nova matriz a_{k+1} de a_{k}, e autovetores 
!      DO k=1,n
!        cs =  c*a(k,i)+s*a(k,j)
!        sc = -s*a(k,i)+c*a(k,j)
!        a(k,i) = cs
!        a(k,j) = sc
!        cs =  c*x(k,i)+s*x(k,j)
!        sc = -s*x(k,i)+c*x(k,j)
!        x(k,i) = cs
!        x(k,j) = sc
!      ENDDO
!    ENDDO
!  ENDDO
!ENDDO
!
!END SUBROUTINE JACOBI

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

real*8 function JKIntegrals(a,b,c,d)

    REAL*8, intent(in) :: a, b, c, d
    REAL*8 :: PI
    PI = 4.D0*DATAN(1.D0)

    JKIntegrals = 16.d0/SQRT(PI)
    JKIntegrals = JKIntegrals*(a*b*c*d)**(3.d0/4.d0)
!   JKIntegrals = JKIntegrals/SQRT((a + c)*(b + d)*(a+b+c+d))
    JKIntegrals = JKIntegrals/((a + c)*(b + d)*SQRT(a+b+c+d))

end function JKIntegrals

END MODULE MathLib
