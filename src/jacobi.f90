MODULE JACOBI

CONTAINS

! Subrotina Jacobi
!
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

                if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip)))  
     *          .and.(abs(d(iq))+g.eq.abs(d(iq))))then                                       
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
       pause 'muitas interacoes no JACOBI'

! A subrotina abaixo ordena os autovalores encontrados na Jacobi. 

102    CONTINUE
       Call EIGSRT(d,v,n,np)

       END  

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
       END

END MODULE JACOBI
