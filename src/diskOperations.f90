!----------------------------------------------------------------------!
! MÓDULO VOLTADO PARA OPERAÇÕES NO DISCO, COMO ABERTURA, LEITURA E     !
! ESCRITA DE/EM ARQUIVOS                                               !
!                                                                      !
! SUBROTINAS PRESENTES:                                                !
!                                                                      !
! -> openFiles()                                                       !
! -> closeFiles()                                                      !
! -> readInput()                                                       !
!----------------------------------------------------------------------!



SUBROUTINE openFiles()

!   Subrotina para abrir arquivos de leitura e escrita
!   Matheus Campos Quintão
!   Última edição: 08/10/17

    IMPLICIT NONE
    logical :: fileExist

    ! Checar se o arquivo "input" existe!
    
    INQUIRE(FILE="input", EXIST=fileExist)
    
    if(fileExist) then  ! Se existir...
    
        OPEN(unit=1, file="input")  ! Abrir input em unit 1
        OPEN(unit=2, file="output") ! Abrir input em unit 2
        OPEN(unit=99, file="debug") ! Abrir input em unit 99
    
    else    ! se não existir ...
    
        PRINT *, "O ARQUIVO INPUT NÃO FOI ENCONTRADO"
        STOP
    
    end if

END SUBROUTINE openFiles


SUBROUTINE closeFiles()

!   Subrotina para fechar arquivos
!   Matheus Campos Quintão
!   Última edição: 08/10/17


close(1)
close(2)
close(99)

END SUBROUTINE closeFiles

SUBROUTINE readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

!   Subrotina para ler input
!   Matheus Campos Quintão
!   Última edição: 08/10/17


! dbg => opção de debug
! ZATOM, CHARGE, NBASIS => Número atômico, Carga e Número de bases


!INPUT EXEMPLO (SEM EXCLAMAÇÃO!!):

! 3 0 4        (ZATOM, CHARGE, NBASIS)
! 10           (BASE 1)
! 20           (BASE 2)
! 30           (BASE 3)
! 40           (BASE 4)
!



    IMPLICIT NONE
    logical, intent(out) :: dbg
    integer, intent(out) :: ZATOM, CHARGE, NBASIS
    integer :: i, j
    real*8, allocatable, intent(INOUT) :: BASIS(:)




! ler a primeira linha do input
    read(1,*) ZATOM, CHARGE, NBASIS

    
! Alocação dinâmica de memória para o vetor de bases
    ALLOCATE(BASIS(1:NBASIS))
    
! Se ZATOM for negativo, o DEBUG será ativado!
    if(ZATOM.lt.0) then
        ZATOM = ABS(ZATOM)
        dbg = .true.
    end if

! Lendo os elementos da base linha a linha até NBASIS    
    do i=1,NBASIS
        read(1,*) BASIS(i)
    end do


! SE DEBUG for VERDADEIRO
    if(dbg) then
    PRINT *, "DBG - readInput"    
    write(99,*) "------DADOS DO INPUT------"
    
    write(99,"(A18, I2)") "Número atômico: ", ZATOM
    write(99,*) "Carga: ", CHARGE
    write(99,"(A18, I2)") "Número de bases: ", NBASIS
    write(99,*) "BASE: "
    
    do i=1,NBASIS
        write(99,"(F9.4)") BASIS(i)
    end do
    write(99,*) "--------------------------"
    write(99,*)
    
    end if



END SUBROUTINE readInput
