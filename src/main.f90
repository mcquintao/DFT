PROGRAM main

    USE DiskOperations
    IMPLICIT NONE
    
    logical :: dbg
    integer :: ZATOM, CHARGE, NBASIS
    real*8, allocatable :: BASIS(:)


! Inicializar variáveis
    ZATOM = 0
    CHARGE = 0
    NBASIS = 0
    BASIS = 0.d0
    dbg = .false.

! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

if(dbg) then
    PRINT *, "DEBUG ON!"
    call dbgInput(ZATOM,CHARGE,NBASIS,BASIS)
end if
! Abrir arquivos de escrita


! --- PROGRAMA PRINCIPAL --- !

! Calcular OVERLAP (SMAT)


! Calcular X (XMAT)


! Calcular E. Cinética (TMAT)



! Calcular E. Potencial (UMAT)


! Fechar arquivos
call closeFiles()

END PROGRAM main
