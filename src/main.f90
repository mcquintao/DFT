PROGRAM main


    IMPLICIT NONE
    
    logical :: dbg
    integer :: ZATOM, CHARGE, NBASIS
    real*8, allocatable :: BASIS(:)


! Inicializar variáveis


! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

! Abrir arquivos de escrita


! --- PROGRAMA PRINCIPAL --- !

! Calcular OVERLAP (SMAT)


! Calcular X (XMAT)


! Calcular E. Cinética (TMAT)



! Calcular E. Potencial (UMAT)


! Fechar arquivos
call closeFiles()

END PROGRAM main
