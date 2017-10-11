PROGRAM main

    USE DiskOperations
    USE DftOperations
    IMPLICIT NONE
    
    logical :: dbg
    integer :: ZATOM, CHARGE, NBASIS
    real*8, allocatable :: BASIS(:)
    real*8, allocatable :: SMAT(:,:)

! Inicializar variáveis
    ZATOM = 0
    CHARGE = 0
    NBASIS = 0
    BASIS = 0.d0
    dbg = .false.

! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

ALLOCATE(SMAT(1:NBASIS,1:NBASIS))

if(dbg) then
    PRINT *, "DEBUG ON!"
    call dbgInput(ZATOM,CHARGE,NBASIS,BASIS)
end if

! --- PROGRAMA PRINCIPAL --- !

! Calcular OVERLAP (SMAT)
call overlapMatrix(BASIS,NBASIS,SMAT)

if(dbg) then
    call dbgOverlap(SMAT,NBASIS)
end if

! Calcular X (XMAT)

! INPUT: SMAT(NBASIS,NBASIS), NBASIS
! OUTPUT: XMAT_VEC(NBASIS,NBASIS),XMAT_VAL(NBASIS,NBASIS)

! WALBER: call xOverlapMatrix(SMAT,NBASIS,SMAT_VEC,SMAT_VAL)


! Calcular E. Cinética (TMAT)



! Calcular E. Potencial (UMAT)
! GABRIEL: call potentialMatrix(BASIS,NBASIS,UMAT)

! Fechar arquivos
call closeFiles()

END PROGRAM main
