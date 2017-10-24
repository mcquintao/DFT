PROGRAM main

    USE DiskOperations
    USE DftOperations
    IMPLICIT NONE
    
    logical :: dbg
    integer :: ZATOM, CHARGE, NBASIS
    real*8, allocatable :: BASIS(:)
    real*8, allocatable :: SMAT(:,:), XMAT(:,:), TMAT(:,:), UMAT(:,:)
    real*8, allocatable :: HCORE(:,:), PMAT(:,:)

! Inicializar variáveis
    ZATOM = 0
    CHARGE = 0
    NBASIS = 0
    dbg = .false.

! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

ALLOCATE(SMAT(1:NBASIS,1:NBASIS))
ALLOCATE(XMAT(1:NBASIS,1:NBASIS))
ALLOCATE(TMAT(1:NBASIS,1:NBASIS))
ALLOCATE(UMAT(1:NBASIS,1:NBASIS))
ALLOCATE(HCORE(1:NBASIS,1:NBASIS))
ALLOCATE(PMAT(1:NBASIS,1:NBASIS))



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
CALL TransferenceMatrix(NBASIS,SMAT,XMAT)

IF(dbg) THEN
    CALL dbgTransference(XMAT,NBASIS)
ENDIF

! Calcular E. Cinética (TMAT)
call kineticMatrix(BASIS,NBASIS,TMAT)

if(dbg) then
    call dbgKinetic(TMAT,NBASIS)
end if

! Calcular E. Potencial (UMAT)
call potentialMatrix(ZATOM,BASIS,NBASIS,UMAT)

if(dbg) then
    call dbgPotencial(UMAT,NBASIS)
end if

! One Electron Matrix
call HCor(TMAT,UMAT,HCORE,NBASIS)

if(dbg) then
    call dbgHCor(HCORE,NBASIS)
end if

call PGuess(PMAT,NBASIS)

if(dbg) then
    call dbgPMat(PMAT,NBASIS)
end if


! Fechar arquivos
call closeFiles()

END PROGRAM main
