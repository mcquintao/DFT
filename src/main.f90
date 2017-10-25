PROGRAM main

    USE DiskOperations
    USE DftOperations
    IMPLICIT NONE
    
    logical :: dbg
    integer :: ZATOM, CHARGE, NBASIS, NE
    real*8, allocatable :: BASIS(:)
    real*8, allocatable :: SMAT(:,:), XMAT(:,:), TMAT(:,:), UMAT(:,:)
    real*8, allocatable :: HCORE(:,:), PMAT(:,:), CMAT(:,:)

! Inicializar variáveis
    ZATOM = 0
    NE = 0
    CHARGE = 0
    NBASIS = 0
    dbg = .false.

! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)

NE = ZATOM - CHARGE
ALLOCATE(SMAT(1:NBASIS,1:NBASIS))
ALLOCATE(XMAT(1:NBASIS,1:NBASIS))
ALLOCATE(TMAT(1:NBASIS,1:NBASIS))
ALLOCATE(UMAT(1:NBASIS,1:NBASIS))
ALLOCATE(HCORE(1:NBASIS,1:NBASIS))
ALLOCATE(PMAT(1:NBASIS,1:NBASIS))
ALLOCATE(CMAT(1:NBASIS,1:NBASIS))



if(dbg) then
    PRINT *, "DEBUG ON!"
    call dbgInput(ZATOM,CHARGE,NBASIS,BASIS)
end if

! --- PROGRAMA PRINCIPAL --- !

! Calcular OVERLAP (SMAT)
call overlapMatrix(BASIS,NBASIS,SMAT)

if(dbg) then
    !call dbgOverlap(SMAT,NBASIS,"TESTE TITULO",12)
    call dbgMatrix(SMAT,NBASIS,"OVERLAP MATRIX",14)
end if

! Calcular X (XMAT)
CALL TransferenceMatrix(NBASIS,SMAT,XMAT)

IF(dbg) THEN
    !CALL dbgTransference(XMAT,NBASIS)
    call dbgMatrix(XMAT,NBASIS,"TRANSFERENCE MATRIX",19)
ENDIF

! Calcular E. Cinética (TMAT)
call kineticMatrix(BASIS,NBASIS,TMAT)

if(dbg) then
    !call dbgKinetic(TMAT,NBASIS)
    call dbgMatrix(TMAT,NBASIS,"KINETIC MATRIX",14)
end if

! Calcular E. Potencial (UMAT)
call potentialMatrix(ZATOM,BASIS,NBASIS,UMAT)

if(dbg) then
    !call dbgPotencial(UMAT,NBASIS)
    call dbgMatrix(UMAT,NBASIS,"POTENTIAL MATRIX",16)
end if

! One Electron Matrix
call HCor(TMAT,UMAT,HCORE,NBASIS)

if(dbg) then
    !call dbgHCor(HCORE,NBASIS)
    call dbgMatrix(HCORE,NBASIS,"H_CORE MATRIX",13)
end if

call PGuess(PMAT,NBASIS)

if(dbg) then
    !call dbgPMat(PMAT,NBASIS)
    call dbgMatrix(PMAT,NBASIS,"DENSITY MATRIX",14)
end if

call scf(XMAT,HCORE,BASIS,NBASIS,NE,PMAT,CMAT)

! Fechar arquivos
call closeFiles()

END PROGRAM main
