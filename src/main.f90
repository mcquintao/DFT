PROGRAM main

    USE DiskOperations
    USE DftOperations
    IMPLICIT NONE
    
    logical :: dbg, openShell
    integer :: ZATOM, CHARGE, NBASIS, NE, NALFA, NBETA
    real*8 :: TOTENERGY
    real*8, allocatable :: BASIS(:)
    real*8, allocatable :: SMAT(:,:), XMAT(:,:), TMAT(:,:), UMAT(:,:)
    real*8, allocatable :: HCORE(:,:), PMAT(:,:), CMAT(:,:), CMAT_ALFA(:,:), CMAT_BETA(:,:)

! Inicializar variáveis
    ZATOM = 0
    NE = 0
    CHARGE = 0
    NBASIS = 0
    dbg = .false.
    openShell = .false.
    TOTENERGY = 0.d0

! Carregar input
call openFiles()
call readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg,openShell)

NE = ZATOM - CHARGE
ALLOCATE(SMAT(1:NBASIS,1:NBASIS))
ALLOCATE(XMAT(1:NBASIS,1:NBASIS))
ALLOCATE(TMAT(1:NBASIS,1:NBASIS))
ALLOCATE(UMAT(1:NBASIS,1:NBASIS))
ALLOCATE(HCORE(1:NBASIS,1:NBASIS))
ALLOCATE(PMAT(1:NBASIS,1:NBASIS))
ALLOCATE(CMAT(1:NBASIS,1:NBASIS))

IF(openShell) THEN
    PRINT *, "** SCF - CAMADA ABERTA **"
    NALFA = NE/2
    NBETA = NALFA + MOD(NE,2)
    ALLOCATE(CMAT_ALFA(1:NBASIS,1:NBASIS))
    ALLOCATE(CMAT_BETA(1:NBASIS,1:NBASIS))
ELSE
    PRINT *, "** SCF - CAMADA FECHADA **"
END IF

if(dbg) then
    PRINT *, "DEBUG ON!"
    call dbgInput(ZATOM,CHARGE,NBASIS,BASIS,openShell,NE,NALFA,NBETA)
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

IF(openShell) THEN
    call scfOpen(XMAT,HCORE,BASIS,NBASIS,NE,PMAT,CMAT_ALFA,CMAT_BETA,TOTENERGY)
    PRINT *, "----------------------------------------------------"
    PRINT *, "OPEN SHELL TOTAL ENERGY (HARTREE): ", TOTENERGY
    PRINT *, "OPEN SHELL TOTAL ENERGY (KJ/MOL): ", TOTENERGY*2625.5
    PRINT *, "----------------------------------------------------"
ELSE
    call scfClose(XMAT,HCORE,BASIS,NBASIS,NE,PMAT,CMAT,TOTENERGY)
    PRINT *, "----------------------------------------------------"
    PRINT *, "CLOSE SHELL TOTAL ENERGY: ", TOTENERGY
    PRINT *, "CLOSE SHELL TOTAL ENERGY (KJ/MOL): ", TOTENERGY*2625.5
    PRINT *, "----------------------------------------------------"
END IF

! Fechar arquivos
call closeFiles()

END PROGRAM main

