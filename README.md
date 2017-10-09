# DFT: Programa FORTRAN para cálculo da estrutura eletrônica de um átomo

em construção ... 

PROGRAMA DFT:

1. VARIÁVEIS e UNIT
    1.1 INTEGER:
        ZATOM,CHARGE,NBASIS
    1.2. REAL:
        BASIS(NBASIS)
    1.3. LOGICAL:
        dbg

2. LEITURA DO INPUT:
    2.1. openFiles()
    2.2. readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg)
3. Cálculos Iniciais:
    3.1. overlapMatrix(BASIS,NBASIS,SMAT)
    3.2. kineticMatrix(BASIS,NBASIS,TMAT)
    3.3. potentialMatrix(BASIS,NBASIS,ZATOM,UMAT)
    3.4. diagOvelap(SMAT,NBASIS,XMAT)



SUBROTINAS:
    openFiles(): src/diskOperations.f90
        * Abre os arquivos input(1), output(2) e debug(99)

    readInput(ZATOM,CHARGE,NBASIS,BASIS,dbg): src/diskOperations.f90
        * Lê os dados em input(1)


