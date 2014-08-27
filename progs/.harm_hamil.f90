!
! Use to generate the hamiltonian, mu, and other matrices required
! for rho_prop.exe. Output is saved in specified folder.
!
! To call:: harm_hamil.exe [N] [output_directory]
!
! N is the number of energy levels used.
!
! output_directory [default=./data/harmonic] is where the matrix
! text files are saved.
!
PROGRAM harm_hamil
    USE double
    USE print_mat_mod
    USE params, ONLY : harm_k
    IMPLICIT NONE

    CHARACTER(LEN=64) :: data_pref='./DATA/data_harmonic/'
    CHARACTER(LEN=64) :: N_in='5'
    
    REAL(KIND=DP), ALLOCATABLE :: ham_mat(:)
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_mat(:, :)
    COMPLEX(KIND=DP), ALLOCATABLE :: mu_mat(:, :)

    REAL(KIND=DP), PARAMETER :: omega = SQRT(harm_k)
    REAL(KIND=DP) :: c
    INTEGER       :: N
    INTEGER       :: i, j


    SELECT CASE ( IARGC() )
        CASE ( 2 )
            CALL GETARG(1, N_in)
            CALL GETARG(2, data_pref)
        CASE ( 1 )
            CALL GETARG(1, N_in)
    END SELECT

    READ(N_in,*) N
    ALLOCATE(rho_mat(N,N), ham_mat(N), mu_mat(N,N))

    rho_mat = 0.0_DP
    ham_mat = 0.0_DP

    rho_mat(1,1) = (1.0_DP, 0.0_DP)

    DO i = 0, N-1
        ham_mat(i+1) = omega * (0.5_DP + REAL(i,KIND=DP))
    ENDDO

    c = 1.0_DP / SQRT(2.0_DP * omega)

    DO i = 0, N-1
        DO j = 0, N-1
            mu_mat(i+1, j+1) = c * (SQRT(REAL(i,KIND=DP)) * kronecker(j, i-1)+&
                                    SQRT(1.0_DP+i)        * kronecker(j, i+1) )
        ENDDO
    ENDDO

    OPEN(UNIT=10, FILE=TRIM(data_pref)//'ham.txt', STATUS='REPLACE')
        DO i = 1, N
            WRITE(10,'(ES16.8)') ham_mat(i)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=TRIM(data_pref)//'mu.txt', STATUS='REPLACE')
    CALL print_mat(mu_mat,10)
    CLOSE(10)

    OPEN(UNIT=10, FILE=TRIM(data_pref)//'rho.txt', STATUS='REPLACE')
    CALL print_mat(rho_mat,10)
    CLOSE(10)

    OPEN(UNIT=10, FILE=TRIM(data_pref)//'gma.txt', STATUS='REPLACE')
        DO i = 1,N
            DO j=1,N
                WRITE(10,'(ES16.8)',ADVANCE='NO') 0.0_DP
            ENDDO
            WRITE(10,*) 
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=TRIM(data_pref)//'big_gma.txt', STATUS='REPLACE')
        DO i = 1,N
            DO j=1,N
                WRITE(10,'(ES16.8)',ADVANCE='NO') 0.0_DP
            ENDDO
            WRITE(10,*) 
        ENDDO
    CLOSE(10)
CONTAINS

FUNCTION kronecker(i, j)
    USE double
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    REAL(KIND=DP)       :: kronecker

    IF ( i == j ) THEN
        kronecker = 1.0_DP
    ELSE
        kronecker = 0.0_DP
    ENDIF

END FUNCTION kronecker

END PROGRAM harm_hamil
