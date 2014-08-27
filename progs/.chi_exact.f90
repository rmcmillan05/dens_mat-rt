PROGRAM chi_exact
    USE double
    USE params , ONLY : omega_au
    USE num_lines
    IMPLICIT NONE

    REAL(KIND=DP) :: wp
    REAL(KIND=DP) :: wq
    COMPLEX(KIND=DP) :: X1 = (0.0_DP, 0.0_DP)
    COMPLEX(KIND=DP) :: X2 = (0.0_DP, 0.0_DP)
    INTEGER :: i, j, m, n, v, levels

    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:)   :: rho, ham, mu, gma
    CHARACTER(LEN=64) :: rho_in, ham_in, mu_in, gma_in, data_pref

    wp = omega_au
    wq = omega_au

    SELECT CASE ( IARGC() )
        CASE ( 0 )
            WRITE(*,*) 'No input folder detected.'
            CALL EXIT(0)
        CASE ( 1 )
            CALL GETARG(1, data_pref)
        CASE DEFAULT
            WRITE(*,*) 'Too many input arguments.'
            CALL EXIT(0)
    END SELECT

    rho_in=TRIM(data_pref)//'rho.txt'
    ham_in=TRIM(data_pref)//'ham.txt'
    gma_in=TRIM(data_pref)//'gma.txt'
    mu_in =TRIM(data_pref)//'mu.txt'

    levels = numlines(rho_in)
    ALLOCATE(rho(levels,levels), ham(levels,levels), mu(levels,levels), &
             gma(levels,levels))
    OPEN(UNIT=10, FILE=rho_in, STATUS='OLD')
        DO i = 1, levels
            READ(10, *) (rho(i, j), j=1,levels)
        ENDDO
    CLOSE(10)
    OPEN(UNIT=10, FILE=ham_in, STATUS='OLD')
        DO i = 1, levels
            READ(10, *) (ham(i, j), j=1,levels)
        ENDDO
    CLOSE(10)
    OPEN(UNIT=10, FILE=gma_in, STATUS='OLD')
        DO i = 1, levels
            READ(10, *) (gma(i, j), j=1,levels)
        ENDDO
    CLOSE(10)
    OPEN(UNIT=10, FILE=mu_in, STATUS='OLD')
        DO i = 1, levels
            READ(10, *) (mu(i, j), j=1,levels)
        ENDDO
    CLOSE(10)

    DO m = 1,levels
        DO n = 1,levels
            X1 = X1 + (rho(m,m)-rho(n,n))*mu(m,n)*mu(n,m) / &
                       (ham(n,n)-ham(m,m) - wp - ci*gma(n,m))
        ENDDO
    ENDDO 

    DO m = 1,levels
        DO n = 1,levels
            DO v = 1,levels
                X2 = X2 + (rho(m,m) - rho(v,v)) *         &
                        mu(m,n) * mu(n,v) * mu(v,m) / &
                        ( (ham(n,n)-ham(m,m) - wp - wq - ci*gma(n,m)) * &
                          (ham(v,v)-ham(m,m) - wp      - ci*gma(v,m)) ) &

                      - (rho(v,v) - rho(n,n)) *  &
                        mu(m,n) * mu(v,m) * mu(n,v) / &
                        ( (ham(n,n)-ham(m,m) - wp - wq - ci*gma(n,m)) * &
                           (ham(n,n)-ham(v,v) - wp     - ci*gma(n,v)) )
            ENDDO
        ENDDO
    ENDDO

    WRITE(*,*) 'chi_1 = ', X1
    WRITE(*,*) 'chi_2 = ', X2

END PROGRAM chi_exact


