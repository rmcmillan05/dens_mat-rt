MODULE runge_mod
    USE double
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: runge

CONTAINS

SUBROUTINE runge
    USE double
    USE params , ONLY : num_lev, field, npts, rk_step, positions, mu, rho_0,  &
                        out_file, npos
    USE fields
    IMPLICIT NONE

    ! rho(t)
    COMPLEX(KIND=DP), ALLOCATABLE                 :: rho(:,:)
    ! RK variables
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: k1, k2, k3, k4
    ! Dummy index variables
    INTEGER                                       :: i, j, n, m 
    ! Time
    REAL(KIND=DP)                                 :: t
    ! Commuted matrix
    COMPLEX(KIND=DP), ALLOCATABLE                 :: comm(:,:)
    ! Dipole at each time-step
    REAL(KIND=DP)                                 :: dipole
    !
    ! Title format specifier
    CHARACTER(LEN=64)                             :: charfmat = '(ES16.8)'
    ! Title for rho element column
    CHARACTER(LEN=2)                              :: poschar
    ! Title for field column
    CHARACTER(LEN=16)                             :: fieldchar
    ! File id for output
    INTEGER                                       :: out_id=50
    !
    ! Percentage increment in screen update
    REAL, PARAMETER                               :: pfraco=5.0 
    ! Percentage increment in screen update
    REAL                                          :: pfrac=0.0 
    ! Percentage complete
    REAL                                          :: pcomp

    ! Write output to file
    OPEN(UNIT=out_id, FILE=out_file, STATUS='REPLACE')

    ALLOCATE(                                                                 &
             rho(num_lev,num_lev),                                            &
             comm(num_lev,num_lev),                                           &
             k1(num_lev,num_lev),                                             &
             k2(num_lev,num_lev),                                             &
             k3(num_lev,num_lev),                                             &
             k4(num_lev,num_lev)                                              &
             )
    
    ! Initializing
    rho    = rho_0
    t      = 0.0_DP
    dipole = 0.0_DP

    ! Printing column headers
    WRITE(out_id, '(A16)', ADVANCE='NO') ' t              '
    fieldchar=' '//TRIM(field)//'(t)'
    WRITE(out_id, '(A16)', ADVANCE='NO') fieldchar
    WRITE(out_id, '(A16)', ADVANCE='NO') ' dipole         '
    DO j = 1, npos 
        WRITE(out_id, '(A5)', ADVANCE='NO') ' rho_'
        WRITE(poschar, '(I2)') positions(j,1)
        WRITE(out_id, '(A3)', ADVANCE='NO') poschar//','
        WRITE(poschar, '(I2)') positions(j,2)
        WRITE(out_id, '(A2)', ADVANCE='NO') poschar
        WRITE(out_id, '(A22)', ADVANCE='NO') ''
    ENDDO
    WRITE(out_id,*)

    ! Printing values at t=0
    WRITE(out_id, charfmat, ADVANCE='NO') t
    WRITE(out_id, charfmat, ADVANCE='NO') efield(field, t)
    WRITE(out_id, charfmat, ADVANCE='NO') dipole
    DO j = 1, npos
            WRITE(out_id, charfmat, ADVANCE='NO')                             &
            REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
            WRITE(out_id, charfmat, ADVANCE='NO')                             &
            IMAG(rho(positions(j,1),positions(j,2)))
    ENDDO                                 
    WRITE(out_id,*)

    ! Percentage complete
    pcomp = REAL(npts)/pfraco

    DO i = 0, npts-1

        ! Initialize RK variables
        k1=0.0_DP
        k2=0.0_DP
        k3=0.0_DP
        k4=0.0_DP

        ! Calculating RK variables
        comm = commute(mu,rho)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k1(n,m) = runge_k_nm(n, m, t, rho, comm)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+0.5_DP*k1)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k2(n,m) = runge_k_nm(n, m, t+0.5_DP*rk_step,              &
                                         rho+0.5_DP*k1, comm)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+0.5_DP*k2)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k3(n,m) = runge_k_nm(n, m, t+0.5_DP*rk_step,              &
                                         rho+0.5_DP*k2, comm)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+k3)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k4(n,m) = runge_k_nm(n, m, t+rk_step, rho+k3, comm)
                ENDIF
            ENDDO
        ENDDO

        rho = rho + (k1 + (2.0_DP *k2) + (2.0_DP *k3) + k4) / 6.0_DP

        ! Exploiting that p_nm=p_mn*
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m > n ) THEN
                    rho(m,n) = CONJG(rho(n,m))
                ENDIF
            ENDDO
        ENDDO

        ! Calculating dipole
        dipole = 0.0_DP
        DO n = 1,num_lev
            DO m = 1,num_lev
                dipole = dipole + REAL(rho(n,m) * mu(m,n), KIND=DP)
            ENDDO
        ENDDO

        ! Next time-step
        t   = t + rk_step

        !! PRINTING VALUES !!
        WRITE(out_id, charfmat, ADVANCE='NO') t  
        WRITE(out_id, charfmat, ADVANCE='NO') efield(field, t)
        WRITE(out_id, charfmat, ADVANCE='NO') dipole
        DO j = 1, npos
                WRITE(out_id, charfmat, ADVANCE='NO')                         &
                REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
                WRITE(out_id, charfmat, ADVANCE='NO')                         &
                IMAG(rho(positions(j,1),positions(j,2)))
        ENDDO                                 
        WRITE(out_id,*)
       
        ! Percentage complete
        pcomp = 100.0*REAL(i)/REAL(npts-1)
        IF (ABS(pcomp-pfrac) <= 1.0E-2) THEN
            WRITE(*,'(A1,A12,I3,A13)',ADVANCE='NO') char(13),                 &
                  '||------->  ', NINT(pcomp), '%  <-------||'
              pfrac = pfrac+pfraco
        ENDIF

    ENDDO

    CLOSE(out_id)

END SUBROUTINE runge

FUNCTION runge_k_nm(n, m, t, rho, comm)
        USE double
        USE fields
        USE params , ONLY : en, num_lev, gma, big_gma, field, rk_step, field, &
                            gma, big_gma, ci, rho_eq
        IMPLICIT NONE
        COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: rho, comm
        INTEGER, INTENT(IN) :: n, m
        REAL(KIND=DP), INTENT(IN) :: t
        COMPLEX(KIND=DP) :: runge_k_nm
        INTEGER :: j

        runge_k_nm = -ci * (en(n) - en(m)) * rho(n,m)                         &
                     +ci * efield(field,t) * comm(n,m)                        &
                     - gma(n,m) * rho(n,m)                                    &
                     + gma(n,m) * rho_eq(n,m)

        IF ( n == m ) THEN
            DO j = 1,num_lev
                IF ( j > n) THEN
                    runge_k_nm = runge_k_nm + big_gma(n,j)*rho(j,j)
                ELSEIF ( j < n) THEN
                    runge_k_nm = runge_k_nm - big_gma(j,n)*rho(n,n)
                ENDIF
            ENDDO
        ENDIF

        runge_k_nm = rk_step * runge_k_nm

END FUNCTION

FUNCTION commute(A, B)
    USE double
    USE omp_lib
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: A, B
    COMPLEX(KIND=DP), ALLOCATABLE                :: commute(:,:)
    COMPLEX(KIND=DP)                             :: y
    INTEGER                                      :: s, nu, n, m

    s = UBOUND(A,1)
    ALLOCATE(commute(s,s))

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n,m,nu,y) SCHEDULE(DYNAMIC)
    DO n = 1,s
        DO m = 1,s
            y = (0.0_DP, 0.0_DP)
            DO nu = 1,s
                y = y + A(n, nu)*B(nu, m) - B(n, nu)*A(nu, m)
            ENDDO
            commute(n,m) = y
        ENDDO
    ENDDO
!$OMP  END PARALLEL DO

END FUNCTION commute

END MODULE runge_mod
