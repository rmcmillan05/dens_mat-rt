MODULE runge_mod
    USE double
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: runge

CONTAINS

SUBROUTINE runge
    USE double
    USE params
    USE fields
!    USE params, ONLY : time_par 
    IMPLICIT NONE
    INTEGER :: out_id=50
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: rho
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:)              :: k1, k2, k3, k4
!    COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:)               :: rho_0, mu
!    REAL(KIND=DP), INTENT(IN), DIMENSION(:,:) :: gma, big_gma
    INTEGER :: npts ! total npts
!    INTEGER, INTENT(IN), DIMENSION(:,:)                        :: positions
    INTEGER                                                    :: i, j
    INTEGER :: n,m
    REAL(KIND=DP)                                              :: h, t
    CHARACTER(LEN=2)        :: poschar
    CHARACTER(LEN=16)       :: fieldchar
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:) :: comm, rho_eq
    CHARACTER(LEN=64) :: charfmat = '(ES16.8)'
    REAL(KIND=DP)   :: dipole
    REAL, PARAMETER :: pfraco=5.0 ! Percentage increment in screen update
    REAL :: pfrac=0.0 ! Percentage increment in screen update
    REAL :: pcomp

    h = 1.0_DP/REAL(nptspau,KIND=DP)
    npts = NINT(nptspau * trange_au)

    OPEN(UNIT=out_id, FILE=out_file, STATUS='REPLACE')

    ALLOCATE(                         &
             rho(num_lev,num_lev),    &
             rho_eq(num_lev,num_lev), &
             comm(num_lev,num_lev),   &
             k1(num_lev,num_lev),     &
             k2(num_lev,num_lev),     &
             k3(num_lev,num_lev),     &
             k4(num_lev,num_lev)      &
             )
    
    rho_eq = 0.0_DP !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rho = rho_0
    t   = 0.0_DP
    dipole = 0.0_DP

    !! PRINTING TITLES FOR COLUMNS !!
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

    !! PRINTING VALUES AT t=0 !!
    WRITE(out_id, charfmat, ADVANCE='NO') t
    WRITE(out_id, charfmat, ADVANCE='NO') efield(field, t)
    WRITE(out_id, charfmat, ADVANCE='NO') dipole
    DO j = 1, npos
            WRITE(out_id, charfmat, ADVANCE='NO') &
            REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
            WRITE(out_id, charfmat, ADVANCE='NO') &
            IMAG(rho(positions(j,1),positions(j,2)))
    ENDDO                                 
    WRITE(out_id,*)

    pcomp = REAL(npts)/pfraco

    DO i = 0, npts-1

        k1=0.0_DP
        k2=0.0_DP
        k3=0.0_DP
        k4=0.0_DP

        comm = commute(mu,rho)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k1(n,m) = runge_k_nm(n,m,t,rho,  h,field,gma,big_gma,comm,rho_eq)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+0.5_DP*k1)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k2(n,m) = runge_k_nm(n,m,t+0.5_DP*h,rho+0.5_DP*k1,  h,field,gma,big_gma,comm,rho_eq)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+0.5_DP*k2)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k3(n,m) = runge_k_nm(n,m,t+0.5_DP*h,rho+0.5_DP*k2,  h,field,gma,big_gma,comm,rho_eq)
                ENDIF
            ENDDO
        ENDDO

        comm = commute(mu,rho+k3)
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m >= n ) THEN
                    k4(n,m) = runge_k_nm(n,m,t+h,rho+k3,  h,field,gma,big_gma,comm,rho_eq)
                ENDIF
            ENDDO
        ENDDO

        rho = rho + (k1 + (2.0_DP *k2) + (2.0_DP *k3) + k4) / 6.0_DP
        DO n = 1,num_lev
            DO m =1,num_lev
                IF ( m > n ) THEN
                    rho(m,n) = CONJG(rho(n,m))
                ENDIF
            ENDDO
        ENDDO
        t   = t + h

        dipole = 0.0_DP
        DO n = 1,num_lev
            DO m = 1,num_lev
                dipole = dipole + REAL(rho(n,m) * mu(m,n), KIND=DP)
            ENDDO
        ENDDO

        !! PRINTING VALUES !!
        WRITE(out_id, charfmat, ADVANCE='NO') t  ! t in a.u.
        WRITE(out_id, charfmat, ADVANCE='NO') efield(field, t)
        WRITE(out_id, charfmat, ADVANCE='NO') dipole
        DO j = 1, npos
                WRITE(out_id, charfmat, ADVANCE='NO') &
                REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
                WRITE(out_id, charfmat, ADVANCE='NO') &
                IMAG(rho(positions(j,1),positions(j,2)))
        ENDDO                                 
        WRITE(out_id,*)
       
        pcomp = 100.0*REAL(i)/REAL(npts-1)
        IF (ABS(pcomp-pfrac) <= 1.0E-2) THEN
            WRITE(*,'(A1,A12,I3,A13)',ADVANCE='NO') char(13), &
                  '||------->  ', NINT(pcomp), '%  <-------||'
              pfrac = pfrac+pfraco
        ENDIF


    ENDDO

    CLOSE(out_id)

END SUBROUTINE runge

FUNCTION commute(A, B)
    USE double
    USE omp_lib
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: A, B
    COMPLEX(KIND=DP), ALLOCATABLE                :: commute(:,:)
    COMPLEX(KIND=DP) :: y
    INTEGER                                   :: s, nu, n, m

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

FUNCTION runge_k_nm(n,m,t,rho,  h,field,gma,big_gma,comm,rho_eq)
        USE double
        USE fields
        USE params , ONLY : en, num_lev
        IMPLICIT NONE
        COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: rho, comm, rho_eq
        REAL(KIND=DP), INTENT(IN), DIMENSION(:,:) :: gma, big_gma
        INTEGER, INTENT(IN) :: n, m
        REAL(KIND=DP), INTENT(IN) :: t, h
        CHARACTER(LEN=*), INTENT(IN) :: field
        COMPLEX(KIND=DP), PARAMETER :: ci=(0.0_DP,1.0_DP)
        COMPLEX(KIND=DP) :: runge_k_nm
        INTEGER :: j

        runge_k_nm = -ci * (en(n) - en(m)) * rho(n,m)   &
                     +ci * efield(field,t) * comm(n,m)  &
                     - gma(n,m) * rho(n,m)    &
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

        runge_k_nm = h * runge_k_nm

END FUNCTION

END MODULE runge_mod
