MODULE runge_mod
    USE double
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: runge

CONTAINS

SUBROUTINE runge
    USE double
    USE params , ONLY : num_lev, field, npts, rk_step, positions, mu, rho_0,  &
                        out_file, npos, E0, omega_au, gma, trange_au, omega_ev,&
                        rad, eps_eff1
    USE fields
    USE global_params , ONLY : pi, ci, power_par, std_out
    USE params , ONLY : omega_g, gamma_g, theta, en, s_alpha, eps_eff2
    USE print_mat_mod
    IMPLICIT NONE

    ! rho(t)
    COMPLEX(KIND=DP), ALLOCATABLE                 :: rho(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_tmp(:,:)
    COMPLEX(KIND=DP) :: s
    COMPLEX(KIND=DP) :: s_tmp
    ! RK variables
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: k_rho
    COMPLEX(KIND=DP) :: k_s(4)
    REAL(KIND=DP) :: Pe
    REAL(KIND=DP) :: fac(3)
    ! Dummy index variables
    INTEGER                                       :: i, j, n, m 
    ! Time
    REAL(KIND=DP)                                 :: t
    REAL(KIND=DP)                                 :: t_tmp
    REAL(KIND=DP) :: field_t
    ! Commuted matrix
    COMPLEX(KIND=DP), ALLOCATABLE                 :: comm(:,:)
    ! Dipole at each time-step
    REAL(KIND=DP)                                 :: dipole
    REAL(KIND=DP)  :: Q_mnp
    REAL(KIND=DP)  :: abs_rate_0
    REAL(KIND=DP)  :: abs_rate_approx
    REAL(KIND=DP)  :: dpdt
    REAL(KIND=DP)  :: Pe_old
    REAL(KIND=DP)  :: E_mnp
    !
    ! Title format specifier
    CHARACTER(LEN=64)                             :: charfmat = '(ES22.14)'
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
             rho_tmp(num_lev,num_lev),                                        &
             comm(num_lev,num_lev),                                           &
             k_rho(4,num_lev,num_lev)                                             &
             )
    
    ! Initializing
    rho    = rho_0
    s = 0.0_DP
    t      = 0.0_DP
    dipole = 0.0_DP
    Pe = 0.0_DP

    ! Printing column headers
    WRITE(out_id, '(A22)', ADVANCE='NO') ' t              '
    fieldchar=' '//TRIM(field)//'(t)'
    WRITE(out_id, '(A22)', ADVANCE='NO') fieldchar
!    WRITE(out_id, '(A22)', ADVANCE='NO') 'effective field '
    WRITE(out_id, '(A22)', ADVANCE='NO') ' P_SQD          '
    WRITE(out_id, '(A22)', ADVANCE='NO') ' P_MNP          '
!    WRITE(out_id, '(A22)', ADVANCE='NO') ' d/dt P_MNP     '
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
    WRITE(out_id, charfmat, ADVANCE='NO') Pe
!    WRITE(out_id, charfmat, ADVANCE='NO') dpdt
    DO j = 1, npos
            WRITE(out_id, charfmat, ADVANCE='NO')                             &
            REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
            WRITE(out_id, charfmat, ADVANCE='NO')                             &
            AIMAG(rho(positions(j,1),positions(j,2)))
    ENDDO                                 

    WRITE(out_id,*)

!    ! Percentage complete
!    pcomp = REAL(npts)/pfraco

    ! Parameters in RK algorithm
    fac(1) = 0.5_DP
    fac(2) = 0.5_DP
    fac(3) = 1.0_DP

!    abs_rate = 0.0_DP
    Q_mnp = 0.0_DP

!theta = ci*2.0_DP*gamma_g*theta*(REAL(mu(1,2)*mu(1,2)))/REAL(gma(1,2))
!theta = ci*2.0_DP*gamma_g*theta*3.40824295834338e+04

    DO i = 0, npts-1

        CALL rk_de(rho, s, t, k_rho(1,:,:), k_s(1))
        k_rho(1,:,:) = rk_step*k_rho(1,:,:)
        k_s(1) = rk_step*k_s(1)

        CALL rk_de(rho+0.5_DP*k_rho(1,:,:), s+0.5_DP*k_s(1), t+0.5_DP*rk_step, k_rho(2,:,:), k_s(2))
        k_rho(2,:,:) = rk_step*k_rho(2,:,:)
        k_s(2) = rk_step*k_s(2)

        CALL rk_de(rho+0.5_DP*k_rho(2,:,:), s+0.5_DP*k_s(2), t+0.5_DP*rk_step, k_rho(3,:,:), k_s(3))
        k_rho(3,:,:) = rk_step*k_rho(3,:,:)
        k_s(3) = rk_step*k_s(3)

        CALL rk_de(rho+k_rho(3,:,:), s+k_s(3), t+rk_step, k_rho(4,:,:), k_s(4))
        k_rho(4,:,:) = rk_step*k_rho(4,:,:)
        k_s(4) = rk_step*k_s(4)

        rho = rho + (k_rho(1,:,:) + 2.0_DP*k_rho(2,:,:) + 2.0_DP*k_rho(3,:,:) + k_rho(4,:,:))/6.0_DP
        s = s + (k_s(1) + 2.0_DP*k_s(2) + 2.0_DP*k_s(3) + k_s(4))/6.0_DP
        t = t + rk_step

        ! Calculating dipole
!        dipole = 0.0_DP
!        DO n = 1,num_lev
!            DO m = 1,num_lev
!                dipole = dipole + REAL(rho(n,m) * mu(m,n))
!            ENDDO
!        ENDDO

        dipole = trace(MATMUL(rho,mu))

        Pe_old = Pe
        ! Calculating gold dipole
        Pe = rad**3*theta*2.0_DP*REAL(s)
!        Pe = theta*2.0_DP*REAL(s)

        dpdt = (Pe - Pe_old)/rk_step
        E_mnp = efield(field,t) + s_alpha*dipole/eps_eff2/rad**3

        Q_mnp = Q_mnp + dpdt*E_mnp

!        abs_rate = abs_rate + SIN(omega_au*t)*(Pe + dipole)
!        abs_rate = abs_rate + SIN(omega_au*t)*(Pe)

        !! PRINTING VALUES !!
        WRITE(out_id, charfmat, ADVANCE='NO') t  
        WRITE(out_id, charfmat, ADVANCE='NO') efield(field, t)
        WRITE(out_id, charfmat, ADVANCE='NO') dipole
        WRITE(out_id, charfmat, ADVANCE='NO') Pe
!        WRITE(out_id, charfmat, ADVANCE='NO') dpdt
        DO j = 1, npos
                WRITE(out_id, charfmat, ADVANCE='NO')                         &
                REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
                WRITE(out_id, charfmat, ADVANCE='NO')                         &
                AIMAG(rho(positions(j,1),positions(j,2)))
        ENDDO                                 
        WRITE(out_id,*)
       
        ! Percentage complete
!        pcomp = 100.0*REAL(i)/REAL(npts-1)
!        IF (ABS(pcomp-pfrac) <= 1.0E-2) THEN
!            WRITE(*,'(A1,A12,I3,A13)',ADVANCE='NO') char(13),                 &
!                  '||------->  ', NINT(pcomp), '%  <-------||'
!              pfrac = pfrac+pfraco
!        ENDIF

    ENDDO

!    abs_rate = 4.0_DP/3.0_DP*pi*rad**3*E0/trange_au*(Pe*COS(omega_au*trange_au) + omega_au*rk_step*abs_rate)
!    abs_rate_0 = 1.0_DP/3.0_DP*E0**2 * rad**6 *omega_au * theta *gamma_g/((omega_g-omega_au)**2 + gamma_g**2)

    Q_mnp = 4.0_DP/3.0_DP*pi*rad**3*Q_mnp*rk_step/trange_au
    Q_mnp = Q_mnp*0.5_DP/pi

    WRITE(std_out,charfmat, ADVANCE='NO') omega_ev
    WRITE(std_out,charfmat) Q_mnp
!    WRITE(0,charfmat, ADVANCE='NO') abs_rate*0.5_DP/pi
!    WRITE(0,charfmat, ADVANCE='NO') abs_rate_0
!    WRITE(0,charfmat, ADVANCE='NO') abs_rate_approx*0.5_DP/pi
!    WRITE(0,charfmat) 4.0_DP/3.0_DP*pi*151.0_DP*REAL(omega_au*mu(1,2)**2*E0**2*gma(1,2)/(2.0_DP*eps_eff1*((en(2)-omega_au)**2+gma(1,2)**2)))
!    WRITE(0,charfmat) REAL(omega_au*mu(1,2)**2*E0**2*gma(1,2)/(2.0_DP*eps_eff1*((en(2)-omega_au)**2+gma(1,2)**2)))

    CLOSE(out_id)

END SUBROUTINE runge

SUBROUTINE rk_de(rho_in, s_in, t_in, rho_out, s_out)
    USE double
    USE params , ONLY : s_alpha, eps_eff1, eps_eff2, dist, theta, omega_g, &
                            gamma_g, rad
    USE params , ONLY : en, gma, field, mu, num_lev, rho_eq
    USE global_params , ONLY : ci
    USE fields , ONLY : efield
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: rho_in(:,:)
    COMPLEX(KIND=DP), INTENT(IN) :: s_in
    REAL(KIND=DP), INTENT(IN) :: t_in
    COMPLEX(KIND=DP), INTENT(OUT) :: rho_out(num_lev,num_lev)
    COMPLEX(KIND=DP) :: comm(num_lev,num_lev)
    COMPLEX(KIND=DP), INTENT(OUT) :: s_out
    REAL(KIND=DP) :: P_mnp, E_sqd, E_mnp

    INTEGER :: n, m

    P_mnp = rad**3*theta*2.0_DP*REAL(s_in)
!    P_mnp = theta*2.0_DP*REAL(s_in)
    E_sqd = efield(field,t_in) + s_alpha*P_mnp/(eps_eff1*dist**3)
    E_mnp = efield(field,t_in) + s_alpha*trace(MATMUL(mu,rho_in))/(eps_eff2*dist**3)

    comm = commute(mu, rho_in)
    DO n=1,num_lev
        DO m=1,num_lev
           rho_out(n,m) = -ci*(en(n)-en(m))*rho_in(n,m)                 &
                          - gma(n,m)*(rho_in(n,m) - rho_eq(n,m))        &
                          + ci*E_sqd*comm(n,m)
        ENDDO
    ENDDO

    s_out = -(gamma_g + ci*omega_g)*s_in + ci*E_mnp

END SUBROUTINE rk_de

FUNCTION runge_k_nm(n, m, field_t, rho, comm)
        USE double
        USE fields
        USE params , ONLY : en, num_lev, gma, big_gma, field, rk_step, field, &
                            gma, big_gma, rho_eq
        USE global_params , ONLY : ci
        IMPLICIT NONE
        COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: rho, comm
        INTEGER, INTENT(IN) :: n, m
        REAL(KIND=DP), INTENT(IN) :: field_t
        COMPLEX(KIND=DP) :: runge_k_nm
        INTEGER :: j

        runge_k_nm = -ci * (en(n) - en(m)) * rho(n,m)                         &
                     +ci * field_t * comm(n,m)                        &
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

FUNCTION trace(A)
    USE double
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
    COMPLEX(KIND=DP) :: trace
    INTEGER :: i

    trace = 0.0_DP
    DO i = 1,SIZE(A,1)
        trace = trace + A(i,i)
    ENDDO

END FUNCTION trace

FUNCTION commute(A, B)
    USE double
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN), DIMENSION(:,:) :: A, B
    COMPLEX(KIND=DP), ALLOCATABLE                :: commute(:,:)
    COMPLEX(KIND=DP)                             :: y
    INTEGER                                      :: s, nu, n, m

    s = UBOUND(A,1)
    ALLOCATE(commute(s,s))

    DO n = 1,s
        DO m = 1,s
            y = (0.0_DP, 0.0_DP)
            DO nu = 1,s
                y = y + A(n, nu)*B(nu, m) - B(n, nu)*A(nu, m)
            ENDDO
            commute(n,m) = y
        ENDDO
    ENDDO

END FUNCTION commute

END MODULE runge_mod
