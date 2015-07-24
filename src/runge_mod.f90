MODULE runge_mod
    USE double
    USE params , ONLY : rk_step, nk
    USE params , ONLY : en, gma, field, mu, num_lev, rho_eq
    USE params , ONLY : gamma_g, omega_g, coupled, &
                        theta, s_alpha, eps_eff2, nk, dist, eps_eff1, eps_0
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: runge

    REAL(KIND=DP), DIMENSION(3) :: P_sqd, P_mnp, E_sqd, E_mnp
CONTAINS

SUBROUTINE runge
    USE double
    USE params , ONLY : npts, rk_step, positions, rho_0,  &
                        out_file, npos, trange, &
                        omega, Q_sqd_start, Q_mnp_start, &
                        check_pt, proc_id, field_width, Q_sqd_end, big_gma, E0
    USE fields
    USE global_params , ONLY : std_out, power_par, au_to_ev, real_fmt
    USE print_mod
    IMPLICIT NONE

    ! rho(t)
    COMPLEX(KIND=DP), ALLOCATABLE                 :: rho(:,:)
    COMPLEX(KIND=DP) :: s(nk,3)
    ! RK variables
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: k_rho
    COMPLEX(KIND=DP) :: k_s(4,nk,3)
!    REAL(KIND=DP) :: P_mnp
    REAL(KIND=DP) :: fac(3)
    ! Dummy index variables
    INTEGER                                       :: i, j, k
    ! Time
    REAL(KIND=DP)                                 :: t
    ! Commuted matrix
    COMPLEX(KIND=DP), ALLOCATABLE                 :: comm(:,:,:)
    ! Dipole at each time-step
!    REAL(KIND=DP)                                 :: P_sqd
    REAL(KIND=DP)  :: Q_mnp
    REAL(KIND=DP)  :: Q_sqd
    REAL(KIND=DP)  :: Q
    REAL(KIND=DP)  :: dpdt_mnp(3)
    REAL(KIND=DP)  :: tmp_e(3)
!    REAL(KIND=DP)  :: E_mnp
    !
    ! Title for rho element column
    CHARACTER(LEN=2)                              :: poschar
    ! File id for output
    INTEGER                                       :: out_id

    ! Write output to file
    out_id = 30 + proc_id
    OPEN(UNIT=out_id, FILE=out_file, STATUS='REPLACE')

    ALLOCATE(                                                                 &
             rho(num_lev,num_lev),                                            &
             comm(num_lev,num_lev,3),                                           &
             k_rho(4,num_lev,num_lev)                                             &
             )
    
    ! Initializing
    rho    = rho_0
    t      = 0.0_DP
    P_sqd = 0.0_DP

    IF ( coupled ) THEN
        s = 0.0_DP
        Q_mnp = 0.0_DP
        Q_sqd = 0.0_DP
        E_sqd = efield(field, t)
        E_mnp = efield(field, t)
        P_mnp = 0.0_DP
        dpdt_mnp = 0.0_DP
    ENDIF


    ! Printing column headers
    WRITE(out_id, '(A22)', ADVANCE='NO') ' 1.t                  '
    WRITE(out_id, '(A22)', ADVANCE='NO') ' 2.field_x            '
    WRITE(out_id, '(A22)', ADVANCE='NO') ' 3.field_y            '
    WRITE(out_id, '(A22)', ADVANCE='NO') ' 4.field_z            '
    IF ( coupled ) THEN
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 5.E_SQD_x              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 6.E_SQD_y              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 7.E_SQD_z              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 8.E_MNP_x              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 9.E_MNP_y              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 10.E_MNP_z             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 11.P_SQD_x             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 12.P_SQD_z             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 13.P_SQD_y             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 14.P_MNP_x             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 15.P_MNP_y             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 16.P_MNP_z             '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 17.d/dt(P_MNP)_x       '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 18.d/dt(P_MNP)_y       '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 19.d/dt(P_MNP)_z       '
    ELSE
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 3.P_SQD_x              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 4.P_SQD_y              '
        WRITE(out_id, '(A22)', ADVANCE='NO') ' 5.P_SQD_x              '
    ENDIF
    k=19
    DO j = 1, npos 
        k=k+1
        WRITE(poschar, '(I2)') k
        WRITE(out_id, '(A8)', ADVANCE='NO') ' '//poschar//'.rho_'
        WRITE(poschar, '(I2)') positions(j,1)
        WRITE(out_id, '(A3)', ADVANCE='NO') poschar//','
        WRITE(poschar, '(I2)') positions(j,2)
        WRITE(out_id, '(A2)', ADVANCE='NO') poschar
        IF ( positions(j,1) == positions(j,2) ) THEN
            WRITE(out_id, '(A9)', ADVANCE='NO') ' '
        ELSE
            k=k+1
            WRITE(out_id, '(A31)', ADVANCE='NO') ' '
        ENDIF
    ENDDO
    WRITE(out_id,*)

    ! Parameters in RK algorithm
    fac(1) = 0.5_DP
    fac(2) = 0.5_DP
    fac(3) = 1.0_DP

    DO i = 0, npts

        IF ( MOD(i, check_pt) == 0 ) THEN
            !! PRINTING VALUES !!
            WRITE(out_id, real_fmt, ADVANCE='NO') t  
            IF ( coupled ) THEN
                tmp_e = efield(field,t)
                WRITE(out_id, real_fmt, ADVANCE='NO') tmp_e(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') tmp_e(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') tmp_e(3)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_sqd(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_sqd(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_sqd(3)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_mnp(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_mnp(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') E_mnp(3)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_sqd(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_sqd(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_sqd(3)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_mnp(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_mnp(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_mnp(3)
                WRITE(out_id, real_fmt, ADVANCE='NO') dpdt_mnp(1)
                WRITE(out_id, real_fmt, ADVANCE='NO') dpdt_mnp(2)
                WRITE(out_id, real_fmt, ADVANCE='NO') dpdt_mnp(3)
            ELSE
                WRITE(out_id, real_fmt, ADVANCE='NO') efield(field, t)
                WRITE(out_id, real_fmt, ADVANCE='NO') P_sqd
            ENDIF
            DO j = 1, npos
                    WRITE(out_id, real_fmt, ADVANCE='NO')                         &
                    REAL(rho(positions(j,1),positions(j,2)),KIND=DP)
                    IF (positions(j,1) /= positions(j,2)) THEN
                        WRITE(out_id, real_fmt, ADVANCE='NO')                         &
                        AIMAG(rho(positions(j,1),positions(j,2)))
                    ENDIF
            ENDDO                                 
            WRITE(out_id,*)
        ENDIF
      
        IF ( i == npts ) THEN
           EXIT
        ENDIF

        CALL rk_de(t,                rho,                     k_rho(1,:,:), s,                 k_s(1,:,:))
        CALL rk_de(t+0.5_DP*rk_step, rho+0.5_DP*k_rho(1,:,:), k_rho(2,:,:), s+0.5_DP*k_s(1,:,:), k_s(2,:,:))
        CALL rk_de(t+0.5_DP*rk_step, rho+0.5_DP*k_rho(2,:,:), k_rho(3,:,:), s+0.5_DP*k_s(2,:,:), k_s(3,:,:))
        CALL rk_de(t+rk_step,        rho+k_rho(3,:,:),        k_rho(4,:,:), s+k_s(3,:,:),        k_s(4,:,:))

        rho = rho + (k_rho(1,:,:) + 2.0_DP*k_rho(2,:,:) + 2.0_DP*k_rho(3,:,:) + k_rho(4,:,:))/6.0_DP
        t = t + rk_step
        P_sqd = REAL(trace(matmul_vec(rho,mu)))

        IF ( coupled ) THEN
            s = s + (k_s(1,:,:) + 2.0_DP*k_s(2,:,:) + 2.0_DP*k_s(3,:,:) + k_s(4,:,:))/6.0_DP

            ! Calculating gold dipole
            P_mnp = 0.0_DP
            dpdt_mnp = 0.0_DP
            DO j=1,nk
                P_mnp = P_mnp + theta(j,:)*2.0_DP*REAL(s(j,:))
                dpdt_mnp = dpdt_mnp + 2.0_DP*theta(j,:)*( omega_g(j,:)*AIMAG(s(j,:)) - &
                                                  gamma_g(j,:)*REAL(s(j,:)) )
            ENDDO

            E_mnp = efield(field,t) + s_alpha*P_sqd/eps_0/dist**3
            E_sqd = efield(field,t) + s_alpha*P_mnp/eps_0/dist**3

            IF ( t >= Q_mnp_start ) THEN
                ! integrating E_mnp*dp/dt
                Q_mnp = Q_mnp + dot_prod_real_real(dpdt_mnp,E_mnp)
            ENDIF

        ENDIF

        IF ( t >= Q_sqd_start .AND. t <= Q_sqd_end) THEN
            ! integrating rho(2,2)
            Q_sqd = Q_sqd + REAL(rho(2,2))
        ENDIF

    ENDDO

    CALL print_param_change


    IF ( coupled ) THEN

!        Q_sqd = power_par*(en(2)-en(1))*gma(1,1)*Q_sqd*rk_step/(Q_sqd_end-Q_sqd_start)
        Q_sqd = power_par*(en(2)-en(1))*big_gma(2,1)*Q_sqd*rk_step/(Q_sqd_end-Q_sqd_start)
        Q_mnp = power_par*Q_mnp*rk_step/(trange-Q_mnp_start)
        Q = Q_mnp + Q_sqd

        WRITE(std_out,real_fmt, ADVANCE='NO') Q_sqd
        WRITE(std_out,real_fmt, ADVANCE='NO') Q_mnp
        WRITE(std_out,real_fmt) Q

    ELSE

        Q_sqd = Q_sqd*rk_step/(Q_sqd_end-Q_sqd_start)

        WRITE(std_out,real_fmt, ADVANCE='NO') Q_sqd
        WRITE(std_out, *)

    ENDIF

    CLOSE(out_id)

END SUBROUTINE runge

SUBROUTINE print_param_change
    USE params , ONLY : field_change_param, field_width, field_height, & 
                        omega, pulse_area, I0
    USE global_params , ONLY : std_out, real_fmt, au_to_ev

    IMPLICIT NONE
    
!    pulse_area = REAL(mu(1,2), KIND=DP)*1.064467019431226_DP*field_height*field_width

    SELECTCASE ( field_change_param )

        CASE( 'field_height' )
            WRITE(std_out,real_fmt, ADVANCE='NO') field_height

        CASE( 'pulse_area' )
            WRITE(std_out,real_fmt, ADVANCE='NO') pulse_area

        CASE ( 'field_width' )
            WRITE(std_out,real_fmt, ADVANCE='NO') field_width

        CASE ( 'omega' )
            WRITE(std_out,real_fmt, ADVANCE='NO') omega*au_to_ev

        CASE( 'I0' )
            WRITE(std_out,real_fmt, ADVANCE='NO') I0


    END SELECT

END SUBROUTINE print_param_change

SUBROUTINE rk_de(t_in, rho_in, rho_out, s_in, s_out)
    USE double
    USE global_params , ONLY : ci
    USE params, ONLY : big_gma, omega, E0
    USE fields !, ONLY : efield
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: rho_in(:,:)
    COMPLEX(KIND=DP), INTENT(IN) :: s_in(:,:)
    REAL(KIND=DP), INTENT(IN) :: t_in
    COMPLEX(KIND=DP), INTENT(OUT) :: rho_out(num_lev,num_lev)
    COMPLEX(KIND=DP) :: comm(num_lev,num_lev,3)
    COMPLEX(KIND=DP), INTENT(OUT) :: s_out(nk,3)
    REAL(KIND=DP) :: ext_field(3)
    REAL(KIND=DP) :: Gamma_n, Gamma_m
    COMPLEX(KIND=DP) :: eta

    INTEGER :: n, m, nu

    IF ( coupled ) THEN

        P_sqd = REAL(trace(matmul_vec(rho_in,mu)))
        P_mnp = 0.0_DP
        DO n=1,nk
            P_mnp = P_mnp + theta(n,:)*2.0_DP*REAL(s_in(n,:))
        ENDDO

        E_mnp = efield(field,t_in) + s_alpha*P_sqd/eps_0/dist**3
        E_sqd = efield(field,t_in) + s_alpha*P_mnp/eps_0/dist**3

        DO n=1,nk
            s_out(n,:) = -(gamma_g(n,:) + ci*omega_g(n,:))*s_in(n,:) + ci*E_mnp
        ENDDO

        s_out = s_out*rk_step

        ext_field = E_sqd
    ELSE
        ext_field = efield(field, t_in)
    ENDIF

    comm = commute(mu, rho_in)
    DO n=1,num_lev
        DO m=1,num_lev
!           rho_out(n,m) = -ci*(en(n)-en(m))*rho_in(n,m)                 &
!                          - gma(n,m)*(rho_in(n,m) - rho_eq(n,m))        &
!                          + ci*ext_field*comm(n,m)

           IF ( n == m) THEN
               eta = 0.0_DP
               DO nu = 1, num_lev
                   IF ( nu > n ) THEN
                       eta = eta + big_gma(n,nu)*rho_in(nu,nu)
                   ELSEIF ( nu < n ) THEN
                       eta = eta - big_gma(nu,n)*rho_in(n,n)
                   ENDIF
               ENDDO
           ELSE
               Gamma_n = 0.0_DP
               Gamma_m = 0.0_DP
               DO nu = 1, num_lev
                   IF ( nu < n ) THEN
                       Gamma_n = Gamma_n + big_gma(nu,n)
                   ENDIF
                   IF ( nu < m) THEN
                       Gamma_m = Gamma_m + big_gma(nu,m)
                   ENDIF
               ENDDO
               eta = -0.5_DP*(Gamma_n + Gamma_m)*rho_in(n,m)
           ENDIF

           rho_out(n,m) = -ci*(en(n)-en(m))*rho_in(n,m)   &
               +ci*dot_prod_real_complex(ext_field,comm(n,m,:))         &
                          + eta
        ENDDO
    ENDDO

    rho_out = rho_out*rk_step

END SUBROUTINE rk_de

FUNCTION matmul_vec(A,B)
    USE double
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
    REAL(KIND=DP), INTENT(IN) :: B(:,:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: matmul_vec(:,:,:)
    COMPLEX(KIND=DP) :: y

    INTEGER :: s, n, m, nu, l

    s = UBOUND(A,1)
    ALLOCATE(matmul_vec(s,s,3))

    DO n = 1,s
        DO m = 1,s
            DO l = 1,3
                y = (0.0_DP, 0.0_DP)
                DO nu = 1,s
                    y = y + A(n, nu)*B(nu, m, l)
                ENDDO
                matmul_vec(n,m,l) = y
            ENDDO
        ENDDO
    ENDDO

END FUNCTION matmul_vec

FUNCTION trace(A)
    USE double
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:,:)
    COMPLEX(KIND=DP) :: trace(3)
    INTEGER :: i, l

    DO l = 1,3
        trace(l) = 0.0_DP
        DO i = 1,SIZE(A,1)
            trace(l) = trace(l) + A(i,i,l)
        ENDDO
    ENDDO

END FUNCTION trace

FUNCTION commute(A, B)
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: A(:,:,:)
    COMPLEX(KIND=DP), INTENT(IN) :: B(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE                :: commute(:,:,:)
    COMPLEX(KIND=DP)                             :: y
    INTEGER                                      :: s, nu, n, m, l

    s = UBOUND(A,1)
    ALLOCATE(commute(s,s,3))

    DO n = 1,s
        DO m = 1,s
            DO l = 1,3
                y = (0.0_DP, 0.0_DP)
                DO nu = 1,s
                    y = y + A(n, nu, l)*B(nu, m) - B(n, nu)*A(nu, m, l)
                ENDDO
                commute(n,m,l) = y
            ENDDO
        ENDDO
    ENDDO

END FUNCTION commute

FUNCTION dot_prod_real_real(A,B)
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: A(3)
    REAL(KIND=DP), INTENT(IN) :: B(3)
    REAL(KIND=DP) :: dot_prod_real_real
    INTEGER :: i

    dot_prod_real_real = 0.0_DP
    DO i = 1,3
        dot_prod_real_real = dot_prod_real_real + A(i)*B(i)
    ENDDO

END FUNCTION dot_prod_real_real

FUNCTION dot_prod_real_complex(A,B)
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: A(3)
    COMPLEX(KIND=DP), INTENT(IN) :: B(3)
    COMPLEX(KIND=DP) :: dot_prod_real_complex
    INTEGER :: i

    dot_prod_real_complex = (0.0_DP,0.0_DP)
    DO i = 1,3
        dot_prod_real_complex = dot_prod_real_complex + A(i)*B(i)
    ENDDO

END FUNCTION dot_prod_real_complex

END MODULE runge_mod
