MODULE print_params_mod
    USE double
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: print_params

CONTAINS

SUBROUTINE print_params
    USE double
    USE params , ONLY : out_folder, omega_ev, omega_au, lambda, I0, E0,       &
                        jname, in_folder, timestamp, in_file, en, npts,       &
                        params_file, trange_au, field, nptspau, rk_step,      &
                        step_centre, step_width, step_height,                 &
                        pulse_start, pulse_cycles, pulse_phase
    USE print_mat_mod
    USE params , ONLY : dist, s_alpha, theta, eps_eff1, omega_g, gamma_g, rad, &
                        eps_eff2, nk
!    USE global_params , ONLY : ci
    IMPLICIT NONE
!    COMPLEX(KIND=DP) :: G
    INTEGER :: i
    CHARACTER(LEN=3) :: istr
    
    ! Current working directory
    CHARACTER(LEN=256) :: cwd
    ! File handle
    INTEGER            :: fid = 10

    ! Get current working directory
    CALL GETCWD(cwd)

    ! Writing parameters to file
    OPEN(fid, FILE=params_file, STATUS='REPLACE')

        CALL print_break(fid)
        CALL print_str(ADJUSTR(TRIM(timestamp)), fid)
        CALL print_break(fid)
        CALL print_break(fid)
        CALL print_str('Laser Properties', fid)
        CALL print_break(fid)
        WRITE(fid,*)
        CALL print_str_str('Field type', field, fid)
        IF (field == 'step') THEN
            CALL print_str_num_real('> Step height', step_height, fid)
            CALL print_str_num_real('> Step width', step_width, fid)
            CALL print_str_num_real('> Step centre', step_centre, fid)
        ELSEIF (field == 'pulse') THEN
            CALL print_str_num_real('> Pulse start', pulse_start, fid)
            CALL print_str_num_real('> Pulse phase', pulse_phase, fid)
            CALL print_str_num_real('> Number of cycles', pulse_cycles, fid)
        ENDIF
        IF (field /= 'step' .AND. field /= 'zero') THEN
            CALL print_str_num_real('Laser Wavelength (nm)', lambda, fid)
            CALL print_str_num_real('Laser Frequency (a.u.)', omega_au, fid)
            CALL print_str_num_real('Laser Energy (eV)', omega_ev, fid)
            CALL print_str_num_real('Laser Intenstiy (W/cm^2)', I0, fid)
            CALL print_str_num_real('Laser Amplitude (E0) (a.u.)', E0, fid)
        ENDIF

        WRITE(fid,*)
        DO i = 1,nk
            WRITE(istr,'(I3)') i
            CALL print_str_num_real('theta_g_'//ADJUSTL(istr)//' (a.u.)', theta(i), fid)
            CALL print_str_num_real('omega_g_'//ADJUSTL(istr)//' (a.u.)', omega_g(i), fid)
            CALL print_str_num_real('gamma_g_'//ADJUSTL(istr)//' (a.u.)', gamma_g(i), fid)
        ENDDO
        CALL print_str_num_real('Separation (a.u.)', dist, fid)
        CALL print_str_num_real('MNP Diameter (a.u.)', rad, fid)
        CALL print_str_num_real('s_alpha (a.u.)', s_alpha, fid)
        CALL print_str_num_real('eps_eff1 (a.u.)', eps_eff1, fid)
        CALL print_str_num_real('eps_eff2 (a.u.)', eps_eff2, fid)
!        G = theta*(omega_au-omega_g+ci*gamma_g)/((omega_au-omega_g)**2+gamma_g**2)
!        G = s_alpha**2*G*rad**3*mu(1,2)**2/(eps_eff1*eps_eff2*dist**6)
!        CALL print_str_num_complex('G', G, fid)
!        CALL print_str_num_real('strength of eff. field', s_alpha*theta &
!                                                         /(eps_eff1*dist**3),  &
!                                                        fid)

        WRITE(fid, *)
        CALL print_break(fid)
        CALL print_str('Runge-Kutta Parameters', fid)
        CALL print_break(fid)
        WRITE(fid, *)
        CALL print_str_num_real('Number of Points per a.u.', nptspau, fid)
        CALL print_str_num_real('Total number of points used',                 &
                                 REAL(npts, KIND=DP), fid)
        CALL print_str_num_real('RK time-step (a.u.)', rk_step, fid)
        CALL print_str_num_real('Total Propagation Time (a.u)', trange_au, fid)
        WRITE(fid, *)
        CALL print_break(fid)
        CALL print_str('Energy Levels', fid)
        CALL print_break(fid)
        WRITE(fid, *)
        CALL print_vec_real(en, fid)
        WRITE(fid, *)
        CALL print_break(fid)
        CALL print_str('Directories', fid)
        CALL print_break(fid)
        WRITE(fid, *)
        CALL print_str('Input folder', fid)
        CALL print_str(TRIM(cwd)//'/'//TRIM(in_folder), fid)
        WRITE(fid, *)
        CALL print_str('Output folder', fid)
        CALL print_str(TRIM(cwd)//'/'//TRIM(out_folder), fid)
        WRITE(fid, *)
        CALL print_str('Output file prefix', fid)
        CALL print_str(jname, fid)

        WRITE(fid, *)
        CALL print_break(fid)
        CALL print_str('Input file used', fid)
        CALL print_break(fid)
        WRITE(fid, *)
        CALL print_str(in_file, fid)
        WRITE(fid, *)

    CLOSE(10)

    CALL SYSTEM('cat '//TRIM(in_file)//' >> '//params_file)
    OPEN(UNIT=10, FILE=params_file, STATUS='OLD', POSITION='APPEND')

        WRITE(fid, *)
        CALL print_break(fid)
        WRITE(fid,'(A80)') '                                    - EOF -    &
                           &                                 '
        CALL print_break(fid)

    CLOSE(10)

END SUBROUTINE print_params

END MODULE print_params_mod
