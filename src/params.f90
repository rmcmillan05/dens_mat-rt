MODULE params
    USE double
    IMPLICIT NONE
    
    ! Input field
    CHARACTER(LEN=256) :: field
    CHARACTER(LEN=256) :: field_change_param
    REAL(KIND=DP)      :: field_param_from
    REAL(KIND=DP)      :: field_param_to
    INTEGER            :: field_param_npts
    REAL(KIND=DP)      :: field_param_step
    REAL(KIND=DP)      :: omega, omega_control
    REAL(KIND=DP)      :: I0, I0_control
    REAL(KIND=DP)      :: E0, E0_control
    !
    ! SQD_MNP
    !
    LOGICAL       :: coupled
    REAL(KIND=DP) :: s_alpha
    REAL(KIND=DP) :: eps_0
    REAL(KIND=DP) :: eps_s
    REAL(KIND=DP) :: ratio
    REAL(KIND=DP) :: dist
    REAL(KIND=DP) :: surf_dist
    REAL(KIND=DP) :: rad
    REAL(KIND=DP) :: s_rad
    REAL(KIND=DP) :: eps_eff1
    REAL(KIND=DP) :: eps_eff2

    CHARACTER(LEN=256)         :: mnp_chi_in_file
    INTEGER                    :: nk
    REAL(KIND=DP), ALLOCATABLE :: theta(:)
    REAL(KIND=DP), ALLOCATABLE :: omega_g(:)
    REAL(KIND=DP), ALLOCATABLE :: gamma_g(:)
    !
    ! Step field parameters
    !
    REAL(KIND=DP)      :: field_height
    REAL(KIND=DP)      :: field_centre
    REAL(KIND=DP)      :: field_width
    REAL(KIND=DP)      :: pulse_area
    !
    ! Pulse field parameters
    ! 
    REAL(KIND=DP)      :: pulse_start
    REAL(KIND=DP)      :: pulse_stop
    REAL(KIND=DP)      :: pulse_phase
    REAL(KIND=DP)      :: pulse_cycles
    !
    ! Runge-Kutta
    !
    REAL(KIND=DP) :: trange
    REAL(KIND=DP) :: Q_sqd_start
    REAL(KIND=DP) :: Q_sqd_end
    REAL(KIND=DP) :: Q_mnp_start
    REAL(KIND=DP) :: rk_step
    INTEGER       :: npts
    REAL(KIND=DP) :: out_pts
    INTEGER       :: check_pt
    !
    ! Directories
    !
    CHARACTER(LEN=256) :: log_file='tmp.log'
    CHARACTER(LEN=256) :: in_folder
    CHARACTER(LEN=256) :: in_file
    CHARACTER(LEN=256) :: out_folder
    CHARACTER(LEN=256) :: out_file
    CHARACTER(LEN=256) :: params_file
    CHARACTER(LEN=256) :: jname
    CHARACTER(LEN=256) :: timestamp
    !
    INTEGER                       :: num_lev
    INTEGER                       :: npos
    REAL(KIND=DP), ALLOCATABLE    :: en(:)
    REAL(KIND=DP), ALLOCATABLE    :: gma(:,:)
    REAL(KIND=DP), ALLOCATABLE    :: big_gma(:,:)
    INTEGER, ALLOCATABLE          :: positions(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_0(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_eq(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: mu(:,:)

    !
    ! MPI Params
    !
    INTEGER          :: npts_per_proc
    INTEGER          :: nprocs
    INTEGER          :: proc_id
    INTEGER          :: remainder
    CHARACTER(LEN=4) :: proc_name
    INTEGER          :: vers_mpi1
    INTEGER          :: vers_mpi2
    
    ! Post-Proc Params
    INTEGER                    :: num_freqs
    INTEGER                    :: read_first_npts
    REAL(KIND=DP), ALLOCATABLE :: freqs_in(:)
    REAL(KIND=DP)              :: start_from
    REAL(KIND=DP)              :: go_to
    REAL(KIND=DP)              :: probe_freq
    INTEGER                    :: max_order
    LOGICAL                    :: c_diffs
    LOGICAL                    :: use_max_freq
    INTEGER                    :: read_col
    CHARACTER(LEN=256)         :: field_in_file
    CHARACTER(LEN=256)         :: field_out_file
    CHARACTER(LEN=256)         :: freqs_out_file

CONTAINS

SUBROUTINE get_params
    USE double
    USE global_params , ONLY : intens_par, length_par, au_to_ev, std_err, pi
    USE print_mod
    IMPLICIT NONE

    CALL read_in_file_rho
    IF (in_folder == '-#error') THEN
        WRITE(std_err,*) 'No input folder given. Exiting...'
        CALL EXIT(0)
    ENDIF
    CALL read_matrices
    CALL read_chi_in_file

    IF ( Q_sqd_end == 0.0_DP ) THEN
        Q_sqd_end = trange
    ENDIF


!    IF ( coupled ) THEN
        IF ( dist < 0.0_DP) THEN
            IF ( s_alpha == 2 ) THEN
                dist = s_rad + surf_dist + ratio*rad
            ELSE
                dist = s_rad + surf_dist + rad
            ENDIF
        ENDIF
        en = en/au_to_ev

        dist     = dist/length_par
        rad      = rad/length_par
        s_rad      = s_rad/length_par
        eps_eff1 = (2.0_DP*eps_0 + eps_s)/(3.0_DP*eps_0)

        mu = mu/eps_eff1

        omega_g    = omega_g/au_to_ev
        gamma_g    = gamma_g/au_to_ev
        theta      = ratio*(rad**3)*theta/au_to_ev
!    ELSE
!        eps_eff1 = 1.0_DP
!    ENDIF

    IF ( pulse_stop < 0.0_DP ) THEN
        pulse_stop = trange
    ENDIF
    omega      = omega/au_to_ev
    omega_control = omega_control/au_to_ev
    E0         = SQRT(I0 / intens_par)
    E0_control = SQRT(I0_control / intens_par)

    IF ( field_change_param == 'omega' ) THEN
        field_param_from = field_param_from/au_to_ev
        field_param_to   = field_param_to/au_to_ev
    ENDIF

    IF ( field_param_npts < 1 ) THEN
        field_param_npts = 1
    ENDIF

    pulse_area = pulse_area * pi
    IF ( field_change_param == 'pulse_area' ) THEN
        field_param_from = field_param_from*pi
        field_param_to   = field_param_to*pi
    ENDIF

    npts_per_proc = field_param_npts/nprocs
    remainder     = MOD(field_param_npts, nprocs)
    IF ( field_param_npts == 1 ) THEN
        field_param_step = 0.0_DP
    ELSE
        field_param_step = (field_param_to - field_param_from)/REAL(field_param_npts-1)
    ENDIF

    npts        = NINT(trange/rk_step)

    IF ( out_pts < 0 ) THEN
        check_pt = 1
    ELSE
        check_pt    = NINT(REAL(npts)/out_pts)
    ENDIF


    IF ( proc_id == 0 ) THEN
        CALL write_log
    ENDIF


END SUBROUTINE get_params

SUBROUTINE write_log
    USE double
    USE print_mod
    USE global_params , ONLY : length_par
    IMPLICIT NONE

    CHARACTER(LEN=256) :: cwd
    CHARACTER(LEN=3) :: istr
    CHARACTER(LEN=12)  :: today
    CHARACTER(LEN=12)  :: now
    CHARACTER(LEN=80) :: str_mpi_info
    INTEGER, PARAMETER :: fh_log = 20
    INTEGER :: i

    ! Get current working directory
    CALL GETCWD(cwd)

    CALL DATE_AND_TIME(DATE=today, TIME=now)
    today = today(7:8)//'/'//today(5:6)//'/'//today(1:4)
    now = now(1:2)//':'//now(3:4)//':'//now(5:6)
    timestamp = TRIM(today)//', '//TRIM(now)

    WRITE(str_mpi_info, "('Running MPI Version: ',I1,'.',I1,   &
           &  ' on ',I4,' processor(s).')") vers_mpi1, vers_mpi2, nprocs

    OPEN(fh_log, FILE=log_file, STATUS='REPLACE', ACTION='WRITE')

    CALL print_break(fh_log)
    CALL print_str("Job '"//TRIM(jname)//"' executed on "//TRIM(today)//' at '//TRIM(now)//'.', fh_log)
    CALL print_str(str_mpi_info, fh_log)
    CALL print_break(fh_log)
    WRITE(fh_log, *) 

    CALL print_title('Runge-Kutta Parameters', fh_log)
        CALL print_str_num_real('Total number of points used', REAL(npts, KIND=DP), fh_log)
        CALL print_str_num_real('RK time-step (a.u.)', rk_step, fh_log)
        CALL print_str_num_real('Total Propagation Time (a.u)', trange, fh_log)
        WRITE(fh_log, *)

    IF ( coupled ) THEN
    CALL print_title('SQD-MNP Properties', fh_log)
        DO i = 1,nk
            WRITE(istr,'(I3)') i
            CALL print_str_num_real('theta_g_'//ADJUSTL(istr)//' (a.u.)', theta(i), fh_log)
            CALL print_str_num_real('omega_g_'//ADJUSTL(istr)//' (a.u.)', omega_g(i), fh_log)
            CALL print_str_num_real('gamma_g_'//ADJUSTL(istr)//' (a.u.)', gamma_g(i), fh_log)
        ENDDO
        CALL print_str_num_real('Separation (a.u.)', dist, fh_log)
        CALL print_str_num_real('MNP Radius (a.u.)', rad, fh_log)
        CALL print_str_num_real('SQD Radius (a.u.)', s_rad, fh_log)
        CALL print_str_num_real('s_alpha', s_alpha, fh_log)
        CALL print_str_num_real('eps_0', eps_0, fh_log)
        CALL print_str_num_real('ratio', ratio, fh_log)
        CALL print_str_num_real('eps_eff1', eps_eff1, fh_log)
        CALL print_str_num_real('eps_eff2', eps_eff2, fh_log)
        WRITE(fh_log, *)
    ENDIF

    CALL print_title('Input System Variables', fh_log)
        CALL print_str('Energy Levels:', fh_log)
        WRITE(fh_log, *)
        CALL print_vec_real(en, fh_log)
        WRITE(fh_log, *)

        CALL print_str('Gamma:', fh_log)
        WRITE(fh_log, *)
        CALL print_mat_real(gma, fh_log)
        WRITE(fh_log, *)

        CALL print_str('mu:', fh_log)
        WRITE(fh_log, *)
        CALL print_mat_complex(mu, fh_log)
        WRITE(fh_log, *)

        CALL print_str('rho_init:', fh_log)
        WRITE(fh_log, *)
        CALL print_mat_complex(rho_0, fh_log)
        WRITE(fh_log, *)

        CALL print_str('rho_eq:', fh_log)
        WRITE(fh_log, *)
        CALL print_mat_complex(rho_eq, fh_log)
        WRITE(fh_log, *)

    CALL print_title('Input file used:', fh_log)

    CLOSE(fh_log)

    CALL SYSTEM('cat '//TRIM(in_file)//' >> '//log_file)

    OPEN(fh_log, FILE=log_file, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
        WRITE(fh_log, *)
        CALL print_eof(fh_log)
    CLOSE(fh_log)

    log_file    = TRIM(out_folder)//'/'//TRIM(jname)//'.log'
    ! Moving the log file to the output directory.
    CALL SYSTEM('mv tmp.log '//log_file)

END SUBROUTINE write_log

SUBROUTINE read_chi_in_file
    USE double
    USE print_mod , ONLY : print_str
    IMPLICIT NONE
    INTEGER :: i
    CHARACTER(LEN=3) :: tmpstr

    OPEN(UNIT=10, FILE=mnp_chi_in_file, STATUS='OLD', ACTION='READ')
        DO i = 1, 1000
            READ(10, '(3A)') tmpstr
            IF ( tmpstr == 'n_k') THEN
                EXIT
            ELSEIF ( i == 1000 ) THEN
                CALL print_str('Error: Problem reading chi input file "'//mnp_chi_in_file//'". Exiting...')
                CALL EXIT(1)
            ENDIF
        ENDDO
!        READ(10, *)
        READ(10, *) nk
            ALLOCATE(theta(nk))
            ALLOCATE(gamma_g(nk))
            ALLOCATE(omega_g(nk))
        READ(10, *) 
        READ(10, *) 
        READ(10, *) (omega_g(i) , i = 1,nk)
        READ(10, *) 
        READ(10, *) 
        READ(10, *) (gamma_g(i) , i = 1,nk)
        READ(10, *) 
        READ(10, *) 
        READ(10, *) (theta(i) , i = 1,nk)
    CLOSE(10)

END SUBROUTINE read_chi_in_file


SUBROUTINE read_matrices
    USE double
    USE num_lines , ONLY : numlines
    IMPLICIT NONE
    ! Names of input matrices.
    CHARACTER(LEN=256) :: rho_in, en_in, gma_in, big_gma_in, mu_in, pos_in,   &
                          rho_eq_in
    ! Dummy sum variables.
    INTEGER            :: i, j
    
    ! Reading in matrices
    rho_in      = TRIM(in_folder)//'/rho.txt'
    CALL check_file(rho_in)
    en_in       = TRIM(in_folder)//'/en.txt'
    CALL check_file(en_in)
    gma_in      = TRIM(in_folder)//'/gma.txt'
    CALL check_file(gma_in)
    big_gma_in  = TRIM(in_folder)//'/big_gma.txt'
    CALL check_file(big_gma_in)
    mu_in       = TRIM(in_folder)//'/mu.txt'
    CALL check_file(mu_in)
    pos_in      = TRIM(in_folder)//'/positions.txt'
    CALL check_file(pos_in)
    rho_eq_in   = TRIM(in_folder)//'/rho_eq.txt'
    CALL check_file(rho_eq_in)

    ! Reading number of levels in system by the number of lines in the rho_in
    ! matrix. The number of rho elements to be read is npos.
    num_lev = numlines(rho_in)
    npos = numlines(pos_in)

    ! Reading in matrices.
    ALLOCATE(                                                                 &
            rho_0(num_lev,num_lev),                                           &
            rho_eq(num_lev,num_lev),                                          &
            en(num_lev),                                                      &
            big_gma(num_lev,num_lev),                                         &
            gma(num_lev,num_lev),                                             &
            mu(num_lev,num_lev),                                              &
            positions(npos,2)                                                 &
            )

    OPEN(UNIT=10, FILE=rho_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) (rho_0(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=rho_eq_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) (rho_eq(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=en_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) en(i)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=gma_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) (gma(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=big_gma_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) (big_gma(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=mu_in, STATUS='OLD', ACTION='READ')
        DO i = 1, num_lev
            READ(10, *) (mu(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=pos_in, STATUS='OLD', ACTION='READ')
        DO i = 1, npos
            READ(10, *) (positions(i, j), j=1,2)
        ENDDO
    CLOSE(10)

END SUBROUTINE read_matrices

SUBROUTINE read_in_file_rho
    USE double
    IMPLICIT NONE

    ! INPUT-RELATED VARIABLES
    CHARACTER(LEN=256) :: buffer, label
    INTEGER            :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER            :: ios = 0
    INTEGER            :: line = 0
    CHARACTER(LEN=6)   :: line_out

    ! SET DEFAULTS
    I0           = 1.0_DP
    I0_control   = 1.0_DP
    field        = 'cosfield'
    trange    = 0.0_DP
    coupled = .FALSE.
    out_pts      = -1.0_DP
    in_folder    = '-#error'
    jname        = 'job'
    out_folder   = '.'
    pulse_area  = 0.0_DP
    field_centre  = 50.0_DP
    field_width   = 20.0_DP
    field_height  = 1.0E-8_DP
    pulse_phase  = 0.0_DP
    pulse_start  = 0.0_DP
    pulse_stop   = -1.0_DP
    pulse_cycles = 6.0_DP
    Q_sqd_start    = 0.0_DP
    Q_sqd_end    = 0.0_DP
    Q_mnp_start    = 0.0_DP

    field_param_npts = -1

    dist = -1.0_DP
    surf_dist = 10.0_DP
    rad = 7.5_DP
    s_rad = 1.0_DP
    s_alpha = 2.0_DP
    eps_0 = 1.0_DP
    ratio = 1.0_DP
    eps_s = 6.0_DP

    CALL check_file(in_file)

    OPEN(fh, FILE=in_file, STATUS='OLD', ACTION='READ')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.
    DO WHILE (ios == 0)
        READ(fh, '(A)', IOSTAT=ios) buffer
        IF (ios == 0) THEN
            line = line + 1

            ! Find the first instance of whitespace. Split label and data.
            pos = SCAN(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            SELECTCASE (label)

            CASE ('mnp_chi_in_file')
                READ(buffer, *, IOSTAT=ios) mnp_chi_in_file

            CASE ('eps_0')
                READ(buffer, *, IOSTAT=ios) eps_0

            CASE ('ratio')
                READ(buffer, *, IOSTAT=ios) ratio

            CASE ('eps_s')
                READ(buffer, *, IOSTAT=ios) eps_s

            CASE ('s_alpha')
                READ(buffer, *, IOSTAT=ios) s_alpha

            CASE ('s_rad')
                READ(buffer, *, IOSTAT=ios) s_rad

            CASE ('rad')
                READ(buffer, *, IOSTAT=ios) rad

            CASE ('dist')
                READ(buffer, *, IOSTAT=ios) dist

            CASE ('surf_dist')
                READ(buffer, *, IOSTAT=ios) surf_dist
                
            CASE ('field_change_param')
                READ(buffer, *, IOSTAT=ios) field_change_param

            CASE ('field_param_from')
                READ(buffer, *, IOSTAT=ios) field_param_from

            CASE ('field_param_to')
                READ(buffer, *, IOSTAT=ios) field_param_to

            CASE ('field_param_npts')
                READ(buffer, *, IOSTAT=ios) field_param_npts

            CASE ('omega')
                READ(buffer, *, IOSTAT=ios) omega

            CASE ('omega_control')
                READ(buffer, *, IOSTAT=ios) omega_control

            CASE ('I0')
                READ(buffer, *, IOSTAT=ios) I0

            CASE ('I0_control')
                READ(buffer, *, IOSTAT=ios) I0_control

            CASE ('field')
                READ(buffer, *, IOSTAT=ios) field
                field = TRIM(field)

            CASE ('trange')
                READ(buffer, *, IOSTAT=ios) trange

            CASE ('coupled')
                READ(buffer, *, IOSTAT=ios) coupled

            CASE ('Q_sqd_end')
                READ(buffer, *, IOSTAT=ios) Q_sqd_end

            CASE ('Q_sqd_start')
                READ(buffer, *, IOSTAT=ios) Q_sqd_start

            CASE ('Q_mnp_start')
                READ(buffer, *, IOSTAT=ios) Q_mnp_start

            CASE ('rk_step')
                READ(buffer, *, IOSTAT=ios) rk_step

            CASE ('out_pts')
                READ(buffer, *, IOSTAT=ios) out_pts

            CASE ('in_folder')
                READ(buffer, *, IOSTAT=ios) in_folder

            CASE ('out_folder')
                READ(buffer, *, IOSTAT=ios) out_folder
                CALL SYSTEM('mkdir -p '//TRIM(out_folder))

            CASE ('name')
                READ(buffer, *, IOSTAT=ios) jname

            CASE ('pulse_area')
                READ(buffer, *, IOSTAT=ios) pulse_area

            CASE ('field_height')
                READ(buffer, *, IOSTAT=ios) field_height

            CASE ('field_centre')
                READ(buffer, *, IOSTAT=ios) field_centre

            CASE ('field_width')
                READ(buffer, *, IOSTAT=ios) field_width

            CASE ('pulse_stop')
                READ(buffer, *, IOSTAT=ios) pulse_stop

            CASE ('pulse_start')
                READ(buffer, *, IOSTAT=ios) pulse_start

            CASE ('pulse_cycles')
                READ(buffer, *, IOSTAT=ios) pulse_cycles

            CASE ('pulse_phase')
                READ(buffer, *, IOSTAT=ios) pulse_phase

            CASE DEFAULT
                WRITE(line_out,'(I3)') line
                IF ( label(1:1) /= '#' .AND. label(1:1) /= '') THEN
                    CALL param_read_fail(label, line_out)
                ENDIF
            END SELECT
        END IF
    END DO

    CLOSE(fh)

END SUBROUTINE read_in_file_rho

SUBROUTINE get_params_pp
    USE double
    USE print_mod , ONLY : print_str
    USE global_params , ONLY : au_to_ev
    IMPLICIT NONE
    
    ! INPUT-RELATED VARIABLES
    CHARACTER(LEN=256) :: buffer, label
    INTEGER            :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER, PARAMETER :: logid=20
    INTEGER            :: ios = 0
    INTEGER            :: line = 0
    CHARACTER(LEN=8)   :: line_out
    CHARACTER(LEN=12)  :: now
    CHARACTER(LEN=12)  :: today
    CHARACTER(LEN=256) :: ignore(23)

    ! SET DEFAULTS
    start_from = 0.0_DP
    probe_freq = 0.0_DP
    go_to      = -1.0_DP
    max_order  = 1
    c_diffs    = .FALSE.
    use_max_freq = .FALSE.
    field_in_file = '-#error'
    out_folder = './out/'
    jname = 'job'
    num_freqs = -1
    read_col = 2
    read_first_npts = 500

    ! SETTING PARAMETERS
    
    ! Inputs to ignore in input file
    ignore  = (/'freq_max           ',                                                    &
                'freq_min           ',                                                    &
                'step               ',                                                        &
                'in_folder          ',                                                   &
                'field              ',                                                       &
                'trange             ',                                                   &
                'I0                 ',                                                          &
                'omega              ',                                                    &
                'chi_out            ',                                                     &
                'dist               ',                                                     &
                'rad                ',                                                      &
                's_alpha            ',                                                     &
                'eps_s              ',                                                       &
                'eps_0              ',                                                       &
                'theta              ',                                                       &
                'omega_g            ',                                                     &
                'm_h                ',                                              &
                'nk                 ',                                              &
                'p22_start          ',                                              &
                'rk_step            ',                                              &
                'out_pts            ',                                              &
                'freq_mult          ',                                              &
                'gamma_g            '                                                      &
              /)

    CALL DATE_AND_TIME(DATE=today, TIME=now)
    timestamp = TRIM(today)//', '//TRIM(now(1:6))

    CALL check_file(in_file)
    OPEN(fh, FILE=in_file, STATUS='OLD', ACTION='READ')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    OPEN(logid, FILE=log_file, STATUS='REPLACE')

    DO WHILE (ios == 0)
        READ(fh, '(A)', IOSTAT=ios) buffer
        IF (ios == 0) THEN
            line = line + 1

            ! Find the first instance of whitespace. Split label and data.
            pos = SCAN(buffer, ' ')
            label = buffer(1:pos)
            buffer = buffer(pos+1:)

            SELECT CASE (label)

            CASE ('read_first_npts')
                READ(buffer, *, IOSTAT=ios) read_first_npts
                CALL param_read_success('read_first_npts',logid)

            CASE ('read_col')
                READ(buffer, *, IOSTAT=ios) read_col
                CALL param_read_success('read_col',logid)

            CASE ('probe_freq')
                READ(buffer, *, IOSTAT=ios) probe_freq
                probe_freq = probe_freq/au_to_ev
                CALL param_read_success('probe_freq',logid)

            CASE ('freqs')
                READ(buffer, *, IOSTAT=ios) freqs_in
                freqs_in = freqs_in/au_to_ev
                CALL param_read_success('freqs',logid)

            CASE ('start_from')
                READ(buffer, *, IOSTAT=ios) start_from
                CALL param_read_success('start_from',logid)

            CASE ('num_freqs')
                READ(buffer, *, IOSTAT=ios) num_freqs
                CALL param_read_success('num_freqs',logid)
                ALLOCATE( freqs_in(num_freqs) )
                freqs_in = 0.0_DP

            CASE ('go_to')
                READ(buffer, *, IOSTAT=ios) go_to
                CALL param_read_success('go_to',logid)

            CASE ('max_order')
                READ(buffer, *, IOSTAT=ios) max_order
                CALL param_read_success('max_order',logid)

            CASE ('calc_diffs')
                READ(buffer, *, IOSTAT=ios) c_diffs
                CALL param_read_success('calc_diffs',logid)

            CASE ('use_max_freq')
                READ(buffer, *, IOSTAT=ios) use_max_freq
                CALL param_read_success('use_max_freq',logid)

            CASE ('out_folder')
                READ(buffer, *, IOSTAT=ios) out_folder
                CALL param_read_success('out_folder',logid)
                CALL SYSTEM('mkdir -p '//TRIM(out_folder))

            CASE ('name')
                READ(buffer, *, IOSTAT=ios) jname
                CALL param_read_success('name',logid)

            CASE ('in_file')
                READ(buffer, *, IOSTAT=ios) field_in_file
                CALL param_read_success('in_file',logid)
                CALL check_file(field_in_file)

            CASE DEFAULT
                WRITE(line_out,'(I3)') line
                IF ( label(1:1) /= '#' .AND. label(1:1) /= '') THEN
                    IF ( ALL(ignore .NE. label) ) THEN
                        WRITE(logid, *)
                        CALL param_read_fail(label, line_out)
                    ENDIF
                ENDIF
            END SELECT
        END IF
    END DO

    CLOSE(logid)

    CLOSE(fh)

    IF ( num_freqs == -1 ) THEN
        CALL print_str('Error: "num_freqs" not given. Exiting...')
        CALL EXIT(1)
    ENDIF

    IF ( use_max_freq .EQV. .FALSE. ) THEN
        IF ( ABS(go_to+1.0_DP) <= 1.0E-12_DP ) THEN
            CALL print_str('Error: either "go_to" must be set or &
                           &"use_max_freq" must be used. Exiting...')
            CALL EXIT(1)
        ENDIF
    ENDIF

    IF ( field_in_file == '-#error' ) THEN
        CALL print_str('Error: "in_file" not specified in '//            &
                        TRIM(in_file)//'.')
        CALL EXIT(1)
    ENDIF

    ! Getting file names
    field_out_file = TRIM(out_folder)//'/'//TRIM(jname)//'-pp.field'
    freqs_out_file = TRIM(out_folder)//'/'//TRIM(jname)//'-pp.freqs'
    log_file    = TRIM(out_folder)//'/'//TRIM(jname)//'-pp.log'

    ! Moving the log file to the output directory.
    CALL SYSTEM('mv tmp.log '//log_file)

END SUBROUTINE get_params_pp

SUBROUTINE param_read_success(pname, fid)
    USE print_mod , ONLY : print_str
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: pname
    INTEGER, INTENT(IN) :: fid

    CALL print_str('Read "'//TRIM(pname)//'" successfully.',fid)

END SUBROUTINE param_read_success

SUBROUTINE param_read_fail(label, line)
    USE global_params , ONLY : std_err
    USE print_mod , ONLY : print_str
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: label
    CHARACTER(LEN=*), INTENT(IN) :: line

    CALL print_str( 'Warning: bad label name in input file at line ' &
                    //TRIM(ADJUSTL(line))//'.')
    CALL print_str( '>> Skipping invalid label "'//                            &
              TRIM(label)//'".', std_err)

END SUBROUTINE param_read_fail

SUBROUTINE check_file(the_file)
    USE print_mod , ONLY : print_str
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: the_file
    LOGICAL :: in_exist

    INQUIRE(FILE=the_file, EXIST=in_exist)
    IF ( in_exist .EQV. .FALSE.) THEN
        CALL print_str('Error: file '//TRIM(the_file)//' does not exist. &
                       &Exiting...')
        CALL EXIT(1)
    ENDIF

END SUBROUTINE check_file

SUBROUTINE print_field_params
    USE double
    USE print_mod
    USE global_params , ONLY : au_to_ev, pi
    IMPLICIT NONE
    
    ! File handle
    INTEGER            :: fid

    ! Writing parameters to file
    fid = 11 + proc_id
    OPEN(fid, FILE=out_file, STATUS='OLD', POSITION='APPEND', ACTION='WRITE')

        CALL print_title('Field Properties', fid)
        WRITE(fid,*)
        CALL print_str_str('Field type', field, fid)

        SELECTCASE ( field )

        CASE ( 'gauss_pulse' )
            CALL print_str_num_real('> Pulse Area', pulse_area, fid)
            CALL print_str_num_real('> Height', 0.939437278699651_DP*pulse_area/REAL(mu(1,2),KIND=DP)/field_width, fid)
            CALL print_str_num_real('> Full Width at Half Maximum', field_width, fid)
            CALL print_str_num_real('> Centre', field_centre, fid)
            CALL print_str_num_real('> Omega (eV)', omega*au_to_ev, fid)

        CASE ( 'sech_pulse' )
            CALL print_str_num_real('> Height', field_height, fid)
            CALL print_str_num_real('> Width', field_width, fid)
            CALL print_str_num_real('> Centre', field_centre, fid)
            CALL print_str_num_real('> Omega (eV)', omega*au_to_ev, fid)

        CASE ( 'gauss' )
            CALL print_str_num_real('> Pulse Area', pulse_area, fid)
            CALL print_str_num_real('> Height', 0.939437278699651_DP*pulse_area/REAL(mu(1,2),KIND=DP)/field_width, fid)
            CALL print_str_num_real('> Full Width at Half Maximum', field_width, fid)
            CALL print_str_num_real('> Centre', field_centre, fid)

        CASE ( 'step' )
            CALL print_str_num_real('> Height', field_height, fid)
            CALL print_str_num_real('> Width', field_width, fid)
            CALL print_str_num_real('> Centre', field_centre, fid)

        CASE ( 'pulse' )
            CALL print_str_num_real('> Pulse start', pulse_start, fid)
            CALL print_str_num_real('> Pulse phase', pulse_phase, fid)
            CALL print_str_num_real('> Number of cycles', pulse_cycles, fid)

        CASE DEFAULT
            CALL print_str_num_real('Laser Frequency (a.u.)', omega, fid)
            CALL print_str_num_real('Laser Energy (eV)', omega*au_to_ev, fid)
            CALL print_str_num_real('Laser Intenstiy (W/cm^2)', I0, fid)
            CALL print_str_num_real('Laser Amplitude (E0) (a.u.)', E0, fid)

        END SELECT

        WRITE(fid, *)
        CALL print_eof(fid)

    CLOSE(fid)

END SUBROUTINE print_field_params

END MODULE params
