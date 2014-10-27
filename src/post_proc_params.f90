MODULE post_proc_params
    USE double
    USE global_params , ONLY : full_wd
    IMPLICIT NONE
    ! USER INPUT VARIABLES FROM FILE WHICH GET SET WHEN get_params IS CALLED
    !
    INTEGER :: num_freqs
    REAL(KIND=DP), ALLOCATABLE :: freqs_in(:)
    REAL(KIND=DP):: start_from
    REAL(KIND=DP) :: go_to
    REAL(KIND=DP) :: probe_freq
    INTEGER :: max_order
    LOGICAL :: c_diffs
    LOGICAL :: use_max_freq
    INTEGER :: read_col
    !
    ! Directories
    !
    ! Input file
    CHARACTER(LEN=256) :: in_file
    ! Input file containing t and field(t)
    CHARACTER(LEN=256) :: field_in_file
    ! Output file containing t, field(t), and the approximate field
    CHARACTER(LEN=256) :: field_out_file
    ! Output file containing the frequencies and real/imag amplitudes
    CHARACTER(LEN=256) :: freqs_out_file
    ! Output directory
    CHARACTER(LEN=256) :: out_folder
    ! Job name
    CHARACTER(LEN=256) :: jname
    ! Log file
    CHARACTER(LEN=256) :: log_file='tmp.log'
    ! Time stamp
    CHARACTER(LEN=256) :: timestamp

CONTAINS

SUBROUTINE get_params_pp
    USE double
    USE print_mat_mod , ONLY : print_str
    USE global_params , ONLY : length_par
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
    CHARACTER(LEN=256) :: ignore(19)

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

    ! SETTING PARAMETERS
    
    ! Inputs to ignore in input file
    ignore  = (/'freq_max           ',                                                    &
                'freq_min           ',                                                    &
                'step               ',                                                        &
                'in_folder          ',                                                   &
                'field              ',                                                       &
                'trange_au          ',                                                   &
                'nptspau            ',                                                     &
                'I0                 ',                                                          &
                'omega_ev           ',                                                    &
                'chi_out            ',                                                     &
                'dist_nm            ',                                                     &
                'rad_nm             ',                                                      &
                's_alpha            ',                                                     &
                'eps_s              ',                                                       &
                'eps_0              ',                                                       &
                'theta              ',                                                       &
                'omega_g            ',                                                     &
                'm_h                ',                                              &
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

            CASE ('read_col')
                READ(buffer, *, IOSTAT=ios) read_col
                CALL param_read_success('read_col',logid)

            CASE ('probe_freq')
                READ(buffer, *, IOSTAT=ios) probe_freq
                CALL param_read_success('probe_freq',logid)

            CASE ('freqs')
                READ(buffer, *, IOSTAT=ios) freqs_in
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
    USE print_mat_mod , ONLY : print_str
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: pname
    INTEGER, INTENT(IN) :: fid

    CALL print_str('Read "'//TRIM(pname)//'" successfully from file '//        &
                    TRIM(in_file)//'.', fid)

END SUBROUTINE param_read_success

SUBROUTINE param_read_fail(label, line)
    USE global_params , ONLY : std_err
    USE print_mat_mod , ONLY : print_str
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: label
    CHARACTER(LEN=*), INTENT(IN) :: line

    CALL print_str( 'Warning: bad label name in "'//TRIM(in_file)//            &
              '" at line '//TRIM(ADJUSTL(line))//'.')
    CALL print_str( '>> Skipping invalid label "'//                            &
              TRIM(label)//'".', std_err)

END SUBROUTINE param_read_fail

SUBROUTINE check_file(the_file)
    USE print_mat_mod , ONLY : print_str
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

END MODULE post_proc_params
