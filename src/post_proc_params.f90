MODULE post_proc_params
    USE double
    IMPLICIT NONE
    ! USER INPUT VARIABLES FROM FILE WHICH GET SET WHEN get_params IS CALLED
    !
    INTEGER :: num_freqs
    REAL(KIND=DP), ALLOCATABLE :: freqs_in(:)
    REAL(KIND=DP):: start_from
    REAL(KIND=DP) :: go_to
    INTEGER :: max_order
    LOGICAL :: c_diffs
    LOGICAL :: use_max_freq
    ! Directories
    !
    ! Input file
    CHARACTER(LEN=256) :: in_file
    ! Input file containing t and field(t)
    CHARACTER(LEN=256) :: field_in_file
    ! Output directory
    CHARACTER(LEN=256) :: out_folder
    ! Output file
    CHARACTER(LEN=256) :: out_file
    ! Parameters file
    CHARACTER(LEN=256) :: params_file
    ! Job name
    CHARACTER(LEN=256) :: jname
    ! Log file
    CHARACTER(LEN=256) :: log_file='tmp.log'
    ! Time stamp
    CHARACTER(LEN=256) :: timestamp
    
CONTAINS

SUBROUTINE get_params_pp
    USE double
    USE global_params , ONLY : std_err
    IMPLICIT NONE
    
    ! INPUT-RELATED VARIABLES
    CHARACTER(LEN=256) :: buffer, label
    INTEGER            :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER            :: ios = 0
    INTEGER            :: line = 0
    CHARACTER(LEN=6)   :: line_out
    CHARACTER(LEN=12)  :: now
    CHARACTER(LEN=12)  :: today

    ! SET DEFAULTS
    start_from = 0.0_DP
    go_to      = -1.0_DP
    max_order  = 1
    c_diffs    = .FALSE.
    use_max_freq = .FALSE.
    in_file = 'pp.in'
    field_in_file = '-#error'
    out_folder = './out/'
    jname = 'job'
    num_freqs = -1

    ! SETTING PARAMETERS

    CALL DATE_AND_TIME(DATE=today, TIME=now)
    timestamp = TRIM(today)//', '//TRIM(now(1:6))

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

            SELECT CASE (label)

            CASE ('num_freqs')
                READ(buffer, *, IOSTAT=ios) num_freqs
                CALL param_read_success('num_freqs',std_err)
                ALLOCATE( freqs_in(num_freqs) )
                freqs_in = 0.0_DP

            CASE ('freqs')
                READ(buffer, *, IOSTAT=ios) freqs_in
                CALL param_read_success('freqs',std_err)

            CASE ('start_from')
                READ(buffer, *, IOSTAT=ios) start_from
                CALL param_read_success('start_from',std_err)

            CASE ('go_to')
                READ(buffer, *, IOSTAT=ios) go_to
                CALL param_read_success('go_to',std_err)

            CASE ('max_order')
                READ(buffer, *, IOSTAT=ios) max_order
                CALL param_read_success('max_order',std_err)

            CASE ('calc_diffs')
                READ(buffer, *, IOSTAT=ios) c_diffs
                CALL param_read_success('calc_diffs',std_err)

            CASE ('use_max_freq')
                READ(buffer, *, IOSTAT=ios) use_max_freq
                CALL param_read_success('use_max_freq',std_err)

            CASE ('out_folder')
                READ(buffer, *, IOSTAT=ios) out_folder
                CALL param_read_success('out_folder',std_err)
                CALL SYSTEM('mkdir -p '//TRIM(out_folder))

            CASE ('name')
                READ(buffer, *, IOSTAT=ios) jname
                CALL param_read_success('name',std_err)

            CASE ('in_file')
                READ(buffer, *, IOSTAT=ios) field_in_file
                CALL param_read_success('in_file',std_err)
                CALL check_file(field_in_file)

            CASE DEFAULT
                WRITE(line_out,'(I5)') line
                IF ( label(1:1) /= '#') THEN
                WRITE(std_err,*) 'Warning: bad label name in "'//TRIM(in_file)//'" at line '&
                                //ADJUSTL(line_out)                     &
                                //'. Skipping invalid label "'                &
                                //ADJUSTL(TRIM(label)//'".')
                ENDIF
            END SELECT
        END IF
    END DO

    CLOSE(fh)

    IF ( num_freqs == -1 ) THEN
        WRITE(std_err,*) 'Error: number of frequencies not given. Exiting...' 
        CALL EXIT(1)
    ENDIF

    ! Getting file names
    params_file = TRIM(out_folder)//'/'//TRIM(jname)//'.params'
    out_file    = TRIM(out_folder)//'/'//TRIM(jname)//'.out'
    log_file    = TRIM(out_folder)//'/'//TRIM(jname)//'.log'

END SUBROUTINE get_params_pp

SUBROUTINE param_read_success(pname, output)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: pname
    INTEGER, INTENT(IN) :: output

    WRITE(output,*) 'Read "'//TRIM(pname)//'" successfully from file '//      &
                    TRIM(in_file)//'.'

END SUBROUTINE param_read_success

SUBROUTINE check_file(the_file)
    USE global_params , ONLY : std_err
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: the_file
    LOGICAL :: in_exist

    INQUIRE(FILE=the_file, EXIST=in_exist)
    IF ( in_exist .EQV. .FALSE.) THEN
        WRITE(std_err,*) 'Error: file '//TRIM(the_file)//' does not exist. Exiting...'
        CALL EXIT(1)
    ENDIF

END SUBROUTINE check_file

END MODULE post_proc_params
