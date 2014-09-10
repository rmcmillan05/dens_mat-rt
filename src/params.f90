MODULE params
    USE double
    IMPLICIT NONE
    
    ! UNIVERSAL PARAMETERS
    !
    ! Imaginary number
    COMPLEX(KIND=DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
    ! pi
    REAL(KIND=DP), PARAMETER    :: pi = 3.1415926535897932
    ! A.u. of wavelength when given in nm
    REAL(KIND=DP), PARAMETER    :: wave_par = 45.5633526_DP                            
    ! A.u. of intensity                                                        
    REAL(KIND=DP), PARAMETER    :: intens_par = 3.50944758E16_DP                  
    ! Convert laser energy in eV to wavelength in nm                           
    REAL(KIND=DP), PARAMETER    :: energy_par = 1.239841E3_DP                     
    ! A.u. conversion for time (fs)                                           
    REAL(KIND=DP), PARAMETER    :: time_par = 2.418884326505E-2_DP               
    ! A.u. conversion for length given in nm                                   
    REAL(KIND=DP), PARAMETER    :: length_par = 0.052917721092_DP

    ! USER INPUT VARIABLES FROM FILE WHICH GET SET WHEN get_params IS CALLED
    !
    ! Field parameters
    !
    ! Input field
    CHARACTER(LEN=256) :: field
    ! Laser freq in eV
    REAL(KIND=DP)      :: omega_ev
    ! Laser freq in a.u.
    REAL(KIND=DP)      :: omega_au
    ! Laser wavelength in nm
    REAL(KIND=DP)      :: lambda
    ! Laser intensity in W/cm2
    REAL(KIND=DP)      :: I0
    ! Laser amplitude in a.u.
    REAL(KIND=DP)      :: E0
    !
    ! Step field parameters
    !
    ! Height
    REAL(KIND=DP)      :: step_height
    ! Centre
    REAL(KIND=DP)      :: step_centre
    ! Width
    REAL(KIND=DP)      :: step_width
    !
    ! Pulse field parameters
    !
    ! Start
    REAL(KIND=DP)      :: pulse_start
    ! Phase
    REAL(KIND=DP)      :: pulse_phase
    ! Number of cycles
    REAL(KIND=DP)      :: pulse_cycles
    ! Limit
    REAL(KIND=DP)      :: pulse_lim
    !
    ! Runge-Kutta
    !
    ! Max propagation time in a.u.
    REAL(KIND=DP)      :: trange_au
    ! No of pts (nptspau+1) used in RK method per atomic unit
    REAL               :: nptspau
    ! RK step-size
    REAL(KIND=DP)      :: rk_step
    ! Total number of points used in propagation
    INTEGER            :: npts 
    !
    ! Directories
    !
    ! Input directory
    CHARACTER(LEN=256) :: in_folder
    ! Input file
    CHARACTER(LEN=256) :: in_file
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
    !
    ! Number of levels in system
    INTEGER                       :: num_lev
    ! Number of elements of rho to be output
    INTEGER                       :: npos
    ! Energy level vector
    REAL(KIND=DP), ALLOCATABLE    :: en(:)
    ! Small gamma
    REAL(KIND=DP), ALLOCATABLE    :: gma(:,:)
    ! Big gamma
    REAL(KIND=DP), ALLOCATABLE    :: big_gma(:,:)
    ! Matrix of rho elements to be output
    INTEGER, ALLOCATABLE          :: positions(:,:)
    ! Input rho_0 matrix
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_0(:,:)
    ! Input rho_eq matrix
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_eq(:,:)
    ! Input mut matrix
    COMPLEX(KIND=DP), ALLOCATABLE :: mu(:,:)
    
CONTAINS

SUBROUTINE get_params
    USE double
    IMPLICIT NONE
    
    ! INPUT-RELATED VARIABLES
    CHARACTER(LEN=256) :: buffer, label
    INTEGER            :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER, PARAMETER :: fh_log = 20
    INTEGER            :: ios = 0
    INTEGER            :: line = 0
    CHARACTER(LEN=6)   :: line_out
    CHARACTER(LEN=12)  :: now
    CHARACTER(LEN=12)  :: today
    LOGICAL            :: ev=.FALSE., nm=.FALSE., au=.FALSE.

    ! SET DEFAULTS
    I0           = 1.0_DP
    field        = 'cosfield'
    trange_au    = 1200.0_DP
    nptspau      = 100.0
    in_folder    = '-#error'
    jname        = 'job'
    out_folder   = '.'
    step_centre  = 50.0_DP
    step_width   = 20.0_DP
    step_height  = 1.0E-8_DP
    pulse_phase  = 0.0_DP
    pulse_start  = 0.0_DP
    pulse_cycles = 6.0_DP

    ! SETTING PARAMETERS

    CALL DATE_AND_TIME(DATE=today, TIME=now)
    timestamp = TRIM(today)//', '//TRIM(now(1:6))

    CALL check_file(in_file)

    OPEN(fh, FILE=in_file, STATUS='OLD', ACTION='READ')
    OPEN(fh_log, FILE=log_file, STATUS='REPLACE', ACTION='WRITE')

    WRITE(fh_log,*) TRIM(timestamp)

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

            CASE ('omega_ev')
                READ(buffer, *, IOSTAT=ios) omega_ev
                WRITE(fh_log,*) 'Read "omega_ev" successfully.'
                ev = .TRUE.; au = .FALSE.; nm=.FALSE.

            CASE ('omega_au')
                READ(buffer, *, IOSTAT=ios) omega_au
                WRITE(fh_log,*) 'Read "omega_au" successfully.'
                ev = .FALSE.; au = .TRUE.; nm=.FALSE.

            CASE ('lambda')
                READ(buffer, *, IOSTAT=ios) lambda
                WRITE(fh_log,*) 'Read "lambda" successfully.'
                ev = .FALSE.; au = .FALSE.; nm=.TRUE.

            CASE ('I0')
                READ(buffer, *, IOSTAT=ios) I0
                WRITE(fh_log,*) 'Read "I0" successfully.'

            CASE ('field')
                READ(buffer, *, IOSTAT=ios) field
                field = TRIM(field)
                WRITE(fh_log,*) 'Read "field" successfully.'

            CASE ('trange_au')
                READ(buffer, *, IOSTAT=ios) trange_au
                WRITE(fh_log,*) 'Read "trange_au" successfully.'

            CASE ('nptspau')
                READ(buffer, *, IOSTAT=ios) nptspau
                WRITE(fh_log,*) 'Read "nptspau" successfully.'

            CASE ('in_folder')
                READ(buffer, *, IOSTAT=ios) in_folder
                WRITE(fh_log,*) 'Read "in_folder" successfully.'

            CASE ('out_folder')
                READ(buffer, *, IOSTAT=ios) out_folder
                WRITE(fh_log,*) 'Read "out_folder" successfully.'
                CALL SYSTEM('mkdir -p '//TRIM(out_folder))

            CASE ('name')
                READ(buffer, *, IOSTAT=ios) jname
                WRITE(fh_log,*) 'Read "name" successfully.'

            CASE ('step_height')
                READ(buffer, *, IOSTAT=ios) step_height
                WRITE(fh_log,*) 'Read "step_height" successfully.'

            CASE ('step_centre')
                READ(buffer, *, IOSTAT=ios) step_centre
                WRITE(fh_log,*) 'Read "step_centre" successfully.'

            CASE ('step_width')
                READ(buffer, *, IOSTAT=ios) step_width
                WRITE(fh_log,*) 'Read "step_width" successfully.'

            CASE ('pulse_start')
                READ(buffer, *, IOSTAT=ios) pulse_start
                WRITE(fh_log,*) 'Read "pulse_start" successfully.'

            CASE ('pulse_cycles')
                READ(buffer, *, IOSTAT=ios) pulse_cycles
                WRITE(fh_log,*) 'Read "pulse_cycles" successfully.'

            CASE ('pulse_phase')
                READ(buffer, *, IOSTAT=ios) pulse_phase
                WRITE(fh_log,*) 'Read "pulse_phase" successfully.'

            CASE DEFAULT
                WRITE(line_out,'(I5)') line
                IF ( label(1:1) /= '#') THEN
                WRITE(fh_log,*) 'Error in file "'//TRIM(in_file)//'" at line '&
                                //ADJUSTL(TRIM(line_out))                     &
                                //'. Skipping invalid label "'                &
                                //ADJUSTL(TRIM(label)//'".')
                ENDIF
            END SELECT
        END IF
    END DO

    CLOSE(fh)

    IF (ev) THEN
        lambda = energy_par/omega_ev
        omega_au = wave_par/lambda
    ELSEIF (au) THEN
        lambda = wave_par/omega_au
        omega_ev = energy_par/lambda
    ELSEIF (nm) THEN
        omega_ev = energy_par/lambda
        omega_au = wave_par/lambda
    ELSE
        omega_ev = 2.0_DP
        lambda = energy_par/omega_ev
        omega_au = wave_par/lambda
    ENDIF

    CLOSE(fh_log)

    IF (in_folder == '-#error') THEN
        WRITE(*,*) 'No input folder given. Exiting...'
        CALL EXIT(0)
    ENDIF

    ! Field calculations
    E0 = SQRT(I0 / intens_par)
    pulse_lim = 2.0_DP * pi * pulse_cycles / omega_au
    
    ! Defining RK step size and npts from npstpau
    rk_step     = 1.0_DP/REAL(nptspau,KIND=DP)
    npts        = NINT(nptspau * trange_au)

    ! Getting file names
    params_file = TRIM(out_folder)//'/'//TRIM(jname)//'.params'
    out_file    = TRIM(out_folder)//'/'//TRIM(jname)//'.out'
    log_file    = TRIM(out_folder)//'/'//TRIM(jname)//'.log'

    ! Moving the log file to the output directory.
    CALL SYSTEM('mv tmp.log '//log_file)

END SUBROUTINE get_params

SUBROUTINE check_file(the_file)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: the_file
    LOGICAL :: in_exist

    INQUIRE(FILE=the_file, EXIST=in_exist)
    IF ( in_exist .EQV. .FALSE.) THEN
        WRITE(*,*) 'Error: file '//TRIM(the_file)//' does not exist. Exiting...'
        CALL EXIT(1)
    ENDIF

END SUBROUTINE check_file

END MODULE params
