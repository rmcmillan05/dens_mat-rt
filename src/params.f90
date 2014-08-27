MODULE params
    USE double
    IMPLICIT NONE
    
    ! imaginary number
    COMPLEX(KIND=DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
    ! pi
    REAL(KIND=DP), PARAMETER :: pi = 3.1415926535897932
    ! a.u. of wavelength when given in nm
    REAL(KIND=DP), PARAMETER :: chi = 45.5633526_DP                            
    ! a.u. of intensity                                                        
    REAL(KIND=DP), PARAMETER :: intens_par = 3.50944758E16_DP                  
    ! Convert laser energy in eV to wavelength in nm                           
    REAL(KIND=DP), PARAMETER :: energy_par = 1.239841E3_DP                     
    ! a.u. conversion for time (fs)                                           
    REAL(KIND=DP), PARAMETER :: time_par = 2.418884326505E-2_DP               
    ! A.U. CONVERSION FOR LENGTH GIVEN IN NM                                   
    REAL(KIND=DP), PARAMETER :: length_par = 0.052917721092_DP


    ! USER INPUT VARIABLES FROM FILE

    ! Input field
    CHARACTER(LEN=256) :: field
    ! Laser freq in eV
    REAL(KIND=DP) :: omega_ev
    ! Laser freq in a.u.
    REAL(KIND=DP) :: omega_au
    ! Laser wavelength in nm
    REAL(KIND=DP) :: lambda
    ! Laser intensity in W/cm2
    REAL(KIND=DP) :: I0
    ! Laser amplitude in a.u.
    REAL(KIND=DP) :: E0

    ! Max propagation time in a.u.
    REAL(KIND=DP) :: trange_au
    ! No of pts (nptspau+1) used in RK method per atomic unit
    REAL :: nptspau

    ! Input directory
    CHARACTER(LEN=256) :: in_folder
    ! Input file
    CHARACTER(LEN=256) :: in_file='rho_prop.in'
    ! Output directory
    CHARACTER(LEN=256) :: out_folder
    ! Output file
    CHARACTER(LEN=256) :: out_file
    ! Params file
    CHARACTER(LEN=256) :: params_file
    ! Job name
    CHARACTER(LEN=256) :: jname
    ! Error file
    CHARACTER(LEN=256) :: log_file='tmp.log'
    ! Time stamp
    CHARACTER(LEN=256) :: timestamp

    INTEGER :: num_lev
    INTEGER :: npos
    REAL(KIND=DP), ALLOCATABLE :: en(:)
    REAL(KIND=DP), ALLOCATABLE :: gma(:,:)
    REAL(KIND=DP), ALLOCATABLE :: big_gma(:,:)
    INTEGER, ALLOCATABLE :: positions(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: rho_0(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: mu(:,:)
    
CONTAINS

SUBROUTINE get_params(omega_ev,   &
                      omega_au,   &
                      lambda,     &
                      I0,         &
                      E0,         &
                      field,      &
                      trange_au,  &
                      nptspau,       &
                      in_folder,  &
                      out_folder,  &
                      timestamp,  &
                      jname)
    USE double
    IMPLICIT NONE
    
    ! INPUT-RELATED VARIABLES
    CHARACTER(LEN=256) :: buffer, label
    INTEGER :: pos
    INTEGER, PARAMETER :: fh = 15
    INTEGER, PARAMETER :: fh_log = 20
    INTEGER :: ios = 0
    INTEGER :: line = 0
    CHARACTER(LEN=6) :: line_out
    INTEGER :: today(3), now(3)

    ! CONTROL FILE VARIABLES
    REAL(KIND=DP), INTENT(OUT) :: omega_ev
    REAL(KIND=DP), INTENT(OUT) :: omega_au
    REAL(KIND=DP), INTENT(OUT) :: lambda
    REAL(KIND=DP), INTENT(OUT) :: I0
    REAL(KIND=DP), INTENT(OUT) :: E0
    CHARACTER(LEN=256), INTENT(OUT) :: field
    REAL(KIND=DP), INTENT(OUT) :: trange_au
    REAL, INTENT(OUT) :: nptspau
    CHARACTER(LEN=256), INTENT(OUT) :: in_folder
    CHARACTER(LEN=256), INTENT(OUT) :: out_folder
    CHARACTER(LEN=256), INTENT(OUT) :: jname
    CHARACTER(LEN=256), INTENT(OUT) :: timestamp

    LOGICAL :: ev=.FALSE., nm=.FALSE., au=.FALSE.

    ! SET DEFAULTS
    I0 = 1.0_DP
    field = 'cosfield'
    trange_au = 1200.0_DP
    nptspau = 100.0
    in_folder = '-#error'
    jname = 'job'
    out_folder = '.'

    !! THE PROGRAM !!

    CALL IDATE(today)
    CALL ITIME(now)

    WRITE(timestamp,'("DATE: ",I2,"/",I2,"/",I4,", ",&
                     &"TIME: ",I2,":",I2,":",I2)') today, now

    OPEN(fh, FILE=in_file)
    OPEN(fh_log, FILE=log_file, STATUS='REPLACE')

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
        omega_au = chi/lambda
    ELSEIF (au) THEN
        lambda = chi/omega_au
        omega_ev = energy_par/lambda
    ELSEIF (nm) THEN
        omega_ev = energy_par/lambda
        omega_au = chi/lambda
    ELSE
        omega_ev = 2.0_DP
        lambda = energy_par/omega_ev
        omega_au = chi/lambda
    ENDIF

    CLOSE(fh_log)

    IF (in_folder == '-#error') THEN
        WRITE(*,*) 'No input folder given. Exiting...'
        CALL EXIT(0)
    ENDIF

    E0 = SQRT(I0 / intens_par)
    
END SUBROUTINE get_params

END MODULE params
