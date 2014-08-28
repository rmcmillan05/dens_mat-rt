MODULE print_params_mod
    USE double
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: print_params

CONTAINS

SUBROUTINE print_params
    USE double
    USE params , ONLY : out_folder, omega_ev, omega_au, lambda, I0, E0,       &
                        jname, in_folder, timestamp, in_file, en, num_lev,    &
                        params_file, trange_au, field, nptspau, rk_step, npts,&
                        step_centre, step_width, step_height,                 &
                        pulse_start, pulse_cycles, pulse_phase
    IMPLICIT NONE
    
    ! Format specifiers used
    CHARACTER(LEN=256) :: char_fmat  = '(A30)'
    CHARACTER(LEN=256) :: title_fmat = '(A45)'
    CHARACTER(LEN=256) :: dir_fmat   = '(A35)'
    CHARACTER(LEN=256) :: real_fmat  = '(ES16.8)'
    CHARACTER(LEN=256) :: int_fmat   = '(I16)'
    ! Current working directory
    CHARACTER(LEN=256) :: cwd
    ! Dummy index variables
    INTEGER            :: i, j, s

    ! Get current working directory
    CALL GETCWD(cwd)

    ! Writing parameters to file
    OPEN(UNIT=10, FILE=params_file, STATUS='REPLACE')

        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) TRIM(ADJUSTL(timestamp))
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 'Laser Properties                              '
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat)
        WRITE(10,char_fmat,ADVANCE='NO') 'Field type                  : '
        WRITE(10,char_fmat) field
        IF (field == 'step') THEN
            WRITE(10,char_fmat,ADVANCE='NO') '> Step height               : '
            WRITE(10,real_fmat) step_height
            WRITE(10,char_fmat,ADVANCE='NO') '> Step width                : '
            WRITE(10,real_fmat) step_width
            WRITE(10,char_fmat,ADVANCE='NO') '> Step centre               : '
            WRITE(10,real_fmat) step_centre
        ELSEIF (field == 'pulse') THEN
            WRITE(10,char_fmat,ADVANCE='NO') '> Pulse start               : '
            WRITE(10,real_fmat) pulse_start
            WRITE(10,char_fmat,ADVANCE='NO') '> Pulse phase               : '
            WRITE(10,real_fmat) pulse_phase
            WRITE(10,char_fmat,ADVANCE='NO') '> Number of cycles          : '
            WRITE(10,real_fmat) pulse_cycles
        ENDIF
        IF (field /= 'step' .AND. field /= 'zero') THEN
            WRITE(10,char_fmat,ADVANCE='NO') 'Laser Wavelength (nm)       : '
            WRITE(10,real_fmat) lambda
            WRITE(10,char_fmat,ADVANCE='NO') 'Laser Frequency (a.u.)      : '
            WRITE(10,real_fmat) omega_au
            WRITE(10,char_fmat,ADVANCE='NO') 'Laser Energy (eV)           : '
            WRITE(10,real_fmat) omega_ev
            WRITE(10,char_fmat,ADVANCE='NO') 'Laser Intensity (W/cm^2)    : '
            WRITE(10,real_fmat) I0
            WRITE(10,char_fmat,ADVANCE='NO') 'Laser Amplitude (E0) (a.u.) : '
            WRITE(10,real_fmat) E0
        ENDIF

        WRITE(10,title_fmat)
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 'Runge-Kutta Parameters                        '
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat)
        WRITE(10,char_fmat,ADVANCE='NO') 'Number of Points per a.u.   : '
        WRITE(10,real_fmat) nptspau
        WRITE(10,char_fmat,ADVANCE='NO') 'Total number of points used : '
        WRITE(10,int_fmat) npts
        WRITE(10,char_fmat,ADVANCE='NO') 'RK time-step (a.u.)         : '
        WRITE(10,real_fmat) rk_step
        WRITE(10,char_fmat,ADVANCE='NO') 'Total Propagation Time (a.u): '
        WRITE(10,real_fmat) trange_au
        WRITE(10,title_fmat) 
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 'Energy Levels and Differences                 '
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 

        IF ( ABS(en(1) - 0.0_DP) < 1.0E-16 ) THEN
            s = 2
        ELSE
            s = 1
        ENDIF
        
        DO i=s,num_lev
            DO j=s,i
                IF ( i == j ) THEN
                    WRITE(10,real_fmat) REAL(en(i),KIND=DP)
                ELSE
                    WRITE(10,real_fmat) REAL(ABS(en(i) - en(j)),KIND=DP)
                ENDIF
            ENDDO
        ENDDO

        WRITE(10,title_fmat)
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 'Directories                                   '
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat)
        WRITE(10,char_fmat)'Root folder:                  '
        WRITE(10,dir_fmat) cwd
        WRITE(10,*)
        WRITE(10,char_fmat)'Input folder:                 '
        WRITE(10,dir_fmat) in_folder
        WRITE(10,*)
        WRITE(10,char_fmat)'Output folder:                '
        WRITE(10,dir_fmat) out_folder
        WRITE(10,*)
        WRITE(10,char_fmat)'Output file prefix:           '
        WRITE(10,dir_fmat) jname

        WRITE(10,title_fmat)
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) 'Input file used                               '
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat)
        WRITE(10,dir_fmat) in_file
        WRITE(10,title_fmat)

    CLOSE(10)

    CALL SYSTEM('cat '//TRIM(in_file)//' >> '//params_file)
    OPEN(UNIT=10, FILE=params_file, STATUS='OLD', POSITION='APPEND')

        WRITE(10,title_fmat)
        WRITE(10,title_fmat) '----------------------------------------------'
        WRITE(10,title_fmat) '                   - EOF -                    '
        WRITE(10,title_fmat) '----------------------------------------------'

    CLOSE(10)

END SUBROUTINE print_params

END MODULE print_params_mod
