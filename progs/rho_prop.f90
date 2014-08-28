PROGRAM rho_prop
    !
    ! Propogates an initial rho matrix in time with the matrices given in
    ! the specified input_directory, using the fourth-order Runge-Kutta method.
    !
    ! N.B. The runge routine only prints out the rho_i,j elements specified in
    ! the file 'positions.txt' located in input_directory.
    !
    USE double
    USE num_lines , ONLY : numlines
    USE runge_mod , ONLY : runge
    USE print_params_mod , ONLY : print_params
    USE params , ONLY : num_lev, npos, in_folder, in_file, rho_0, en, gma,    &
                        big_gma, mu, positions, rho_eq, get_params
    IMPLICIT NONE

    ! Names of input matrices.
    CHARACTER(LEN=256) :: rho_in, en_in, gma_in, big_gma_in, mu_in, pos_in,   &
                          rho_eq_in
    ! Used to calculate the total computation time.
    INTEGER            :: count1, count2, count_rate
    ! Dummy sum variables.
    INTEGER            :: i, j
    
    CALL SYSTEM_CLOCK(count1,count_rate)

    ! Take first input argument as the file from which variables are read 
    ! (default rho_prop.in).
    SELECT CASE ( IARGC() )
        CASE ( 1 )
            CALL GETARG(1, in_file)
        CASE DEFAULT
            in_file = 'rho_prop.in'
    END SELECT

    ! Calling get_params to read in variables from file and store them.
    CALL   get_params

    ! Reading in matrices
    rho_in      = TRIM(in_folder)//'/rho.txt'
    en_in       = TRIM(in_folder)//'/en.txt'
    gma_in      = TRIM(in_folder)//'/gma.txt'
    big_gma_in  = TRIM(in_folder)//'/big_gma.txt'
    mu_in       = TRIM(in_folder)//'/mu.txt'
    pos_in      = TRIM(in_folder)//'/positions.txt'
    rho_eq_in   = TRIM(in_folder)//'/rho_eq.txt'

    ! Reading number of levels in system by the number of lines in the rho_in
    ! matrix. The number of rho elements to be read is npos.
    num_lev = numlines(rho_in)
    npos = numlines(pos_in)

    ! Reading in matrices.
    ALLOCATE(                                                                 &
            rho_0(num_lev,num_lev),                                           &
            rho_eq(num_lev,num_lev),                                           &
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

    ! Performing the RK calculations.
    CALL runge
    ! Outputing all parameters used to file.
    CALL print_params

    ! Printing information about job to screen.
    WRITE(*,'(A1)',ADVANCE='NO') char(10)
    WRITE(*,'(A28)') '**** Completed. Yay :-) ****'
    CALL SYSTEM_CLOCK(count2,count_rate)
    WRITE(*,*)
    WRITE(*,'(A17,F16.2,A9)') 'Total time taken: ' ,                          &
                               (REAL(count2-count1)/REAL(count_rate)),        &
                              ' seconds.'

END PROGRAM rho_prop
