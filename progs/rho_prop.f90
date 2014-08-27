PROGRAM rho_prop
    !
    ! Propogates an initial rho matrix in time with the matrices given in
    ! the specified input_directory, using the fourth-order Runge-Kutta method.
    !
    ! N.B. The runge routine only prints out the rho_i,j elements specified in
    ! the file 'positions.txt' located in input_directory.
    !
    USE double
    USE num_lines
    USE runge_mod
    USE print_params_mod
    USE params
    IMPLICIT NONE

    INTEGER :: count1, count2, count_rate
    INTEGER           :: i, j
    CHARACTER(LEN=64) :: rho_in, en_in, gma_in, big_gma_in, mu_in, pos_in
    
    CALL SYSTEM_CLOCK(count1,count_rate)

    SELECT CASE ( IARGC() )
        CASE ( 1 )
            CALL GETARG(1, in_file)
    END SELECT

    CALL get_params(omega_ev,   &
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

    params_file = TRIM(out_folder)//'/'//TRIM(jname)//'.params'
    out_file = TRIM(out_folder)//'/'//TRIM(jname)//'.out'

    rho_in=TRIM(in_folder)//'/rho.txt'
    en_in=TRIM(in_folder)//'/en.txt'
    gma_in=TRIM(in_folder)//'/gma.txt'
    big_gma_in=TRIM(in_folder)//'/big_gma.txt'
    mu_in =TRIM(in_folder)//'/mu.txt'
    pos_in=TRIM(in_folder)//'/positions.txt'

    num_lev = numlines(rho_in)
    ALLOCATE(                         &
            rho_0(num_lev,num_lev),   &
            en(num_lev),              &
            big_gma(num_lev,num_lev), &
            gma(num_lev,num_lev),     &
            mu(num_lev,num_lev)       &
            )

    OPEN(UNIT=10, FILE=rho_in, STATUS='OLD')
        DO i = 1, num_lev
            READ(10, *) (rho_0(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=en_in, STATUS='OLD')
        DO i = 1, num_lev
            READ(10, *) en(i)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=gma_in, STATUS='OLD')
        DO i = 1, num_lev
            READ(10, *) (gma(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=big_gma_in, STATUS='OLD')
        DO i = 1, num_lev
            READ(10, *) (big_gma(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=mu_in, STATUS='OLD')
        DO i = 1, num_lev
            READ(10, *) (mu(i, j), j=1,num_lev)
        ENDDO
    CLOSE(10)
    
    npos = numlines(pos_in)
    ALLOCATE(positions(npos, 2))
    
    OPEN(UNIT=10, FILE=pos_in, STATUS='OLD')
        DO i = 1, npos
            READ(10, *) (positions(i, j), j=1,2)
        ENDDO
    CLOSE(10)

    CALL runge
    CALL print_params
    CALL SYSTEM('mv tmp.log '//TRIM(out_folder)//TRIM(jname)//'.log')

    WRITE(*,'(A1)',ADVANCE='NO') char(10)
    WRITE(*,'(A28)') '**** Completed. Yay :-) ****'

    CALL SYSTEM_CLOCK(count2,count_rate)

    WRITE(*,*)
    WRITE(*,'(A17,F16.2,A9)') 'Total time taken: ' ,                       &
                           (REAL(count2-count1)/REAL(count_rate)),' seconds.'

END PROGRAM rho_prop
