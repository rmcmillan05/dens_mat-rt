PROGRAM rho_prop
    !
    ! Propogates an initial rho matrix in time with the matrices given in
    ! the specified input_directory, using the fourth-order Runge-Kutta method.
    !
    ! N.B. The runge routine only prints out the rho_i,j elements specified in
    ! the file 'positions.txt' located in input_directory.
    !
    USE double
    USE mpi
    USE communications
    USE runge_mod , ONLY : runge
    USE params , ONLY : npts_per_proc, nprocs, omega,  &
                        omega_from, remainder, proc_id, &
                        get_params, out_file, params_file, out_folder, jname, in_file, &
                        omega_step, field, print_field_params
    USE print_mod
    USE params , ONLY : check_file
    IMPLICIT NONE

    CHARACTER(LEN=8) :: suffix
    INTEGER :: j

    CALL set_up_mpi

    ! Take first input argument as the file from which variables are read 
    ! (default rho_prop.in).
    SELECT CASE ( IARGC() )
        CASE ( 1 )
            CALL GETARG(1, in_file)
        CASE DEFAULT
            in_file = 'rho_prop.in'
    END SELECT

    CALL get_params

    SELECTCASE ( field )

    CASE ( 'cosfield' )
        DO j = 1, npts_per_proc
            IF ( j == npts_per_proc .AND. proc_id < remainder ) THEN
                omega = REAL(npts_per_proc*nprocs + proc_id)*omega_step + omega_from
                WRITE(suffix,'(I8.8)') npts_per_proc*nprocs + proc_id + 1
            ELSE
                omega = REAL(proc_id*npts_per_proc + j - 1)*omega_step + omega_from
                WRITE(suffix,'(I8.8)') proc_id*npts_per_proc + j
            ENDIF

            out_file    = TRIM(out_folder)//'/'//TRIM(jname)//'_'//suffix//'.out'
            params_file = TRIM(out_folder)//'/'//TRIM(jname)//'_'//suffix//'.field_params'

            CALL runge
            CALL print_field_params
        ENDDO

    CASE DEFAULT
        IF ( proc_id == 0 ) THEN
            out_file = TRIM(out_folder)//'/'//TRIM(jname)//'.out'
            params_file = TRIM(out_folder)//'/'//TRIM(jname)//'.params'

            CALL runge
            CALL print_field_params
        ENDIF

    END SELECT

END PROGRAM rho_prop
