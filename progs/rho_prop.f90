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
                        field_param_from, remainder, proc_id, &
                        get_params, out_file, params_file, out_folder, jname, in_file, &
                        field_param_step, field, print_field_params, field_width
    USE print_mod
    USE params , ONLY : check_file
    IMPLICIT NONE

    CHARACTER(LEN=8) :: suffix
    INTEGER :: suff_int
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

        DO j = 1, npts_per_proc + 1

            IF ( j == npts_per_proc + 1 .AND. proc_id >= remainder ) THEN
                EXIT
            ENDIF

            CALL change_parameter(j, suff_int)

            WRITE(suffix,'(I8.8)') suff_int + 1

            out_file    = TRIM(out_folder)//'/'//TRIM(jname)//'_'//suffix//'.out'
            params_file = TRIM(out_folder)//'/'//TRIM(jname)//'_'//suffix//'.field_params'

            CALL runge
            CALL print_field_params
        ENDDO

!    CASE DEFAULT
!        IF ( proc_id == 0 ) THEN
!            out_file = TRIM(out_folder)//'/'//TRIM(jname)//'.out'
!            params_file = TRIM(out_folder)//'/'//TRIM(jname)//'.params'
!
!            CALL runge
!            CALL print_field_params
!        ENDIF
!
!    END SELECT

CONTAINS

SUBROUTINE change_parameter(i, suffix_int)
    USE double
    USE mpi
    USE communications
    USE params , ONLY : field_change_param, field_height, field_width, &
                        field_param_step, field_param_from, omega
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(OUT) :: suffix_int
    REAL(KIND=DP) :: out_val

    IF ( i == npts_per_proc + 1 ) THEN
        suffix_int = npts_per_proc*nprocs + proc_id
        out_val = REAL(suffix_int)*field_param_step + field_param_from 
    ELSE
        suffix_int = proc_id*npts_per_proc + i - 1
        out_val = REAL(suffix_int)*field_param_step + field_param_from
    ENDIF

    SELECTCASE ( field_change_param )
        CASE ( 'field_height' )
            field_height  = out_val
        CASE ( 'field_width' )
            field_width  = out_val
        CASE ( 'omega' )
            omega        = out_val
    END SELECT

END SUBROUTINE change_parameter

END PROGRAM rho_prop
