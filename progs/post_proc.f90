PROGRAM post_proc
    USE double
    USE print_mod
    USE num_lines , ONLY : numlines
    USE freq_analysis , ONLY : get_freq_amp, find_closest, remove_duplicates
    USE global_params , ONLY : std_out, energy_par, wave_par
    USE params , ONLY : get_params_pp, freqs_in, field_in_file,      &
                                  probe_freq, field_out_file, freqs_out_file,  &
                                  read_col, in_file
    IMPLICIT NONE

    REAL(KIND=DP), ALLOCATABLE :: field(:)
    REAL(KIND=DP), ALLOCATABLE :: field_approx(:)
    COMPLEX(KIND=DP), ALLOCATABLE :: exp_coeff(:)
    REAL(KIND=DP), ALLOCATABLE :: t(:)
    REAL(KIND=DP), ALLOCATABLE :: freqs_out(:)
    CHARACTER(LEN=256), ALLOCATABLE :: headers(:)
    REAL(KIND=DP), ALLOCATABLE :: table(:,:)
    REAL(KIND=DP), ALLOCATABLE :: cos_coeff(:)
    INTEGER, ALLOCATABLE :: t_out_ids(:)
    REAL(KIND=DP), ALLOCATABLE :: sin_coeff(:)
    REAL(KIND=DP) :: max_err
    REAL(KIND=DP), ALLOCATABLE :: pos_freqs(:)
    INTEGER :: nf
    INTEGER :: i,j,nt
    INTEGER :: fid=20
    INTEGER :: probe_freq_id
    REAL(KIND=DP) :: tmp

    ! Take first input argument as the file from which variables are read 
    ! (default pp.in).
    SELECT CASE ( IARGC() )
        CASE ( 1 )
            CALL GETARG(1, in_file)
        CASE DEFAULT
            in_file = 'pp.in'
    END SELECT

    CALL get_params_pp

    nt = numlines(field_in_file)-1

    ALLOCATE( t(nt) )
    ALLOCATE( field(nt) )
    ALLOCATE( field_approx(nt) )

    OPEN(UNIT=10, FILE=field_in_file, STATUS='OLD', ACTION='READ')
        READ(10, *)
        DO i = 1, nt
! WARNING: INPUT TABLE *MUST* BE FORMATTED IN ES22.14 FMT
            READ(10, '(ES22.14)', ADVANCE='NO') t(i)
            DO j = 2,read_col
                IF ( j == read_col ) THEN
                    READ(10, *) field(i)
                ELSE
                    READ(10, '(ES22.14)', ADVANCE='NO') tmp
                ENDIF
            ENDDO
        ENDDO
    CLOSE(10)

    CALL get_freq_amp(field,                                                   &
                      t,                                                       &
                      freqs_in,                                                &
                      freqs_out,                                               &
                      exp_coeff,                                               &
                      cos_coeff,                                               &
                      sin_coeff,                                               &
                        t_out_ids,                                             &
                      field_approx,                                            &
                      max_err)
    OPEN(fid, FILE=freqs_out_file, STATUS='REPLACE')

    ALLOCATE(headers(3))
    ALLOCATE(table(SIZE(exp_coeff),3))
        headers(1) = 'freqs'
        table(:,1) = freqs_out

        headers(2) = 'exp. coeff. (real pt)'
        table(:,2) = REAL(exp_coeff)

        headers(3) = 'exp. coeff. (imag pt)'
        table(:,3) = AIMAG(exp_coeff)
        CALL print_str('## These are the actual coefficients of the &
                       &exponentials', fid)
        CALL print_table(table, headers, fid)
        WRITE(fid,*)
    DEALLOCATE( headers )
    DEALLOCATE( table )

    nf = SIZE(sin_coeff)
    ALLOCATE( pos_freqs(nf) )
    pos_freqs = freqs_out(nf:2*nf-1)

    ALLOCATE(headers(3))
    ALLOCATE(table(nf,3))
        headers(1) = 'freqs'
        table(:,1) = pos_freqs

        headers(2) = 'real amp. (cosine)'
        table(:,2) = cos_coeff

        headers(3) = 'imaginary amp. (sine)'
        table(:,3) = sin_coeff

        CALL print_str('## These are the coefficients of the cosine/sine &
                       &parts (divide by the field to get chi)', fid)
        CALL print_table(table, headers, fid)
        WRITE(fid,*)
    DEALLOCATE( headers )
    DEALLOCATE( table )

    CALL find_closest(pos_freqs, probe_freq, tmp, probe_freq_id)
    CALL print_str('## Selecting "freq_probe" for comparison', fid)
    WRITE(fid, '(A2)', ADVANCE='NO') '! ' 
    WRITE(std_out, '(ES22.14)', ADVANCE='NO') pos_freqs(probe_freq_id)
    WRITE(fid, '(ES22.14)', ADVANCE='NO') pos_freqs(probe_freq_id)
    ! WRITE FREQUENCY ALSO IN eV
    WRITE(std_out, '(ES22.14)', ADVANCE='NO') pos_freqs(probe_freq_id)*        &
                                              energy_par/wave_par
    WRITE(fid, '(ES22.14)', ADVANCE='NO') pos_freqs(probe_freq_id)*            &
                                          energy_par/wave_par
    ! WRITING LASER FREQ - THIS NEEDS TO BE REMOVED
    WRITE(std_out, '(ES22.14)', ADVANCE='NO') pos_freqs(2)*        &
                                              energy_par/wave_par
    WRITE(fid, '(ES22.14)', ADVANCE='NO') pos_freqs(2)*            &
                                          energy_par/wave_par
    WRITE(std_out, '(ES22.14)', ADVANCE='NO') cos_coeff(probe_freq_id)
    WRITE(fid, '(ES22.14)', ADVANCE='NO') cos_coeff(probe_freq_id)
    WRITE(std_out, '(ES22.14)', ADVANCE='NO') sin_coeff(probe_freq_id)
    WRITE(fid, '(ES22.14)', ADVANCE='NO') sin_coeff(probe_freq_id)
    WRITE(std_out, '(F12.5)') max_err
    WRITE(fid, '(F12.5)') max_err
    WRITE(fid, *) 

    ALLOCATE(headers(2))
    ALLOCATE(table(SIZE(freqs_out),2))
        headers(1) = 't_out'
        table(:,1) = t(t_out_ids)

        headers(2) = 'field_out'
        table(:,2) = field(t_out_ids)

        CALL print_str('## These are the points used in the matrix inversion', &
                       fid)
        CALL print_table(table, headers, fid)
    DEALLOCATE( headers )
    DEALLOCATE( table )

        WRITE(fid,*)
        CALL print_str_num_real('Estimated relative error (%)', max_err, fid)

    CLOSE(fid)


    ALLOCATE(headers(3))
    ALLOCATE(table(nt,3))
        headers(1) = 't'
        table(:,1) = t

        headers(2) = 'field'
        table(:,2) = field

        headers(3) = 'field_approx'
        table(:,3) = field_approx

    OPEN(fid, FILE=field_out_file, STATUS='REPLACE')
        CALL print_table(table, headers, fid)
    CLOSE(fid)
    DEALLOCATE( headers )
    DEALLOCATE( table )

END PROGRAM post_proc
