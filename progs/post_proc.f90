PROGRAM post_proc
    USE double
    USE print_mat_mod
    USE num_lines , ONLY : numlines
    USE freq_analysis , ONLY : get_freq_amp
    USE print_output_mod
!    USE global_params , ONLY : std_err
    USE post_proc_params , ONLY : get_params_pp, freqs_in, field_in_file
    IMPLICIT NONE

    REAL(KIND=DP), ALLOCATABLE :: field(:)
    REAL(KIND=DP), ALLOCATABLE :: field_approx(:)
    REAL(KIND=DP), ALLOCATABLE :: t(:)
    REAL(KIND=DP), ALLOCATABLE :: freqs_out(:)
    REAL(KIND=DP), ALLOCATABLE :: t_out_mat(:,:)
    CHARACTER(LEN=256), ALLOCATABLE :: t_out_headers(:)
    REAL(KIND=DP), ALLOCATABLE :: coeff_output_mat(:,:)
    CHARACTER(LEN=256), ALLOCATABLE :: coeff_headers(:)
    REAL(KIND=DP), ALLOCATABLE :: output_mat(:,:)
    CHARACTER(LEN=256), ALLOCATABLE :: headers(:)
    REAL(KIND=DP), ALLOCATABLE :: cos_coeff(:)
    INTEGER, ALLOCATABLE :: t_out_ids(:)
    REAL(KIND=DP), ALLOCATABLE :: sin_coeff(:)
    REAL(KIND=DP) :: max_err
    REAL(KIND=DP), ALLOCATABLE :: pos_freqs(:)
    INTEGER :: nf
    INTEGER :: i,nt

    CALL get_params_pp

    nt = numlines(field_in_file)-1

    ALLOCATE( t(nt) )
    ALLOCATE( field(nt) )
    ALLOCATE( field_approx(nt) )

    OPEN(UNIT=10, FILE=field_in_file, STATUS='OLD', ACTION='READ')
        READ(10, *)
        DO i = 1, nt
            READ(10, *) t(i), field(i)
        ENDDO
    CLOSE(10)


    CALL get_freq_amp(field,                                               &
                      t,                                                   &
                      freqs_in,                                               &
                      freqs_out,                                              &
                      cos_coeff,                                            &
                      sin_coeff,                                            &
                        t_out_ids,                                                &
                      field_approx,                                           &
                      max_err)

    nf = SIZE(sin_coeff)
    ALLOCATE( pos_freqs(nf) )
    pos_freqs = freqs_out(nf:2*nf-1)

    ALLOCATE(t_out_headers(2))
    ALLOCATE(t_out_mat(SIZE(freqs_out),2))

    t_out_headers(1) = 't_out'
    t_out_headers(2) = 'field_out'

    t_out_mat(:,1) = t(t_out_ids)
    t_out_mat(:,2) = field(t_out_ids)

    CALL print_output(t_out_mat, t_out_headers)

    ALLOCATE(coeff_headers(3))
    ALLOCATE(coeff_output_mat(nf,3))

    coeff_headers(1) = 'freqs'
    coeff_headers(2) = 'cosine coeff'
    coeff_headers(3) = 'sine coeff'

    coeff_output_mat(:,1) = pos_freqs
    coeff_output_mat(:,2) = cos_coeff
    coeff_output_mat(:,3) = sin_coeff

    CALL print_output(coeff_output_mat, coeff_headers, 18, 8, 'F')

    ALLOCATE(output_mat(nt,3))
    ALLOCATE(headers(3))

    headers(1) = 't'
    headers(2) = 'field'
    headers(3) = 'field_approx'

    output_mat(:,1) = t
    output_mat(:,2) = field
    output_mat(:,3) = field_approx

    CALL print_output(output_mat, headers, 18, 8, 'ES')

END PROGRAM post_proc
