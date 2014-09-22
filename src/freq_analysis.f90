MODULE freq_analysis
    USE double
    IMPLICIT NONE

CONTAINS

SUBROUTINE get_freq_amp(field_in,                                             &
                        t_in,                                                 &
                        freqs_in,                                             &
                        freqs_out,                                            &
                        exp_coeff,                                            &
                        cos_coeff,                                            &
                        sin_coeff,                                            &
                        t_out_ids,                                                &
                        field_approx,                                         &
                        max_err)
    USE double
    USE print_mat_mod
    USE global_params , ONLY : ci, pi, std_err
    USE post_proc_params , ONLY : start_from, go_to, use_max_freq
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: field_in(:)
    REAL(KIND=DP), INTENT(IN) :: t_in(:)
    REAL(KIND=DP), INTENT(INOUT) :: freqs_in(:)

    COMPLEX(KIND=DP), ALLOCATABLE, INTENT(OUT) :: exp_coeff(:)
    REAL(KIND=DP), INTENT(OUT) :: field_approx(SIZE(t_in))
    REAL(KIND=DP), INTENT(OUT) :: max_err
    COMPLEX(KIND=DP), ALLOCATABLE :: field (:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: t (:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: A(:,:)
    COMPLEX(KIND=DP), ALLOCATABLE :: coeff(:,:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: freqs_out(:)
    REAL(KIND=DP), ALLOCATABLE :: t_out(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: t_out_ids(:)
    REAL(KIND=DP), ALLOCATABLE :: t_tmp(:)
    REAL(KIND=DP) :: dt
    REAL(KIND=DP) :: min_period
    CHARACTER(LEN=20) :: gt_char_o, gt_char_n, max_err_char
    REAL(KIND=DP) :: dt_in
    REAL(KIND=DP), PARAMETER :: max_err_tol = 1.0_DP

    INTEGER :: nf ! Number of frequencies
    INTEGER :: i, j
    INTEGER :: s_id, e_id

    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: cos_coeff(:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: sin_coeff(:)

    CALL add_freqs(freqs_in, freqs_out)
    nf = SIZE(freqs_out)
    ALLOCATE( coeff(nf,1) )
    ALLOCATE( A(nf,nf) )
    ALLOCATE( t(nf,1) )
    ALLOCATE( field(nf,1) )
    ALLOCATE( t_out(nf) )
    ALLOCATE( t_out_ids(nf) )

    min_period= 2.0_DP*pi/MAXVAL(freqs_out)
    dt_in = t_in(2) - t_in(1)

    IF ( use_max_freq ) THEN
        go_to = start_from + 0.95_DP*min_period
    ENDIF

    CALL linspace(start_from, go_to, nf, t_tmp, dt)

    IF ( dt < dt_in ) THEN
        WRITE(gt_char_o, '(F16.12)') go_to
        CALL print_str('Warning: dt is too large for optimal frequency &
                       &resolution. Increasing "go_to"...')
        DO WHILE ( dt < dt_in ) 
            go_to = go_to + min_period
            CALL linspace(start_from, go_to, nf, t_tmp, dt)
        ENDDO
        WRITE(gt_char_n, '(F16.12)') go_to
        CALL print_str('"go_to" changed from '//gt_char_o//&
                       &' to '//TRIM(ADJUSTL(gt_char_n)), std_err)
    ENDIF

    IF ( go_to > t_in(SIZE(t_in)) ) THEN
        CALL print_str('Warning: "go_to" is greater than the maximum value of &
                       &t.')
    ENDIF

    CALL find_closest(t_in, start_from, t_out(1), t_out_ids(1))

    DO i = 2, nf
        CALL find_closest(t_in, t_tmp(i), t_out(i), t_out_ids(i))
        IF ( t_out_ids(i) == t_out_ids(i-1) ) THEN
            CALL print_str('Error: singular matrix inversion. Exiting...')
            CALL EXIT(1)
        ENDIF
    ENDDO

    DO i = 1,nf
        DO j = 1,nf
            A(i,j) = EXP(ci*freqs_out(j)*t_out(i))
        ENDDO
    ENDDO

    t(:,1) = t_out
    field(:,1) = field_in(t_out_ids)

    CALL lin_sys(A, field, coeff)

    field_approx = 0.0_DP

    DO i = 1,SIZE(t_in)
        DO j = 1,SIZE(freqs_out)
            field_approx(i) = field_approx(i)                                 &
                              + REAL(coeff(j,1)*EXP(ci*freqs_out(j)*t_in(i)))
        ENDDO
    ENDDO

    ALLOCATE( exp_coeff(SIZE(coeff)) )
    exp_coeff = coeff(:,1)

    s_id = t_out_ids(1)
    e_id = SIZE(t_in)

    max_err = MAXVAL(ABS(field_approx(s_id:e_id) - field_in(s_id:e_id))) /                          &
              (MAXVAL(ABS(field_in(s_id:e_id)))                               &
              -MINVAL(ABS(field_in(s_id:e_id))))
    max_err = max_err*100.0_DP

    IF ( max_err >= max_err_tol ) THEN
        WRITE(max_err_char, '(F8.2)') max_err
        CALL print_str('Warning: maximum error very large (~'//               &
                       TRIM(ADJUSTL(max_err_char))//'%)')
    ENDIF

    nf = (nf-1)/2+1

    ALLOCATE( cos_coeff(nf) )
    ALLOCATE( sin_coeff(nf) )

    cos_coeff(1) = REAL(coeff(nf,1))
    sin_coeff(1) = 0.0_DP

    DO i = 1,nf-1
        cos_coeff(i+1) = REAL(coeff(i,1) + coeff(i+nf,1))
        sin_coeff(i+1) = REAL(ci * ( coeff(i+nf,1) - coeff(i,1) ))
    ENDDO

END SUBROUTINE get_freq_amp

SUBROUTINE add_freqs(freqs_in, freqs_out)
    USE double
    USE post_proc_params , ONLY : c_diffs, max_order
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: freqs_in(:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: freqs_out(:)

    REAL(KIND=DP), ALLOCATABLE :: tmp_freqs(:)
    INTEGER :: nf
    INTEGER :: i, j, k
    INTEGER :: sort_ok

    EXTERNAL DLASRT

    CALL remove_duplicates(freqs_in, freqs_out)
    nf = SIZE(freqs_out)

    ALLOCATE( tmp_freqs(max_order*nf) )

    k = 0
    DO i = 1, max_order
        DO j = 1, nf
            k = k+1
            tmp_freqs(k) = REAL(i,KIND=DP)*freqs_out(j)
        ENDDO
    ENDDO

    DEALLOCATE( freqs_out )
    IF ( c_diffs ) THEN
        CALL calc_diffs(tmp_freqs, freqs_out)
    ELSE 
        ALLOCATE(freqs_out(SIZE(tmp_freqs)))
        freqs_out = tmp_freqs
    ENDIF

    DEALLOCATE( tmp_freqs )
    CALL remove_duplicates(freqs_out, tmp_freqs)

    CALL DLASRT('I', SIZE(tmp_freqs), tmp_freqs,sort_ok)

    DEALLOCATE( freqs_out )
    CALL add_neg(tmp_freqs, freqs_out)

END SUBROUTINE add_freqs

SUBROUTINE remove_duplicates(a_in, a_out)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: a_in(:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: a_out(:)
    REAL(KIND=DP) :: a_tmp(SIZE(a_in))
    REAL(KIND=DP) :: tol = 1.0E-12
    INTEGER :: k
    INTEGER :: i, j

    a_tmp = 0.0_DP
    k = 1

outer: DO i = 1, SIZE(a_in)
        DO j = 1,k
            IF ( a_in(i) <= tol ) THEN ! Also remove 0
                CYCLE outer
            ELSEIF ( ABS(a_in(i)-a_tmp(j)) <= tol ) THEN
                CYCLE outer
            ENDIF
        ENDDO
        a_tmp(k) = a_in(i)
        k = k + 1
    ENDDO outer
    
    k = k-1

    ALLOCATE( a_out(k) )
    a_out = a_tmp(1:k)

END SUBROUTINE remove_duplicates

SUBROUTINE calc_diffs(a, diffs)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: a(:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: diffs(:)

    INTEGER :: i, j, k, n

    n = SIZE(a)

    ALLOCATE( diffs(n*(n+1)/2) )

    k = 0
    DO i = 1,n
        DO j = 1,n
            IF ( i == j ) THEN
                k = k+1
                diffs(k) = a(i)
            ELSEIF ( i > j ) THEN
                k = k+1
                diffs(k) = ABS(a(i)-a(j))
            ENDIF
        ENDDO
    ENDDO

END SUBROUTINE calc_diffs

SUBROUTINE add_neg(a_in, a_out)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: a_in(:)
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: a_out(:)
    
    ALLOCATE( a_out(2*SIZE(a_in)+1) )

    a_out(1:SIZE(a_in)) = -1.0_DP * a_in
    a_out(SIZE(a_in)+1) = 0.0_DP
    a_out(SIZE(a_in)+2:2*SIZE(a_in)+1) = a_in

END SUBROUTINE add_neg

SUBROUTINE lin_sys(A, B, X)
    USE double
    USE print_mat_mod , ONLY : print_str
    IMPLICIT NONE

    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
    COMPLEX(KIND=DP), INTENT(IN) :: B(:,:)
    COMPLEX(KIND=DP), INTENT(OUT) :: X(:,:)

    COMPLEX(KIND=DP) :: A_new(SIZE(A,1),SIZE(A,2))
    INTEGER :: n
    INTEGER :: nrhs
    INTEGER :: info
    INTEGER :: ipiv(SIZE(A,1))

    EXTERNAL ZGESV

    n = SIZE(B,1)
    nrhs = SIZE(B,2)
    A_new = A
    X = B

    CALL ZGESV(n, nrhs, A_new, n, ipiv, X, n, info)

    IF ( info /= 0 ) THEN
        CALL print_str('Error: no solution to linear system of equations.')
    ENDIF

END SUBROUTINE lin_sys

SUBROUTINE find_closest(t_in, find_t, found_t, found_t_id)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t_in(:)
    REAL(KIND=DP), INTENT(IN) :: find_t
    REAL(KIND=DP), INTENT(OUT) :: found_t
    INTEGER, INTENT(OUT) :: found_t_id

    INTEGER :: i, n
    REAL(KIND=DP) :: err, err_tmp
    REAL(KIND=DP) :: found_t_tmp
    INTEGER :: found_t_id_tmp

    n = SIZE(t_in)

    err_tmp = 0.0_DP
    found_t = t_in(1)
    found_t_id = 1

    DO i = 2, n+1
        found_t_tmp = t_in(i-1)
        found_t_id_tmp = i-1

        err = ABS(found_t_tmp - find_t)

        IF ( err < err_tmp ) THEN
            found_t = found_t_tmp
            found_t_id = found_t_id_tmp
        ENDIF

        err_tmp = err

    ENDDO

END SUBROUTINE find_closest

SUBROUTINE linspace(a, b, n, t, dt)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: a
    REAL(KIND=DP), INTENT(IN) :: b
    INTEGER, INTENT(IN) :: n
    REAL(KIND=DP), ALLOCATABLE, INTENT(OUT) :: t(:)
    REAL(KIND=DP), INTENT(OUT) :: dt

    INTEGER :: i

    ALLOCATE( t(n) )

    dt = (b-a)/REAL(n-1,KIND=DP)

    t(1) = a
    DO i = 2,n
        t(i) = t(i-1) + dt
    ENDDO

END SUBROUTINE linspace

END MODULE freq_analysis
