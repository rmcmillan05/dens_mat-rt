MODULE print_mat_mod
    USE double
    IMPLICIT NONE

CONTAINS

SUBROUTINE print_str(str_in, out_id_in)
    USE double
    USE global_params , ONLY : full_wd, char_fmt, std_out, std_err
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: str_in
    INTEGER, INTENT(IN), OPTIONAL :: out_id_in

    CHARACTER(LEN=full_wd) :: str_out
    INTEGER :: out_id

    str_out = TRIM(ADJUSTL(str_in))

    IF ( PRESENT(out_id_in) ) THEN
        out_id = out_id_in
    ELSEIF ( str_out(1:5) == 'Error' .OR. str_out(1:7) == 'Warning') THEN
        out_id = std_err
    ELSE
        out_id = std_out
    ENDIF

    WRITE(out_id,char_fmt) str_out

END SUBROUTINE print_str

SUBROUTINE print_str_str(str1_in, str2_in, out_id)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: str1_in
    CHARACTER(LEN=*), INTENT(IN) :: str2_in
    INTEGER, INTENT(IN) :: out_id

    CHARACTER(LEN=29) :: str1_out

    str1_out = TRIM(ADJUSTL(str1_in))

    WRITE(out_id, '(A29)', ADVANCE='NO') str1_out
    WRITE(out_id, '(A3)', ADVANCE='NO') ' : '
    WRITE(out_id, '(A48)') TRIM(str2_in)

END SUBROUTINE print_str_str

SUBROUTINE print_str_num_real(str_in, num_in, out_id)
    USE double
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: str_in
    CHARACTER(LEN=29) :: str_out
    REAL(KIND=DP), INTENT(IN) :: num_in
    INTEGER, INTENT(IN) :: out_id

    str_out = TRIM(ADJUSTL(str_in))

    WRITE(out_id, '(A29)', ADVANCE='NO') str_out 
    WRITE(out_id, '(A3)', ADVANCE='NO') ' : '
    WRITE(out_id, '(A26)', ADVANCE='NO') ''
    WRITE(out_id, '(ES22.14)') num_in

END SUBROUTINE print_str_num_real

SUBROUTINE print_str_num_complex(str_in, num_in, out_id)
    USE double
    IMPLICIT NONE
    
    CHARACTER(LEN=*), INTENT(IN) :: str_in
    CHARACTER(LEN=29) :: str_out
    COMPLEX(KIND=DP), INTENT(IN) :: num_in
    INTEGER, INTENT(IN) :: out_id

    str_out = TRIM(ADJUSTL(str_in))

    WRITE(out_id, '(A29)', ADVANCE='NO') str_out
    WRITE(out_id, '(A3)', ADVANCE='NO') ' : '
    WRITE(out_id, '(A1)', ADVANCE='NO') '('
    WRITE(out_id, '(ES22.14)', ADVANCE='NO') REAL(num_in)
    WRITE(out_id, '(A2)', ADVANCE='NO') ', '
    WRITE(out_id, '(ES22.14)', ADVANCE='NO') AIMAG(num_in)
    WRITE(out_id, '(A1)') ')'

END SUBROUTINE print_str_num_complex

SUBROUTINE print_mat_complex(A, fid, fspec)
    USE double
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: fid_o
    INTEGER :: n, m, i, j

    fid_o = 6
    fspec_o = '(ES22.14)'
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = SIZE(A, 1)
    m = SIZE(A, 2)

    DO i = 1,n
        DO j = 1,m
            WRITE(fid_o, FMT='(A1)' , ADVANCE='NO') '('
            WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') REAL(A(i,j))
            WRITE(fid_o, FMT='(A2)' , ADVANCE='NO') ', '
            WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') AIMAG(A(i,j))
            WRITE(fid_o, FMT='(A3)' , ADVANCE='NO') ')  '
        ENDDO
        WRITE(fid_o,*)
    ENDDO

END SUBROUTINE print_mat_complex

SUBROUTINE print_mat_real(A, fid, fspec)
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: A(:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    INTEGER :: fid_o
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: n, m, i, j

    fspec_o = '(ES22.14)'
    fid_o = 6
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = UBOUND(A, 1)
    m = UBOUND(A, 2)

    DO i = 1,n
        DO j = 1,m
                WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') A(i,j)
                WRITE(fid_o, FMT='(A2)', ADVANCE='NO') ''
        ENDDO
        WRITE(fid,*)
    ENDDO

END SUBROUTINE print_mat_real

SUBROUTINE print_vec_complex(A, fid, fspec)
    USE double
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN) :: A(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: fid_o
    INTEGER :: n, i

    fid_o = 6
    fspec_o = '(ES22.14)'
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = SIZE(A)

    DO i = 1,n
        WRITE(fid_o, FMT='(A1)' , ADVANCE='NO') '('
        WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') REAL(A(i))
        WRITE(fid_o, FMT='(A2)' , ADVANCE='NO') ', '
        WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') AIMAG(A(i))
        WRITE(fid_o, FMT='(A3)' , ADVANCE='NO') ')  '
        WRITE(fid_o,*)
    ENDDO

END SUBROUTINE print_vec_complex

SUBROUTINE print_vec_real(A, fid, fspec)
    USE double
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: A(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    INTEGER :: fid_o
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: n, i

    fspec_o = '(ES22.14)'
    fid_o = 6
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = SIZE(A)

    DO i = 1,n
        WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') A(i)
        WRITE(fid_o,*)
    ENDDO

END SUBROUTINE print_vec_real

SUBROUTINE print_table(A, headers, fid_in, F1_in, F2_in, ftype_in)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: A(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: headers(:)
    INTEGER, OPTIONAL, INTENT(IN) :: F1_in, F2_in
    INTEGER, OPTIONAL, INTENT(IN) :: fid_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ftype_in

    INTEGER :: F1, F2
    INTEGER :: fid
    CHARACTER(LEN=512) :: F1_char, F2_char
    CHARACTER(LEN=512) :: fspec_real
    CHARACTER(LEN=512) :: fspec_char
    CHARACTER(LEN=512) :: ftype

    INTEGER :: i, j, numcol, numrow

    F1 = 22
    F2 = 14
    fid = 6
    ftype = 'ES'

    IF ( PRESENT(F1_in) ) F1 = F1_in
    IF ( PRESENT(F2_in) ) F2 = F2_in
    IF ( PRESENT(fid_in) ) fid = fid_in
    IF ( PRESENT(ftype_in) ) ftype = ftype_in

    WRITE(F1_char, '(I2)') F1
    WRITE(F2_char, '(I2)') F2
    fspec_real = '('//TRIM(ftype)//TRIM(F1_char)//'.'//TRIM(F2_char)//')'
    fspec_char = '(A'//TRIM(F1_char)//')'

    numcol = SIZE(A,2)
    numrow = SIZE(A,1)

    DO i = 1,numcol
        WRITE(fid,fspec_char,ADVANCE='NO') ADJUSTL(headers(i))
        WRITE(fid,'(A2)', ADVANCE='NO') ''
    ENDDO
    WRITE(fid,*)

    DO i = 1,numrow
        DO j = 1,numcol
            WRITE(fid,fspec_real,ADVANCE='NO') A(i,j)
            WRITE(fid,'(A2)', ADVANCE='NO') ''
        ENDDO
        WRITE(fid,*)
    ENDDO

END SUBROUTINE print_table

SUBROUTINE print_break(fid)
    USE global_params , ONLY : tmp_msg
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: fid

    tmp_msg = '--------------------------------------------------------------&
              &------------------'

    CALL print_str(tmp_msg, fid)

END SUBROUTINE print_break

END MODULE print_mat_mod
