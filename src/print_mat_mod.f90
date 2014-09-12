MODULE print_mat_mod
    USE double
    IMPLICIT NONE

CONTAINS

SUBROUTINE print_mat_complex(A, fspec, fid)
    USE double
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: fid_o
    INTEGER :: n, m, i, j

    fid_o = 6
    fspec_o = '(F16.5)'
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

    fspec_o = '(F16.5)'
    fid_o = 6
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = UBOUND(A, 1)
    m = UBOUND(A, 2)

    DO i = 1,n
        DO j = 1,m
                WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') A(i,j)
        ENDDO
        WRITE(fid,*)
    ENDDO

END SUBROUTINE print_mat_real

SUBROUTINE print_vec_complex(A, fspec, fid)
    USE double
    IMPLICIT NONE
    COMPLEX(KIND=DP), INTENT(IN) :: A(:)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
    INTEGER, INTENT(IN), OPTIONAL :: fid
    CHARACTER(LEN=26) :: fspec_o
    INTEGER :: fid_o
    INTEGER :: n, i

    fid_o = 6
    fspec_o = '(F16.5)'
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

    fspec_o = '(F16.5)'
    fid_o = 6
    IF (PRESENT(fspec)) fspec_o = fspec
    IF (PRESENT(fid)) fid_o = fid

    n = SIZE(A)

    DO i = 1,n
        WRITE(fid_o, FMT=fspec_o, ADVANCE='NO') A(i)
        WRITE(fid_o,*)
    ENDDO

END SUBROUTINE print_vec_real

END MODULE print_mat_mod
