MODULE print_mat_mod
    USE double
    IMPLICIT NONE
    PRIVATE
    PUBLIC print_mat

CONTAINS
    SUBROUTINE print_mat(A, fid, fspec)
        USE double
        IMPLICIT NONE
        COMPLEX(KIND=DP), INTENT(IN) :: A(:,:)
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: fspec
        INTEGER, INTENT(IN), OPTIONAL :: fid
        CHARACTER(LEN=26) :: fspec_o
        INTEGER :: n, m, i, j

        fspec_o = '(A1, F16.5, A2, F16.5, A3)'
        IF (PRESENT(fspec)) fspec_o = fspec


        n = UBOUND(A, 1)
        m = UBOUND(A, 2)

        DO i = 1,n
            DO j = 1,m
                IF (PRESENT(fid)) THEN
                    WRITE(fid, FMT=fspec_o, ADVANCE='NO') &
                    '(', REAL(A(i,j),KIND=DP), ', ', IMAG(A(i,j)), ')  '
                ELSE
                    WRITE(*, FMT=fspec_o, ADVANCE='NO') &
                    '(', REAL(A(i,j),KIND=DP), ', ', IMAG(A(i,j)), ')  '
                ENDIF
            ENDDO
            IF (PRESENT(fid)) THEN
                WRITE(fid,*)
            ELSE
                WRITE(*,*)
            ENDIF
        ENDDO

    END SUBROUTINE print_mat
END MODULE print_mat_mod

