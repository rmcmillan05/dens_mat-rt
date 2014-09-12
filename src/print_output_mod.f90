MODULE print_output_mod
    USE double
    IMPLICIT NONE

CONTAINS

SUBROUTINE print_output(A, headers, F1_in, F2_in, ftype_in, fid_in, fname_in)
    USE double
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: A(:,:)
    CHARACTER(LEN=*), INTENT(IN) :: headers(:)
    INTEGER, OPTIONAL, INTENT(IN) :: F1_in, F2_in
    INTEGER, OPTIONAL, INTENT(IN) :: fid_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: fname_in
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ftype_in

    INTEGER :: F1, F2
    INTEGER :: fid
    CHARACTER(LEN=512) :: fname
    CHARACTER(LEN=512) :: F1_char, F2_char
    CHARACTER(LEN=512) :: fspec_real
    CHARACTER(LEN=512) :: fspec_char
    CHARACTER(LEN=512) :: ftype

    INTEGER :: i, j, numcol, numrow

    F1 = 18
    F2 = 8
    fid = 6
    fname = 'out'
    ftype = 'ES'

    IF ( PRESENT(F1_in) ) F1 = F1_in
    IF ( PRESENT(F2_in) ) F2 = F2_in
    IF ( PRESENT(fid_in) ) fid = fid_in
    IF ( PRESENT(fname_in) ) fname = fname_in
    IF ( PRESENT(ftype_in) ) ftype = ftype_in

    WRITE(F1_char, '(I2)') F1
    WRITE(F2_char, '(I2)') F2

    fspec_real = '('//TRIM(ftype)//TRIM(F1_char)//'.'//TRIM(F2_char)//')'
    fspec_char = '(A'//TRIM(F1_char)//')'

    numcol = SIZE(A,2)
    numrow = SIZE(A,1)

    DO i = 1,numcol
        WRITE(fid,fspec_char,ADVANCE='NO') ADJUSTL(headers(i))
    ENDDO
    WRITE(fid,*)

    DO i = 1,numrow
        DO j = 1,numcol
            WRITE(fid,fspec_real,ADVANCE='NO') A(i,j)
        ENDDO
        WRITE(fid,*)
    ENDDO

END SUBROUTINE print_output

END MODULE print_output_mod

