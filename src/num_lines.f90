MODULE num_lines
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: numlines

CONTAINS

FUNCTION numlines(fname)
        ! From http://physics.bu.edu/py502/lectures1/examples/read.f90
        IMPLICIT NONE
        INTEGER :: numlines
        CHARACTER :: c
        CHARACTER(LEN=*), INTENT(IN) :: fname

        OPEN(1, FILE=fname, STATUS='OLD', ACTION='READ')
        numlines = 0
        DO
            READ(1, *, END=10)c 
            numlines = numlines + 1
        ENDDO
10      CLOSE(1)

END FUNCTION numlines

END MODULE num_lines
