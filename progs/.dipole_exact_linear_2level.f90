PROGRAM exact_dipole_2level
    USE double
    USE params
    USE fields
    IMPLICIT NONE
    REAL(KIND=DP), PARAMETER :: wp = omega
    REAL(KIND=DP), PARAMETER :: w21 = 0.8_DP
    REAL(KIND=DP), PARAMETER :: g21 = 0.0_DP
    COMPLEX(KIND=DP), PARAMETER :: mu21 = (0.005_DP,0.0_DP)
    INTEGER, PARAMETER :: N=1024
    REAL(KIND=DP) :: tmax = 200.0_DP
    REAL(KIND=DP) :: t(N)
    REAL(KIND=DP) :: chi1(N)
    REAL(KIND=DP) :: k
    INTEGER :: i

    DO i=1,N
        t(i) = tmax * (i - 1.0_DP) / ( N - 1.0_DP)
    ENDDO

    k = mu21 * CONJG(mu21)

    chi1 = E0 * k * REAL( &
             ci* ( (w21+wp-g21*ci)*EXP(-ci*wp*t)        &
                  -(w21-wp-g21*ci)*EXP( ci*wp*t)        &
                  -2.0_DP*wp*EXP(-(ci*w21+g21)*t) ) / &
           ( (w21+wp-g21*ci)*(w21-wp-g21*ci) )          &
                        , KIND=DP)

    WRITE(*,*) 'wp = ', wp, 'w21 = ', w21, 'mu21 = ', mu21, 'exact = ', &
               2.0_DP*k*w21/( (w21+wp+ci*g21)*(w21-wp-ci*g21)) 

STOP

    WRITE(*,'(3A16)') 'tvec', 'field', 'x'
    DO i=1,N
        WRITE(*,'(4F16.12)') t(i), efield('sinfield',t(i)), chi1(i)
    ENDDO


END PROGRAM exact_dipole_2level

