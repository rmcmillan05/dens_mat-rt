MODULE gold_dielectric
    USE double
    IMPLICIT NONE

CONTAINS

FUNCTION au_di_function(omega)
    USE double
    USE params , ONLY : wave_par
    IMPLICIT NONE

    COMPLEX(KIND=DP)            :: au_di_function
    REAL(KIND=DP),INTENT(IN)    :: omega
    REAL(KIND=DP)               :: lambda
    COMPLEX(KIND=DP), PARAMETER :: c_i = (0.0_DP, 1.0_DP)
    REAL(KIND=DP), PARAMETER    :: pi = 4.0_DP * ATAN(1.0_DP)
    REAL(KIND=DP), PARAMETER    :: e_inf = 1.53_DP
    REAL(KIND=DP), PARAMETER    :: lambda_p = 145.0_DP
    REAL(KIND=DP), PARAMETER    :: gamma_p = 17000.0_DP
    REAL(KIND=DP), PARAMETER    :: A_1 = 0.94_DP
    REAL(KIND=DP), PARAMETER    :: phi_1 = -pi/4.0_DP
    REAL(KIND=DP), PARAMETER    :: lambda_1 = 468.0_DP
    REAL(KIND=DP), PARAMETER    :: gamma_1  = 2300.0_DP
    REAL(KIND=DP), PARAMETER    :: A_2 = 1.36_DP
    REAL(KIND=DP), PARAMETER    :: phi_2 = -pi/4.0_DP
    REAL(KIND=DP), PARAMETER    :: lambda_2 = 331.0_DP
    REAL(KIND=DP), PARAMETER    :: gamma_2  = 940.0_DP

    lambda = wave_par/omega

    au_di_function = e_inf - 1.0_DP / (lambda_p**2 *                          &
                               (1.0_DP/(lambda**2) +                          &
                                c_i/(gamma_p*lambda)                          &
                               )                                              &
                            )                                                 &
                 + A_1/lambda_1 * ( EXP(c_i*phi_1) /                          &
                                    (1.0_DP/lambda_1                          &
                                     - 1.0_DP/lambda                          &
                                     - c_i/gamma_1                            &
                                    )                                         &
                                   +EXP(-c_i*phi_1) /                         &
                                    (1.0_DP/lambda_1                          &
                                     + 1.0_DP/lambda                          &
                                     + c_i/gamma_1                            &
                                    )                                         &
                                  )                                           &
                 + A_2/lambda_2 * ( EXP(c_i*phi_2) /                          &
                                    (1.0_DP/lambda_2                          &
                                     - 1.0_DP/lambda                          &
                                     - c_i/gamma_2                            &
                                    )                                         &
                                   +EXP(-c_i*phi_2) /                         &
                                    (1.0_DP/lambda_2                          &
                                     + 1.0_DP/lambda                          &
                                     + c_i/gamma_2                            &
                                    )                                         &
                                  )                                  
END FUNCTION au_di_function

FUNCTION au_gamma(omega)
    USE double
    USE sqd_mnp_params , ONLY : eps_0
    IMPLICIT NONE

    COMPLEX(KIND=DP)          :: au_gamma
    REAL(KIND=DP), INTENT(IN) :: omega

    au_gamma = (au_di_function(omega) - eps_0)/(2.0_DP*eps_0 +          &
                                                 au_di_function(omega)) 

END FUNCTION au_gamma

END MODULE gold_dielectric
