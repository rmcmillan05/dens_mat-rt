MODULE fields
    USE double
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: efield

CONTAINS

FUNCTION sinfield(t)
    USE double
    USE params , ONLY : E0, omega_au
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: sinfield

    sinfield = E0 * SIN(omega_au * t)

END FUNCTION sinfield

FUNCTION cosfield(t)
    USE double
    USE params , ONLY : E0, omega_au
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: cosfield

    cosfield = E0 * COS(omega_au * t)

END FUNCTION cosfield

FUNCTION step(t)
    USE double
    USE params, ONLY : step_height, step_centre, step_width
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: step

    IF (t > step_centre .AND. t <= step_centre+step_width) THEN
       step = step_height
    ELSE
       step = 0.0_DP
    ENDIF 

END FUNCTION step

FUNCTION pulse(t_in)
    USE double
    USE params , ONLY : omega_au, E0, pulse_phase, pulse_lim, omega_au,   &
                        pulse_start
    USE global_params , ONLY : pi
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t_in
    REAL(KIND=DP)             :: t
    REAL(KIND=DP)             :: pulse
    REAL(KIND=DP)             :: f, df, w
    
    t = t_in - pulse_start
    IF ( t >= 0.0_DP .AND. t <= pulse_lim ) THEN
        f     = SIN(pi * t / pulse_lim) ** 2
        df    = 2.0_DP * SIN(pi * t / pulse_lim) *                            &
                         COS(pi * t / pulse_lim) * pi / pulse_lim
        w     = omega_au * t + pulse_phase
        pulse = E0 * f * SIN(w) - E0 * df * COS(w) / omega_au
    ELSE
        pulse = 0.0_DP
    ENDIF

END FUNCTION pulse

FUNCTION e_sqd(t)
    USE double
    USE gold_dielectric , ONLY : au_gamma
    USE params , ONLY : omega_au, E0
    USE params , ONLY : s_alpha, eps_eff1, dist, rad
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: e_sqd
    REAL(KIND=DP)             :: p_mnp
    COMPLEX(KIND=DP)          :: gma

    gma = au_gamma(omega_au)
    p_mnp = 0.5_DP*E0*rad**3 * ( REAL(gma)*COS(omega_au*t)                    &
                                +AIMAG(gma)*SIN(omega_au*t) )
    e_sqd = E0*COS(omega_au*t) + s_alpha*p_mnp/eps_eff1/dist**3

END FUNCTION e_sqd

FUNCTION efield(field, t)
    USE double
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(INOUT) :: field
    REAL(KIND=DP), INTENT(IN)       :: t
    REAL(KIND=DP)                   :: efield

    SELECT CASE ( TRIM(field) )
        CASE ('step')
            efield = step(t)
        CASE ('pulse')
            efield = pulse(t)
        CASE ('sinfield')
            efield = sinfield(t)
        CASE ('cosfield')
            efield = cosfield(t)
        CASE ('e_sqd')
            efield = e_sqd(t)
        CASE ('zero')
            efield = 0.0_DP * t
        CASE DEFAULT
            field = 'zero'
            efield = 0.0_DP * t
    END SELECT

END FUNCTION efield

END MODULE fields
