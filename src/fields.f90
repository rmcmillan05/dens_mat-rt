MODULE fields
    USE double
    IMPLICIT NONE

CONTAINS

FUNCTION sinfield(t)
    USE double
    USE params , ONLY : E0, omega
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: sinfield

    sinfield = E0 * SIN(omega * t)

END FUNCTION sinfield

FUNCTION cosfield(t)
    USE double
    USE params , ONLY : E0, omega
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: cosfield

    cosfield = E0 * COS(omega * t)

END FUNCTION cosfield

FUNCTION step(t)
    USE double
    USE params, ONLY : field_height, field_centre, field_width
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: step

    IF (t > field_centre .AND. t <= field_centre+field_width) THEN
       step = field_height
    ELSE
       step = 0.0_DP
    ENDIF 

END FUNCTION step

FUNCTION gauss(t)
    USE double
    USE params, ONLY : field_height, field_centre, field_width
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: gauss

    gauss = field_height * EXP(-2.772589_DP*(t-field_centre)**2/field_width**2)

END FUNCTION gauss


FUNCTION pulse(t_in)
    USE double
    USE params ,        ONLY : omega, E0, pulse_phase, pulse_start, pulse_cycles
    USE global_params , ONLY : pi
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t_in
    REAL(KIND=DP)             :: pulse_lim
    REAL(KIND=DP)             :: t
    REAL(KIND=DP)             :: pulse
    REAL(KIND=DP)             :: f, df, w
    
    pulse_lim = 2.0_DP * pi * pulse_cycles / omega
    t = t_in - pulse_start

    IF ( t >= 0.0_DP .AND. t <= pulse_lim ) THEN
        f     = SIN(pi * t / pulse_lim) ** 2
        df    = 2.0_DP * SIN(pi * t / pulse_lim) *                             &
                         COS(pi * t / pulse_lim) * pi / pulse_lim
        w     = omega * t + pulse_phase
        pulse = E0 * f * SIN(w) - E0 * df * COS(w) / omega
    ELSE
        pulse = 0.0_DP
    ENDIF

END FUNCTION pulse

FUNCTION e_sqd_RWA(t)
    USE double
    USE gold_dielectric , ONLY : au_gamma
    USE params ,          ONLY : omega, E0
    USE params ,          ONLY : s_alpha, eps_eff1, dist, rad
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: e_sqd_RWA
    REAL(KIND=DP)             :: p_mnp
    COMPLEX(KIND=DP)          :: gma

    gma   = au_gamma(omega)
    p_mnp = 0.5_DP*E0*rad**3 * ( REAL(gma)*COS(omega*t)                        &
                                +AIMAG(gma)*SIN(omega*t) )
    e_sqd_RWA = E0*COS(omega*t) + s_alpha*p_mnp/eps_eff1/dist**3

END FUNCTION e_sqd_RWA

FUNCTION efield(field, t)
    USE double
    USE print_mod , ONLY : print_str
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
        CASE ('e_sqd_RWA')
            efield = e_sqd_RWA(t)
        CASE ('zero')
            efield = 0.0_DP * t
        CASE ('gauss')
            efield = gauss(t)
        CASE DEFAULT
            CALL print_str("Error: field '"//TRIM(field)//                     &
                           "' not recognized.  Exiting...")
            CALL EXIT(2)
    END SELECT

END FUNCTION efield

END MODULE fields
