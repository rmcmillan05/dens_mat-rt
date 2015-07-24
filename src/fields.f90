MODULE fields
    USE double
    IMPLICIT NONE

CONTAINS

FUNCTION sinfield(t)
    USE double
    USE params , ONLY : E0, omega, pulse_stop
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: sinfield

    IF ( t <= pulse_stop ) THEN
        sinfield = E0 * SIN(omega * t)
    ELSE
        sinfield = 0.0_DP
    ENDIF

END FUNCTION sinfield

FUNCTION pump_probe(t)
    USE double
    USE params , ONLY : E0, omega, E0_control, omega_control

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP) :: pump_probe

    pump_probe = E0*COS(omega*t) + E0_control*COS(omega_control*t)

END FUNCTION pump_probe

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
    USE params, ONLY : field_width, field_centre, pulse_area, mu, eps_eff1, field_height
    USE global_params , ONLY : pi
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: gauss

!    gauss = field_height * EXP(-2.772589_DP*(t-field_centre)**2/field_width**2)
    gauss = eps_eff1*0.939437278699651_DP*pulse_area/REAL(mu(1,2),KIND=DP)/field_width *&
            EXP(-2.772588722239781_DP*((t-field_centre)/field_width)**2)
!    gauss = 0.939437278699651_DP*pulse_area/REAL(mu(1,2),KIND=DP)/field_width *&
!            EXP(-2.772588722239781_DP*((t-field_centre)/field_width)**2)

END FUNCTION gauss

FUNCTION sech_pulse(t)
    USE double
    USE params, ONLY : field_height, field_centre, field_width, omega
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: sech_pulse

    sech_pulse = field_height*COS(omega*t)/COSH((t-field_centre)/field_width)

END FUNCTION sech_pulse

FUNCTION gauss_pulse(t)
    USE double
    USE params, ONLY : omega
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP)             :: gauss_pulse

    gauss_pulse = SIN(omega*t)*gauss(t)

END FUNCTION gauss_pulse


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

FUNCTION e_sqd_rwa(t, mu, rho21)
    USE double
    USE gold_dielectric , ONLY : au_gamma
    USE global_params ,   ONLY : ci
    USE params ,          ONLY : omega, E0
    USE params ,          ONLY : s_alpha, eps_eff1, dist, rad
    USE params , ONLY : field_height, field_centre, field_width
    IMPLICIT NONE

    REAL(KIND=DP), INTENT(IN) :: t, mu
    COMPLEX(KIND=DP), INTENT(IN) :: rho21
    REAL(KIND=DP)             :: e_sqd_rwa
    REAL(KIND=DP)             :: p_mnp
    COMPLEX(KIND=DP)          :: gma
    COMPLEX(KIND=DP)          :: omega_eff, G

    gma   = au_gamma(omega)
!    p_mnp = 0.5_DP*E0*rad**3 * ( REAL(gma)*COS(omega*t)                        &
!                                +AIMAG(gma)*SIN(omega*t) )
!    e_sqd_rwa = E0*COS(omega*t) + s_alpha*p_mnp/eps_eff1/dist**3

    omega_eff = 0.5_DP*field_height*mu/eps_eff1 * (1.0_DP + s_alpha*rad**3*gma/dist**3)
    omega_eff = omega_eff /COSH((t-field_centre)/field_width)
    G = s_alpha**2*mu**2*rad**3*gma/eps_eff1**2/dist**6

    e_sqd_rwa = 1.0_DP/mu*2.0_DP*REAL((omega_eff + G*rho21*EXP(ci*omega*t))*EXP(-ci*omega*t))

END FUNCTION e_sqd_rwa

FUNCTION e_sqd_chi_const(t, psqd)
    USE double
    USE params , ONLY : field, rad, s_alpha, dist, eps_eff1, eps_eff2, omega, field_centre, field_width, field_height
    USE global_params , ONLY : ci
    USE gold_dielectric , ONLY : au_gamma

    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP), INTENT(IN) :: psqd
    REAL(KIND=DP) :: e_sqd_chi_const

    REAL(KIND=DP) :: f, pmnp
    COMPLEX(KIND=DP) :: gma
!    COMPLEX(KIND=DP) :: omega_eff
!    COMPLEX(KIND=DP) :: G


    gma = au_gamma(omega)

!    omega_eff = (1.0_DP + rad**3*s_alpha*gma/dist**3)/eps_eff1
!    G = s_alpha**2*rad**3*gma/eps_eff1/dist**6
!
!    e_sqd_chi_const = omega_eff*efield(field, t) + G*psqd

    f = 1.0_DP/(COSH((t-field_centre)/field_width))
    pmnp = 2.0_DP*REAL(rad**3*gma*(0.5_DP*field_height*f + s_alpha*psqd/dist**3)*EXP(ci*omega*t))

    e_sqd_chi_const = (field_height*f*COS(omega*t) + s_alpha*pmnp/dist**3)/eps_eff1

END FUNCTION e_sqd_chi_const


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
!        CASE ('e_sqd_rwa')
!            efield = e_sqd_rwa(t)
        CASE ('zero')
            efield = 0.0_DP * t
        CASE ('gauss')
            efield = gauss(t)
        CASE ('gauss_pulse')
            efield = gauss_pulse(t)
        CASE ('sech_pulse')
            efield = sech_pulse(t)
        CASE ('pump_probe')
            efield = pump_probe(t)
        CASE DEFAULT
            CALL print_str("Error: field '"//TRIM(field)//                     &
                           "' not recognized.  Exiting...")
            CALL EXIT(2)
    END SELECT

END FUNCTION efield

END MODULE fields
