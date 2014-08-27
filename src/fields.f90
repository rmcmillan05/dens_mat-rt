MODULE fields
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: efield

CONTAINS

FUNCTION sinfield(t)
    USE double
    USE params , ONLY : E0, omega_au
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP) :: sinfield

    sinfield = E0 * SIN(omega_au * t)

END FUNCTION sinfield

FUNCTION cosfield(t)
    USE double
    USE params , ONLY : E0, omega_au
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP) :: cosfield

    cosfield = E0 * COS(omega_au * t)
END FUNCTION cosfield

!FUNCTION step(omega,t)
!    USE double
!    IMPLICIT NONE
!    REAL(KIND=DP), INTENT(IN) :: omega
!    REAL(KIND=DP), INTENT(IN) :: t
!    REAL(KIND=DP)             :: step
!    REAL(KIND=DP), PARAMETER  :: c=5.0_DP, w=20.0_DP, h=1.8_DP
!
!!    IF (t > c-w .AND. t < c+w) THEN
!    IF (t > c .AND. t <= c+w) THEN
!       step = h
!    ELSE
!       step = 0.0_DP
!    ENDIF 
!
!END FUNCTION step

FUNCTION pulse(t_in)
    USE double
    USE params , ONLY : omega_au, E0, pi
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: t_in
    REAL(KIND=DP) :: t
    REAL(KIND=DP)             :: lim
    REAL(KIND=DP)             :: pulse
    REAL(KIND=DP)             :: f, df, w
    REAL(KIND=DP) :: phi = 0.0_DP
    REAL(KIND=DP) :: num_cycles = 6.0_DP
    REAL(KIND=DP) :: start = 0.0_DP

    lim = 2.0_DP * pi * num_cycles / omega_au
    
    t = t_in - start
    IF ( t >= 0.0_DP .AND. t <= lim ) THEN
        f     = SIN(pi * t / lim) ** 2
        df    = 2.0_DP * SIN(pi * t / lim) * COS(pi * t / lim) * pi / lim
        w     = omega_au * t + phi
        pulse = E0 * f * SIN(w) - E0 * df * COS(w) / omega_au
    ELSE
        pulse = 0.0_DP
    ENDIF

END FUNCTION pulse

FUNCTION e_sqd(t)
    USE double
    USE gold_dielectric
    USE params , ONLY : omega_au, E0
    USE sqd_mnp_params , ONLY : s_alpha, eps_eff1, dist, rad
    IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN) :: t
    REAL(KIND=DP) :: e_sqd
    REAL(KIND=DP) :: p_mnp
    COMPLEX(KIND=DP) :: gma

    gma = au_gamma(omega_au)
    p_mnp = 0.5_DP*E0*rad**3 * ( REALPART(gma)*COS(omega_au*t)   &
                                +IMAGPART(gma)*SIN(omega_au*t) )
    e_sqd = E0*COS(omega_au*t) + s_alpha*p_mnp/eps_eff1/dist**3

END FUNCTION e_sqd

FUNCTION efield(field, t)
    USE double
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: field
    REAL(KIND=DP), INTENT(IN)    :: t
    REAL(KIND=DP)                :: efield

    SELECT CASE ( TRIM(field) )
!        CASE ('step')
!            efield = step(omega,t)
        CASE ('pulse')
            efield = pulse(t)
        CASE ('sinfield')
            efield = sinfield(t)
        CASE ('cosfield')
            efield = cosfield(t)
        CASE ('e_sqd')
            efield = e_sqd(t)
        CASE DEFAULT
            efield = 0.0_DP * t
    END SELECT

END FUNCTION efield


END MODULE fields
