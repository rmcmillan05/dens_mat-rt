MODULE sqd_mnp_params
    USE double
    USE params , ONLY : length_par
    IMPLICIT NONE

    ! SQD_MNP
    !
    ! Distance of centres between MNP and SQD in nm
    REAL(KIND=DP), PARAMETER :: dist_nm = 13.0_DP 
    ! Radius of MNP in nm
    REAL(KIND=DP), PARAMETER :: rad_nm = 7.5_DP
    ! s_alpha = 2 for z-axis, -1 for x,y (z is axis of molecule)
    REAL(KIND=DP), PARAMETER :: s_alpha = 2.0_DP
    ! Dielectric constant of background medium
    REAL(KIND=DP), PARAMETER :: eps_0 = 1.0_DP
    ! Dielectric constant of SQD
    REAL(KIND=DP), PARAMETER :: eps_s = 6.0_DP

    ! CALCULATIONS AND UNIT CONVERSIONS
    !
    REAL(KIND=DP), PARAMETER :: dist = dist_nm/length_par
    REAL(KIND=DP), PARAMETER :: rad = rad_nm/length_par
    REAL(KIND=DP), PARAMETER :: eps_eff1 = (2.0_DP*eps_0 + eps_s)/(3.0_DP*eps_0)
    REAL(KIND=DP), PARAMETER :: eps_eff2 = (2.0_DP*eps_0 + eps_s)/3.0_DP

END MODULE sqd_mnp_params
