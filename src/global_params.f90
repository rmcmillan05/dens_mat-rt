MODULE global_params
    USE double
    IMPLICIT NONE
    
    COMPLEX(KIND=DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
    REAL(KIND=DP), PARAMETER    :: pi = 3.1415926535897932
    REAL(KIND=DP), PARAMETER    :: wave_par = 45.5633526_DP
    REAL(KIND=DP), PARAMETER    :: intens_par = 3.50944758E16_DP
    REAL(KIND=DP), PARAMETER    :: au_to_ev = 27.211_DP
    REAL(KIND=DP), PARAMETER    :: energy_par = 1.239841E3_DP
    REAL(KIND=DP), PARAMETER    :: time_par = 2.418884326505E-2_DP
    REAL(KIND=DP), PARAMETER    :: length_par = 0.052917721092_DP
    REAL(KIND=DP), PARAMETER    :: mass_par = 9.10938291E-13_DP
    REAL(KIND=DP), PARAMETER    :: power_par = 0.18023781169_DP
    REAL(KIND=DP), PARAMETER    :: nm_par = 18.897261245650668_DP

    INTEGER, PARAMETER          :: std_err = 0
    INTEGER, PARAMETER          :: std_out = 6
    INTEGER, PARAMETER          :: full_wd=80
    CHARACTER(LEN=5), PARAMETER :: char_fmt = '(A80)'
    CHARACTER(LEN=64)           :: real_fmt = '(ES22.13E3)'
    CHARACTER(LEN=full_wd)      :: tmp_msg
    
END MODULE global_params
