MODULE global_params
    USE double
    IMPLICIT NONE
    
    ! UNIVERSAL PARAMETERS
    !
    ! Imaginary number
    COMPLEX(KIND=DP), PARAMETER :: ci = (0.0_DP, 1.0_DP)
    ! pi
    REAL(KIND=DP), PARAMETER    :: pi = 3.1415926535897932
    ! A.u. of wavelength when given in nm
    REAL(KIND=DP), PARAMETER    :: wave_par = 45.5633526_DP                            
    ! A.u. of intensity                                                        
    REAL(KIND=DP), PARAMETER    :: intens_par = 3.50944758E16_DP                  
    ! Multiply by this to convert energy in eV to a.u.
    REAL(KIND=DP), PARAMETER    :: ev_to_au = 0.036749309_DP
    ! Convert laser energy in eV to wavelength in nm                           
    REAL(KIND=DP), PARAMETER    :: energy_par = 1.239841E3_DP                     
    ! A.u. conversion for time (fs)                                           
    REAL(KIND=DP), PARAMETER    :: time_par = 2.418884326505E-2_DP               
    ! A.u. conversion for length given in nm                                   
    REAL(KIND=DP), PARAMETER    :: length_par = 0.052917721092_DP
    ! A.u. conversion for mass given in kg
    REAL(KIND=DP), PARAMETER    :: mass_par = 9.10938291E-13_DP
    ! Multiplu a.u. by watt_par to convert a.u. of power to watts
    REAL(KIND=DP), PARAMETER :: power_par = 0.18023781169_DP
    ! Multiply length in nm by this to get length in a.u.
    REAL(KIND=DP), PARAMETER :: nm_par = 18.897261245650668_DP

    ! FOR TEXT OUTPUT
    ! Standard error number
    INTEGER, PARAMETER :: std_err = 0
    ! Standard output number
    INTEGER, PARAMETER :: std_out = 6
    !
    INTEGER, PARAMETER :: full_wd=80
    !
    CHARACTER(LEN=5), PARAMETER :: char_fmt = '(A80)'
    ! Message string holder
    CHARACTER(LEN=full_wd) :: tmp_msg
    

END MODULE global_params
