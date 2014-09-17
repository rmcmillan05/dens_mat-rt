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
    ! Convert laser energy in eV to wavelength in nm                           
    REAL(KIND=DP), PARAMETER    :: energy_par = 1.239841E3_DP                     
    ! A.u. conversion for time (fs)                                           
    REAL(KIND=DP), PARAMETER    :: time_par = 2.418884326505E-2_DP               
    ! A.u. conversion for length given in nm                                   
    REAL(KIND=DP), PARAMETER    :: length_par = 0.052917721092_DP

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
