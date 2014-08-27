PROGRAM rho_perturb
    !
    ! Calculates the rho matrices 
    !           p^(0), p^(1), p^(2), ...
    ! in perturbation theory along with the dipole
    !           d^(i) = trace(p^(i) * mu)
    ! for each order, at each time-step.
    ! The sum of the d^(i) is also calculated.
    !
    ! To call:: rho_perturb.exe input_directory [trange] [electric_field] ...
    !                           [npts] [max_order]
    !
    ! input_directory is the path containing the rho_0, mu, ham, gma matrices.
    !
    ! trange is the range in time over which to approximate, split into npts
    ! evenly spaced points (starting at t=0 up to trange).
    !
    ! electric_field at the moment can be either 'step', 'pulse' or 'sinfield'
    ! see src/fields.f90 to modify or add fields.
    !
    ! max_order is the maximum order in the perturbation expansion to which
    ! the calculation is performed.
    !
    USE double
    USE fields
    USE num_lines
    USE omp_lib
    IMPLICIT NONE

    CHARACTER(LEN=64)  :: field='sinfield'
    CHARACTER(LEN=64)  :: L_in='200'
    CHARACTER(LEN=64)  :: npts_in = '10000'
    CHARACTER(LEN=64)  :: out_fmat='(ES16.8)'
    CHARACTER(LEN=64)  :: mo_in='7'

    INTEGER                                     :: num_levels
    INTEGER                                     :: npts
    INTEGER                                     :: max_order

    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:,:) :: p0, p1
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:,:)   :: rho_0, ham, gma, rho,&
                                                       mu, dipole
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:)     :: dip_current
    COMPLEX(KIND=DP), ALLOCATABLE, DIMENSION(:)     :: dip_all
    REAL(KIND=DP),    ALLOCATABLE, DIMENSION(:)     :: t
    REAL(KIND=DP),    ALLOCATABLE, DIMENSION(:)     :: E

    REAL(KIND=DP)                               :: delta_t
    REAL(KIND=DP)                               :: w_nm
    REAL(KIND=DP)                               :: t_max
    COMPLEX(KIND=DP)                            :: c
    COMPLEX(KIND=DP), PARAMETER                 :: ci = (0.0_DP, 1.0_DP)
    COMPLEX(KIND=DP), PARAMETER                 :: c0 = (0.0_DP, 0.0_DP)
    INTEGER                                     :: i, j, n, m, v

    CHARACTER(LEN=1)  :: dip_num_char
    CHARACTER(LEN=64) :: rho_in, ham_in, gma_in, data_pref, mu_in


    SELECT CASE ( IARGC() )
        CASE ( 5)
            CALL GETARG(1, data_pref)
            CALL GETARG(2, L_in)
            CALL GETARG(3, field)
            CALL GETARG(4, npts_in)
            CALL GETARG(5, mo_in)
        CASE ( 4 )
            CALL GETARG(1, data_pref)
            CALL GETARG(2, L_in)
            CALL GETARG(3, field)
            CALL GETARG(4, npts_in)
        CASE ( 3 )
            CALL GETARG(1, data_pref)
            CALL GETARG(2, L_in)
            CALL GETARG(3, field)
        CASE ( 2 )
            CALL GETARG(1, data_pref)
            CALL GETARG(2, L_in)
        CASE ( 1 )
            CALL GETARG(1, data_pref)
        CASE DEFAULT
            WRITE(*,*) 'Specify directory containing the matrix files:'
            WRITE(*,*) 'rho.txt, ham.txt, mu.txtm and gma.txt'
            CALL EXIT(0)
    END SELECT

    rho_in=TRIM(data_pref)//'rho.txt'
    ham_in=TRIM(data_pref)//'ham.txt'
    gma_in=TRIM(data_pref)//'gma.txt'
    mu_in =TRIM(data_pref)//'mu.txt'
    
    READ(L_in,*) t_max
    READ(npts_in,*) npts
    READ(mo_in,*) max_order

    ALLOCATE(t(npts+1),E(npts+1),dip_current(npts+1))
    ALLOCATE(dip_all(npts+1))
    dip_all = 0.0_DP
    
    num_levels = numlines(rho_in)

    ALLOCATE(rho_0(num_levels,num_levels), ham(num_levels,num_levels),  &
             gma  (num_levels,num_levels), rho(num_levels,num_levels),  &
             mu   (num_levels,num_levels))

    ALLOCATE(p0(num_levels,num_levels,npts+1), p1(num_levels,num_levels,npts+1))

    ALLOCATE(dipole(max_order,npts+1))

    OPEN(UNIT=10, FILE=rho_in, STATUS='OLD')
        DO i = 1, num_levels
            READ(10, *) (rho_0(i, j), j=1,num_levels)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=ham_in, STATUS='OLD')
        DO i = 1, num_levels
            READ(10, *) (ham(i, j), j=1,num_levels)
        ENDDO
    CLOSE(10)

    OPEN(UNIT=10, FILE=gma_in, STATUS='OLD')
        DO i = 1, num_levels
            READ(10, *) (gma(i, j), j=1,num_levels)
        ENDDO
    CLOSE(10)
    
    OPEN(UNIT=10, FILE=mu_in, STATUS='OLD')
        DO i = 1, num_levels
            READ(10, *) (mu(i, j), j=1,num_levels)
        ENDDO
    CLOSE(10)
    
    delta_t = t_max/npts

    DO i=1,npts+1
        p0(:,:,i) = rho_0
        t(i) = delta_t*(i-1.0_DP) 
        E(i) = efield(field,t(i))
    ENDDO

    p1 = c0

    DO j=1,max_order

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(n,m,v,i,w_nm,c) SCHEDULE(DYNAMIC)        
        DO n=1,num_levels
            DO m=1,num_levels
                w_nm = REAL(ham(n,n) - ham(m,m),KIND=DP)
                c = 0.0_DP
                DO v=1,num_levels
                    c = c + mu(n,v)*p0(v,m,1) - p0(n,v,1)*mu(v,m)
                ENDDO
                p1(n,m,1) = c*E(1)*EXP((ci*w_nm + gma(n,m))*t(1))
                DO i=2,npts+1
                    c = 0.0_DP
                    DO v=1,num_levels
                        c = c + mu(n,v)*p0(v,m,i) - p0(n,v,i)*mu(v,m)
                    ENDDO
                    p1(n,m,i) = p1(n,m,i-1) + c*E(i)*EXP((ci*w_nm + gma(n,m))*t(i))
                ENDDO
                p1(n,m,:) = EXP(-(ci*w_nm + gma(n,m))*t) * ci * delta_t * p1(n,m,:)
            ENDDO
        ENDDO
!$OMP END PARALLEL DO        

        dip_current = 0.0_DP
        DO n=1,num_levels
            DO m=1,num_levels
                dip_current = dip_current + p1(n,m,:)*mu(m,n)
            ENDDO
        ENDDO

        p0 = p1

        dipole(j,:) = dip_current
        dip_all = dip_all + dip_current

    ENDDO

    WRITE(*,'(A16)',ADVANCE='NO') ' t              '
    WRITE(*,'(A16)',ADVANCE='NO') field
    DO j=1,max_order
        WRITE(dip_num_char,'(I1)') j
        WRITE(*,'(A16)',ADVANCE='NO') 'dip'//dip_num_char//'_num'
    ENDDO
    WRITE(*,'(A16)',ADVANCE='NO') 'dip_sum'
    DO j=1,max_order
        WRITE(*,'(F12.8)',ADVANCE='NO') SUM(ABS(IMAG(dipole(j,:))))
    ENDDO
    WRITE(*,*)

    DO i=1,npts+1
        WRITE(*,out_fmat,ADVANCE='NO') t(i)
        WRITE(*,out_fmat,ADVANCE='NO') E(i)
        DO j =1,max_order
            WRITE(*,out_fmat,ADVANCE='NO') REAL(dipole(j,i),KIND=DP)
        ENDDO
        WRITE(*,out_fmat,ADVANCE='NO') REAL(dip_all(i),KIND=DP)
        WRITE(*,*)
    ENDDO

END PROGRAM rho_perturb
