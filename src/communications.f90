MODULE communications
    USE double
    USE mpi
    IMPLICIT NONE

CONTAINS

SUBROUTINE set_up_mpi
    USE double
    USE mpi
    USE params , ONLY : nprocs, proc_id, proc_name, vers_mpi1, vers_mpi2
    IMPLICIT NONE

    INTEGER :: error
    INTEGER :: namelen

    CALL MPI_Init(error)
    CALL MPI_Comm_size(MPI_COMM_WORLD, nprocs, error)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, proc_id, error)
    CALL MPI_Get_processor_name(proc_name, namelen, error)
    CALL MPI_Finalize(error)

    vers_mpi1 = MPI_VERSION
    vers_mpi2 = MPI_SUBVERSION

END SUBROUTINE set_up_mpi

END MODULE communications
