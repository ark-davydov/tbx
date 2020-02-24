program tbx
#ifdef MPI
  use mpi
#endif
use modcom
use parameters
use tasksclass
implicit none
character(len=20) arg
type(CLpars) pars
type(CLtasks) tasks
#ifdef _OPENMP
  integer omp_get_max_threads
#endif

#ifdef MPI
  call mpi_init(mpi_err)
  call mpi_comm_dup(mpi_comm_world,mpi_com,mpi_err)
  call mpi_comm_size(mpi_com,np_mpi,mpi_err)
  call mpi_comm_rank(mpi_com,lp_mpi,mpi_err)
  if (lp_mpi.eq.0) then
    mp_mpi=.true.
    write(*,*)
    if (np_mpi.gt.1) then
      write(*,*)
      write(*,'("Using MPI, number of processes : ",I8)') np_mpi
    end if
  else
    mp_mpi=.false.
  end if
#else
  mp_mpi=.true.
  lp_mpi=0
  np_mpi=1
#endif

#ifdef _OPENMP
  if (mp_mpi) write(*,'("Number of OMP threads at each process : ",I8)') omp_get_max_threads()
#endif

#ifdef PARDI
  call mkl_disable_fast_mm()
#endif

call get_command_argument(1, arg)
call pars%init(trim(adjustl(arg)))
call tasks%run(pars)

call message("Finished")

#ifdef MPI
  call mpi_finalize(mpi_err)
#endif

end
