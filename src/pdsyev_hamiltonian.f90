program exact_diagonalization

  use m_npy  
  implicit none
  
  ! -- file names and IO
  character(len=1023) :: project_dir
  character(len=1023) :: group_specs_fname
  character(len=1023) :: spec_info_fname
  character(len=1023) :: hamiltonian_fname
  character(len=1023) :: eigen_fname
  character(len=1023) :: fname_prefix
  character(1) :: char_buffer
  integer :: iostatus
  
  ! -- double precision kind
  integer, parameter :: dp = selected_real_kind(15, 307)
  
  ! -- group specs
  integer :: p,q,N
  integer :: hdim
  logical :: periodic_boundary

  ! -- hamiltonian
  real(dp), dimension(:,:), allocatable :: hamiltonian
  integer :: i,j
  real(dp) :: val

  !-- BLACKS
  integer :: IAM, NPROCS
  integer :: NPROW, NPCOL

  !-- ScaLAPACK
  integer :: LDA, LWORK, INFO, MAXN, MAXPROCS
  real(dp) :: MONE

  character :: JOBZ, UPLO 
  real(dp), dimension(:), allocatable :: W, WORK, Z
  integer, dimension(:), allocatable :: DESCA, DESCZ
  real(dp), dimension(1) :: WORK_QUERY
  !-- location of files is relative to this path
  call get_environment_variable("HYPERBOLIC_DIR", project_dir)

  !#####################################################################
  ! Reading of group_specs.inp                                         # 
  !#####################################################################
  
  group_specs_fname = trim(trim(project_dir) // trim("/group_specs.inp"))

  write(*,*) "Loading: ", trim(group_specs_fname)
  open(unit=10, file=trim(group_specs_fname))
  read(10,*,iostat=iostatus) char_buffer, p
  read(10,*,iostat=iostatus) char_buffer, q
  read(10,*,iostat=iostatus) char_buffer, N

  if (iostatus>0) then
    write(*,*) "Error occured while trying to read the group_specs.inp file"
  end if
  
  close(10)

  periodic_boundary = .TRUE.
  if (N<0) then
    periodic_boundary = .FALSE.
  endif

  if(periodic_boundary) then
    write(fname_prefix, "(A,A,I0,A,I0,A,I0)") trim(project_dir), "/", p,"_",q,"_modulo_", N
  else
    write(fname_prefix, "(A,A,I0,A,I0,A,I0)") trim(project_dir), "/", p,"_",q,"_open_", -N
  endif

  !#####################################################################
  ! Reading of spectral information file *.info                        # 
  !#####################################################################
  
  spec_info_fname = trim(fname_prefix) // ".info"
  
  write(*,*) "Loading: ", trim(spec_info_fname)
  open(unit=10, file=trim(spec_info_fname))
  read(10,*,iostat=iostatus) hdim
  write(*, *) "Hilbert space dimension: ", hdim 
  close(10)

  !#####################################################################
  ! Setup of the Hamiltonian                                           # 
  !#####################################################################
  
  hamiltonian_fname = trim(fname_prefix) // ".hamiltonian"
  write(*,*) "Loading: ", trim(hamiltonian_fname)

  allocate(hamiltonian(hdim,hdim))
  
  !-- loop over file until EOF
  
  iostatus = 0
  
  open(unit=10, file=trim(hamiltonian_fname))
  do while(.TRUE.)
    read(10, '(I10, A1, I10, A1, EN15.6)', iostat=iostatus) i, char_buffer, j, char_buffer, val
    if(iostatus==0) then
      hamiltonian(i+1,j+1) = val
    else if(iostatus>0) then
      write(*, *) "A horrible error has occured!"
    else
      write(*, *) "Hamiltonian loaded succesfully!"
      exit
    endif    
  enddo 
  close(10)

  !#####################################################################
  ! Initialize BLACS                                                   # 
  !#####################################################################

  ! -- Returns the number of processes available for use
  call BLACS_PINFO(IAM, NPROCS)
  if ( NPROCS.LT.1 ) then
    ! -- Allocates virtual machine and spawns processes
    call BLACS_SETUP( IAM, NPROW*NPCOL )
  endif

  ! -- No idea what this does
  call BLACS_GET( -1, 0, CONTEXT )

  ! -- Assigns available processes into BLACS process grid.
  call BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )

  ! -- Assigns a grid position to the calling rank
  call BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )


  bailing_out = MYROW.EQ.-1

  if(not(bailing_out)) then

    ! -- LAPACK vars
    LDA = max(1, hdim)
    JOBZ = 'N'
    UPLO = 'U'
    
    ! -- array descriptors
    CALL DESCINIT( DESCA, hdim, hdim, NB, NB, 0, 0, CONTEXT, LDA, INFO )
    CALL DESCINIT( DESCZ, hdim, hdim, NB, NB, 0, 0, CONTEXT, LDA, INFO )

    ! --  distribute matrix to process grid

    !CALL PDLAMODHILB( N, A, 1, 1, DESCA, INFO )

    CALL PDSYEV( JOBZ, UPLO, hdim, hamiltonian, 1, 1, DESCA, W, Z, 1, 1,
                 DESCZ, WORK, LWORK, INFO )

    CALL BLACS_GRIDEXIT( CONTEXT )
  endif

  ! --  release all BLACS context and memory allocated by the BLACS
  call BLACS_EXIT( 0 )

  

  ! -- allocating space for eigenvalues
  allocate(W(hdim))

  

  ! -- automatically determine optimal size of work memory
  LWORK = -1
  call DSYEV(JOBZ, UPLO, hdim, hamiltonian, LDA, W, WORK_QUERY, LWORK, INFO)
  LWORK = max(1, 3*hdim-1)
  LWORK = max(LWORK, WORK_QUERY(1))

  ! write(*,*) "LWORK: ", LWORK, WORK_QUERY(1), hdim

  ! -- allocate work memory
  allocate(WORK(LWORK))

  ! -- calculate all eigenvalues
  call DSYEV(JOBZ, UPLO, hdim, hamiltonian, LDA, W, WORK, LWORK, INFO)

  deallocate(WORK)
  deallocate(hamiltonian)

  if(INFO==0) then
    write(*,*) "Diagonalization succesful!"
  endif
 
  !#####################################################################
  ! Writing the eigenvalues to file
  !#####################################################################

  eigen_fname = trim(fname_prefix) // ".eig"
  
  ! -- store as binary 
  call save_npy(trim(eigen_fname), W)

  deallocate(W)
end program
