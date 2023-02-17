program exact_diagonalization
  implicit none;
  integer dim;
  integer, dimension(:,:), allocatable :: hamiltonian;

  dim = 10;

  allocate(hamiltonian(dim,dim));
  print *, "Hello World!"

  deallocate(hamiltonian);
end program
