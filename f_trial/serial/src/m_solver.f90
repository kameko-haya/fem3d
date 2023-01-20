module m_solver

  integer,allocatable,dimension(:) :: index, item
  integer,allocatable,dimension(:) :: icelnode
  real,allocatable,dimension(:) :: U, X, Diag, AMat, Rhs
  real,allocatable,dimension(:,:) :: W
  real,dimension(2,2) :: Kmat

end module m_solver
