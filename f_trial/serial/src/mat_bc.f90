subroutine mat_bc

  use m_util
  use m_solver

  integer :: i, j 
  

  !C at X = 0, u=0
  i=1
  js = item(i-1)+1
  je = item(i  )
  do j = js, je
    Amat(j) =0.0d
  enddo
  Diag(i)=1.0d
  Rhs(i)=0

  !C at X = Xmax, F

end subroutine mat_bc
