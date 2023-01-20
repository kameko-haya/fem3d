subroutine input_ctrl
  use m_util
  use m_solver

  implicit none
  integer :: i, icel

  open (11,file= 'input.dat', status='unknown')
  read (11,*) NE
  read (11,*) dX, F, Area, Young
  read (11,*) IterMax
  read (11,*) Resid
  close (11)
  N = NE+1

  allocate(U(N))
  allocate(X(N))
  allocate(Diag(N))
  allocate(Amat(2*N-2))
  allocate(Rhs(N))

  allocate(index(0:N))
  allocate(item(2*N-2))
  allocate(icelnode(2*NE))

  allocate(W(N,4))

  do i=1, N
     U(i)=0.d0
     Diag(i)=0.d0
     Rhs(i)=0.d0
     X(i) = i*dX
  end do 
  do i=1, 2*N-2
     Amat(i)=0.d0
  end do 
  do icel=1, NE
     icelnode(2*icel-1)=icel
     icelnode(2*icel  )=icel+1
  end do

  Kmat(1,1) = 1.d0
  Kmat(2,1) =-1.d0
  Kmat(1,2) =-1.d0
  Kmat(2,2) = 1.d0

end subroutine input_ctrl
