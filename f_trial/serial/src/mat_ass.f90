subroutine mat_ass_main

    use m_util
    use m_solver
    implicit none
    
    integer :: i, j, icel, k1, k2, in1, in2
    real, dimension(2,2):: Emat


    do icel = 1, NE

       in1 = icelnode(icel*2-1)
       in2 = icelnode(icel*2  )

       X1 = X(in1)
       X2 = X(in2)

       DL = abs(X2-X1)
       Ck = Area * Young / DL
       !c Kmat(col, row)
       Emat(1,1)=Ck*Kmat(1,1)
       Emat(2,1)=Ck*Kmat(2,1)
       Emat(1,2)=Ck*Kmat(1,2)
       Emat(2,2)=Ck*Kmat(2,2)

       Diag(in1) = Diag(in1) + Emat(1,1)
       Diag(in2) = Diag(in2) + Emat(2,2)

       k1 = index(in1  )
       k2 = index(in2-1)+1

       print *, "icel =",  icel, "Amat_index=", k1, k2
       Amat(k1) = Amat(k1)+Emat(2,1)
       Amat(k2) = Amat(k2)+Emat(1,2)
    end do
    print *, "Amat" 
    do i =1, 2*N-2
     print *, Amat(i)
    enddo
    print *, "Diag" 
    do i =1, N
     print *, Diag(i)
    enddo
end subroutine mat_ass_main
