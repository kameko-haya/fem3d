subroutine mat_con

   use m_solver
   use m_util

   implicit none
   integer :: i, j, js, je
   do i=1, N+1
     index(i) =2
   end do
   index(0) =0
   index(1) =1
   index(N) =1
!   print *, "N=", N
!   print *, "index"
!   do i =0, N
!     print *, index(i)
!   enddo

   do i =1, N
     index(i)= index(i) + index(i-1)
   end do

!   print *, "index"
!   do i =0, N
!     print *, index(i)
!   enddo

   NPLU = index(N)

   do i= 1, N
      js= index(i-1)+1
      je= index(i  )
      !print *, "col=", i, "js, je=", js, je
      if (i .eq. 1) then
         item(js)=i+1
         print *, "col=",i, "item(",js,")=", item(js) 
      else if (i .eq. N) then
         item(js)=i-1
         print *, "col=",i, "item(",js,")=", item(js)
      else
         item(js)   = i-1
         item(je)   = i+1
         print *, "col=",i, "item(",js,")=", item(js), "item(",je,")=",item(je)
      endif
   enddo
   do i = 1, 2*N-2
      print *, item(i)
   enddo
end subroutine mat_con
