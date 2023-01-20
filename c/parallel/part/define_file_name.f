!C
!C***
!C*** DEFINE_FILE_NAME
!C***
!C          
      subroutine DEFINE_FILE_NAME (HEADER, my_rank, filname)

      character (len=80) ::  filname
      character (len=80) ::  HEADER
      character (len= 1) ::  SUBindex1
      character (len= 2) ::  SUBindex2
      character (len= 3) ::  SUBindex3
      character (len= 4) ::  SUBindex4
      character (len= 5) ::  SUBindex5
      character (len= 6) ::  SUBindex6
    
      HEADER= adjustL (HEADER)
      LENGTH= len_trim(HEADER)

      if (my_rank.le.9) then
        ID= 1
        write(SUBindex1 ,'(i1.1)') my_rank
       else if (my_rank.le.99) then
        ID= 2
        write(SUBindex2 ,'(i2.2)') my_rank
       else if (my_rank.le.999) then
        ID= 3
        write(SUBindex3 ,'(i3.3)') my_rank
       else if (my_rank.le.9999) then
        ID= 4
        write(SUBindex5 ,'(i4.4)') my_rank
       else if (my_rank.le.99999) then
        ID= 5
        write(SUBindex6 ,'(i5.5)') my_rank
       else if (my_rank.le.99999) then
        ID= 6
        write(SUBindex4 ,'(i6.6)') my_rank
      endif

      if (ID.eq.1) filname= HEADER(1:LENGTH)//'.'//SUBindex1
      if (ID.eq.2) filname= HEADER(1:LENGTH)//'.'//SUBindex2
      if (ID.eq.3) filname= HEADER(1:LENGTH)//'.'//SUBindex3
      if (ID.eq.4) filname= HEADER(1:LENGTH)//'.'//SUBindex4
      if (ID.eq.5) filname= HEADER(1:LENGTH)//'.'//SUBindex5
      if (ID.eq.6) filname= HEADER(1:LENGTH)//'.'//SUBindex6

      return
      end
