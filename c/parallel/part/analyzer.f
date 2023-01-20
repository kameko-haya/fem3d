      module  analyzer

      use  geofem_util
      use  partitioner

      type (local_mesh) :: local_mesh_buf
      type (grp_data)   ::   grp_data_buf

      contains

      subroutine  init_analyzer ( errno )
      integer(kind=kint)    :: errno
      character*80 GRIDFIL      
      character*1  LINE(132)

!C
!C +--------------------+
!C | Original GRID FILE |
!C +--------------------+
!C===
      write (*,*      ) 'Original GRID-FILE ?'
      read  (*,'(a80)')  GRIDFIL

      IUNIT= 11
      open (IUNIT,file=GRIDFIL,status='unknown',form='formatted')

!C
!C-- NODE info.
        read (IUNIT,'(10i10)') local_mesh_buf%n_node
        local_mesh_buf%n_internal= local_mesh_buf%n_node
        if (local_mesh_buf%n_node.ne.local_mesh_buf%n_internal)         &
     &      call ERROR_EXIT (2)

        IELMDMY= 8

        allocate (local_mesh_buf%node(local_mesh_buf%n_node,3))

        write(* ,'(  " * INODTOT =",i8)') local_mesh_buf%n_node
        write(* ,'(  " * GRID")')

        do i= 1, local_mesh_buf%n_node
          read (IUNIT,*) ii,(local_mesh_buf%node(i,k),k =1,3)            
        enddo

!C
!C-- ELEMENT Info.
        read (IUNIT,*) local_mesh_buf%n_elem

        write(* ,'(  " * IELMTOT =",i8)') local_mesh_buf%n_elem
        write(* ,'(  " * ELM")')

        allocate (local_mesh_buf%elem     (local_mesh_buf%n_elem,8))
        allocate (local_mesh_buf%elem_type(local_mesh_buf%n_elem  ))
        allocate (local_mesh_buf%elem_mat (local_mesh_buf%n_elem  ))
        allocate (local_mesh_buf%elem_cond(local_mesh_buf%n_elem  ))

        read (IUNIT,'(10i10)')                                          &
     &       (local_mesh_buf%elem_type(i),i=1,local_mesh_buf%n_elem)

        do icel= 1, local_mesh_buf%n_elem
          ityp= local_mesh_buf%elem_type(icel)

          if (ityp.eq.341) IELMDMY= 4
          if (ityp.eq.351) IELMDMY= 6
          if (ityp.eq.361) IELMDMY= 8
          read (IUNIT,'(10i10)')                                         &
     &          icc,      local_mesh_buf%elem_mat(icel),                 &
     &                   (local_mesh_buf%elem(icel,k),k=1,IELMDMY)          
        enddo

!C
!C-- IMPORT/EXPORT

!C
!C-- BOUNDARY Info. : NODE group
      N2= max (local_mesh_buf%n_elem*2, local_mesh_buf%n_node*2)

      allocate (grp_data_buf%node_grp%enum_grp_index(0:100))
      allocate (grp_data_buf%node_grp%enum_grp_name (0:100))
      allocate (grp_data_buf%node_grp%enum_grp_node(N2))

      allocate (grp_data_buf%elem_grp%enum_grp_index(0:100))
      allocate (grp_data_buf%elem_grp%enum_grp_name (0:100))
      allocate (grp_data_buf%elem_grp%enum_grp_node(N2))

      allocate (grp_data_buf%surf_grp_index(0:100))
      allocate (grp_data_buf%surf_grp_name (0:100))
      allocate (grp_data_buf%surf_grp_node (N2,2))

        write(* ,'(  " * BOUNDARY : NODE group")')
        grp_data_buf%node_grp%enum_grp_index(0)= 0
          read  (IUNIT,*) grp_data_buf%node_grp%n_enum_grp
          read  (IUNIT,'(10i10)')                                       & 
     &                  (grp_data_buf%node_grp%enum_grp_index(ig),      &
     &                   ig= 1,grp_data_buf%node_grp%n_enum_grp)
          do ig= 1, grp_data_buf%node_grp%n_enum_grp
              read  (IUNIT,'(a80)')                                     &
     &               grp_data_buf%node_grp%enum_grp_name(ig)
              write  (*, '(a64)')                                       &
     &               grp_data_buf%node_grp%enum_grp_name(ig)
              read (IUNIT,'(10i10)')                                    &
     &             (grp_data_buf%node_grp%enum_grp_node(is),            &
     &              is= grp_data_buf%node_grp%enum_grp_index(ig-1)+1,   &
     &                  grp_data_buf%node_grp%enum_grp_index(ig))
          enddo
!C===

!C
!C-- POINTER copy and allocation
      XYZ     => local_mesh_buf%node
      ICELNOD => local_mesh_buf%elem
      IELMTYP => local_mesh_buf%elem_type
      IELMMAT => local_mesh_buf%elem_mat

      CONDW   => local_mesh_buf%elem_cond

      N      = local_mesh_buf%n_node
      IELMTOT= local_mesh_buf%n_elem

      allocate (ICELNODG(IELMTOT,8))
      ICELNODG= ICELNOD

      if (      N.le.0) call ERROR_EXIT(1001,0)
      if (IELMTOT.le.0) call ERROR_EXIT(1001,0)

      NODGRPNAME => grp_data_buf%node_grp%enum_grp_name
      NODGRPITEMG => grp_data_buf%node_grp%enum_grp_node
      NODGRPSTACKG => grp_data_buf%node_grp%enum_grp_index
      NODGRPTOT = grp_data_buf%node_grp%n_enum_grp

      if (NODGRPTOT.lt.0) call ERROR_EXIT(1002,1)
 
      if (NODGRPTOT.gt.0) then
        do is= 1, NODGRPSTACKG(NODGRPTOT)
          in= NODGRPITEMG(is)
          if (in.le.0) call ERROR_EXIT(1003,1)
          if (in.gt.N) call ERROR_EXIT(2002,1)
        enddo
      endif

      allocate (RHO(local_mesh_buf%n_node))
      allocate (NODELM(local_mesh_buf%n_elem))

      ELMGRPTPT= 0
      SUFGRPTPT= 0

      allocate (NODGRPSTACK(0:NODGRPTOT))
      allocate (NODGRPITEM (NODGRPSTACKG(NODGRPTOT)))


!C
!C +--------------+
!C | ELEMENT-TYPE |
!C +--------------+
!C
!C   ROD :             111  1-2
!C                     112  1-3:2
!C   2D  : triangle    211  1-2-3
!C                     212  1-2-3:4-5-6
!C   2D  : quad.       221  1-2-3-4
!C                     222  1-2-3-4:5-6-7-8
!C   3D  : tet.        311  1-2-3-4
!C                     312  1-2-3-4:5-6-7:8-9-10
!C   3D  : prism       321  1-2-3-4-5-6
!C                     322  1-2-3-4-5-6:7-8-9:10-11-12:13-14-15
!C   3D  : hexa.       331  1-2-3-4-5-6-7-8
!C                     332  1-2-3-4-5-6-7-8:9-10-11-12:13-14-15-16:17-18-19-20
!C
!C   master-slave(tri) 411  1:2-3-4 
!C                     412  1:2-3-4:5-6-7
!C   master-slave(quad)421  1:2-3-4-5 
!C                     422  1:2-3-4-6:6-7-8-9
!C
!C   joint (tri)       511 (1-2-3)*(4-5-6)
!C                     512 (1-2-3:4-5-6)*(7-8-9:10-11-12)
!C   joint (quad)      521 (1-2-3-4)*(5-6-7-8)
!C                     522 (1-2-3-4:5-6-7-8)*(9-10-11-12:13-14-15-16)
!C
!C   beam              611  1-2
!C                     612  1-2:3
!C
!C   shell: triangle   711  1-2-3
!C                     712  1-2-3:4-5-6
!C   shell: quad.      721  1-2-3-4
!C                     722  1-2-3-4:5-6-7-8
!C===   
        NODELM= 0
        do icel= 1, local_mesh_buf%n_elem
          ityp= local_mesh_buf%elem_type(icel)
          if (ityp.le.  0) call ERROR_EXIT(1004, icel)
          if (ityp.eq.111) NODELM(icel)=  2
          if (ityp.eq.112) NODELM(icel)=  3
          if (ityp.eq.211) NODELM(icel)=  3
          if (ityp.eq.212) NODELM(icel)=  6
          if (ityp.eq.221) NODELM(icel)=  4
          if (ityp.eq.222) NODELM(icel)=  8
          if (ityp.eq.311) NODELM(icel)=  4
          if (ityp.eq.312) NODELM(icel)= 10
          if (ityp.eq.321) NODELM(icel)=  6
          if (ityp.eq.322) NODELM(icel)= 15
          if (ityp.eq.361) NODELM(icel)=  8
          if (ityp.eq.332) NODELM(icel)= 20
          if (ityp.eq.411) NODELM(icel)=  4
          if (ityp.eq.412) NODELM(icel)=  7
          if (ityp.eq.421) NODELM(icel)=  5
          if (ityp.eq.422) NODELM(icel)=  9
          if (ityp.eq.511) NODELM(icel)=  6
          if (ityp.eq.512) NODELM(icel)= 12
          if (ityp.eq.521) NODELM(icel)=  8
          if (ityp.eq.522) NODELM(icel)= 16
          if (ityp.eq.611) NODELM(icel)=  2
          if (ityp.eq.612) NODELM(icel)=  3
          if (ityp.eq.711) NODELM(icel)=  3
          if (ityp.eq.712) NODELM(icel)=  6
          if (ityp.eq.721) NODELM(icel)=  4
          if (ityp.eq.722) NODELM(icel)=  8

          if (NODELM(icel).eq.0) call ERROR_EXIT(33,icel)

          do k= 1, NODELM(icel)
            in= ICELNOD(icel,k)
            if (in.le.0) call ERROR_EXIT(1005,icel)
            if (in.gt.N) call ERROR_EXIT(2001,icel)
          enddo
        enddo
!C
!C== EDGE information
      NE= 4*N

  100 continue
      allocate (IEDGNOD(NE,2)) 
      IEDGTOT= 0
      IEDGNOD= 0

      do icel= 1, IELMTOT
        ityp= local_mesh_buf%elem_type(icel)
!C
!C-- ROD or BEAM
        if (ityp.eq.111 .or. ityp.eq.611) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          call EDGE_INFO (in1,in2, iedge, 0)
        endif        
        if (ityp.eq.112 .or. ityp.eq.612) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in3,in2, iedge, 0)
        endif        
!C
!C-- 2D or Shell : triangle
        if (ityp.eq.211 .or. ityp.eq.711) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in1, iedge, 0)
        endif        
        if (ityp.eq.212 .or. ityp.eq.712) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in4,in2, iedge, 0)
          call EDGE_INFO (in2,in5, iedge, 0)
          call EDGE_INFO (in5,in3, iedge, 0)
          call EDGE_INFO (in3,in6, iedge, 0)
          call EDGE_INFO (in6,in1, iedge, 0)
        endif        
!C
!C-- 2D or Shell : quadrilateral
        if (ityp.eq.221 .or. ityp.eq.721) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in1, iedge, 0)
        endif        

        if (ityp.eq.222 .or. ityp.eq.722) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          in6= ICELNOD(icel,7)
          in6= ICELNOD(icel,8)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in5,in2, iedge, 0)
          call EDGE_INFO (in2,in6, iedge, 0)
          call EDGE_INFO (in6,in3, iedge, 0)
          call EDGE_INFO (in3,in7, iedge, 0)
          call EDGE_INFO (in7,in4, iedge, 0)
          call EDGE_INFO (in4,in8, iedge, 0)
          call EDGE_INFO (in8,in1, iedge, 0)
        endif        
!C
!C-- 3D : tetrahedron
        if (ityp.eq.311) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in2, iedge, 0)
        endif        

        if (ityp.eq.312) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in8= ICELNOD(icel, 8)
          in9= ICELNOD(icel, 9)
          in0= ICELNOD(icel,10)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in5,in2, iedge, 0)
          call EDGE_INFO (in1,in6, iedge, 0)
          call EDGE_INFO (in6,in3, iedge, 0)
          call EDGE_INFO (in1,in7, iedge, 0)
          call EDGE_INFO (in7,in4, iedge, 0)
          call EDGE_INFO (in2,in8, iedge, 0)
          call EDGE_INFO (in8,in3, iedge, 0)
          call EDGE_INFO (in3,in9, iedge, 0)
          call EDGE_INFO (in9,in4, iedge, 0)
          call EDGE_INFO (in4,in0, iedge, 0)
          call EDGE_INFO (in0,in2, iedge, 0)
        endif        
!C
!C-- 3D : prisms
        if (ityp.eq.321) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in1, iedge, 0)
          call EDGE_INFO (in4,in5, iedge, 0)
          call EDGE_INFO (in5,in6, iedge, 0)
          call EDGE_INFO (in6,in4, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in5, iedge, 0)
          call EDGE_INFO (in3,in6, iedge, 0)
        endif        

        if (ityp.eq.322) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in8= ICELNOD(icel, 8)
          in9= ICELNOD(icel, 9)
          in0= ICELNOD(icel,10)
          ina= ICELNOD(icel,11)
          inb= ICELNOD(icel,12)
          inc= ICELNOD(icel,13)
          ind= ICELNOD(icel,14)
          ine= ICELNOD(icel,15)
          call EDGE_INFO (in1,in7, iedge, 0)
          call EDGE_INFO (in7,in2, iedge, 0)
          call EDGE_INFO (in2,in8, iedge, 0)
          call EDGE_INFO (in8,in3, iedge, 0)
          call EDGE_INFO (in3,in9, iedge, 0)
          call EDGE_INFO (in9,in1, iedge, 0)
          call EDGE_INFO (in4,in0, iedge, 0)
          call EDGE_INFO (in0,in5, iedge, 0)
          call EDGE_INFO (in5,ina, iedge, 0)
          call EDGE_INFO (ina,in6, iedge, 0)
          call EDGE_INFO (in6,inb, iedge, 0)
          call EDGE_INFO (inb,in4, iedge, 0)
          call EDGE_INFO (in1,inc, iedge, 0)
          call EDGE_INFO (inc,in4, iedge, 0)
          call EDGE_INFO (in2,ind, iedge, 0)
          call EDGE_INFO (ind,in5, iedge, 0)
          call EDGE_INFO (in3,ine, iedge, 0)
          call EDGE_INFO (ine,in6, iedge, 0)
        endif        
!C
!C-- 3D : hexahedron
        if (ityp.eq.361) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          in7= ICELNOD(icel,7)
          in8= ICELNOD(icel,8)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in1, iedge, 0)
          call EDGE_INFO (in5,in6, iedge, 0)
          call EDGE_INFO (in6,in7, iedge, 0)
          call EDGE_INFO (in7,in8, iedge, 0)
          call EDGE_INFO (in8,in5, iedge, 0)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in2,in6, iedge, 0)
          call EDGE_INFO (in3,in7, iedge, 0)
          call EDGE_INFO (in4,in8, iedge, 0)
        endif        

        if (ityp.eq.332) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in8= ICELNOD(icel, 8)
          in9= ICELNOD(icel, 9)
          in0= ICELNOD(icel,10)
          ina= ICELNOD(icel,11)
          inb= ICELNOD(icel,12)
          inc= ICELNOD(icel,13)
          ind= ICELNOD(icel,14)
          ine= ICELNOD(icel,15)
          ind= ICELNOD(icel,16)
          ing= ICELNOD(icel,17)
          inh= ICELNOD(icel,18)
          ini= ICELNOD(icel,19)
          inj= ICELNOD(icel,20)
          call EDGE_INFO (in1,in9, iedge, 0)
          call EDGE_INFO (in9,in2, iedge, 0)
          call EDGE_INFO (in2,in0, iedge, 0)
          call EDGE_INFO (in0,in3, iedge, 0)
          call EDGE_INFO (in3,ina, iedge, 0)
          call EDGE_INFO (ina,in4, iedge, 0)
          call EDGE_INFO (in4,inb, iedge, 0)
          call EDGE_INFO (inb,in1, iedge, 0)
          call EDGE_INFO (in5,inc, iedge, 0)
          call EDGE_INFO (inc,in6, iedge, 0)
          call EDGE_INFO (in6,ind, iedge, 0)
          call EDGE_INFO (ind,in7, iedge, 0)
          call EDGE_INFO (in7,ine, iedge, 0)
          call EDGE_INFO (ine,in8, iedge, 0)
          call EDGE_INFO (in8,inf, iedge, 0)
          call EDGE_INFO (inf,in5, iedge, 0)
          call EDGE_INFO (in1,ing, iedge, 0)
          call EDGE_INFO (ing,in5, iedge, 0)
          call EDGE_INFO (in2,inh, iedge, 0)
          call EDGE_INFO (inh,in6, iedge, 0)
          call EDGE_INFO (in3,ini, iedge, 0)
          call EDGE_INFO (ini,in7, iedge, 0)
          call EDGE_INFO (in4,inj, iedge, 0)
          call EDGE_INFO (inj,in8, iedge, 0)
        endif        

!C
!C-- Master-Slave : triangle
        if (ityp.eq.411) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in2, iedge, 0)
        endif        

        if (ityp.eq.412) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in5, iedge, 0)
          call EDGE_INFO (in5,in3, iedge, 0)
          call EDGE_INFO (in3,in6, iedge, 0)
          call EDGE_INFO (in6,in4, iedge, 0)
          call EDGE_INFO (in4,in7, iedge, 0)
          call EDGE_INFO (in7,in2, iedge, 0)
        endif        

!C
!C-- Master-Slave : quadrilateral
        if (ityp.eq.421) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in4= ICELNOD(icel,5)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in5, iedge, 0)
          call EDGE_INFO (in5,in2, iedge, 0)
        endif        

        if (ityp.eq.422) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in7= ICELNOD(icel, 8)
          in7= ICELNOD(icel, 9)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in1,in3, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in2,in6, iedge, 0)
          call EDGE_INFO (in6,in3, iedge, 0)
          call EDGE_INFO (in3,in7, iedge, 0)
          call EDGE_INFO (in7,in4, iedge, 0)
          call EDGE_INFO (in4,in8, iedge, 0)
          call EDGE_INFO (in8,in5, iedge, 0)
          call EDGE_INFO (in5,in9, iedge, 0)
          call EDGE_INFO (in9,in2, iedge, 0)
        endif        

!C
!C-- Joint : prisms
        if (ityp.eq.511) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in1, iedge, 0)
          call EDGE_INFO (in4,in5, iedge, 0)
          call EDGE_INFO (in5,in6, iedge, 0)
          call EDGE_INFO (in6,in4, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in5, iedge, 0)
          call EDGE_INFO (in3,in6, iedge, 0)
        endif        

        if (ityp.eq.512) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in8= ICELNOD(icel, 8)
          in9= ICELNOD(icel, 9)
          in0= ICELNOD(icel,10)
          ina= ICELNOD(icel,11)
          inb= ICELNOD(icel,12)
          call EDGE_INFO (in1,in7, iedge, 0)
          call EDGE_INFO (in7,in2, iedge, 0)
          call EDGE_INFO (in2,in8, iedge, 0)
          call EDGE_INFO (in8,in3, iedge, 0)
          call EDGE_INFO (in3,in9, iedge, 0)
          call EDGE_INFO (in9,in1, iedge, 0)
          call EDGE_INFO (in4,in0, iedge, 0)
          call EDGE_INFO (in0,in5, iedge, 0)
          call EDGE_INFO (in5,ina, iedge, 0)
          call EDGE_INFO (ina,in6, iedge, 0)
          call EDGE_INFO (in6,inb, iedge, 0)
          call EDGE_INFO (inb,in4, iedge, 0)
          call EDGE_INFO (in1,in4, iedge, 0)
          call EDGE_INFO (in2,in5, iedge, 0)
          call EDGE_INFO (in3,in6, iedge, 0)
        endif        
!C
!C-- 3D : hexahedron
        if (ityp.eq.521) then
          in1= ICELNOD(icel,1)
          in2= ICELNOD(icel,2)
          in3= ICELNOD(icel,3)
          in4= ICELNOD(icel,4)
          in5= ICELNOD(icel,5)
          in6= ICELNOD(icel,6)
          in7= ICELNOD(icel,7)
          in8= ICELNOD(icel,8)
          call EDGE_INFO (in1,in2, iedge, 0)
          call EDGE_INFO (in2,in3, iedge, 0)
          call EDGE_INFO (in3,in4, iedge, 0)
          call EDGE_INFO (in4,in1, iedge, 0)
          call EDGE_INFO (in5,in6, iedge, 0)
          call EDGE_INFO (in6,in7, iedge, 0)
          call EDGE_INFO (in7,in8, iedge, 0)
          call EDGE_INFO (in8,in5, iedge, 0)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in2,in6, iedge, 0)
          call EDGE_INFO (in3,in7, iedge, 0)
          call EDGE_INFO (in4,in8, iedge, 0)
        endif        

        if (ityp.eq.522) then
          in1= ICELNOD(icel, 1)
          in2= ICELNOD(icel, 2)
          in3= ICELNOD(icel, 3)
          in4= ICELNOD(icel, 4)
          in5= ICELNOD(icel, 5)
          in6= ICELNOD(icel, 6)
          in7= ICELNOD(icel, 7)
          in8= ICELNOD(icel, 8)
          in9= ICELNOD(icel, 9)
          in0= ICELNOD(icel,10)
          ina= ICELNOD(icel,11)
          inb= ICELNOD(icel,12)
          inc= ICELNOD(icel,13)
          ind= ICELNOD(icel,14)
          ine= ICELNOD(icel,15)
          ind= ICELNOD(icel,16)
          call EDGE_INFO (in1,in9, iedge, 0)
          call EDGE_INFO (in9,in2, iedge, 0)
          call EDGE_INFO (in2,in0, iedge, 0)
          call EDGE_INFO (in0,in3, iedge, 0)
          call EDGE_INFO (in3,ina, iedge, 0)
          call EDGE_INFO (ina,in4, iedge, 0)
          call EDGE_INFO (in4,inb, iedge, 0)
          call EDGE_INFO (inb,in1, iedge, 0)
          call EDGE_INFO (in5,inc, iedge, 0)
          call EDGE_INFO (inc,in6, iedge, 0)
          call EDGE_INFO (in6,ind, iedge, 0)
          call EDGE_INFO (ind,in7, iedge, 0)
          call EDGE_INFO (in7,ine, iedge, 0)
          call EDGE_INFO (ine,in8, iedge, 0)
          call EDGE_INFO (in8,inf, iedge, 0)
          call EDGE_INFO (inf,in5, iedge, 0)
          call EDGE_INFO (in1,in5, iedge, 0)
          call EDGE_INFO (in2,in6, iedge, 0)
          call EDGE_INFO (in3,in7, iedge, 0)
          call EDGE_INFO (in4,in8, iedge, 0)
        endif        

        if (IEDGTOT.ge.NE-24 .and. icel.lt.IELMTOT) then
          nn= IELMTOT/icel + 1
          NE= nn * NE + 1
          deallocate (IEDGNOD)
          goto 100
        endif
      enddo

      write(* ,'(  " * IEDGTOT =",2i8)') IEDGTOT, NE

      allocate (IACTEDG (IEDGTOT)) 
      allocate (IEDGFLAG(IEDGTOT))
      do ie= 1, IEDGTOT
        IACTEDG (ie)= ie
        IEDGFLAG(ie)= 0
      enddo

      IACTEDGTOT= IEDGTOT

      return

 998  continue
      call ERROR_EXIT (11,0)

 999  continue
      call ERROR_EXIT (12,0)

      end subroutine init_analyzer
      end module analyzer

!C
!C***
!C*** EDGE_INFO
!C***      

      subroutine EDGE_INFO (nod1,nod2,iedge,NFLAG)
      use partitioner
      parameter (nbuckets= 30000000)
      integer(kind=kint) ieaddrs (-nbuckets:nbuckets)

!C
!C     NFALG= 0 : CREATE NEW EDGEs
!C     NFLAG= 1 : REFER  the EDGE INFORMATION
!C

      iedge= 0
       
      nn1= mod(nod1, nbuckets) * mod(nod2, nbuckets)
      iarg= mod (nn1, nbuckets)

      if (NFLAG.eq.0) then
      if (ieaddrs (iarg).gt.IEDGTOT) then
        ieaddrs (iarg)= 0
      endif
      endif

   50 continue


!C
!C-- NEW EDGE

!      if (iarg.gt.nbuckets) write (*,*) nod1,nod2,iarg
      if (ieaddrs (iarg).eq.0) then
        iedgtot= iedgtot + 1
        iedge= iedgtot
        iedgnod (iedge,1)= nod1
        iedgnod (iedge,2)= nod2

!      if (iarg.gt.nbuckets) write (*,*) nod1,nod2,iarg
        ieaddrs (iarg)= iedgtot
        return
      else

!      if (iarg.gt.nbuckets) write (*,*) nod1,nod2,iarg
        iedge= ieaddrs (iarg)
        in1= iedgnod (iedge,1)
        in2= iedgnod (iedge,2)

!C
!C-- EXISTING EDGE
        if (in1.eq.nod1 .and. in2.eq.nod2  .or.                         &
     &      in1.eq.nod2 .and. in2.eq.nod1) return

        incr= 1
        ioldadd= iarg
  100   continue
        inewadd= mod (ioldadd + incr**3, nbuckets)

        if (inewadd .eq. ioldadd) then
          icount= icount+ 1
          ioldadd= ioldadd + 1
          inewadd= ioldadd
        endif

        if (NFLAG .eq. 0) then
        if (ieaddrs (inewadd).gt.IEDGTOT) then
            ieaddrs (inewadd)= 0
            goto 50
        endif
        endif

        if (ieaddrs (inewadd) .ne. 0) then
          iedge= ieaddrs (inewadd)
          in1= iedgnod (iedge,1)
          in2= iedgnod (iedge,2)
!C
!C-- EXISTING EDGE
          if (in1.eq.nod1 .and. in2.eq.nod2  .or.                       &
     &        in1.eq.nod2 .and. in2.eq.nod1) return
          incr= incr + 1
          go to 100

         else
!C
!C-- NEW EDGE
          iedgtot= iedgtot + 1
          iedge= iedgtot
          iedgnod (iedge,1)= nod1
          iedgnod (iedge,2)= nod2

!      if (inewadd.gt.nbuckets) write (*,*) nod1,nod2,inewadd
          ieaddrs (inewadd)= iedge
          return
        endif
      endif

      return
      end


      subroutine CHKCOM(LINE,IUNIT)
      character*1 LINE(132)

      do
        read (IUNIT,'(132a1)',err=998, end=999) LINE
!        iflag= 1
!        do i= 1, 30
!          if (LINE(i).ne.' ') iflag= 1          
!        enddo        
!        write (*,'(80a1)') (LINE(k),k=1,80)
!        write (*,*) iflag

        if   (LINE(1).ne.'!'.and.LINE(1).ne.'#') exit
!        if   (LINE(1).ne.'!'.and.LINE(1).ne.'#'.and.iflag.eq.1) exit
      enddo

      backspace (IUNIT)

      return

 998  continue
      call ERROR_EXIT (11)

 999  continue
      call ERROR_EXIT (12)

      end




