      module geofem_util
        implicit none
        public
        integer,parameter:: geofem_name_len = 80
        integer:: kreal
        integer:: kint
!        parameter (kreal = selected_real_kind(10))
        parameter (kreal= 8)
        parameter (kint = 4) 

        integer:: input_datum_initiated = 0

!! LOCAL MESH info.
        type local_mesh
! NODE info
          integer n_node
          real(kind=kreal),pointer:: node(:,:)
          integer n_elem
! ELEMENT info
          integer,pointer:: elem_type(:)
          integer,pointer:: elem_mat(:)
          real(kind=kreal),pointer:: elem_cond(:)
          integer,pointer:: elem(:,:)
! PE info
          integer n_neighbor_pe
          integer,pointer:: neighbor_pe(:)
          integer           n_internal
          integer,pointer:: import_index(:)
          integer,pointer:: import_node(:)
          integer,pointer:: export_index(:)
          integer,pointer:: export_node(:)
          integer,pointer:: global_node_id(:)
          integer,pointer:: global_elem_id(:)
        end type local_mesh

!! interface structures for input of Analysis
        type node_elem_grp
          integer n_enum_grp
          character(geofem_name_len),pointer:: enum_grp_name(:)
          integer,pointer::                    enum_grp_index(:)
          integer,pointer::                    enum_grp_node(:)
        end type node_elem_grp

        type grp_data
          type(node_elem_grp) node_grp
          type(node_elem_grp) elem_grp
! surface group
          integer n_surf_grp
          character(geofem_name_len),pointer:: surf_grp_name(:)
          integer,pointer:: surf_grp_index(:)
          integer,pointer:: surf_grp_node(:,:)
        end type grp_data

!
! grp_data_buf%node_grp
! grp_data_buf%elem_grp
! grp_data_buf%
!
      end module geofem_util
