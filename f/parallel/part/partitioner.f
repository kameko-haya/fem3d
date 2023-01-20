      module partitioner
        use geofem_util
        implicit REAL*8 (A-H,O-Z)
        integer(kind=kint), dimension (:)  , allocatable :: RHO
        integer(kind=kint), dimension (:,:), allocatable :: STACK_EXPORT
        integer(kind=kint), dimension (:,:), allocatable :: STACK_IMPORT

        integer(kind=kint), dimension (:,:), allocatable ::             &
     &                      HOME_NODE, HOME_ELEM, HOME_EDGE
        integer(kind=kint), dimension (:)  , allocatable ::             &
     &                      ELEM_internal_LIST, nELEM_internal
        integer(kind=kint), dimension (:)  , allocatable ::             &
     &                      EDGE_internal_LIST, nEDGE_internal

        integer(kind=kint), dimension (:)  , allocatable ::             &
     &           ELMGRPITEM, NODGRPITEM,                                &
     &           NEIBNODTOT, IWORK, IACTEDG, IEDGFLAG, IMASK,           &
     &           IDEAD, ISTACK, IGROUP, ICOND1, ICOND2

        integer(kind=kint), dimension (:,:), allocatable ::             &
     &           SUFGRPITEM, IEDGNOD, NEIBNOD, NEIBPE

        integer(kind=kint), pointer::                                   &
     &                              NODGRPITEMG(:), ELMGRPITEMG(:),     &
     &                              SUFGRPITEMG(:,:),                   &
     &                              ELMGRPSTACKG(:), NODGRPSTACKG(:),   &
     &                              SUFGRPSTACKG(:)

        integer(kind=kint ), pointer::  ICELNOD(:,:), IELMTYP(:)
        integer(kind=kint ), dimension(:,:), allocatable:: ICELNODG
        integer(kind=kint ), pointer::  IELMMAT(:)
        real   (kind=kreal), pointer::  XYZ(:,:), CONDW(:)

        integer(kind=kint ), dimension (:) , allocatable ::             &
     &           NPN, NPNID, ISTACKN, NPC, NPCID, ISTACKC, NEIBPETOT,   &
     &           NODTOT, INTNODTOT, INODLOCAL, NOD_EXPORT, NOD_IMPORT,  &
     &           ELMGRPSTACK, NODGRPSTACK , SUFGRPSTACK, NPE

        integer(kind=kint ), dimension (:) , allocatable ::  NODELM

        integer(kind=kint )                   ::                        &
     &           RHOMAX, RHOMIN, ELMGRPTOT, NODGRPTOT, SUFGRPTOT,       &
     &           N, NP, NPOWER, NTYP

        character(geofem_name_len),pointer::                            &
     &                NODGRPNAME(:), ELMGRPNAME(:), SUFGRPNAME(:)
      

        character (len=80) :: GRIDFIL, METISFIL, HEADER, HEADW
        character (len=80), dimension(:), allocatable, save :: FILNAME
        character (len=80), dimension(:), allocatable, save :: WORKFIL

        integer(kind=kint )                   ::                        &
     &         IOPTFLAG, ITERMAX, IEDGCUT, IEDGTOT, IELMTOT, IACTEDGTOT,&
     &         HOWmanyADAPTATIONs, CoarseGridLevels,                    &
     &         WhenIwasRefinedN, WhenIwasRefinedE,                      &
     &         WhereIwas, adapt_type, adapt_level,                      &
     &         adapt_par1, adapt_par2, adapt_chi1, adapt_chi2,          &
     &         adapt_par_type, IELMTOTG

      end module partitioner
