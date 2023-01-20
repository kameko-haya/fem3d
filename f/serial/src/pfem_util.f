!C
!C***
!C*** pfem_util
!C***

      module pfem_util
      include 'precision.inc'
!C
!C +--------------+
!C | file settings |
!C +--------------+
!C===
      character(len= 80) ::  fname
!C===

!C
!C +-------------------------+
!C | DISTRIBUTED MESH FILE's |
!C +-------------------------+
!C===

!C
!C-- CONNECTIVITIES & BOUNDARY nodes
      integer:: ICELTOT, NODGRPtot

      real(kind=kreal)  , dimension(:,:),allocatable :: XYZ
      integer(kind=kint), dimension(:,:),allocatable :: iCELNOD

      integer(kind=kint), dimension(:),  allocatable ::                 &
     &                                   NODGRP_INDEX, NODGRP_ITEM

      character(len=80), dimension(:),  allocatable :: NODGRP_NAME
!C===

!C
!C +-----------------+
!C | MATRIX & SOLVER |
!C +-----------------+
!C===

!C
!C-- MATRIX SCALARs
      integer(kind=kint) :: N, NP, N2, NL, NU, NPL, NPU

!C
!C-- MATRIX arrays
      real(kind=kreal), dimension(:), allocatable ::  D, B, X, ALUG
      real(kind=kreal), dimension(:), allocatable ::  AL,AU

      integer(kind=kint), dimension(:), allocatable :: indexL, indexU
      integer(kind=kint), dimension(:), allocatable ::  itemL,  itemU

      integer(kind=kint), dimension(:),  allocatable :: INL, INU
      integer(kind=kint), dimension(:,:),allocatable :: IAL, IAU

      integer(kind=kint), dimension(:,:),allocatable :: IWKX
!C
!C-- PARAMETER's for LINEAR SOLVER
      integer(kind=kint) :: METHOD, PRECOND, NSET, iterPREmax
      integer(kind=kint) :: NREST
      integer(kind=kint) :: ITER, ITERactual

      real   (kind=kreal) :: RESID, SIGMA_DIAG, SIGMA
!C===

!C
!C +-------------+
!C | PARAMETER's |
!C +-------------+
!C===

!C
!C-- GENERAL PARAMETER's
      integer(kind=kint ), dimension(100) :: pfemIarray
      real   (kind=kreal), dimension(100) :: pfemRarray

      real   (kind=kreal), parameter :: O8th= 0.125d0
!C
!C-- PARAMETER's for FEM
      real(kind=kreal), dimension(2,2,8) :: PNQ, PNE, PNT
      real(kind=kreal), dimension(2)     :: WEI, POS
      integer(kind=kint), dimension(100) :: NCOL1, NCOL2

      real(kind=kreal), dimension(2,2,2,8) :: SHAPE
      real(kind=kreal), dimension(2,2,2,8) :: PNX, PNY, PNZ
      real(kind=kreal), dimension(2,2,2  ) :: DETJ

!C
!C-- PROBLEM PARAMETER's
      real(kind=kreal):: ELAST, POISSON
      real(kind=kreal), dimension(:),allocatable :: SIGMA_N, TAU_N
!C===
      end module pfem_util
