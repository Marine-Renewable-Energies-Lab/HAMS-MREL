!  ------------------------------------------------------------------------------------------------------
!
!    Program HAMS-MREL for the diffraction and radiation of waves
!    by 3D structures.
!
!  License:
!
!    This routine is part of HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Vaibhav Raghavan
!
!  ------------------------------------------------------------------------------------------------------
MODULE SecondOrderTypes
   ! Options for the optional second-order (mean wave drift) calculation, read from SecondOrder.in.

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: so_options, panel_mesh, first_order_state, surface_geom

   TYPE :: so_options
      LOGICAL :: active = .FALSE.
      INTEGER :: mean_drift = 0          ! 0 off, 1 unidirectional, 2 all heading pairs
      LOGICAL :: do_far = .FALSE., do_near = .FALSE.
      INTEGER :: velocity_method = 1     ! 1 differencing, 2 source
      INTEGER :: waterline_method = 1    ! 1 extrapolate, 2 field
      REAL(8) :: waterline_offset = 0.01D0
      LOGICAL :: motion_fixed = .FALSE.
      INTEGER :: n_azimuth = 360
      REAL(8) :: far_tolerance = 1D-4
      INTEGER :: write_terms = 1         ! 0 none, 1 term file, 2 also panel diagnostics
      REAL(8) :: wl_z_tolerance = 1D-6   ! relative to body size
   END TYPE so_options

   TYPE :: panel_mesh
      ! Full panel set with symmetry images expanded. Panels 1..nbase are image 1, the next nbase are image 2.
      INTEGER :: n = 0, nnode = 0, nbody = 1, nbase = 0, nsys = 1
      REAL(8), ALLOCATABLE :: xc(:,:), nrm(:,:), area(:), psize(:)   ! (3,n), (3,n), (n), (n)
      REAL(8), ALLOCATABLE :: nj(:,:)       ! (6,n), generalized normals as used in the radiation boundary condition
      REAL(8), ALLOCATABLE :: xnode(:,:)    ! (3,nnode)
      INTEGER, ALLOCATABLE :: ncn(:), ncon(:,:), body(:), base(:), image(:)   ! (n), (4,n), (n), (n), (n)
      REAL(8), ALLOCATABLE :: xr(:,:)       ! (3,nbody), body reference points
   END TYPE panel_mesh

   TYPE :: first_order_state
      ! First-order solution for one frequency and heading, in physical quantities with time factor e^{-i omega t}.
      REAL(8) :: omega = 0D0, k = 0D0, nu = 0D0, h = -1D0, amp = 1D0, beta = 0D0, g = 0D0, rho = 0D0   ! h < 0: deep water
      COMPLEX(8), ALLOCATABLE :: phi_i(:), phi_p(:), dphi_p_dn(:), dphi_dn(:)   ! (n)
      COMPLEX(8), ALLOCATABLE :: grad_i(:,:), grad_p(:,:)   ! (3,n); grad_p is filled later by the surface velocity code
      COMPLEX(8), ALLOCATABLE :: xi(:,:)                    ! (6,nbody), translations then rotations
   END TYPE first_order_state

   TYPE :: surface_geom
      ! Panel neighbours, waterline segments and per-body waterplane and volume properties, built by build_surface_geom.
      ! Neighbours are in CSR form: the neighbours of panel i are nbr(nbr_start(i):nbr_start(i+1)-1).
      ! Waterline segments are oriented counter-clockwise around the waterplane (body on the left). Per-body
      ! properties are relative to the body reference point xr: xf, yf the centre of flotation, l11 = int x^2 dA,
      ! l12 = int x y dA, l22 = int y^2 dA over the waterplane, xb, yb, zb the centre of buoyancy.
      INTEGER, ALLOCATABLE :: nbr_start(:), nbr(:)
      INTEGER :: nwl = 0
      REAL(8), ALLOCATABLE :: wl_x1(:,:), wl_x2(:,:), wl_xm(:,:), wl_len(:)   ! (3,nwl), (3,nwl), (3,nwl), (nwl)
      REAL(8), ALLOCATABLE :: wl_nh(:,:), wl_nz(:)    ! unit horizontal normal into the body (3,nwl); n_z of the panel
      INTEGER, ALLOCATABLE :: wl_panel(:), wl_body(:)
      REAL(8), ALLOCATABLE :: awp(:), xf(:), yf(:), l11(:), l12(:), l22(:), vol(:), xb(:), yb(:), zb(:)   ! (nbody)
      INTEGER, ALLOCATABLE :: n_wl_gaps(:)            ! (nbody)
   END TYPE surface_geom

END MODULE SecondOrderTypes
