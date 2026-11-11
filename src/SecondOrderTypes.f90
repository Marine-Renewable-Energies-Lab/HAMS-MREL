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

   PUBLIC :: so_options, panel_mesh, first_order_state, surface_geom, cs_points, cs_fields
   PUBLIC :: CS_WALL, CS_LINE, CS_FREE, QTF_TERMS_QUADRATIC, QTF_TERMS_BODY, QTF_B2_BODY, QTF_B2_CONTROL

   INTEGER, PARAMETER :: CS_WALL = 1, CS_LINE = 2, CS_FREE = 3   ! C panel, Gamma_c segment midpoint, F panel
   ! qtf_terms: the quadratic part only, or also the potential part with the body forcing (indirect method)
   INTEGER, PARAMETER :: QTF_TERMS_QUADRATIC = 0, QTF_TERMS_BODY = 1
   ! qtf_b2: the surface piece of B2 integrated on the hull, or moved to a control surface C and the free surface F
   INTEGER, PARAMETER :: QTF_B2_BODY = 1, QTF_B2_CONTROL = 2

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
      ! Control-surface method (drift_methods = control)
      LOGICAL :: do_control = .FALSE.
      INTEGER :: cs_velocity = 1         ! 1 green, 2 source
      INTEGER :: cs_shape = 1            ! 1 cylinder, 2 box, 3 file
      REAL(8) :: cs_margin = 0D0, cs_bottom_margin = 0D0, cs_panel_size = 0D0   ! m; 0 sets the default
      REAL(8) :: cs_fd_step = 0.02D0     ! fraction of the local control-surface panel size
      INTEGER :: cs_fs_rows = 0          ! free-surface rows; 0 sets the default
      ! Quadratic part of the QTFs (near-field method)
      INTEGER :: qtf = 0                 ! 0 off, 1 difference frequency, 2 sum frequency, 3 both
      REAL(8) :: qtf_max_dw = 0D0        ! rad/s, largest |omega_i - omega_j| of QTF_Quadratic.10d; 0: all pairs
      INTEGER :: qtf_terms = QTF_TERMS_QUADRATIC
      INTEGER :: qtf_b2 = QTF_B2_CONTROL ! used with qtf_terms = body only
      ! Panel weights of the B2 hull term at sharp edges (SharpEdge), with qtf_b2 = control; edges where the normals differ
      ! by more than qtf_edge_angle degrees. Off by default: on the sharp-cylinder fixture they moved B2 F_x and F_z away
      ! from the exact reference (15.1% to 17.2%, 30.0% to 37.4%) for an M_y gain of 11.9% to 11.1%
      LOGICAL :: qtf_edge_weights = .FALSE.
      REAL(8) :: qtf_edge_angle = 45D0
      ! The keys were in SecondOrder.in, so that the log can say when they have no effect
      LOGICAL :: qtf_b2_given = .FALSE., qtf_edge_given = .FALSE.
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
      ! Neighbours along the waterline loop: wl_next starts at the end node of the segment, wl_prev ends at its start
      ! node (nwl); 0 at a gap or where more than two segments meet.
      INTEGER, ALLOCATABLE :: wl_prev(:), wl_next(:)
      REAL(8), ALLOCATABLE :: awp(:), xf(:), yf(:), l11(:), l12(:), l22(:), vol(:), xb(:), yb(:), zb(:)   ! (nbody)
      INTEGER, ALLOCATABLE :: n_wl_gaps(:)            ! (nbody)
   END TYPE surface_geom

   TYPE :: cs_points
      ! Evaluation points of the control-surface method, all bodies together. w is the panel area or segment length.
      ! On CS_LINE points nrm holds n-bar (horizontal, into F); on CS_FREE points it holds (0, 0, -1).
      INTEGER :: n = 0
      REAL(8), ALLOCATABLE :: x(:,:), nrm(:,:), w(:), psize(:)   ! (3,n), (3,n), (n), (n)
      INTEGER, ALLOCATABLE :: kind(:), body(:)                   ! (n), (n)
   END TYPE cs_points

   TYPE :: cs_fields
      ! Total first-order potential and its gradient at the points of a cs_points set, for one state. Filled by fields_green
      ! (FieldVelocity). phi is zero at CS_WALL points below the stencil reach of z = 0, where it is not computed.
      COMPLEX(8), ALLOCATABLE :: phi(:), grad(:,:)   ! (n), (3,n)
   END TYPE cs_fields

END MODULE SecondOrderTypes
