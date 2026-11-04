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
MODULE FirstOrderState
   ! Converts HAMS's first-order solution into physical quantities for the second-order loads, with time factor
   ! e^{-i omega t}. Conventions: docs/superpowers/specs/2026-11-04-second-order-waves-design.md, section 5.2.
   !
   !   incident potential   phi_I = -i (g amp / omega) Z(z) e^{i k (x cos(beta) + y sin(beta))}
   !   HAMS MXPOT(:,6N+1,:) holds the scattered potential divided by i, MXPOT(:,j,:) the unit-velocity radiation
   !   potential psi_j (no factor of i), and DSPL holds the motions in e^{+i omega t}, so xi_j = conj(DSPL(j)).

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: depth_factor, wave_gravity, incident_phi, incident_grad, expand_images, build_state

CONTAINS

   PURE SUBROUTINE depth_factor(k, z, h, zf, dzf)
      ! Vertical profile Z(z) and its derivative dZ/dz. Deep water (h < 0): Z = e^{kz}. Finite depth:
      ! Z = cosh k(z+h) / cosh kh, written with decaying exponentials so that it cannot overflow for large kh.
      REAL(8), INTENT(IN) :: k, z, h
      REAL(8), INTENT(OUT) :: zf, dzf
      REAL(8) :: ea, eb, den

      IF (h < 0D0) THEN
         zf = EXP(k*z)
         dzf = k*zf
      ELSE
         ea = EXP(k*z)
         eb = EXP(-k*(z + 2D0*h))
         den = 1D0 + EXP(-2D0*k*h)
         zf = (ea + eb)/den
         dzf = k*(ea - eb)/den
      END IF
   END SUBROUTINE depth_factor

   PURE REAL(8) FUNCTION wave_gravity(omega, k, h)
      ! The gravity of the wave the solver used. VINP and DINP in finite depth have the factor omega/(k sinh kh), and k
      ! comes from DISPERSION, which has its own g = 9.807, so their wave has omega^2/(k tanh kh), not the g of Const_mod.
      REAL(8), INTENT(IN) :: omega, k, h

      IF (h < 0D0) THEN
         wave_gravity = omega*omega/k
      ELSE
         wave_gravity = omega*omega/(k*TANH(k*h))
      END IF
   END FUNCTION wave_gravity

   PURE COMPLEX(8) FUNCTION incident_phi(x, omega, k, h, amp, beta, g)
      ! Incident wave potential at x(3), for an elevation amp e^{i k (x cos(beta) + y sin(beta))} at the origin.
      REAL(8), INTENT(IN) :: x(3), omega, k, h, amp, beta, g
      REAL(8) :: zf, dzf

      CALL depth_factor(k, x(3), h, zf, dzf)
      incident_phi = CMPLX(0D0, -g*amp/omega*zf, KIND=8) &
                     *EXP(CMPLX(0D0, k*(x(1)*COS(beta) + x(2)*SIN(beta)), KIND=8))
   END FUNCTION incident_phi

   PURE SUBROUTINE incident_grad(x, omega, k, h, amp, beta, g, grad)
      ! Gradient of incident_phi at x(3).
      REAL(8), INTENT(IN) :: x(3), omega, k, h, amp, beta, g
      COMPLEX(8), INTENT(OUT) :: grad(3)
      REAL(8) :: zf, dzf
      COMPLEX(8) :: pref

      CALL depth_factor(k, x(3), h, zf, dzf)
      pref = CMPLX(0D0, -g*amp/omega, KIND=8) &
             *EXP(CMPLX(0D0, k*(x(1)*COS(beta) + x(2)*SIN(beta)), KIND=8))
      grad(1) = CMPLX(0D0, k*COS(beta), KIND=8)*pref*zf
      grad(2) = CMPLX(0D0, k*SIN(beta), KIND=8)*pref*zf
      grad(3) = pref*dzf
   END SUBROUTINE incident_grad

   SUBROUTINE expand_images(xc_b, dxyz_b, area_b, psize_b, xnode_b, ncn_b, ncon_b, body_b, xr, isx, isy, nsys, mesh)
      ! Builds the full panel set from the base panels by adding the symmetry image. Base arrays use the layout of the
      ! HAMS globals (panel index first), the result uses the layout of panel_mesh.
      !   xc_b(nbase,3)      panel centroids              dxyz_b(nbase,6)  normal and moment-arm components, as DXYZ
      !   area_b, psize_b    (nbase)                      xnode_b(nnode_b,3), ncn_b(nbase), ncon_b(nbase,4)
      !   body_b(nbase)      body of each panel           xr(3,nbody)      body reference points
      !   isx, isy, nsys     symmetry flags of HAMS (nsys = 1: no image, 2: one image)
      ! Image rule, as in EFORCE_MULTI: with ISX = 1 and ISY = 0 the image is x -> -x, otherwise y -> -y. The image
      ! panels follow the base panels, and image nodes follow the base nodes, so ncon stays valid. Components of nj
      ! take the factors (fx, fy, 1, fy, fx, fx*fy), which are the RX/RY signs of the force integration.
      REAL(8), INTENT(IN) :: xc_b(:,:), dxyz_b(:,:), area_b(:), psize_b(:), xnode_b(:,:), xr(:,:)
      INTEGER, INTENT(IN) :: ncn_b(:), ncon_b(:,:), body_b(:), isx, isy, nsys
      TYPE(panel_mesh), INTENT(OUT) :: mesh
      INTEGER :: nbase, nnode_b, is, off, node_off, i
      REAL(8) :: fx, fy, f3(3), fj(6)

      IF (nsys /= 1 .AND. nsys /= 2) THEN
         WRITE (*, *) 'expand_images: nsys must be 1 or 2, got ', nsys
         ERROR STOP 1
      END IF
      nbase = SIZE(xc_b, 1)
      nnode_b = SIZE(xnode_b, 1)
      mesh%nbase = nbase
      mesh%nsys = nsys
      mesh%nbody = SIZE(xr, 2)
      mesh%n = nbase*nsys
      mesh%nnode = nnode_b*nsys
      ALLOCATE (mesh%xc(3, mesh%n), mesh%nrm(3, mesh%n), mesh%area(mesh%n), mesh%psize(mesh%n), mesh%nj(6, mesh%n))
      ALLOCATE (mesh%xnode(3, mesh%nnode), mesh%ncn(mesh%n), mesh%ncon(4, mesh%n))
      ALLOCATE (mesh%body(mesh%n), mesh%base(mesh%n), mesh%image(mesh%n), mesh%xr(3, mesh%nbody))
      mesh%xr = xr

      DO is = 1, nsys
         fx = 1D0
         fy = 1D0
         IF (is == 2) THEN
            IF (isx == 1 .AND. isy == 0) THEN
               fx = -1D0
            ELSE
               fy = -1D0
            END IF
         END IF
         f3 = [fx, fy, 1D0]
         fj = [fx, fy, 1D0, fy, fx, fx*fy]
         off = (is - 1)*nbase
         node_off = (is - 1)*nnode_b
         DO i = 1, nbase
            mesh%xc(:, off + i) = f3*xc_b(i, :)
            mesh%nrm(:, off + i) = f3*dxyz_b(i, 1:3)
            mesh%nj(:, off + i) = fj*dxyz_b(i, 1:6)
            mesh%area(off + i) = area_b(i)
            mesh%psize(off + i) = psize_b(i)
            mesh%ncn(off + i) = ncn_b(i)
            mesh%ncon(:, off + i) = ncon_b(i, :) + node_off
            mesh%body(off + i) = body_b(i)
            mesh%base(off + i) = i
            mesh%image(off + i) = is
         END DO
         DO i = 1, nnode_b
            mesh%xnode(:, node_off + i) = f3*xnode_b(i, :)
         END DO
      END DO
   END SUBROUTINE expand_images

   SUBROUTINE build_state(mesh, omega, k, nu, h, amp, beta, g, rho, mxpot, dspl, motion_fixed, state)
      ! Converts the HAMS solution for one frequency and heading into the physical first-order state.
      !   mxpot(nbase, 6*nbody+1, nsys)  HAMS MXPOT: modes 1..6*nbody, then the scattered potential over i
      !   dspl(6*nbody)                  HAMS motions in e^{+i omega t}
      !   motion_fixed                   true: body held fixed, xi = 0
      ! phi_p is the perturbation potential (scattered plus radiated). grad_p is allocated here and filled later.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: omega, k, nu, h, amp, beta, g, rho
      COMPLEX(8), INTENT(IN) :: mxpot(:,:,:), dspl(:)
      LOGICAL, INTENT(IN) :: motion_fixed
      TYPE(first_order_state), INTENT(OUT) :: state
      COMPLEX(8), PARAMETER :: ci = (0D0, 1D0)
      INTEGER :: nmode, md, i, b, im, j, jb
      COMPLEX(8) :: acc
      REAL(8) :: g_wave

      nmode = 6*mesh%nbody
      md = nmode + 1
      state%omega = omega
      state%k = k
      state%nu = nu
      state%h = h
      state%amp = amp
      state%beta = beta
      state%g = g
      state%rho = rho
      ! The incident field must be the one the solver used. state%g stays the given g, for the rho g terms.
      g_wave = wave_gravity(omega, k, h)
      ALLOCATE (state%phi_i(mesh%n), state%phi_p(mesh%n), state%dphi_p_dn(mesh%n), state%dphi_dn(mesh%n))
      ALLOCATE (state%grad_i(3, mesh%n), state%grad_p(3, mesh%n), state%xi(6, mesh%nbody))
      state%grad_p = (0D0, 0D0)

      IF (motion_fixed) THEN
         state%xi = (0D0, 0D0)
      ELSE
         state%xi = RESHAPE(CONJG(dspl(1:nmode)), [6, mesh%nbody])
      END IF

      DO i = 1, mesh%n
         b = mesh%base(i)
         im = mesh%image(i)
         state%phi_i(i) = incident_phi(mesh%xc(:, i), omega, k, h, amp, beta, g_wave)
         CALL incident_grad(mesh%xc(:, i), omega, k, h, amp, beta, g_wave, state%grad_i(:, i))

         ! Radiation modes of every body act on every panel
         acc = ci*mxpot(b, md, im)
         DO jb = 1, mesh%nbody
            DO j = 1, 6
               acc = acc - ci*omega*state%xi(j, jb)*mxpot(b, 6*(jb - 1) + j, im)
            END DO
         END DO
         state%phi_p(i) = acc

         ! Normal velocity of the panel's own body only
         jb = mesh%body(i)
         acc = (0D0, 0D0)
         DO j = 1, 6
            acc = acc - ci*omega*state%xi(j, jb)*mesh%nj(j, i)
         END DO
         state%dphi_dn(i) = acc
         state%dphi_p_dn(i) = acc - SUM(state%grad_i(:, i)*mesh%nrm(:, i))
      END DO
   END SUBROUTINE build_state

END MODULE FirstOrderState
