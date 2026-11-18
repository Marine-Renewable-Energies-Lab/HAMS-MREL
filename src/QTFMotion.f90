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
MODULE QTFMotion
   ! Second-order motion of free bodies (Phase 6c): the added mass and damping at a second-order frequency Omega from the
   ! stored radiation potentials psi, and the matrix of the equation of motion at Omega. On explicit arrays only (no
   ! globals of WavDynMods), for any number of bodies.
   !
   ! Added mass and damping. psi(:, k) of radiation_psi (FirstOrderState) is the unit-velocity radiation potential of
   ! mode k on the full panel set (symmetry images expanded), time factor e^{-i Omega t}, dpsi_k/dn = n_k on the body of
   ! mode k, n out of the fluid (into the body). RFORCE (PotentWavForce.f90) takes
   !    A_lk + i B_lk/Omega = rho int psi_k n_l dS                                                         (1)
   ! over the panels of the body of mode l, l the force and k the motion: AMAS(l, k), BDMP(l, k). radiation_coefficients
   ! is the same sum on the expanded mesh. The image signs RX, RY of RFORCE are in mesh%nj of the image panels already,
   ! so no factor is needed here; the half-mesh test checks it against the whole mesh.
   !
   ! Sign of the damping term. Take the motion xi(t) = Re(xi e^{-i Omega t}). Its velocity is Re(-i Omega xi e^{-i Omega t}),
   ! so the radiation potential is phi = Re(psi_k (-i Omega xi_k) e^{-i Omega t}). The linear pressure p = -rho dphi/dt
   ! has the amplitude -rho (-i Omega)(-i Omega xi_k) psi_k = rho Omega^2 xi_k psi_k. The force on the body is
   ! F_l = int p n_l dS with n into the body, so by (1)
   !    F_l = Omega^2 xi_k rho int psi_k n_l dS = (Omega^2 A_lk + i Omega B_lk) xi_k,
   ! that is F = -A xi'' - B xi' with xi'' = -Omega^2 xi and xi' = -i Omega xi, the radiation force of Newman (1977),
   ! chapter 6. Newton's law -Omega^2 M xi = F - C xi - K_ext xi - B_ext xi' + F_ext, with the hydrostatic restoring C,
   ! an external spring K_ext and linear damper B_ext and the exciting force F_ext, gives
   !    [-Omega^2 (M + A) - i Omega (B + B_ext) + C + K_ext] xi = F_ext                                    (2)
   ! motion_matrix returns the matrix of (2). SolveMotion writes +i omega (BDMP + BLNR) because HAMS's first-order motions
   ! are in e^{+i omega t} (build_state of FirstOrderState conjugates DSPL); here everything is e^{-i Omega t}.
   ! motion_solve_omega solves (2) for all force columns at one Omega and flags an ill-conditioned matrix.
   !
   ! Source: Newman, J.N. (1977). Marine Hydrodynamics. MIT Press. Chapter 6: the radiation force and the equations of
   ! motion, and section 6.15, reciprocity (A and B symmetric).
   USE SecondOrderTypes, ONLY: panel_mesh
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: radiation_coefficients, motion_matrix, motion_solve_omega

CONTAINS

   SUBROUTINE radiation_coefficients(mesh, psi, omega, rho, a, b)
      ! A and B of (1) at omega from psi: a(l, k) = rho Re S_lk, b(l, k) = rho omega Im S_lk, S_lk = sum psi_k n_l dS over
      ! the panels of the body of mode l, as RFORCE (AMAS(l, k), BDMP(l, k)) and write_psi_rows of SecondOrderDriver
      !   psi(n, 6*nbody)  radiation_psi at omega (> 0)
      !   a, b(6*nbody, 6*nbody)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      REAL(8), INTENT(IN) :: omega, rho
      REAL(8), INTENT(OUT) :: a(:, :), b(:, :)
      COMPLEX(8) :: s
      INTEGER :: nmode, k, l, lb, lm, i

      nmode = 6*mesh%nbody
      IF (.NOT. (omega > 0D0)) THEN
         WRITE (*, *) 'radiation_coefficients: omega must be positive, got ', omega
         ERROR STOP 1
      END IF
      IF (SIZE(psi, 1) /= mesh%n .OR. SIZE(psi, 2) /= nmode .OR. ANY(SHAPE(a) /= [nmode, nmode]) .OR. &
          ANY(SHAPE(b) /= [nmode, nmode])) THEN
         WRITE (*, *) 'radiation_coefficients: psi, a or b does not have the size of the mesh'
         ERROR STOP 1
      END IF
      DO k = 1, nmode
         DO l = 1, nmode
            lb = (l - 1)/6 + 1
            lm = l - 6*(lb - 1)
            s = (0D0, 0D0)
            DO i = 1, mesh%n
               IF (mesh%body(i) == lb) s = s + psi(i, k)*mesh%nj(lm, i)*mesh%area(i)
            END DO
            a(l, k) = rho*REAL(s, 8)
            b(l, k) = rho*omega*AIMAG(s)
         END DO
      END DO
   END SUBROUTINE radiation_coefficients

   SUBROUTINE motion_matrix(omega, mass, a, b, b_ext, c, k_ext, lhs)
      ! The matrix of (2): lhs = -omega^2 (mass + a) - i omega (b + b_ext) + c + k_ext, time factor e^{-i omega t}.
      ! All arguments (n, n), n = 6*nbody: mass the body mass matrix, a and b of radiation_coefficients at omega, b_ext
      ! the external linear damping (BLNR), c the hydrostatic restoring, k_ext the external stiffness (KSTF)
      REAL(8), INTENT(IN) :: omega
      REAL(8), INTENT(IN) :: mass(:, :), a(:, :), b(:, :), b_ext(:, :), c(:, :), k_ext(:, :)
      COMPLEX(8), INTENT(OUT) :: lhs(:, :)
      INTEGER :: n

      n = SIZE(lhs, 1)
      IF (SIZE(lhs, 2) /= n .OR. ANY(SHAPE(mass) /= [n, n]) .OR. ANY(SHAPE(a) /= [n, n]) .OR. &
          ANY(SHAPE(b) /= [n, n]) .OR. ANY(SHAPE(b_ext) /= [n, n]) .OR. ANY(SHAPE(c) /= [n, n]) .OR. &
          ANY(SHAPE(k_ext) /= [n, n])) THEN
         WRITE (*, *) 'motion_matrix: the matrices do not all have the size of lhs'
         ERROR STOP 1
      END IF
      lhs = CMPLX(-omega**2*(mass + a) + c + k_ext, -omega*(b + b_ext), KIND=8)
   END SUBROUTINE motion_matrix

   SUBROUTINE motion_solve_omega(lhs, forces, xi, ok, rcond)
      ! The motion xi of (2) at one Omega for all npair force columns: lhs xi = forces, lhs of motion_matrix, forces
      ! the e^{-i Omega t} exciting forces (unscaled). One LU factorisation (ZGETRF) of a copy of lhs, then ZGETRS on all
      ! columns at once. The reciprocal condition number in the 1-norm, rcond = 1/(||lhs||_1 ||lhs^-1||_1), is estimated
      ! by ZGECON from the LU factors and ||lhs||_1 (ZLANGE '1' before factoring), LAPACK Users' Guide, 3rd ed. (Anderson
      ! et al. 1999), section 4.4 (error bounds for linear equation solving), and the ZGECON reference.
      ! ok = rcond >= RCOND_MIN. When not ok, or when ZGETRF finds an exactly zero pivot (info > 0, rcond = 0), xi is
      ! zero (R6): a free body with no restoring in some mode at Omega -> 0 has no bounded motion there, and the caller
      ! flags the frequency instead of writing a meaningless value.
      !   lhs(n, n), forces(n, npair), xi(n, npair)
      COMPLEX(8), INTENT(IN) :: lhs(:, :), forces(:, :)
      COMPLEX(8), INTENT(OUT) :: xi(:, :)
      LOGICAL, INTENT(OUT) :: ok
      REAL(8), INTENT(OUT) :: rcond
      ! Below this the solve loses more than 12 of the 16 digits of the forces
      REAL(8), PARAMETER :: RCOND_MIN = 1D-12
      REAL(8), EXTERNAL :: ZLANGE
      EXTERNAL :: ZGETRF, ZGETRS, ZGECON
      COMPLEX(8), ALLOCATABLE :: lu(:, :), work(:)
      REAL(8), ALLOCATABLE :: rwork(:)
      INTEGER, ALLOCATABLE :: ipiv(:)
      REAL(8) :: anorm
      INTEGER :: n, npair, info

      n = SIZE(lhs, 1)
      npair = SIZE(forces, 2)
      IF (n < 1 .OR. SIZE(lhs, 2) /= n .OR. SIZE(forces, 1) /= n .OR. ANY(SHAPE(xi) /= [n, npair])) THEN
         WRITE (*, *) 'motion_solve_omega: lhs is not square, or forces or xi does not have its size'
         ERROR STOP 1
      END IF
      ALLOCATE (lu(n, n), ipiv(n), work(2*n), rwork(2*n))
      lu = lhs
      anorm = ZLANGE('1', n, n, lu, n, rwork)
      CALL ZGETRF(n, n, lu, n, ipiv, info)
      IF (info < 0) THEN
         WRITE (*, *) 'motion_solve_omega: ZGETRF argument ', -info, ' is illegal'
         ERROR STOP 1
      ELSE IF (info > 0) THEN
         rcond = 0D0
      ELSE
         CALL ZGECON('1', n, lu, n, anorm, rcond, work, rwork, info)
         IF (info /= 0) THEN
            WRITE (*, *) 'motion_solve_omega: ZGECON argument ', -info, ' is illegal'
            ERROR STOP 1
         END IF
      END IF
      ok = rcond >= RCOND_MIN
      IF (.NOT. ok .OR. npair == 0) THEN
         xi = (0D0, 0D0)
         RETURN
      END IF
      xi = forces
      CALL ZGETRS('N', n, npair, lu, n, ipiv, xi, n, info)
      IF (info /= 0) THEN
         WRITE (*, *) 'motion_solve_omega: ZGETRS argument ', -info, ' is illegal'
         ERROR STOP 1
      END IF
   END SUBROUTINE motion_solve_omega

END MODULE QTFMotion
