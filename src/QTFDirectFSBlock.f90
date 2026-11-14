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
SUBMODULE (QTFDirect) QTFDirectFSBlock_sm
   ! The free-surface part of the right sides of a block of pairs of direct_at_omega (Phase 7b, qtf_terms = full).
   ! A submodule of QTFDirect: direct_at_omega calls it, and it calls direct_fs_rhs of QTFDirectFS, which uses QTFDirect, so
   ! as a separate module the two would use each other. QTFDirect holds the interface; direct_fs_data, first_order_state
   ! and pw_set come by host association.
   !
   ! Per pair (header of QTFDirect): the ordered pair (a, b), omega_a > omega_b, is the pair itself when Omega > 0 and
   ! the swapped one when Omega < 0. The right side at |Omega| of a pair with Omega < 0 is that of the conjugated forcing,
   ! and conj Q_D(1, 2) = Q_D(2, 1) [derived in fs_pair_terms of QTFFreeSurface: pair_product of KERNEL_DIFF], so the
   ! swapped pair gives it, Q_D at the inner points and the outer moments alike. The moments are not conjugated: they are
   ! the integrals against the B_mn of G at |Omega|, the Green function of the right side at |Omega|.
   !   q(:, col)       Q_D(a, b) at the inner points (fs_forcing_points), with the g of fs (= Omega^2 / V, the g of the
   !                   Green function; QTFDirectFS header)
   !   mom(:, :, col)  fs_outer_moments of (a, b) with k0, kev, mg and nev of fs, to r_far = direct_fs_rfar of
   !                   fs_outer_rfar (the r_far of the indirect outer integral, with the basis B_mn in place of psi)
   ! then direct_fs_rhs adds r_FS = s_D / g [gfs (w q) + a^T mom] of every column.
   USE QTFDirectFS, ONLY: direct_fs_rhs, direct_fs_rfar
   USE QTFFreeSurface, ONLY: fs_forcing_points, fs_outer_moments, fs_outer_rfar
   IMPLICIT NONE

CONTAINS

   MODULE SUBROUTINE direct_fs_block(fs, store, pairs, dws, rhs, seconds)
      TYPE(direct_fs_data), INTENT(IN) :: fs
      TYPE(first_order_state), INTENT(IN) :: store(:, :)
      INTEGER, INTENT(IN) :: pairs(:, :)
      REAL(8), INTENT(IN) :: dws(:)
      COMPLEX(8), INTENT(INOUT) :: rhs(:, :, :)
      REAL(8), INTENT(INOUT) :: seconds(2)
      COMPLEX(8), ALLOCATABLE :: q(:, :), mom(:, :, :)
      TYPE(pw_set) :: basis(1)
      REAL(8) :: rf
      INTEGER :: nc, col, fa, ba, fb, bb
      INTEGER(8) :: c0, c1, rate

      nc = SIZE(pairs, 2)
      IF (SIZE(pairs, 1) /= 4 .OR. SIZE(dws) /= nc .OR. SIZE(rhs, 2) /= nc) THEN
         WRITE (*, '(a)') ' direct_fs_block: pairs, dws and rhs do not match.'
         ERROR STOP 1
      END IF
      IF (nc == 0) RETURN
      CALL SYSTEM_CLOCK(c0, rate)
      ALLOCATE (q(fs%nfs, nc), mom(-fs%mg:fs%mg, 0:fs%nev, nc))
      ! The basis B_mn as a pw_set of zero amplitudes: its wave numbers enter r_far as those of psi do
      basis(1)%mmax = fs%mg
      basis(1)%nev = fs%nev
      basis(1)%k0 = fs%k0
      basis(1)%h = fs%h
      basis(1)%radius = fs%radius
      basis(1)%centre = fs%centre
      ALLOCATE (basis(1)%kev(fs%nev), basis(1)%a(-fs%mg:fs%mg, 0:fs%nev))
      basis(1)%kev = fs%kev(1:fs%nev)
      basis(1)%a = (0D0, 0D0)
      DO col = 1, nc
         IF (dws(col) > 0D0) THEN
            fa = pairs(1, col)
            ba = pairs(2, col)
            fb = pairs(3, col)
            bb = pairs(4, col)
         ELSE
            fa = pairs(3, col)
            ba = pairs(4, col)
            fb = pairs(1, col)
            bb = pairs(2, col)
         END IF
         ASSOCIATE (sa => store(fa, ba), sb => store(fb, bb))
            CALL fs_forcing_points(sa%omega, fs%f(:, fa, ba), fs%inc(:, fa, ba), sb%omega, fs%f(:, fb, bb), &
                                   fs%inc(:, fb, bb), fs%g, q(:, col))
            rf = direct_fs_rfar(fs, fs_outer_rfar(fs%pw(fa, ba), sa%k, fs%pw(fb, bb), sb%k, basis))
            CALL fs_outer_moments(sa%omega, fs%pw(fa, ba), sa%beta, sa%k, sb%omega, fs%pw(fb, bb), sb%beta, sb%k, fs%k0, &
                                  fs%kev, fs%h, fs%centre, fs%radius, fs%mg, fs%nev, sa%omega - sb%omega, fs%g, fs%np, rf, &
                                  mom(:, :, col))
         END ASSOCIATE
      END DO
      CALL SYSTEM_CLOCK(c1)
      seconds(1) = seconds(1) + REAL(c1 - c0, 8)/REAL(rate, 8)
      CALL direct_fs_rhs(fs, q, mom, rhs)
      CALL SYSTEM_CLOCK(c0)
      seconds(2) = seconds(2) + REAL(c0 - c1, 8)/REAL(rate, 8)
   END SUBROUTINE direct_fs_block

END SUBMODULE QTFDirectFSBlock_sm
