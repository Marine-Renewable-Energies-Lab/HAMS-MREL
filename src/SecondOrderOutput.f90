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
MODULE SecondOrderOutput
   ! Writers for the second-order result files. The loads are computed with time factor e^{-i omega t}. This module is
   ! the one place where they are conjugated, to the e^{+i omega t} convention of the WAMIT-format outputs.

   USE Const_mod, ONLY: G, RHO
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: write_drift_row, write_terms_row, write_qtf_row, write_qtf_terms_row, write_qtf_terms_header, write_field_row
   PUBLIC :: write_direct_terms_header, write_direct_sum_terms_header, write_potential_sum_terms_header
   PUBLIC :: write_motion_row, write_motion_header

CONTAINS

   SUBROUTINE write_drift_row(unit, oufr, beta1_deg, beta2_deg, index, value, refl, moment)
      ! Writes one row PER BETA1 BETA2 I Mod Pha Re Im, with the phase in degrees.
      !   oufr     output frequency in the form chosen by OUFT in ControlFile.in
      !   value    load per unit amplitude pair, e^{-i omega t}
      !   moment   true: scale by rho g refl^2, false: scale by rho g refl
      INTEGER, INTENT(IN) :: unit, index
      REAL(8), INTENT(IN) :: oufr, beta1_deg, beta2_deg, refl
      COMPLEX(8), INTENT(IN) :: value
      LOGICAL, INTENT(IN) :: moment
      REAL(8), PARAMETER :: RAD_TO_DEG = 180D0/3.141592653589793D0
      COMPLEX(8) :: scaled

      scaled = output_value(value, refl, moment)
      WRITE (unit, '(3ES14.6,I6,4ES14.6)') oufr, beta1_deg, beta2_deg, index, ABS(scaled), &
         ATAN2(AIMAG(scaled), REAL(scaled, 8))*RAD_TO_DEG, REAL(scaled, 8), AIMAG(scaled)
   END SUBROUTINE write_drift_row

   SUBROUTINE write_terms_row(unit, oufr, beta1_deg, beta2_deg, body, index, term, value, refl, moment)
      ! Writes one row PER BETA1 BETA2 body I term Re Im of MeanDriftTerms.txt. I = 1..6 and the moments are about XR of
      ! the body; term is the term index of NearFieldDrift. The other arguments are as for write_drift_row.
      INTEGER, INTENT(IN) :: unit, body, index, term
      REAL(8), INTENT(IN) :: oufr, beta1_deg, beta2_deg, refl
      COMPLEX(8), INTENT(IN) :: value
      LOGICAL, INTENT(IN) :: moment
      COMPLEX(8) :: scaled

      scaled = output_value(value, refl, moment)
      WRITE (unit, '(3ES14.6,3I6,2ES14.6)') oufr, beta1_deg, beta2_deg, body, index, term, REAL(scaled, 8), AIMAG(scaled)
   END SUBROUTINE write_terms_row

   SUBROUTINE write_qtf_row(unit, oufr1, oufr2, beta1_deg, beta2_deg, index, value, refl, moment)
      ! Writes one row PER1 PER2 BETA1 BETA2 I Mod Pha Re Im of QTF_Quadratic.10d or .10s, scaled and conjugated as
      ! write_drift_row. oufr1, beta1 belong to the first component of the pair, oufr2, beta2 to the second.
      INTEGER, INTENT(IN) :: unit, index
      REAL(8), INTENT(IN) :: oufr1, oufr2, beta1_deg, beta2_deg, refl
      COMPLEX(8), INTENT(IN) :: value
      LOGICAL, INTENT(IN) :: moment
      REAL(8), PARAMETER :: RAD_TO_DEG = 180D0/3.141592653589793D0
      COMPLEX(8) :: scaled

      scaled = output_value(value, refl, moment)
      WRITE (unit, '(4ES14.6,I6,4ES14.6)') oufr1, oufr2, beta1_deg, beta2_deg, index, ABS(scaled), &
         ATAN2(AIMAG(scaled), REAL(scaled, 8))*RAD_TO_DEG, REAL(scaled, 8), AIMAG(scaled)
   END SUBROUTINE write_qtf_row

   SUBROUTINE write_motion_row(unit, oufr1, oufr2, beta1_deg, beta2_deg, index, xi, refl, rotation)
      ! Writes one row PER1 PER2 BETA1 BETA2 I Mod Pha Re Im of QTF_Motion.11d/.11s/.12d/.12s (Phase 6c, ruling R4): the
      ! second-order motion xi per unit amplitude pair, e^{-i Omega t}, conjugated as write_qtf_row and scaled as WAMIT
      ! scales motions. The force files hold F/(rho g A^2 L^k), k = 1 for a force and 2 for a moment. Dividing by a
      ! stiffness of the size rho g L^2 (translations) or rho g L^4 (rotations) gives a translation of the size A^2/L and
      ! a rotation of the size A^2/L^2 [derived], so with A = 1 the row holds xi L (translation) or xi L^2 (rotation,
      ! radians). oufr1, beta1 belong to the first component of the pair, oufr2, beta2 to the second.
      INTEGER, INTENT(IN) :: unit, index
      REAL(8), INTENT(IN) :: oufr1, oufr2, beta1_deg, beta2_deg, refl
      COMPLEX(8), INTENT(IN) :: xi
      LOGICAL, INTENT(IN) :: rotation
      REAL(8), PARAMETER :: RAD_TO_DEG = 180D0/3.141592653589793D0
      COMPLEX(8) :: scaled

      IF (rotation) THEN
         scaled = CONJG(xi)*refl**2
      ELSE
         scaled = CONJG(xi)*refl
      END IF
      WRITE (unit, '(4ES14.6,I6,4ES14.6)') oufr1, oufr2, beta1_deg, beta2_deg, index, ABS(scaled), &
         ATAN2(AIMAG(scaled), REAL(scaled, 8))*RAD_TO_DEG, REAL(scaled, 8), AIMAG(scaled)
   END SUBROUTINE write_motion_row

   SUBROUTINE write_motion_header(unit, force_file, route, left_out)
      ! The comment lines of a QTF_Motion file (Phase 6c): the force file it follows and its route (indirect or direct),
      ! the columns and scaling, the matrix of the equation of motion, and one line per Omega whose pairs are left out
      ! (left_out, from QTFDriver: Omega, rcond and the modes without restoring), or "none"
      INTEGER, INTENT(IN) :: unit
      CHARACTER(LEN=*), INTENT(IN) :: force_file, route, left_out(:)
      INTEGER :: m

      WRITE (unit, '(a)') '# Second-order motion of the free bodies (qtf_motion = 1) from the force of '//force_file// &
         ' ('//TRIM(route)//' route): its rows and order, less those left out below.'
      WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I Mod Pha Re Im; translations over A^2/L, rotations over A^2/L^2 (rad),'// &
         ' conjugated as QTF_Quadratic.10d; I = 6 (body - 1) + mode, about the reference point of the body.'
      WRITE (unit, '(a)') '# Solved per Omega: [-Omega^2 (M + A) - i Omega (B + BLNR) + C + KSTF] xi = F, A and B at |Omega|'// &
         ' from psi; BQDR is not used. At Omega = 0: (C + KSTF) xi = F.'
      IF (SIZE(left_out) == 0) THEN
         WRITE (unit, '(a)') '# Left out: none (every Omega has rcond >= 1e-12).'
      END IF
      DO m = 1, SIZE(left_out)
         WRITE (unit, '(a)') '# Left out: the pairs of '//TRIM(left_out(m))//'.'
      END DO
   END SUBROUTINE write_motion_header

   SUBROUTINE write_qtf_terms_header(unit, full)
      ! The four comment lines of QTF_PotentialTerms.txt: the 5 body-forcing terms (qtf_terms = body), or with full
      ! (qtf_terms = full) the 7 terms, FS_inner and FS_outer after B3
      INTEGER, INTENT(IN) :: unit
      LOGICAL, INTENT(IN) :: full

      IF (full) THEN
         WRITE (unit, '(a)') '# Terms of the potential part of QTF_Indirect.11d (qtf_terms = full), scaled and conjugated'// &
            ' as QTF_Quadratic.10d.'
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each term: 1 Froude-Krylov of phi_I2, 2 incident'// &
            ' diffraction, 3 B1, 4 B2, 5 B3, 6 FS_inner, 7 FS_outer.'
         WRITE (unit, '(a)') '# .11d - .10d is their sum. FS_inner: the free-surface integral from the waterlines to the'// &
            ' circle R (inner mesh); FS_outer: beyond the inner mesh (partial waves). One row per pair with omega_i /='// &
            ' omega_j, by ascending |omega_i - omega_j|.'
      ELSE
         WRITE (unit, '(a)') '# Body-forcing terms of the potential part of QTF_Indirect.11d (qtf_terms = body), scaled'// &
            ' and conjugated as QTF_Quadratic.10d.'
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each term: 1 Froude-Krylov of phi_I2, 2 incident'// &
            ' diffraction, 3 B1, 4 B2, 5 B3.'
         WRITE (unit, '(a)') '# .11d - .10d is their sum; the free-surface integral is omitted. One row per pair with'// &
            ' omega_i /= omega_j, by ascending |omega_i - omega_j|.'
      END IF
      WRITE (unit, '(a)') '# The rows are not in the order of .11d: match them by PER1 PER2 BETA1 BETA2 I.'
   END SUBROUTINE write_qtf_terms_header

   SUBROUTINE write_direct_terms_header(unit, full)
      ! The four comment lines of QTF_DirectTerms.txt: the parts D_I2 and D_D of the direct load (qtf_terms = body), or with
      ! full (qtf_terms = full) D_I2, D_D and D_FS, the part from phi_D split into that of the body forcing and that of the
      ! free-surface forcing (Phase 7b)
      INTEGER, INTENT(IN) :: unit
      LOGICAL, INTENT(IN) :: full

      WRITE (unit, '(a)') '# Parts of the direct load of QTF_Direct.12d (qtf_direct = 1), scaled and conjugated as'// &
         ' QTF_Quadratic.10d.'
      IF (full) THEN
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each part: 1 D_I2, the second-order incident'// &
            ' potential phi_I2; 2 D_D, the second-order diffracted and radiated potential phi_D of the body forcing;'// &
            ' 3 D_FS, phi_D of the free-surface forcing (qtf_terms = full).'
      ELSE
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each part: 1 the second-order incident potential'// &
            ' phi_I2, 2 the second-order diffracted and radiated potential phi_D.'
      END IF
      WRITE (unit, '(a)') '# .12d - .10d is their sum. One row per pair with omega_i /= omega_j, by ascending'// &
         ' |omega_i - omega_j|.'
      WRITE (unit, '(a)') '# The rows are not in the order of .12d: match them by PER1 PER2 BETA1 BETA2 I.'
   END SUBROUTINE write_direct_terms_header

   SUBROUTINE write_direct_sum_terms_header(unit, full)
      ! The four comment lines of QTF_DirectTermsSum.txt (Phase 7c-2, ruling R1 of its plan): the parts of the direct sum
      ! load of QTF_Direct.12s beside the indirect sum load of the same run, from the same psi at Omega = omega_i +
      ! omega_j, so that the file checks .12s pair by pair. Columns: the 5 keys, then Re Im of D_I2, D_D and P_BODY (11 in
      ! all, qtf_terms = body), or with full of D_I2, D_D, D_FS, P_BODY and P_FS (15 in all, qtf_terms = full)
      INTEGER, INTENT(IN) :: unit
      LOGICAL, INTENT(IN) :: full

      WRITE (unit, '(a)') '# Parts of the direct sum-frequency load of QTF_Direct.12s (qtf_direct = 1) and the indirect'// &
         ' load of the same pairs, scaled and conjugated as QTF_Quadratic.10s.'
      IF (full) THEN
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of: 1 D_I2, the second-order incident potential'// &
            ' phi_I2; 2 D_D, phi_D of the body forcing; 3 D_FS, phi_D of the free-surface forcing; 4 P_BODY, the five'// &
            ' body-forcing terms of the indirect method; 5 P_FS, its free-surface integral (FS_inner + FS_outer).'
         WRITE (unit, '(a)') '# .12s - .10s = D_I2 + D_D + D_FS; the indirect method gives P_BODY + P_FS for the same'// &
            ' pair. One row per pair of .12s (the diagonal included), by ascending omega_i + omega_j.'
      ELSE
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of: 1 D_I2, the second-order incident potential'// &
            ' phi_I2; 2 D_D, phi_D of the body forcing; 3 P_BODY, the five body-forcing terms of the indirect method.'
         WRITE (unit, '(a)') '# .12s - .10s = D_I2 + D_D; the indirect method gives P_BODY for the same pair (the'// &
            ' free-surface integral is omitted). One row per pair of .12s (the diagonal included), by ascending'// &
            ' omega_i + omega_j.'
      END IF
      WRITE (unit, '(a)') '# The rows are not in the order of .12s: match them by PER1 PER2 BETA1 BETA2 I.'
   END SUBROUTINE write_direct_sum_terms_header

   SUBROUTINE write_potential_sum_terms_header(unit, full)
      ! The four comment lines of QTF_PotentialTermsSum.txt (Phase 6b): the terms of T+_pot, the potential part of
      ! QTF_Indirect.11s, as the header of QTF_PotentialTerms.txt gives those of QTF_Indirect.11d. Columns: the 5 keys, then
      ! Re Im of the 5 body-forcing terms (qtf_terms = body, 15 in all), or with full (qtf_terms = full) of 7 terms, FS_inner
      ! and FS_outer after B3 (19 in all)
      INTEGER, INTENT(IN) :: unit
      LOGICAL, INTENT(IN) :: full

      IF (full) THEN
         WRITE (unit, '(a)') '# Terms of the potential part of QTF_Indirect.11s (qtf_terms = full), scaled and conjugated'// &
            ' as QTF_Quadratic.10s.'
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each term: 1 Froude-Krylov of phi_I2, 2 incident'// &
            ' diffraction, 3 B1, 4 B2, 5 B3, 6 FS_inner, 7 FS_outer.'
         WRITE (unit, '(a)') '# .11s - .10s is their sum. FS_inner: the free-surface integral from the waterlines to the'// &
            ' circle R (inner mesh); FS_outer: beyond the inner mesh (partial waves). One row per pair of .11s (the'// &
            ' diagonal included), by ascending omega_i + omega_j.'
      ELSE
         WRITE (unit, '(a)') '# Body-forcing terms of the potential part of QTF_Indirect.11s (qtf_terms = body), scaled'// &
            ' and conjugated as QTF_Quadratic.10s.'
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each term: 1 Froude-Krylov of phi_I2, 2 incident'// &
            ' diffraction, 3 B1, 4 B2, 5 B3.'
         WRITE (unit, '(a)') '# .11s - .10s is their sum; the free-surface integral is omitted. One row per pair of .11s'// &
            ' (the diagonal included), by ascending omega_i + omega_j.'
      END IF
      WRITE (unit, '(a)') '# The rows are not in the order of .11s: match them by PER1 PER2 BETA1 BETA2 I.'
   END SUBROUTINE write_potential_sum_terms_header

   SUBROUTINE write_qtf_terms_row(unit, oufr1, oufr2, beta1_deg, beta2_deg, index, values, refl, moment)
      ! Writes one row PER1 PER2 BETA1 BETA2 I, then Re Im of each of values, of QTF_PotentialTerms.txt. The keys are as in
      ! write_qtf_row; each value is scaled and conjugated as there, and written with 16 digits, so that the sum of the
      ! terms can be checked against QTF_Indirect.11d - QTF_Quadratic.10d. QTFDriver writes the rows in ascending
      ! |omega_1 - omega_2|, not in the order of .11d: match them by their keys.
      INTEGER, INTENT(IN) :: unit, index
      REAL(8), INTENT(IN) :: oufr1, oufr2, beta1_deg, beta2_deg, refl
      COMPLEX(8), INTENT(IN) :: values(:)
      LOGICAL, INTENT(IN) :: moment
      COMPLEX(8) :: scaled(SIZE(values))
      INTEGER :: t

      DO t = 1, SIZE(values)
         scaled(t) = output_value(values(t), refl, moment)
      END DO
      WRITE (unit, '(4ES14.6,I6,*(ES24.15))') oufr1, oufr2, beta1_deg, beta2_deg, index, &
         (REAL(scaled(t), 8), AIMAG(scaled(t)), t = 1, SIZE(values))
   END SUBROUTINE write_qtf_terms_row

   SUBROUTINE write_field_row(unit, oufr1, oufr2, beta1_deg, beta2_deg, kernel, ipt, surface, quad, pot, refl, direct, &
                              has_direct)
      ! Writes one row PER1 PER2 BETA1 BETA2 K IPT TYPE Mod Pha Re Im Re_quad Im_quad Re_inc Im_inc of QTF_FieldPartial.txt.
      ! quad and pot are the quadratic part and the second-order incident part, per unit amplitude pair, e^{-i omega t}, of
      ! the pressure (surface false, TYPE 1) or the elevation (surface true, TYPE 2). The row holds their sum and each part,
      ! conjugated, the pressure scaled by refl/(rho g) and the elevation by refl, so that they are pressure/(rho g A^2/L)
      ! and elevation/(A^2/L). kernel is 1 for the difference frequency and 2 for the sum frequency; oufr1, beta1 belong to
      ! the first component of the pair, oufr2, beta2 to the second.
      ! With direct (QTF_Field.txt), the row goes on with Re_dir Im_dir has_direct after those 15 columns, which are the
      ! columns of the row without direct to the byte: direct is the direct part (phi_D of the direct solve, per unit
      ! amplitude pair, e^{-i omega t}, as quad and pot), conjugated and scaled as they are, and not added to the sum
      ! columns. has_direct (default .TRUE.) = .FALSE. marks a row the direct solve does not cover (mean row, sum row): it
      ! gets 0 0 0 whatever direct holds.
      INTEGER, INTENT(IN) :: unit, kernel, ipt
      REAL(8), INTENT(IN) :: oufr1, oufr2, beta1_deg, beta2_deg, refl
      LOGICAL, INTENT(IN) :: surface
      COMPLEX(8), INTENT(IN) :: quad, pot
      COMPLEX(8), INTENT(IN), OPTIONAL :: direct
      LOGICAL, INTENT(IN), OPTIONAL :: has_direct
      REAL(8), PARAMETER :: RAD_TO_DEG = 180D0/3.141592653589793D0
      COMPLEX(8) :: total, q, p, d
      LOGICAL :: have

      q = field_value(quad, refl, surface)
      p = field_value(pot, refl, surface)
      total = q + p
      IF (.NOT. PRESENT(direct)) THEN
         WRITE (unit, '(4ES14.6,3I6,8ES14.6)') oufr1, oufr2, beta1_deg, beta2_deg, kernel, ipt, MERGE(2, 1, surface), &
            ABS(total), ATAN2(AIMAG(total), REAL(total, 8))*RAD_TO_DEG, REAL(total, 8), AIMAG(total), REAL(q, 8), AIMAG(q), &
            REAL(p, 8), AIMAG(p)
         RETURN
      END IF
      have = .TRUE.
      IF (PRESENT(has_direct)) have = has_direct
      d = (0D0, 0D0)
      IF (have) d = field_value(direct, refl, surface)
      WRITE (unit, '(4ES14.6,3I6,10ES14.6,I6)') oufr1, oufr2, beta1_deg, beta2_deg, kernel, ipt, MERGE(2, 1, surface), &
         ABS(total), ATAN2(AIMAG(total), REAL(total, 8))*RAD_TO_DEG, REAL(total, 8), AIMAG(total), REAL(q, 8), AIMAG(q), &
         REAL(p, 8), AIMAG(p), REAL(d, 8), AIMAG(d), MERGE(1, 0, have)
   END SUBROUTINE write_field_row

   PURE COMPLEX(8) FUNCTION field_value(value, refl, surface)
      COMPLEX(8), INTENT(IN) :: value
      REAL(8), INTENT(IN) :: refl
      LOGICAL, INTENT(IN) :: surface

      IF (surface) THEN
         field_value = CONJG(value)*refl
      ELSE
         field_value = CONJG(value)*refl/(RHO*G)
      END IF
   END FUNCTION field_value

   PURE COMPLEX(8) FUNCTION output_value(value, refl, moment)
      COMPLEX(8), INTENT(IN) :: value
      REAL(8), INTENT(IN) :: refl
      LOGICAL, INTENT(IN) :: moment

      IF (moment) THEN
         output_value = CONJG(value)/(RHO*G*refl**2)
      ELSE
         output_value = CONJG(value)/(RHO*G*refl)
      END IF
   END FUNCTION output_value

END MODULE SecondOrderOutput
