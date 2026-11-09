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

   PUBLIC :: write_drift_row, write_terms_row, write_qtf_row, write_qtf_terms_row

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

   SUBROUTINE write_qtf_terms_row(unit, oufr1, oufr2, beta1_deg, beta2_deg, index, values, refl, moment)
      ! Writes one row PER1 PER2 BETA1 BETA2 I, then Re Im of each of values, of QTF_PotentialTerms.txt. The keys are as in
      ! write_qtf_row; each value is scaled and conjugated as there, and written with 16 digits, so that the sum of the
      ! terms can be checked against QTF_Indirect.11d - QTF_Quadratic.10d.
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
