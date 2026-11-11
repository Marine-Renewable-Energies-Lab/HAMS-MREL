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
MODULE QTFFreeSurface
   ! Free-surface integral of the second-order potential in the difference-frequency QTF (Phase 6a). The forcing of the
   ! second-order potential phi_D (diffracted and radiated part) on z = 0 is
   !   -Omega^2 phi_D + g d phi_D/dz = Q_D,   Q_D = Q[tot_1, tot_2] - Q[inc_1, inc_2],
   ! with Q of SecondOrderIncident (free_surface_forcing): the incident-incident part is carried by the incident
   ! second-order potential, so only the difference of the total and the incident forcing enters here.

   USE SecondOrderIncident, ONLY: free_surface_forcing
   USE NearFieldDrift, ONLY: KERNEL_DIFF
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fs_point_field, fs_forcing_diff

   TYPE :: fs_point_field
      ! One first-order state at one point of z = 0: potential, gradient, and d^2 phi / dz^2.
      COMPLEX(8) :: phi = (0D0, 0D0)
      COMPLEX(8) :: grad(3) = (0D0, 0D0)
      COMPLEX(8) :: phizz = (0D0, 0D0)
   END TYPE fs_point_field

CONTAINS

   PURE COMPLEX(8) FUNCTION fs_forcing_diff(w1, tot1, inc1, w2, tot2, inc2, g)
      ! Q_D = Q[tot1, tot2] - Q[inc1, inc2] at the difference frequency (KERNEL_DIFF), for the states 1 (omega w1) and 2
      ! (omega w2), total and incident, at one point of z = 0.
      REAL(8), INTENT(IN) :: w1, w2, g
      TYPE(fs_point_field), INTENT(IN) :: tot1, inc1, tot2, inc2

      fs_forcing_diff = free_surface_forcing(KERNEL_DIFF, w1, tot1%phi, tot1%grad, tot1%phizz, &
                                             w2, tot2%phi, tot2%grad, tot2%phizz, g) &
                        - free_surface_forcing(KERNEL_DIFF, w1, inc1%phi, inc1%grad, inc1%phizz, &
                                               w2, inc2%phi, inc2%grad, inc2%phizz, g)
   END FUNCTION fs_forcing_diff

END MODULE QTFFreeSurface
