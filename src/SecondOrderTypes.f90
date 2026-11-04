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

   PUBLIC :: so_options

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

END MODULE SecondOrderTypes
