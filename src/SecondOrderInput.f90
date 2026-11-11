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
MODULE SecondOrderInput
   ! Reader and validator for the optional SecondOrder.in file (mean wave drift options).

   USE SecondOrderTypes
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: read_second_order_options
   PUBLIC :: validate_second_order

CONTAINS

   ! Reads <inputdir>/SecondOrder.in. Lines are "key = value", "!" starts a comment, key order is free.
   ! An absent file is not an error: it leaves opts%active false.
   SUBROUTINE read_second_order_options(inputdir, opts, ok, msg)
      CHARACTER(LEN=*), INTENT(IN)  :: inputdir
      TYPE(so_options), INTENT(OUT) :: opts
      LOGICAL,          INTENT(OUT) :: ok
      CHARACTER(LEN=*), INTENT(OUT) :: msg

      CHARACTER(LEN=512) :: line
      CHARACTER(LEN=:), ALLOCATABLE :: key, value
      LOGICAL :: exists
      INTEGER :: unit, ios, ic, ie, iline

      ok = .TRUE.
      msg = ''
      INQUIRE(FILE=TRIM(inputdir)//'/SecondOrder.in', EXIST=exists)
      IF (.NOT. exists) RETURN

      OPEN(NEWUNIT=unit, FILE=TRIM(inputdir)//'/SecondOrder.in', STATUS='OLD', ACTION='READ', IOSTAT=ios)
      IF (ios /= 0) THEN
         ok = .FALSE.
         msg = 'SecondOrder.in could not be opened.'
         RETURN
      END IF

      iline = 0
      DO
         READ(unit, '(A)', IOSTAT=ios) line
         IF (ios /= 0) EXIT
         iline = iline + 1
         ic = INDEX(line, '!')
         IF (ic > 0) line(ic:) = ' '
         IF (LEN_TRIM(line) == 0) CYCLE
         ie = INDEX(line, '=')
         IF (ie == 0) THEN
            ok = .FALSE.
            msg = 'Line '//itoa(iline)//' of SecondOrder.in has no "=".'
            EXIT
         END IF
         key = lower(TRIM(ADJUSTL(line(:ie-1))))
         value = lower(TRIM(ADJUSTL(line(ie+1:))))
         CALL set_option(key, value, opts, ok, msg)
         IF (.NOT. ok) EXIT
      END DO
      CLOSE(unit)
      IF (.NOT. ok) RETURN

      IF (opts%mean_drift < 0 .OR. opts%mean_drift > 2) THEN
         ok = .FALSE.
         msg = 'mean_drift must be 0, 1 or 2 in SecondOrder.in.'
      ELSE IF (opts%n_azimuth < 1) THEN
         ok = .FALSE.
         msg = 'n_azimuth must be positive in SecondOrder.in.'
      ELSE IF (opts%write_terms < 0 .OR. opts%write_terms > 2) THEN
         ok = .FALSE.
         msg = 'write_terms must be 0, 1 or 2 in SecondOrder.in.'
      ELSE IF (opts%waterline_offset <= 0D0) THEN
         ok = .FALSE.
         msg = 'waterline_offset must be positive in SecondOrder.in.'
      ELSE IF (opts%far_tolerance <= 0D0) THEN
         ok = .FALSE.
         msg = 'far_tolerance must be positive in SecondOrder.in.'
      ELSE IF (opts%wl_z_tolerance <= 0D0) THEN
         ok = .FALSE.
         msg = 'wl_z_tolerance must be positive in SecondOrder.in.'
      ELSE IF (opts%qtf < 0 .OR. opts%qtf > 3) THEN
         ok = .FALSE.
         msg = 'qtf must be 0, 1, 2 or 3 in SecondOrder.in.'
      ELSE IF (opts%qtf_max_dw < 0D0) THEN
         ok = .FALSE.
         msg = 'qtf_max_dw must not be negative in SecondOrder.in.'
      ELSE IF (.NOT. (opts%qtf_edge_angle > 0D0 .AND. opts%qtf_edge_angle < 180D0)) THEN
         ok = .FALSE.
         msg = 'qtf_edge_angle must be between 0 and 180 degrees in SecondOrder.in.'
      ELSE IF (opts%qtf > 0 .AND. opts%mean_drift == 0) THEN
         ok = .FALSE.
         msg = 'qtf needs mean_drift = 1 or 2.'
      ELSE IF (opts%qtf_field < 0 .OR. opts%qtf_field > 1) THEN
         ok = .FALSE.
         msg = 'qtf_field must be 0 or 1 in SecondOrder.in.'
      ELSE IF (opts%qtf_field == 1 .AND. opts%qtf == 0) THEN
         ok = .FALSE.
         msg = 'qtf_field needs qtf = 1, 2 or 3 in SecondOrder.in.'
      ELSE IF (opts%qtf_terms == QTF_TERMS_BODY .AND. opts%qtf == 2) THEN
         ok = .FALSE.
         msg = 'qtf_terms = body needs the difference QTF (qtf = 1 or 3); sum-frequency body forcing needs the free-surface'// &
               ' integrals (Phase 6).'
      END IF
      IF (.NOT. ok) RETURN

      opts%active = opts%mean_drift > 0
      ! Without drift_methods the far-field method is used
      IF (opts%active .AND. .NOT. (opts%do_far .OR. opts%do_near .OR. opts%do_control)) opts%do_far = .TRUE.
   END SUBROUTINE read_second_order_options

   SUBROUTINE set_option(key, value, opts, ok, msg)
      CHARACTER(LEN=*), INTENT(IN)    :: key, value
      TYPE(so_options), INTENT(INOUT) :: opts
      LOGICAL,          INTENT(OUT)   :: ok
      CHARACTER(LEN=*), INTENT(OUT)   :: msg

      INTEGER :: ios

      ok = .TRUE.
      msg = ''
      ios = 0
      SELECT CASE (key)
      CASE ('mean_drift')
         READ(value, *, IOSTAT=ios) opts%mean_drift
      CASE ('drift_methods')
         CALL parse_methods(value, opts, ios)
      CASE ('velocity_method')
         CALL choose(value, 'differencing', 'source', opts%velocity_method, ios)
      CASE ('waterline_method')
         CALL choose(value, 'extrapolate', 'field', opts%waterline_method, ios)
      CASE ('waterline_offset')
         READ(value, *, IOSTAT=ios) opts%waterline_offset
      CASE ('motion_state')
         IF (value == 'free') THEN
            opts%motion_fixed = .FALSE.
         ELSE IF (value == 'fixed') THEN
            opts%motion_fixed = .TRUE.
         ELSE
            ios = 1
         END IF
      CASE ('n_azimuth')
         READ(value, *, IOSTAT=ios) opts%n_azimuth
      CASE ('far_tolerance')
         READ(value, *, IOSTAT=ios) opts%far_tolerance
      CASE ('write_terms')
         READ(value, *, IOSTAT=ios) opts%write_terms
      CASE ('wl_z_tolerance')
         READ(value, *, IOSTAT=ios) opts%wl_z_tolerance
      CASE ('cs_velocity')
         CALL choose(value, 'green', 'source', opts%cs_velocity, ios)
      CASE ('cs_shape')
         IF (value == 'cylinder') THEN
            opts%cs_shape = 1
         ELSE IF (value == 'box') THEN
            opts%cs_shape = 2
         ELSE IF (value == 'file') THEN
            opts%cs_shape = 3
         ELSE
            ios = 1
         END IF
      CASE ('cs_margin')
         READ(value, *, IOSTAT=ios) opts%cs_margin
      CASE ('cs_bottom_margin')
         READ(value, *, IOSTAT=ios) opts%cs_bottom_margin
      CASE ('cs_panel_size')
         READ(value, *, IOSTAT=ios) opts%cs_panel_size
      CASE ('cs_fd_step')
         READ(value, *, IOSTAT=ios) opts%cs_fd_step
      CASE ('cs_fs_rows')
         READ(value, *, IOSTAT=ios) opts%cs_fs_rows
      CASE ('qtf')
         READ(value, *, IOSTAT=ios) opts%qtf
      CASE ('qtf_field')
         READ(value, *, IOSTAT=ios) opts%qtf_field
      CASE ('qtf_max_dw')
         READ(value, *, IOSTAT=ios) opts%qtf_max_dw
      CASE ('qtf_terms')
         IF (value == 'quadratic') THEN
            opts%qtf_terms = QTF_TERMS_QUADRATIC
         ELSE IF (value == 'body') THEN
            opts%qtf_terms = QTF_TERMS_BODY
         ELSE IF (value == 'full') THEN
            ok = .FALSE.
            msg = 'qtf_terms = full is not implemented (Phase 6); use quadratic or body.'
            RETURN
         ELSE
            ios = 1
         END IF
      CASE ('qtf_b2')
         ! Read always, used with qtf_terms = body only
         CALL choose(value, 'body', 'control', opts%qtf_b2, ios)
         opts%qtf_b2_given = .TRUE.
      CASE ('qtf_edge_weights')
         ! Read always, used with qtf_terms = body and qtf_b2 = control only; off by default (SecondOrderTypes)
         IF (value == 'on') THEN
            opts%qtf_edge_weights = .TRUE.
         ELSE IF (value == 'off') THEN
            opts%qtf_edge_weights = .FALSE.
         ELSE
            ios = 1
         END IF
         opts%qtf_edge_given = .TRUE.
      CASE ('qtf_edge_angle')
         READ(value, *, IOSTAT=ios) opts%qtf_edge_angle
         opts%qtf_edge_given = .TRUE.
      CASE DEFAULT
         ok = .FALSE.
         msg = 'Unknown key "'//key//'" in SecondOrder.in.'
         RETURN
      END SELECT
      IF (ios /= 0) THEN
         ok = .FALSE.
         msg = 'Invalid value "'//value//'" for key "'//key//'" in SecondOrder.in.'
      END IF
   END SUBROUTINE set_option

   ! drift_methods is a list of the words "far", "near" and "control", separated by blanks or commas.
   SUBROUTINE parse_methods(value, opts, ios)
      CHARACTER(LEN=*), INTENT(IN)    :: value
      TYPE(so_options), INTENT(INOUT) :: opts
      INTEGER,          INTENT(OUT)   :: ios

      CHARACTER(LEN=LEN(value)) :: rest
      INTEGER :: n, i

      ios = 0
      opts%do_far = .FALSE.
      opts%do_near = .FALSE.
      opts%do_control = .FALSE.
      rest = value
      DO i = 1, LEN(rest)
         IF (rest(i:i) == ',') rest(i:i) = ' '
      END DO
      rest = ADJUSTL(rest)
      DO WHILE (LEN_TRIM(rest) > 0)
         n = INDEX(TRIM(rest), ' ')
         IF (n == 0) n = LEN_TRIM(rest) + 1
         SELECT CASE (rest(:n-1))
         CASE ('far')
            opts%do_far = .TRUE.
         CASE ('near')
            opts%do_near = .TRUE.
         CASE ('control')
            opts%do_control = .TRUE.
         CASE DEFAULT
            ios = 1
            RETURN
         END SELECT
         rest = ADJUSTL(rest(n:))
      END DO
      IF (.NOT. (opts%do_far .OR. opts%do_near .OR. opts%do_control)) ios = 1
   END SUBROUTINE parse_methods

   SUBROUTINE choose(value, word1, word2, code, ios)
      CHARACTER(LEN=*), INTENT(IN)  :: value, word1, word2
      INTEGER,          INTENT(OUT) :: code, ios

      ios = 0
      code = 1
      IF (value == word1) THEN
         code = 1
      ELSE IF (value == word2) THEN
         code = 2
      ELSE
         ios = 1
      END IF
   END SUBROUTINE choose

   ! Checks the options against the first-order set-up. Depths above 500 m are rejected because the incident
   ! wave velocity routines VINP and DINP in PatclVelct.f90 treat them differently, and first-order results
   ! must stay unchanged.
   SUBROUTINE validate_second_order(opts, isol, h, isolv, irsp, has_gen_modes, ok, msg)
      TYPE(so_options), INTENT(IN)  :: opts
      INTEGER,          INTENT(IN)  :: isol, isolv, irsp
      REAL(8),          INTENT(IN)  :: h
      LOGICAL,          INTENT(IN)  :: has_gen_modes
      LOGICAL,          INTENT(OUT) :: ok
      CHARACTER(LEN=*), INTENT(OUT) :: msg

      ok = .FALSE.
      msg = ''
      IF (isol == 2) THEN
         msg = 'Mean drift needs ISOL=1 (incident and scattering potentials combined) in ControlFile.in.'
      ELSE IF (h > 500.0D0) THEN
         msg = 'Mean drift does not support a water depth above 500 m. Set H<0 in ControlFile.in for infinite depth.'
      ELSE IF (opts%do_near .AND. has_gen_modes) THEN
         msg = 'The near-field drift method does not support generalized modes (gen_mod files). Use far only.'
      ELSE IF (opts%velocity_method == 2 .AND. isolv == 2) THEN
         msg = 'velocity_method = source needs the direct solver (ISOLV=1) in ControlFile.in.'
      ELSE IF (opts%velocity_method == 2 .AND. irsp == 1) THEN
         msg = 'velocity_method = source cannot be used with irregular frequency removal (IRSP=1) in ControlFile.in.'
      ELSE IF (opts%do_control .AND. has_gen_modes) THEN
         msg = 'drift_methods = control cannot be used with generalized modes.'
      ELSE IF (opts%cs_velocity == 2 .AND. isolv == 2) THEN
         msg = 'cs_velocity = source needs the direct solver (ISOLV=1) in ControlFile.in.'
      ELSE IF (opts%cs_velocity == 2 .AND. irsp == 1) THEN
         msg = 'cs_velocity = source cannot be used with irregular frequency removal (IRSP=1) in ControlFile.in.'
      ELSE IF (opts%cs_margin < 0D0) THEN
         msg = 'cs_margin must not be negative in SecondOrder.in.'
      ELSE IF (opts%cs_bottom_margin < 0D0) THEN
         msg = 'cs_bottom_margin must not be negative in SecondOrder.in.'
      ELSE IF (opts%cs_panel_size < 0D0) THEN
         msg = 'cs_panel_size must not be negative in SecondOrder.in.'
      ELSE IF (opts%cs_fd_step <= 0D0 .OR. opts%cs_fd_step > 0.5D0) THEN
         msg = 'cs_fd_step must be above 0 and at most 0.5 in SecondOrder.in.'
      ELSE IF (opts%cs_fs_rows < 0) THEN
         msg = 'cs_fs_rows must not be negative in SecondOrder.in.'
      ELSE IF (opts%qtf > 0 .AND. .NOT. opts%do_near) THEN
         msg = 'qtf needs drift_methods to include near.'
      ELSE
         ok = .TRUE.
      END IF
   END SUBROUTINE validate_second_order

   FUNCTION lower(s) RESULT(r)
      CHARACTER(LEN=*), INTENT(IN) :: s
      CHARACTER(LEN=LEN(s)) :: r
      INTEGER :: i, c

      r = s
      DO i = 1, LEN(s)
         c = IACHAR(s(i:i))
         IF (c >= IACHAR('A') .AND. c <= IACHAR('Z')) r(i:i) = ACHAR(c + 32)
      END DO
   END FUNCTION lower

   FUNCTION itoa(i) RESULT(r)
      INTEGER, INTENT(IN) :: i
      CHARACTER(LEN=:), ALLOCATABLE :: r
      CHARACTER(LEN=16) :: buf

      WRITE(buf, '(I0)') i
      r = TRIM(buf)
   END FUNCTION itoa

END MODULE SecondOrderInput
