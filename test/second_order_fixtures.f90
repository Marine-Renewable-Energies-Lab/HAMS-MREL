!> Geometry fixtures shared by the second-order unit tests.
module second_order_fixtures
  use SecondOrderTypes, only : panel_mesh
  implicit none
  private

  public :: make_cylinder, set_panel

contains

!> Fills panel i from its nodes (already stored in mesh%xnode and listed in mesh%ncon) and a given unit normal.
!> The centroid is the vertex average, the area comes from the diagonals, and nj uses mesh%xr(:,1) as moment origin.
subroutine set_panel(mesh, i, nrm)
  type(panel_mesh), intent(inout) :: mesh
  integer, intent(in) :: i
  real(8), intent(in) :: nrm(3)
  real(8) :: v(3,4), d1(3), d2(3), c(3), x(3)
  integer :: m

  do m = 1, mesh%ncn(i)
    v(:,m) = mesh%xnode(:, mesh%ncon(m,i))
  end do
  mesh%xc(:,i) = sum(v(:,1:mesh%ncn(i)), dim=2) / real(mesh%ncn(i), 8)
  if (mesh%ncn(i) == 4) then
    d1 = v(:,3) - v(:,1)
    d2 = v(:,4) - v(:,2)
  else
    d1 = v(:,2) - v(:,1)
    d2 = v(:,3) - v(:,1)
  end if
  c = [d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1)]
  mesh%area(i) = 0.5d0 * sqrt(sum(c**2))
  mesh%psize(i) = sqrt(mesh%area(i))
  mesh%nrm(:,i) = nrm
  x = mesh%xc(:,i) - mesh%xr(:,1)
  mesh%nj(1:3,i) = nrm
  mesh%nj(4,i) = x(2)*nrm(3) - x(3)*nrm(2)
  mesh%nj(5,i) = x(3)*nrm(1) - x(1)*nrm(3)
  mesh%nj(6,i) = x(1)*nrm(2) - x(2)*nrm(1)
end subroutine set_panel

!> Vertical cylinder about the z axis with its free surface at z = 0 (open top), single body, no symmetry.
!> Normals point into the body. Side wall: ntheta x nz panels. Bottom (optional): nr rings of ntheta panels,
!> the innermost ring being triangles.
subroutine make_cylinder(radius, draft, ntheta, nz, nr, mesh, with_bottom)
  real(8), intent(in) :: radius, draft
  integer, intent(in) :: ntheta, nz, nr
  type(panel_mesh), intent(out) :: mesh
  logical, intent(in) :: with_bottom
  real(8), parameter :: pi = 3.141592653589793d0
  integer :: nside, nnode, npan, it, iz, ir, ip, itn, i
  real(8) :: th, thm, r

  nside = ntheta * (nz + 1)
  nnode = nside
  npan = ntheta * nz
  if (with_bottom) then
    nnode = nside + ntheta * (nr - 1) + 1
    npan = npan + ntheta * nr
  end if
  mesh%n = npan
  mesh%nbase = npan
  mesh%nnode = nnode
  mesh%nbody = 1
  mesh%nsys = 1
  allocate(mesh%xc(3,npan), mesh%nrm(3,npan), mesh%area(npan), mesh%psize(npan), mesh%nj(6,npan))
  allocate(mesh%xnode(3,nnode), mesh%ncn(npan), mesh%ncon(4,npan), mesh%body(npan), mesh%base(npan), mesh%image(npan))
  allocate(mesh%xr(3,1))
  mesh%xr = 0d0
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(i, i = 1, npan)]
  mesh%ncon = 0

  ! Side nodes: index it + ntheta*iz, iz = 0..nz
  do iz = 0, nz
    do it = 1, ntheta
      th = 2d0 * pi * real(it - 1, 8) / real(ntheta, 8)
      mesh%xnode(:, it + ntheta*iz) = [radius*cos(th), radius*sin(th), -draft*real(iz, 8)/real(nz, 8)]
    end do
  end do
  if (with_bottom) then
    do ir = 1, nr - 1
      r = radius * real(ir, 8) / real(nr, 8)
      do it = 1, ntheta
        th = 2d0 * pi * real(it - 1, 8) / real(ntheta, 8)
        mesh%xnode(:, nside + (ir-1)*ntheta + it) = [r*cos(th), r*sin(th), -draft]
      end do
    end do
    mesh%xnode(:, nnode) = [0d0, 0d0, -draft]
  end if

  ip = 0
  do iz = 1, nz
    do it = 1, ntheta
      itn = mod(it, ntheta) + 1
      ip = ip + 1
      mesh%ncn(ip) = 4
      mesh%ncon(1:4, ip) = [it + ntheta*(iz-1), it + ntheta*iz, itn + ntheta*iz, itn + ntheta*(iz-1)]
      thm = 2d0 * pi * (real(it - 1, 8) + 0.5d0) / real(ntheta, 8)
      call set_panel(mesh, ip, [-cos(thm), -sin(thm), 0d0])
    end do
  end do
  if (with_bottom) then
    do ir = 1, nr
      do it = 1, ntheta
        itn = mod(it, ntheta) + 1
        ip = ip + 1
        if (ir == 1) then
          mesh%ncn(ip) = 3
          mesh%ncon(1:3, ip) = [nnode, bottom_node(1, it), bottom_node(1, itn)]
        else
          mesh%ncn(ip) = 4
          mesh%ncon(1:4, ip) = [bottom_node(ir-1, it), bottom_node(ir, it), bottom_node(ir, itn), bottom_node(ir-1, itn)]
        end if
        call set_panel(mesh, ip, [0d0, 0d0, 1d0])
      end do
    end do
  end if

contains

  integer function bottom_node(ring, j)
    integer, intent(in) :: ring, j
    if (ring == nr) then
      bottom_node = j + ntheta*nz
    else
      bottom_node = nside + (ring-1)*ntheta + j
    end if
  end function bottom_node

end subroutine make_cylinder

end module second_order_fixtures
