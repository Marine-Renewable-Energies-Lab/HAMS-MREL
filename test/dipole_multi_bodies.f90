!> Test geometry for dipole panels on several bodies: copies of one body of dipole_bodies (column, plate as a sheet or
!> as a closed plate, optional lid), placed at x = XS(b), written straight into the multi-body arrays the program fills
!> from HullMesh_N.pnl (the combined hull, lid and node arrays of PanelMesh_mod and Inerfs_mod) and into the raw dipole
!> block of DipolePanelMulti_mod, as DIPOLE_READ_PANELS_MULTI leaves it. Each body's node list is its hull nodes, then
!> its sheet nodes; the hull, sheet and lid connectivity index the combined node lists. XR_MULTI(b) = (XS(b), 0, 0).
!> The single-body arrays of dipole_bodies keep the body at x = 0, so a test can run the single path on the same body.
module dipole_multi_bodies

  use PanelMesh_mod
  use Inerfs_mod, only : iNELEM, iNTND, iXYZ, iNCN, iNCON, iXYZ_P, iDXYZ_P, iDS, iPNSZ, tNELEM
  use Potentials_mod, only : CGRN_MULTI_COMB, RKBN_MULTI_COMB, DGRN_MULTI_COMB, PKBN_MULTI_COMB, MXPOT_MULTI_COMB
  use DipolePanel_mod, only : DNELEM, DPXYZ, DNCN, DNCON
  use DipolePanelMulti_mod
  use Body_mod, only : XR_MULTI
  use WaveDyn_mod, only : NBODY
  use dipole_bodies, only : build_body, release_body

  implicit none
  private
  public :: build_bodies, release_bodies, hull_start

contains

  !> First hull panel of body b in the combined hull arrays.
  integer function hull_start(b)
    integer, intent(in) :: b
    hull_start = sum(NELEM_MULTI(1:b-1)) + 1
  end function hull_start

  !> kind, t, b, th, nt, nz, nr, lid as build_body; xs(nb) the x positions of the nb copies.
  subroutine build_bodies(kind, t, bp, th, nt, nz, nr, lid, xs)
    integer, intent(in) :: kind, nt, nz, nr
    real(8), intent(in) :: t, bp, th, xs(:)
    logical, intent(in) :: lid
    integer :: nb, b, nh, nd, nl, nn, ndn, nln, k, e0, n0, l0, ln0, d0
    real(8) :: sh(3)

    call release_bodies
    call build_body(kind, t, bp, th, nt, nz, nr, lid)
    nb = size(xs)
    NBODY = nb
    nh = NELEM
    nd = DNELEM
    nl = iNELEM
    nn = NTND
    ndn = 4*nd
    nln = iNTND
    allocate(NELEM_MULTI(nb), NTND_MULTI(nb), XR_MULTI(nb,3))
    NELEM_MULTI = nh
    NTND_MULTI = nn + ndn
    XR_MULTI = 0d0
    XR_MULTI(:,1) = xs
    tNELEM = nb*nh
    allocate(XYZ_GLOBAL_MULTI_COMB(nb*(nn + ndn),3), XYZ_GLOBAL_MULTI_COMB_P(nb*nh,3), PNSZ_MULTI_COMB(nb*nh))
    allocate(NCN_MULTI_COMB(nb*nh), NCON_MULTI_COMB(nb*nh,4), DS_MULTI_COMB(nb*nh), DXYZ_MULTI_COMB(nb*nh,6))
    allocate(GEONRM_MULTI_COMB(nb*nh,3))
    NCON_MULTI_COMB = 0
    call DIPOLE_MULTI_INIT(nb, [(nd, b = 1, nb)])
    do b = 1, nb
      sh = [xs(b), 0d0, 0d0]
      e0 = (b - 1)*nh
      n0 = (b - 1)*(nn + ndn)
      do k = 1, nn
        XYZ_GLOBAL_MULTI_COMB(n0 + k,:) = XYZ(k,:) + sh
      end do
      do k = 1, ndn
        XYZ_GLOBAL_MULTI_COMB(n0 + nn + k,:) = DPXYZ(k,:) + sh
      end do
      do k = 1, nh
        XYZ_GLOBAL_MULTI_COMB_P(e0 + k,:) = XYZ_P(k,:) + sh
        NCN_MULTI_COMB(e0 + k) = NCN(k)
        NCON_MULTI_COMB(e0 + k,1:NCN(k)) = NCON(k,1:NCN(k)) + n0
      end do
      PNSZ_MULTI_COMB(e0+1:e0+nh) = PNSZ
      DS_MULTI_COMB(e0+1:e0+nh) = DS
      DXYZ_MULTI_COMB(e0+1:e0+nh,:) = DXYZ_P          ! (x - XR_MULTI(b)) x n: the same as about 0 for the body at 0
      GEONRM_MULTI_COMB(e0+1:e0+nh,:) = GEONRM_P
      if (nd > 0) then
        d0 = DSTART(b) - 1
        do k = 1, nd
          DNCN_COMB(d0 + k) = DNCN(k)
          DNCON_COMB(d0 + k,1:DNCN(k)) = DNCON(k,1:DNCN(k)) + n0 + nn
        end do
      end if
    end do
    allocate(RKBN_MULTI_COMB(nb*nh,nb*nh,1,4), CGRN_MULTI_COMB(nb*nh,nb*nh,1,4))

    if (lid) then
      allocate(iNELEM_MULTI(nb), iNTND_MULTI(nb))
      iNELEM_MULTI = nl
      iNTND_MULTI = nln
      NELEM_TOTAL = nb*nh
      iNELEM_TOTAL = nb*nl
      tNELEM = nb*(nh + nl)
      allocate(iXYZ_GLOBAL_MULTI_COMB(nb*nln,3), iXYZ_GLOBAL_MULTI_COMB_P(nb*nl,3), iPNSZ_MULTI_COMB(nb*nl))
      allocate(iNCN_MULTI_COMB(nb*nl), iNCON_MULTI_COMB(nb*nl,4), iDS_MULTI_COMB(nb*nl), iDXYZ_MULTI_COMB(nb*nl,6))
      iNCON_MULTI_COMB = 0
      do b = 1, nb
        sh = [xs(b), 0d0, 0d0]
        l0 = (b - 1)*nl
        ln0 = (b - 1)*nln
        do k = 1, nln
          iXYZ_GLOBAL_MULTI_COMB(ln0 + k,:) = iXYZ(k,:) + sh
        end do
        do k = 1, nl
          iXYZ_GLOBAL_MULTI_COMB_P(l0 + k,:) = iXYZ_P(k,:) + sh
          iNCN_MULTI_COMB(l0 + k) = iNCN(k)
          iNCON_MULTI_COMB(l0 + k,1:iNCN(k)) = iNCON(k,1:iNCN(k)) + ln0
        end do
        iPNSZ_MULTI_COMB(l0+1:l0+nl) = iPNSZ
        iDS_MULTI_COMB(l0+1:l0+nl) = iDS
        iDXYZ_MULTI_COMB(l0+1:l0+nl,:) = iDXYZ_P
      end do
      allocate(DGRN_MULTI_COMB(nb*nl,nb*nh,1,4), PKBN_MULTI_COMB(nb*nl,nb*nh,1,4))
    end if
  end subroutine build_bodies

  !> Free what build_bodies allocated, and the single-body arrays of dipole_bodies.
  subroutine release_bodies
    if (allocated(NELEM_MULTI)) deallocate(NELEM_MULTI)
    if (allocated(NTND_MULTI)) deallocate(NTND_MULTI)
    if (allocated(XR_MULTI)) deallocate(XR_MULTI)
    if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate(XYZ_GLOBAL_MULTI_COMB)
    if (allocated(XYZ_GLOBAL_MULTI_COMB_P)) deallocate(XYZ_GLOBAL_MULTI_COMB_P)
    if (allocated(PNSZ_MULTI_COMB)) deallocate(PNSZ_MULTI_COMB)
    if (allocated(NCN_MULTI_COMB)) deallocate(NCN_MULTI_COMB)
    if (allocated(NCON_MULTI_COMB)) deallocate(NCON_MULTI_COMB)
    if (allocated(DS_MULTI_COMB)) deallocate(DS_MULTI_COMB)
    if (allocated(DXYZ_MULTI_COMB)) deallocate(DXYZ_MULTI_COMB)
    if (allocated(GEONRM_MULTI_COMB)) deallocate(GEONRM_MULTI_COMB)
    if (allocated(RKBN_MULTI_COMB)) deallocate(RKBN_MULTI_COMB)
    if (allocated(CGRN_MULTI_COMB)) deallocate(CGRN_MULTI_COMB)
    if (allocated(DGRN_MULTI_COMB)) deallocate(DGRN_MULTI_COMB)
    if (allocated(PKBN_MULTI_COMB)) deallocate(PKBN_MULTI_COMB)
    if (allocated(MXPOT_MULTI_COMB)) deallocate(MXPOT_MULTI_COMB)
    if (allocated(iNELEM_MULTI)) deallocate(iNELEM_MULTI)
    if (allocated(iNTND_MULTI)) deallocate(iNTND_MULTI)
    if (allocated(iXYZ_GLOBAL_MULTI_COMB)) deallocate(iXYZ_GLOBAL_MULTI_COMB)
    if (allocated(iXYZ_GLOBAL_MULTI_COMB_P)) deallocate(iXYZ_GLOBAL_MULTI_COMB_P)
    if (allocated(iPNSZ_MULTI_COMB)) deallocate(iPNSZ_MULTI_COMB)
    if (allocated(iNCN_MULTI_COMB)) deallocate(iNCN_MULTI_COMB)
    if (allocated(iNCON_MULTI_COMB)) deallocate(iNCON_MULTI_COMB)
    if (allocated(iDS_MULTI_COMB)) deallocate(iDS_MULTI_COMB)
    if (allocated(iDXYZ_MULTI_COMB)) deallocate(iDXYZ_MULTI_COMB)
    call DIPOLE_MULTI_RELEASE
    call release_body
    NBODY = 1
  end subroutine release_bodies

end module dipole_multi_bodies
