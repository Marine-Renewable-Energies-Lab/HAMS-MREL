!> Unit tests of PartialWaves (Phase 6a, Task 5): the partial-wave decomposition of a first-order field on a vertical
!> cylinder r = R, and its reconstruction on z = 0 outside R. pw_hams_field sets module globals (the mesh arrays and the wave
!> of WaveDyn_mod), so the suite runs serially (see unit_tester.f90).
module test_partial_waves
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD
  use WaveDyn_mod, only : BETA, H, NK, W1, WK, WVN
  use SecondOrderTypes, only : panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL, CS_LINE, CS_FREE
  use FirstOrderState, only : wave_gravity
  use FieldVelocity, only : fields_green
  use PartialWaves, only : pw_set, pw_radial, pw_sample_points, pw_decompose, pw_eval_surface, pw_eval_point, pw_radial_at, &
                           pw_bessel_k, pw_evanescent_from_wvn, pw_bessel_k_scaled, pw_bessel_k_scaled_c, pw_hankel_scaled
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state
  use test_waterline_field, only : set_wave, load_panels
  implicit none
  private

  public :: collect_partial_waves

  real(8), parameter :: pi = 3.14159265358979323846d0
  complex(8), parameter :: ci = (0d0, 1d0)
  ! The g of DISPERSION (WavDynSubs.f90), which fills WVN
  real(8), parameter :: g_dispersion = 9.807d0

  ! K_m(x), m = 0..10 (rows), x = 0.1, 1, 10 (columns). Generated with mpmath 1.2.1 at 30 digits:
  !   python -I -c "import mpmath as mp; mp.mp.dps=30; [print([mp.nstr(mp.besselk(m, mp.mpf(x)), 17) for m in range(11)])
  !                 for x in ['0.1', '1', '10']]"
  real(8), parameter :: k_table(0:10, 3) = reshape([ &
    2.4270690247020166d0, 9.8538447808706061d0, 1.9950396464211414d+2, 7.9900124304654362d+3, 4.7960024979256828d+5, &
    3.8376009995835928d+7, 3.8380805998333854d+9, 4.6060804799000208d+11, 6.4488964799200125d+13, 1.031869497592001d+16, &
    1.857429584630401d+18, &
    4.2102443824070833d-1, 6.0190723019723457d-1, 1.6248388986351775d0, 7.1012628247379445d0, 4.4232415847062845d+1, &
    3.609605896012407d+2, 3.6538383118594699d+3, 4.4207020331914879d+4, 6.2255212295866777d+5, 1.0005040987670599d+7, &
    1.8071328990102945d+8, &
    1.7780062316167652d-5, 1.8648773453825585d-5, 2.1509817006932769d-5, 2.7252700256598692d-5, 3.7861437160891984d-5, &
    5.7541849985312279d-5, 9.5403287146204263d-5, 1.720257945607574d-4, 3.3623939953126462d-4, 7.1000883381078078d-4, &
    1.61425530039067d-3], [11, 3])
  real(8), parameter :: k_args(3) = [0.1d0, 1d0, 10d0]

  ! e^x K_m(x) far out (the range of the outer radial integral of QTFFreeSurface, Task 7; review M1 of Task 5): m = 0, 1, 5,
  ! 10, 40 (rows) at x = 600, 2000, 1e4 (columns). mpmath 1.2.1 at 60 digits (scipy kve agrees to 4e-16):
  !   python -I -c "import mpmath as mp; mp.mp.dps=60; [print([mp.nstr(mp.besselk(m, x)*mp.exp(x), 18) for m in
  !                 (0, 1, 5, 10, 40)]) for x in (600, 2000, 10000)]"
  real(8), parameter :: kfar_table(5, 3) = reshape([ &
    5.11556857202359605d-02, 5.11982977254724428d-02, 5.22316967228415585d-02, 5.55973611433654077d-02, 1.93757455427221670d-01, &
    2.80232050146043202d-02, 2.80302099405703373d-02, 2.81988543876665367d-02, 2.87324348517065224d-02, 4.18009747997934861d-02, &
    1.25329847176992859d-02, 1.25336113512705058d-02, 1.25486599595387317d-02, 1.25958034111704649d-02, 1.35767645150143869d-02], &
    [5, 3])
  real(8), parameter :: kfar_args(3) = [600d0, 2000d0, 1d4]
  integer, parameter :: kfar_orders(5) = [0, 1, 5, 10, 40]
  ! e^w K_m(w) at complex w (the Hankel functions of the tail contours, H_m(z) = (2/pi) i^-(m+1) K_m(-i z)): m = 0, 1, 12
  ! (rows) at w = 2i, 3i, 0.5 + 30i, 2 + 0.1i, 1 + 1.5i (columns), re and im. mpmath at 60 digits (scipy kve agrees to 2e-15).
  real(8), parameter :: kc_table(2, 3, 5) = reshape([ &
    6.53411265314505774d-01, -5.82626965863348500d-01, 5.29871256808186808d-01, -7.53783022009367643d-01, &
    -9.09987791671622358d+06, 1.98835959988381416d+07, &
    5.28384701332361040d-01, -4.87937343470068785d-01, 4.55291808365371309d-01, -5.80053060508276652d-01, &
    -1.87254778536004247d+05, 2.66925213548079664d+04, &
    1.63776863865670336d-01, -1.59734902200358242d-01, 1.61184958988351695d-01, -1.62528699625233142d-01, &
    -2.48622112564824332d-01, 1.51493704878725675d-02, &
    8.40904669107718217d-01, -1.91650846726027710d-02, 1.03204084563524634d+00, -3.24148464703462036d-02, &
    1.15944065635207564d+08, -6.46278295345855802d+07, &
    8.08489914254441611d-01, -3.88176581561545719d-01, 8.57072866860195681d-01, -6.13068093047916096d-01, &
    -1.14851330700899497d+08, 1.56311577571564078d+08], [2, 3, 5])
  complex(8), parameter :: kc_args(5) = [(0d0, 2d0), (0d0, 3d0), (0.5d0, 30d0), (2d0, 0.1d0), (1d0, 1.5d0)]
  integer, parameter :: kc_orders(3) = [0, 1, 12]
  ! The scaled Hankel functions of pw_hankel_scaled: kind 1 H^(1)_m(z) e^{-iz}, kind 2 H^(2)_m(z) e^{iz}, at z in both
  ! half-planes (the tail contours of QTFFreeSurface go up for a phase e^{i kappa r} with kappa > 0, down otherwise, and
  ! each Hankel function is needed on both sides). Columns: kind, m, Re z, Im z, Re, Im. scipy 1.7.3 hankel1e/hankel2e
  ! (Amos), cross-checked against mpmath at 120 digits to 3e-14 (scratch script ref_hankel.py of Task 7).
  real(8), parameter :: hankel_table(6, 144) = reshape([ &
    1d0, 0d0, 3.0000d0, 40.0000d0, 4.67131316757484188d-03, -1.25506764899253498d-01, &
    1d0, 1d0, 3.0000d0, 40.0000d0, -1.27053108545090054d-01, -4.84462919537343838d-03, &
    1d0, 7d0, 3.0000d0, 40.0000d0, 2.28165398525886076d-01, 1.86324664639270950d-02, &
    1d0, 20d0, 3.0000d0, 40.0000d0, 5.83187708700602059d+00, -1.45026782936039886d+01, &
    2d0, 0d0, 3.0000d0, 40.0000d0, 1.26287265488215700d-01, 4.75936871238591668d-03, &
    2d0, 1d0, 3.0000d0, 40.0000d0, -4.58058145460838850d-03, 1.24712049553174137d-01, &
    2d0, 7d0, 3.0000d0, 40.0000d0, -6.02432476824806711d-04, -6.83084739644215233d-02, &
    2d0, 20d0, 3.0000d0, 40.0000d0, 8.62612247494891457d-04, -2.88121938736356188d-04, &
    1d0, 0d0, 3.0000d0, -40.0000d0, 1.26287265488215700d-01, -4.75936871238591668d-03, &
    1d0, 1d0, 3.0000d0, -40.0000d0, -4.58058145460838850d-03, -1.24712049553174137d-01, &
    1d0, 7d0, 3.0000d0, -40.0000d0, -6.02432476824806711d-04, 6.83084739644215233d-02, &
    1d0, 20d0, 3.0000d0, -40.0000d0, 8.62612247494891457d-04, 2.88121938736356188d-04, &
    2d0, 0d0, 3.0000d0, -40.0000d0, 4.67131316757484188d-03, 1.25506764899253498d-01, &
    2d0, 1d0, 3.0000d0, -40.0000d0, -1.27053108545090054d-01, 4.84462919537343838d-03, &
    2d0, 7d0, 3.0000d0, -40.0000d0, 2.28165398525886076d-01, -1.86324664639270950d-02, &
    2d0, 20d0, 3.0000d0, -40.0000d0, 5.83187708700602059d+00, 1.45026782936039886d+01, &
    1d0, 0d0, 14.0000d0, 2.0000d0, 1.37411581631751573d-01, -1.61231389289476779d-01, &
    1d0, 1d0, 14.0000d0, 2.0000d0, -1.57350515781715855d-01, -1.43786973000265716d-01, &
    1d0, 7d0, 14.0000d0, 2.0000d0, -2.21338496326563078d-01, 1.91668252525199256d-01, &
    1d0, 20d0, 14.0000d0, 2.0000d0, 1.05970041722183321d+01, 4.31320065863250406d+01, &
    2d0, 0d0, 14.0000d0, 2.0000d0, 1.59202486507210250d-01, 1.40553719561578527d-01, &
    2d0, 1d0, 14.0000d0, 2.0000d0, -1.34343997919724567d-01, 1.63442226196832130d-01, &
    2d0, 7d0, 14.0000d0, 2.0000d0, -1.51192608129932099d-01, -8.83114080123570400d-02, &
    2d0, 20d0, 14.0000d0, 2.0000d0, 3.99936344037381208d-01, 7.07509063843815289d-01, &
    1d0, 0d0, 14.0000d0, -2.0000d0, 1.59202486507210250d-01, -1.40553719561578527d-01, &
    1d0, 1d0, 14.0000d0, -2.0000d0, -1.34343997919724567d-01, -1.63442226196832130d-01, &
    1d0, 7d0, 14.0000d0, -2.0000d0, -1.51192608129932099d-01, 8.83114080123570400d-02, &
    1d0, 20d0, 14.0000d0, -2.0000d0, 3.99936344037381208d-01, -7.07509063843815289d-01, &
    2d0, 0d0, 14.0000d0, -2.0000d0, 1.37411581631751573d-01, 1.61231389289476779d-01, &
    2d0, 1d0, 14.0000d0, -2.0000d0, -1.57350515781715855d-01, 1.43786973000265716d-01, &
    2d0, 7d0, 14.0000d0, -2.0000d0, -2.21338496326563078d-01, -1.91668252525199256d-01, &
    2d0, 20d0, 14.0000d0, -2.0000d0, 1.05970041722183321d+01, -4.31320065863250406d+01, &
    1d0, 0d0, 125.0000d0, 400.0000d0, 5.87494104461924477d-03, -3.85191289792354424d-02, &
    1d0, 1d0, 125.0000d0, 400.0000d0, -3.85608833288417011d-02, -5.89532025369762618d-03, &
    1d0, 7d0, 125.0000d0, 400.0000d0, 4.06130355606319493d-02, 6.91928511862245459d-03, &
    1d0, 20d0, 125.0000d0, 400.0000d0, 1.77562337358795723d-02, -5.87916633101210548d-02, &
    2d0, 0d0, 125.0000d0, 400.0000d0, 3.85400212897233641d-02, 5.88514496765310188d-03, &
    2d0, 1d0, 125.0000d0, 400.0000d0, -5.86470859287702411d-03, 3.84982064580437211d-02, &
    2d0, 7d0, 125.0000d0, 400.0000d0, 4.92765768385133775d-03, -3.65379916364645188d-02, &
    2d0, 20d0, 125.0000d0, 400.0000d0, 2.47111209315936352d-02, 2.20662851315717976d-04, &
    1d0, 0d0, 125.0000d0, -400.0000d0, 3.85400212897233641d-02, -5.88514496765310188d-03, &
    1d0, 1d0, 125.0000d0, -400.0000d0, -5.86470859287702411d-03, -3.84982064580437211d-02, &
    1d0, 7d0, 125.0000d0, -400.0000d0, 4.92765768385133775d-03, 3.65379916364645188d-02, &
    1d0, 20d0, 125.0000d0, -400.0000d0, 2.47111209315936352d-02, -2.20662851315717976d-04, &
    2d0, 0d0, 125.0000d0, -400.0000d0, 5.87494104461924477d-03, 3.85191289792354424d-02, &
    2d0, 1d0, 125.0000d0, -400.0000d0, -3.85608833288417011d-02, 5.89532025369762618d-03, &
    2d0, 7d0, 125.0000d0, -400.0000d0, 4.06130355606319493d-02, -6.91928511862245459d-03, &
    2d0, 20d0, 125.0000d0, -400.0000d0, 1.77562337358795723d-02, 5.87916633101210548d-02, &
    1d0, 0d0, 0.6000d0, 20.0000d0, 2.62675985846297019d-03, -1.77268767816415623d-01, &
    1d0, 1d0, 0.6000d0, 20.0000d0, -1.81641892604171323d-01, -2.82132890638303992d-03, &
    1d0, 7d0, 0.6000d0, 20.0000d0, 5.78810118584229816d-01, 2.83537698640646402d-02, &
    1d0, 20d0, 0.6000d0, 20.0000d0, 4.30715483188164683d+02, -1.64579015298423974d+03, &
    2d0, 0d0, 0.6000d0, 20.0000d0, 1.79498433156658754d-01, 2.72733127461467816d-03, &
    2d0, 1d0, 0.6000d0, 20.0000d0, -2.52016903374127285d-03, 1.74958131759319374d-01, &
    2d0, 7d0, 0.6000d0, 20.0000d0, -1.15202417236245241d-03, -5.18210133638597747d-02, &
    2d0, 20d0, 0.6000d0, 20.0000d0, 1.28466949798701280d-05, -3.15701550762322698d-06, &
    1d0, 0d0, 0.6000d0, -20.0000d0, 1.79498433156658754d-01, -2.72733127461467816d-03, &
    1d0, 1d0, 0.6000d0, -20.0000d0, -2.52016903374127285d-03, -1.74958131759319374d-01, &
    1d0, 7d0, 0.6000d0, -20.0000d0, -1.15202417236245241d-03, 5.18210133638597747d-02, &
    1d0, 20d0, 0.6000d0, -20.0000d0, 1.28466949798701280d-05, 3.15701550762322698d-06, &
    2d0, 0d0, 0.6000d0, -20.0000d0, 2.62675985846297019d-03, 1.77268767816415623d-01, &
    2d0, 1d0, 0.6000d0, -20.0000d0, -1.81641892604171323d-01, 2.82132890638303992d-03, &
    2d0, 7d0, 0.6000d0, -20.0000d0, 5.78810118584229816d-01, -2.83537698640646402d-02, &
    2d0, 20d0, 0.6000d0, -20.0000d0, 4.30715483188164683d+02, 1.64579015298423974d+03, &
    1d0, 0d0, 2.1000d0, 0.0000d0, 3.63285123877580496d-01, -4.05475280393981619d-01, &
    1d0, 1d0, 2.1000d0, 0.0000d0, -3.31509533058069905d-01, -4.64465348553197710d-01, &
    1d0, 7d0, 2.1000d0, 0.0000d0, -1.69568111961925780d+02, 9.91712568381099970d+01, &
    1d0, 20d0, 2.1000d0, 0.0000d0, -1.33510746885919440d+16, 7.80834674549948100d+15, &
    2d0, 0d0, 2.1000d0, 0.0000d0, 3.63285123877580496d-01, 4.05475280393981619d-01, &
    2d0, 1d0, 2.1000d0, 0.0000d0, -3.31509533058069905d-01, 4.64465348553197710d-01, &
    2d0, 7d0, 2.1000d0, 0.0000d0, -1.69568111961925780d+02, -9.91712568381099970d+01, &
    2d0, 20d0, 2.1000d0, 0.0000d0, -1.33510746885919440d+16, -7.80834674549948100d+15, &
    1d0, 0d0, 2.1000d0, 0.0000d0, 3.63285123877580496d-01, -4.05475280393981619d-01, &
    1d0, 1d0, 2.1000d0, 0.0000d0, -3.31509533058069905d-01, -4.64465348553197710d-01, &
    1d0, 7d0, 2.1000d0, 0.0000d0, -1.69568111961925780d+02, 9.91712568381099970d+01, &
    1d0, 20d0, 2.1000d0, 0.0000d0, -1.33510746885919440d+16, 7.80834674549948100d+15, &
    2d0, 0d0, 2.1000d0, 0.0000d0, 3.63285123877580496d-01, 4.05475280393981619d-01, &
    2d0, 1d0, 2.1000d0, 0.0000d0, -3.31509533058069905d-01, 4.64465348553197710d-01, &
    2d0, 7d0, 2.1000d0, 0.0000d0, -1.69568111961925780d+02, -9.91712568381099970d+01, &
    2d0, 20d0, 2.1000d0, 0.0000d0, -1.33510746885919440d+16, -7.80834674549948100d+15, &
    1d0, 0d0, 1.0000d0, 3.0000d0, 6.52862933829046477d-02, -4.28995556176446913d-01, &
    1d0, 1d0, 1.0000d0, 3.0000d0, -4.86902665718281336d-01, -9.36811015738065189d-02, &
    1d0, 7d0, 1.0000d0, 3.0000d0, 1.11621651413778977d+01, 1.33762409046700469d+02, &
    1d0, 20d0, 1.0000d0, 3.0000d0, -5.10886141722504062d+13, -5.27200219361865234d+13, &
    2d0, 0d0, 1.0000d0, 3.0000d0, 4.60748431932595293d-01, 8.40726319421193496d-02, &
    2d0, 1d0, 1.0000d0, 3.0000d0, -4.10944372581844714d-02, 3.92058009751650893d-01, &
    2d0, 7d0, 1.0000d0, 3.0000d0, 3.12384266126978782d-01, 1.12735725940798337d-01, &
    2d0, 20d0, 1.0000d0, 3.0000d0, -1.71526045181436707d+11, 6.07677871557743149d+10, &
    1d0, 0d0, 1.0000d0, -3.0000d0, 4.60748431932595293d-01, -8.40726319421193496d-02, &
    1d0, 1d0, 1.0000d0, -3.0000d0, -4.10944372581844714d-02, -3.92058009751650893d-01, &
    1d0, 7d0, 1.0000d0, -3.0000d0, 3.12384266126978782d-01, -1.12735725940798337d-01, &
    1d0, 20d0, 1.0000d0, -3.0000d0, -1.71526045181436707d+11, -6.07677871557743149d+10, &
    2d0, 0d0, 1.0000d0, -3.0000d0, 6.52862933829046477d-02, 4.28995556176446913d-01, &
    2d0, 1d0, 1.0000d0, -3.0000d0, -4.86902665718281336d-01, 9.36811015738065189d-02, &
    2d0, 7d0, 1.0000d0, -3.0000d0, 1.11621651413778977d+01, -1.33762409046700469d+02, &
    2d0, 20d0, 1.0000d0, -3.0000d0, -5.10886141722504062d+13, 5.27200219361865234d+13, &
    1d0, 0d0, 0.5000d0, 0.5000d0, 2.74080367005704484d-01, -8.11098157030117783d-01, &
    1d0, 1d0, 0.5000d0, 0.5000d0, -1.07772806491427087d+00, -7.06357352085897117d-01, &
    1d0, 7d0, 0.5000d0, 0.5000d0, 1.64955275839513313d+05, -5.21718178624222521d+05, &
    1d0, 20d0, 0.5000d0, 0.5000d0, 3.24669288858237112d+25, 6.03709546174175549d+25, &
    2d0, 0d0, 0.5000d0, 0.5000d0, 8.27507156287965562d-01, 5.22663167845433430d-01, &
    2d0, 1d0, 0.5000d0, 0.5000d0, 1.41861192379814677d-01, 8.77473898711614853d-01, &
    2d0, 7d0, 0.5000d0, 0.5000d0, -1.94290533060463611d+05, 5.26363583848979979d+04, &
    2d0, 20d0, 0.5000d0, 0.5000d0, 1.22351000348986417d+25, -2.20501582959378641d+25, &
    1d0, 0d0, 0.5000d0, -0.5000d0, 8.27507156287965562d-01, -5.22663167845433430d-01, &
    1d0, 1d0, 0.5000d0, -0.5000d0, 1.41861192379814677d-01, -8.77473898711614853d-01, &
    1d0, 7d0, 0.5000d0, -0.5000d0, -1.94290533060463611d+05, -5.26363583848979979d+04, &
    1d0, 20d0, 0.5000d0, -0.5000d0, 1.22351000348986417d+25, 2.20501582959378641d+25, &
    2d0, 0d0, 0.5000d0, -0.5000d0, 2.74080367005704484d-01, 8.11098157030117783d-01, &
    2d0, 1d0, 0.5000d0, -0.5000d0, -1.07772806491427087d+00, 7.06357352085897117d-01, &
    2d0, 7d0, 0.5000d0, -0.5000d0, 1.64955275839513313d+05, 5.21718178624222521d+05, &
    2d0, 20d0, 0.5000d0, -0.5000d0, 3.24669288858237112d+25, -6.03709546174175549d+25, &
    1d0, 0d0, 300.0000d0, 0.5000d0, 3.25326909217107663d-02, -3.26141240492938828d-02, &
    1d0, 1d0, 300.0000d0, 0.5000d0, -3.25600392395993601d-02, -3.25871828928806645d-02, &
    1d0, 7d0, 300.0000d0, 0.5000d0, 2.98595812204826734d-02, 3.50944203564648449d-02, &
    1d0, 20d0, 300.0000d0, 0.5000d0, 4.58378419136295998d-02, -5.51518528848620772d-03, &
    2d0, 0d0, 300.0000d0, 0.5000d0, 3.25870476435432196d-02, 3.25599030477690954d-02, &
    2d0, 1d0, 300.0000d0, 0.5000d0, -3.25055462303151868d-02, 3.26412687379388639d-02, &
    2d0, 7d0, 300.0000d0, 0.5000d0, 2.97929464343476835d-02, -3.51346464050526358d-02, &
    2d0, 20d0, 300.0000d0, 0.5000d0, 4.57451939922714143d-02, 5.42636279236420691d-03, &
    1d0, 0d0, 300.0000d0, -0.5000d0, 3.25870476435432196d-02, -3.25599030477690954d-02, &
    1d0, 1d0, 300.0000d0, -0.5000d0, -3.25055462303151868d-02, -3.26412687379388639d-02, &
    1d0, 7d0, 300.0000d0, -0.5000d0, 2.97929464343476835d-02, 3.51346464050526358d-02, &
    1d0, 20d0, 300.0000d0, -0.5000d0, 4.57451939922714143d-02, -5.42636279236420691d-03, &
    2d0, 0d0, 300.0000d0, -0.5000d0, 3.25326909217107663d-02, 3.26141240492938828d-02, &
    2d0, 1d0, 300.0000d0, -0.5000d0, -3.25600392395993601d-02, 3.25871828928806645d-02, &
    2d0, 7d0, 300.0000d0, -0.5000d0, 2.98595812204826734d-02, -3.50944203564648449d-02, &
    2d0, 20d0, 300.0000d0, -0.5000d0, 4.58378419136295998d-02, 5.51518528848620772d-03, &
    1d0, 0d0, 40.0000d0, 3000.0000d0, 9.70924944708208949d-05, -1.45657347242112122d-02, &
    1d0, 1d0, 40.0000d0, 3000.0000d0, -1.45681614973511356d-02, -9.71410294855317569d-05, &
    1d0, 7d0, 40.0000d0, 3000.0000d0, 1.46851231943003601d-02, 9.94866159308629153d-05, &
    1d0, 20d0, 40.0000d0, 3000.0000d0, 1.17615866801277461d-04, -1.55694239267619855d-02, &
    2d0, 0d0, 40.0000d0, 3000.0000d0, 1.45669482624376269d-02, 9.71167670353872422d-05, &
    2d0, 1d0, 40.0000d0, 3000.0000d0, -9.70682117977319538d-05, 1.45645208827981758d-02, &
    2d0, 7d0, 40.0000d0, 3000.0000d0, 9.47533644215628240d-05, -1.44484810683359974d-02, &
    2d0, 20d0, 40.0000d0, 3000.0000d0, 1.36275706915300384d-02, 7.87383894280300968d-05, &
    1d0, 0d0, 40.0000d0, -3000.0000d0, 1.45669482624376269d-02, -9.71167670353872422d-05, &
    1d0, 1d0, 40.0000d0, -3000.0000d0, -9.70682117977319538d-05, -1.45645208827981758d-02, &
    1d0, 7d0, 40.0000d0, -3000.0000d0, 9.47533644215628240d-05, 1.44484810683359974d-02, &
    1d0, 20d0, 40.0000d0, -3000.0000d0, 1.36275706915300384d-02, -7.87383894280300968d-05, &
    2d0, 0d0, 40.0000d0, -3000.0000d0, 9.70924944708208949d-05, 1.45657347242112122d-02, &
    2d0, 1d0, 40.0000d0, -3000.0000d0, -1.45681614973511356d-02, 9.71410294855317569d-05, &
    2d0, 7d0, 40.0000d0, -3000.0000d0, 1.46851231943003601d-02, -9.94866159308629153d-05, &
    2d0, 20d0, 40.0000d0, -3000.0000d0, 1.17615866801277461d-04, 1.55694239267619855d-02], [6, 144])

  ! pw_hams_field: the 448-panel sharp cylinder (radius 10, draft 10, sharp bottom edge) of the QTF fixtures, in depth 15
  ! (finite depth, so that the field has evanescent modes; 5 m of water under the keel). Modes chosen so the gate passes.
  real(8), parameter :: hf_depth = 15d0, hf_radius = 15d0
  integer, parameter :: hf_mmax = 16, hf_nev = 8
  ! The gate of the brief, and a limit for the horizontal gradient fixed before the measurement
  real(8), parameter :: hf_phi_limit = 1d-4, hf_grad_limit = 1d-3
  ! With the near-panel integrals of fields_green (see test_hams_field): measured 4.9e-4, the Task 6 warning level
  real(8), parameter :: hf_stock_limit = 1d-3

contains

subroutine collect_partial_waves(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("pw_bessel_k_values", test_bessel_k_values), &
    new_unittest("pw_bessel_k_far", test_bessel_k_far), &
    new_unittest("pw_hankel_scaled_values", test_hankel_scaled), &
    new_unittest("pw_round_trip", test_round_trip), &
    new_unittest("pw_eval_point_matches", test_eval_point), &
    new_unittest("pw_hams_field", test_hams_field) &
    ]
end subroutine collect_partial_waves

! ---------------------------------------------------------------------------------------------- K_m

!> K_m(x) for m = 0..10 at x = 0.1, 1 and 10 against mpmath, 1e-12 relative.
subroutine test_bessel_k_values(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: kv(0:10), worst
  integer :: ix, m
  character(len=120) :: info

  worst = 0d0
  do ix = 1, 3
    call pw_bessel_k(10, k_args(ix), kv)
    do m = 0, 10
      worst = max(worst, abs(kv(m) - k_table(m, ix)) / k_table(m, ix))
    end do
  end do
  print '(4x,a,es10.2)', "K_m against mpmath, worst relative error ", worst
  write (info, '(a,es10.3)') "K_m against mpmath, worst relative error ", worst
  call check(error, worst <= 1d-12, trim(info))
end subroutine test_bessel_k_values


!> e^x K_m(x) (pw_bessel_k_scaled) far out, where K_m itself underflows: x = 600, 2000 and 1e4, m up to 40, against mpmath to
!> 1e-12; and e^w K_m(w) at complex w (pw_bessel_k_scaled_c) to 1e-12. Below 600 the scaled values times e^-x are those
!> of pw_bessel_k (x = 0.1, 1, 10, 300, to 1e-14).
subroutine test_bessel_k_far(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: kv(0:40), ku(0:40), worst, worst_c, worst_u
  complex(8) :: kc(0:12)
  real(8), parameter :: xs(4) = [0.1d0, 1d0, 10d0, 300d0]
  integer :: ix, j, m
  character(len=160) :: info

  worst = 0d0
  do ix = 1, 3
    call pw_bessel_k_scaled(40, kfar_args(ix), kv)
    do j = 1, 5
      worst = max(worst, abs(kv(kfar_orders(j)) - kfar_table(j, ix)) / kfar_table(j, ix))
    end do
  end do
  worst_c = 0d0
  do ix = 1, 5
    call pw_bessel_k_scaled_c(12, kc_args(ix), kc)
    do j = 1, 3
      worst_c = max(worst_c, abs(kc(kc_orders(j)) - cmplx(kc_table(1, j, ix), kc_table(2, j, ix), 8)) &
                / abs(cmplx(kc_table(1, j, ix), kc_table(2, j, ix), 8)))
    end do
  end do
  worst_u = 0d0
  do ix = 1, 4
    call pw_bessel_k_scaled(10, xs(ix), kv(0:10))
    call pw_bessel_k(10, xs(ix), ku(0:10))
    do m = 0, 10
      worst_u = max(worst_u, abs(kv(m) * exp(-xs(ix)) - ku(m)) / ku(m))
    end do
  end do
  print '(4x,a,3es10.2)', "e^x K_m far out, e^w K_m complex, against pw_bessel_k: worst relative ", worst, worst_c, worst_u
  write (info, '(a,3es10.3)') "scaled K_m far out, complex, against pw_bessel_k ", worst, worst_c, worst_u
  call check(error, worst <= 1d-12 .and. worst_c <= 1d-12 .and. worst_u <= 1d-14, trim(info))
end subroutine test_bessel_k_far

!> pw_hankel_scaled against scipy (Amos) at 18 points in both half-planes, |z| from 0.7 to 3000, m = 0, 1, 7, 20, both
!> kinds, 1e-12 relative; and on the real axis against the intrinsics, H^(1) = J + iY (z = 2.1 and 300, m = 0..20).
subroutine test_hankel_scaled(error)
  type(error_type), allocatable, intent(out) :: error
  complex(8) :: hs(0:20), ref, z
  real(8) :: worst, worst_r, x
  real(8), parameter :: xr(2) = [2.1d0, 300d0]
  integer :: i, m, kind
  character(len=160) :: info

  worst = 0d0
  do i = 1, size(hankel_table, 2)
    kind = nint(hankel_table(1, i))
    m = nint(hankel_table(2, i))
    z = cmplx(hankel_table(3, i), hankel_table(4, i), 8)
    call pw_hankel_scaled(kind, 20, z, hs)
    ref = cmplx(hankel_table(5, i), hankel_table(6, i), 8)
    worst = max(worst, abs(hs(m) - ref) / abs(ref))
  end do
  worst_r = 0d0
  do i = 1, 2
    x = xr(i)
    do kind = 1, 2
      call pw_hankel_scaled(kind, 20, cmplx(x, 0d0, 8), hs)
      do m = 0, 20
        ref = cmplx(bessel_jn(m, x), merge(1d0, -1d0, kind == 1) * bessel_yn(m, x), 8) &
              * exp(cmplx(0d0, merge(-x, x, kind == 1), 8))
        worst_r = max(worst_r, abs(hs(m) - ref) / abs(ref))
      end do
    end do
  end do
  print '(4x,a,2es10.2)', "scaled Hankel functions against scipy, and on the real axis against J + iY: worst ", worst, worst_r
  write (info, '(a,2es10.3)') "scaled Hankel functions against scipy, real axis ", worst, worst_r
  call check(error, worst <= 1d-12 .and. worst_r <= 1d-12, trim(info))
end subroutine test_hankel_scaled

! ---------------------------------------------------------------------------------------------- round trip

!> H_m^(1)(x) for m = 0..mmax + 1, by the intrinsics.
subroutine hankel_all(mmax, x, hk)
  integer, intent(in) :: mmax
  real(8), intent(in) :: x
  complex(8), intent(out) :: hk(0:mmax + 1)
  integer :: m
  do m = 0, mmax + 1
    hk(m) = cmplx(bessel_jn(m, x), bessel_yn(m, x), 8)
  end do
end subroutine hankel_all

!> A field of known modes, sampled on the cylinder and decomposed, gives back its coefficients to 1e-12; and pw_eval_surface
!> reproduces phi, grad phi, phi_z and phi_zz of the same modes at r = 1.5 R on z = 0 to 1e-12. The reference gradient uses
!> other derivative identities (H_m' = H_{m-1} - m H_m / x, K_m' = -K_{m-1} - m K_m / x) than the module.
subroutine test_round_trip(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: mmax = 6, nev = 4, npt = 7
  real(8), parameter :: h = 20d0, radius = 12d0, omega = 0.9d0, centre(2) = [3d0, -2d0]
  real(8) :: wvn_local(nev + 1), k0, kev(nev), nu, zf, r, th, xk, worst_root, worst_a, x(3)
  real(8) :: kvals(0:mmax + 1), kr(0:mmax + 1, nev)
  complex(8) :: b(-mmax:mmax, 0:nev), a_true(-mmax:mmax, 0:nev), hk(0:mmax + 1), hr(0:mmax + 1)
  complex(8), allocatable :: samples(:)
  complex(8) :: phi, grad(3), phizz, phi_ref, gr_ref(3), phizz_ref, dr, dth, e, term, dterm, sgn
  real(8) :: err_phi, err_grad, err_zz, max_phi, max_grad, max_zz
  type(cs_points) :: pts
  type(pw_set) :: pw
  integer :: m, n, i, ip, am
  character(len=160) :: info
  external :: DISPERSION

  call DISPERSION(wvn_local, nev + 1, omega, h)
  k0 = wvn_local(1)
  kev = wvn_local(2:)
  nu = omega**2 / g_dispersion
  ! The roots of WVN satisfy the dispersion relations, so they can be used as the vertical modes
  worst_root = abs(k0 * tanh(k0 * h) - nu)
  do n = 1, nev
    worst_root = max(worst_root, abs(kev(n) * tan(kev(n) * h) + nu))
  end do
  write (info, '(a,es10.3)') "WVN roots off the dispersion relation, relative to omega^2/g ", worst_root / nu
  call check(error, worst_root <= 1d-12 * nu, trim(info))
  if (allocated(error)) return

  ! Coefficients at r = R, O(1), and the stored coefficients a that multiply H_m(k0 r) and K_m(kn r)
  call hankel_all(mmax, k0 * radius, hk)
  do n = 0, nev
    if (n > 0) call pw_bessel_k(mmax + 1, kev(n) * radius, kr(:, n))
    do m = -mmax, mmax
      b(m, n) = cmplx(sin(1.3d0 * (m + 7 * n) + 0.2d0), cos(0.7d0 * m - 1.1d0 * n), 8)
      am = abs(m)
      if (n == 0) then
        sgn = merge(-1d0, 1d0, m < 0 .and. mod(am, 2) == 1)
        a_true(m, n) = b(m, n) / (sgn * hk(am))
      else
        a_true(m, n) = b(m, n) / kr(am, n)
      end if
    end do
  end do

  ! Samples: phi(R, theta, z) = sum b_mn e^{i m theta} Z_n(z)
  call pw_sample_points(centre, radius, h, mmax, nev, pts)
  call check(error, pts%n == (4 * mmax + 4) * 2 * (nev + 1), "wrong number of sample points")
  if (allocated(error)) return
  call check(error, all(pts%kind == CS_WALL), "sample points are not CS_WALL")
  if (allocated(error)) return
  allocate (samples(pts%n))
  do i = 1, pts%n
    x = pts%x(:, i)
    call check(error, abs(norm2(x(1:2) - centre) - radius) <= 1d-12 * radius .and. x(3) < 0d0 .and. x(3) > -h, &
               "sample point not on the cylinder")
    if (allocated(error)) return
    call check(error, abs(sum(pts%nrm(1:2, i) * (x(1:2) - centre)) / radius - 1d0) <= 1d-12, "normal not outward")
    if (allocated(error)) return
    th = atan2(x(2) - centre(2), x(1) - centre(1))
    samples(i) = (0d0, 0d0)
    do n = 0, nev
      if (n == 0) then
        zf = (exp(k0 * x(3)) + exp(-k0 * (x(3) + 2d0 * h))) / (1d0 + exp(-2d0 * k0 * h))
      else
        zf = cos(kev(n) * (x(3) + h)) / cos(kev(n) * h)
      end if
      do m = -mmax, mmax
        samples(i) = samples(i) + b(m, n) * exp(ci * m * th) * zf
      end do
    end do
  end do

  call pw_decompose(samples, centre, radius, h, k0, kev, mmax, nev, pw)
  worst_a = maxval(abs(pw%a - a_true) / abs(a_true))
  print '(4x,a,es10.2)', "round trip: worst relative error of the coefficients ", worst_a
  write (info, '(a,es10.3)') "coefficients not recovered, worst relative error ", worst_a
  call check(error, worst_a <= 1d-12, trim(info))
  if (allocated(error)) return

  ! Reconstruction at r = 1.5 R on z = 0
  err_phi = 0d0; err_grad = 0d0; err_zz = 0d0; max_phi = 0d0; max_grad = 0d0; max_zz = 0d0
  r = 1.5d0 * radius
  do ip = 1, npt
    th = 2d0 * pi * (ip - 0.37d0) / npt
    x = [centre(1) + r * cos(th), centre(2) + r * sin(th), 0d0]
    call pw_eval_surface(pw, x, phi, grad, phizz)
    call hankel_all(mmax, k0 * r, hr)
    phi_ref = (0d0, 0d0); dr = (0d0, 0d0); dth = (0d0, 0d0); phizz_ref = (0d0, 0d0); gr_ref(3) = (0d0, 0d0)
    do n = 0, nev
      if (n > 0) then
        xk = kev(n) * r
        call pw_bessel_k(mmax + 1, xk, kvals)
      end if
      do m = -mmax, mmax
        am = abs(m)
        e = exp(ci * m * th)
        if (n == 0) then
          sgn = merge(-1d0, 1d0, m < 0 .and. mod(am, 2) == 1)
          term = sgn * hr(am)
          if (am == 0) then
            dterm = -sgn * hr(1) * k0
          else
            dterm = sgn * (hr(am - 1) - am * hr(am) / (k0 * r)) * k0
          end if
          gr_ref(3) = gr_ref(3) + a_true(m, n) * e * term * k0 * tanh(k0 * h)
          phizz_ref = phizz_ref + a_true(m, n) * e * term * k0**2
        else
          term = kvals(am)
          if (am == 0) then
            dterm = -kvals(1) * kev(n)
          else
            dterm = (-kvals(am - 1) - am * kvals(am) / xk) * kev(n)
          end if
          gr_ref(3) = gr_ref(3) - a_true(m, n) * e * term * kev(n) * tan(kev(n) * h)
          phizz_ref = phizz_ref - a_true(m, n) * e * term * kev(n)**2
        end if
        phi_ref = phi_ref + a_true(m, n) * e * term
        dr = dr + a_true(m, n) * e * dterm
        dth = dth + a_true(m, n) * e * term * ci * m
      end do
    end do
    gr_ref(1) = cos(th) * dr - sin(th) * dth / r
    gr_ref(2) = sin(th) * dr + cos(th) * dth / r
    max_phi = max(max_phi, abs(phi_ref))
    max_grad = max(max_grad, norm2(abs(gr_ref)))
    max_zz = max(max_zz, abs(phizz_ref))
    err_phi = max(err_phi, abs(phi - phi_ref))
    err_grad = max(err_grad, norm2(abs(grad - gr_ref)))
    err_zz = max(err_zz, abs(phizz - phizz_ref))
  end do
  err_phi = err_phi / max_phi
  err_grad = err_grad / max_grad
  err_zz = err_zz / max_zz
  print '(4x,a,3es10.2)', "round trip at r = 1.5 R: phi, grad phi, phi_zz errors relative to the largest ", &
    err_phi, err_grad, err_zz
  write (info, '(a,3es10.3)') "reconstruction at r = 1.5 R off (phi, grad, phi_zz) ", err_phi, err_grad, err_zz
  call check(error, err_phi <= 1d-12 .and. err_grad <= 1d-12 .and. err_zz <= 1d-12, trim(info))
end subroutine test_round_trip

! ---------------------------------------------------------------------------------------------- pw_eval_point

!> pw_eval_point against pw_eval_surface on z = 0 (1e-14 relative), and against an independent sum of the modes at depth,
!> at the circle samples (with and without the radial functions of pw_radial_at) and at r = 1.4 R (1e-12). Also
!> pw_eval_surface with the radial functions of pw_radial_at equals it without them, to the bit (the same arithmetic).
!> A pw_set of known coefficients: mmax 5, nev 3, h 18, R 9, centre (-1, 2), omega 1.0.
subroutine test_eval_point(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: mmax = 5, nev = 3, npt = 9
  real(8), parameter :: h = 18d0, radius = 9d0, omega = 1.0d0, centre(2) = [-1d0, 2d0]
  real(8) :: wvn_local(nev + 1), x(3), th, r, zf, kvals(0:mmax + 1), worst_s, worst_d, worst_c, scale
  complex(8) :: hk(0:mmax + 1), phi, phi_s, grad(3), phizz, phi_c, grad_c(3), phizz_c, ref, sgn
  type(pw_set) :: pw
  type(pw_radial) :: rad
  type(cs_points) :: pts
  integer :: m, n, ip, i, am
  character(len=160) :: info
  external :: DISPERSION

  call DISPERSION(wvn_local, nev + 1, omega, h)
  pw%mmax = mmax
  pw%nev = nev
  pw%k0 = wvn_local(1)
  pw%h = h
  pw%radius = radius
  pw%centre = centre
  allocate (pw%kev(nev), pw%a(-mmax:mmax, 0:nev))
  pw%kev = wvn_local(2:)
  do n = 0, nev
    do m = -mmax, mmax
      pw%a(m, n) = cmplx(cos(0.8d0 * m + 2.1d0 * n + 0.4d0), sin(1.9d0 * m - 0.6d0 * n), 8)
    end do
  end do

  ! z = 0: pw_eval_point is the phi of pw_eval_surface; with the radial functions given, pw_eval_surface is unchanged
  worst_s = 0d0; worst_c = 0d0; scale = 0d0
  do ip = 1, npt
    th = 2d0 * pi * (ip - 0.41d0) / npt
    r = radius * (1d0 + 0.3d0 * ip)
    x = [centre(1) + r * cos(th), centre(2) + r * sin(th), 0d0]
    call pw_eval_point(pw, x, phi)
    call pw_eval_surface(pw, x, phi_s, grad, phizz)
    call pw_radial_at(pw, r, rad)
    call pw_eval_surface(pw, x, phi_c, grad_c, phizz_c, rad=rad)
    scale = max(scale, abs(phi_s))
    worst_s = max(worst_s, abs(phi - phi_s))
    if (phi_c /= phi_s .or. any(grad_c /= grad) .or. phizz_c /= phizz) worst_c = 1d0
  end do
  print '(4x,a,es10.2)', "pw_eval_point against pw_eval_surface on z = 0, relative ", worst_s / scale
  write (info, '(a,es10.3)') "pw_eval_point differs from pw_eval_surface on z = 0, relative ", worst_s / scale
  call check(error, worst_s <= 1d-14 * scale, trim(info))
  if (allocated(error)) return
  call check(error, worst_c == 0d0, "pw_eval_surface with the radial functions of pw_radial_at differs from without")
  if (allocated(error)) return

  ! At depth: the circle samples (r = R, with and without the radial functions) and r = 1.4 R, against a direct sum
  call pw_sample_points(centre, radius, h, mmax, nev, pts)
  call pw_radial_at(pw, radius, rad)
  worst_d = 0d0; scale = 0d0
  do i = 1, pts%n + npt
    if (i <= pts%n) then
      x = pts%x(:, i)
    else
      th = 2d0 * pi * (i - pts%n - 0.2d0) / npt
      x = [centre(1) + 1.4d0 * radius * cos(th), centre(2) + 1.4d0 * radius * sin(th), -h * (i - pts%n) / (npt + 1d0)]
    end if
    r = norm2(x(1:2) - centre)
    th = atan2(x(2) - centre(2), x(1) - centre(1))
    call hankel_all(mmax, pw%k0 * r, hk)
    ref = (0d0, 0d0)
    do n = 0, nev
      if (n == 0) then
        zf = cosh(pw%k0 * (x(3) + h)) / cosh(pw%k0 * h)
      else
        zf = cos(pw%kev(n) * (x(3) + h)) / cos(pw%kev(n) * h)
        call pw_bessel_k(mmax + 1, pw%kev(n) * r, kvals)
      end if
      do m = -mmax, mmax
        am = abs(m)
        if (n == 0) then
          sgn = merge(-1d0, 1d0, m < 0 .and. mod(am, 2) == 1)
          ref = ref + pw%a(m, n) * exp(ci * m * th) * sgn * hk(am) * zf
        else
          ref = ref + pw%a(m, n) * exp(ci * m * th) * kvals(am) * zf
        end if
      end do
    end do
    call pw_eval_point(pw, x, phi)
    scale = max(scale, abs(ref))
    worst_d = max(worst_d, abs(phi - ref))
    if (i <= pts%n) then
      call pw_eval_point(pw, x, phi_c, rad=rad)
      worst_d = max(worst_d, abs(phi_c - ref))
    end if
  end do
  print '(4x,a,es10.2)', "pw_eval_point at depth against the direct sum of the modes, relative ", worst_d / scale
  write (info, '(a,es10.3)') "pw_eval_point at depth differs from the direct sum of the modes, relative ", worst_d / scale
  call check(error, worst_d <= 1d-12 * scale, trim(info))
end subroutine test_eval_point

! ---------------------------------------------------------------------------------------------- HAMS field

!> The scattered field of one state on the 448-panel sharp cylinder, by fields_green: sampled on r = R = 15 m, decomposed,
!> and evaluated at r = 20 m and 30 m on z = 0, against fields_green there. Hull data: the MacCamy-Fuchs scattered wave of
!> depth 15 m on the hull panels (maccamy_fuchs_state, with amp = 0 so that only the panel integrals enter, as in
!> green_perturbation_only). The unit tests have no BEM solve; whatever the hull data, a sum of Green functions over the
!> panels satisfies Laplace, the free-surface and bottom conditions and the radiation condition outside the hull, which is
!> all the decomposition assumes.
!> The gate (1e-4) is on that Green-function sum: fields_green is called with centroid_only, so that every panel takes its
!> centroid value (near_panel = .false. of field_influence_row), the route by which the driver samples the circle.
!> With the stock route, panels within 50 psize take the panel integrals of 1/r and of its image in z = 0 in place of
!> their centroid values (field_influence_row). That correction satisfies phi_z = 0 on z = 0, not the free-surface
!> condition, so it is not a sum of the modes; at r = R it is extrapolated by the propagating modes. Measured with the
!> stock route: 4.9e-4 (phi) and 3.8e-4 (gradient), the same for mmax = 16..32 and nev = 8..16; the correction alone is
!> 2.5e-4 at the rings. That route is checked against 1e-3, the level at which the mode check warns.
subroutine test_hams_field(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, mesh
  type(first_order_state) :: states(1)
  real(8) :: err_phi, err_grad
  integer :: nthread_saved
  character(len=160) :: info

  call set_wave(hf_depth)
  nthread_saved = NTHREAD
  NTHREAD = 4
  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call load_panels(full, 0, 0, mesh)
  call check(error, mesh%n == 448, "the fixture is not the 448-panel cylinder")
  if (allocated(error)) return
  call maccamy_fuchs_state(mesh, 10d0, hf_depth, WK, W1, 1d0, BETA, wave_gravity(W1, WK, hf_depth), 1000d0, states(1))
  states(1)%amp = 0d0

  call hams_field_errors(mesh, states, .true., err_phi, err_grad, error)
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if
  print '(4x,a,i0,a,i0,a,2es10.2)', "partial waves (mmax = ", hf_mmax, ", nev = ", hf_nev, &
    ") against the Green-function sum at r = 20, 30 m: phi and horizontal gradient errors, relative to the largest ", &
    err_phi, err_grad
  write (info, '(a,es10.3)') "phi against fields_green (centroid values), relative to max |phi_P| ", err_phi
  call check(error, err_phi <= hf_phi_limit, trim(info))
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if
  write (info, '(a,es10.3)') "horizontal gradient against fields_green (centroid values), relative to the largest ", err_grad
  call check(error, err_grad <= hf_grad_limit, trim(info))
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if

  call hams_field_errors(mesh, states, .false., err_phi, err_grad, error)
  NTHREAD = nthread_saved
  if (allocated(error)) return
  print '(4x,a,2es10.2)', "the same with the near-panel integrals of fields_green (stock route): phi and gradient errors ", &
    err_phi, err_grad
  write (info, '(a,2es10.3)') "partial waves against fields_green with near-panel integrals (phi, gradient) ", &
    err_phi, err_grad
  call check(error, err_phi <= hf_stock_limit .and. err_grad <= hf_stock_limit, trim(info))
end subroutine test_hams_field

!> pw_hams_field for one route: sample by fields_green on r = R, decompose, and compare pw_eval_surface with fields_green at
!> 16 points on each of r = 20 m and 30 m, z = 0. Errors relative to the largest |phi_P| and |grad_h phi_P| there.
!> centroid: both fields_green calls with centroid_only (every panel at its centroid value), else the stock route.
subroutine hams_field_errors(mesh, states, centroid, err_phi, err_grad, error)
  type(panel_mesh), intent(in) :: mesh
  type(first_order_state), intent(in) :: states(1)
  logical, intent(in) :: centroid
  real(8), intent(out) :: err_phi, err_grad
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nring = 16
  real(8), parameter :: rings(2) = [20d0, 30d0], step = 0.02d0
  type(cs_points) :: pts, ring
  type(cs_fields) :: fields(1), ref(1)
  type(pw_set) :: pw
  real(8) :: kev(hf_nev), th, max_phi, max_grad
  complex(8) :: phi, grad(3), phizz
  integer :: i, l, j

  ! Sample on the cylinder; fields_green leaves phi zero at CS_WALL points below the surface, so it is given the points as
  ! CS_FREE (the kind only decides whether phi is computed)
  call pw_sample_points([0d0, 0d0], hf_radius, hf_depth, hf_mmax, hf_nev, pts)
  pts%kind = CS_FREE
  call fields_green(mesh, pts, states, step, .false., fields, centroid_only=centroid)
  call pw_evanescent_from_wvn(hf_nev, kev)
  call check(error, all(kev == WVN(2:hf_nev + 1)), "pw_evanescent_from_wvn does not return WVN(2:nev+1)")
  if (allocated(error)) return
  call pw_decompose(fields(1)%phi, [0d0, 0d0], hf_radius, hf_depth, WK, kev, hf_mmax, hf_nev, pw)

  ! The reference at r = 20 m and 30 m on z = 0
  ring%n = 2 * nring
  allocate (ring%x(3, ring%n), ring%nrm(3, ring%n), ring%w(ring%n), ring%psize(ring%n), ring%kind(ring%n), ring%body(ring%n))
  ring%w = 1d0
  ring%psize = 0.5d0
  ring%kind = CS_LINE
  ring%body = 0
  i = 0
  do l = 1, 2
    do j = 1, nring
      i = i + 1
      th = 2d0 * pi * (j - 0.3d0) / nring
      ring%x(:, i) = [rings(l) * cos(th), rings(l) * sin(th), 0d0]
      ring%nrm(:, i) = [cos(th), sin(th), 0d0]
    end do
  end do
  call fields_green(mesh, ring, states, step, .false., ref, centroid_only=centroid)

  max_phi = 0d0; max_grad = 0d0; err_phi = 0d0; err_grad = 0d0
  do i = 1, ring%n
    call pw_eval_surface(pw, ring%x(:, i), phi, grad, phizz)
    max_phi = max(max_phi, abs(ref(1)%phi(i)))
    max_grad = max(max_grad, norm2(abs(ref(1)%grad(1:2, i))))
    err_phi = max(err_phi, abs(phi - ref(1)%phi(i)))
    err_grad = max(err_grad, norm2(abs(grad(1:2) - ref(1)%grad(1:2, i))))
  end do
  err_phi = err_phi / max_phi
  err_grad = err_grad / max_grad
end subroutine hams_field_errors

end module test_partial_waves
