module dcpm
use input
use constants
use ieee_arithmetic
use cm1libs, only: rslf, rsif

#ifdef MPI
    use mpi
#endif
implicit none

    private
    ! constants
    real (kind=dp), parameter :: C_p = 1004.0d0 !water specific heat at constant pressure J K-1 kg-1ke
    real (kind=dp), parameter :: L = 2.25d6 !latent heat of water condensation J kg-1
    real (kind=dp), parameter :: L_ice = 2.8d6 !latent heat of ice sublimation J kg-1
    real (kind=dp), parameter :: rho_w = 1.0d3 !water density kg m-3
    real (kind=dp), parameter :: R_d = 287.d0 !air specific gas constant J K-1 kg-1
    real (kind=dp), parameter :: R_v = 461.5d0 !water vapor gas constant J K-1 kg-1
    real (kind=dp), parameter :: R = 8.314d0 !molar gas constant J K-1 mol-1
    real (kind=dp), parameter :: M_w = 0.018d0 !water molecular weight kg mol-1
    real (kind=dp), parameter :: M_a = 0.0289d0 !air molecular weight kg mol-1
    real (kind=dp), parameter :: coeff_thermal_conductivity_air = 0.96d0 !thermal accommodation coefficient of air
    real (kind=dp), parameter :: coeff_condensation_diffusivity = 0.042d0 ! condensation coefficient to calculate diffusivity of water vapor in air
    real (kind=dp), parameter :: MwMa = 0.622d0 !ratio of molecular weight between water and air Mw/Ma
    real (kind=dp), parameter :: kappa = 0.286d0 !Rd/Cp in equation T2 = T1(P2/P1)^(Rd/Cp)
    real (kind=dp), parameter :: gamma = 0.5d0
    logical, parameter :: if_dynamic = .false. ! if_dynamic; whether to use dynamic bin (bin radius can change)
    logical, parameter :: if_coll = .true. ! if_coll: whether to consider collision 
    logical, parameter :: if_turb = .true. ! if_turb: whether to consider turbulence
    logical, parameter :: if_square_radius_scheme = .true. ! if_square_radius_scheme: whether use derivative of r^2 or r to solve condensation/evaporation radius change
    logical, parameter :: if_limit_t_step_cond_by_radius = .true. ! if_limit_t_step_cond_by_radius: whether to limit condensation time step by radius
    logical, parameter :: if_limit_t_step_cond_by_bin = .false. ! more rigorous limitation: the bin radius is not changed beyond nearby bins, require logical, parameter :: if_limit_t_step_cond_by_radius = .true.
    logical, parameter :: if_limit_t_step_cond_by_temp = .true. ! if_limit_t_step_cond_by_temp: whether to limit condensation time step by temperature
    logical, parameter :: if_limit_t_step_cond_by_moisture = .true. ! if dwl > 0, dwl should be no larger than wl*percent_change_maximum; if dwl < 0, -dwl should not be larger than wv*percent_change_maximum
    logical, parameter :: if_limit_t_step_cond_by_supersaturation = .true. ! in addition to limit by moisture, dwl can't exceeds (wv-wv_sat), require if_limit_t_step_cond_by_moisture = .true. 
    ! logical, parameter :: if_redistribute_aerosol = .false. ! If redistribute dry aerosol larger than max CCN radius to the initial aerosol spectrum
    ! logical, parameter :: if_correct_aerosol = (if_redistribute_aerosol .eqv. .false.) ! If correct dry radius that caused by numerical error
    logical, parameter :: if_aerosol_replenish = .true. ! if replenish aerosol at bottom boundary ! For short period simulation
    logical, parameter :: if_warning_tolerance = .false. ! Whether output warning message for moisture
    logical, parameter :: if_warning_tolerance_aerosol = .false. ! Whether output warning message for aerosol
    logical, parameter :: if_warning_wv_wl_balance = .false. ! Whether output warning message for wl, wv imbalance
    logical, parameter :: if_warning_timestep = .false. ! Whether output warning message for time step breaking limitation rules
    logical, parameter :: if_warning_general = .false. !General warning switch
    integer, parameter :: iv_fall = 1 ! 1=> Best 1950 2=> Lhermitte 1990
    ! iv_fall: selection of fall velocity
        ! iv_fall = 1 => Best falling velocities (Bs50)
    integer, parameter :: ieff_coal = 4 
    ! ieff_coal: selection of coalescence efficiency
        ! ieff_coal = 1 => constant
        ! ieff_coal = 2 => E=(BO1984)
        ! ieff_coal = 3 => E=a[1+ds/dl]^(-2).exp[-b.sig.Et^2/Sc](LL82a)
        ! ieff_coal = 4 => E=f(LL82a,BO95) (Seifert at al. 2005)
          
    real (kind=dp), parameter :: t_step_cond_min = 1d-2 ! minimum time step for condensation (s)
    real (kind=dp), parameter :: t_step_cond_max = 1d-1 ! maximum time step for condensation (s)
    real (kind=dp), parameter :: t_step_coll_min = 1d-1 ! time step for coalescence (s)
    real (kind=dp), parameter :: t_step_coll_max = 1d6 ! time step for coalescence (s) (not a limit anymore)
    real (kind=dp), parameter :: percent_radius_change_maximum = 0.5 ! maximum change percent of radius during a time step (from 0 to 1)
    real (kind=dp), parameter :: percent_wl_change_maximum = 0.5 ! maximum change percent of wl during a time step (from 0 to 1)
    real (kind=dp), parameter :: temperature_change_maximum = 0.1 ! maximum change of temperature during a time step
    integer, parameter :: max_iter_cond = 2 ! maximum iteration of rk during condensation

   
    ! aerosol and drop distrubution parameters
    real (kind=dp) :: aer_kappa ! aer_kappa: hygroscopicity (will define according to aerosol type)
    real (kind=dp), parameter :: a_cond = 0.042d0 ! condensation coefficient
    real (kind=dp), parameter :: ccn_radius_min = 0.02d-6 ! minimum cloud condensation nuclei (CCN) radius [m]
    real (kind=dp), parameter :: ccn_radius_max = 5d-6 ! maximum cloud condensation nuclei (CCN) radius [m]
    real (kind=dp), parameter :: drop_radius_min = 50d-6 ! minumum drop radius [m]
    real (kind=dp), parameter :: drop_radius_max = 4000d-6 ! maximum drop radius [m]
    real (kind=dp), parameter :: drop_brkp_parent_min = 250d-6 ! minumum parent drop radius during breakup [m]
    real (kind=dp), parameter :: drop_brkp_fragment_min = 10d-6 ! minimum fragment drop radius during breakup [m]
    ! coalescence threshold deprecated
    real (kind=dp), parameter :: drop_coal_large_min = 10d-6 ! minimum coalescence large drop radius [m]
    real (kind=dp), parameter :: drop_coal_small_min = 10d-6 ! minimum coalescence small drop radius [m]
    real (kind=dp) :: ratio_wet_coll_min = 1 ! minimum radius ratio (r_wet/r_dry) for collision (reinitialized in dcpm_init)
    real (kind=dp), dimension(:), allocatable :: aer_num_concs, aer_gmds, aer_gsds
    ! real (kind=dp), dimension(4) :: aer_num_concs = (/ 1401.9d6, 415.7d6, 0.3d6, 0.3d6 /)
    ! real (kind=dp), dimension(4) :: aer_gmds = (/ 0.076d-6, 0.195d-6, 0.75d-6, 2.2d-6 /)
    ! real (kind=dp), dimension(4) :: aer_gsds = (/ 1.63d0, 1.35d0, 1.30d0, 1.40d0 /)
    real (kind=dp) :: ratio_bins
    real (kind=dp), dimension(11) :: qc_keys = (/ 1e-7, 5e-7, 1e-06, 5e-06, 1e-05, 5e-05, 0.0001, 0.0005, 0.001,  0.005, 0.010 /)
    real (kind=dp), dimension(9) :: hdf_keys = (/ 1e-4, 1e-3,1e-2,5e-2,0.1,0.25,0.5,0.75,1.0 /)

    real (kind=sp), dimension(:,:,:), allocatable :: qcbin_lookup_tbl, ncbin_lookup_tbl

    ! coalescence paramters
    real (kind=dp), parameter :: const_eff_coal = 0.157d0 ! coalescence efficiency as a constant (used only when selecting constant option)
    real (kind=dp), parameter :: effi_etlmt = 5d-6 ! indicater limit to switch coalescence efficient calculation

    ! ice process
    real (kind=dp), parameter :: temp_ice = -20+273.15
    real (kind=dp), parameter :: temp_water = -5+273.15

    ! tolerance in numerical comparsion (i.e. if very small number is zero)
    real (kind=dp), parameter :: tolerance_t = 1d-6 ! time tolerance (for comparsion only)
    real (kind=dp), parameter :: tolerance_temp = 0.01 ! temperature tolerance
    ! real (kind=dp), parameter :: tolerance_wv = 1d-7 ! water vapor mixing ratio tolerance
    ! real (kind=dp), parameter :: tolerance_wl = 1d-7 ! liquid water mixing ratio tolerance
    ! real (kind=dp), parameter :: tolerance_wi = 1d-7 ! ice water mixing ratio tolerance
    real (kind=dp), parameter :: tolerance_bin_num_conc = 1 ! bin number concentration tolerance
    real (kind=dp), parameter :: threshold_bin_num_conc_coll = tolerance_bin_num_conc*10! bin number concentration threshold
    real (kind=dp), parameter :: tolerance_bin_r = 1e-12 ! bin radius concentration tolerance (compare with 0 only)
    real (kind=dp), parameter :: tolerance_r_ratio = 1e-2 ! radius 0.01 off is fine (compare two radius)
    real (kind=dp), parameter :: scale_bin_wl = 1e9
    ! real (kind=dp), parameter :: threshold_wl_coll = 1e-4 ! threshold of total liquid water mixing ratio for collision
    real (kind=dp) :: tolerance_bin_wl ! tolerance of bin liquid water mixing ratio
    real (kind=dp) :: threshold_bin_wl_coll ! threshold of bin liquid water mixing ratio for collision
    integer, parameter :: nk_kernel = 10
    integer, parameter :: nk_brkp_kernel = 1
    real (kind=dp), dimension(21,21):: turbulence_kernel
    ! turbulence_kernel: turbulence kernels for collision
    real (kind=dp), dimension(40) :: rd1
    real (kind=dp), dimension(22) :: rd2
    real (kind=dp), dimension(40,40):: collision_kernel_0000_1, collision_kernel_2500_1, collision_kernel_5000_1
    real (kind=dp), dimension(40,22):: collision_kernel_0000_2, collision_kernel_2500_2, collision_kernel_5000_2
    real (kind=dp), dimension(22,22):: collision_kernel_0000_3, collision_kernel_2500_3, collision_kernel_5000_3
    real (kind=dp), dimension(:), allocatable :: bin_r_wet, bin_r_wet_check, bin_r_dry_init_vol_weight, bin_num_conc_base
    real (kind=dp), dimension(:), allocatable ::  bin_r_wet_grid, bin_r_wet_grid_check
    ! real (kind=dp), dimension(num_bins) :: bin_r_wet
    ! real (kind=dp), dimension(num_bins+1) :: bin_r_wet_grid

    ! Coalescence kernels and weights
    ! integer :: nk_kernel
    real (kind=dp), dimension(:), allocatable :: z_kernels, z_brkp_kernels
    real (kind=dp), dimension(:,:,:), allocatable :: kernels_coal, kernels_brkp
    real (kind=dp), dimension(:,:,:,:), allocatable :: weights_brkp
    real (kind=dp), dimension(:,:,:), allocatable :: weights_coal
    logical :: if_dcpm_initialized = .false.
    logical :: if_dcpm_activated_myid = .false.
    logical :: if_dcpm_activated = .false.
    integer :: ind_bin_base_empty_start = 0


    ! timer
    integer :: t_clock_init, t_clock_curr, clock_rate_dcpm, clock_max_dcpm
    real (kind=sp) :: t_cpu_init, t_cpu_curr

    real (kind=dp), dimension(:,:,:,:), allocatable :: q3d_dcpm
 
    public :: dcpm_calculate, dcpm_init, dcpm_fallout, apply_DSD_bdy, correct_wl_dcpm
    public :: if_dcpm_activated, bin_r_wet, scale_bin_wl

    contains

    character(len=20) function str(k)
    !   "Convert an integer to string."
        integer, intent(in) :: k
        write (str, *) k
        str = adjustl(str)
    end function str

    subroutine correct_wl_dcpm(q3d)
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(inout) :: q3d
        real (kind=dp), dimension(num_bins) :: bin_wl
        real (kind=dp) :: wl
        integer :: i3d, j3d, k3d
        do i3d=1,ni
        do j3d=1,nj
        do k3d=1,nk
            bin_wl = q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl
            wl = sum(bin_wl)
            q3d(i3d,j3d,k3d,nql1) = wl
        end do
        end do 
        end do
    end subroutine correct_wl_dcpm

    subroutine dcpm_calculate (dt, q3d, dcpm3d, adddiag3d, p3d, th0, th3d, pi0, pp3d, zh, zf, temp3d, rho3d, dorad, effc, effi)
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke), intent(in) :: pi0,th0,zh
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke), intent(inout) :: pp3d,th3d,p3d,temp3d, rho3d
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(inout) :: q3d
        real (kind=sp), dimension(ibr:ier,jbr:jer,kbr:ker), intent(inout) :: effc ! EFFC - DROPLET EFFECTIVE RADIUS
        real (kind=sp), dimension(ibr:ier,jbr:jer,kbr:ker), intent(inout) :: effi ! EFFI - ICE EFFECTIVE RADIUS
        real (kind=dp) :: effi_curr, effc_curr
        real, dimension(ibm:iem,jbm:jem,kbm:kem,num_bins+1), intent(inout) :: dcpm3d
        real, dimension(ibm:iem,jbm:jem,kbm:kem,nadddiagvars3d), intent(inout) :: adddiag3d
        logical, intent(in) :: dorad
        real (kind=dp), dimension(num_bins) :: bin_r_dry, bin_wl, bin_r_dry_temp, bin_r_wet_temp, bin_r_wet_predict, bin_shift_curr, bin_num_conc, bin_dnum_conc_dt, bin_dwl_dt, bin_deriv, k1, k2, k3, k4
        real (kind=dp), dimension(num_bins) :: bin_r_dry_bc, bin_num_conc_bc
        logical, dimension(num_bins) :: bin_if_act ! if the bin aerosol is activated
        real (kind=dp), dimension(num_bins) :: bin_weight, bin_num_conc_base_curr, bin_num_conc_preadv, bin_num_conc_temp, bin_r_dry_vol_defecit, bin_r_dry_vol_excess
        real (kind=sp), intent(in) :: dt
        real (kind=dp), dimension(num_bins, num_bins) :: kernels_coal_curr, kernels_brkp_curr
        real (kind=dp), dimension(num_bins, num_bins, num_bins) :: weights_brkp_curr, weights_coal_curr
        integer :: i3d, j3d, k3d, ind_bin, ind_bin_i, ind_bin_j, ind_bin_k, ind_zk, num_cond_bins, num_coal_bins, ind_bin_radius_control
        real (kind=dp) :: r_dry_cubic, wl_bin_defecit, wl_bin_excess, r_wet_cubic
        real (kind=dp) :: wl_bin, supersat, supersat_crit, p, r_wet_temp, r_crit, r_shift, z, mass_change_ratio, radius_change_max
        real (kind=dp) :: wv, wv_init, wv_predict, wv_sat, wv_sat_predict, wl, wl_init, wl_predict, temp, temp_init, temp_predict
        logical :: if_do_cond, if_did_cond_curr, if_did_cond_stat, if_act, if_did_coll_stat, if_grid_counted, if_debug_cond
        real (kind=dp) :: t_dcpm, t_cond, t_step_dcpm, t_step_cond, t_step_cond_max_curr, t_step_cond_radius, t_step_cond_temperature, t_step_cond_moisture, t_step_coll, t_step_coll_suggest, t_step_cond_min_stat, t_coll,   t_step_coll_min_negc
        real (kind=dp) :: t_step_cond_ratio_radius, t_step_cond_ratio_temperature, t_step_cond_ratio_moisture
        real (kind=dp) :: thermal_conductivity_air, diffusivity_water_vapor, G_a, G_b, coeff_growth, dr_dt, dwl_dt, r_dry, r_wet, rho_d, dwl, dnum_conc_dt
        real (kind=dp) :: numerator, denominator
        real (kind=dp) :: wi, wv_sat_i, wv_sat_w, portion_ice
        real (kind=dp) :: wl_max, wl_min, wi_max, wi_min, supersat_max, supersat_min, num_conc_max, num_conc_min, r_drop_max
        real (kind=dp) :: wl_max_all, wl_min_all, wi_max_all, wi_min_all, supersat_max_all, supersat_min_all, num_conc_max_all, num_conc_min_all, r_drop_max_all, t_step_coll_min_negc_all
        integer :: num_grid_active_stat, num_rk_stat, num_dist_stat, num_rk_stat_all, num_dist_stat_all
        integer :: count, num_iter
        real (kind=dp) :: temp_guess, temp_result, wv_result, wl_result, wi_result, wtot
        real (kind=dp) :: temp_error
        integer :: count_loop, ind_bin_start_redistribution
        real (kind=dp) :: hdf
        logical :: if_grid_wet, if_r_dry_vol_defecit, if_r_dry_vol_excess
        logical :: cond_wl, cond_dry, cond_radius, cond_supersat, if_start_redistribution
        real (kind=dp) :: vol_loss, num_conc_increase, r_dry_vol_defecit_tot, r_dry_vol_excess_tot, scale_r_dry_vol_tot, r_dry_vol_tolerance, r_dry_vol_diff
        integer, dimension(num_bins) :: bin_ind_coal
        integer, dimension(num_bins) :: bin_ind_brkp
        integer :: num_coal_bins_curr, num_brkp_bins_curr, ind_ind_bin_i, ind_ind_bin_j
        real (kind=dp) :: num_conc_coal, num_conc_brkp, wl_coal_i, wl_coal_j, wl_brkp_i, wl_brkp_j
        real :: r_dry_max_possible
        logical :: if_r_dry_adjusted = .false.
        logical :: if_evaporation, if_condensation, if_coal, if_brkp
        
        ! real (kind=dp), dimension(num_bins, num_bins, num_bins) :: vol_transfer_ijk
        ! real (kind=dp), dimension(num_bins, num_bins) :: vol_loss_ij


        ! call system ("touch -a "//'Output'//trim(str(myid))//'.txt')
        ! open(unit=240311, file='Output'//trim(str(myid))//'.txt', Access = 'append',Status='old')

        call correct_wl_dcpm(q3d)

        if ((ib.ne.ibm) .or. (jb.ne.jbm) .or. (kb.ne.kbm) ) then
            print *, "ERROR: ib, ibm,  not equal"
            print *, ib, ibm,  jb, jbm, kb, kbm
            call stopcm1
        end if
        if ((ie.ne.iem) .or.  (je.ne.jem) .or.  (ke.ne.kem) ) then
            print *, "ERROR: ib, ibm, not equal"
            print *, ie, iem, ier, je, jem, jer, ke, kem, ker
            call stopcm1
        end if
        
        ! if (myid .eq. 266) then
        !     i3d = 10
        !     j3d = 3
        !     print *, p3d(i3d,j3d,1:nk)
        !     print *, "1st level height:", zh(i3d,j3d,1), "west:", zh(i3d-1,j3d,1), "east:", zh(i3d+1,j3d,1), "north:", zh(i3d,j3d+1,1), "south:", zh(i3d,j3d-1,1)
        ! end if
        call dcpm_preprocess (q3d, dcpm3d, th0, th3d, pi0, pp3d)

        wl_max = maxval(q3d(1:ni,1:nj,1:nk,nql1))
        wl_min = minval(q3d(1:ni,1:nj,1:nk,nql1))
        wi_max = maxval(q3d(1:ni,1:nj,1:nk,nqs1))
        wi_min = minval(q3d(1:ni,1:nj,1:nk,nqs1))
        ! supersat_max = maxval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1))
        ! supersat_min = minval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1))
        num_conc_max = maxval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))
        num_conc_min = minval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))
        wl_max_all = wl_max
        wl_min_all = wl_min
        ! supersat_max_all = supersat_max
        ! supersat_min_all = supersat_min
        num_conc_max_all = num_conc_max
        num_conc_min_all = num_conc_min
#ifdef MPI
    call MPI_REDUCE(wl_max,wl_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wl_min,wl_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wi_max,wi_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wi_min,wi_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    ! call MPI_REDUCE(supersat_max,supersat_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    ! call MPI_REDUCE(supersat_min,supersat_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_conc_max,num_conc_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_conc_min,num_conc_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
#endif    
        if (myid .eq. 0) then
            print *, "Enter DCPM dt = ", dt, "s"    
            print *, "Before microphysics:" 
            print *, "begin variable", " wl (g/kg) ", " min ", wl_min_all*1e3," max ", wl_max_all*1e3, " end variable"
            print *, "begin variable", " number concentration (cm^{-3}) ", " min ", num_conc_min_all*1e-6," max ", num_conc_max_all*1e-6, " end variable"


        end if
        ! print *,    "extreme values before bin microphysics:", "temp", minval(temp3d(1:ni,1:nj,1:nk)), maxval(temp3d(1:ni,1:nj,1:nk)),&
        !             "wv", minval(q3d(1:ni,1:nj,1:nk,nqv)), maxval(q3d(1:ni,1:nj,1:nk,nqv)), "wl", minval(q3d(1:ni,1:nj,1:nk,nql1)), maxval(q3d(1:ni,1:nj,1:nk,nql1)), &
        !             "bin_wl", minval(q3d(1:ni,1:nj,1:nk,nnc1:nnc1+num_bins-1)/scale_bin_wl), maxval(q3d(1:ni,1:nj,1:nk,nnc1:nnc1+num_bins-1)/scale_bin_wl), &
        !             "p", minval(p3d(1:ni,1:nj,1:nk)), maxval(p3d(1:ni,1:nj,1:nk)), "theta", minval(th3d(1:ni,1:nj,1:nk)), maxval(th3d(1:ni,1:nj,1:nk)), &
        !             "rho", minval(rho3d(1:ni,1:nj,1:nk)), maxval(rho3d(1:ni,1:nj,1:nk)), "supersat", minval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1)), maxval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1)),&
        !             "number concentraton", minval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2)), maxval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))     

        ! print *, "run dcpm ", dt, "s"
        ! do i3d=1,ni
        ! do j3d=1,nj
        ! do k3d=1,nk
        !     z = zh(i3d,j3d,k3d)
        !     t_dcpm = 0
        !     wv = q3d(i3d,j3d,k3d,nqv)
        !     bin_r_wet_temp = bin_r_wet
        !     bin_r_dry = dcpm3d(i3d,j3d,k3d,1:num_bins)
        !     bin_num_conc = q3d(i3d,j3d,k3d,nnc1:nnc2)
        !     temp = (th0(i3d,j3d,k3d)+th3d(i3d,j3d,k3d))*(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d))
        !     p = p3d(i3d,j3d,k3d)
        !     rho_d = calculate_rho_d (temp, p, wv)
        !     wl = sum((4d0*pi/3d0)*bin_num_conc*(bin_r_wet_temp**3d0 - bin_r_dry**3d0) *rho_w/rho_d) 
        !     q3d(i3d,j3d,k3d,nql1) = wl
        !     supersat = calculate_supersat (temp, p, wv) 
        ! end do
        ! end do
        ! end do
        ! print *,    "extreme wl after updating wl", minval(q3d(1:ni,1:nj,1:nk,nql1)), maxval(q3d(1:ni,1:nj,1:nk,nql1))
        
        if_did_cond_stat = .false.
        if_did_coll_stat = .false.
        t_step_cond_min_stat = t_step_cond_max
        t_step_coll_min_negc = dt
        num_grid_active_stat = 0
        num_rk_stat = 0
        num_dist_stat = 0

        !$omp parallel do default(private)  &
        !$omp shared(pi0,th0,zh,pp3d,th3d,p3d,temp3d, rho3d,q3d,dcpm3d,dt) &
        !omp shared(z_kernels, kernels_coal, weights_coal, bin_r_wet, bin_r_wet_grid) &
        !omp shared(if_did_cond_stat,if_did_coll_stat,t_step_cond_min_stat)


        do i3d=1,ni
        do j3d=1,nj
        do k3d=1,nk

            if( radopt.eq.1 )then
                effc (i3d, j3d, k3d) = 25 !initialize with minimum value 25 microns
                ! effc (i3d, j3d, k3d) = 0
            elseif( radopt.eq.2 )then
                effc (i3d, j3d, k3d) = 25 *1e-6
                ! effc (i3d, j3d, k3d) = 0 *1e-6
            end if
 
            z = zh(i3d,j3d,k3d)
            t_dcpm = 0
            wv = q3d(i3d,j3d,k3d,nqv)
            wi = q3d(i3d,j3d,k3d,nqs1)
            bin_r_wet_temp = bin_r_wet
            bin_wl = q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl
            wl = sum(bin_wl)
            q3d(i3d,j3d,k3d,nql1) = wl
            bin_num_conc = q3d(i3d,j3d,k3d,nnc1+num_bins:nnc2)
            bin_num_conc_preadv = q3d_dcpm(i3d,j3d,k3d,nnc1+num_bins:nnc2)
            temp = (th0(i3d,j3d,k3d)+th3d(i3d,j3d,k3d))*(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d))
            p = p3d(i3d,j3d,k3d)
            hdf = exp((-(zh(i3d,j3d,k3d)-zf(i3d,j3d,1))/Hs_aerosol))  
            bin_num_conc_base_curr = bin_num_conc_base*hdf*scale_num_conc

            ! if (myid.eq.240 .and. i3d.eq.1 .and. j3d.eq.2 .and. k3d.eq.1) then
            !     print *, "At sample point before dcpm calculation:", " wl:", wl, "wv:", wv, "temp:",temp, "p:", p
            ! end if

            ! if (temp < 220) then
            !     print *, "ERROR:temperature is too low:", temp
            !     print *, "myid", myid, i3d, j3d, k3d, p, wv, q3d(i3d,j3d,k3d,nqs1)
            !     call stopcm1
            ! end if
            if (ieee_is_nan(temp)) then
                print *, "temperature is nan before microphysics."
                print *, "i3d, j3d, k3d:", i3d, j3d, k3d
                print *, "th0",th0(i3d,j3d,k3d),"th3d",th3d(i3d,j3d,k3d),"pi0",pi0(i3d,j3d,k3d),"pp3d",pp3d(i3d,j3d,k3d)
            end if
             
            ! Begin ice process following Grabowski 1998
            if (temp<temp_water) then
                wtot = wi+wl+wv 
                if (wtot < -tolerance_wv) then
                    print *, "ERROR: wtot is negative before ice process", "wtot, wi, wl, wv", wtot, wi, wl, wv
                    print *, "myid, i3d, j3d, k3d", myid, i3d, j3d, k3d
                    do ind_bin = 1, num_bins
                        print *, ind_bin, bin_r_wet(ind_bin), bin_wl(ind_bin)
                    end do
                    call stopcm1
                end if 
                count_loop = 0
                temp_guess = temp
                do while (.true.)
                    count_loop = count_loop+1
                    if (count_loop > 100) then
                        print *, "ERROR: ice preocess iterates over 100 times."
                        call stopcm1
                    end if
                    wv_sat_i = rsif(real(p),real(temp_guess))
                    wv_sat_w = rslf(real(p),real(temp_guess))
                    ! if (wv_sat_w <= wv_sat_i) then
                    !     print *, "ERROR:wv_sat_w<=wv_sat_i", wv_sat_w, wv_sat_i
                    !     call stopcm1
                    ! end if
                    if (temp_guess >= temp_water) then
                        portion_ice = 0d0
                    else if (temp_guess <= temp_ice) then
                        portion_ice = 1d0
                    else
                        portion_ice = (temp_water-temp_guess)/(temp_water-temp_ice)
                    end if
                    wv_sat = wv_sat_i*portion_ice+wv_sat_w*(1-portion_ice)
                    if (wtot > wv_sat) then
                        wi_result = (wtot-wv_sat)*portion_ice
                        wl_result = (wtot-wv_sat)*(1-portion_ice)         
                    else
                        wi_result = 0
                        wl_result = 0
                    end if
                    wv_result = wtot-wi_result-wl_result
                    temp_result = temp + L_ice*(wi_result-wi)/C_p + L*(wl_result-wl)/C_p
                    if (abs(temp_guess-temp_result) < tolerance_temp) then
                        exit
                    else
                        temp_guess = 0.5d0*(temp_guess+temp_result)
                    end if
                end do
                
                ! Restart the iteration of as final step
                temp_result = temp
                wv_result = wv
                ! Lose water immediately if wl larger than equlibrium (freezing) this is only condition to update wl in ice process
                if (wl_result < wl-tolerance_wl .and. wl >tolerance_wl) then
                    do ind_bin = 1, num_bins
                        if (bin_wl(ind_bin) > tolerance_bin_wl) then
                            bin_num_conc(ind_bin) = bin_num_conc(ind_bin)*wl_result/wl
                            bin_wl(ind_bin) = bin_wl(ind_bin)*wl_result/wl
                        end if
                    end do
                    wl_result = sum(bin_wl)
                    wv_result = wv_result-(wl_result-wl)
                    temp_result = temp_result+L*(wl_result-wl)/C_p
                    wl = wl_result
                end if

                wv_result = wv_result-(wi_result-wi)
                temp_result = temp_result + L_ice*(wi_result-wi)/C_p 
                wi = wi_result
                wv = wv_result
                temp = temp_result
                ! if (wi > 1e-4) then
                !     print *, "example calculation:", count_loop, i3d,j3d,k3d,temp, L_ice*(wi-0)/C_p, wv, wv_sat, wv_sat_w, wv_sat_i, wi, wv
                !     call stopcm1
                ! end if
                if (wi<-tolerance_wi .or. wv<-tolerance_wv .or. wl<-tolerance_wl) then
                    print *, "ERROR: got negative values in ice processes"
                    print *, "wv, wl, wi", wv, wl, wi
                    print *, "temp, wv_sat_w, wv_sat_i, wv_sat, portion_ice, wv_result, wl_result, wi_result", temp, wv_sat_w, wv_sat_i, wv_sat, portion_ice, wv_result, wl_result, wi_result
                    call stopcm1
                end if
            else
                wv = wv+wi
                temp = temp-L_ice*wi/C_p
                wi = 0
            end if 
            temp_error = abs((wl - q3d(i3d,j3d,k3d,nql1))*L/Cp+(wi - q3d(i3d,j3d,k3d,nqs1))*L_ice/Cp+temp3d(i3d,j3d,k3d)-temp)
            if (temp_error > 0.01 .and. temp_error/0.01>100) then
                print *, "myid", myid, "myi", myi, "myj", myj
                print *, "ERROR: the energy is imbalanced after ice scheme of grid i j k at time(s) ", i3d, j3d, k3d
                print *, "z", z, "wv", wv, "wl", wl, "wi", wi, "temp", temp, "p", p, "rho_d", rho_d, "supersat", supersat, "wv_init", wv_init
                print *, "was", temp3d(i3d,j3d,k3d), q3d(i3d,j3d,k3d,nqv) , q3d(i3d,j3d,k3d,nql1), q3d(i3d,j3d,k3d,nqs1), p3d(i3d,j3d,k3d), rho3d(i3d,j3d,k3d), dcpm3d(i3d,j3d,k3d,num_bins+1)
                print *, "wv change:", wv-q3d(i3d,j3d,k3d,nqv), "wl change:", wl-q3d(i3d,j3d,k3d,nql1), "wi change:", wi-q3d(i3d,j3d,k3d,nqs1)
                print *, "expected temperature change:", (wl - q3d(i3d,j3d,k3d,nql1))*L/Cp+(wi - q3d(i3d,j3d,k3d,nqs1))*L_ice/Cp, "actual temperature change:", temp-temp3d(i3d,j3d,k3d)
                call stopcm1
            end if

            ! End ice process
           
            if (wi >= tolerance_wi) then
                effi_curr = 2*(1.2433*(temp3d(i3d,j3d,k3d)-273.15)+8.6629*log(wi*rho3d(i3d,j3d,k3d)*1e6)+132.23)
                ! if (effi_curr > 360d0) then
                !     print *, "ice effective radius too large:", temp,wi, wv, "was:", (th0(i3d,j3d,k3d)+th3d(i3d,j3d,k3d))*(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d)), q3d(i3d,j3d,k3d,nqs1), q3d(i3d,j3d,k3d,nqv)
                ! end if
            else
                effi_curr = 25
            end if

            if( radopt.eq.1 )then
                effi (i3d, j3d, k3d) = max(25d0,min(effi_curr,360d0))
            elseif( radopt.eq.2 )then
                effi (i3d, j3d, k3d) = effi_curr*1e-6
            end if

            rho_d = calculate_rho_d (temp, p, wv)
            supersat = calculate_supersat (temp, p, wv) 
  
            ! bin_r_dry = dcpm3d(i3d,j3d,k3d,1:num_bins)
            ! if (supersat > 0.3) then
            !     print *, "WARNING: supersaturation too high"
            !     print *, i3d, j3d, k3d, "wv:", wv, "wl:", wl, "wi:", wi, "temp:",temp, "p:", p, "rho_d:", rho_d, "supersaturation:",supersat
            !     print *, "wv_sat_i:", rsif(real(p),real(temp)), "wv_sat_w:", rslf(real(p),real(temp))
            !     ! print *, bin_wl
            !     ! print *, bin_num_conc
            !     ! stop
            ! end if

            if (k3d .eq. 1) then
                if (if_aerosol_replenish .eqv. .true.) then
                    call fill_bin_r_dry (bin_r_wet, bin_r_dry_bc, supersat, temp) ! supersat and temp not used anymore, naive filling
                    call fill_bin_num_conc (bin_r_dry_bc, bin_num_conc_bc)
                    do ind_bin = 1, num_bins
                        if (bin_wl(ind_bin) < tolerance_bin_wl) then
                            bin_num_conc(ind_bin) = bin_num_conc_bc(ind_bin)
                        ! else
                        !     bin_num_conc(ind_bin) = max(bin_num_conc(ind_bin), bin_num_conc_bc(ind_bin))
                        end if
                    end do
                end if
            end if


            call calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl)
            rho_d = calculate_rho_d (temp, p, wv)
            supersat = calculate_supersat (temp, p, wv)  
            bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet_temp**3d0 - bin_r_dry**3d0) *rho_w/rho_d
            wl = sum(bin_wl)      

            if (wl > 0.02) then
                print *, "wl too large before condensation", wl, "myid", myid, "i3d,j3d,k3d", i3d,j3d,k3d
                print *, "sum of bin_wl:", sum(bin_wl), "wv:", wv, "temp:", temp, "supersat:", supersat, "pressure:", p3d(i3d,j3d,k3d)
                do ind_bin =1, num_bins
                    print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                end do
                call stopcm1
            end if

            ! Correct dry radius (the dry radius cannot exceed the maximum previous radius of neighboring bins)
            if_r_dry_adjusted = .false.
            if (if_correct_aerosol) then
                if (wl>tolerance_wl .or. maxval(bin_num_conc(ind_bin_base_empty_start:num_bins))>tolerance_bin_num_conc) then
                    if ((myi.ne.1 .or. i3d.ne.1) .and. (myi.ne.nodex .or. i3d.ne.ni) .and.(myj.ne.1 .or. j3d.ne.1) .and. (myj.ne.nodey .or. j3d.ne.nj)) then
                        call calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl)
                        do ind_bin = 1, num_bins
                            r_dry_max_possible = max(ccn_radius_min,maxval(dcpm3d(i3d-1:i3d+1,j3d,k3d,ind_bin)))
                            r_dry_max_possible = max(r_dry_max_possible,maxval(dcpm3d(i3d,j3d-1:j3d+1,k3d,ind_bin)))
                            r_dry_max_possible = max(r_dry_max_possible,maxval(dcpm3d(i3d,j3d,max(1,k3d-1):min(nk,k3d+1),ind_bin)))
                            if (bin_r_dry(ind_bin) > r_dry_max_possible) then
                                if_r_dry_adjusted = .true.
                                if (abs(bin_r_wet_temp(ind_bin)-bin_r_dry(ind_bin))<tolerance_bin_r) then
                                    bin_r_wet_temp(ind_bin) = bin_r_wet_temp(ind_bin)-bin_r_dry(ind_bin)+r_dry_max_possible
                                else
                                    r_wet_cubic = (bin_r_wet_temp(ind_bin)-bin_r_dry(ind_bin))*(bin_r_wet_temp(ind_bin)**2d0+bin_r_wet_temp(ind_bin)*bin_r_dry(ind_bin)+bin_r_dry(ind_bin)**2d0)+r_dry_max_possible**3d0
                                    bin_r_wet_temp(ind_bin) = r_wet_cubic**(1d0/3d0)  
                                end if 
                                bin_r_wet_temp(ind_bin) = max(bin_r_wet_temp(ind_bin), bin_r_wet(1))         
                                bin_r_dry(ind_bin) = r_dry_max_possible
                            end if
                        end do
                    end if
                end if
            end if

            if (if_redistribute_aerosol) then
                if (if_correct_aerosol) then
                    print *, "ERROR: if_redistribute_aerosol and if_correct_aerosol cant be both .true."
                    call stopcm1
                end if
                !------------------------------------------------------------
                ! Redistribute Dry Aerosols
                !
                ! Purpose:
                !   Adjust dry aerosol number concentrations to compensate for
                !   numerical drift caused by droplet coalescence and evaporation,
                !   while conserving total dry-aerosol volume.   
                !
                ! Workflow:
                ! 1. Copy original aerosol arrays to working (_temp) arrays.
                ! 2. For effectively dry bins (wl < tolerance_bin_wl):
                !      - Remove residual liquid water (update wv and temp).
                !      - Compute dry-volume difference: r_dry_vol_diff = V_temp - V_base
                !      - Classify bins:
                !          • Excess if r_dry_vol_diff > 0
                !          • Deficit if r_dry_vol_diff < 0
                ! 3. Sum total excess and deficit volumes.
                ! 4. If both excess and deficit exist:
                !      - Scale deficit volumes proportionally to total available excess.
                !      - Update N_temp for each deficit bin, capped at N_base.
                ! 5. Apply _temp arrays back to main arrays.
                ! 6. Optional diagnostics if test flags are enabled.
                !
                ! Notes:
                !   - Redistribution is volume-conservative: added volume per deficit
                !     bin is proportional to its deficit.
                !   - Excess is clipped to base (initial number concentration) to avoid creating extra aerosols.
                !   - This ensures dry-aerosol mass conservation but is **not true damping**.
                !------------------------------------------------------------

                if (wl < tolerance_wl) then
                    bin_wl = 0d0
                    bin_num_conc = bin_num_conc_base_curr
                    call fill_bin_r_dry (bin_r_wet, bin_r_dry, supersat, temp) ! supersat and temp not used anymore, naive filling
                    temp = temp-L*(wl)/C_p
                    wv = wv + wl
                    wl = 0
                else
                    bin_num_conc_temp = bin_num_conc
                    bin_r_dry_temp = bin_r_dry
                    bin_r_dry_vol_defecit = 0
                    bin_r_dry_vol_excess = 0
                    scale_r_dry_vol_tot = 0
                    if_r_dry_vol_defecit = .false.
                    if_r_dry_vol_excess = .false.
                    do ind_bin = 1, num_bins
                        if  (bin_wl(ind_bin) < tolerance_bin_wl) then
                            wv = wv+(bin_wl(ind_bin))
                            wl = wl-(bin_wl(ind_bin))
                            temp = temp-L*(bin_wl(ind_bin))/C_p
                            bin_wl(ind_bin) = 0
                            if (bin_r_wet(ind_bin) .gt. ccn_radius_max) then
                                r_dry_vol_diff = (4.0*pi/3.0) * bin_r_dry_temp(ind_bin)**3 * bin_num_conc_temp(ind_bin)
                            else
                                r_dry_vol_diff = (4.0*pi/3.0) * bin_r_dry_temp(ind_bin)**3 * bin_num_conc_temp(ind_bin) - &
                                                (4.0*pi/3.0) * bin_r_wet(ind_bin)**3 * bin_num_conc_base_curr(ind_bin)
                            end if
                            if (bin_num_conc_temp(ind_bin) .gt. bin_num_conc_base_curr(ind_bin)) then
                                bin_num_conc_temp(ind_bin) = bin_num_conc_base_curr(ind_bin)
                            end if
 
                            if (bin_r_wet(ind_bin) > ccn_radius_max) then
                                bin_r_dry_temp(ind_bin) = 0
                            else
                                bin_r_dry_temp(ind_bin) = bin_r_wet(ind_bin)
                            end if

                            if (r_dry_vol_diff .gt. 0) then
                                if_r_dry_vol_excess = .true.
                                bin_r_dry_vol_excess(ind_bin) = r_dry_vol_diff 
                            else if (r_dry_vol_diff .lt. 0) then
                                if_r_dry_vol_defecit = .true.
                                bin_r_dry_vol_defecit(ind_bin) = -r_dry_vol_diff 
                            end if
                        end if
                    end do
                    r_dry_vol_excess_tot = sum(bin_r_dry_vol_excess)
                    r_dry_vol_defecit_tot = sum(bin_r_dry_vol_defecit)

                    ! Compute tolerance: 0.1% of the smallest non-zero bin volume
                    r_dry_vol_tolerance = huge(1.0)  ! start with a very large number
                    do ind_bin = 1, num_bins
                        if (bin_num_conc_base_curr(ind_bin) > tolerance_bin_num_conc) then
                            r_dry_vol_tolerance = min(r_dry_vol_tolerance,(4.0*pi/3.0) * bin_r_wet(ind_bin)**3 * bin_num_conc_base_curr(ind_bin))
                        end if
                    end do
                    r_dry_vol_tolerance = r_dry_vol_tolerance * 0.001

                    if ((r_dry_vol_excess_tot .gt. r_dry_vol_tolerance) .and. (r_dry_vol_defecit_tot .gt. r_dry_vol_tolerance)) then
                        scale_r_dry_vol_tot = r_dry_vol_excess_tot/r_dry_vol_defecit_tot
                        do ind_bin = 1, num_bins
                            if ((bin_r_dry_temp(ind_bin) .gt. tolerance_bin_r) .and. (bin_wl(ind_bin) .lt. tolerance_bin_wl) .and. &
                            bin_num_conc_temp(ind_bin) .lt. (bin_num_conc_base_curr(ind_bin)-tolerance_bin_num_conc)) then
                                num_conc_increase = bin_r_dry_vol_defecit(ind_bin)*scale_r_dry_vol_tot/(bin_r_dry_temp(ind_bin)**3*4*pi/3)
                                bin_num_conc_temp(ind_bin) = min(bin_num_conc_base_curr(ind_bin),bin_num_conc_temp(ind_bin)+num_conc_increase)
                            end if
                        end do
                    end if
                    

                    if (if_test .and. (maxval(abs(bin_num_conc-bin_num_conc_temp)) .gt. tolerance_bin_num_conc) .and. (myid.eq.7) .and. .false.) then
                        print *, "output of test case: myid, myi, myj, i3d, j3d, k3d", myid, myi, myj, i3d, j3d, k3d
                        print *, "r_dry_vol_excess_tot",r_dry_vol_excess_tot,"r_dry_vol_defecit_tot",r_dry_vol_defecit_tot,"scale_r_dry_vol_tot",scale_r_dry_vol_tot, "r_dry_vol_tolerance",r_dry_vol_tolerance
                        print *, "supersat", supersat, "wl", wl
                        print *, "wv", wv, "temp", temp
                        print *, "ind_bin, bin_r_dry(ind_bin), bin_r_dry_temp(ind_bin),bin_r_wet(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), bin_num_conc_temp(ind_bin), bin_num_conc_base_curr(ind_bin), bin_r_dry_vol_excess(ind_bin), bin_r_dry_vol_defecit(ind_bin)"
                        do ind_bin = 1, num_bins
                            print *, ind_bin, bin_r_dry(ind_bin), bin_r_dry_temp(ind_bin), bin_r_wet(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), bin_num_conc_temp(ind_bin), bin_num_conc_base_curr(ind_bin), bin_r_dry_vol_excess(ind_bin), bin_r_dry_vol_defecit(ind_bin)
                        end do
                        call stopcm1
                    end if
                    bin_num_conc = bin_num_conc_temp
                    bin_r_dry = bin_r_dry_temp
                end if
            end if

            supersat = calculate_supersat (temp, p, wv)  
            q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1) = bin_wl*scale_bin_wl
            q3d(i3d,j3d,k3d,nnc1+num_bins:nnc2) = bin_num_conc
            dcpm3d(i3d,j3d,k3d,1:num_bins) = bin_r_dry
            q3d(i3d,j3d,k3d,nql1) = wl
            q3d(i3d,j3d,k3d,nqs1) = wi
            q3d(i3d,j3d,k3d,nqv) = wv
            temp3d(i3d,j3d,k3d) = temp
            dcpm3d(i3d,j3d,k3d,num_bins+1) = supersat
            adddiag3d(i3d,j3d,k3d,1) = supersat ! supersaturation before microphysics
            th3d(i3d,j3d,k3d) = temp/(p/p00)**(1.0/cpdrd)-th0(i3d,j3d,k3d)
            rho3d(i3d,j3d,k3d) = rho_d

            if (supersat < 0 .and. wl < tolerance_wl) then
                adddiag3d(i3d,j3d,k3d,2) = supersat ! supersaturation after microphysics (no change here since skip)
                if (if_r_dry_adjusted) then
                    call update_bins_static(bin_r_dry, bin_r_wet_temp, bin_num_conc, "after aerosol correction")         
                end if       
                cycle
            end if
            
            if_grid_counted = .false.
            wv_init = wv
            wl_init = wl
            temp_init = temp
            t_step_dcpm = dt

            do while (t_dcpm - dt < -tolerance_t)

                t_cond = 0
                if_did_cond_curr = .false.
                do while(t_cond-t_step_dcpm < -tolerance_t)
                ! Can't break out loop here: need to run collision coalescence & breakup
                    t_step_cond = min(t_step_cond_max,t_step_dcpm-t_cond)
                    t_step_cond_max_curr = t_step_cond
                    rho_d = calculate_rho_d (temp, p, wv)               
                    wv_sat = calculate_wv (temp, 0d0, p)
                    supersat = calculate_supersat (temp, p, wv)   

                    num_cond_bins = 0
                    do ind_bin = 1, num_bins
                        if (bin_num_conc(ind_bin) > tolerance_bin_num_conc) then
                            num_cond_bins = max(num_cond_bins, ind_bin)
                        end if
                    end do

                    ! update if bin is activated
                    bin_if_act = .false.
                    if_do_cond = .false.

                    do ind_bin = 1, num_cond_bins
                        if (abs(bin_num_conc(ind_bin)) < tolerance_bin_num_conc) then
                            bin_if_act(ind_bin) = .false.
                            cycle
                        end if  

                        supersat_crit = calculate_supersat_crit (temp, bin_r_dry(ind_bin), aer_kappa)

                        r_crit = calculate_r_crit (temp, bin_r_dry(ind_bin), aer_kappa)

                        cond_wl         = (bin_wl(ind_bin) >= tolerance_bin_wl)
                        cond_dry        = (supersat < calculate_equilibrium_supersat (bin_r_wet_temp(ind_bin), bin_r_dry(ind_bin), temp, aer_kappa))
                        cond_radius     = (bin_r_wet_temp(ind_bin) > r_crit*(1+tolerance_r_ratio))
                        cond_supersat   = (supersat_crit <= supersat)

                        if_evaporation = (cond_wl .and. cond_dry)
                        if_condensation = (cond_radius .or. cond_supersat) .and. (supersat>0)

                        if (if_evaporation .or. if_condensation) then
                            bin_if_act(ind_bin) = .true.
                            if_do_cond = .true.
                        else
                            bin_if_act(ind_bin) = .false.
                        end if
                   
                    end do        


                    ! if (abs(wv - wv_init) > 1e-4) then
                    !     print *, "wv not equal to wv_init", wv, wv_init
                    !     print *, "t_dcpm", t_dcpm, "dt", dt
                    !     print *, "wl, wl_init", wl, wl_init
                    !     call stopcm1
                    ! end if
                    ! if (abs(wl - wl_init) > 1e-4) then
                    !     print *, "wl not equal to wl_init", wl, wl_init
                    !     call stopcm1
                    ! end if
                    ! condensation
                    if(if_do_cond) then
                        if_did_cond_curr = .true.

                        if (if_grid_counted .eqv. .false.) then
                            !$omp atomic update
                            if_did_cond_stat = .true.
                            !$omp atomic update
                            num_grid_active_stat = num_grid_active_stat + 1
                            if_grid_counted = .true.
                        end if

                        num_iter = 0
                        do while(.true.)
                            num_iter = num_iter + 1
                            ! print *, "iteration #", num_iter, "t_step_cond:", t_step_cond
                            bin_shift_curr = bin_r_wet_temp-bin_r_wet_temp
                            call calculate_bin_deriv (k1, bin_r_wet_temp, bin_shift_curr, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins)
                            bin_shift_curr = k1*0.5*t_step_cond
                            call calculate_bin_deriv (k2, bin_r_wet_temp, bin_shift_curr, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins)
                            bin_shift_curr = k2*0.5*t_step_cond
                            call calculate_bin_deriv (k3, bin_r_wet_temp, bin_shift_curr, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins)
                            bin_shift_curr = k3*1.0*t_step_cond
                            call calculate_bin_deriv (k4, bin_r_wet_temp, bin_shift_curr, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins)
 
                            if (if_square_radius_scheme) then
                                bin_r_wet_predict = bin_r_wet_temp**2d0+t_step_cond/6.0*(k1+2.0*k2+2.0*k3+k4) 
                                do ind_bin = 1, num_bins
                                    if (bin_r_wet_predict(ind_bin) <= bin_r_dry(ind_bin)**2d0) then
                                        bin_r_wet_predict(ind_bin) = bin_r_dry(ind_bin)
                                    else 
                                        bin_r_wet_predict(ind_bin) = sqrt(bin_r_wet_predict(ind_bin))
                                    end if
                                end do                       
                            else
                                bin_r_wet_predict = bin_r_wet_temp+t_step_cond/6d0*(k1+2.0*k2+2.0*k3+k4) 
                                do ind_bin = 1, num_bins
                                    if (bin_r_wet_predict(ind_bin) <= bin_r_dry(ind_bin)) then
                                        bin_r_wet_predict(ind_bin) = bin_r_dry(ind_bin)
                                    end if
                                end do
                            end if

                            ! wl_predict = wl + sum((4d0*pi/3d0)*bin_num_conc*(bin_r_wet_predict- bin_r_wet_temp)*(bin_r_wet_predict**2d0+bin_r_wet_predict*bin_r_wet_temp+bin_r_wet_temp**2d0) *rho_w/rho_d) 
                            wl_predict = sum((4d0*pi/3d0)*bin_num_conc*(bin_r_wet_predict**3d0-bin_r_dry**3d0) *rho_w/rho_d) 
                            if (abs(wl_predict-wl) < tolerance_wl*0.01) then
                                exit
                            end if
                            wv_predict = wv_init - (wl_predict-wl_init)
                            temp_predict = temp_init + L*(wl_predict-wl_init)/C_p
                            if (ieee_is_nan(temp_predict)) then
                                print *, "got nan temp_predict"
                                print *, "temp_init", temp_init, "L", L, "wl_predict", wl_predict, "wl_init", wl_init, "C_p", C_p
                                do ind_bin = 1, num_bins
                                    print *, bin_r_wet_temp(ind_bin), k1(ind_bin), k2(ind_bin), k3(ind_bin), k4(ind_bin), bin_r_wet_predict(ind_bin), bin_r_dry(ind_bin), bin_num_conc(ind_bin)
                                end do
                                call stopcm1
                            end if
                            wv_sat_predict = calculate_wv (temp_predict, 0d0, p)
 
                            if (abs(sum((bin_r_wet_predict**3d0-bin_r_dry**3d0)*4d0/3d0*pi*bin_num_conc*rho_w/rho_d) - wl_predict) > max(tolerance_wl,0.001 *wl_predict) ) then
                                print *, "WARNING mass imbalance occurs at radius prediction"
                                print *, sum((bin_r_wet_predict**3d0-bin_r_dry**3d0)*4d0/3d0*pi*bin_num_conc*rho_w/rho_d), wl_predict, sum((bin_r_wet_predict**3d0-bin_r_dry**3d0)*4d0/3d0*pi*bin_num_conc*rho_w/rho_d) - wl_predict
                            end if
                            ! print *, "wv:", wv
                            t_step_cond_ratio_radius = 1.0
                            t_step_cond_ratio_temperature = 1.0
                            t_step_cond_ratio_moisture = 1.0
                            if (t_step_cond > t_step_cond_min+tolerance_t .or. num_iter < max_iter_cond) then
                                call determine_t_step_cond_ratio (t_step_cond_ratio_radius, t_step_cond_ratio_temperature, t_step_cond_ratio_moisture, ind_bin_radius_control, bin_r_wet_predict, bin_r_wet_temp, bin_r_dry, bin_num_conc, wv, wv_predict, wv_sat, wv_sat_predict, wl, wl_predict, temp, temp_predict, num_cond_bins) 
                                if (t_step_cond*t_step_cond_ratio_radius < t_step_cond_min+tolerance_t) then
                                    if (if_warning_timestep) then
                                        print *, "WARNING condensation time step limited by radius is smaller than limit:", t_step_cond*t_step_cond_ratio_radius
                                        print *, "Use minimum condensation time step ", t_step_cond_min, " instead."
                                        print *, t_step_cond_ratio_radius, t_step_cond_ratio_temperature, t_step_cond_ratio_moisture
                                        print *, "ind_bin_radius_control:", ind_bin_radius_control
                                        print *, "wet radius of bin:", bin_r_wet_temp(ind_bin_radius_control), "deriv of bin:", (1d0/6d0*(k1(ind_bin_radius_control)+2.0*k2(ind_bin_radius_control)+2.0*k3(ind_bin_radius_control)+k4(ind_bin_radius_control)))
                                        print *, "predict wet radius of bin:", bin_r_wet_predict(ind_bin_radius_control), "neighboring grids:", bin_r_wet_grid(ind_bin_radius_control), bin_r_wet_grid(ind_bin_radius_control+1)
                                        print *, "if_square_radius_scheme:", if_square_radius_scheme
                                    end if
                                    t_step_cond_ratio_radius = t_step_cond_min/t_step_cond
                                end if
                            end if
                            if (abs(t_step_cond_ratio_radius-1)<1d-6 .and. abs(t_step_cond_ratio_temperature-1)<1d-6 .and. abs(t_step_cond_ratio_moisture-1)<1d-6) then
                                ! if (abs(wv-wv_predict) < tolerance_wv) then
                                !     print *, "condensation did nothing", wv, wv_predict
                                !     do ind_bin=1, num_bins
                                !         print *, bin_r_wet_temp(ind_bin), bin_r_wet_predict(ind_bin), k1(ind_bin), k2(ind_bin), k3(ind_bin), k4(ind_bin), t_step_cond
                                !     end do
                                ! end if
                                t_step_cond = t_step_cond_max_curr
                                num_rk_stat = num_rk_stat + 1   
                                bin_r_wet_temp = bin_r_wet_predict
                                wl = wl_predict
                                wv = wv_predict
                                wv_sat = wv_sat_predict
                                temp = temp_predict
                                if (wv < -tolerance_wv) then
                                    print *, "natural growth ERROR: negative water vapor mixing ratio!", "myid", myid, "index", i3d, j3d, k3d, "wv:", wv, "wl:", wl, "wv_init:", wv_init, "wl_init:", wl_init
                                    print *, "wv:", wv, "temp:", temp, "supersat:", supersat
                                    do ind_bin =1, num_bins
                                        print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                                    end do
                                    call stopcm1
                                end if
                                exit
                            end if

                            t_step_cond_radius = t_step_cond*t_step_cond_ratio_radius
                            t_step_cond_temperature = t_step_cond*t_step_cond_ratio_temperature
                            t_step_cond_moisture = t_step_cond*t_step_cond_ratio_moisture

                            ! if (t_step_cond_temperature < t_step_cond_min*0.5 .and. abs(t_step_cond_ratio_radius-1)>1d-6 .and. myid.eq.25) then
                            !     print *, "ERROR condensation time step limited by temperature is too small:", t_step_cond_temperature
                            !     print *, "myid:",myid, "wv:", wv, "wv_predict:", wv_predict, "wl:", wl, "wl_predict:", wl_predict, "temp:", temp, "temp_predict", temp_predict, "supersat:", supersat, "wl:", wl
                            !     call calculate_bin_deriv (bin_deriv, bin_r_wet_temp, bin_r_wet_temp-bin_r_wet_temp, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins) 
                            !     do ind_bin =1, num_bins
                            !         print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d, bin_deriv(ind_bin)

                            !     end do
                            !     call stopcm1
                            ! end if
                            if (t_step_cond_moisture < t_step_cond_min) then
                                t_step_cond = min(t_step_cond_radius, t_step_cond_temperature)
                                ! if (wv <= wv_sat) then  
                                !     print *, "ERROR: perform direct volume distribution at subsaturation"
                                !     stop
                                ! end if
                                !$omp atomic update
                                dwl = calculate_equilibrium_dwl(wv, temp, p)
                                if (ieee_is_nan(dwl)) then
                                    print *, "dwl is nan", "wv", wv,"temp", temp,"p", p
                                    call stopcm1
                                end if                                
                                if (abs(dwl) .le. tolerance_wl*10) then
                                    exit
                                end if                                
                                call calculate_bin_deriv (bin_deriv, bin_r_wet_temp, bin_r_wet_temp-bin_r_wet_temp, bin_r_dry, bin_if_act, bin_num_conc, wv, temp, p, num_cond_bins)         
                                if (if_square_radius_scheme) then
                                    bin_weight = bin_deriv*1.5d0*bin_r_wet_temp*bin_num_conc
                                else
                                    bin_weight = bin_deriv*3.0d0*bin_r_wet_temp**2*bin_num_conc
                                end if

                                do ind_bin = 1, num_bins
                                    if (ieee_is_nan(bin_weight(ind_bin)/sum(bin_weight))) then
                                        print *, "bin_weight/sum(bin_weight)", bin_weight/sum(bin_weight)
                                        print *, "bin_weight", bin_weight
                                        print *, "bin_deriv", bin_deriv
                                        call stopcm1
                                    end if
                                end do
                                
                                bin_weight = bin_weight/sum(bin_weight)
                                do ind_bin = 1, num_bins
                                    if ((bin_weight(ind_bin) .le. 1e-3) .or. (abs(dwl*bin_weight(ind_bin)) .lt. tolerance_bin_wl)) then
                                        bin_weight(ind_bin) = 0
                                    end if
                                end do
                                if (sum(bin_weight) < 1e-2) then
                                    exit
                                end if                                
                                bin_weight = bin_weight/sum(bin_weight) 

                                num_dist_stat = num_dist_stat + 1 
                                do ind_bin = 1, num_bins
                                    if (bin_weight(ind_bin) .ge. 1e-4) then
                                        bin_r_wet_temp(ind_bin) = (bin_r_wet_temp(ind_bin)**3d00 + dwl*rho_d/rho_w*bin_weight(ind_bin)*3d0/4d0/pi/bin_num_conc(ind_bin))**(1d0/3d0)
                                    end if
                                    if (ieee_is_nan(bin_r_wet_temp(ind_bin))) then
                                        print *, "ind_bin", ind_bin, "bin_r_wet_temp(ind_bin)", bin_r_wet_temp(ind_bin), bin_r_wet_temp(ind_bin), dwl, rho_d, bin_num_conc(ind_bin)
                                        call stopcm1
                                    end if
                                end do
                                wl_predict = sum((4d0*pi/3d0)*bin_num_conc*(bin_r_wet_temp**3d0 - bin_r_dry**3d0) *rho_w/rho_d) 
                                if (abs(wl_predict-wl-dwl) .ge. tolerance_wl*10) then
                                    print *, "mass distribution WARNING: mass not conserved!", "index", i3d, j3d, k3d, "wl_predict:", wl_predict, "wl:", wl, "dwl:", dwl, "abs(wl_predict-wl-dwl)", tolerance_wl
                                    do ind_bin = 1, num_bins
                                            print *, "ind_bin:", ind_bin, "bin weight:", bin_weight(ind_bin), "bin deriv:", bin_deriv(ind_bin)
                                    end do
                                    ! call stopcm1
                                end if
                                wl = wl_predict
                                wv = wv_init - (wl-wl_init)
                                if (wv < -tolerance_wv) then
                                    print *, "mass distribution ERROR: negative water vapor mixing ratio!", "index", i3d, j3d, k3d, "wv:", wv, "wl:", wl, "wv_init:", wv_init, "wl_init:", wl_init
                                    print *, "wv:", wv, "temp:", temp, "supersat:", supersat
                                    do ind_bin =1, num_bins
                                        print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                                    end do
                                    call stopcm1
                                end if
                                temp = temp + L*dwl/C_p

                                wv_sat = calculate_wv (temp, 0d0, p)     
                                ! if (abs(wv-wv_sat) < tolerance_wv) then
                                !     print *, "ERROR: the supersaturation is not 0 after distributing water vapor."
                                !     print *, "after distribution, wv:", wv, "wv_sat", wv_sat, "wl:", wl, "temp", temp, "dwl:", dwl
                                ! end if                 
                                exit
                            else
                                ! prepare for next iteration
                                t_step_cond = min(min(min(t_step_cond_radius, t_step_cond_temperature), t_step_cond_moisture),t_step_cond)
                                t_step_cond = t_step_cond*0.9 ! Force the time step decrease and incrase success probability for next iteration
                                cycle
                            end if      
                                       
                        end do

 
                        if (wl > 0.02) then
                            print *, "wl too large", wl, "myid", myid, "i3d,j3d,k3d", i3d,j3d,k3d
                            print *,"wv:", wv, "temp:", temp, "supersat:", supersat, "pressure:", p3d(i3d,j3d,k3d)
                            do ind_bin =1, num_bins
                                print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                            end do
                            call stopcm1
                        end if

                        do ind_bin=1, num_bins
                            if (bin_r_wet_temp(ind_bin)*(1+tolerance_r_ratio) < bin_r_dry(ind_bin)) then
                                print *, "ERROR: resulting wet radius smaller than dry radius", "ind_bin", ind_bin,"bin_r_wet_temp(ind_bin)",bin_r_wet_temp(ind_bin),"bin_r_dry(ind_bin)",bin_r_dry(ind_bin)
                                call stopcm1 
                            end if
                            bin_r_wet_temp(ind_bin) = max(bin_r_wet_temp(ind_bin), bin_r_dry(ind_bin))
                        end do                      

                        if (wv < -tolerance_wv) then
                            print *, "line 683 ERROR: negative water vapor mixing ratio!", "index", i3d, j3d, k3d, "wv:", wv, "wl:", wl, "wv_init:", wv_init, "wl_init:", wl_init
                            print *, "wv:", wv, "temp:", temp, "supersat:", supersat
                            do ind_bin =1, num_bins
                                print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                            end do
                            call stopcm1
                        end if

                        wv = max(wv,0.0)
                        wl = max(wl,0.0)
                    end if
 
                    !$omp atomic update
                    t_step_cond_min_stat = min(t_step_cond_min_stat, t_step_cond)

                    if (minval(bin_num_conc) < -tolerance_bin_num_conc) then
                        print *, "ERROR: negative number concentration after condensation:", minval(bin_num_conc)
                        call stopcm1
                    end if
                
                    t_cond = t_cond + t_step_cond

                end do

                if ( ((supersat > 0.1 .and. wl .le. 1e-6) .or. (supersat < -0.1 .and. wl .ge. 1e-6)) .and. (if_test .and. myid.eq.16 .and. .false.)) then
                    print *, "WARNING: wet rate or dry rate so slow at myid, i3d, j3d, k3d", myid, i3d, j3d, k3d
                    print *, "wv:", wv, "temp:", temp, "supersat:", supersat, "pressure:", p, "rho_d", rho_d
                    print *, "change of wl by condensation:", wl-wl_init, "t_step_cond", t_step_cond
                    do ind_bin =1, num_bins
                        print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), k1(ind_bin), k2(ind_bin), k3(ind_bin), k4(ind_bin), bin_if_act(ind_bin)
                    end do
                    call stopcm1
                end if

                if (minval(bin_num_conc) < -tolerance_bin_num_conc) then
                    print *, "ERROR: negative number concentration just before assignment."
                    call stopcm1
                end if

                if (if_did_cond_curr .eqv. .true.) then
                    bin_wl = (bin_r_wet_temp**3d0-bin_r_dry**3d0)*4d0/3d0*pi*bin_num_conc*rho_w/rho_d
                    do ind_bin = 1, num_bins
                        if (bin_wl(ind_bin) < tolerance_bin_wl) then
                            wv = wv + bin_wl(ind_bin)
                            temp = temp-L*(bin_wl(ind_bin))/C_p
                            wl = wl-bin_wl(ind_bin)
                            bin_wl(ind_bin) = 0
                            bin_r_wet_temp(ind_bin) = bin_r_dry(ind_bin)
                        end if
                    end do
                    call update_bins_static(bin_r_dry, bin_r_wet_temp, bin_num_conc, "after condensation")                
                end if

                bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet**3d0 - bin_r_dry**3d0) *rho_w/rho_d 
                call calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl) ! This will also update bin_wl if necesary
                
                ! Collision
                ! if (dt-t_dcpm < tolerance_t .or. t_cond > t_step_coll_min) then
                t_coll = 0
                if (wl > 0.02) then
                    print *, "wl too large after condensation but before collision", wl, "myid", myid, "i3d,j3d,k3d", i3d,j3d,k3d
                    print *, "wv:", wv, "temp:", temp, "supersat:", supersat, "pressure:", p3d(i3d,j3d,k3d)
                    do ind_bin =1, num_bins
                        print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                    end do
                    call stopcm1
                end if
                if (if_coll .and. (wl .ge. threshold_wl_coll)) then
                    bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet_temp**3d0 - bin_r_dry**3d0) *rho_w/rho_d 
                    wl = sum(bin_wl) 
                    bin_ind_coal(:) = 0
                    bin_ind_brkp(:) = 0
                    num_coal_bins_curr = 0
                    num_brkp_bins_curr = 0
                    do ind_bin = 1, num_bins
                        if ((bin_num_conc(ind_bin) >= threshold_bin_num_conc_coll) .and. &
                            (bin_wl(ind_bin) >= threshold_bin_wl_coll) .and. &
                            (bin_r_wet(ind_bin) >= bin_r_dry(ind_bin)*ratio_wet_coll_min)) then
                            if (bin_r_wet(ind_bin) >= drop_coal_small_min) then
                                num_coal_bins_curr  = num_coal_bins_curr +1
                                bin_ind_coal(num_coal_bins_curr ) = ind_bin
                            end if
                            if (bin_r_wet(ind_bin) >= drop_brkp_parent_min) then
                                num_brkp_bins_curr = num_brkp_bins_curr+1
                                bin_ind_brkp(num_brkp_bins_curr) = ind_bin
                            end if
                        end if
                    end do
                    if (num_coal_bins_curr >= 2) then
                        if ((bin_r_wet(bin_ind_coal(num_coal_bins_curr)) >= drop_coal_large_min)) then
                            if_coal = .true.
                        end if
                    else
                        if_coal = .false.
                    end if
                    if (num_brkp_bins_curr >=2) then
                        if_brkp = .true.
                    else
                        if_brkp = .false.
                    end if
                    
                    if (if_coal .or. if_brkp) then
                        if_did_coll_stat = .true.
                        
                        kernels_coal_curr = 0
                        kernels_brkp_curr = 0
                        weights_coal_curr = 0
                        weights_brkp_curr = 0
                        if (z < z_kernels(1) .or. z>z_kernels(nk_kernel)) then
                            print *, "ERROR: the height is not within kernel height range."
                            call stopcm1             
                        end if           
                        ! prepare coalescence kernel
                        if (if_coal) then
                            do ind_zk = 1, nk_kernel-1
                                if (z >= z_kernels(ind_zk) .and. z <= z_kernels(ind_zk+1)) then
                                    kernels_coal_curr = kernels_coal(ind_zk,:,:) + (kernels_coal(ind_zk+1,:,:)-kernels_coal(ind_zk,:,:))*(z-z_kernels(ind_zk))/(z_kernels(ind_zk+1)-z_kernels(ind_zk))
                                    exit
                                end if
                            end do
                            weights_coal_curr = weights_coal(:,:,:)
                            if (maxval(abs(kernels_coal_curr)) <1d-10) then
                                print *, "ERROR: fail to interpolate coalescence kernel"
                                call stopcm1
                            end if  
                        end if
                        ! prepare breakup kernel
                        if (if_brkp) then
                            do ind_zk = 1, nk_kernel-1
                                if (z >= z_kernels(ind_zk) .and. z <= z_kernels(ind_zk+1)) then
                                    kernels_brkp_curr = kernels_brkp(ind_zk,:,:) + (kernels_brkp(ind_zk+1,:,:)-kernels_brkp(ind_zk,:,:))*(z-z_kernels(ind_zk))/(z_kernels(ind_zk+1)-z_kernels(ind_zk))
                                    exit
                                end if
                            end do
                            if (nk_brkp_kernel .eq. 1) then
                                weights_brkp_curr = weights_brkp(1,:,:,:)
                            else
                                if (z < z_brkp_kernels(1) .or. z>z_brkp_kernels(nk_brkp_kernel)) then
                                    print *, "ERROR: the height is not within kernel height range."
                                    call stopcm1
                                else 
                                    do ind_zk = 1, nk_brkp_kernel-1
                                        if (z >= z_brkp_kernels(ind_zk) .and. z <= z_brkp_kernels(ind_zk+1)) then
                                            weights_brkp_curr = weights_brkp(ind_zk,:,:,:) + (weights_brkp(ind_zk+1,:,:,:)-weights_brkp(ind_zk,:,:,:))*(z-z_brkp_kernels(ind_zk))/(z_brkp_kernels(ind_zk+1)-z_brkp_kernels(ind_zk))
                                            exit
                                        end if
                                    end do
                                end if
                            end if  
                            if (maxval(abs(kernels_brkp_curr)) <1d-10) then
                                print *, "ERROR: fail to interpolate breakup kernel"
                                call stopcm1
                            end if                                                        
                        end if
      
                        do while (t_coll-t_step_dcpm < -tolerance_t) 
                            t_step_coll = min(t_step_coll_max,t_step_dcpm-t_coll)
                            
                            do ind_bin = 1, num_bins
                                if (bin_num_conc(ind_bin) < tolerance_bin_num_conc) then
                                    bin_num_conc(ind_bin) = 0.0
                                    bin_wl(ind_bin) = 0
                                    if (bin_r_wet_temp(ind_bin) <= ccn_radius_max) then
                                        bin_r_dry(ind_bin) = bin_r_wet_temp(ind_bin)
                                    else
                                        bin_r_dry(ind_bin) = 0
                                    end if
                                end if 
                            end do

                            bin_dnum_conc_dt = 0
                            bin_dwl_dt = 0
                            if (if_coal) then
                                do ind_ind_bin_i = 1, num_coal_bins_curr
                                    do ind_ind_bin_j = ind_ind_bin_i, num_coal_bins_curr
                                        ind_bin_i = bin_ind_coal(ind_ind_bin_i)
                                        ind_bin_j = bin_ind_coal(ind_ind_bin_j)
                                        if (bin_r_wet(ind_bin_j) < drop_coal_large_min) then
                                            cycle
                                        end if
                                        ! Handles number concentration (conc) and liquid water content (wl) during coalescence between ind_bin_i and ind_bin_j
                                        num_conc_coal = kernels_coal_curr (ind_bin_i, ind_bin_j) * bin_num_conc(ind_bin_i) * bin_num_conc(ind_bin_j)
                                        wl_coal_i = kernels_coal_curr (ind_bin_i, ind_bin_j) * bin_wl(ind_bin_i) * bin_num_conc(ind_bin_j)
                                        wl_coal_j = kernels_coal_curr (ind_bin_i, ind_bin_j) * bin_wl(ind_bin_j) * bin_num_conc(ind_bin_i)
                                        bin_dnum_conc_dt(ind_bin_i) = bin_dnum_conc_dt(ind_bin_i) - num_conc_coal
                                        bin_dwl_dt(ind_bin_i) = bin_dwl_dt(ind_bin_i) - wl_coal_i
                                        if (ind_bin_j .ne. ind_bin_i) then
                                            bin_dnum_conc_dt(ind_bin_j) = bin_dnum_conc_dt(ind_bin_j) - num_conc_coal
                                            bin_dwl_dt(ind_bin_j) = bin_dwl_dt(ind_bin_j) - wl_coal_j
                                        end if
                                        
                                        ! do ind_bin = max(ind_bin_i,ind_bin_j), num_bins
                                        !     bin_dnum_conc_dt(ind_bin) = bin_dnum_conc_dt(ind_bin) + &
                                        !     (1-0.5*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                        !     weights_coal_curr(ind_bin, ind_bin_i, ind_bin_j) * &
                                        !     num_conc_coal

                                        !     bin_dwl_dt(ind_bin) = bin_dwl_dt(ind_bin) + &
                                        !     (1-0.5*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                        !     weights_coal_curr(ind_bin, ind_bin_i, ind_bin_j) * (bin_r_wet(ind_bin)**3)/(bin_r_wet(ind_bin_i)**3+bin_r_wet(ind_bin_j)**3) * &
                                        !     (wl_coal_i + wl_coal_j)
                                        ! end do

                                        ! --- Vectorized redistribution to larger bins ---
                                        bin_dnum_conc_dt(max(ind_bin_i, ind_bin_j):num_bins) = &
                                            bin_dnum_conc_dt(max(ind_bin_i, ind_bin_j):num_bins) + &
                                            (1.0_dp - 0.5_dp*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                            weights_coal_curr(max(ind_bin_i, ind_bin_j):num_bins, ind_bin_i, ind_bin_j) * num_conc_coal

                                        bin_dwl_dt(max(ind_bin_i, ind_bin_j):num_bins) = &
                                            bin_dwl_dt(max(ind_bin_i, ind_bin_j):num_bins) + &
                                            (1.0_dp - 0.5_dp*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                            weights_coal_curr(max(ind_bin_i, ind_bin_j):num_bins, ind_bin_i, ind_bin_j) * &
                                            (bin_r_wet(max(ind_bin_i, ind_bin_j):num_bins)**3) / &
                                            (bin_r_wet(ind_bin_i)**3 + bin_r_wet(ind_bin_j)**3) * &
                                            (wl_coal_i + wl_coal_j)
                                    end do
                                end do
                            end if
                            if (if_brkp) then
                                do ind_ind_bin_i = 1, num_brkp_bins_curr
                                    do ind_ind_bin_j = ind_ind_bin_i, num_brkp_bins_curr
                                        ind_bin_i = bin_ind_brkp(ind_ind_bin_i)
                                        ind_bin_j = bin_ind_brkp(ind_ind_bin_j)
                                        ! Handles number concentration (conc) and liquid water content (wl) during breakup between ind_bin_i and ind_bin_j
                                        num_conc_brkp = kernels_brkp_curr (ind_bin_i, ind_bin_j) * bin_num_conc(ind_bin_i) * bin_num_conc(ind_bin_j)
                                        wl_brkp_i = kernels_brkp_curr (ind_bin_i, ind_bin_j) * bin_wl(ind_bin_i) * bin_num_conc(ind_bin_j)
                                        wl_brkp_j = kernels_brkp_curr (ind_bin_i, ind_bin_j) * bin_wl(ind_bin_j) * bin_num_conc(ind_bin_i)
                                        bin_dnum_conc_dt(ind_bin_i) = bin_dnum_conc_dt(ind_bin_i) - num_conc_brkp
                                        bin_dwl_dt(ind_bin_i) = bin_dwl_dt(ind_bin_i) - wl_brkp_i
                                        if (ind_bin_j .ne. ind_bin_i) then
                                            bin_dnum_conc_dt(ind_bin_j) = bin_dnum_conc_dt(ind_bin_j) - num_conc_brkp
                                            bin_dwl_dt(ind_bin_j) = bin_dwl_dt(ind_bin_j) - wl_brkp_j
                                        end if

                                        ! do ind_bin = 1, num_bins
                                        !     bin_dnum_conc_dt(ind_bin) = bin_dnum_conc_dt(ind_bin) + &
                                        !     (1-0.5*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                        !     weights_brkp_curr(ind_bin, ind_bin_i, ind_bin_j) * &
                                        !     num_conc_brkp

                                        !     bin_dwl_dt(ind_bin) = bin_dwl_dt(ind_bin) + &   
                                        !     (1-0.5*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                        !     weights_brkp_curr(ind_bin, ind_bin_i, ind_bin_j) * (bin_r_wet(ind_bin)**3)/(bin_r_wet(ind_bin_i)**3+bin_r_wet(ind_bin_j)**3) * &
                                        !     (wl_brkp_i + wl_brkp_j)                
                                        ! end do

                                        ! --- Vectorized redistribution to all bins ---
                                        bin_dnum_conc_dt(1:num_bins) = bin_dnum_conc_dt(1:num_bins) + &
                                            (1.0_dp - 0.5_dp*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                            weights_brkp_curr(1:num_bins, ind_bin_i, ind_bin_j) * num_conc_brkp

                                        bin_dwl_dt(1:num_bins) = bin_dwl_dt(1:num_bins) + &
                                            (1.0_dp - 0.5_dp*real(kronecker_delta(ind_bin_i, ind_bin_j), kind=dp)) * &
                                            weights_brkp_curr(1:num_bins, ind_bin_i, ind_bin_j) * &
                                            (bin_r_wet(1:num_bins)**3) / (bin_r_wet(ind_bin_i)**3 + bin_r_wet(ind_bin_j)**3) * &
                                            (wl_brkp_i + wl_brkp_j)
                                            
                                    end do
                                end do
                            end if
                            
                            if (max(abs(sum(bin_dnum_conc_dt*bin_r_wet_temp**3)),abs(sum(bin_dwl_dt)))>tolerance_wl/t_step_dcpm) then
                                print *, "ERROR: volume or mass imbalanced after collision"
                                print *, "sum(bin_dnum_conc_dt*bin_r_wet_temp**3)",sum(bin_dnum_conc_dt*bin_r_wet_temp**3)
                                print *, "sum(bin_dwl_dt)",sum(bin_dwl_dt)
                                print *, "ind_bin, bin_r_dry(ind_bin),bin_r_wet(ind_bin),bin_num_conc(ind_bin),bin_dnum_conc_dt(ind_bin), bin_wl(ind_bin), bin_dwl_dt(ind_bin)"
                                do ind_bin = 1, num_bins
                                    print *, ind_bin, bin_r_dry(ind_bin),bin_r_wet(ind_bin),bin_num_conc(ind_bin),bin_dnum_conc_dt(ind_bin), bin_wl(ind_bin), bin_dwl_dt(ind_bin)
                                end do                            
                                call stopcm1
                            end if                       
                            
                            do ind_bin = 1, num_bins
                                if (bin_num_conc(ind_bin)+bin_dnum_conc_dt(ind_bin)*t_step_coll <-tolerance_bin_num_conc) then
                                    t_step_coll_suggest = (-bin_num_conc(ind_bin)-tolerance_bin_num_conc)/bin_dnum_conc_dt(ind_bin)
                                    if (t_step_coll_suggest < t_step_coll_min) then
                                        print *, "ERROR: collision time step is too small:", t_step_coll_suggest, "minimum:", t_step_coll_min
                                        print *, "ind_bin","wet radius","number concentration", "num conc derivative", "resulting min time step"
                                        print *, ind_bin, bin_r_wet(ind_bin), bin_num_conc(ind_bin), bin_dnum_conc_dt(ind_bin), abs((-bin_num_conc(ind_bin)-tolerance_bin_num_conc)/bin_dnum_conc_dt(ind_bin))
                                        call stopcm1
                                    end if
                                    t_step_coll = min(t_step_coll, t_step_coll_suggest)
                                !$omp atomic update
                                t_step_coll_min_negc = min(t_step_coll,   t_step_coll_min_negc)
                                end if
                            end do
                            if (t_step_coll < 0) then
                                print *, "ERROR: negative collision time step"
                                call stopcm1
                            end if
                    
                            bin_num_conc = bin_num_conc+bin_dnum_conc_dt*t_step_coll
                            bin_wl = bin_wl+bin_dwl_dt*t_step_coll
                            call calculate_bin_r_dry(bin_wl, bin_r_dry_temp, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl)
                            if (if_test .and. (sum(abs(bin_dnum_conc_dt)) .ge. 100 .and. myid.eq.31) .and. .false.) then
                                print *, "output of test case: myid, myi, myj, i3d, j3d, k3d", myid, myi, myj, i3d, j3d, k3d
                                print *, "sum(bin_dnum_conc_dt*bin_r_wet_temp**3)",sum(bin_dnum_conc_dt*bin_r_wet_temp**3)
                                print *, "sum(bin_dwl_dt)", sum(bin_dwl_dt)
                                print *, "t_step_coll:", t_step_coll
                                print *, "ind_bin, bin_r_dry_old,bin_r_dry_new,bin_r_wet,bin_num_conc_new,bin_dnum_conc_dt, bin_wl_new, bin_dwl_dt"
                                print *, "num_bins", num_bins
                                do ind_bin = 1, num_bins
                                    print *, ind_bin, bin_r_dry(ind_bin),bin_r_dry_temp(ind_bin),bin_r_wet(ind_bin),bin_num_conc(ind_bin),bin_dnum_conc_dt(ind_bin), bin_wl(ind_bin), bin_dwl_dt(ind_bin)
                                end do 
                                call stopcm1
                            end if
                            bin_r_dry = bin_r_dry_temp
                            if (minval(bin_num_conc) < -tolerance_bin_num_conc*1.01) then
                                print *, "ERROR: negative number concentration after collision."
                                print *, minval(bin_num_conc), "tolerance:", tolerance_bin_num_conc
                                call stopcm1
                            end if                
                            if (minval(bin_wl) < -tolerance_bin_wl) then
                                print *, "ERROR: negative bin liquid water mixing ratio after collision."
                                print *, minval(bin_wl), "tolerance:", tolerance_bin_wl
                                call stopcm1
                            end if
                            t_coll = t_coll + t_step_coll 
                            ! print *, maxval(bin_dnum_conc_dt*t_step_coll)
                        end do
                    end if
                end if
                t_dcpm = t_dcpm + t_step_dcpm

            end do
            call calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl) ! This will also update bin_wl if necesary
            supersat = calculate_supersat (temp, p, wv)
            bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet**3d0 - bin_r_dry**3d0) *rho_w/rho_d ! Now bin_r_wet_temp = bin_r_wet
            wl = sum(bin_wl)
            if (wl > 0.02) then
                print *, "wl too large after collision: ", wl, "myid", myid, "i3d,j3d,k3d", i3d,j3d,k3d
                print *, "sum of bin_wl:", sum(bin_wl), "wv:", wv, "temp:", temp, "supersat:", supersat, "pressure:", p3d(i3d,j3d,k3d)
                do ind_bin =1, num_bins
                    print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), (4d0*pi/3d0)*bin_num_conc(ind_bin)*(bin_r_wet_temp(ind_bin)**3d0 - bin_r_dry(ind_bin)**3d0) *rho_w/rho_d
                end do
                call stopcm1
            end if
            if (minval(bin_wl) < -tolerance_bin_wl) then
                print *, "ERROR: negative bin liquid water content found after condensation and collision coalescence & breakup."
                print *, "num conc", "r wet", "r_dry", "bin wl"
                do ind_bin = 1, num_bins
                    print *, bin_num_conc(ind_bin), bin_r_wet(ind_bin), bin_r_dry(ind_bin), bin_wl(ind_bin)
                end do
                call stopcm1 
            end if
            if (minval(bin_num_conc) < -tolerance_bin_num_conc) then
                print *, "ERROR: negative bin number concentration after condensation and collision coalescence & breakup."
                call stopcm1
            end if

            if (abs(wl + wv - q3d(i3d,j3d,k3d,nqv) - q3d(i3d,j3d,k3d,nql1)) > max(tolerance_wl,0.001 *wl) .and. if_warning_wv_wl_balance) then
                print *, "WARNING: the water content is imbalanced after condensaton and collision coalescence & breakup of grid i j k at time(s) ", i3d, j3d, k3d
                print *, "z", z, "wv", wv, "wl", wl, "temp", temp, "p", p, "rho_d", rho_d, "supersat", supersat, "wv_init", wv_init
                print *, "was", temp3d(i3d,j3d,k3d), q3d(i3d,j3d,k3d,nqv) , q3d(i3d,j3d,k3d,nql1), p3d(i3d,j3d,k3d), rho3d(i3d,j3d,k3d), dcpm3d(i3d,j3d,k3d,num_bins+1)
                print *, "wv change:", wv-q3d(i3d,j3d,k3d,nqv), "wl change:", wl-q3d(i3d,j3d,k3d,nql1)
            end if
            temp_error = abs((wl - q3d(i3d,j3d,k3d,nql1))*L/Cp+temp3d(i3d,j3d,k3d)-temp)
            if (temp_error > 0.01 .and. temp_error/0.01>100) then
                print *, "myid", myid
                print *, "ERROR: the energy is imbalanced after condensaton and collision coalescence & breakup of grid i j k at time(s) ", i3d, j3d, k3d
                print *, "z", z, "wv", wv, "wl", wl, "temp", temp, "p", p, "rho_d", rho_d, "supersat", supersat, "wv_init", wv_init
                print *, "was", temp3d(i3d,j3d,k3d), q3d(i3d,j3d,k3d,nqv) , q3d(i3d,j3d,k3d,nql1), p3d(i3d,j3d,k3d), rho3d(i3d,j3d,k3d), dcpm3d(i3d,j3d,k3d,num_bins+1)
                print *, "wv change:", wv-q3d(i3d,j3d,k3d,nqv), "wl change:", wl-q3d(i3d,j3d,k3d,nql1)
                print *, "expected temperature change:", (wl - q3d(i3d,j3d,k3d,nql1))*L/Cp, "actual temperature change:", temp-temp3d(i3d,j3d,k3d)
                call stopcm1
            end if

            ! calculate effective radius
            if( dorad ) then
                if_grid_wet = .false.
                numerator = 0d0
                denominator = 0d0
                do ind_bin = 1, num_bins
                    if (bin_wl(ind_bin) > tolerance_bin_wl) then
                        numerator = numerator + (bin_r_wet(ind_bin)*1d6)**3d0 * bin_num_conc(ind_bin)
                        denominator = denominator + (bin_r_wet(ind_bin)*1d6)**2d0 * bin_num_conc(ind_bin)
                        if_grid_wet = .true.
                    end if
                end do
                if (if_grid_wet .eqv. .true.) then
                    effc_curr = numerator/denominator
                    if( radopt.eq.1 )then
                        effc(i3d,j3d,k3d) = max(effc_curr,25d0)
                    elseif( radopt.eq.2 )then
                        effc(i3d,j3d,k3d) = max(effc_curr,25d0)*1e-6
                    end if
                end if
            end if

            wv = wv_init - (wl-wl_init)
            supersat = calculate_supersat (temp, p, wv)
            temp3d(i3d,j3d,k3d) = temp
            ! if (myid.eq.240 .and. i3d.eq.1 .and. j3d.eq.2 .and. k3d.eq.1) then
            !     print *, "At sample point after dcpm calculation:", " wl:", wl, "wv:", wv, "temp:",temp, "p:", p
            ! end if

            q3d(i3d,j3d,k3d,nqv) = wv
            q3d(i3d,j3d,k3d,nql1) = wl
            q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1) = bin_wl*scale_bin_wl
            q3d(i3d,j3d,k3d,nnc1+num_bins:nnc2) = bin_num_conc
            dcpm3d(i3d,j3d,k3d,1:num_bins) = bin_r_dry
            dcpm3d(i3d,j3d,k3d,num_bins+1) = supersat
            adddiag3d(i3d,j3d,k3d,2) = supersat
            p3d(i3d,j3d,k3d) = p
            pp3d(i3d,j3d,k3d) = (p/p00)**(1.0/cpdrd)-pi0(i3d,j3d,k3d)
            th3d(i3d,j3d,k3d) = temp/(p/p00)**(1.0/cpdrd)-th0(i3d,j3d,k3d)
            rho3d(i3d,j3d,k3d) = calculate_rho_d (temp, p, wv)

        enddo
        enddo
        enddo
        ! close (240311)
        !$omp end parallel do

        if (if_did_cond_stat .eqv. .false.) then
            ! print *, "No condensation performed in dcpm this time (no bin is activated)"
            if_dcpm_activated_myid = .false.
        else
            if_dcpm_activated_myid = .true.
            ! print *, "the minimum condensation step among all bins at dcpm_calculate:", t_step_cond_min_stat
            ! print *, "the active grid number:", num_grid_active_stat, "total grid number:", ni*nj*nk, "times to update by distrbuting mass:", num_dist_stat, "times to update by Runge-Kutta:", num_rk_stat
            ! print *, "suggested LES time step:", t_step_LES_suggest
        end if
        
        q3d_dcpm = q3d

        wl_max = maxval(q3d(1:ni,1:nj,1:nk,nql1))
        wl_min = minval(q3d(1:ni,1:nj,1:nk,nql1))
        wi_max = maxval(q3d(1:ni,1:nj,1:nk,nqs1))
        wi_min = minval(q3d(1:ni,1:nj,1:nk,nqs1))
        supersat_max = maxval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1))
        supersat_min = minval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1))
        num_conc_max = maxval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))
        num_conc_min = minval(q3d(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))
        r_drop_max = 0d0
        do ind_bin = num_bins, 1, -1
            if (maxval(q3d(1:ni,1:nj,1:nk,nnc1+ind_bin-1)) > tolerance_bin_wl*scale_bin_wl) then
                r_drop_max = bin_r_wet(ind_bin)
                exit
            end if
        end do
        t_step_coll_min_negc_all = t_step_coll_min_negc
        wl_max_all = wl_max
        wl_min_all = wl_min
        wi_max_all = wi_max
        wi_min_all = wi_min
        supersat_max_all = supersat_max
        supersat_min_all = supersat_min
        num_conc_max_all = num_conc_max
        num_conc_min_all = num_conc_min
        r_drop_max_all = r_drop_max
        num_rk_stat_all = num_rk_stat
        num_dist_stat_all = num_dist_stat
        if_dcpm_activated = if_dcpm_activated_myid
        
#ifdef MPI
    call MPI_REDUCE(wl_max,wl_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wl_min,wl_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wi_max,wi_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(wi_min,wi_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(supersat_max,supersat_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(supersat_min,supersat_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_conc_max,num_conc_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_conc_min,num_conc_min_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(r_drop_max,r_drop_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(if_dcpm_activated_myid,if_dcpm_activated,1,MPI_LOGICAL, MPI_LOR,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(t_step_coll_min_negc,t_step_coll_min_negc_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_rk_stat,num_rk_stat_all,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(num_dist_stat,num_dist_stat_all,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD,ierr)
#endif    
        if (myid .eq. 0) then
            if (if_dcpm_activated .eqv. .true.) then
                print *, "After bin microphysics:"
                print *, "begin variable", " wl (g/kg) ", " min ", wl_min_all*1e3," max ", wl_max_all*1e3, " end variable"
                print *, "begin variable", " wi (g/kg) ", " min ", wi_min_all*1e3," max ", wi_max_all*1e3, " end variable"
                print *, "begin variable", " number concentration (cm^{-3}) ", " min ", num_conc_min_all*1e-6," max ", num_conc_max_all*1e-6, " end variable"
                print *, "begin variable", " drop radius (mm) ", " min ", 0," max ", r_drop_max_all * 1e3, " end variable"
                print *, "begin variable", " total rk time steps: ", num_rk_stat_all, " end variable"
                print *, "begin variable", " total distribution time steps: ", num_dist_stat_all, " end variable"
            else
                print *, "No condensation performed in dcpm this time (no bin is activated)"
                print *, "begin variable", " supersaturation ", " min ", supersat_min_all," max ", supersat_max_all, " end variable"
            end if
            print *, "Exit DCPM"
        end if

        call MPI_BARRIER (MPI_COMM_WORLD, ierr)
        call dcpm_postcheck (q3d, dcpm3d,pp3d,th3d,p3d,temp3d,rho3d)

        ! print *,    "extreme values in q3d_dcpm (1:ni, 1:nk, :1:nj, :) after bin microphysics:", "temp", minval(temp3d(1:ni,1:nj,1:nk)), maxval(temp3d(1:ni,1:nj,1:nk)),&
        !             "wv", minval(q3d_dcpm(1:n i,1:nj,1:nk,nqv)), maxval(q3d_dcpm(1:ni,1:nj,1:nk,nqv)), "wl", minval(q3d_dcpm(1:ni,1:nj,1:nk,nql1)), maxval(q3d_dcpm(1:ni,1:nj,1:nk,nql1)), &
        !             "bin_wl", minval(q3d_dcpm(1:ni,1:nj,1:nk,nnc1:nnc1+num_bins-1)/scale_bin_wl), maxval(q3d_dcpm(1:ni,1:nj,1:nk,nnc1:nnc1+num_bins-1)/scale_bin_wl), &
        !             "p", minval(p3d(1:ni,1:nj,1:nk)), maxval(p3d(1:ni,1:nj,1:nk)), "theta", minval(th3d(1:ni,1:nj,1:nk)), maxval(th3d(1:ni,1:nj,1:nk)), &
        !             "rho", minval(rho3d(1:ni,1:nj,1:nk)), maxval(rho3d(1:ni,1:nj,1:nk)), "supersat", minval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1)), maxval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1)),&
        !             "number concentraton", minval(q3d_dcpm(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2)), maxval(q3d_dcpm(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2))   
        ! print *,    "extreme values in q3d_dcpm (:, :, ::, :) after bin microphysics:", "temp", minval(temp3d(:,:,:)), maxval(temp3d(:,:,:)),&
        !             "wv", minval(q3d_dcpm(:,:,:,nqv)), maxval(q3d_dcpm(:,:,:,nqv)), "wl", minval(q3d_dcpm(:,:,:,nql1)), maxval(q3d_dcpm(:,:,:,nql1)), &
        !             "bin_wl", minval(q3d_dcpm(:,:,:,nnc1:nnc1+num_bins-1)/scale_bin_wl), maxval(q3d_dcpm(:,:,:,nnc1:nnc1+num_bins-1)/scale_bin_wl), &
        !             "p", minval(p3d(:,:,:)), maxval(p3d(:,:,:)), "theta", minval(th3d(:,:,:)), maxval(th3d(:,:,:)), &
        !             "rho", minval(rho3d(:,:,:)), maxval(rho3d(:,:,:)), "supersat", minval(dcpm3d(:,:,:,num_bins+1)), maxval(dcpm3d(:,:,:,num_bins+1)),&
        !             "number concentraton", minval(q3d_dcpm(:,:,:,nnc1+num_bins:nnc2)), maxval(q3d_dcpm(:,:,:,nnc1+num_bins:nnc2))  

        ! call cpu_time ( t_cpu_curr )
        ! write ( *, * ) 'Elapsed CPU time (count from dcpm initialization) = ', t_cpu_curr - t_cpu_init
        ! call system_clock ( t_clock_curr, clock_rate_dcpm, clock_max_dcpm )
        ! write ( *, * ) 'Elapsed real time (count from dcpm initialization) = ', real ( t_clock_curr - t_clock_init ) / real ( clock_rate_dcpm )

        ! print *, "exit dcpm_calculate"

    end subroutine dcpm_calculate

    subroutine dcpm_fallout (dt, q3d, mh, rain, prate, rho, temp3d)
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke,numq), intent(inout) :: q3d
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke), intent(in) :: rho, mh, temp3d
        real (kind=sp), dimension(ib:ie,jb:je,nrain), intent(inout) :: rain
        real (kind=sp), intent(inout), dimension(ib:ie,jb:je) :: prate
        real (kind=sp), intent(in) :: dt
        real (kind=sp) :: dtfall
        real (kind=dp), dimension(1:nk+1):: rq, rn
        real (kind=dp), dimension(1:nk) :: vice_1d
        real (kind=dp) :: crmax, v_term
        real (kind=dp), dimension(1:nk) :: v_term_column
        real (kind=dp) :: v_term_max, prate_max, rain_avg, v_term_max_all, prate_max_all, rain_avg_all
        integer :: n, nfall, ind_bin_drop_start, i3d, j3d, k3d, nr, ind_bin

        prate = 0d0
        ! Fallout for ice


        do i3d=1,ni
        do j3d=1,nj
            do k3d=1,nk
                ! Fall velocity by Mitchell et al. 2007
                if (q3d(i3d,j3d,k3d,nqs1)>=tolerance_wi) then
                    vice_1d(k3d) = 0.01*(1.0121*(temp3d(i3d,j3d,k3d)-273.15)+6.6303*log(q3d(i3d,j3d,k3d,nqs1)*rho(i3d,j3d,k3d)*1e6)+82.082)
                else
                    vice_1d(k3d) = 0
                end if
            end do
            crmax = 0.0
            do k3d=1,nk
                crmax = max( crmax , vice_1d(k3d)*dt*rdz*mh(i3d,j3d,k3d) )
            enddo
            nfall = max( 1 , int(crmax+1.0) )
            do k3d=(nk-1),1,-1
                if( q3d(i3d,j3d,k3d,nqs1).lt.tolerance_wi ) then 
                    vice_1d(k3d) = vice_1d(k3d+1)
                end if
            enddo
            dtfall=dt/nfall
            do n=1,nfall
                do k3d=1,nk
                rq(k3d)=rho(i3d,j3d,k3d)*vice_1d(k3d)*max(0.0,q3d(i3d,j3d,k3d,nqs1))
                enddo
                rq(nk+1)=0.0
                do k3d=1,nk
                    q3d(i3d,j3d,k3d,nqs1)=q3d(i3d,j3d,k3d,nqs1)+dtfall*(rq(k3d+1)-rq(k3d))   &
                                                *rdz*mh(i3d,j3d,k3d)/rho(i3d,j3d,k3d)
                enddo
                do nr=1,nrain
                rain(i3d,j3d,nr)=rain(i3d,j3d,nr)+0.1*dtfall*rq(1)
                                        ! mult by 100 cm/m  (convert from m to cm)
                                        ! div by 1000 kg/m^3  (density of water)
                enddo
                prate(i3d,j3d) = prate(i3d,j3d)+rq(1)
            enddo
            if (q3d(i3d,j3d,1,nqs1) > tolerance_wi) then
                print *, "ERROR: Ice shows up at ground level."
                print *, "k3d ", "rq ", "vice ", "qice"
                do k3d=1, nk
                    print *, k3d, rq(k3d), q3d(i3d,j3d,k3d,nqs1)
                end do 
                call stopcm1
            end if
        enddo
        enddo

        ! Fallout for water

        do ind_bin =1, num_bins
            if (bin_r_wet(ind_bin) > drop_radius_min) then
                ind_bin_drop_start = ind_bin
                exit
            end if
        end do
        v_term_max = 0
        do i3d=1,ni
        do j3d=1,nj
            do ind_bin = ind_bin_drop_start, num_bins
                ! Enable sediments
                if (maxval(q3d(i3d,j3d,1:nk,nnc1+num_bins+ind_bin-1)) < tolerance_bin_num_conc) then
                    cycle
                end if            
                do k3d=1,nk
                    if (iv_fall .eq. 1) then
                        v_term = calculate_terminal_velocity_Best (bin_r_wet(ind_bin))
                    else if (iv_fall .eq. 2) then
                        v_term = calculate_terminal_velocity_Lhermitte (bin_r_wet(ind_bin),(rho(i3d,j3d,1)*1.5-rho(i3d,j3d,2)*0.5)/rho(i3d,j3d,k3d))
                    else
                        print *, "ERROR: iv_fall option unknown! iv_fall=", iv_fall
                        call stopcm1
                    end if                    
                    v_term_column(k3d) = v_term
                end do
                v_term_max = max(v_term_max, maxval(v_term_column))
                crmax = max(0d0, maxval(v_term_column)*dt*rdz*maxval(mh(i3d,j3d,1:nk)))
                nfall = max( 1 , int(crmax+1.0) )
                dtfall = dt/nfall
                do n=1,nfall
                    do k3d=1,nk
                        v_term = v_term_column(k3d)
                        rq(k3d)=rho(i3d,j3d,k3d)*v_term*max(0.0,q3d(i3d,j3d,k3d,nnc1+ind_bin-1))
                        rn(k3d)=rho(i3d,j3d,k3d)*v_term*max(0.0,q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1))
                    end do
                    rq(nk+1)=0.0
                    do k3d=1,nk
                        q3d(i3d,j3d,k3d,nnc1+ind_bin-1)=q3d(i3d,j3d,k3d,nnc1+ind_bin-1) &
                            +dtfall*(rq(k3d+1)-rq(k3d))*rdz*mh(i3d,j3d,k3d)/rho(i3d,j3d,k3d)
                        q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1)=q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1) &
                            +dtfall*(rn(k3d+1)-rn(k3d))*rdz*mh(i3d,j3d,k3d)/rho(i3d,j3d,k3d)
                    end do
                    do nr=1,nrain
                        rain(i3d,j3d,nr)=rain(i3d,j3d,nr)+0.1*dtfall*rq(1)/scale_bin_wl
                                            ! mult by 100 cm/m  (convert from m to cm)
                                            ! div by 1000 kg/m^3  (density of water)
                    end do  
                end do
                prate(i3d,j3d) = prate(i3d,j3d) + rq(1)/scale_bin_wl
                if (prate(i3d,j3d) > 1e-3*10 .and. if_warning_tolerance) then
                    print *, "WARNING: precipitation rate too large"                 
                    print *, "rate:", prate(i3d,j3d), rq(1), v_term, sum(q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1))
                    ! call stopcm1
                end if
            end do
        end do
        end do

    

        prate_max = maxval(prate(1:ni,1:nj))

        rain_avg = sum(rain(1:ni,1:nj,1))/ni/nj
        prate_max_all = prate_max
        rain_avg_all = rain_avg
#ifdef MPI
    call MPI_REDUCE(prate_max,prate_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(v_term_max,v_term_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(rain_avg,rain_avg_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    rain_avg_all = rain_avg_all/numprocs
#endif   
        if (myid .eq. 0) then
            print *, "begin variable", " maximum precipitation rate (kg/m2/s):", prate_max_all, " end variable"
            print *, "begin variable", " accumulated surface rainfall (cm):", rain_avg_all, " end variable"
            print *, "begin variable", " maximum terminal velocity (m/s):", v_term_max_all, " end variable"
        end if      


    end subroutine dcpm_fallout

    subroutine calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl)
        real (kind=dp), dimension(num_bins), intent(inout) :: bin_wl, bin_r_dry, bin_num_conc   
        real (kind=dp), dimension(num_bins), intent(in) :: bin_r_wet_temp  
        real (kind=dp), intent(in) :: rho_d
        real (kind=dp), intent(inout) :: wv
        real (kind=dp), intent(out) :: wl    
        real (kind=dp) :: wl_bin_defecit, wl_bin_excess, r_dry_cubic
        integer :: ind_bin
        wl = sum(bin_wl)
        do ind_bin = 1, num_bins

            if (bin_num_conc(ind_bin) < tolerance_bin_num_conc) then
                bin_num_conc(ind_bin) = 0
                bin_wl(ind_bin) = 0
                if (bin_r_wet_temp(ind_bin) > ccn_radius_max) then
                    bin_r_dry(ind_bin) = 0
                else
                    bin_r_dry(ind_bin) = bin_r_wet_temp(ind_bin)
                end if
            else
                r_dry_cubic = (bin_r_wet_temp(ind_bin)**3d0-bin_wl(ind_bin)*3d0/4d0/pi/bin_num_conc(ind_bin)*rho_d/rho_w)
                if (r_dry_cubic <= (ccn_radius_min-tolerance_bin_r)**3d0) then
                    bin_r_dry(ind_bin) = ccn_radius_min
                    wl_bin_excess = bin_wl(ind_bin) - (bin_r_wet_temp(ind_bin)**3d0-bin_r_dry(ind_bin)**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc(ind_bin)
                    if (abs(wl_bin_excess) > tolerance_bin_wl) then
                        if (if_warning_tolerance) then
                            print *, "WARNING: wl excess exceeds tolerance."
                            print *, "wl excess:", wl_bin_excess, "bin wl tolerance:", tolerance_bin_wl
                            print *, ind_bin, bin_r_wet_temp(ind_bin), bin_wl(ind_bin), bin_num_conc(ind_bin), bin_r_dry(ind_bin), rho_d, rho_w
                        end if
                    end if
                    bin_wl(ind_bin) = (bin_r_wet_temp(ind_bin)**3d0-bin_r_dry(ind_bin)**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc(ind_bin)
                else
                    bin_r_dry(ind_bin) = r_dry_cubic**(1d0/3d0)  
                    if (abs((bin_r_wet_temp(ind_bin)**3d0-bin_r_dry(ind_bin)**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc(ind_bin) - bin_wl(ind_bin)) > tolerance_bin_wl) then
                        print *, "ERROR: the updating didn't work!"
                        print *, bin_wl(ind_bin), abs((bin_r_wet_temp(ind_bin)**3d0-bin_r_dry(ind_bin)**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc(ind_bin) - bin_wl(ind_bin)) , bin_r_wet_temp(ind_bin), bin_r_dry(ind_bin), bin_num_conc(ind_bin)     
                        call stopcm1
                    end if                                
                end if
            end if


            if (ieee_is_nan(bin_r_dry(ind_bin)) .or. bin_r_dry(ind_bin)<0) then
                print *, "ERROR: dry radius is nan or negative"
                print *, "ind_bin:", ind_bin
                print *, bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin)
                print *, "wl after advection & diffusion:", bin_wl(ind_bin), "maximum allowance:", bin_r_wet_temp(ind_bin)**3d0*4d0/3d0*pi*rho_w/rho_d*bin_num_conc(ind_bin)
                print *, ind_bin, bin_r_wet_temp(ind_bin), bin_wl(ind_bin), bin_num_conc(ind_bin), bin_r_dry(ind_bin), rho_d, rho_w
                print *,  bin_r_wet_temp(ind_bin)**3d0-bin_wl(ind_bin)*3d0/4d0/pi/bin_num_conc(ind_bin)*rho_d/rho_w
                call stopcm1
            end if               
            if (bin_r_dry(ind_bin) >= bin_r_wet_temp(ind_bin)+tolerance_bin_r) then
                print *, "ERROR: dry radius is larger than wet radius"
                print *, bin_r_dry(ind_bin), bin_r_wet_temp(ind_bin), bin_num_conc(ind_bin), bin_wl(ind_bin), rho_d, rho_w
                call stopcm1
            end if

        end do
        wv = wv - (sum((bin_r_wet_temp**3d0-bin_r_dry**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc)-wl)
        wl = sum((bin_r_wet_temp**3d0-bin_r_dry**3d0)*4d0/3d0*pi*rho_w/rho_d*bin_num_conc)
        ! if (minval(bin_r_dry) < 0) then
        !     print *, "ERROR: negative dry radius after advection/difusion"
        !     do ind_bin = 1, num_bins
    end subroutine calculate_bin_r_dry

    subroutine dcpm_preprocess (q3d, dcpm3d, th0, th3d, pi0, pp3d)
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(inout) :: q3d
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,num_bins+1), intent(inout) :: dcpm3d
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem), intent(in) :: th0, pi0, pp3d
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem), intent(inout) :: th3d
        integer :: i3d, j3d, k3d, ind_bin, ind_bin_err
        real (kind=dp) :: wl_old, wl_new
        real (kind=dp) :: temp
        do i3d=1,ni
        do j3d=1,nj
        do k3d=1,nk
            temp = (th0(i3d,j3d,k3d)+th3d(i3d,j3d,k3d))*(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d))
            if (ieee_is_nan(q3d(i3d,j3d,k3d,nqv))) then
                print *, "got nan for water vapor mixing ratio after advection and diffusion"
                print *, q3d(i3d,j3d,k3d,nqv)
                print *, q3d(i3d,j3d,k3d,:)
                call stopcm1
            end if
            if (ieee_is_nan(q3d(i3d,j3d,k3d,nql1))) then
                print *, "got nan for liquid water mixing ratio after advection and diffusion"
                print *, q3d(i3d,j3d,k3d,nql1)
                print *, q3d(i3d,j3d,k3d,:)
                call stopcm1
            end if
            if (ieee_is_nan(q3d(i3d,j3d,k3d,nqs1))) then
                print *, "got nan for ice water mixing ratio after advection and diffusion"
                print *, q3d(i3d,j3d,k3d,nqs1)
                print *, q3d(i3d,j3d,k3d,:)
                call stopcm1
            end if

            wl_old = sum(q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl)
            ! bin liquid water mixing ratio
            do ind_bin = 1, num_bins
                if (ieee_is_nan(q3d(i3d,j3d,k3d,nnc1+ind_bin-1))) then
                    print *, "got nan for bin liquid wter mixing ratio after advection and diffusion"
                    do ind_bin_err = 1, num_bins
                        print *, ind_bin_err, q3d(i3d,j3d,k3d,nnc1+ind_bin_err-1), "was:", q3d_dcpm(i3d,j3d,k3d,nnc1+ind_bin_err-1)/scale_bin_wl
                        print *, ind_bin_err, q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin_err-1), "was:", q3d_dcpm(i3d,j3d,k3d,nnc1+num_bins+ind_bin_err-1)
                    end do
                    call stopcm1
                end if            
                if ( q3d(i3d,j3d,k3d,nnc1+ind_bin-1)<-tolerance_bin_wl*scale_bin_wl) then
                    if (if_warning_tolerance) then
                        print *, "WARNING: got negative bin liquid wter mixing ratio after advection and diffusion:", q3d(i3d,j3d,k3d,nnc1+ind_bin-1)/scale_bin_wl, "was ", q3d_dcpm(i3d,j3d,k3d,nnc1+ind_bin-1)/scale_bin_wl
                        print *, "bin tolerance:", -tolerance_bin_wl
                        Print *, "cell index:", i3d, j3d, k3d
                        ! do ind_bin_err = 1, num_bins
                        !     print *, ind_bin_err, q3d(i3d,j3d,k3d,nnc1+ind_bin_err-1), "was:", q3d_dcpm(i3d,j3d,k3d,nnc1+ind_bin_err-1)/scale_bin_wl
                        !     print *, ind_bin_err, q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin_err-1), "was:", q3d_dcpm(i3d,j3d,k3d,nnc1+num_bins+ind_bin_err-1)
                        ! end do
                        ! stop
                    end if
                end if
                if ( q3d(i3d,j3d,k3d,nnc1+ind_bin-1)>1d-2*scale_bin_wl) then
                    if (if_warning_tolerance) then
                        print *, "WARNING: got too large bin liquid water mixing ratio after advection and diffusion"
                        print *, "processor", myid
                        print *, i3d, j3d, k3d
                        print *, q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl
                    end if
                end if
                
                !! Need to remove tiny values (otherwise the bin computation would be super heavy)
                if (q3d(i3d,j3d,k3d,nnc1+ind_bin-1)<tolerance_wl*scale_bin_wl) then
                    q3d(i3d,j3d,k3d,nnc1+ind_bin-1) = 0.0
                    ! if (bin_r_wet(ind_bin)>ccn_radius_max) then
                    !     q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1) = 0d0
                    ! end if
                end if
            end do          
       
            ! bin number concentration
            do ind_bin = 1, num_bins
                if (ieee_is_nan(q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1))) then
                    print *, "got nan for number concentration after advection and diffusion"
                    print *, q3d(i3d,j3d,k3d,nnc1:nnc2)
                    call stopcm1
                end if            
                if ( ((q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1)<-tolerance_bin_num_conc) .and. q3d(i3d,j3d,k3d,nnc1+ind_bin-1)>tolerance_bin_wl*scale_bin_wl) .and. if_warning_tolerance) then 
                    print *, "WARNING: got negative number concentration after advection and diffusion:", q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1), "was", q3d_dcpm(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1)
                    print *, "bin liquid wter mixing ratio after advection and diffusion:", q3d(i3d,j3d,k3d,nnc1+ind_bin-1)/scale_bin_wl, "was ", q3d_dcpm(i3d,j3d,k3d,nnc1+ind_bin-1)/scale_bin_wl
                    print *, "myid", myid, "grid index:", i3d, j3d, k3d
                    ! print *, "the neighboring grid change", " left:", q3d_dcpm(i3d-1,j3d,k3d,nnc1+num_bins+ind_bin-1), "to", q3d(i3d-1,j3d,k3d,nnc1+num_bins+ind_bin-1)
                    ! print *, "the neighboring grid change", " right:", q3d_dcpm(i3d+1,j3d,k3d,nnc1+num_bins+ind_bin-1), "to", q3d(i3d+1,j3d,k3d,nnc1+num_bins+ind_bin-1)
                    ! print *, "the neighboring grid change", " front:", q3d_dcpm(i3d,j3d-1,k3d,nnc1+num_bins+ind_bin-1), "to", q3d(i3d,j3d-1,k3d,nnc1+num_bins+ind_bin-1)                    
                    ! print *, "the neighboring grid change", " back:", q3d_dcpm(i3d,j3d+1,k3d,nnc1+num_bins+ind_bin-1), "to", q3d(i3d,j3d+1,k3d,nnc1+num_bins+ind_bin-1)
                    ! print *, "the neighboring grid change", " bot:", q3d_dcpm(i3d,j3d,k3d-1,nnc1+num_bins+ind_bin-1), "to", q3d(i3d,j3d,k3d-1,nnc1+num_bins+ind_bin-1)                    
                    ! print *, "the neighboring grid change", " top:", q3d_dcpm(i3d,j3d,k3d+1,nnc1+num_bins+ind_bin-1), "to", q3d(i3d,j3d,k3d+1,nnc1+num_bins+ind_bin-1)
                    ! stop
                end if

                if (q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1) < tolerance_bin_num_conc) then
                    q3d(i3d,j3d,k3d,nnc1+num_bins+ind_bin-1) = 0d0
                    q3d(i3d,j3d,k3d,nnc1+ind_bin-1) = 0.0
                end if
            end do          

            ! liquid water mixing ratio and temperature
            wl_new = sum(q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl)
            q3d(i3d,j3d,k3d,nql1) = wl_new
            q3d(i3d,j3d,k3d,nqv) = q3d(i3d,j3d,k3d,nqv) - (wl_new-wl_old)
            temp = temp+L*(wl_new-wl_old)/C_p
            th3d(i3d,j3d,k3d) = temp/(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d))-th0(i3d,j3d,k3d)
            ! water vapor mixing ratio
            if (q3d(i3d,j3d,k3d,nqv) < min(0d0,tolerance_wv)) then
                q3d(i3d,j3d,k3d,nqv) = 0d0
            end if

            ! ice water mixing ratio
            if (q3d(i3d,j3d,k3d,nqs1) < 0d0) then
                q3d(i3d,j3d,k3d,nqs1) = 0d0
            end if

            ! liquid water mixing ratio
            if (q3d(i3d,j3d,k3d,nql1) < 0d0) then
                q3d(i3d,j3d,k3d,nql1) = 0d0
            end if 

        end do 
        end do
        end do

    end subroutine

    subroutine dcpm_postcheck (q3d, dcpm3d,pp3d,th3d,p3d,temp3d,rho3d)
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(in) :: q3d
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,num_bins+1), intent(in) :: dcpm3d
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke), intent(in) :: pp3d,th3d,p3d,temp3d,rho3d
        integer :: i3d, j3d, k3d, ind_bin

        if ( maxval(abs(bin_r_wet - bin_r_wet_check))      >= tolerance_bin_r  .or. &
            maxval(abs(bin_r_wet_grid - bin_r_wet_grid_check)) >= tolerance_bin_r ) then
        print *, "ERROR: bin_r_wet or bin_r_wet_grid was modified somewhere in the calculation"
        print *, "myid, myi, myj", myid, myi, myj
        print *, "ind_bin, r_wet, left grid, right grid"
        do ind_bin = 1, num_bins
            print *, ind_bin, bin_r_wet(ind_bin), bin_r_wet_grid(ind_bin), bin_r_wet_grid(ind_bin+1)
        end do
        call stopcm1
        end if

        do i3d=1,ni
        do j3d=1,nj
        do k3d=1,nk
            if (ieee_is_nan(temp3d(i3d,j3d,k3d))) then
                print *, "temp is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(p3d(i3d,j3d,k3d))) then
                print *, "p is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(pp3d(i3d,j3d,k3d))) then
                print *, "pp is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(th3d(i3d,j3d,k3d))) then
                print *, "th is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(rho3d(i3d,j3d,k3d))) then
                print *, "rho is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(q3d(i3d,j3d,k3d,nqv))) then
                print *, "qv is nan after bin microphysics"
                call stopcm1
            end if
            if (ieee_is_nan(q3d(i3d,j3d,k3d,nql1))) then
                print *, i3d, j3d, k3d
                print *, "ql is nan after bin microphysics"
                call stopcm1
            end if
            do ind_bin = 1, num_bins
                if ( q3d(i3d,j3d,k3d,nnc1+ind_bin-1)<-tolerance_bin_num_conc) then
                    if (q3d(i3d,j3d,k3d,nnc1+ind_bin-1)<-tolerance_bin_num_conc*10) then
                        print *, "ERROR: got negative number concentration after bin microphysics (<-10)"
                        print *, q3d(i3d,j3d,k3d,nnc1:nnc2)
                        call stopcm1
                    else
                        print *, "WARNING: got negative number concentration after bin microphysics"
                    end if
                end if
            end do 
            do ind_bin = 1, num_bins
                if (ieee_is_nan(q3d(i3d,j3d,k3d,nnc1+ind_bin-1))) then
                    print *, "got nan for number concentration after bin microphysics"
                    print *, q3d(i3d,j3d,k3d,nnc1:nnc2)
                    call stopcm1
                end if
            end do
            do ind_bin = 1, num_bins
                if ( q3d(i3d,j3d,k3d,nnc1+ind_bin-1)<-tolerance_bin_wl*scale_bin_wl) then
                    print *, "got negative bin liquid water mixing ratio after bin microphysics"
                    print *, q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)
                    call stopcm1
                end if
                if ( q3d(i3d,j3d,k3d,nnc1+ind_bin-1)>1d-2*scale_bin_wl) then
                    print *, "processor", myid
                    print *, "got too large bin liquid water mixing ratio after bin microphysics"
                    print *, i3d, j3d, k3d
                    print *, q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)/scale_bin_wl
                    call stopcm1
                end if
            end do 
            do ind_bin = 1, num_bins
                if (ieee_is_nan(q3d(i3d,j3d,k3d,nnc1+ind_bin-1))) then
                    print *, "got nan for bin liquid water mixing ratio after bin microphysics"
                    print *, q3d(i3d,j3d,k3d,nnc1:nnc1+num_bins-1)
                    call stopcm1
                end if
            end do
        end do 
        end do
        end do
    end subroutine
!#####################     Condensation/Evaporation      ###################

    subroutine calculate_bin_deriv (bin_deriv, bin_r_wet_curr, bin_shift, bin_r_dry, bin_if_act, bin_num_conc, wv_curr, temp_curr, p, num_cond_bins) 
        ! function to calculate radius derivatives of rk4
        real (kind=dp), dimension(num_bins), intent(in) :: bin_r_wet_curr, bin_r_dry, bin_num_conc, bin_shift
        logical, dimension(num_bins), intent(in) :: bin_if_act
        real (kind=dp), dimension(num_bins), intent(out) :: bin_deriv
        real (kind=dp), dimension(num_bins) :: bin_r_wet_predict
        real (kind=dp), intent(in) :: wv_curr, temp_curr, p
        real (kind=dp) :: rho_d_curr, wl_shift, wv, temp, supersat, rho_d, r_dry, r_wet, e_s
        logical :: if_act
        real (kind=dp) :: thermal_conductivity_air, diffusivity_water_vapor, G_a, G_b, coeff_growth, dr_dt
        integer, intent(in) :: num_cond_bins
        integer :: ind_bin
        if (if_square_radius_scheme) then
            bin_r_wet_predict = bin_r_wet_curr**2+bin_shift
            do ind_bin = 1, num_bins
                if (bin_r_wet_predict(ind_bin) <= bin_r_dry(ind_bin)**2d0) then
                    bin_r_wet_predict(ind_bin) = bin_r_dry(ind_bin)
                else 
                    bin_r_wet_predict(ind_bin) = sqrt(bin_r_wet_predict(ind_bin))
                end if
            end do                       
        else
            bin_r_wet_predict = bin_r_wet_curr+bin_shift
            do ind_bin = 1, num_bins
                if (bin_r_wet_predict(ind_bin) <= bin_r_dry(ind_bin)) then
                    bin_r_wet_predict(ind_bin) = bin_r_dry(ind_bin)
                end if
            end do
        end if

        rho_d_curr = calculate_rho_d (temp_curr, p, wv_curr)        
        wl_shift = sum((4d0*pi/3d0)*bin_num_conc*(bin_r_wet_predict **3d0 - bin_r_wet_curr**3d0) *rho_w/rho_d_curr) 
        wv = wv_curr-wl_shift
        temp = temp_curr + L*wl_shift/C_p
        supersat = calculate_supersat (temp, p, wv) 
        rho_d = calculate_rho_d (temp, p, wv)      
        e_s = calculate_e_s(temp)
        ! print * , wv, temp, supersat, rho_d, e_s

        bin_deriv = 0
        do ind_bin=1, num_cond_bins
            r_dry = bin_r_dry(ind_bin)       
            r_wet = bin_r_wet_predict(ind_bin)   
            if_act = bin_if_act(ind_bin)
            if (if_act) then
                thermal_conductivity_air = calculate_thermal_conductivity_air(temp, rho_d, r_wet)
                diffusivity_water_vapor = calculate_diffusivity_water_vapor(temp, r_wet, p)
                G_a = (rho_w*R*temp)/e_s/diffusivity_water_vapor/M_w
                G_b = (L*rho_w*((L*M_w/R/temp)-1d0))/thermal_conductivity_air/temp
                coeff_growth = 1d0/(G_a + G_b)
                if (if_square_radius_scheme) then
                    bin_deriv(ind_bin) = 2d0*coeff_growth*(supersat - calculate_equilibrium_supersat (r_wet, r_dry, temp, aer_kappa))
                else
                    bin_deriv(ind_bin) = (coeff_growth/r_wet)*(supersat - calculate_equilibrium_supersat (r_wet, r_dry, temp, aer_kappa))
                end if
                if (ieee_is_nan(bin_deriv(ind_bin))) then
                    print *, "ERROR: got nan bin radius derivative."
                    print *, thermal_conductivity_air, diffusivity_water_vapor, G_a, G_b, coeff_growth, supersat, calculate_equilibrium_supersat (r_wet, r_dry, temp, aer_kappa), r_wet
                    call stopcm1
                end if
            else
                bin_deriv(ind_bin) = 0
            end if                  
        end do    
        ! print *, bin_if_act
        ! stop
    end subroutine calculate_bin_deriv

    subroutine determine_t_step_cond_ratio (t_step_cond_ratio_radius, t_step_cond_ratio_temperature, t_step_cond_ratio_moisture, ind_bin_radius_control, bin_r_wet_predict, bin_r_wet_temp, bin_r_dry, bin_num_conc, wv, wv_predict, wv_sat, wv_sat_predict, wl, wl_predict, temp, temp_predict, num_cond_bins)                          
        real(kind=dp), dimension(num_bins), intent(in) :: bin_r_wet_predict, bin_r_wet_temp, bin_r_dry, bin_num_conc
        real(kind=dp), intent(in) :: wv, wv_predict, wv_sat, wv_sat_predict, wl, wl_predict, temp, temp_predict
        real(kind=dp) :: radius_change_max, dwl, radius_change_predict
        real(kind=dp), intent(out) :: t_step_cond_ratio_radius, t_step_cond_ratio_temperature, t_step_cond_ratio_moisture
        integer, intent(out) :: ind_bin_radius_control
        integer, intent(in) :: num_cond_bins
        integer :: ind_bin
        ! Initialize cpondensation time step
        t_step_cond_ratio_radius = 1.0
        t_step_cond_ratio_temperature = 1.0
        t_step_cond_ratio_moisture = 1.0
        ind_bin_radius_control = -1
                            
        if (if_limit_t_step_cond_by_radius) then
        ! radias limitation for time steps     
            do ind_bin = 1, num_cond_bins
                if (bin_r_wet_temp(ind_bin)*(1+tolerance_r_ratio) < bin_r_dry(ind_bin)) then
                    print *, "ERROR: wet radius < dry radius" 
                    print *, "r_wet", bin_r_wet_temp(ind_bin), "r_dry", bin_r_dry(ind_bin)
                    call stopcm1
                end if
                radius_change_predict = bin_r_wet_predict(ind_bin) - bin_r_wet_temp(ind_bin)
                if (abs(radius_change_predict ) < tolerance_bin_r*1) then
                    cycle
                else if (radius_change_predict  < 0) then
                    if (if_limit_t_step_cond_by_bin) then
                        radius_change_max = bin_r_wet_temp(ind_bin) - max(bin_r_wet_grid(ind_bin), bin_r_dry(ind_bin))
                    else
                        radius_change_max = bin_r_wet_temp(ind_bin) - bin_r_dry(ind_bin)
                    end if
                    if (radius_change_max*percent_radius_change_maximum < -radius_change_predict) then
                        if (abs(radius_change_max*percent_radius_change_maximum/radius_change_predict) < t_step_cond_ratio_radius) then
                            t_step_cond_ratio_radius = abs(radius_change_max*percent_radius_change_maximum/radius_change_predict)
                            ind_bin_radius_control = ind_bin
                        end if
                    end if

                else if (radius_change_predict > 0) then
                    if (if_limit_t_step_cond_by_bin) then
                        radius_change_max = min(bin_r_wet_grid(ind_bin+1)-bin_r_wet_temp(ind_bin), bin_r_wet_temp(ind_bin))
                    else
                        radius_change_max = bin_r_wet_temp(ind_bin)
                        ! radius_change_max = bin_r_wet_temp(num_bins) - bin_r_wet_temp(ind_bin)
                    end if
                    if (radius_change_max*percent_radius_change_maximum < radius_change_predict) then
                        if (abs(radius_change_max*percent_radius_change_maximum/radius_change_predict) < t_step_cond_ratio_radius) then
                            t_step_cond_ratio_radius = abs(radius_change_max*percent_radius_change_maximum/radius_change_predict)
                            ind_bin_radius_control = ind_bin
                        end if
                    end if
                end if
                        
            end do
        end if

        dwl = wl_predict - wl
        if (if_limit_t_step_cond_by_temp) then
        ! temperature limitation
            if (dwl /=0) then
                t_step_cond_ratio_temperature = min(temperature_change_maximum/abs(temp_predict-temp), t_step_cond_ratio_temperature)
            end if
        end if

        if (if_limit_t_step_cond_by_moisture .and. dwl < -tolerance_wl) then   
            ! liquid water limitation
            if (-dwl > percent_wl_change_maximum*wl) then
                t_step_cond_ratio_moisture = min((wl*percent_wl_change_maximum)/(-dwl), t_step_cond_ratio_moisture)
            end if
        end if

        if (if_limit_t_step_cond_by_supersaturation) then
            if ((dwl>tolerance_wl .and. wv_predict+tolerance_wv < wv_sat_predict) .or. (dwl<-tolerance_wl .and. wv_predict > wv_sat_predict+tolerance_wv)) then
                t_step_cond_ratio_moisture = max(0.1d0,(dwl-(wv_sat_predict-wv_predict))/dwl)
            end if
        end if

    end subroutine determine_t_step_cond_ratio

    subroutine update_bins_static(bin_r_dry_curr, bin_r_wet_curr, bin_num_conc_curr, note)
        ! subroutine to change bin dry & wet radius and bin number concentration if using static bins
        real (kind=dp), dimension(:), intent(inout) :: bin_r_dry_curr, bin_r_wet_curr, bin_num_conc_curr
        real (kind=dp), dimension(num_bins+1) :: bin_r_wet_curr_grid
        real (kind=dp), dimension(num_bins) :: bin_num_conc_temp, bin_r_dry_temp, bin_vol_wet_temp, bin_vol_dry_temp
        real (kind=dp):: ratio_mean, radius_overlap, distribute_ratio, distribute_num_conc, vol_bin, distribute_ratio_total, r_wet_curr, num_conc_curr
        integer :: ind_bin, ind_bin_init, ind_bin_shift
        real (kind=dp) :: locator_bin
        integer :: ind_left, ind_right
        real (kind=dp) :: dn_left, dn_right
        real (kind=dp) :: ratio_change_num_tot_stat
        real (kind=dp) :: old_sum_wet, new_sum_wet, old_sum_dry, new_sum_dry, old_sum_conc, new_sum_conc
        character(len=*), intent(in) :: note
        bin_num_conc_temp = 0
        bin_vol_dry_temp = 0
        do ind_bin = 1, num_bins
            if (bin_r_dry_curr(ind_bin) > bin_r_wet_curr(ind_bin)*(1+tolerance_r_ratio)) then
                print *, "ERROR: dry radius larger than wet radius before static bin update"
                print *, "bin ", ind_bin, "dry radius:", bin_r_dry_curr(ind_bin), "wet radius:", bin_r_wet_curr(ind_bin)
                call stopcm1
            end if
            bin_r_dry_curr(ind_bin) = min(bin_r_dry_curr(ind_bin),bin_r_wet_curr(ind_bin))
        end do
        do ind_bin = 1, num_bins
            r_wet_curr = bin_r_wet_curr(ind_bin)
            num_conc_curr = bin_num_conc_curr(ind_bin)
            if (r_wet_curr < bin_r_wet(1)) then
                if (if_warning_general) then
                    print *, "WARNING: The bin is smaller than the first bin", ind_bin, note
                    print *, "r_wet_curr, r_dry_curr, num_conc_curr: ", r_wet_curr, bin_r_dry_curr(ind_bin), num_conc_curr
                    print *, "Modified to: ", r_wet_curr, bin_r_dry_curr(ind_bin), num_conc_curr
                end if
                bin_num_conc_curr(ind_bin) = num_conc_curr * (r_wet_curr**3d0-bin_r_dry_curr(ind_bin)**3d0)/(bin_r_wet(1)**3d0-bin_r_dry_curr(ind_bin)**3d0)
                bin_r_wet_curr(ind_bin) = bin_r_wet(1)
                
            end if
        end do

        do ind_bin = 1, num_bins
            num_conc_curr = bin_num_conc_curr(ind_bin)
            r_wet_curr = bin_r_wet_curr(ind_bin)
            if (num_conc_curr < tolerance_bin_num_conc) then
                cycle
            end if
            
            if (r_wet_curr > bin_r_wet(num_bins)) then
                print *, "ERROR: The bin is larger than the last bin", locator_bin, ">", num_bins, r_wet_curr, ">", bin_r_wet(num_bins), bin_r_dry_curr(ind_bin), note
                print *, ind_bin, "number conc:", num_conc_curr
                call stopcm1
            end if
            
            locator_bin = (log(r_wet_curr) - log(bin_r_wet(1)))/log(ratio_bins) + 1
            ind_left = max(1,floor(locator_bin))
            ind_right = min(num_bins, ceiling(locator_bin))
            
            if (r_wet_curr < bin_r_wet(ind_left)) then
                do ind_bin_shift = 1, ind_left-1
                    if (bin_r_wet(ind_left-ind_bin_shift) <= r_wet_curr) then
                        ind_left = ind_left-ind_bin_shift
                        ind_right = ind_left + 1
                        exit
                    end if
                end do
            end if
            if (r_wet_curr > bin_r_wet(ind_right)) then
                do ind_bin_shift = 1, num_bins-ind_right
                    if (bin_r_wet(ind_right+ind_bin_shift) >= r_wet_curr) then
                        ind_right = ind_right+ind_bin_shift
                        ind_left = ind_right - 1
                        exit
                    end if
                end do
            end if
            if (r_wet_curr < bin_r_wet(ind_left) .or. r_wet_curr > bin_r_wet(ind_right)) then
                print *, "ERROR in finding bins to distribute mass of radius ", r_wet_curr, note
                print *, "left bin index:", ind_left, "left radius:", bin_r_wet(ind_left)
                print *, "right bin index:", ind_right, "right radius:", bin_r_wet(ind_right)
                call stopcm1
            end if
            if (ind_left == ind_right) then
                dn_left = num_conc_curr
                dn_right = 0
            else
                dn_left = (bin_r_wet(ind_right)-r_wet_curr)*(r_wet_curr**2+bin_r_wet(ind_right)**2+bin_r_wet(ind_right)*r_wet_curr)/(bin_r_wet(ind_right)**3-bin_r_wet(ind_left)**3)*num_conc_curr
                dn_right = (r_wet_curr-bin_r_wet(ind_left))*(r_wet_curr**2+bin_r_wet(ind_left)**2+r_wet_curr*bin_r_wet(ind_left))/(bin_r_wet(ind_right)**3-bin_r_wet(ind_left)**3)*num_conc_curr     
                if (dn_left < -tolerance_bin_num_conc .or. dn_right < -tolerance_bin_num_conc) then
                    if (abs(min(dn_left, dn_right))/abs(dn_left+dn_right)>1e-6) then
                        print *, "ERROR: distributed negative number concentraton ", note
                        print *, dn_left, dn_right
                        print *, "left:",bin_r_wet(ind_left), "right:", bin_r_wet(ind_right), "radius:", r_wet_curr
                        call stopcm1
                    end if
                end if
                dn_left = max(0d0,dn_left)
                dn_right = max(0d0, dn_right)
            end if         
            bin_num_conc_temp(ind_left) = bin_num_conc_temp(ind_left) + dn_left
            bin_num_conc_temp(ind_right) = bin_num_conc_temp(ind_right) + dn_right
            bin_vol_dry_temp(ind_left) = bin_vol_dry_temp(ind_left) + dn_left*(bin_r_dry_curr(ind_bin)*bin_r_wet(ind_left)/r_wet_curr)**3
            bin_vol_dry_temp(ind_right) = bin_vol_dry_temp(ind_right) + dn_right*(bin_r_dry_curr(ind_bin)*bin_r_wet(ind_right)/r_wet_curr)**3
        end do

        do ind_bin = 1, num_bins
            if (bin_num_conc_temp(ind_bin)>tolerance_bin_num_conc) then
                bin_r_dry_temp(ind_bin) = (bin_vol_dry_temp(ind_bin)/bin_num_conc_temp(ind_bin))**(1d0/3d0)
            else
                bin_r_dry_temp(ind_bin) = bin_r_dry_curr(ind_bin)
            end if
        end do

        if ((minval(bin_r_dry_temp) < -tolerance_bin_r) .or. (minval(bin_num_conc_temp) < -tolerance_bin_num_conc) ) then
            print *, "negative dry radius or number concentration for bin in redistribution", note
            do ind_bin = 1, num_bins
                print *, ind_bin, bin_r_dry_curr(ind_bin), bin_r_wet(ind_bin), bin_num_conc_temp(ind_bin)
            end do
            call stopcm1
        end if

        do ind_bin = 1, num_bins
            if (bin_r_dry_temp(ind_bin) > bin_r_wet(ind_bin)*(1+tolerance_r_ratio) .and. .true.) then
                print *, "WARNING: dry radius larger than wet radius after static bin update", note
                print *, "bin ", ind_bin, "dry radius:", bin_r_dry_temp(ind_bin), "wet radius:", bin_r_wet(ind_bin), "number concentration:", bin_num_conc_temp(ind_bin)
            end if
            bin_r_dry_temp(ind_bin) = min(bin_r_dry_temp(ind_bin), bin_r_wet(ind_bin))
            
        end do

        old_sum_conc = sum(bin_num_conc_curr)
        new_sum_conc = sum(bin_num_conc_temp)
        if (abs(new_sum_conc-old_sum_conc)/abs(old_sum_conc) > 0.01) then
            print *, "myid, myi, myj", myid, myi, myj
            print *, "ERROR: the number concentration is not balanced upon redistributing mass ", note
            print *, "ind_bin, bin_num_conc_curr(ind_bin), bin_num_conc_temp(ind_bin), bin_r_wet_curr, bin_r_wet"
            do ind_bin = 1, num_bins
                print *, ind_bin, bin_num_conc_curr(ind_bin), bin_num_conc_temp(ind_bin), bin_r_wet_curr(ind_bin), bin_r_wet(ind_bin)
            end do
            print *, "sum before:", old_sum_conc
            print *, "sum after:", new_sum_conc
            call stopcm1
        end if

        old_sum_wet = sum(bin_r_wet_curr**3*bin_num_conc_curr)*4*pi/3*1e3/1
        new_sum_wet = sum(bin_r_wet**3*bin_num_conc_temp)*4*pi/3*1e3/1
        if (abs(new_sum_wet - old_sum_wet)> max(0.01*abs(old_sum_wet), tolerance_bin_wl)) then
            print *, "myid, myi, myj", myid, myi, myj
            print *, "ERROR: the mass is not balanced for wet radius upon redistributing mass ", note
            print *, "ind_bin, bin_num_conc_curr(ind_bin), bin_num_conc_temp(ind_bin), bin_r_wet_curr, bin_r_wet"
            do ind_bin = 1, num_bins
                print *, ind_bin, bin_num_conc_curr(ind_bin), bin_num_conc_temp(ind_bin), bin_r_wet_curr(ind_bin), bin_r_wet(ind_bin)
            end do
            print *, "sum before:", old_sum_wet
            print *, "sum after:", new_sum_wet
            call stopcm1
        end if
        ! bin_num_conc_temp = bin_num_conc_temp * (old_sum_wet/new_sum_wet)

        old_sum_dry = sum(bin_r_dry_curr**3*bin_num_conc_curr) *4*pi/3*1e3/1
        new_sum_dry = sum(bin_r_dry_temp**3*bin_num_conc_temp) *4*pi/3*1e3/1
        if ((abs(new_sum_dry - old_sum_dry)> max(0.01*abs(old_sum_dry), tolerance_wl*0.001)).and.(if_warning_tolerance_aerosol .eqv. .true.)) then
            print *, "WARNING: the mass is not balanced for dry radius upon redistributing mass ", note
            print *, "before:", old_sum_dry
            print *, "after:", new_sum_dry
            ! print *, bin_r_dry_temp
            ! call stopcm1
        end if
        ! bin_r_dry_temp = bin_r_dry_temp * (old_sum_dry/new_sum_dry) ** (1d0/3d0) ! This code causes instability



        bin_num_conc_curr = bin_num_conc_temp
        bin_r_wet_curr = bin_r_wet
        bin_r_dry_curr = bin_r_dry_temp

    end subroutine update_bins_static

    function calculate_equilibrium_dwl(wv, temp, p) result(dwl)
        real (kind=dp), intent(in) :: wv, temp, P
        real (kind=dp) :: wv_sat, dwl, dwl_left, dwl_right, dwl_mid, error_left, error_right, error_mid
        integer :: num_iter
        integer :: max_num_iter = 10000
        integer :: warning_max_num_iter = 1000
        wv_sat = calculate_wv (temp, 0d0, p)
        if (wv <= wv_sat ) then
            dwl = 0
        else
            ! bisection method
            dwl_left = 0
            dwl_right = wv-wv_sat
            num_iter = 0

            do while (.true.)
                num_iter = num_iter + 1
                if (num_iter > max_num_iter) then
                    print *, "ERROR: reached maximum iteration time!", " maximum iteration time is ", max_num_iter
                    print *, "temp:", temp, "wv:", wv, "p:", p, "dwl_Left:", dwl_left, "dwl_right:", dwl_right, "dwl_mid:", dwl_mid, "error_left:", error_left, "error_right:", error_right, "error_mid:", error_mid
                    call stopcm1
                end if

                dwl_mid = 0.5d0*(dwl_left+dwl_right)
                error_left = wv - dwl_left - calculate_wv (temp + L*dwl_left/C_p, 0d0, p)
                error_right = wv - dwl_right - calculate_wv (temp + L*dwl_right/C_p, 0d0, p)
                error_mid = wv - dwl_mid - calculate_wv (temp + L*dwl_mid/C_p, 0d0, p)


                if (abs(error_mid)  < tolerance_wv) then
                    if (num_iter > warning_max_num_iter) then
                        print *, "Warning: the iteration of bisection method to find equilibrium wl is too large:", num_iter
                    end if
                    dwl = dwl_mid
                    exit
                end if

                if (error_left*error_right > 0) then
                    print *, "ERROR: same sign for dwl_left and dwl_right. Bisection Failed."
                    print *, "temp:", temp, "wv:", wv, "p:", p, "dwl_Left:", dwl_left, "dwl_right:", dwl_right, "dwl_mid:", dwl_mid, "error_left:", error_left, "error_right:", error_right, "error_mid:", error_mid
                    call stopcm1
                else 
                    if (error_mid*error_left < 0) then
                        dwl_right = dwl_mid
                    else if (error_mid*error_right < 0) then
                        dwl_left = dwl_mid
                    else
                        ! beyond computation limit
                        dwl = dwl_mid
                        exit
                    end if
                end if
                if (abs(dwl_left-dwl_right)<tolerance_wv) then
                    exit
                end if
            end do
        end if

    end function calculate_equilibrium_dwl


!#####################           Coalescence & Breakup            ###################

    function calculate_terminal_velocity_Best (r_drop) result(v_term)
        ! Best falling velocities (Bs50)
        real (kind=dp), intent(in) :: r_drop
        real(kind=dp) :: v_term, d1
        d1=r_drop*2d0
        v_term = 943d0*(1d0-exp(-(d1*1000d0/1.77d0)**(1.147d0)))/100d0 ! unit[m/s]
    end function calculate_terminal_velocity_Best

    function calculate_terminal_velocity_Lhermitte (r_drop, rho0rhoz) result(v_term)
        ! Lhermitte (1990) falling velocities
        real (kind=dp), intent(in) :: r_drop
        real (kind=sp), intent(in) :: rho0rhoz ! surface density/z-level density
        real(kind=dp) :: v_term, d1
        d1=r_drop*2d0
        v_term = 923d0*(rho0rhoz)**0.5d0*(1d0-exp(-6.8d0*(d1*1000d0/10d0)**2d0-4.88d0*(d1/10d0)))/100d0 ! unit[m/s]
    end function calculate_terminal_velocity_Lhermitte

    function calculate_weight(vols, i, j, k) result(weight)
        real (kind=dp), dimension(:), intent(in) :: vols
        integer, intent(in) :: i, j, k ! the weight is the volume contribution for bin ind i from (bin ind j and bin ind k)        
        integer num_coal_bins ! size of coalescence bins (not the same with total bins)
        real (kind=dp) :: vol, weight
        weight = ieee_value(weight, ieee_quiet_nan) !Initialize with NaN
        num_coal_bins = size(vols)

        vol = vols(j) + vols(k)
        if (i == 1) then
            if ((vol >= vols(i)) .and. (vol <= vols(i+1))) then
                weight = (vols(i+1)-vol)/(vols(i+1)-vols(i))
            else if (vol < vols(i)) then
                weight = vol/vols(i)
            else
                weight = 0d0
            end if
        end if
        if ((i > 1) .and. (i < num_coal_bins)) then
            if ((vol >= vols(i-1)) .and. (vol <= vols(i))) then
                weight = (vol-vols(i-1))/(vols(i)-vols(i-1))
            else if((vol >= vols(i)) .and. (vol <= vols(i+1))) then
                weight = (vols(i+1)-vol)/(vols(i+1)-vols(i))
            else
                weight = 0d0
            end if
        end if
        if(i == num_coal_bins) then
            if((vol >= vols(i-1)).and.(vol <= vols(i))) then
                weight = (vol-vols(i-1))/(vols(i)-vols(i-1))
            else if (vol > vols(i)) then
                ! weight = 1d0
                weight = vol/vols(i)
            else
                weight = 0d0
            end if
        end if
        if (ieee_is_nan(weight)) then
            print *, "ERROR: The weight is NaN!"
            call stopcm1
        end if
    end function

    subroutine calculate_weights(weights, bin_coal_vol)
        real (kind=dp), dimension(:,:,:), intent(out) :: weights
        real (kind=dp), dimension(:), intent(in) :: bin_coal_vol
        integer num_coal_bins ! size of coalescence bins (not the same with total bins)
        integer i, j, k ! the weight is the volume contribution for bin ind i from (bin ind j and bin ind k)
        num_coal_bins = size(bin_coal_vol)

        do i = 1, num_coal_bins
            do j = 1, num_coal_bins
                do k = 1, num_coal_bins               
                    weights(i,j,k) = calculate_weight(bin_coal_vol, i, j, k)
                end do
            end do
        end do
        do j = 1, num_coal_bins
            do k = 1, num_coal_bins
                if (abs( sum( weights(1:num_coal_bins,j,k)*bin_coal_vol(1:num_coal_bins) ) / (bin_coal_vol(j)+bin_coal_vol(k)) -1) > 1e-6) then
                    print *, "ERROR in weight calcualtion:", abs(sum(weights(1:num_coal_bins,j,k)*bin_coal_vol(1:num_coal_bins))/(bin_coal_vol(j)+bin_coal_vol(k))-1)
                    print *, j, k
                    call stopcm1
                end if
            end do
        end do

    end subroutine calculate_weights

    subroutine calculate_kernels_coal (kernels_coal_curr, bin_coal_r_wet, bin_coal_v_term, num_coal_bins, num_cond_coal_bins, z)
        integer, intent(in) :: num_coal_bins, num_cond_coal_bins ! size of coalescence bins (not the same with total bins)    
        real(kind=dp), dimension(num_coal_bins, num_coal_bins), intent(out) :: kernels_coal_curr
        real(kind=dp), dimension(num_coal_bins), intent(in) :: bin_coal_r_wet, bin_coal_v_term
        real(kind=dp), intent(in) :: z
        real(kind=dp) :: eff_coll, kernel_grav_coll, eff_coal, r_j, r_k, v_term_j, v_term_k
        integer:: j, k
        real (kind=dp), dimension(size(rd1),size(rd1)):: EF1
        real (kind=dp), dimension(size(rd1),size(rd2)):: EF2
        real (kind=dp), dimension(size(rd2),size(rd2)):: EF3
        call interpolate_collision_kernel (z, EF1, EF2, EF3)
        ! initialize coalescence kernel
        kernels_coal_curr = 0d0
        do j=1, num_cond_coal_bins
            do k=1, num_cond_coal_bins
                r_j = bin_coal_r_wet(j)
                r_k = bin_coal_r_wet(k)
                v_term_j = bin_coal_v_term(j)
                v_term_k = bin_coal_v_term(k)
                eff_coll = calculate_eff_coll (r_j, r_k, EF1, EF2, EF3)
                kernel_grav_coll = pi*(r_j+r_k)**2d0*abs(v_term_j-v_term_k) * eff_coll
                eff_coal = calculate_eff_coal (r_j, r_k, v_term_j, v_term_k)
                ! print *, j, k, eff_coll, kernel_grav_coll
                kernels_coal_curr(j,k) = kernel_grav_coll * eff_coal
                ! print *, 'small r: ',min(r_j,r_k), 'large r: ',max(r_j,r_k), 'eff_coll:', eff_coll
            end do
        end do
        
    end subroutine calculate_kernels_coal




    function calculate_eff_coal (r_drop_1, r_drop_2, v_term_1, v_term_2) result(eff_coal)
        real(kind=dp), intent(in) :: r_drop_1, r_drop_2, v_term_1, v_term_2
        real(kind=dp) :: r_l, r_s, v_term_l, v_term_s, d_l, d_s
        real(kind=dp) :: beta, bb, aa, xx, cke, sc, st, et, eff_coal_ll82, wenbr, qbo95, eff_coal_b95a, eff_coal_b95b, eff_coal_bo95, bo1, bo2
        real(kind=dp) :: eff_coal
        integer :: ind_iter
        eff_coal = ieee_value(eff_coal, ieee_quiet_nan) !Initialize with NaN
        ! get radius, velocity and diameter for large and small drops
        if (r_drop_1 > r_drop_2) then
            r_l = r_drop_1
            r_s = r_drop_2
            v_term_l = v_term_1
            v_term_s = v_term_2
        else
            r_l = r_drop_2
            r_s = r_drop_1
            v_term_l = v_term_2
            v_term_s = v_term_1
        end if 
        d_l = r_l*2
        d_s = r_s*2

        if(r_s > 14d-6) then
            ! ieffi = 1 => constant
            if (ieff_coal == 1) then
                eff_coal = const_eff_coal
            ! ieffi = 2 => E=(BO1984)
            else if (ieff_coal == 2) then
                if((r_s >= 1d-6).and.(r_s <= 30d-6).and.(r_l >= 50d-6).and.(r_l <= 500d-6)) then
                    beta = log(r_s*1d6)+0.44d0*log(r_l/200d0)
                    bb = 0.0946d0*beta-0.319d0
                    aa = sqrt(bb**2d0+0.00441d0)
                    xx = (aa-bb)**(1d0/3d0)-(aa+bb)**(1d0/3d0)
                    eff_coal=xx+0.459d0
                else
                    eff_coal=0d0
                end if
            ! ieffi = 3 => E=a[1+d_s/d_l]^(-2).exp[-b.sig.Et^2/Sc](LL82a)
            else if (ieff_coal == 3) then
                cke=(1000d0*pi/12d0)*(d_l*d_s)**3d0/ (d_l**3d0+d_s**3d0)*(v_term_l-v_term_s)**2d0
                sc=7.28d-2*pi*(d_l**3d0+d_s**3d0)**(2d0/3d0)
                st=7.28d-2*pi*(d_l**2d0+d_s**2d0)
                et=cke+(st-sc)
                if(et < effi_etlmt) then
                    eff_coal = 0.778d0*(1d0+d_s/d_l)**(-2d0)*exp(-2.61d6*7.28d-2*et**2d0/sc)
                else
                    eff_coal = 0d0
                end if
            ! ieffi = 4 => E=f(LL82a,BO95) (Seifert at al. 2005)
            else if(ieff_coal == 4) then
                ! Coalescence efficiency from LL82a : eff_coal_ll82        
                cke=(1000d0*pi/12d0)*(d_l*d_s)**3d0/ (d_l**3d0+d_s**3d0)*(v_term_l-v_term_s)**2d0
                sc=7.28d-2*pi*(d_l**3d0+d_s**3d0)**(2d0/3d0)
                st=7.28d-2*pi*(d_l**2d0+d_s**2d0)
                et=cke+(st-sc)
                if(et < effi_etlmt) then
                    eff_coal_ll82=0.778d0*(1d0+d_s/d_l)**(-2d0)* exp(-2.61d6*7.28d-2*et**2d0/sc)
                else
                    eff_coal_ll82=0d0
                end if

                ! Coalescence efficiency from BO95 : eff_coal_bo95
                wenbr=1000d0*r_s*(v_term_l-v_term_s)**2d0/7.28d-2
                qbo95=d_s/d_l
                ! Coalescence efficiency from BO95 : eff_coal_b95a
                eff_coal_b95a=0.767d0-10.14d0*(2d0**(3d0/2d0)/(6d0*pi))*sqrt(wenbr)*(qbo95**4d0*(1d0+qbo95)/((1d0+qbo95**2d0)*(1d0+qbo95**3d0)))
                if(eff_coal_b95a.lt.0d0) then
                    eff_coal_b95a= 0d0
                end if
                ! Coalescence efficiency from BO95 : eff_coal_b95b
                bo1=1d0
                do ind_iter = 1, 10
                    bo2=bo1
                    bo1=bo2-(-5.29d0*bo2**3d0+7.27d0*bo2**2d0-5.94d0*bo2+5.07d0 -log(r_s*1d6)-log(r_l*1d6/200d0))/(3d0*(-5.29d0)*bo2**2d0+2d0*7.27d0*bo2-5.94d0)
                    if(abs((bo2-bo1)/bo2).le.1d-6) then
                        exit
                    end if
                end do
                if(bo1.lt.0.) then
                    bo1=0d0
                end if
                eff_coal_b95b=bo1
                ! Coalescence efficiency from BO95 : eff_coal_bo95=max(eff_coal_b95a,eff_coal_b95b)
                eff_coal_bo95=max(eff_coal_b95a,eff_coal_b95b)

                !	Coalescence efficiency combined : eff_coal
                if(d_s < 3d-4) then
                    eff_coal=eff_coal_bo95
                else if (d_s > 6d-4) then
                    eff_coal=eff_coal_ll82
                else if ((d_s >= 3d-4).and.(d_s <= 6d-4)) then
                    eff_coal=cos(pi/2d0*(d_s-3d-4)/3d-4)**2d0*eff_coal_bo95+sin(pi/2d0*(d_s-3d-4)/3d-4)**2d0*eff_coal_ll82
                end if
            else
                print *, "Error: coalescence efficiency option not specified"
                call stopcm1
            end if
        else if(r_s <= 0d0) then
                eff_coal=0d0
        else if((r_s <= 14d-6).and.(r_s > 0d0)) then
                eff_coal=1d0
        end if
        if(eff_coal > 1d0) then 
            eff_coal=1d0
        end if 
        if(eff_coal <= 0d0) then
            eff_coal=0d0
        end if
        if (ieee_is_nan(eff_coal)) then
            print *, "ERROR: The coalescence efficiency is NaN!"
            call stopcm1
        end if
        ! print *, 'r_drop_1: ', r_drop_1, 'r_drop_2: ', r_drop_2, 'eff_coal: ', eff_coal
    end function calculate_eff_coal

    function calculate_eff_coll (r1, r2, EF1, EF2, EF3) result(eff_coll)
        real (kind=dp), intent(in) :: r1, r2
        real (kind=dp), dimension(:,:), intent(in) :: EF1, EF2, EF3        
        real (kind=dp) :: eff_coll

        if ((r1 > 1d-6) .and. (r2 > 1d-6)) then
            eff_coll = precalculate_eff_coll (r1, r2, EF1, EF2, EF3)
        else
            eff_coll = 0d0
        end if
        if (ieee_is_nan(eff_coll)) then
            print *, "ERROR: nan eff_coll."
            call stopcm1
        end if
        if (eff_coll < 0) then
            eff_coll = 0d0
        end if
    end function calculate_eff_coll

    function kronecker_delta (j, k) result(delta_jk)
        integer, intent(in) :: j, k
        integer :: delta_jk
        if (j == k) then
            delta_jk = 1
        else
            delta_jk = 0
        end if
    end function kronecker_delta

    function precalculate_eff_coll(R1, R2, EF1, EF2, EF3) result(COLK)
        !!----(PINKSY AND KHAIN, 2002)---------------------------------------
        ! Interpolate data/EFF to each drop pair.
        !
        !    Parameters
        !    ----------
        !    Zh: Parcel Height [m]
        !    R1: Droplet radius [m]
        !    R2: Droplet radius [m]
        !    if_turb: Turn on/off turbulence (1/0) defined in parameters_and_properties.f90
        !    Returns
        !    ----------
        !    COLK: Ecoll for different drop pairs
        !!-------------------------------------------

        real (kind=dp), intent(in) :: R1, R2
        real (kind=dp), dimension(:,:), intent(in) :: EF1, EF2, EF3
        real (kind=dp) :: COLK
        real (kind=dp) :: EE, E1, E2, RR1, RR2, FCT, T1, T2, TT
        integer :: I, J, II, JJ, ND1, ND2
        
        II = 0
        JJ = 0
        ND1=size(RD1)
        ND2=size(RD2)
        COLK = ieee_value(COLK, ieee_quiet_nan)
        !----------------------------------------------
        ! For big drops, EE =1
        ! r1> 250 micron or r2> 250 micron
        !----------------------------------------------
        IF(R1 >= RD2(ND2).or.R2 >= RD2(ND2)) THEN
            EE = 1d0
            COLK = EE
            RETURN
        END IF

        RR1 = min(R1, R2)
        RR2 = max(R1, R2)
        !----------------------------------------------
        ! INTERPOLATION OF THE 1st TABLE,
        ! r1 < 40 micron & r2 < 40 micron
        !----------------------------------------------

        IF(R1 <= RD1(ND1).and.R2 <= RD1(ND1)) THEN

            DO I=1,ND1
                IF(R1 <= RD1(I)) THEN
                    II = I
                    exit
                END IF
            END DO

            DO J=1,ND1
                IF(R2 <= RD1(J)) THEN
                    JJ = J
                    exit
                END IF
            END DO

            ! Using Bilinear interpolation
            E1 = (EF1(II,JJ-1)-EF1(II-1,JJ-1))/(RD1(II)-RD1(II-1))
            E1 = EF1(II-1,JJ-1) + E1*(R1-RD1(II-1))
            E2 = (EF1(II,JJ)-EF1(II-1,JJ))/(RD1(II)-RD1(II-1))
            E2 = EF1(II-1,JJ) + E2*(R1-RD1(II-1))
            EE = (E2-E1)/(RD1(JJ)-RD1(JJ-1))
            EE = E1 + EE*(R2-RD1(JJ-1))

            COLK = EE
            if (ieee_is_nan(COLK)) then
                print *, "2323 ERROR: The collision efficiency is NaN!", RD1(JJ), RD1(JJ-1), E2, E1,RD1(II),RD1(II-1),EF1(II,JJ),EF1(II-1,JJ),R1
                call stopcm1
            end if
            !Pinsky, M., Khain, A., and Krugliak, H.: Collisions of Cloud Droplets in a Turbulent Flow.
            !Part V: Application of Detailed Tables of Turbulent Collision Rate Enhancement to Simulation
            !of Droplet Spectra Evolution, Journal of the Atmospheric Sciences, 65, 357-374, 2008.
            IF((if_turb).and.(RR2 <= 21d-6)) THEN

                T1 = (turbulence_kernel(II,JJ-1)-turbulence_kernel(II-1,JJ-1))/(RD1(II)-RD1(II-1))
                T1 = turbulence_kernel(II-1,JJ-1) + T1*(R1-RD1(II-1))
                T2 = (turbulence_kernel(II,JJ)-turbulence_kernel(II-1,JJ))/(RD1(II)-RD1(II-1))
                T2 = turbulence_kernel(II-1,JJ) + T2*(R1-RD1(II-1))
                TT = (T2-T1)/(RD1(JJ)-RD1(JJ-1))
                TT = T1 + TT*(R2-RD1(JJ-1))
                FCT = TT

                COLK = FCT*EE

            END IF
            if (ieee_is_nan(COLK)) then
                print *, "2340 ERROR: The collision efficiency is NaN!"
                call stopcm1
            end if
            RETURN
        END IF

        !----------------------------------------------
        ! INTERPOLATION OF THE 2nd TABLE,
        ! min(R1, R2) < 40 micron & max(R1, R2) > 50 micron
        !----------------------------------------------
        IF((RR1.lt.RD1(ND1)).and.(RR2.gt.RD2(1))) THEN

            DO I=1,ND1
                IF(RR1 <= RD1(I)) THEN
                    II = I
                    exit
                END IF
            END DO

            DO J=1,ND2
                IF(RR2 <= RD2(J)) THEN
                    JJ = J
                    exit
                END IF
            END DO

            E1 = (EF2(II,JJ-1)-EF2(II-1,JJ-1))/(RD1(II)-RD1(II-1))
            E1 = EF2(II-1,JJ-1) + E1*(RR1-RD1(II-1))
            E2 = (EF2(II,JJ)-EF2(II-1,JJ))/(RD1(II)-RD1(II-1))
            E2 = EF2(II-1,JJ) + E2*(RR1-RD1(II-1))
            EE = (E2-E1)/(RD2(JJ)-RD2(JJ-1))
            EE = E1 + EE*(RR2-RD2(JJ-1))

            COLK = EE
            if (ieee_is_nan(COLK)) then
                print *, "2371 ERROR: The collision efficiency is NaN!"
                call stopcm1
            end if
            RETURN

        END IF

        !----------------------------------------------
        ! INTERPOLATION OF THE 3rd TABLE,
        ! r1 > 50 micron & r2 > 50 micron
        !----------------------------------------------
        IF((R1 >= RD2(1)).and.(R2 >= RD2(1))) THEN

            DO I=1,ND2
                IF(R1 <= RD2(I)) THEN
                    II = I
                    exit
                END IF
            END DO

            DO J=1,ND2
                IF(R2 <= RD2(J)) THEN
                    JJ = J
                    exit
                END IF
            END DO

            E1 = (EF3(II,JJ-1)-EF3(II-1,JJ-1))/(RD2(II)-RD2(II-1))
            E1 = EF3(II-1,JJ-1) + E1*(R1-RD2(II-1))
            E2 = (EF3(II,JJ)-EF3(II-1,JJ))/(RD2(II)-RD2(II-1))
            E2 = EF3(II-1,JJ) + E2*(R1-RD2(II-1))
            EE = (E2-E1)/(RD2(JJ)-RD2(JJ-1))
            EE = E1 + EE*(R2-RD2(JJ-1))

            COLK = EE
        if (ieee_is_nan(COLK)) then
            print *, "2403 ERROR: The collision efficiency is NaN!"
            call stopcm1
        end if
            RETURN
        END IF

        if (max(R1, R2) >= maxval(RD2)) then
            COLK = 1
        end if

        if (ieee_is_nan(COLK)) then
            print *, "ERROR: The collision efficiency is NaN!"
            call stopcm1
        end if
        return
    end function precalculate_eff_coll

    subroutine interpolate_collision_kernel (Ze, EF1, EF2, EF3)
        !!----(PINKSY ET AL., 2000)---------------------------------------
        ! Interpolate Ecoll (1000, 750, 500 mb levels) to a given height.
        !
        !    Parameters
        !    ----------
        !    Ze: Parcel Height [m]
        !    rd1: Droplet radius [m]
        !    rd2: Droplet radius [m]
        !    Returns
        !    ----------
        !    EF1: Ecoll for different drop pairs
        !    EF2: Ecoll for different drop pairs
        !    EF3: Ecoll for different drop pairs
        !!-------------------------------------------
        real (kind=dp), dimension(3) :: Z_TAB
        real (kind=dp), intent(in) :: Ze
        real (kind=dp), dimension(:,:), intent(out) :: EF1
        real (kind=dp), dimension(:,:), intent(out) :: EF2
        real (kind=dp), dimension(:,:), intent(out) :: EF3
        real (kind=dp) :: CC
        integer :: nd1, nd2, i, j
        nd1 = size(rd1)
        nd2 = size(rd2)
        !----------------------------------------------
        ! INTERPOLATE TO HEIGHT, Z
        !----------------------------------------------
        Z_TAB = (/0d0, 2500d0, 5000d0/)

        if (Ze.lt.Z_TAB(1)) then
            print *, "ERROR: negative Ze.", Ze
            call stopcm1
        end if

        !-----------------------
        !  Z < 2500 m
        !-----------------------
        IF ((Ze.ge.Z_TAB(1)).and.(Ze <= Z_TAB(2))) THEN
            CC = (Ze-Z_TAB(1))/(Z_TAB(2)-Z_TAB(1))

            DO i=1,nd1
                DO j=1,nd1
                    EF1(i,j) = collision_kernel_0000_1(i,j)+CC*(collision_kernel_2500_1(i,j)-collision_kernel_0000_1(i,j))
                END DO
            END DO
            DO i=1,nd1
                DO j=1,nd2
                    EF2(i,j) = collision_kernel_0000_2(i,j)+CC*(collision_kernel_2500_2(i,j)-collision_kernel_0000_2(i,j))
                END DO
            END DO
            DO i=1,nd2
                DO j=1,nd2
                    EF3(i,j) = collision_kernel_0000_3(i,j)+CC*(collision_kernel_2500_3(i,j)-collision_kernel_0000_3(i,j))
                END DO
            END DO
        RETURN
        END IF

        !-----------------------
        !  2500 m < Z < 5000 m
        !-----------------------
        IF((Ze.gt.Z_TAB(2)).and.(Ze <= Z_TAB(3))) THEN
            CC = (Ze-Z_TAB(2))/(Z_TAB(3)-Z_TAB(2))
            DO i=1,nd1
                DO j=1,nd1
                    EF1(i,j)=collision_kernel_2500_1(i,j)+CC*(collision_kernel_5000_1(i,j)-collision_kernel_2500_1(i,j))
                END DO
            END DO
            DO i=1,nd1
                DO j=1,nd2
                    EF2(i,j)=collision_kernel_2500_2(i,j)+CC*(collision_kernel_5000_2(i,j)-collision_kernel_2500_2(i,j))
                END DO
            END DO
            DO i=1,nd2
                DO j=1,nd2
                    EF3(i,j)=collision_kernel_2500_3(i,j)+CC*(collision_kernel_5000_3(i,j)-collision_kernel_2500_3(i,j))
                END DO
            END DO
        RETURN
        END IF

        !-----------------------
        !  Z > 5000 m
        !-----------------------
        IF(Ze.gt.Z_TAB(3)) THEN
            DO i=1,nd1
                DO j=1,nd1
                    EF1(i,j) = collision_kernel_5000_1(i,j)
                END DO
            END DO
            DO i=1,nd1
                DO j=1,nd2
                    EF2(i,j) = collision_kernel_5000_2(i,j)
                END DO
            END DO
            DO i=1,nd2
                DO j=1,nd2
                    EF3(i,j) = collision_kernel_5000_3(i,j)
                END DO
            END DO
        RETURN
        END IF

        return
    end subroutine interpolate_collision_kernel


subroutine calculate_kernels_coal_brkp (kernels_coal_curr, kernels_brkp_curr, bin_curr_r_wet, bin_curr_v_term, num_curr_bins, z)
    integer, intent(in) :: num_curr_bins ! size of coalescence/breakup bins (not necessarily the same with total bins)    
    real(kind=dp), dimension(num_curr_bins, num_curr_bins), intent(out) :: kernels_coal_curr, kernels_brkp_curr
    real(kind=dp), dimension(num_curr_bins), intent(in) :: bin_curr_r_wet, bin_curr_v_term
    real(kind=dp), intent(in) :: z
    real(kind=dp) :: eff_coll, kernel_grav_coll, eff_coal, eff_brkp, r_j, r_k, v_term_j, v_term_k
    integer:: j, k
    real (kind=dp), dimension(size(rd1),size(rd1)):: EF1
    real (kind=dp), dimension(size(rd1),size(rd2)):: EF2
    real (kind=dp), dimension(size(rd2),size(rd2)):: EF3
    call interpolate_collision_kernel (z, EF1, EF2, EF3)
    ! initialize coalescence and breakup kernels
    kernels_coal_curr = 0d0
    kernels_brkp_curr = 0d0
    do j=1, num_curr_bins
        do k=1, num_curr_bins
            r_j = bin_curr_r_wet(j)
            r_k = bin_curr_r_wet(k)
            v_term_j = bin_curr_v_term(j)
            v_term_k = bin_curr_v_term(k)
            eff_coll = calculate_eff_coll (r_j, r_k, EF1, EF2, EF3)
            kernel_grav_coll = pi*(r_j+r_k)**2d0*abs(v_term_j-v_term_k) * eff_coll
            eff_coal = calculate_eff_coal (r_j, r_k, v_term_j, v_term_k)
            if (min(r_j,r_k) < drop_brkp_parent_min) then
                eff_coal = 1d0
            end if
            eff_brkp = 1d0 - eff_coal
            ! print *, j, k, eff_coll, kernel_grav_coll
            kernels_coal_curr(j,k) = kernel_grav_coll * eff_coal
            kernels_brkp_curr(j,k) = kernel_grav_coll * eff_brkp
            ! print *, 'small r: ',min(r_j,r_k), 'large r: ',max(r_j,r_k), 'eff_coll:', eff_coll
        end do
    end do    
end subroutine calculate_kernels_coal_brkp


subroutine calculate_kernels_coal_brkp_nominal (kernels_coal_curr, kernels_brkp_curr, bin_curr_r_wet, bin_curr_v_term, num_curr_bins)
    ! The kernel with absense of collision efficiency
    integer, intent(in) :: num_curr_bins ! size of coalescence/breakup bins (not necessarily the same with total bins)    
    real(kind=dp), dimension(num_curr_bins, num_curr_bins), intent(out) :: kernels_coal_curr, kernels_brkp_curr
    real(kind=dp), dimension(num_curr_bins), intent(in) :: bin_curr_r_wet, bin_curr_v_term
    real(kind=dp) :: eff_coll, kernel_grav_coll, eff_coal, eff_brkp, r_j, r_k, v_term_j, v_term_k
    integer:: j, k
    ! initialize coalescence and breakup kernels
    kernels_coal_curr = 0d0
    kernels_brkp_curr = 0d0
    do j=1, num_curr_bins
        do k=1, num_curr_bins
            r_j = bin_curr_r_wet(j)
            r_k = bin_curr_r_wet(k)
            v_term_j = bin_curr_v_term(j)
            v_term_k = bin_curr_v_term(k)
            kernel_grav_coll = pi*(r_j+r_k)**2d0*abs(v_term_j-v_term_k)
            eff_coal = calculate_eff_coal (r_j, r_k, v_term_j, v_term_k)
            if (min(r_j,r_k) < drop_brkp_parent_min) then
                eff_coal = 1d0
            end if
            eff_brkp = 1d0 - eff_coal
            ! print *, j, k, eff_coll, kernel_grav_coll
            kernels_coal_curr(j,k) = kernel_grav_coll * eff_coal
            kernels_brkp_curr(j,k) = kernel_grav_coll * eff_brkp
            ! print *, 'small r: ',min(r_j,r_k), 'large r: ',max(r_j,r_k), 'eff_coll:', eff_coll
        end do
    end do    
end subroutine calculate_kernels_coal_brkp_nominal

! kappa is the contribution to droplet population located at the ith interval resulting from the collisional breakup of two drops of volume xj and xk. 
subroutine calculate_kappa_brkp (kappa_brkp, bin_curr_r_wet, bin_curr_r_wet_grid, bin_curr_v_term, num_curr_bins)
    !   NOTE: ONLY IN THIS SUBROUTINE
    !   All dimensions are in cm, volume in cm^3, density in cm^-3
    !   All mass are in grams
    !       => Equivalence mass/volume : rho_water=1gr/cm^3
    integer, intent(in) :: num_curr_bins ! size of breakup bins (not necessarily the same with total bins)    
    real (kind=dp), dimension(1:num_curr_bins, 1:num_curr_bins, 1:num_curr_bins), intent(inout) :: kappa_brkp
    real (kind=dp), dimension(1:num_curr_bins), intent(in) :: bin_curr_r_wet, bin_curr_v_term
    real (kind=dp), dimension(1:num_curr_bins+1), intent(in) :: bin_curr_r_wet_grid
    real (kind=dp), dimension(1:num_curr_bins, 1:num_curr_bins) :: kernel_mass_conservation
    real (kind=dp) :: dl, ds, xl, xs, ml, ms, mcoal, xcoal, dcoal, dcoal_real, vl, vs, cke, sc, st, w1, w2, ds0
    real (kind=dp) :: scdi, scfi, scs2, scsh, scf2, rfi, rsh, rdi, ffi, fsh, fdi, ff1, ff2, fs1, fd1, rcorr, pfi1, pfi2, pfi3, pdi1, pdi2, psh1, psh2, ra2di1, ra2di2, ra2sh1, ra2sh2, ra2sh3, ra2fi1, ra2fi2, ra2fi3
    real (kind=dp) :: sgsh1, sgsh1min, sgsh2, sgsh2max, sgfi1, sgfi2, sgfi3, sgdi1, sgd2, sgdi1min, sgdi1max, sgdi2, sgdi2min, sgf3, sgs2, sgfi3min, sd1i, sd2i, sdi1, sf3i, sfi3
    real (kind=dp) :: massfil,masshee,massdis,cnbrfi1,cmasdi1,cmasdi2, cmasfi1, cmasfi2, cmasfi3, cmassh1, cmassh2, cnbrdi1, cnbrdi2, cnbrfi2, cnbrfi3, cnbrsh1, cnbrsh2
    real (kind=dp) :: ps2i, pfi301, pfi302, pfi303, pf3i, pfi3m, pd2i, mush2, mush1, mus2, mufi3, mufi2, mufi1, muf3, mudi1, mudi2, mudi2min, mudi2max, mud2, ms2i, mf3i, md2i
    real (kind=dp) :: hdi1, hdi2, hd2i, hd2, hf3, hf3i, hfi1, hfi2, hfi2max, hfi3, hs2, hs2i, hsh1, hsh2, ddd1, ddd2, dff3, dss2
    integer :: id1, id2, idi1, iff1, ish1, iexfi1,iexfi2,iexfi3, iexsh1,iexsh2, iexdi1,iexdi2 
    integer :: i, j, k
    real (kind=dp) :: coeff_ms, rate_ms
    real (kind=dp), dimension(num_curr_bins) :: d, x
    real (kind=dp), dimension(num_curr_bins+1) :: d1
    real (kind=dp) :: w0 = 0.86
    real (kind=dp) :: cke_fdlmt = 3.65E-06
    real (kind=dp) :: cke_rflmt = 8.93E-07
    real (kind=dp) :: dcut
    integer :: igaus = 1
    logical :: if_scaling_VY85 = .true.
    d = bin_curr_r_wet*2d0*100d0 ! convert to cm
    d1 = bin_curr_r_wet_grid*2d0*100d0 ! convert to cm
    x = pi/6d0*d**3d0
    dcut = d1(1)

    do j=1, num_curr_bins
        do k=1, num_curr_bins
            if (bin_curr_r_wet(min(k,j)) < drop_brkp_parent_min) then
                cycle
            end if
            if(k.eq.j) then
                cycle ! In this case, gravity collision kernel is 0 (no chance to colide)
            else if (d(j)>d(k)) then
                dl = d(j)
                vl = bin_curr_v_term(j)*100d0
                ds = d(k)
                vs = bin_curr_v_term(k)*100d0
            else
                dl = d(k)
                vl = bin_curr_v_term(k)*100d0
                ds = d(j)
                vs = bin_curr_v_term(j)*100d0
            end if
            
            xl=pi/6.*dl**3.
            xs=pi/6.*ds**3.
            ml = xl * 1d0 ! equivalence mass/volume : rho_water=1gr/cm^3
            ms = xs * 1d0 ! equivalence mass/volume : rho_water=1gr/cm^3
            xcoal = xl+xs ! total volume
            mcoal = ml+ms ! total mass
            dcoal_real=(dl**3.+ds**3.)**(1/3.) ! real coalescence diameter
            dcoal = dcoal_real ! The following numerical treatment is errorness
            ! dcoal = 0 ! numerical bin diameter (largest bin diameter to receive fragments)
            ! do i = 1, num_curr_bins
            !     dcoal = d(i)
            !     if (dcoal > d(i)) then
            !         exit
            !     end if
            ! end do
            ! dcut = d(1) ! This is a modification of the rainshaft model

            ! scaling criteria (VY85)
            scfi=.15*(ml+ms)
            scsh=.15*ml+ms

            cke=(1000.*pi/12d0)*((0.01*dl)*(0.01*ds))**3./((0.01*dl)**3.+(0.01*ds)**3.)*(0.01*(vl-vs))**2.
            sc=7.28e-06*pi*(dl**3.+ds**3.)**(2/3.)
            st=7.28e-06*pi*(dl**2.+ds**2d0)
            w1=cke/sc
            w2=cke/st

            ! rfi : Filament Breakup (LL82)
            if (cke >= cke_rflmt) then
                rfi=1.11e-04*cke**(-0.654)
            else
                rfi=1.
            endif
            ! rsh : Sheet Breakup (LL82)
            if (w2 >= w0) then
                rsh = 0.685*(1.-exp(-1.63*(w2-w0)))
            else 
                rsh = 0
            end if
            ! rsh and rfi adjustment in case rsh + rfi > 1 (List, Donaldson and Stewart, 1987)
            ! and rdi : Disc Breakup (LL82)
            if (rfi+rsh > 1.) then
                rcorr=rfi+rsh
                rfi=rfi/rcorr
                rsh=rsh/rcorr
                rdi=0.
            else
                rdi=1.-(rfi+rsh)
            endif

            ! fdi : Disk Breakup (LL82)
            if(rdi.ne.0.) then
                if(cke < cke_fdlmt) then
                    fd1=2.
                else
                    fd1=297.5+23.7*log(cke)
                endif
                if(fd1 < 2d0) then
                    fdi=2.
                else
                    fdi=fd1
                endif
            else
                fdi=0.
            endif
            ! ffi : Filament Breakup (LL82)
            if(rfi.ne.0.) then
                ff1=(-2.25e+04*(dl-0.403)**2.-37.9)*ds**2.5+9.67*(dl+0.170)**2.+4.95
                ff2=1.02e+04*ds**2.83+2.
                if(ff1 < 2d0) then
                    ff1=2.
                    ds0=((ff1-2d0)/1.02e+04)**(1/2.83)
                else
                    ds0=((ff1-2d0)/1.02e+04)**(1/2.83)
                endif
                if(ds >= ds0) then
                    if(ff1 >= 2d0) then
                        ffi=ff1
                    else
                        ffi=2.
                    endif
                else 
                    ffi=ff2
                endif
            else
                ffi = 0
            end if
            ! fsh : Sheet Breakup (LL82)
            if(rsh.ne.0.) then
                fs1=5.*calculate_errfunct((st-2.53e-06)/1.85e-06)+6.
                if(fs1 < 2d0) then
                    fsh=2.
                else
                    fsh=fs1
                endif
            else
                fsh=0.
            endif
            ! fto : total drop number (LL82), not used in this code
            ! fto = rdi*fdi+rfi*ffi+rsh*fsh

            ! if (rdi .ne. 0.) then
            !     massdis = mcoal
            ! else
            !     massdis = 0
            ! end if

            if (rfi .ne. 0.) then
                massfil = mcoal
            else
                massfil = 0
            end if

            if (rsh .ne. 0.) then
                masshee = mcoal
            else
                masshee = 0
            end if
            ! mat1 = massfil
            ! mat2 = massdis
            ! mat3 = masshee

            cmasfi1=0.
            cmasfi2=0.
            cmasfi3=0.
            cmasdi1=0.
            cmasdi2=0.
            cmassh1=0.
            cmassh2=0.
            cnbrfi1=0.
            cnbrfi2=0.
            cnbrfi3=0.
            cnbrdi1=0.
            cnbrdi2=0.
            cnbrsh1=0.
            cnbrsh2=0.

            hfi1=0.
            hfi2=0.
            hfi3=0.
            sgfi1=0.
            sgfi2=0.
            sgfi3=0.
            mufi1=0.
            mufi2=0.
            mufi3=0.
            hsh1=0.
            hsh2=0.
            sgsh1=0.
            sgsh2=0.
            mush1=0.
            mush2=0.
            hdi1=0.
            hdi2=0.
            sgdi1=0.
            sgdi2=0.
            mudi1=0.
            mudi2=0.

            dff3=0.
            hf3i=0.
            sfi3=0.
            sf3i=0.
            mf3i=0.
            pf3i=0.
            pfi301=0.
            pfi302=0.
            pfi303=0.
            dss2=0.
            hs2i=0.
            ! ss2i=0.
            ! ssh2=0.
            ms2i=0.
            ps2i=0.
            ddd1=0.
            sd1i=0.
            sdi1=0.
            ms2i=0.
            ddd2=0.
            hd2i=0.
            sd2i=0.
            sdi1=0.
            md2i=0.
            pd2i=0.

            hf3=0.
            hd2=0.
            hs2=0.
            muf3=0.
            mud2=0.
            mus2=0.
            sgf3=0.
            sgd2=0.
            sgs2=0.
            pfi3m=0.

            iexfi1=0
            iexfi2=0
            iexfi3=0
            iexsh1=0
            iexsh2=0
            iexdi1=0
            iexdi2=0

            ! pfi : Filament Breakup (MF04)
            if (rfi .ne. 0) then
                if(((j >= 2).and.(j <= num_curr_bins-1)).and.((k >= 2).and.(k <= num_curr_bins-1))) then
                    if(ffi > 2d0) then ! ffilimit=2.
                        iff1=1
                        ! filament lognormal function
                        mufi3=-3.884+1.115*(1.-exp(-6.60e+06*(st-sc)))
                        sgfi3min=0.408-0.105*(3.-1.546)**2.
                        sgfi3=0.408-0.105*(w2-1.546)**2.
                        if(sgfi3.le.sgfi3min) then
                            sgfi3=sgfi3min
                        else
                            sgfi3=sgfi3
                        endif
                        hfi3=sqrt(2./pi)*(ffi-2d0)/sgfi3/(1.+calculate_errfunct((mufi3-log(dcut))/(sqrt(2d0)*sgfi3)))
                        ! filament gaussian function
                        hfi2max=54.1*0.15**(-0.851)+10.4
                        if (w2<=0) then
                            print*, "ERROR: negative w2!", " w2",w2,"cke",cke,d(j),d(k)
                            call stopcm1
                        end if
                        hfi2=54.1*w2**(-0.851)+10.4
                        if(hfi2.gt.hfi2max) then
                            hfi2=hfi2max
                        else
                            hfi2=hfi2
                        endif
                        sgfi2=1./(sqrt(2.*pi)*hfi2)
                        ! filament dirac function
                        hfi1=1.
                        ! mass conservation verification + scaling
                        cmasfi1=hfi1*ml
                        cmasfi2=hfi2*pi/6.*calculate_gausint(dcut,dcoal,ds,sqrt(2d0)*sgfi2,3)
                        cmasfi3=hfi3*pi/6.*calculate_lgnrint(dcut,dcoal,mufi3,sqrt(2d0)*sgfi3,2d0)
                        if(cmasfi3 >= scfi .and. if_scaling_VY85) then
                            hfi3=hfi3*scfi/cmasfi3
                            cmasfi3=scfi
                        endif
                        scf2=(mcoal-cmasfi3)-(cmasfi2+cmasfi1)
                        if(scf2.ne.0.) then
                            hfi1=hfi1*(1.+scf2/(cmasfi1+cmasfi2))
                            hfi2=hfi2*(1.+scf2/(cmasfi1+cmasfi2))
                        endif
                        cmasfi1=hfi1*ml
                        cmasfi2=hfi2*pi/6.*calculate_gausint(dcut,dcoal,ds,sqrt(2d0)*sgfi2,3)
                        cmasfi3=hfi3*pi/6.*calculate_lgnrint(dcut,dcoal,mufi3,sqrt(2d0)*sgfi3,2d0)
                        cnbrfi2=hfi2*calculate_gausint(dcut,dcoal,ds,sqrt(2d0)*sgfi2,0)
                        cnbrfi3=hfi3*calculate_lgnrint(dcut,dcoal,mufi3,sqrt(2d0)*sgfi3,-1d0)
                        iexfi3=1
                        iexfi2=1
                        iexfi1=1
                    else ! ffi <= ffilimit
                        iff1=2
                        cmasfi3=0.
                        cmasfi2=ms
                        cmasfi1=ml
                        sgfi3=0.
                        mufi3=0.
                        hfi3=0.
                        hfi2=cmasfi2/ms
                        hfi1=cmasfi1/ml
                        iexfi3=2
                        iexfi2=3
                        iexfi1=1
                    endif
                else ! At least ne of the interactive bins is at the boundary
                    iff1=3
                    cmasfi3=0.
                    cmasfi2=ms
                    cmasfi1=ml
                    sgfi3=0.
                    mufi3=0.
                    hfi3=0.
                    hfi2=cmasfi2/ms
                    hfi1=cmasfi1/ml
                    iexfi3=2
                    iexfi2=3
                    iexfi1=1
                endif
                hf3=hfi3
                sgf3=sgfi3
                muf3=mufi3
            end if
            ! psh : Sheet Breakup (MF04)
            if(rsh.ne.0.) then
                if(((j >= 2).and.(j <= num_curr_bins-1)).and.((k >= 2).and.(k <= num_curr_bins-1))) then
                    if(fsh > 2d0) then ! fshlimit=2.
                        ish1=1
                        ! sheet lognormal function
                        mush2=-3.65+0.51*w2
                        sgsh2max=0.279+1.91e+04*(1.5e-05)
                        sgsh2=0.279+1.91e+04*cke
                        if(sgsh2.ge.sgsh2max) then
                            sgsh2=sgsh2max
                        else
                            sgsh2=sgsh2
                        endif
                        hsh2=sqrt(2./pi)*(fsh-1d0)/sgsh2/(1.+calculate_errfunct((mush2-log(dcut))/(sqrt(2d0)*sgsh2)))
                        if(igaus.eq.1) then ! A - Uses MF04 Gaussian function for the remaining parent drop
                            ! sheet gaussian function
                            sgsh1min=0.002
                            sgsh1=-0.00612+8.65e-04*dl/ds+3.84e+04*(st-sc)
                            if(sgsh1.le.sgsh1min) then
                                sgsh1=sgsh1min
                            else
                                sgsh1=sgsh1
                            endif
                            hsh1=1./calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,0)
                            if (ieee_is_nan(hsh1)) then
                                print *, "hsh1 is nan line 2131", calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,0)
                                call stopcm1
                            end if
                            ! mass conservation verification + scaling
                            cmassh1=hsh1*pi/6.*calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,3)
                            cmassh2=hsh2*pi/6.*calculate_lgnrint(dcut,dcoal,mush2,sqrt(2d0)*sgsh2,2d0)
                            if(cmassh2.ge.scsh .and. if_scaling_VY85) then
                                hsh2=hsh2*scsh/cmassh2
                                cmassh2=scsh
                            endif
                            scs2=(mcoal-cmassh2)-cmassh1
                            if(scs2.ne.0.) then
                                cmassh1=cmassh1+scs2
                                hsh1=cmassh1*(6./pi)/calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,3)
                                if (ieee_is_nan(hsh1)) then
                                    print *, "hsh1 is nan line 2145", calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,3), dcut, dcoal, dl, sgsh1, ds
                                    call stopcm1
                                end if
                            endif
                            cmassh1=hsh1*pi/6.*calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,3)
                            cmassh2=hsh2*pi/6.*calculate_lgnrint(dcut,dcoal,mush2,sqrt(2d0)*sgsh2,2d0)
                            cnbrsh1=hsh1*calculate_gausint(dcut,dcoal,dl,sqrt(2d0)*sgsh1,0)
                            cnbrsh2=hsh2*calculate_lgnrint(dcut,dcoal,mush2,sqrt(2d0)*sgsh2,-1d0)
                            iexsh2=1
                            iexsh1=1
                        else if(igaus.eq.9) then ! B - Uses a Delta function for the remaining parent drop
                            cmassh2=hsh2*pi/6.*calculate_lgnrint(dcut,dcoal,mush2,sqrt(2d0)*sgsh2,2d0)
                            ! scaling
                            if(cmassh2.ge.scsh .and. if_scaling_VY85) then
                                ish1=3
                                hsh2=hsh2*scsh/cmassh2
                                cmassh2=scsh
                            endif
                            cnbrsh2=hsh2*calculate_lgnrint(dcut,dcoal,mush2,sqrt(2d0)*sgsh2,-1d0)
                            ! sheet delta function
                            cmassh1=mcoal-cmassh2
                            hsh1=cmassh1/ml
                            iexsh2=1
                            iexsh1=3
                        endif
                    else ! fsh <= fshlimit
                        ish1=2
                        cmassh2=ms
                        cmassh1=ml
                        sgsh2=0.
                        mush2=0.
                        hsh2=cmassh2/ms
                        hsh1=cmassh1/ml
                        iexsh2=3
                        iexsh1=3
                    endif
                else
                    ish1=3
                    cmassh2=ms
                    cmassh1=ml
                    sgsh2=0.
                    mush2=0.
                    hsh2=cmassh2/ms
                    hsh1=cmassh1/ml
                    iexsh2=3
                    iexsh1=3
                endif
                sgs2=sgsh2
                mus2=mush2
                if(ish1.eq.1) then
                    hs2=hsh2
                else
                    hs2=0.
                end if
            endif

            ! pdi : Disk Breakup (MF04)
            if(rdi.ne.0d0) then
                if(((j >= 2).and.(j <= num_curr_bins-1)).and.((k >= 2).and.(k <= num_curr_bins-1))) then
                    if(fdi > 2d0) then ! fdilimit=2
                        idi1=1
                        ! disc lognormal function
                        mudi2min=-2.905+1.224*(0.2)
                        mudi2max=-2.905+1.224*(0.55)
                        mudi2=-2.905+1.224*ds/dl
                        if(mudi2.ge.mudi2max) then
                            mudi2=mudi2max
                        else if(mudi2.le.mudi2min) then
                            mudi2=mudi2min
                        else
                            mudi2=mudi2
                        endif
                        sgdi2min=0.1
                        sgdi2=0.0949+2.131*dl-2.986*ds
                        if(sgdi2.le.sgdi2min) then
                            sgdi2=sgdi2min
                        else
                            sgdi2=sgdi2
                        endif
                        hdi2=sqrt(2./pi)*(fdi-1d0)/sgdi2/(1.+calculate_errfunct((mudi2-log(dcut))/(sqrt(2d0)*sgdi2)))
                        ! disc gaussian function
                        ddd1=dl*(1.-exp(-3.7*(3.1-w1)))
                        sgdi1min=-0.0233+6.76e+03*(4.1e-06)
                        sgdi1max=-0.0233+6.76e+03*(1.5e-05)
                        sgdi1=-0.0233+6.76e+03*cke
                        if(cke >= 1.5e-05) then
                            sgdi1=sgdi1max
                        else if(cke <= 4.1e-06) then
                            sgdi1=sgdi1min
                        else
                            sgdi1=sgdi1
                        endif
                        if (calculate_gausint(dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,0) == 0) then
                            print *, "ERROR: gausint is 0", dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,0
                            print *, "ds, dl, vs, vl, cke, sc, st", ds, dl, vs, vl, cke, sc, st
                            call stopcm1
                        end if
                        hdi1=1./calculate_gausint(dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,0)
                        ! mass conservation verification + scaling
                        cmasdi1=hdi1*pi/6.*calculate_gausint(dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,3)
                        cmasdi2=hdi2*pi/6.*calculate_lgnrint(dcut,dcoal,mudi2,sqrt(2d0)*sgdi2,2d0)
                        scdi=mcoal-(cmasdi2+cmasdi1)
                        if(scdi.ne.0.) then
                            hdi1=hdi1*(1.+scdi/(mcoal-scdi))
                            hdi2=hdi2*(1.+scdi/(mcoal-scdi))
                        endif
                        cmasdi1=hdi1*pi/6.*calculate_gausint(dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,3)
                        cmasdi2=hdi2*pi/6.*calculate_lgnrint(dcut,dcoal,mudi2,sqrt(2d0)*sgdi2,2d0)
                        cnbrdi1=hdi1*calculate_gausint(dcut,dcoal,ddd1,sqrt(2d0)*sgdi1,0)
                        cnbrdi2=hdi2*calculate_lgnrint(dcut,dcoal,mudi2,sqrt(2d0)*sgdi2,-1d0)
                        iexdi1=1
                        iexdi2=1
                    else ! fdi <= fdilimit
                        idi1=2
                        cmasdi2=ms
                        cmasdi1=ml
                        sgdi2=0.
                        mudi2=0.
                        hdi2=cmasdi2/ms
                        hdi1=cmasdi1/ml
                        iexdi2=3
                        iexdi1=3
                    endif
                else
                    idi1=3
                    cmasdi2=ms
                    cmasdi1=ml
                    sgdi2=0.
                    mudi2=0.
                    hdi2=cmasdi2/ms
                    hdi1=cmasdi1/ml
                    iexdi2=3
                    iexdi1=3
                endif
                sgd2=sgdi2
                mud2=mudi2
                if(idi1.eq.1) then
                    hd2=hdi2
                else
                    hd2=0d0
                endif
            end if

            ra2di1 = 0 ! LIHUI added
            ra2di2 = 0 ! LIHUI added
            ! Breakup part function (explicit) for number and mass preservation (KR96) : kappa_brkp
            do i=1, num_curr_bins
                if (x(i) > xcoal) then
                    kappa_brkp(i,j,k) = 0d0
                else
                    ! Function Ra2 in the Ramkrishna Scheme for Disc function 1
                    if(iexdi1.eq.1) then
                        pdi1=hdi1*exp(-((d(i)-ddd1)/(sqrt(2d0)*sgdi1))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2di1=hdi1/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),ddd1,sqrt(2d0)*sgdi1,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),ddd1,sqrt(2d0)*sgdi1,3))
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2di1=hdi1/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),ddd1,sqrt(2d0)*sgdi1,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),ddd1,sqrt(2d0)*sgdi1,3)) &
                            +hdi1/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),ddd1,sqrt(2d0)*sgdi1,3)-x(i-1)*calculate_gausint(d(i-1),d(i),ddd1,sqrt(2d0)*sgdi1,0))
                        else if(i.eq.num_curr_bins) then
                            ra2di1=hdi1/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),ddd1,sqrt(2d0)*sgdi1,3)-x(i-1)*calculate_gausint(d(i-1),d(i),ddd1,sqrt(2d0)*sgdi1,0))
                        endif
                    else if(iexdi1.eq.3) then
                        if(x(i).eq.xl) then
                            pdi1=hdi1
                            ra2di1=hdi1
                        else if(x(i).ne.xl) then
                            pdi1=0.
                            ra2di1=0.
                        endif
                    endif

                    ! Function Ra2 in the Ramkrishna Scheme for Disc function 2
                    if(iexdi2.eq.1) then
                        pdi2=hdi2/d(i)*exp(-((log(d(i))-mudi2)/(sqrt(2d0)*sgdi2))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2di2=hdi2/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mudi2,sqrt(2d0)*sgdi2,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mudi2,sqrt(2d0)*sgdi2,2d0))
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2di2=hdi2/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mudi2,sqrt(2d0)*sgdi2,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mudi2,sqrt(2d0)*sgdi2,2d0)) &
                            +hdi2/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mudi2,sqrt(2d0)*sgdi2,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mudi2,sqrt(2d0)*sgdi2,-1d0))
                        else if(i.eq.num_curr_bins) then
                            ra2di2=hdi2/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mudi2,sqrt(2d0)*sgdi2,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mudi2,sqrt(2d0)*sgdi2,-1d0))
                        endif
                    else if(iexdi2.eq.3) then
                        if(x(i).eq.xs) then
                            pdi2=hdi2
                            ra2di2=hdi2
                        else if(x(i).ne.xs) then
                            pdi2=0.
                            ra2di2=0.
                        endif
                    endif
                    ra2fi1 = 0 ! LIHUI added
                    ra2fi2 = 0 ! LIHUI added
                    ra2fi3 = 0 ! LIHUI added
                    ! Function Ra2 in the Ramkrishna Scheme for Fila function 1
                    if(iexfi1.eq.1) then 
                        if(x(i).eq.xl) then
                            pfi1=hfi1
                            ra2fi1=hfi1
                        else if(x(i).ne.xl) then
                            pfi1=0.
                            ra2fi1=0.
                        endif
                    else if(iexfi1.ne.1) then 
                        pfi1=0.
                        ra2fi1=0.
                    endif

                    ! Function Ra2 in the Ramkrishna Scheme for Fila function 2
                    if(iexfi2.eq.1) then
                        pfi2=hfi2*exp(-((d(i)-ds)/(sqrt(2d0)*sgfi2))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2fi2=hfi2/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,3))
                            if (ra2fi2 < -1e-1) then
                                print *, d(i), d(i+1), ds, sqrt(2d0)*sgfi2
                                print *, calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,0), calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,1)/d(i+1), calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,2)/d(i+1)**2, calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,3)/d(i+1)**3
                                print *, "negative ra2fi2", ra2fi2, i, j, k, x(i), x(i+1), d(i), d(i+1), (pi/6.)*d(i)**3, (pi/6.)*d(i+1)**3, x(i+1)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,0), (pi/6.)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,3), (x(i+1)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,3))
                                call stopcm1
                            end if
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2fi2=hfi2/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),ds,sqrt(2d0)*sgfi2,3)) &
                            +hfi2/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),ds,sqrt(2d0)*sgfi2,3)-x(i-1)*calculate_gausint(d(i-1),d(i),ds,sqrt(2d0)*sgfi2,0))
                        else if(i.eq.num_curr_bins) then
                            ra2fi2=hfi2/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),ds,sqrt(2d0)*sgfi2,3)-x(i-1)*calculate_gausint(d(i-1),d(i),ds,sqrt(2d0)*sgfi2,0))
                        endif
                    else if(iexfi2.eq.3) then
                        if(x(i).eq.xs) then
                            pfi2=hfi2
                            ra2fi2=hfi2
                        else if(x(i).ne.xs) then
                            pfi2=0.
                            ra2fi2=0.
                        endif
                    endif

                    ! Function Ra2 in the Ramkrishna Scheme for Fila function 3
                    if(iexfi3.eq.1) then
                        pfi3=hfi3/d(i)*exp(-((log(d(i))-mufi3)/(sqrt(2d0)*sgfi3))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2fi3=hfi3/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mufi3,sqrt(2d0)*sgfi3,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mufi3,sqrt(2d0)*sgfi3,2d0))
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2fi3=hfi3/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mufi3,sqrt(2d0)*sgfi3,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mufi3,sqrt(2d0)*sgfi3,2d0)) &
                            +hfi3/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mufi3,sqrt(2d0)*sgfi3,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mufi3,sqrt(2d0)*sgfi3,-1d0))
                        else if(i.eq.num_curr_bins) then
                            ra2fi3=hfi3/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mufi3,sqrt(2d0)*sgfi3,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mufi3,sqrt(2d0)*sgfi3,-1d0))
                        endif
                    else if(iexfi3.ne.1) then
                        pfi3=0.
                        ra2fi3=0.
                    endif
                    ra2sh1 = 0 ! LIHUI added
                    ra2sh2 = 0 ! LIHUI added
                    ! Function Ra2 in the Ramkrishna Scheme for Shee function 1
                    if(iexsh1.eq.1) then
                        psh1=hsh1*exp(-((d(i)-dl)/(sqrt(2d0)*sgsh1))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2sh1=hsh1/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,3))
                            if (ieee_is_nan(ra2sh1)) then
                                print *, i, hsh1, x(i+1), x(i), calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,0), (pi/6.)*calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,3)
                                call stopcm1
                            end if
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2sh1=hsh1/(x(i+1)-x(i))*(x(i+1)*calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,0)-(pi/6.)*calculate_gausint(d(i),d(i+1),dl,sqrt(2d0)*sgsh1,3)) &
                            +hsh1/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),dl,sqrt(2d0)*sgsh1,3)-x(i-1)*calculate_gausint(d(i-1),d(i),dl,sqrt(2d0)*sgsh1,0))
                        else if(i.eq.num_curr_bins) then
                            ra2sh1=hsh1/(x(i)-x(i-1))*((pi/6.)*calculate_gausint(d(i-1),d(i),dl,sqrt(2d0)*sgsh1,3)-x(i-1)*calculate_gausint(d(i-1),d(i),dl,sqrt(2d0)*sgsh1,0))
                        endif
                    else if(iexsh1.eq.3) then
                        if(x(i).eq.xl) then
                            psh1=hsh1
                            ra2sh1=hsh1
                        else if(x(i).ne.xl) then
                            psh1=0.
                            ra2sh1=0.
                        endif
                    endif

                    ! Function Ra2 in the Ramkrishna Scheme for Shee function 2
                    if(iexsh2.eq.1) then
                        psh2=hsh2/d(i)*exp(-((log(d(i))-mush2)/(sqrt(2d0)*sgsh2))**2d0)*(d1(i+1)-d1(i))
                        if(i.eq.1) then
                            ra2sh2=hsh2/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mush2,sqrt(2d0)*sgsh2,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mush2,sqrt(2d0)*sgsh2,2d0))
                        else if((i > 1).and.(i < num_curr_bins)) then
                            ra2sh2=hsh2/(x(i+1)-x(i))*(x(i+1)*calculate_lgnrint(d(i),d(i+1),mush2,sqrt(2d0)*sgsh2,-1d0)-(pi/6.)*calculate_lgnrint(d(i),d(i+1),mush2,sqrt(2d0)*sgsh2,2d0)) &
                            +hsh2/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mush2,sqrt(2d0)*sgsh2,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mush2,sqrt(2d0)*sgsh2,-1d0))
                        else if(i.eq.num_curr_bins) then
                            ra2sh2=hsh2/(x(i)-x(i-1))*((pi/6.)*calculate_lgnrint(d(i-1),d(i),mush2,sqrt(2d0)*sgsh2,2d0)-x(i-1)*calculate_lgnrint(d(i-1),d(i),mush2,sqrt(2d0)*sgsh2,-1d0))
                        endif
                    else if(iexsh2.eq.3) then
                        if(x(i).eq.xs) then
                            psh2=hsh2
                            ra2sh2=hsh2
                        else if(x(i).ne.xs) then
                            psh2=0.
                            ra2sh2=0.
                        endif
                    endif
                    kappa_brkp(i,j,k) = max(0d0, rdi*(ra2di1+ra2di2) + rfi*(ra2fi1+ra2fi2+ra2fi3) + rsh*(ra2sh1+ra2sh2))
                    if (ieee_is_nan(kappa_brkp(i,j,k))) then
                        print *, "nan kappa", kappa_brkp(i,j,k)
                        print *, i, j, k
                        print *, rdi, ra2di1, ra2di2, rfi, ra2fi1, ra2fi2, ra2fi3, rsh, ra2sh1, ra2sh2, hfi2
                        print *, "iexfi2:", iexfi2
                        call stopcm1
                    end if
                end if
            end do

            ! mass conservation adjustment
            coeff_ms = xcoal/sum(kappa_brkp(1:num_curr_bins, j, k) * x)
            if (ieee_is_nan(coeff_ms)) then
                print *, "result negative mass conservation coefficients"
                call stopcm1
            end if
            kappa_brkp(1:num_curr_bins, j, k) = kappa_brkp(1:num_curr_bins, j, k) * coeff_ms

            ! check mass conservation
            rate_ms = abs(sum(kappa_brkp(1:num_curr_bins, j, k) * x)/xcoal)           
            
            if (abs(rate_ms-1d0) > 0.01) then
                print *, "ERROR: mass is not conserved for collision between bin ", j, " and ", k
                print *, "mass conservation between bin ", j, d(j), x(j), " and ", k, d(k), x(k), ":", "xcoal:", xcoal, rate_ms
                print *, kappa_brkp(1:num_curr_bins, j, k) 
                call stopcm1
                ! do i=1, num_curr_bins
                !     print *, i, kappa_brkp(i,j,k), d(i), x(i)
                ! end do
                ! stop    
            end if
            kernel_mass_conservation(j,k) = rate_ms
        end do
    end do
    ! print *, "the maximum mass conservation: ", maxval(kernel_mass_conservation), "the minimum mass conservation: ", minval(kernel_mass_conservation), "the mean mass conservation:", sum(kernel_mass_conservation)/num_curr_bins**2

    ! open(unit=10086, file="breakup_kappa_examples.txt", status="new", action="write")
    !     write(10086, *) d
    !     write(10086, *) d(99), d(57), new_line("a"), kappa_brkp(1:num_curr_bins,99,57), new_line("a")
    !     write(10086, *) d(80), d(66), new_line("a"), kappa_brkp(1:num_curr_bins,80,66), new_line("a")
    ! close(10086)

end subroutine calculate_kappa_brkp


function calculate_errfunct(x) result(errfunct)
real (kind=dp) :: errfunct,x
    !   Returns the error function erf(x).
    if(x.lt.0.)then
    errfunct=-calculate_gammp(0.5d0,x**2d0)
    else
    errfunct=calculate_gammp(0.5d0,x**2d0)
    endif
return
end function calculate_errfunct

function calculate_gammp(a,x) result (gammp)
    !   USES gcf,gser
    !   Returns the incomplete gamma function P(a, x). From Press et al. (1992)
    real (kind=dp) :: a,gammp,x
    real (kind=dp) :: gammcf,gamser,gln
    if(x.lt.0..or.a.le.0.) then
        print *, 'bad arguments in gammp'
        call stopcm1
    end if
    if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
    else 
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
    endif

end function calculate_gammp

function calculate_lgnrint(x1,x2,a,b,c) result(lgnrint)
    real (kind=dp) :: lgnrint,a,b,c,x1,x2
    ! USES erf from Press et al. (1992)
    ! Returns the definite integral of generalized lognormal functions between limits x1 (lower) and x2 (upper) Output=Integral[x^c.exp{-{lnx-a)/b}^2}.dx]
    lgnrint=sqrt(pi)/2.*b*exp((c+1.)*(a+(c+1.)/4.*b**2d0))*(calculate_errfunct((log(x2)-a)/b-(c+1.)/2.*b)-calculate_errfunct((log(x1)-a)/b-(c+1.)/2.*b))

end function calculate_lgnrint

function calculate_gausint(x1,x2,a,b,n) result(gausint)
! USES erf,bico from Press et al. (1992) 
! Returns the definite integral of generalized gaussian function between limits x1 (lower) and x2 (upper) Output=Integral[x^n.exp{-{(x-a)/b}^2}.dx]
    real (kind=dp) :: gausint,x1,x2,a,b
    integer :: n
    real (kind=dp) :: t1,t2
    real (kind=dp) :: p0,p1,pn,pm
    real (kind=dp) :: ak,bk
    real (kind=dp) :: bico
    integer :: k,k1,k2,l1,l
    t1=(x1-a)/b
    t2=(x2-a)/b
    k1=0
    l1=0
    bk=0.
    p0=.5*sqrt(pi)*(calculate_errfunct(t2)-calculate_errfunct(t1))
    p1=-.5*(exp(-t2**2d0)-exp(-t1**2d0))
    do 11 k=0,n
        if(k.eq.0) then
            ak=p0
        else if(k.eq.1) then
            ak=p1
        else if(k.ge.2) then
                k1=k
                k1=k1-2
                l1=l1+1
1               if(k1.eq.0)then
                    k2=0
                    goto 2
                else if(k1.eq.1)then
                    k2=1
                    goto 2
                else if((k1.ne.0).or.(k1.ne.1))then 
                    k1=k1-2
                    l1=l1+1
                    goto 1
2               endif
            l1=l1
            if(k2.eq.0)pn=p0
            if(k2.eq.1)pn=p1
            do 12 l=1,l1
                pm=pn
                pn=-.5*(x2**real(2*l+k2-1, kind=dp)*exp(-t2**2d0)-x1**real(2*l+k2-1, kind=dp)*exp(-t1**2d0)-real(2*l+k2-1, kind=dp)*pm)
12          enddo
            ak=pn
        endif
        bk=bk+calculate_bico(n,k)*a**real(n-k, kind=dp)*b**real(k+1, kind=dp)*ak
    l1=0
11  continue
    gausint = bk
end function calculate_gausint

function calculate_bico(n,k) result(bico)
	integer :: k,n
	real (kind=dp) :: bico
	bico=nint(exp(calculate_factln(n)-calculate_factln(k)-calculate_factln(n-k)))
end function calculate_bico

FUNCTION calculate_factln(n) result(factln)
    ! Returns ln(n!).
    integer :: n
    real (kind=dp) :: factln
	real (kind=dp) :: a(100)
	SAVE a
	DATA a/100*-1./
	if (n.lt.0) then
        print *, 'negative factorial in factln'
        call stopcm1
    end if
	if (n.le.99) then
		if (a(n+1).lt.0.) then
            a(n+1)=calculate_gammln(n+1. + 0d0)
        end if
		factln=a(n+1)
	else
		factln=calculate_gammln(n+1. + 0d0)
	endif
end function calculate_factln

function calculate_gammln(xx) result(gammln)
    implicit none
    !  gammln retruns the natural logarithm of Gamma(xx).  From
    !  Section 6.1, p 207 of Numerical Recipes in Fortran, 2nd ed.
    !  adapted from https://github.com/kkadam/ms_hydro/blob/master/gammln.F90 (C) Copr. 1986-92 Numerical Recipes Software .
    real (kind=dp) :: gammln,xx
    integer :: j
    real (kind=dp) :: ser,stp,tmp,x,y,cof(6)
    SAVE cof,stp
    data cof,stp/76.18009172947146e0,-86.50532032941677e0,             &
    24.01409824083091e0,-1.231739572450155e0,.1208650973866179e-2,     &
    -.5395239384953e-5,2.5066282746310005e0/
    x=xx
    y=x
    tmp=x+5.5e0
    tmp=(x+0.5e0)*log(tmp)-tmp
    ser=1.000000000190015e0
    do j= 1, 6
        y=y+1.e0
        ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
end function calculate_gammln

subroutine gcf(gammcf,a,x,gln)
	integer :: ITMAX
	real (kind=dp) :: a,gammcf,gln,x,EPS,FPMIN
	PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
	integer :: i
	real (kind=dp) :: an,b,c,d,del,h,gammln
	gln=calculate_gammln(a)
	b=x+1.-a
	c=1./FPMIN
	d=1./b
	h=d
	do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d) < FPMIN) then
            d=FPMIN
        end if
        c=b+an/c
        if(abs(c) < FPMIN) then
            c=FPMIN
        end if
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.) < EPS) then
            goto 1
        end if
11	end do
	print *, 'a too large, ITMAX too small in gcf'
1	gammcf=exp(-x+a*log(x)-gln)*h
	return
end subroutine gcf

subroutine gser(gamser,a,x,gln)
	integer :: ITMAX
	real (kind=dp) :: a,gamser,gln,x,EPS
	PARAMETER (ITMAX=100,EPS=3.e-7)
	integer :: n
	real (kind=dp) :: ap,del,sum,gammln
	gln=calculate_gammln(a)
	if(x <= 0.)then
	    if(x < 0.) then
            print *, 'x < 0 in gser'
            call stopcm1
        end if
        gamser=0.
    else
        ap=a
        sum=1./a
        del=sum
        do 11 n=1,ITMAX
            ap=ap+1.
            del=del*x/ap
            sum=sum+del
            if(abs(del) < abs(sum)*EPS) then
                goto 1
            end if
    11	end do
        print *, 'a too large, ITMAX too small in gser'
        call stopcm1
    1	gamser=sum*exp(-x+a*log(x)-gln)
	endif
end subroutine gser



!#####################          Initialization           ###################    
    subroutine dcpm_init(qa, dcpm3d, p3d, th0, th3d, pi0, rho0, pp3d, zh, zf, irst) ! qa is then passed to q3d
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke), intent(in) :: pi0,th0,rho0,zh,pp3d,th3d,p3d
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(inout) :: qa
        real (kind=sp), dimension(ibm:iem,jbm:jem,kbm:kem,num_bins+1), intent(inout) :: dcpm3d
        integer, intent(in) :: irst
        real (kind=sp), dimension(ib:ie,jb:je,kb:ke) :: temp3d
        real (kind=dp), dimension(num_bins) :: bin_r_dry, bin_num_conc, bin_vol, bin_v_term, bin_wl
        real (kind=dp), dimension(num_bins, num_bins) :: kernels_coal_z_temp, kernels_brkp_z_temp
        real (kind=dp), dimension(num_bins,num_bins,num_bins) :: weights_coal_z_temp
        real (kind=dp), dimension(:,:,:), allocatable :: weights_brkp_z_temp
        integer :: i3d, j3d, k3d, ind_zk
        real (kind=dp) :: wv, wl, supersat, rho_d, temp, p, zf_min, zf_max, hdf
        integer :: ind_bin, ind_bin_brkp_start, ind_bin_i, ind_bin_j
        real (kind=sp), dimension(:), allocatable :: rho_kernels, rho_brkp_kernels
        integer :: iproc
        tolerance_bin_wl = tolerance_wl/(num_bins*0.1)
        threshold_bin_wl_coll = tolerance_bin_wl*10
        ! ratio_bins = (drop_radius_max/ccn_radius_min)**(1d0/real((num_bins-1))) ! ratio_bins_init: initial wet radius ratio of bins
        if (myid .eq. 0) then
            print *, "Begin initialization of dcpm"
            print *, "number of bins: ", num_bins
            print *, "minimum CCN radius: ", ccn_radius_min
            print *, "maximum drop radius: ", drop_radius_max
            if (irst .eq. 1) then
                print *, "This is a restart"
            end if
        end if


        ! call cpu_time ( t_cpu_init )
        ! call system_clock ( t_clock_init, clock_rate_dcpm, clock_max_dcpm )

        ! ******** Begin of reading/calculating properties ********

        ! First, get aerosol type correct
        if (aerosol_type .eq. 1) then
            ! SAM aerosol
            allocate(aer_num_concs(4))
            allocate(aer_gmds(4))
            allocate(aer_gsds(4))
            aer_num_concs = (/ 1401.9d6, 415.7d6, 0.3d6, 0.3d6 /)
            aer_gmds = (/ 0.076d-6, 0.195d-6, 0.75d-6, 2.2d-6 /)
            aer_gsds = (/ 1.63d0, 1.35d0, 1.30d0, 1.40d0 /)
            aer_kappa = 0.14d0 
        else if (aerosol_type .eq. 2) then
            ! Continental aerosol
            allocate(aer_num_concs(1))
            allocate(aer_gmds(1))
            allocate(aer_gsds(1))
            aer_num_concs = (/842d6/)
            aer_gmds = (/0.0218*2d-6/)
            aer_gsds = (/1.16d0/)
            aer_kappa = 0.3 ! Assuming an organic fraction around 0.4-0.6
        else if (aerosol_type .eq. 3) then
            ! SAM aerosol but without mode 3 and 4
            allocate(aer_num_concs(2))
            allocate(aer_gmds(2))
            allocate(aer_gsds(2))
            aer_num_concs = (/ 1401.9d6, 415.7d6/)
            aer_gmds = (/ 0.076d-6, 0.195d-6/)
            aer_gsds = (/ 1.63d0, 1.35d0/)
            aer_kappa = 0.14d0 
        end if
        ratio_wet_coll_min = (99*aer_kappa)**(1d0/3d0)
        if (myid.eq.0) then
            if (aerosol_type .eq. 1) then
                print *, "aerosol type = SAM"
            else if (aerosol_type .eq. 2) then
                print *, "aerosol type = continental"
            end if
            print *, "N0 (cm^(-3)) = ", aer_num_concs*1e-6
            print *, "Dm (micron) = ", aer_gmds*1e6
            print *, "std = ", aer_num_concs
            print *, "hygroscopicity = ", aer_kappa
            print *, "minimum factor (r_wet/r_dry) for collision = ", ratio_wet_coll_min
        end if
        ! Now, only allocate q3d_dcpm. The array will be updated in dcpm_calculate subroutine
        allocate(bin_r_wet(1:num_bins))
        allocate(bin_r_dry_init_vol_weight(1:num_bins))
        allocate(bin_num_conc_base(1:num_bins))
        allocate(bin_r_wet_check(1:num_bins))
        allocate(bin_r_wet_grid(1:num_bins+1))
        allocate(bin_r_wet_grid_check(1:num_bins+1))
        allocate(q3d_dcpm(ibm:iem,jbm:jem,kbm:kem,numq))
        ratio_bins = (drop_radius_max/ccn_radius_min)**(1d0/real(num_bins))

        call input_turbulence_kernel (turbulence_kernel)
        call input_collision_kernels (rd1, rd2, collision_kernel_0000_1, collision_kernel_2500_1, collision_kernel_5000_1, collision_kernel_0000_2, collision_kernel_2500_2, collision_kernel_5000_2, collision_kernel_0000_3, collision_kernel_2500_3, collision_kernel_5000_3)
        call make_bin_r_wet (bin_r_wet)
        call make_bin_r_wet_grid(bin_r_wet_grid, bin_r_wet)  
        bin_r_wet_check = bin_r_wet
        bin_r_wet_grid_check = bin_r_wet_grid
        call fill_bin_num_conc(bin_r_wet,bin_num_conc_base)
        bin_r_dry_init_vol_weight = bin_r_wet**3*bin_num_conc_base
        bin_r_dry_init_vol_weight = bin_r_dry_init_vol_weight/sum(bin_r_dry_init_vol_weight)
        do ind_bin=1, num_bins
            if (bin_r_wet(ind_bin) > ccn_radius_max) then
                ind_bin_base_empty_start = ind_bin
                exit
            end if
        end do
        ! print *, calculate_equilibrium_dwl(1.69442017E-02, 294.279816 , 94579.3750)
        ! stop

        if ((irst .eq. 0) .or. (if_restart_morrison .and. (irst.eq.1))) then
            if (myid.eq.0) then
                print *, ni, nj, nk
            end if        
            do i3d=1, ni
                do j3d=1, nj
                    do k3d=1, nk
                        wv = qa(i3d,j3d,k3d,nqv)
                        wl = qa(i3d,j3d,k3d,nql1)
                        temp3d(i3d,j3d,k3d) = (th0(i3d,j3d,k3d)+th3d(i3d,j3d,k3d))*(pi0(i3d,j3d,k3d)+pp3d(i3d,j3d,k3d))
                        temp = temp3d(i3d,j3d,k3d)
                        ! print *, temp
                        p = p3d(i3d,j3d,k3d)
                        supersat = calculate_supersat (temp, p, wv)
                        rho_d = calculate_rho_d (temp, p, wv)
                        dcpm3d(i3d,j3d,k3d,num_bins+1) = supersat
                        hdf = exp((-(zh(i3d,j3d,k3d)-zf(i3d,j3d,1))/Hs_aerosol))

                        if ((if_DSD_lookup .eqv. .false.) .or. (wl < tolerance_wl)) then
                            call preprocess_bins (bin_r_dry, bin_r_wet, bin_num_conc, bin_r_wet_grid, supersat, temp)
                            bin_num_conc = bin_num_conc*hdf*scale_num_conc
                            ! print *, bin_num_conc, "hdf:",hdf, scale_num_conc
                            ! call stopcm1
                            dcpm3d(i3d,j3d,k3d,1:num_bins) = bin_r_dry
                            do ind_bin=1, num_bins
                                if (ieee_is_nan(bin_r_dry(ind_bin))) then
                                    print *, "got nan in r dry initialization"
                                    call stopcm1
                                end if
                            end do
                            bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet**3d0 - bin_r_dry**3d0) *rho_w/rho_d
                            qa(i3d,j3d,k3d,nql1) = 0d0
                        else
                            call interpolate_DSD (hdf, wl, wv, temp, p, bin_num_conc, bin_wl)
                        end if
                        ! if (i3d*j3d*k3d .eq. 1) then
                        !     print *, bin_num_conc
                        ! end if
                        qa(i3d,j3d,k3d,nnc1+num_bins:nnc2) = bin_num_conc
                        qa(i3d,j3d,k3d,nnc1:nnc1+num_bins-1) = bin_wl*scale_bin_wl
                        qa(i3d,j3d,k3d,nql1) = sum(bin_wl)
                        if (if_DSD_lookup .eqv. .false.) then
                            qa(i3d,j3d,k3d,nqs1) = 0d0
                        end if
                    enddo
                    ! call stopcm1
                enddo

            enddo

            ! call stopcm1
            ! print *, "At initialization, the maximum supersaturation is ", maxval(dcpm3d(1:ni,1:nj,1:nk,num_bins+1))
            ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            ! call stopcm1
        end if  
        
        q3d_dcpm = qa
        ! Coalescence kernel and weights initialization
        do ind_bin =1, num_bins
            if (bin_r_wet(ind_bin) > drop_brkp_fragment_min) then
                ind_bin_brkp_start = ind_bin
                exit
            end if
        end do
        
        ! nk_kernel = 2*nk
        allocate(z_kernels(1:nk_kernel))
        allocate(z_brkp_kernels(1:nk_brkp_kernel))
        allocate(rho_kernels(1:nk_kernel))
        allocate(rho_brkp_kernels(1:nk_brkp_kernel))
        allocate(kernels_coal(1:nk_kernel, num_bins, num_bins))
        allocate(kernels_brkp(1:nk_kernel, num_bins, num_bins))
        allocate(weights_coal(1:num_bins,1:num_bins,1:num_bins))
        allocate(weights_brkp(1:nk_brkp_kernel, 1:num_bins,1:num_bins,1:num_bins))
        allocate(weights_brkp_z_temp(num_bins-ind_bin_brkp_start+1,num_bins-ind_bin_brkp_start+1,num_bins-ind_bin_brkp_start+1))
        zf_min = minval(zf(1:ni,1:nj,1:nk+1))
        zf_max = maxval(zf(1:ni,1:nj,1:nk+1))
        if (zf_min .lt. 0) then
            print *, "ERROR: negative zf", zf_min
            call stopcm1
        end if
        do ind_zk = 1, nk_kernel
            z_kernels(ind_zk) = zf_min+(zf_max-zf_min)/(nk_kernel-1)*(ind_zk-1)
        end do
        do ind_zk = 2, nk_kernel-1
            rho_kernels(ind_zk) = interpolate_zmean(rho0,zh,z_kernels(ind_zk))
        end do
        call extrapolate_endpoints(z_kernels, rho_kernels, nk_kernel)

        if (nk_brkp_kernel .eq. 1) then
            z_brkp_kernels(1) = zf_min+1000 ! 1 km height
            rho_brkp_kernels(1) = interpolate_zmean(rho0,zh,z_brkp_kernels(1))
        else
            do ind_zk = 1, nk_brkp_kernel
                z_brkp_kernels(ind_zk) = zf_min+(zf_max-zf_min)/(nk_brkp_kernel-1)*(ind_zk-1)
            end do
            do ind_zk = 2, nk_brkp_kernel-1
                rho_brkp_kernels(ind_zk) = interpolate_zmean(rho0,zh,z_brkp_kernels(ind_zk))
            end do
            call extrapolate_endpoints(z_brkp_kernels, rho_brkp_kernels, nk_brkp_kernel)
        end if
        
        do ind_zk = 1, nk_kernel-1
            if (z_kernels(ind_zk) >= z_kernels(ind_zk+1) .or. rho_kernels(ind_zk) <= rho_kernels(ind_zk+1)) then
                print *, "ERROR: z_kernels do not increase or rho_kernels do not decrease with index"
                print *, "myid", myid, "ind_zk", ind_zk, new_line("a"), "z_kernels:", z_kernels, new_line("a"), "rho_kernels:", rho_kernels
                call stopcm1
            end if
        end do                         
        if (myid.eq.0) then
            print *, "Finish initialing lookup height kernels. Kernel number: ", nk_kernel
            print *, "Exemplary kernels", new_line("a"), "z_kernels:", z_kernels, new_line("a"), "rho_kernels:", rho_kernels
            print *, "Finish initialing lookup height weights kernels. Kernel number: ", nk_brkp_kernel
            print *, "Exemplary kernels", new_line("a"), "z_brkp_kernels:", z_brkp_kernels, new_line("a"), "rho_brkp_kernels:", rho_brkp_kernels
        end if
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD, ierr)
#endif
        print *, "myid", myid, "dcpm max num conc:", maxval(qa(1:ni,1:nj,1:nk,nnc1+num_bins:nnc2)), new_line("a"), &
            "rho0 extremes:", minval(rho0(1:ni,1:nj,1:nk)), maxval(rho0(1:ni,1:nj,1:nk)), new_line("a"), &
            "rho kernel extremes:", minval(rho_kernels), maxval(rho_kernels)        
        
        if (myid.eq.0) then                                          
            print *, "myid",myid,"Initializing weights"
        end if
        weights_coal(:,:,:) = 0
        weights_brkp(:,:,:,:) = 0
        do ind_zk = 1, nk_kernel
            do ind_bin = 1, num_bins
                if (iv_fall .eq. 1) then
                    bin_v_term(ind_bin) = calculate_terminal_velocity_Best (bin_r_wet(ind_bin))
                else if (iv_fall .eq. 2) then
                    bin_v_term(ind_bin) = calculate_terminal_velocity_Lhermitte (bin_r_wet(ind_bin), rho_kernels(1)/rho_kernels(ind_zk))
                else
                    print *, "ERROR: iv_fall option unknown! iv_fall=", iv_fall
                    call stopcm1
                end if
            end do        
            call calculate_kernels_coal_brkp (kernels_coal_z_temp, kernels_brkp_z_temp, bin_r_wet, bin_v_term, num_bins, z_kernels(ind_zk ))
            if (any(abs(kernels_coal_z_temp - transpose(kernels_coal_z_temp)) > &
                    maxval(max(abs(kernels_coal_z_temp), abs(transpose(kernels_coal_z_temp)))) * 1e-6)) then
                print *, "ERROR: coalescence kernel is asymmetric"
                call stopcm1
            end if

            if (any(abs(kernels_brkp_z_temp - transpose(kernels_brkp_z_temp)) > &
                    maxval(max(abs(kernels_brkp_z_temp), abs(transpose(kernels_brkp_z_temp)))) * 1e-6)) then
                print *, "ERROR: breakup kernel is asymmetric"
                call stopcm1
            end if
            kernels_coal(ind_zk ,:,:) = kernels_coal_z_temp
            kernels_brkp(ind_zk ,:,:) = kernels_brkp_z_temp
        end do      
        do ind_zk = 1, nk_brkp_kernel
            do ind_bin = 1, num_bins
                if (iv_fall .eq. 1) then
                    bin_v_term(ind_bin) = calculate_terminal_velocity_Best (bin_r_wet(ind_bin))
                else if (iv_fall .eq. 2) then
                    bin_v_term(ind_bin) = calculate_terminal_velocity_Lhermitte (bin_r_wet(ind_bin), rho_kernels(1)/rho_brkp_kernels(ind_zk))
                else
                    print *, "ERROR: iv_fall option unknown! iv_fall=", iv_fall
                    call stopcm1
                end if
            end do   
            weights_brkp_z_temp = 0
            call calculate_kappa_brkp (weights_brkp_z_temp, bin_r_wet(ind_bin_brkp_start:num_bins), bin_r_wet_grid(ind_bin_brkp_start:num_bins+1), &
                                        bin_v_term(ind_bin_brkp_start:num_bins), num_bins-ind_bin_brkp_start+1)
            weights_brkp(ind_zk,ind_bin_brkp_start:num_bins,ind_bin_brkp_start:num_bins,ind_bin_brkp_start:num_bins) = weights_brkp_z_temp
        end do
        do ind_zk = 1, nk_brkp_kernel
            do ind_bin = 1, num_bins
                ! Elementwise max over i,j for slice k
                if (any(abs(weights_brkp(ind_zk,ind_bin,:,:) - transpose(weights_brkp(ind_zk,ind_bin,:,:))) > &
                        maxval(max(weights_brkp(ind_zk,ind_bin,:,:), transpose(weights_brkp(ind_zk,ind_bin,:,:)))) * 1e-6)) then
                    print *, "ERROR: 3D breakup weights are asymmetric at ind_bin =", ind_bin
                    call stopcm1
                end if
            end do
        end do

        do ind_bin = 1, num_bins
            bin_vol(ind_bin) = 4*pi/3 * bin_r_wet(ind_bin)**3
        end do
        weights_coal_z_temp = 0
        call calculate_weights(weights_coal_z_temp, bin_vol(1:num_bins))
        weights_coal(:,:,:) = weights_coal_z_temp        

        do ind_bin = 1, num_bins
            ! Elementwise max over i,j for slice k
            if (any(abs(weights_coal(ind_bin,:,:) - transpose(weights_coal(ind_bin,:,:))) > &
                    maxval(max(weights_coal(ind_bin,:,:), transpose(weights_coal(ind_bin,:,:)))) * 1e-6)) then
                print *, "ERROR: 3D coalescence weights are asymmetric at ind_bin =", ind_bin
                call stopcm1
            end if
        end do

        
        if (if_test .and. myid.eq.0 .and. .false.) then
            print *, "test output of breakup weights"
            do ind_zk = 1, nk_kernel
                print *, "ind_zk", ind_zk, "zk", z_kernels(ind_zk)
                print *, "wet radius of breakup bins:", bin_r_wet(ind_bin_brkp_start:)
                do ind_bin_i = ind_bin_brkp_start, num_bins
                    do ind_bin_j = ind_bin_brkp_start, num_bins
                        print *, "weights between ", ind_bin_i, ind_bin_j, bin_r_wet(ind_bin_i), bin_r_wet(ind_bin_j)
                        print *, weights_brkp(ind_zk,ind_bin_brkp_start:,ind_bin_i,ind_bin_j)
                        print *, new_line("a")
                    end do
                end do
                exit
            end do
            call stopcm1
        end if
        call MPI_BARRIER (MPI_COMM_WORLD, ierr)
        if_dcpm_initialized = .true.
        if (myid .eq. 0) then
            print *, "exit dcpm bin microphysics initialization"
        end if
    end subroutine dcpm_init

    function interpolate_zmean(array3d, z3d, ztarget) result(zmean)
        !==== Arguments ====
        real(kind=sp), intent(in), dimension(ib:ie, jb:je, kb:ke) :: array3d, z3d
        real(kind=dp), intent(in) :: ztarget
        real(kind=sp) :: zmean

        !==== Local variables ====
        integer :: i, j, k
        real(kind=sp) :: z_lower, z_upper, f_lower, f_upper
        real(kind=sp) :: weight, val_interp
        real(kind=sp) :: sum_val, count

        sum_val = 0.0_sp
        count   = 0.0_sp

        do j = 1, nj
            do i = 1, ni
                ! Skip column if ztarget is outside vertical range
                if (ztarget < z3d(i,j,1) .or. ztarget > z3d(i,j,nk)) cycle
                ! Find vertical interval
                do k = 1, nk - 1
                    if (ztarget >= z3d(i,j,k) .and. ztarget <= z3d(i,j,k+1)) then
                        ! Linear interpolation
                        z_lower = z3d(i,j,k)
                        z_upper = z3d(i,j,k+1)
                        f_lower = array3d(i,j,k)
                        f_upper = array3d(i,j,k+1)

                        weight = (ztarget - z_lower) / (z_upper - z_lower)
                        val_interp = f_lower + weight * (f_upper - f_lower)

                        sum_val = sum_val + val_interp
                        count   = count + 1.0_sp
                        exit
                    end if
                end do
            end do
        end do

        !==== Return mean or error ====
        if (count > 0.0_sp) then
            zmean = sum_val / count
        else
            write(*,*) 'Error in interpolate_zmean: no valid columns found for ztarget =', ztarget
            stop
        end if

    end function interpolate_zmean


 

    subroutine interpolate_DSD (hdf, wl, wv, temp, p, bin_num_conc, bin_wl)
        real (kind=dp), intent(in) :: hdf
        real (kind=dp), intent(in) :: wl, wv, temp, p
        real (kind=dp) :: wl_adjust, ratio
        real (kind=dp), dimension(num_bins), intent(out) :: bin_num_conc, bin_wl
        real (kind=dp), dimension(num_bins) :: bin_num_conc_left, bin_wl_left, bin_num_conc_right, bin_wl_right
        integer :: ind_hdf_key, ind_qc_key, ind_bin
        real (kind=dp) :: rho_d
        integer :: num_hdf_keys = size(hdf_keys)
        integer :: num_qc_keys = size(qc_keys)
        integer :: ind_hdf_left, ind_hdf_right, ind_qc_left, ind_qc_right, ind_hdf, ind_qc
        rho_d = calculate_rho_d (temp, p, wv)

        wl_adjust = wl
        if (wl < minval(qc_keys)) then
            print *, "ERROR: liquid water mixing ratio is not in the lookup table", wl
            call stopcm1        
        else if (wl > maxval(qc_keys)) then
            wl_adjust = maxval(qc_keys)-tolerance_wl
        end if
        if (hdf > maxval(hdf_keys) .or. hdf < minval(hdf_keys)) then
            print *, "ERROR: Height decay factor is not in the lookup table", hdf
            call stopcm1        
        end if

        do ind_hdf = 1, num_hdf_keys-1
            if (hdf>=hdf_keys(ind_hdf) .and. hdf<=hdf_keys(ind_hdf+1)) then
                ind_hdf_left = ind_hdf
                ind_hdf_right = ind_hdf+1
                exit
            end if
        end do

        do ind_qc = 1, num_qc_keys-1
            if (wl_adjust>=qc_keys(ind_qc) .and. wl_adjust<=qc_keys(ind_qc+1)) then
                ind_qc_left = ind_qc
                ind_qc_right = ind_qc+1
                exit
            end if
        end do

        bin_num_conc_left = ncbin_lookup_tbl(ind_hdf_left,ind_qc_left,1:num_bins) + (wl_adjust-qc_keys(ind_qc_left))/(qc_keys(ind_qc_right)-qc_keys(ind_qc_left))*(ncbin_lookup_tbl(ind_hdf_left,ind_qc_right,1:num_bins)-ncbin_lookup_tbl(ind_hdf_left,ind_qc_left,1:num_bins))
        bin_num_conc_right = ncbin_lookup_tbl(ind_hdf_right,ind_qc_left,1:num_bins) + (wl_adjust-qc_keys(ind_qc_left))/(qc_keys(ind_qc_right)-qc_keys(ind_qc_left))*(ncbin_lookup_tbl(ind_hdf_right,ind_qc_right,1:num_bins)-ncbin_lookup_tbl(ind_hdf_right,ind_qc_left,1:num_bins))
        bin_num_conc = bin_num_conc_left + (hdf-hdf_keys(ind_hdf_left))/(hdf_keys(ind_hdf_right)-hdf_keys(ind_hdf_left))*(bin_num_conc_right-bin_num_conc_left)

        bin_wl_left = qcbin_lookup_tbl(ind_hdf_left,ind_qc_left,1:num_bins) + (wl_adjust-qc_keys(ind_qc_left))/(qc_keys(ind_qc_right)-qc_keys(ind_qc_left))*(qcbin_lookup_tbl(ind_hdf_left,ind_qc_right,1:num_bins)-qcbin_lookup_tbl(ind_hdf_left,ind_qc_left,1:num_bins))
        bin_wl_right = qcbin_lookup_tbl(ind_hdf_right,ind_qc_left,1:num_bins) + (wl_adjust-qc_keys(ind_qc_left))/(qc_keys(ind_qc_right)-qc_keys(ind_qc_left))*(qcbin_lookup_tbl(ind_hdf_right,ind_qc_right,1:num_bins)-qcbin_lookup_tbl(ind_hdf_right,ind_qc_left,1:num_bins))
        bin_wl = bin_wl_left + (hdf-hdf_keys(ind_hdf_left))/(hdf_keys(ind_hdf_right)-hdf_keys(ind_hdf_left))*(bin_wl_right-bin_wl_left)

        wl_adjust = sum(bin_wl)
        if (wl > wl_adjust) then
            do ind_bin = 1, num_bins
                ratio = wl/wl_adjust
                if (ieee_is_nan(ratio)) then
                    print *, "ERROR: the adjust ratio is nan.", wl_adjust, wl, ratio
                    call stopcm1
                end if
                if (ratio >=1.5) then
                    print *, "ERROR: the adjust ratio exceeds 1.5.", wl_adjust, wl, ratio
                    call stopcm1
                end if
                if (bin_wl(ind_bin)>tolerance_bin_wl) then
                    bin_wl(ind_bin) = bin_wl(ind_bin)*ratio
                    bin_num_conc(ind_bin) = bin_num_conc(ind_bin)*ratio
                end if
            end do
        end if

        do ind_bin = 1, num_bins
            if (bin_wl(ind_bin) < tolerance_bin_wl) then
                cycle
            end if
            bin_num_conc(ind_bin) = max(bin_num_conc(ind_bin), bin_wl(ind_bin)*rho_d/rho_w/((4d0*pi/3d0)*bin_r_wet(ind_bin)**3d0))
            if((bin_r_wet(ind_bin)**3d0-bin_wl(ind_bin)*3d0/4d0/pi/bin_num_conc(ind_bin)*rho_d/rho_w) > max(0.1d0*bin_r_wet(ind_bin),ccn_radius_max)**3d0) then
                bin_num_conc(ind_bin) = bin_wl(ind_bin)*rho_d/rho_w/((4d0*pi/3d0)*(bin_r_wet(ind_bin)**3d0-max(0.1d0*bin_r_wet(ind_bin),ccn_radius_max)**3d0))
            end if
        end do
        if (abs(sum(bin_wl)-wl) > tolerance_wl) then
            print *, "mass imbalance in DSD interpolation"
            call stopcm1
        end if

    end subroutine interpolate_DSD



    subroutine apply_DSD_bdy (hdf, wl_new, temp, p, qpoint, bin_r_dry)
        real (kind=dp), intent(in) :: hdf
        real (kind=sp), dimension(numq), intent(inout) :: qpoint
        real (kind=dp), dimension(num_bins), intent(inout) :: bin_r_dry
        real (kind=dp), intent(in) :: wl_new, temp, p
        real (kind=dp), dimension(num_bins) :: bin_num_conc, bin_wl, bin_r_wet_temp
        real (kind=dp) :: wv, rho_d, wl
        real (kind=dp) :: ratio
        integer :: ind_bin
        
        
        if (if_DSD_lookup .eqv. .true.) then
            bin_num_conc = qpoint(nnc1+num_bins:nnc2)
            bin_wl = qpoint(nnc1:nnc1+num_bins-1)/scale_bin_wl
            wv = qpoint(nqv)
            rho_d = calculate_rho_d (temp, p, wv)
            if (rho_d < 0.1d0) then
                print *, "ERROR: density too small", rho_d, temp, p, wv
                call stopcm1
            end if
            bin_r_wet_temp = bin_r_wet
            call calculate_bin_r_dry(bin_wl, bin_r_dry, bin_num_conc, bin_r_wet_temp, rho_d, wv, wl)
            rho_d = calculate_rho_d (temp, p, wv)
            if (wl_new < tolerance_wl) then
                do ind_bin = 1, num_bins
                    if (bin_r_wet(ind_bin) .le. ccn_radius_max) then
                        bin_r_dry(ind_bin) = bin_r_wet(ind_bin)
                    else
                        bin_r_dry(ind_bin) = 0
                        bin_num_conc(ind_bin) = 0
                    end if
                    bin_wl(ind_bin) = 0d0
                end do
            else
                do ind_bin = 1, num_bins
                    if (bin_wl(ind_bin)<=tolerance_bin_wl) then
                        bin_wl(ind_bin) = 0d0
                        if (bin_r_wet(ind_bin) > ccn_radius_max) then
                            bin_r_dry(ind_bin) = 0d0
                            bin_num_conc(ind_bin) = 0d0
                        else
                            bin_r_dry(ind_bin) = bin_r_wet(ind_bin)
                        end if
                    end if
                end do
                if (sum(bin_wl) > tolerance_wl) then
                    ratio = wl_new/sum(bin_wl)        
                    do ind_bin = 1, num_bins
                        if (bin_wl(ind_bin)>tolerance_bin_wl) then
                            bin_wl(ind_bin) = bin_wl(ind_bin)*ratio
                            bin_r_wet_temp(ind_bin) = ( bin_wl(ind_bin)/bin_num_conc(ind_bin)*rho_d/rho_w/(4d0/3d0*pi)+bin_r_dry(ind_bin)**3d0 )**(1d0/3d0)      
                            if (bin_r_wet_temp(ind_bin) .le. bin_r_wet(1)) then
                                print *, "scale not working:", "ratio:", ratio, "wl_new:", wl_new, bin_r_dry(ind_bin), bin_r_wet(ind_bin), bin_num_conc(ind_bin)
                                call stopcm1
                            end if
                            if (bin_r_wet_temp(ind_bin) .ge. drop_radius_max) then
                                print *, "ERROR: wet radius too large from apply_DSD_bdy!", "ratio:", ratio, "wl_new:", wl_new
                                call stopcm1
                            end if    
                        end if
                    end do
                    if (abs(sum(bin_wl)-wl_new) > tolerance_wl) then
                        print *, "mass imbalance in DSD boundary scaling, sum(bin_wl), wl_new", sum(bin_wl), wl_new
                        call stopcm1
                    end if
                    call update_bins_static(bin_r_dry, bin_r_wet_temp, bin_num_conc, "at apply_DSD_bdy")
                    bin_wl = (4d0*pi/3d0)*bin_num_conc*(bin_r_wet**3d0 - bin_r_dry**3d0) *rho_w/rho_d
                    ! if (abs(ratio-1) > 0.1) then
                    !     if (myid .eq. 0) then
                    !         print *, "myid", myid, "ratio:", ratio, "wl_new", wl_new, "bin_wl", bin_wl, "num conc", bin_num_conc
                    !         call stopcm1
                    !     end if
                    ! end if
                else
                    call interpolate_DSD (hdf, wl_new, wv, temp, p, bin_num_conc, bin_wl)
                end if
            end if
            qpoint(nnc1+num_bins:nnc2) = bin_num_conc
            qpoint(nnc1:nnc1+num_bins-1) = bin_wl*scale_bin_wl
            qpoint(nql1) = wl_new
        end if
    end subroutine apply_DSD_bdy

    subroutine input_turbulence_kernel (turbulence_kernel)
        real (kind=dp), dimension(21,21), intent(out) :: turbulence_kernel
        integer :: i, j
        OPEN(UNIT=104,FILE='data/KERNEL_0500_20000.TXT',STATUS='OLD',FORM='FORMATTED')
            DO i=1,21
                READ(104,*) (turbulence_kernel(i,j), j=1,21)
            END DO
        CLOSE(104)
    end subroutine input_turbulence_kernel

    subroutine input_collision_kernels (rd1, rd2, T0000_1, T2500_1, T5000_1, T0000_2, T2500_2, T5000_2, T0000_3, T2500_3, T5000_3)
        real (kind=dp), dimension(:), intent(out) :: rd1, rd2
        real (kind=dp), dimension(size(rd1),size(rd1)), intent(out) :: T0000_1, T2500_1, T5000_1
        real (kind=dp), dimension(size(rd1),size(rd2)), intent(out) :: T0000_2, T2500_2, T5000_2
        real (kind=dp), dimension(size(rd2),size(rd2)), intent(out) :: T0000_3, T2500_3, T5000_3
        integer :: nd1, nd2, i, j
        nd1 = size(rd1)
        nd2 = size(rd2)
        !-----------------------
        ! INPUT OF 1th TABLE
        !-----------------------
        OPEN(UNIT=100,FILE='data/EFF_1000_1',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd1
                READ(100,*) rd1(i), (T0000_1(i,j), j=1,nd1)
            END DO
        CLOSE(100)

        OPEN(UNIT=101,FILE='data/EFF_750_1',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd1
                READ(101,*) rd1(i), (T2500_1(i,j), j=1,nd1)
            END DO
        CLOSE(101)

        OPEN(UNIT=102,FILE='data/EFF_500_1',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd1
                READ(102,*) rd1(i), (T5000_1(i,j), j=1,nd1)
            END DO
        CLOSE(102)

        !-----------------------
        ! INPUT OF 2nd TABLE
        !-----------------------
        OPEN(UNIT=100,FILE='data/EFF_1000_2',STATUS='OLD',FORM='FORMATTED')
            READ(100,*) (rd2(j), j=1,nd2)
            DO i=1,nd1
                READ(100,*) rd1(i), (T0000_2(i,j), j=1,nd2)
            END DO
        CLOSE(100)

        OPEN(UNIT=101,FILE='data/EFF_750_2',STATUS='OLD',FORM='FORMATTED')
            READ(101,*) (rd2(j), j=1,nd2)
            DO i=1,nd1
                READ(101,*) rd1(i), (T2500_2(i,j), j=1,nd2)
            END DO
        CLOSE(101)

        OPEN(UNIT=102,FILE='data/EFF_500_2',STATUS='OLD',FORM='FORMATTED')
            READ(102,*) (rd2(j), j=1,nd2)
            DO i=1,nd1
                READ(102,*) rd1(i), (T5000_2(i,j), j=1,nd2)
            END DO
        CLOSE(102)

        !-----------------------
        ! INPUT OF 3rd TABLE
        !-----------------------
        OPEN(UNIT=100,FILE='data/EFF_1000_3',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd2
                READ(100,*) rd2(i), (T0000_3(i,j), j=1,nd2)
            END DO
        CLOSE(100)

        OPEN(UNIT=101,FILE='data/EFF_750_3',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd2
                READ(101,*) rd2(i), (T2500_3(i,j), j=1,nd2)
            END DO
        CLOSE(101)

        OPEN(UNIT=102,FILE='data/EFF_500_3',STATUS='OLD',FORM='FORMATTED')
            DO i=1,nd2
                READ(102,*) rd2(i), (T5000_3(i,j), j=1,nd2)
            END DO
        CLOSE(102)

        DO i=1,nd1
            rd1(i) = rd1(i)*1d-6 ! convert micron to m
        END DO
        DO j=1,nd2
            rd2(j) = rd2(j)*1d-6 ! convert micron to m
        END DO

    end subroutine input_collision_kernels

    subroutine preprocess_bins (bin_r_dry, bin_r_wet, bin_num_conc, bin_r_wet_grid, supersat, temp)
        real (kind=dp), dimension(:), intent(in) :: bin_r_wet
        real (kind=dp), dimension(:), intent(out) :: bin_r_dry, bin_num_conc,  bin_r_wet_grid
        real (kind=dp):: supersat, temp
        
        call fill_bin_r_dry (bin_r_wet, bin_r_dry, supersat, temp) ! supersat and temp not used anymore, naive filling
        call fill_bin_num_conc (bin_r_dry, bin_num_conc)
        call make_bin_r_wet_grid(bin_r_wet_grid, bin_r_wet)  
        
    end subroutine preprocess_bins

    subroutine make_bin_r_wet_grid(bin_r_wet_grid, bin_r_wet)
        real (kind=dp), dimension(:), intent(in) :: bin_r_wet
        real (kind=dp), dimension(:), intent(out) :: bin_r_wet_grid
        integer :: ind_bin
        do ind_bin = 1, num_bins + 1
            if (ind_bin == 1) then
                bin_r_wet_grid(ind_bin) = ccn_radius_min
            else if (ind_bin == num_bins + 1) then
                bin_r_wet_grid(num_bins+1) = drop_radius_max
            else
                bin_r_wet_grid(ind_bin) = (bin_r_wet(ind_bin-1)*bin_r_wet(ind_bin))**0.5 
            end if
        end do     
    end subroutine make_bin_r_wet_grid

    subroutine make_bin_r_wet (bin_r_wet)
        real (kind=dp), dimension(:), intent(out) :: bin_r_wet
        integer :: ind_bin
        ! ratio_bins = (drop_radius_max/ccn_radius_min)**(1d0/real((num_bins-1)))
        ! bin_r_wet(1) = ccn_radius_min
        ! do ind_bin = 2, num_bins
        !     bin_r_wet(ind_bin) = bin_r_wet(ind_bin-1) * ratio_bins
        ! end do
        ! bin_r_wet = bin_r_wet/(bin_r_wet(num_bins)/drop_radius_max)
        ! call fill_logspace (ccn_radius_min, drop_radius_max, bin_r_wet)
        real (kind=dp), dimension(size(bin_r_wet)+1) :: bin_points
        real (kind=dp) :: a, b ! left and right values for each bin for convenience
        call fill_logspace (ccn_radius_min, drop_radius_max, bin_points)
        do ind_bin = 1, num_bins
            a = bin_points(ind_bin)
            b = bin_points(ind_bin+1)
            bin_r_wet(ind_bin) = sqrt(a*b) ! geometric mean
        end do
    end subroutine make_bin_r_wet

    subroutine fill_bin_num_conc (bin_r_dry, bin_num_conc)
        real (kind=dp), dimension(:), intent(in) :: bin_r_dry
        real (kind=dp), dimension(:), intent(out) :: bin_num_conc
        integer :: num_bins, ind_bin
        real (kind=dp) :: a, b ! left and right values for each bin for convenience
        num_bins = size(bin_num_conc)
        do ind_bin=1, num_bins
            if (bin_r_dry(ind_bin) > ccn_radius_max) then
                bin_num_conc(ind_bin) = 0
            else if (bin_r_dry(ind_bin) < tolerance_bin_r) then
                bin_num_conc(ind_bin) = 0
            else
                if (ind_bin == 1) then
                    a = ccn_radius_min
                else
                    a = sqrt(bin_r_dry(ind_bin) * bin_r_dry(ind_bin-1))
                end if
                if (ind_bin == num_bins) then
                    b = ccn_radius_max
                else
                    b = sqrt(bin_r_dry(ind_bin) * bin_r_dry(ind_bin+1))
                end if
                bin_num_conc(ind_bin) = calculate_num_conc (a, b)
            end if
        end do
        if (if_fill_left_line) then
            do ind_bin = num_bins, 1, -1
                if (bin_r_dry(ind_bin) .lt. aer_gmds(1)/2 .and. bin_r_dry(ind_bin) .ge. ccn_radius_min/10) then
                    bin_num_conc(ind_bin) = bin_num_conc(ind_bin+1)
                end if
            end do
        end if
    end subroutine fill_bin_num_conc

    subroutine fill_bin_r_dry (bin_r_wet_curr, bin_r_dry, supersat, temp)
        real (kind=dp), dimension(:), intent(out) :: bin_r_dry
        real (kind=dp), dimension(:), intent(in) :: bin_r_wet_curr
        real (kind=dp), intent(in) :: supersat, temp
        integer ind_bin
        do ind_bin = 1, num_bins
            if (bin_r_wet_curr(ind_bin) > ccn_radius_max) then
                bin_r_dry(ind_bin) = 0d0
            else
                bin_r_dry(ind_bin) = bin_r_wet_curr(ind_bin)
                ! The following method is deprecated
                ! bin_r_dry(ind_bin) = calculate_r_dry(bin_r_wet_curr(ind_bin), supersat, temp)
            end if
        end do

    end subroutine fill_bin_r_dry

    function calculate_num_conc (a,b) result(num_conc_tot)
        real(kind=dp), intent(in) :: a,b
        real (kind=dp) :: aer_num_conc, aer_gmd, aer_gsd
        real(kind=dp) :: num_conc_tot
        integer ind_dist
        num_conc_tot = 0
        if (a<b) then
            do ind_dist = 1, size(aer_num_concs)
                aer_num_conc = aer_num_concs(ind_dist)
                aer_gmd = aer_gmds(ind_dist)
                aer_gsd = aer_gsds(ind_dist)
                num_conc_tot = num_conc_tot + aer_num_conc * 0.5*(b*2-a*2)*(calculate_lognorm_prob_density(a*2, aer_gmd, aer_gsd) +  calculate_lognorm_prob_density(b*2, aer_gmd, aer_gsd))
            end do
        end if
    end function calculate_num_conc
    
    function calculate_lognorm_prob_density( x, gm, gsd) result(lognorm_prob_density)
        ! function to calculate probability density of a lognormal distribution given geometric mean and geometric standard deviation
        ! https://en.wikipedia.org/wiki/Log-normal_distribution
        real (kind=dp), intent(in) :: x, gm, gsd
        real (kind=dp) :: sigma, mu, lognorm_prob_density
        ! gm: geometric mean
        ! gsd: geometric standard deviation
        sigma = log (gsd)
        mu = log (gm)
        lognorm_prob_density = 1/(x*sigma*sqrt(2*pi))*exp(-(log(x)-mu)**2/(2*sigma**2))
    end function calculate_lognorm_prob_density

    subroutine fill_linspace(xl, xr, xs)
        ! function to make bins linearly between xl and xr
        real (kind=dp), intent(in) :: xl, xr
        real (kind=dp), dimension(:), intent(out) :: xs
        integer :: i, n ! i is index, n is number of bin points
        n = size(xs)
        if (n == 1) then
            if ( xl /= xr ) then
                print *, "ERROR: n = 1, but xl != xr"
                call stopcm1
            else
                xs = xl
            end if
        else
            do i = 1, n
            xs(i) = (xr - xl)*real(i-1, dp) / real(n-1, dp) + xl
            end do
        end if
        return
    end subroutine fill_linspace

    subroutine fill_logspace(xl, xr, xs)
        ! function to make bin_points between xl and xr with equal space in log scale    
        real (kind=dp), intent(in) :: xl, xr
        real (kind=dp), dimension(:), intent(out) :: xs
        call fill_linspace ( log10(xl), log10(xr), xs)
        xs = 10.d0**xs
        return
    end subroutine fill_logspace

!#####################       Thermal Calculation         ###################

    function calculate_rho_d (temp, p, wv) result (rho_d)
        ! calculate moist air density (rho_d) [kg m-3]
        real (kind=dp), intent(in) :: temp, p, wv
        real (kind=dp) :: rho_d
        rho_d = p/R_d/((1d0 + reps*wv)*temp)
    end function calculate_rho_d

    function calculate_wv (temp, supersat, p) result(wv)
        ! calculate water vapor mixing ratio
        ! supersat: supersaturation
        ! temp: temperature
        real(kind=dp), intent(in) :: supersat, temp, p
        real (kind=dp) :: wv
        real (kind=dp) :: e_s
        e_s = (supersat + 1d0) * calculate_e_s(temp)
        wv = e_s * MwMa /(p-e_s)
    end function calculate_wv

    function calculate_e_s (temp) result(e_s)
        ! calculate the saturation vapor pressure over water for a given temperature.
        ! Uses an empirical fit [Huang et al. 2018, JAMC]
        ! temp: temperature [K]
        real(kind=dp), intent(in) :: temp
        real(kind=dp) :: temp_celsius
        real (kind=dp) :: e_s
        temp_celsius = temp - 273.15
        e_s = exp (34.494 - 4924.99/(temp_celsius + 237.1))/(temp_celsius + 105)**1.57
    end function calculate_e_s

    function calculate_supersat (temp, p, wv) result(supersat)
        ! calculate supersaturation
        real(kind=dp), intent(in) :: temp, p, wv
        real (kind=dp) :: e, e_s, supersat
        if (wv < tolerance_wv) then
            supersat = 0
        else
            e_s = calculate_e_s (temp)
            e = p/(MwMa/wv+1)
            supersat = e/e_s-1
        end if
        ! if (supersat > 0) then
        !     print *, "supersat, temp, p, wv",supersat, temp, p, wv
        ! end if
    end function calculate_supersat

    function calculate_sigma_w (temp) result (sigma_w)
        ! Surface tension of water for a given temperature.
        ! Parameters
            ! temp : float, ambient air temperature [K]
        ! Returns
            ! sigma_w: float, water surface tension coefficient, [J/m^2]
        real (kind=dp), intent(in) :: temp
        real (kind=dp) :: sigma_w
        sigma_w = 0.0761d0 - (1.55d-4)*(temp - 273.15d0)
    end function calculate_sigma_w

    function calculate_equilibrium_r_wet (supersat, temp, r_dry) result (equilibrium_r_wet)
        ! function to calculate equilibrium wet radius with bisection method
        ! input:
        !   supersat: supersaturation
        !   temp: temperature [K]
        !   r_dry: dry radius of core material [m]
        !   aer_kappa: particle hygroscopicity parameter
        ! return:
        !   equilibrium_r_wet: equilibrium wet radius, the wet particle size to reach equilibrium with supersaturation in the environment [m]
        real (kind=dp), intent(in) :: supersat, temp, r_dry
        real (kind=dp) :: equilibrium_r_wet
        integer :: times_iter ! iteration times
        real (kind=dp) :: r_wet_lower, r_wet_upper, r_wet_mid, error_lower, error_upper, error_mid !lower bound, upper bound, mid wet radius and errors in bisection method
        integer :: times_iter_max = 3000 ! maximum iteration times
        real (kind=dp) :: tolerance = 1d-8 ! tolerance to end iteration
        r_wet_lower = r_dry ! initial lower bound: dry radius (core material)
        r_wet_upper = calculate_r_crit (temp, r_dry, aer_kappa) ! initial upper bound

        times_iter = 0
        do while (times_iter < times_iter_max)
            times_iter = times_iter + 1
            r_wet_mid = (r_wet_lower + r_wet_upper)/2d0
            error_mid = calculate_equilibrium_supersat (r_wet_mid, r_dry, temp, aer_kappa) - supersat
            error_lower = calculate_equilibrium_supersat (r_wet_lower, r_dry, temp, aer_kappa) - supersat
            error_upper = calculate_equilibrium_supersat (r_wet_upper, r_dry, temp, aer_kappa) - supersat

            if (abs(error_mid) < tolerance) then 
                equilibrium_r_wet = r_wet_mid
                return
            else if (error_mid * error_lower < 0) then
                r_wet_upper = r_wet_mid
            else if (error_mid * error_upper < 0) then
                r_wet_lower = r_wet_mid
            else
                print *, "Bisection failed to calculate equilibrium wet radius: same sigh for upper and lower wet radius bound."
                print *, "Dry radius is: ", r_dry, " r_crit: ", calculate_r_crit (temp, r_dry, aer_kappa), " temp: ", temp, " kappa: ", aer_kappa, "supersat:", supersat
                print *, "upper supersat:", calculate_equilibrium_supersat (r_wet_upper, r_dry, temp, aer_kappa), "lower supersat:", calculate_equilibrium_supersat (r_wet_lower, r_dry, temp, aer_kappa)
                call stopcm1
            end if
        end do

        print *, "Bisection failed to calculate equilibrium wet radius: iteration number exceeds maximum setting."
        print *, "supersat",supersat, "temp",temp, "r_dry",r_dry
        call stopcm1
    end function calculate_equilibrium_r_wet

    function calculate_equilibrium_r_dry (supersat, temp, r_wet) result (equilibrium_r_dry)
        ! function to calculate equilibrium dry radius given wet radius (to update r_dry for newly added bin)
        ! with bisection method
        ! input:
        !   supersat: supersaturation
        !   temp: temperature [K]
        !   r_wet: wet radius [m]
        !   aer_kappa: particle hygroscopicity parameter
        ! return:
        !   equilibrium_r_dry: equilibrium dry radius, the core material radius to enable the drop (particle) to reach equilibrium with supersaturation in the environment [m]
        real (kind=dp), intent(in) :: supersat, temp, r_wet
        real (kind=dp) :: equilibrium_r_dry
        integer :: times_iter ! iteration times
        real (kind=dp) :: r_dry_lower, r_dry_upper, r_dry_mid, error_lower, error_upper, error_mid !lower bound, upper bound, mid wet radius and errors in bisection method
        integer :: times_iter_max = 3000 ! maximum iteration times
        real (kind=dp) :: tolerance = 1d-8 ! tolerance to end iteration

        r_dry_lower = 0 ! initial lower bound: dry radius (core material)
        r_dry_upper = r_wet ! initial upper bound: surface temp

        times_iter = 0
        do while (times_iter < times_iter_max)
            times_iter = times_iter + 1
            r_dry_mid = (r_dry_lower + r_dry_upper)/2d0
            error_mid = calculate_equilibrium_supersat (r_wet, r_dry_mid, temp, aer_kappa) - supersat
            error_lower = calculate_equilibrium_supersat (r_wet, r_dry_lower, temp, aer_kappa) - supersat
            error_upper = calculate_equilibrium_supersat (r_wet, r_dry_upper, temp, aer_kappa) - supersat

            if (abs(error_mid) < tolerance) then 
                equilibrium_r_dry = r_dry_mid
                return
            else if (error_mid * error_lower < 0) then
                r_dry_upper = r_dry_mid
            else if (error_mid * error_upper < 0) then
                r_dry_lower = r_dry_mid
            else
                print *, "Bisection failed to calculate equilibrium dry radius: same sigh for upper and lower dry radius bound."
                print *, "iteration time:", times_iter, "r_dry_lower:", r_dry_lower, "r_dry_upper", r_dry_upper
                print *, "temp:", temp, "supersat:", supersat, "aer_kappa:", aer_kappa, "r_dry_lower:", r_dry_lower, "r_dry_upper:", r_dry_upper
                print *, "upper supersat:", calculate_equilibrium_supersat (r_wet, r_dry_upper, temp, aer_kappa), "lower supersat:", calculate_equilibrium_supersat (r_wet, r_dry_lower, temp, aer_kappa)
                ! print *, calculate_equilibrium_r_wet (supersat, temp, 5d-6, aer_kappa)
                call stopcm1
            end if
        end do
        if (equilibrium_r_dry > r_wet .or. equilibrium_r_dry < 0 .or. ieee_is_nan(equilibrium_r_dry)) then
            print *, "wrong equilibrium_r_dry value in calculate_equilibrium_r_dry:", equilibrium_r_dry
            print *, "r_wet:", r_wet, "temp:", temp, "supersat:", supersat
            call stopcm1
        end if


        print *, "Bisection failed to calculate equilibrium wet radius: iteration number exceeds maximum setting."
        call stopcm1
    end function calculate_equilibrium_r_dry

    function calculate_equilibrium_supersat (r_wet, r_dry, temp, aer_kappa) result(equilibrium_supersat)
        ! kappa-Kohler theory equilibrium supersaturation over aerosol.
        !
        !    Calculates the equilibrium supersaturation (relative to 100% RH) over an aerosol particle of given dry/wet radius and of specified hygroscopicity bathed in gas at a particular temperature, following [PK2007].
        !
        !    Parameters
        !    ----------
        !    r_wet : float
        !        droplet radius with water [m]
        !    r_dry : float
        !        dry particle radius, i.e. core material [m]
        !    temp : float
        !        ambient air temperature, [K]
        !    aer_kappa: float
        !        particle hygroscopicity parameter
        !
        !    Returns
        !    -------
        !    equilibrium_supersat: float, equilibrium supersaturation
        !
        !    References
        !    ----------
        !
        !    [PK2007] Petters, M. D., and S. M. Kreidenweis. "A Single Parameter
        !        Representation of Hygroscopic Growth and Cloud Condensation Nucleus
        !        Activity." Atmospheric Chemistry and Physics 7.8 (2007): 1961-1971
        real (kind=dp), intent(in) :: r_wet, r_dry, temp, aer_kappa
        real (kind=dp) :: A, B ! intermediate value for calculation
        real (kind=dp) :: sigma_w, equilibrium_supersat ! water surface tension and calculated equilibrium supersaturation
        sigma_w = calculate_sigma_w (temp) ! calculate water surface tension coefficient
        A = (2d0*M_w*sigma_w)/R/temp/rho_w/r_wet
        B = (r_wet**3d0 - r_dry**3d0)/(r_wet**3d0 - (r_dry**3d0)*(1d0-aer_kappa))
        equilibrium_supersat = exp(A)*B - 1d0
    end function calculate_equilibrium_supersat

    function calculate_r_crit (temp, r_dry, aer_kappa) result(r_crit)
        ! This function computes the critical size.
        ! For a given aerosol particle, it will compute the analytical critical point according to the approximated kappa-Kohler equation.
        ! Parameters
            ! temp: float, ambient air temperature [K]
            ! r_dry : float, dry particle radius [m]
            ! aer_kappa : float, particle hygroscopicity parameter
        !Returns
            !r_crit: the critical radius [m] of the aerosol droplet.
        real (kind=dp), intent(in) :: temp, r_dry, aer_kappa
        real (kind=dp) :: A, r_crit
        A = 2d0 * M_w * calculate_sigma_w (temp) / (R * temp * rho_w)
        r_crit = sqrt (3d0*aer_kappa*(r_dry**3d0)/A)
    end function calculate_r_crit

    function calculate_supersat_crit (temp, r_dry, aer_kappa) result (supersat_crit)
        ! This function computes the critical supersaturation.
        ! For a given aerosol particle, it will compute the analytical critical point according to the approximated kappa-Kohler equation.
        ! Parameters
            ! temp: float, ambient air temperature [K]
            ! r_dry : float, dry particle radius [m]
            ! aer_kappa : float, particle hygroscopicity parameter
        !Returns
            !supersat_crit: the critical supersaturation of the aerosol droplet.
        real (kind=dp), intent(in) :: temp, r_dry, aer_kappa
        real (kind=dp) :: A, supersat_crit
        A = 2d0 * M_w * calculate_sigma_w (temp) / (R * temp * rho_w)
        supersat_crit = sqrt( 4d0*(A**3d0)/(27d0*aer_kappa*(r_dry**3d0)))
    end function calculate_supersat_crit

    function calculate_thermal_conductivity_air(temp, rho, r_wet) result(thermal_conductivity_air)
        ! Thermal conductivity of air, modified for non-continuum effects.
        !
        !    Modification to account for non-continuum effects (small aerosol/droplet
        !    size) yields the equation
        !
        !    Parameters
        !    ----------
        !    temp : float
        !        ambient air temperature, K
        !    rho : float
        !        ambient air density (excluding water vapor), kg/m^3
        !    r_wet : float
        !        droplet radius, m
        !
        !    Returns
        !    -------
        !    thermal_conductivity_air: float, J/m/s/K
        !
        !    References
        !    ----------
        !
        !    [SP2006] Seinfeld, John H, and Spyros N Pandis. Atmospheric Chemistry
        !       and Physics: From Air Pollution to Climate Change. Vol. 2nd. Wiley, 2006.
        real(kind=dp), intent(in) :: temp, rho, r_wet
        real(kind=dp) :: thermal_conductivity_air, numerator, denominator
        numerator = (1d-3)*(4.39d0 + 0.071d0*temp)
        denominator = 1d0 + (numerator/coeff_thermal_conductivity_air/r_wet/rho/C_p)*sqrt(2.d0*pi*M_a/R/temp)
        thermal_conductivity_air = numerator/denominator
    end function calculate_thermal_conductivity_air

    function calculate_diffusivity_water_vapor(temp, r_wet, p) result(diffusivity_water_vapor)
        ! Diffusivity of water vapor in air, modified for non-continuum effects.
        !
        !    The diffusivity of water vapor in air as a function of temperature and pressure
        !    is given by [SP2006].
        !
        !    Parameters
        !    ----------
        !    temp: float
        !        ambient temperature of air surrounding droplets, K
        !    r_wet: float
        !        radius of aerosol/droplet, m
        !    p : float
        !        ambient pressure of surrounding air, Pa
        !
        !    Returns
        !    -------
        !    diffusivity_water_vapor: float, m^2/s
        !
        !    References
        !    ----------
        !
        !    [SP2006] Seinfeld, John H, and Spyros N Pandis. Atmospheric Chemistry
        !       and Physics: From Air Pollution to Climate Change. Vol. 2nd. Wiley, 2006.

        real(kind=dp), intent(in) :: temp, r_wet, p
        real(kind=dp) :: diffusivity_water_vapor, numerator, denominator
        numerator = (1.d-4)*(0.211d0/(p*1.01325d-5))*((temp/273.d0)**1.94d0)
        denominator = 1.d0 + (numerator/coeff_condensation_diffusivity/r_wet)*sqrt(2.d0*pi*M_w/R/temp)
        diffusivity_water_vapor = numerator/denominator
    end function calculate_diffusivity_water_vapor

    subroutine extrapolate_endpoints(x, y, n)
        implicit none
        integer, intent(in) :: n
        real(kind=dp), intent(in)    :: x(n)
        real(kind=sp), intent(inout) :: y(n)

        ! Local variables for slope calculations
        real(kind=dp) :: slope_low, slope_high

        if (n < 4) then
            print *, "Error: need at least 4 points for extrapolation."
            call stopcm1
        end if

        ! ---- Lower boundary extrapolation ----
        slope_low = (y(3) - y(2)) / (x(3) - x(2))
        y(1) = y(2) - (x(2) - x(1)) * slope_low

        ! ---- Upper boundary extrapolation ----
        slope_high = (y(n-1) - y(n-2)) / (x(n-1) - x(n-2))
        y(n) = y(n-1) + (x(n) - x(n-1)) * slope_high

    end subroutine extrapolate_endpoints


end module dcpm
