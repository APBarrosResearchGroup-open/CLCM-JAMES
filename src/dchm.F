module dchm
use input, only:    ib, ie, jb, je, kb, ke, ibm, iem, jbm, jem, kbm, kem, ibr, ier, jbr, jer, ibl, iel, jbl, jel, ni, nj, &
                    sp, dp, myid, ierr, numprocs, myi, myj, isfc, irst, days_spinup_dchm, shadingoption, scale_moisture
use constants, only : eps
use cm1libs , only : rslf
use ieee_arithmetic

#ifdef MPI
    use mpi
#endif

implicit none
integer, parameter :: Ntotlayer=12 !! hardcoded: Ntotlayer = Nlayer
integer, parameter :: Ntotv = 2 !! hardcoded: Ntotv = Nv
private

! ====================================== Model Setup  ======================================
integer, parameter :: yrn0 = 2002                           ! start year
integer, parameter :: yrn2 = 2008                           ! end year
integer, parameter :: shft = 117418
integer, parameter :: N_LC = 1                              ! number of vegetation types in each grid
integer, parameter :: Nh = INT(N_LC+2)                      ! number of heights in each grid (= N_LC + 2)
integer, parameter :: Nv = INT((2**N_LC)-1+1)               ! number of distinct vegetation cover patch in each grid (= 2^(N_LC)-1+1)
integer, parameter :: Nlayer = 12                         ! number of soil layers
! integer, parameter :: Ntotv = Nv

integer :: ispin                                            ! flag of spinup 1 ==>> spinup
integer :: itsp                                             ! timestep counter
logical :: if_newday                                        ! if this is a new day
real :: dts                                                 ! timestep

real(kind=dp), parameter :: hcf = 0.5                       ! thickness of capillary fringe [m]

real(kind=dp), dimension(:), allocatable :: dxgrid          ! grid dimension [m]
integer :: Ngrid
! ======================================================================================================


! ====================================== Subgrid Specification  ======================================
! Altitudes of each level aboveground [m] (ascending heights)
! (variables with height dimension are calculated following this order)
real(kind=dp), dimension(:,:), allocatable :: harray

! The index of where input forcing is in harray
integer, dimension(:), allocatable :: hfarray

! Land cover types in each grid sorted by ascending heights
! (if same height: sort by ascending # of shading above)
! (variables with dimension Nv are calculated following this order)
CHARACTER(LEN=10), dimension(:,:), allocatable :: varray 
! The index of where bare soil is in varray
integer, PARAMETER :: vbs = 1 

! Height of each vegetation type (the correspoing index in harray)
integer, dimension(:,:), allocatable :: hvarray

! Shading relationship at each height for each shaded vegetation type
! If unshaded, number is 0; If shaded, the numbers represent which 
! vegetation provide the shade (i.e., the index in varray) 
integer, dimension(:,:,:), allocatable :: shading

! Specify which grid has groundwater layer
! {0: groundwater present in soil column; 1: with groundwater layer}

integer, dimension(:), allocatable :: gwlayer
! Specify which vegetation in a subgrid has access to groundwater
! {0: groundwater inaccessible; 1: vegetation roots reach groundwater}
integer, dimension(:,:), allocatable :: gwaccess

! Specify which vegetation performs C3 photosynthesis 
! {0: C4 photosynthesis; 1: C3 photosynthesis}
integer, dimension(:,:), allocatable :: photoc3
! ======================================================================================================

! ====================================== Vegetation Specification ======================================
! Specify shape of initial root profile
! {'AK': more deep roots [Amenu & Kumar, 2008]; 'Zeng': more shallow roots [Zeng, 2001]} 
CHARACTER(LEN=4), PROTECTED :: RootType = 'Zeng' 
real, dimension(:,:), allocatable :: z50, z95, maxrootdepth
real, PARAMETER :: Grmax=0.1   ! maximum daily growth rate (as a fraction of previous root fraction)
real, PARAMETER :: Srmax=0.0   ! maximum daily reduction rate (as a fraction of previous root fraction)
real, PARAMETER :: wcrit=0.09  ! soil moisture threshold (above which roots grow, and shrink otherwise)

! Stomatal sensitivity to VPD and Temperature
real, dimension(:,:), allocatable :: vdelta, vmm, T_b, T_high
! ======================================================================================================

! Constants
real (kind=dp), PARAMETER :: rw = 1000.0         ! water density [kg/m3]
real (kind=dp), PARAMETER :: sig = 5.6696e-08    ! Stefan-Boltzman constat
real (kind=dp), PARAMETER :: emg_soil = 0.985         ! emissivity of soil
real (kind=dp), PARAMETER :: Ti = 273.15         ! Freezing temperature
real (kind=dp), PARAMETER :: CI = 2.09E+06        ! Heat capacity of snow
real (kind=dp), PARAMETER :: CW=4.18E+06          ! Heat capacity of water
real (kind=dp), PARAMETER :: LHM=3.34E+05         ! Latent heat for melting
real (kind=dp), PARAMETER :: vh = 2.5e+06         ! Latent heat for evaporation
real (kind=dp), parameter :: Cp_air = 1005.70001  ! Air heat capacity
! ======================================================================================================

! Soil Parameters
real (kind=dp), dimension(:), allocatable :: area
integer, dimension(:), allocatable :: soiltext
real (kind=dp), dimension(:), allocatable :: WD, sf, lambda, wscr, wswp, power, WSAT
real (kind=dp), dimension(:,:), allocatable :: WSATi
real (kind=dp), dimension(:), allocatable :: C, TD
real (kind=dp), dimension(:,:), allocatable :: dsoil

! Vegetation Parameters
integer, dimension(:,:), allocatable :: lctype
real (kind=dp), dimension(:,:), allocatable :: rcmin, critleafwpot, zref
real (kind=dp), dimension(:,:), allocatable :: vlai, fpar, cv, cvsep
real (kind=dp), dimension(:), allocatable :: zrefl, zrefs, zreff
real (kind=dp), dimension(:), allocatable :: rzsoilwpot

! Atmospheric Forcing
real (kind=dp), dimension(:), allocatable :: pr0, pr, psn
real (kind=dp), dimension(:), allocatable :: pnots, qs, issrd, istrd
real (kind=dp), dimension(:), allocatable :: u0, uh, u, vsgl
real (kind=dp), dimension(:,:), allocatable :: tair, qsat, vpd
real (kind=dp), PARAMETER :: hu = 2 ! wind input is at height 2m

! Ancillary Variables
real (kind=dp), dimension(:), allocatable :: a0, cossza, wtd, albdeff, emisseff

! Heat Transfer Variables
real (kind=dp), dimension(:), allocatable :: aircoefs, aircoefl
real (kind=dp), dimension(:), allocatable :: emgs
real (kind=dp), dimension(:,:), allocatable :: FO, FLIN, FOT, FSM
real (kind=dp), dimension(:,:,:), allocatable :: Swi
real (kind=dp), dimension(:), allocatable :: FOavg, FSMavg
real (kind=dp), dimension(:), allocatable :: gflux, hflux, sflux, NR, esfcbgt, etotbgt
real (kind=dp), dimension(:,:,:), allocatable :: soilT, soilTp
real (kind=dp), dimension(:,:), allocatable :: soilTavg, soilTavgp
real (kind=dp), dimension(:), allocatable :: gswnettot, glwnettot

! Skin water reservoir and fluxes
real (kind=dp), dimension(:,:), allocatable :: wskmax, wsk, wskp
real (kind=dp), dimension(:,:), allocatable :: intercep
real (kind=dp), dimension(:), allocatable :: wsktotal, interceptotal

! Snow Routine Variables
real (kind=dp), dimension(:,:), allocatable :: dsn, csn, esn, SLW, psm, CM, TDM
real (kind=dp), dimension(:), allocatable :: pres, prespr

! Vegetation Variables
real (kind=dp), dimension(:,:), allocatable :: rc, rootwpot, plwpot
real (kind=dp), dimension(:,:,:), allocatable :: rootfrc, rootup, dailyrootup, dailyrutot
real (kind=dp), dimension(:,:), allocatable :: dailyswtot
real (kind=dp), dimension(:,:), allocatable :: totrootup

! Evapotranspiration Fluxes
real (kind=dp), dimension(:), allocatable :: ebare, wsksoil
real (kind=dp), dimension(:,:), allocatable :: esk, eskprev, ev, trans
real (kind=dp), dimension(:), allocatable :: etotal

! Soil Water Fluxes
real (kind=dp), dimension(:), allocatable :: ptot, filtot, filflux, PERCR
real (kind=dp), dimension(:), allocatable :: baseflow, exfilt, excess, runf, evpond
real (kind=dp), dimension(:,:), allocatable :: soilwpot, soilW, soilWp
real (kind=dp), dimension(:,:), allocatable :: dailysoilW

! Photosynthesis
real (kind=dp), dimension(:,:), allocatable :: wdr, gpp, funcleafwpot

! CLCM 2D parameters to cummulate between turb.F
real, dimension(:,:), allocatable :: lai2d, sai2d, cv2d, fpar2d, hcanopytop2d, hcanopybot2d, hcanopy2d
integer, dimension(:,:), allocatable :: soilC, soilV, soilC_dchm, soilV_dchm
integer, parameter :: soiltype_default=1
integer, parameter :: lctype_default=4
real (kind=sp), dimension(:,:), allocatable :: ts2m_dchm, qv2m_dchm, s2m_dchm, tscanopy_dchm, qvcanopy_dchm, scanopy_dchm,tdiff_dchm, lh_dchm, gflux_dchm,gswnetsfc_dchm,glwnetsfc_dchm,gswnettot_dchm, glwnettot_dchm,esfcbgt_dchm,etotbgt_dchm
real (kind=sp), dimension(:,:,:), allocatable :: dsoil_dchm
real (kind=sp), dimension(:,:), allocatable :: runoff_sfc_dchm, flow_base_dchm, percolation_dchm
real (kind=sp), dimension(:,:), allocatable :: trans_dchm, ebare_dchm
real (kind=sp) :: rtime_daystart = 0
real (kind=sp) :: rtime_prevday = 0


public :: Ntotlayer, Ntotv, Nv, dsoil_dchm, soilC, soilV, soilC_dchm, soilV_dchm, lai2d, sai2d, cv2d, hcanopytop2d, hcanopybot2d, hcanopy2d, fpar2d
public :: dchm_init, dchm_restart_read
public :: ts2m_dchm, qv2m_dchm, s2m_dchm, tscanopy_dchm, qvcanopy_dchm, scanopy_dchm, tdiff_dchm, lh_dchm, gflux_dchm,gswnetsfc_dchm,glwnetsfc_dchm, gswnettot_dchm, glwnettot_dchm, esfcbgt_dchm,etotbgt_dchm,runoff_sfc_dchm, flow_base_dchm, percolation_dchm
public :: trans_dchm, ebare_dchm

contains

subroutine dchm_init(albd,mavail,sfem,znt,tsk,soilT_dchm, soilW_dchm, xf, yf)
    real, intent(in), dimension(ibl:iel,jbl:jel) :: albd,sfem
    real, intent(inout), dimension(ib:ie,jb:je) :: znt, tsk, mavail
    real, dimension(ibl:iel,jbl:jel,Ntotlayer), intent(inout) :: soilW_dchm
    real, dimension(ibl:iel,jbl:jel,Ntotlayer,Nv), intent(inout) :: soilT_dchm
    real, intent(in), dimension(ib:ie+1) :: xf
    real, intent(in), dimension(jb:je+1) :: yf
    integer :: i, j, i2d, j2d, v
    real :: hcanopy, sai, hcanopytop, hcanopybot
    integer :: ii, vv, cvi
    real :: cvb, evi

    Ngrid = ni*nj

    allocate(harray(Ngrid, Nh))
    allocate(hfarray(Ngrid))
    allocate(varray(Ngrid, Nv))
    allocate(hvarray(Ngrid,Nv))
    allocate(dxgrid(Ngrid))
    allocate(shading(Ngrid, Nh, Nv))
    allocate(gwlayer(Ngrid)); gwlayer = 0
    allocate(gwaccess(Ngrid,Nv-1)); gwaccess = 0
    allocate(photoc3(Ngrid, Nv-1)); photoc3 = 0
    allocate(z50(Ngrid, Nv-1)); z50 = 0.75
    allocate(z95(Ngrid, Nv-1)); z95 = 2.5
    allocate(maxrootdepth(Ngrid, Nv-1)); maxrootdepth = 3
    allocate(vdelta(Ngrid, Nv-1)); vdelta = 0.2
    allocate(vmm(Ngrid, Nv-1)); vmm = 0.6
    allocate(T_b(Ngrid, Nv-1)); T_b = 18.35
    allocate(T_high(Ngrid, Nv-1)); T_high = 40.0

    ! Soil Parameters
    allocate(area(Ngrid)); area = 0
    allocate(soiltext(Ngrid)); soiltext = 0
    allocate(WD(Ngrid)); WD = 0
    allocate(sf(Ngrid)); sf = 0
    allocate(lambda(Ngrid)); lambda = 0
    allocate(wscr(Ngrid)); wscr = 0
    allocate(wswp(Ngrid)); wswp = 0
    allocate(power(Ngrid)); power = 0
    allocate(WSAT(Ngrid)); WSAT = 0
    allocate(WSATi(Ngrid,Nlayer)); WSATi = 0
    allocate(C(Ngrid)); C=0
    allocate(TD(Ngrid)); TD=0
    allocate(dsoil(Ngrid,Nlayer)); dsoil = 0

    ! Vegetation Parameters
    allocate(lctype(Ngrid,Nv)); lctype = 0
    allocate(rcmin(Ngrid,Nv)); rcmin = 0
    allocate(critleafwpot(Ngrid,Nv)); critleafwpot = 0
    allocate(zref(Ngrid,Nv)); zref = 0
    allocate(vlai(Ngrid,Nv)); vlai = 0
    allocate(fpar(Ngrid,Nv)); fpar = 0
    allocate(cv(Ngrid,Nv)); cv = 0
    allocate(cvsep(Ngrid,Nv)); cvsep = 0
    allocate(zrefl(Ngrid)); zrefl = 0
    allocate(zrefs(Ngrid)); zrefs = 0
    allocate(zreff(Ngrid)); zreff = 0
    allocate(rzsoilwpot(Ngrid)); zreff = 0

    ! Stmospheric Forcing
    allocate(pr0(Ngrid)); pr0 = 0
    allocate(pr(Ngrid)); pr = 0
    allocate(psn(Ngrid)); psn = 0
    allocate(pnots(Ngrid)); pnots = 0
    allocate(qs(Ngrid)); qs = 0
    allocate(issrd(Ngrid)); issrd = 0
    allocate(istrd(Ngrid)); istrd = 0
    allocate(u0(Ngrid)); u0 = 0
    allocate(uh(Ngrid)); uh = 0
    allocate(u(Ngrid)); u = 0
    allocate(vsgl(Ngrid)); vsgl = 0
    allocate(tair(Ngrid,Nh)); tair = 0
    allocate(qsat(Ngrid,Nh)); qsat = 0
    allocate(vpd(Ngrid,Nh)); vpd = 0

    allocate(a0(Ngrid)); a0 = 0
    allocate(albdeff(Ngrid)); albdeff = 0
    allocate(emisseff(Ngrid)); emisseff = 0
    allocate(cossza(Ngrid)); cossza = 0
    allocate(wtd(Ngrid)); wtd = 0

    ! Heat Transfer Variables
    allocate(aircoefs(Ngrid)); aircoefs = 0
    allocate(aircoefl(Ngrid)); aircoefl = 0
    allocate(emgs(Ngrid)); emgs = 0
    allocate(FO(Ngrid,Nv)); FO = 0
    allocate(FLIN(Ngrid,Nv)); FLIN = 0
    allocate(FOT(Ngrid,Nv)); FOT = 0
    allocate(FSM(Ngrid,Nv)); FSM = 0
    allocate(SWi(Ngrid, Nh, Nv)); SWi = 0
    allocate(FOavg(Ngrid)); FOavg = 0
    allocate(FSMavg(Ngrid)); FSMavg = 0
    allocate(gflux(Ngrid)); gflux = 0
    allocate(hflux(Ngrid)); hflux = 0
    allocate(sflux(Ngrid)); sflux = 0
    allocate(esfcbgt(Ngrid)); esfcbgt = 0
    allocate(etotbgt(Ngrid)); etotbgt = 0
    allocate(NR(Ngrid)); NR = 0
    allocate(soilT(Ngrid,Nlayer,Nv)); soilT = 0
    allocate(soilTp(Ngrid,Nlayer,Nv)); soilTp = 0
    allocate(soilTavg(Ngrid, Nlayer)); soilTavg = 0
    allocate(soilTavgp(Ngrid, Nlayer)); soilTavgp = 0
    allocate(gswnettot(Ngrid)); gswnettot = 0
    allocate(glwnettot(Ngrid)); glwnettot = 0
    ! Skin water reservoir and flux
    allocate(wskmax(Ngrid, Nv)); wskmax = 0
    allocate(wsk(Ngrid, Nv)); wsk = 0
    allocate(wskp(Ngrid, Nv)); wskp = 0
    allocate(intercep(Ngrid, Nv)); intercep = 0
    allocate(wsktotal(Ngrid)); wsktotal = 0
    allocate(interceptotal(Ngrid)); interceptotal = 0

    ! Snow Routine Variables
    allocate(dsn(Ngrid, Nv)); dsn = 0
    allocate(csn(Ngrid, Nv)); csn = 0
    allocate(esn(Ngrid, Nv)); esn = 0
    allocate(SLW(Ngrid, Nv)); SLW = 0
    allocate(psm(Ngrid, Nv)); psm = 0
    allocate(CM(Ngrid, Nv)); CM = 0
    allocate(TDM(Ngrid, Nv)); TDM = 0
    allocate(pres(Ngrid)); pres = 0
    allocate(prespr(Ngrid)); prespr = 0

    ! Vegetation Variables
    allocate(rc(Ngrid, Nv-1)); rc = 0
    allocate(rootwpot(Ngrid, Nv-1)); rootwpot = 0
    allocate(plwpot(Ngrid, Nv-1)); plwpot = 0
    allocate(rootfrc(Ngrid, Nlayer, Nv-1)); rootfrc = 0
    allocate(rootup(Ngrid, Nlayer, Nv-1)); rootup = 0
    allocate(dailyrootup(Ngrid, Nlayer, Nv-1)); dailyrootup = 0
    allocate(dailyrutot(Ngrid, Nlayer, Nv-1)); dailyrutot = 0
    allocate(dailyswtot(Ngrid, Nlayer)); dailyswtot = 0
    allocate(totrootup(Ngrid, Nlayer)); totrootup = 0

    ! Evapotranspiration Fluxes
    allocate(ebare(Ngrid)); ebare = 0
    allocate(wsksoil(Ngrid)); wsksoil = 0
    allocate(esk(Ngrid, Nv-1)); esk = 0
    allocate(eskprev(Ngrid, Nv-1)); eskprev = 0
    allocate(ev(Ngrid, Nv-1)); ev = 0
    allocate(trans(Ngrid, Nv-1)); trans = 0
    allocate(etotal(Ngrid)); etotal = 0

    ! Soil Water Fluxes
    allocate(ptot(Ngrid)); ptot = 0
    allocate(filtot(Ngrid)); filtot = 0
    allocate(filflux(Ngrid)); filflux = 0
    allocate(PERCR(Ngrid)); PERCR = 0
    allocate(baseflow(Ngrid)); baseflow = 0
    allocate(exfilt(Ngrid)); exfilt = 0
    allocate(excess(Ngrid)); excess = 0
    allocate(runf(Ngrid)); runf = 0
    allocate(evpond(Ngrid)); evpond = 0
    allocate(soilwpot(Ngrid,Nlayer)); soilwpot = 0
    allocate(soilW(Ngrid,Nlayer)); soilW = 0
    allocate(soilWp(Ngrid,Nlayer)); soilWp = 0
    allocate(dailysoilW(Ngrid, Nlayer)); dailysoilW = 0

    ! Photosynthesis
    allocate(wdr(Ngrid, Nv-1)); wdr = 0
    allocate(gpp(Ngrid, Nv-1)); gpp = 0
    allocate(funcleafwpot(Ngrid, Nv-1)); funcleafwpot = 0

    ! CLCM related variables
    allocate(runoff_sfc_dchm(ib:ie,jb:je))
    allocate(dsoil_dchm(ib:ie,jb:je,Ntotlayer))
    allocate(flow_base_dchm(ib:ie,jb:je))
    allocate(percolation_dchm(ib:ie,jb:je))
    allocate(trans_dchm(ib:ie,jb:je))
    allocate(ebare_dchm(ib:ie,jb:je))
    allocate(ts2m_dchm(ib:ie,jb:je))
    allocate(qv2m_dchm(ib:ie,jb:je))
    allocate(s2m_dchm(ib:ie,jb:je))
    allocate(tscanopy_dchm(ib:ie,jb:je))
    allocate(qvcanopy_dchm(ib:ie,jb:je))
    allocate(scanopy_dchm(ib:ie,jb:je))
    allocate(lh_dchm(ib:ie,jb:je))
    allocate(gflux_dchm(ib:ie,jb:je))
    allocate(esfcbgt_dchm(ib:ie,jb:je))
    allocate(etotbgt_dchm(ib:ie,jb:je))
    allocate(gswnetsfc_dchm(ib:ie,jb:je))
    allocate(glwnetsfc_dchm(ib:ie,jb:je))
    allocate(gswnettot_dchm(ib:ie,jb:je))
    allocate(glwnettot_dchm(ib:ie,jb:je))
    allocate(tdiff_dchm(ib:ie,jb:je))

    if (isfc .eq. 2023) then
        allocate(soilC(ib:ie,jb:je))
        allocate(soilV(ib:ie,jb:je))
        allocate(soilC_dchm(ib:ie,jb:je))
        allocate(soilV_dchm(ib:ie,jb:je))
        allocate(lai2d(ib:ie,jb:je))
        allocate(sai2d(ib:ie,jb:je))
        allocate(cv2d(ib:ie,jb:je))
        allocate(fpar2d(ib:ie,jb:je))
        call read2Dinteger(soilC(ib:ie,jb:je), "soilC     ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je)
        call read2Dinteger(soilV(ib:ie,jb:je), "soilV     ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je)
        call read2D(lai2d(ib:ie,jb:je), "lai       ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je)
        call read2D(cv2d(ib:ie,jb:je), "cv        ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je)
        call read2D(fpar2d(ib:ie,jb:je), "fpar      ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je)
        if (irst.eq.0) then
            call read3D(soilW_dchm(ib:ie,jb:je,1:Ntotlayer), "soilW_dchm ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je,1,Ntotlayer)
            do v = 1, Nv, 1
                call read3D(soilT_dchm(ib:ie,jb:je,1:Ntotlayer,v), "soilT_dchm ", "DATA_IC/IC_                   ", myi, myj, ib, ie, jb, je,1,Ntotlayer)
            end do
        end if  
        ! print *, "range soilT_dchm:", minval(soilT_dchm), maxval(soilT_dchm)
        ! print *, "range lai:", minval(lai2d(1:ni,1:nj)), maxval(lai2d(1:ni,1:nj))
        ! print *, "range cv:", minval(cv2d(1:ni,1:nj)), maxval(cv2d(1:ni,1:nj))
        ! print *, "range fpar:", minval(fpar2d(1:ni,1:nj)), maxval(fpar2d(1:ni,1:nj))
    end if

    !*************** SOIL PROPERTIES *************** 
    ! Soil depth [dsoil] [m]
    ! Soil texture [soiltext] [-]
    ! Saturated hydraulic conductivity [WD] [mm/s]
    ! Suction front head [sf] [mm]
    ! Pore size distribution [lambda] [-]
    ! field capacity [wscr] [m/m]
    ! wilting point [wswp] [m/m]
    !***********************************************
    allocate(hcanopy2d(ib:ie,jb:je))
    allocate(hcanopytop2d(ib:ie,jb:je))
    allocate(hcanopybot2d(ib:ie,jb:je))
    do i = 1, Ngrid, 1
        call unpack_coord(i, i2d, j2d)
        dxgrid(i) = ((xf(i2d+1)-xf(i2d))*(yf(j2d+1)-yf(j2d)))**0.5
        area(i) = dxgrid(i)*dxgrid(i)
        ! Define soil layer depth
        do j=1,Ntotlayer
            dsoil(i,j)=0.2
            if (j.eq.1) dsoil(i,j)=0.01
            if (j.eq.2) dsoil(i,j)=0.02
            if (j.eq.3) dsoil(i,j)=0.04
            if (j.eq.4) dsoil(i,j)=0.08
            if (j.eq.5) dsoil(i,j)=0.15
            if (j.eq.Ntotlayer) dsoil(i,j)=0.5
        enddo
        
        if (isfc.eq.2023) then
            soiltext(i) = soilC(i2d,j2d)
        else
            soiltext(i) = soiltype_default
        end if
        soilC_dchm(i2d,j2d) = soilC(i2d,j2d)

        ! SET SOIL PROPERTIES      
        if(soiltext(i).eq.1) then
            ! print*,'Soil Texture = Sand'
            WD(i)=23.56*(10.0)*(1./(60.*60.))
            sf(i)=4.95*10. 
            lambda(i)=0.694 
            wscr(i)=0.091  
            wswp(i)=0.033
            WSAT(i)=0.339
        else if(soiltext(i).eq.2) then
            ! print*,'Soil Texture = Loamy sand'
            WD(i)=5.98*(10.)*(1./(60.*60.)) 
            sf(i)=6.13*10. 
            lambda(i)=0.553 
            wscr(i)=0.125  
            wswp(i)=0.055
            WSAT(i)=0.421
        else if(soiltext(i).eq.3) then
            ! print*,'Soil Texture = Sandy Loam'
            WD(i)=2.18*(10.)*(1./(60.*60.)) 
            sf(i)=11.01*10. 
            lambda(i)=0.378 
            wscr(i)=0.207  
            wswp(i)=0.095
            WSAT(i)=0.434
        else if(soiltext(i).eq.4) then
            ! print*,'Soil Texture = Silt Loam'
            WD(i)=0.68*(10.)*(1./(60.*60.)) 
            sf(i)=16.68*10. 
            lambda(i)=0.234 
            wscr(i)=0.333  
            wswp(i)=0.133
            WSAT(i)=0.476
        else if(soiltext(i).eq.5) then
            ! print*,'Soil Texture = Silt'
            WD(i)=0.68*(10.)*(1./(60.*60.)) 
            sf(i)=16.68*10. 
            lambda(i)=0.234 
            wscr(i)=0.333  
            wswp(i)=0.133
            WSAT(i)=0.476
        else if(soiltext(i).eq.6) then
            ! print*,'Soil Texture = Loam'
            WD(i)=1.32*(10.)*(1./(60.*60.)) 
            sf(i)=8.89*10. 
            lambda(i)=0.252 
            wscr(i)=0.27 
            wswp(i)=0.117
            WSAT(i)=0.439
        else if(soiltext(i).eq.7) then
            ! print*,'Soil Texture = Sandy clay loam'
            WD(i)=0.3*(10.)*(1./(60.*60.)) 
            sf(i)=21.85*10.
            lambda(i)=0.319 
            wscr(i)=0.255 
            wswp(i)=0.148
            WSAT(i)=0.404
        else if(soiltext(i).eq.8) then
            ! print*,'Soil Texture = Silty clay loam'
            WD(i)=0.2*(10.)*(1./(60.*60.)) 
            sf(i)=27.3*10.
            lambda(i)=0.177 
            wscr(i)=0.366  
            wswp(i)=0.208
            WSAT(i)=0.464
        else if(soiltext(i).eq.9) then
            ! print*,'Soil Texture = Clay loam'
            WD(i)=0.2*(10.)*(1./(60.*60.)) 
            sf(i)=20.88*10. 
            lambda(i)=0.242 
            wscr(i)=0.318  
            wswp(i)=0.197
            WSAT(i)=0.465
        else if(soiltext(i).eq.10) then
            ! print*,'Soil Texture = Sandy clay'
            WD(i)=0.12*(10.)*(1./(60.*60.)) 
            sf(i)=23.90*10. 
            lambda(i)=0.223 
            wscr(i)=0.339  
            wswp(i)=0.239
            WSAT(i)=0.406
        else if(soiltext(i).eq.11) then
            ! print*,'Soil Texture = Silty clay'
            WD(i)=0.10*(10.)*(1./(60.*60.)) 
            sf(i)=29.22*10. 
            lambda(i)=0.150 
            wscr(i)=0.387 
            wswp(i)=0.250
            WSAT(i)=0.468
        else if(soiltext(i).eq.12) then
            ! print*,'Soil Texture = Clay'
            WD(i)=0.06*(10.)*(1./(60.*60.)) 
            sf(i)=31.63*10. 
            lambda(i)=0.165 
            wscr(i)=0.396  
            wswp(i)=0.272
            WSAT(i)=0.468
        else 
            if (soiltext(i).ne.14) then
                print*,'Soil Texture = NO DATA'            
                call stopcm1    
            end if
            ! place holder for water
            WD(i)=23.56*(10.0)*(1./(60.*60.)) 
            sf(i)=4.95*10. 
            lambda(i)=0.694 
            wscr(i)=0.091  
            wswp(i)=0.033
            WSAT(i)=1
        end if      

        power(i)=nint(3.+(2./lambda(i)))
        WSATi(i,:) = WSAT(i) * dsoil(i,:)

        !************** SOIL THERMAL PROPERTIES **************
        ! Heat capacity [C] [J/m2/K] -> [CDM] for snowmelt
        ! Heat diffusivity [TD] [m2/s] -> [TDM] for snowmelt
        !*****************************************************
        if (soiltext(i).eq.3) then      ! Sandy Loam
            ! [Arkhangelskaya & Lukyashchenko, 2008; Abu-Hamdeh, 2003]
            ! estimated using bulk density≈1.36 𝑔/cm3 and 𝜃≈0.12 m3/m3
            C(i)=1800000.0
            TD(i)=0.0000004
        else 
            ! Clay (used for entire SE US)
            C(i)=2400000.0
            TD(i)=0.00000081
        end if

        !************** LAND COVER **************
        ! Land cover (lctype) 
        !****************************************
        v = 2
        if (isfc.eq.2023) then
            if ((soilV(i2d,j2d) .eq. 28) .or. (soilV(i2d,j2d) .eq. 39))  then
                lctype(i,v) = 4
            elseif (soilV(i2d,j2d) .eq. 29) then
                lctype(i,v) = 1
            elseif (soilV(i2d,j2d) .eq. 30) then
                lctype(i,v) = 5
            elseif (soilV(i2d,j2d) .eq. 32) then
                ! lctype(i,v) = 10
                if (lai2d(i2d,j2d) .ge. 2.5) then
                    lctype(i,v) = 6
                else
                    lctype(i,v) = 7
                end if
            elseif ((soilV(i2d,j2d) .eq. 33) .or. (soilV(i2d,j2d) .eq. 37) .or. (soilV(i2d,j2d) .eq. 40)) then
                lctype(i,v) = 10
            elseif (soilV(i2d,j2d) .eq. 38) then
                lctype(i,v) = 12
                ! lctype(i,v) = 13
            elseif (soilV(i2d,j2d) .eq. 27) then
                lctype(i,v) = 13
            elseif (soilV(i2d,j2d) .eq. 17) then
                lctype(i,v) = 17        
            else
                print *, "ERROR: undefined land cover type.", soilV(i2d,j2d)
                call stopcm1
            end if
        else
            lctype(i,2) = lctype_default
        end if
        lctype(i,1) = 0
        ! Bare soil component
        v = 1
        rcmin(i,v) = 0.0
        zref(i,v) = 0.005
        critleafwpot(i,v) = 0.0
        ! Vegetation component
        v = 2
        soilV_dchm(i2d,j2d) = lctype(i,2)

        if (lctype(i,v).eq.0) then                                  ! Bare Soil
            rcmin(i,v) = 0.0
            zref(i,v) = 0.005
            critleafwpot(i,v) = 0.0
            sai = 0
            hcanopytop = 0
            hcanopybot = 0
        ELSEif (lctype(i,v).eq.1) then                              ! ENF
            rcmin(i,v)=175.0
            zref(i,v)=0.50
            critleafwpot(i,v)=-2500.0
            sai = 0.5
            hcanopytop = 20.0
            hcanopybot = 8.5
        ELSEif (lctype(i,v).eq.2) then                              ! EBF
            rcmin(i,v)=150.0 
            zref(i,v)=0.50 
            critleafwpot(i,v)=-3000.0 
            sai = 0.5
            hcanopytop = 20.0
            hcanopybot = 8.0
        ELSEif (lctype(i,v).eq.3) then                              ! DNF
            rcmin(i,v)=175.0 
            zref(i,v)=0.50 
            critleafwpot(i,v)=-2500.0 
            sai = 0.7
            hcanopytop = 18.0
            hcanopybot = 7.0
        ELSEif (lctype(i,v).eq.4) then                              ! DBF
            rcmin(i,v)=175.0 
            zref(i,v)=0.50 
            critleafwpot(i,v)=-3000.0 
            sai = 0.4
            hcanopytop = 16.0
            hcanopybot = 11.5
        ELSEif (lctype(i,v).eq.5) then                              ! MF
            rcmin(i,v)=175.0 
            zref(i,v)=0.50 
            critleafwpot(i,v)=-2500.0 
            sai = 0.4
            hcanopytop = 16.0
            hcanopybot = 10.0
        ELSEif (lctype(i,v).eq.6) then                              ! Closed Shrubland
            rcmin(i,v)=175.0 
            zref(i,v)=0.10 
            critleafwpot(i,v)=-2500.0 
            sai = 0.3
            hcanopytop = 1.1
            hcanopybot = 0.1
        ELSEif (lctype(i,v).eq.7) then                              ! Open Shrubland
            rcmin(i,v)=178.87 
            zref(i,v)=0.10 
            critleafwpot(i,v)=-3000.0 
            sai = 0.2
            hcanopytop = 1.1
            hcanopybot = 0.1
        ELSEif (lctype(i,v).eq.8) then                              ! Woody Savanna (WSV)
            rcmin(i,v)=169.15 
            zref(i,v)=0.15 
            critleafwpot(i,v)=-2500.0 
            sai = 0.4
            hcanopytop = 13.0
            hcanopybot = 0.1
        ELSEif (lctype(i,v).eq.9) then                              ! Savannas
            rcmin(i,v)=165.0 
            zref(i,v)=0.15 
            critleafwpot(i,v)=-4000.0 
            sai = 0.4
            hcanopytop = 10.0
            hcanopybot = 0.1
        ELSEif (lctype(i,v).eq.10) then                             ! GRA
            rcmin(i,v)=165.0 
            zref(i,v)=0.11 
            critleafwpot(i,v)=-3500.0 
            sai = 0.4
            hcanopytop = 1.0
            hcanopybot = 0.05
        ELSEif (lctype(i,v).eq.12) then                             ! CRO
            rcmin(i,v)=117.12
            zref(i,v)=0.10 
            critleafwpot(i,v)=-1800.0 
            sai = 0.3
            hcanopytop = 2.0
            hcanopybot = 0.1
        elseif (lctype(i,v).eq.13) then                             ! Barren 
            rcmin(i,v)=175
            zref(i,v)=0.01
            critleafwpot(i,v)=-2000.0
            sai = 0
            hcanopytop = 0
            hcanopybot = 0
            cv(i,vbs) = 0
            cv(i,3-vbs) = 1
        else                                                        ! NO DATA                  
            if (WSAT(i).lt.0.99) then
                print *, "vegation type:", lctype(i,v), " ERROR: NO DATA for vegetation type."
                call stopcm1
            end if
            rcmin(i,v)=0.0 
            zref(i,v)=znt(i2d,j2d) 
            critleafwpot(i,v)=0.0 
            sai = 0
            hcanopytop = 0
            hcanopybot = 0
        end if 
        hcanopy = (hcanopytop+hcanopybot)*0.5
        hcanopy2d(i2d,j2d) = hcanopy
        hcanopytop2d(i2d,j2d) = hcanopytop
        hcanopybot2d(i2d,j2d) = hcanopybot
        sai2d(i2d,j2d) = sai

        if (hcanopy .ge. 0.01) then
            a0(i) = 0.16 ! grassland albedo
            emgs(i) = 0.96 ! grassland emissivity
        end if 

        if (hcanopy .ge. 2.) then
            harray(i,1) = 0
            harray(i,2) = 2
            harray(i,3) = hcanopy
            hfarray(i) = 2
            hvarray(i,1) = 2
            hvarray(i,2) = 3
            shading(i,1:Nh,1:Nv) = reshape((/0, 2, 0, 2, 0, 0  /), shape=(/3,2/), order=(/2,1/) )
        else
            harray(i,1) = 0
            harray(i,2) = hcanopy
            harray(i,3) = 2
            hfarray(i) = 3
            hvarray(i,1) = 3
            hvarray(i,2) = 2
            shading(i,1:Nh,1:Nv) = reshape((/0, 2, 0, 0, 0, 0  /), shape=(/3,2/), order=(/2,1/) )
        end if

        ! if (hcanopy .le. 3.) then
        !     shading(i,1:Nh,1:Nv) = reshape((/0, 0, 0, 0, 0, 0  /), shape=(/3,2/), order=(/2,1/) )
        ! end if

        

        !**************** VEGETATION PARAMETERS ****************
        ! LAI [vlai][-]
        ! Vegetation cover fraction [cv] [-]
        ! Fraction of photosynthetically active radiation [fpar] [-] 
        !*******************************************************    
        vlai(i,vbs) = 0.0
        vlai(i,3-vbs) = lai2d(i2d,j2d)
        cv(i,vbs) = 1.0 - cv2d(i2d,j2d)
        cv(i,3-vbs) = cv2d(i2d,j2d)
        fpar(i,vbs) = 0
        fpar(i,3-vbs) = fpar2d(i2d,j2d)
        cvi = 0; cvb = 0.0
        DO v = 1, Nv
            if (vlai(i,v).LE.0.0) vlai(i,v) = 0.0
            if (fpar(i,v).LE.0.0) fpar(i,v) = 0.0
            if (fpar(i,v).GE.1.0) fpar(i,v) = 1.0

            if ((v.NE.vbs).AND.(cv(i,v).GT.cvb)) then
                cvi = v
                cvb = cv(i,v)
            end if
            
        ENDDO ! DO v = 1, Nv
        if (cvi.eq.0) then
            do v = 1, Nv
                if (v.NE.vbs) then
                    cvi = v
                    cvb = cv(i,v)
                end if
            end do
        end if
        zrefs(i)=zref(i,cvi)
        zrefl(i)=zref(i,cvi)
        zreff(i)=zref(i,cvi)
        a0(i) = albd(i2d,j2d)
        emgs(i) = sfem(i2d,j2d)
        if (WSAT(i).lt.0.99) then
            znt(i2d,j2d) = zref(i,cvi)
        end if      
    end do
    cvsep = cv


    ! Initialize soil water content
    if (irst.eq.0) then      
        do i = 1, Ngrid, 1
            call unpack_coord(i, i2d, j2d)      
            do j=1,Ntotlayer
                if (isfc .eq. 2023) then
                    soilW(i,j) = soilW_dchm(i2d,j2d,j)
                    if (soilW(i,j)/dsoil(i,j) > WSAT(i)) then
                        if((soilW(i,j)/dsoil(i,j)/WSAT(i)-1.gt. 0.01) .and. (j < Ntotlayer)) then
                            print *, "WARNING: soil moisture exceeds porosity!", "layer:", j, soilW(i,j)/dsoil(i,j),WSAT(i)
                            ! call stopcm1
                        end if
                        soilW(i,j) = WSAT(i)*dsoil(i,j) 
                    end if
                    do v = 1, Nv, 1
                        soilT(i,j,v) = soilT_dchm(i2d,j2d,j,v)
                    end do
                else
                    soilW(i,:) = 0.75*WSATi(i,:)
                    do v = 1, Nv, 1
                        soilT(i,j,v)=tsk(i2d,j2d)
                    end do
                end if
                soilTavg(i,j) = SUM(cvsep(i,:)*soilT(i,j,:))
            enddo
            if (WSAT(i).lt.0.99) then
                soilW(i,:) = soilW(i,:)*scale_moisture
                mavail(i2d,j2d) = soilW(i,1)/dsoil(i,1)
            else
                if (mavail(i2d,j2d) < 0.99) then
                    print * , "ERROR: WSAT and mavail mismatch", mavail(i2d,j2d), WSAT(i), lctype(i,2), soiltext(i)
                    call stopcm1
                end if
                mavail(i2d,j2d) = 1
                WSAT(i) = 1
            end if
        end do
    end if

    if (myid.eq.0) then
        print *, "dchm_init complete"
    end if
    

end subroutine dchm_init

subroutine dchm_restart_read(soilTavg_dchm, soilW_dchm, soilT_dchm, wsk_dchm, dailyrootup_dchm, dailyrutot_dchm, dailysoilW_dchm, rootfrc_dchm, filtot_dchm, wsksoil_dchm, rtime_daystart_dchm)
    real, dimension(ibl:iel,jbl:jel,Ntotlayer), intent(in) :: soilTavg_dchm, soilW_dchm, dailysoilW_dchm
    real, dimension(ibl:iel,jbl:jel,Ntotlayer,Nv), intent(in) :: soilT_dchm
    real, dimension(ibl:iel,jbl:jel,Ntotlayer,Nv-1), intent(in) :: rootfrc_dchm, dailyrootup_dchm, dailyrutot_dchm
    real, dimension(ibl:iel,jbl:jel,Nv), intent(in) :: wsk_dchm
    real, dimension(ibl:iel,jbl:jel), intent(in) :: filtot_dchm, wsksoil_dchm
    real, intent(in) :: rtime_daystart_dchm
    integer:: i, j, v, vv, i2d, j2d
    rtime_daystart = rtime_daystart_dchm
    do i=1,Ngrid
        call unpack_coord(i, i2d, j2d)
        filtot(i) = filtot_dchm(i2d,j2d)
        wsksoil(i) = wsksoil_dchm(i2d,j2d)
        do v=1, Nv, 1
            wsk(i,v) = wsk_dchm(i2d,j2d,v)
        end do
        do j=1,Ntotlayer
            soilW(i,j) = soilW_dchm(i2d,j2d,j)
            soilTavg(i,j) = soilTavg_dchm(i2d,j2d,j)
            dailysoilW(i,j) = dailysoilW_dchm(i2d,j2d,j)
            do v=1, Nv, 1
                soilT(i,j,v) = soilT_dchm(i2d,j2d,j,v)
            end do
            do vv=1, Nv-1, 1
                rootfrc(i,j,vv) = rootfrc_dchm(i2d,j2d,j,vv)
                dailyrootup(i,j,vv) = dailyrootup_dchm(i2d,j2d,j,vv)
                dailyrutot(i,j,vv) = dailyrutot_dchm(i2d,j2d,j,vv)
            end do
        enddo
    end do
    if (myid.eq.0) then
        print *, "Complete dchm_restart_read."
        ! print *, "maximum rootfrc:", maxval(rootfrc_dchm)
        print *, "rtime_daystart:", rtime_daystart
        if (myid.eq.0) print *, "root frc at restart_read:", maxval(rootfrc), rootfrc(10,:,1)
    end if
end subroutine dchm_restart_read



subroutine unpack_coord(coord_packed, i_unpacked, j_unpacked)
    integer, intent(in) :: coord_packed
    integer, intent(out) :: i_unpacked, j_unpacked
    i_unpacked = mod(coord_packed-1,ni)+1
    j_unpacked = (coord_packed-i_unpacked)/ni+1
    if ((j_unpacked-1)*ni+i_unpacked .ne. coord_packed) then
        print *, "error in calculing coordinates."
        print *, coord_packed, i_unpacked, j_unpacked, ni, nj
        print *, "Ngrid=", Ngrid
        call stopcm1
    end if
end subroutine unpack_coord


end module dchm