!!****if* source/physics/Diffuse/DiffuseMain/Unsplit/diff_advanceMag
!!
!!  NAME 
!!
!!  diff_advanceMag
!!
!!  SYNOPSIS
!!
!!  call diff_advanceMag(integer(IN)                  :: blockCount,
!!                         integer(IN)                  :: blockList(blockCount),
!!                         real(IN)                     :: dt,
!!                         integer, OPTIONAL, intent(IN):: pass)
!!
!!  DESCRIPTION 
!!      This routine advances the magnetic diffusion equation (i.e., resistivity).
!!      An implicit scheme is used. 
!!
!!      Supported boundary conditions are: 
!!                PERIODIC, OUTFLOW (tested).
!!                DIRICHLET (untested).
!!
!!
!! ARGUMENTS
!!
!!  blockCount   - The number of blocks in the list
!!  blockList(:) - The list of blocks on which the solution must be updated
!!   dt           : The time step
!!   pass         : pass=1 directional order of solution sweep X-Y-Z, 
!!                  pass=2 directional order of solution sweep Z-Y-X.
!!  dt           - The time step
!!
!! SIDE EFFECTS
!!
!!  Updates certain variables in permanent UNK storage to contain the
!!  updated B-field and some auxiliaries.  Invokes a solver (of the magnetic
!!  diffusion equation). On return,
!!     MAG(X,Y,Z)_VAR:  contains updated B-fields for the current simulation time.
!!
!!  May modify certain variables used for intermediate results by the solvers
!!  invoked. The list of variables depends on the Diffuse implementation.
!!  The following information is subject to change.
!!     RES1_VAR:  contains parallel resistivity that was passed to Grid_advanceMagDiffusion
!!     RES2_VAR:  contains perpendicular resistivity that was passed to Grid_advanceMagDiffusion
!!  
!! NOTES
!!
!!  The interface of this subroutine must be explicitly known to code that
!!  calls it.  The simplest way to make it so is to have something like
!!    use diff_interface,ONLY: diff_advanceMag
!!  in the calling routine.
!!***

!!REORDER(4): solnVec

#define DEBUG_GRID_GCMASK

subroutine diff_advanceMag(blockCount,blockList,dt,pass)

#include "Flash.h"

  use Diffuse_data, ONLY : diffusion_cutoff_density
  use Diffuse_data, ONLY : useDiffuse, useDiffuseMagneticResistivity, diff_meshMe, diff_meshcomm, &
       diff_resistivitySolver, diff_resistivityForm, diff_useAnisoMagRes, &
       diff_mele, diff_boltz, &
       diff_singleSpeciesA, diff_singleSpeciesZ, diff_avo, diff_mele, &
       diff_magFlMode, diff_magFlCoef, diff_speedlt, &
       diff_geometry
  
  use diff_saData, ONLY : diff_boundary, &
       updateDiffuse, &
       diff_magxDomainBC, diff_magyDomainBC, diff_magzDomainBC, diff_magThetaImplct, &
       diff_xmax, diff_neghLevels
  use Eos_interface, ONLY : Eos_wrapped, Eos, Eos_getAbarZbar
  use Driver_interface, ONLY : Driver_abortFlash
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use MagneticResistivity_interface, ONLY : MagneticResistivity
#ifdef MAGNETIC_RESISTIVITY
  use MagneticResistivity_data, ONLY : res_vacSpecVar, res_vacDens, res_vacFrac
#endif
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getBlkBC, &
      Grid_advanceDiffusion, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
      Grid_getDeltas, GRID_PDE_BND_PERIODIC, GRID_PDE_BND_NEUMANN, &
      GRID_PDE_BND_DIRICHLET, Grid_getCellCoords, Grid_getBlkData

#ifdef DEBUG_GRID_GCMASK
  use Logfile_interface, ONLY : Logfile_stampVarMask
#endif
  use Diffuse_interface, ONLY: Diffuse_solveScalar, &
       Diffuse_solveScalarMag, &
       Diffuse_fluxLimiter, &
       Diffuse_setContextInfo

  use Driver_data, ONLY      : dr_simTime
#ifdef USE_CIRCUIT
  use circ_commonData, ONLY : circ_IloadOld, circ_IloadNew
#endif

  implicit none

#include "constants.h"
#include "Flash_mpi.h"
#include "Eos.h"
#include "Eos_components.h"


  integer,intent(IN)                       :: blockCount
  integer,dimension(blockCount),intent(IN) :: blockList
  real,intent(in)                          :: dt
  integer, OPTIONAL, intent(IN)            :: pass
  
  integer :: i,j,k
  integer :: lb, blockID
  integer :: bcTypes(6,3)
  
  real    :: bcValues(2,6,3)
  real    :: res_zone(2), diff_coeff, xdens, xtemp
  
  integer ::EosMode
  logical :: gcmask(NUNK_VARS)  
  
  integer, dimension(2,MDIM):: blkLimitsGC, blkLimits
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec

  real :: mion
  real :: abar

  logical,save :: gcMaskLogged =.FALSE.

  integer,parameter :: cvarPARRES   = RES1_VAR
  integer,parameter :: cvarPERPRES  = RES2_VAR

  integer :: cvarsRes(2)
  real :: theta, maxBvel, oldMAGP, newMAGP1, newMAGP2, oldENER, newENER, Qohm

  !! variables used for E-field calculation and
  !! subsequent face-centered B-field update
  type :: blockdata
     real, allocatable, dimension(:,:,:) :: oldEx   !! edge-centered E-fields
     real, allocatable, dimension(:,:,:) :: oldEy   !! from time step n
     real, allocatable, dimension(:,:,:) :: oldEz
     real, allocatable, dimension(:,:,:) :: oldQ   !! cell-centered j . E
  end type blockdata
  type(blockdata), allocatable, dimension(:) :: blkdat
  real, allocatable, dimension(:,:,:) :: newEx   !! edge-centered E-fields
  real, allocatable, dimension(:,:,:) :: newEy   !! from time step n+1
  real, allocatable, dimension(:,:,:) :: newEz
  real, allocatable, dimension(:,:,:) :: newQ   !! cell-centered j . E
  real, POINTER, DIMENSION(:,:,:,:) :: facexData, faceyData, facezData 
  real, allocatable, dimension(:,:,:) :: facexEy, facexEz, &  !! face-centered E-fields
                                         faceyEx, faceyEz, &
                                         facezEx, facezEy
  real :: facexBy, facexBz, &  !! face-centered B-fields
          faceyBx, faceyBz, &
          facezBx, facezBy
  integer :: j2D, k3D, ilo, ihi, jlo, jhi, klo, khi, datasize(MDIM), datasizeGC(MDIM)
  real :: dx, dy, dz ,dr
  real, dimension(MDIM) :: del
  integer, dimension(2,MDIM):: faces
  logical :: notLowerIBoundary, notUpperIBoundary, &
             notLowerJBoundary, notUpperJBoundary, &
             notLowerKBoundary, notUpperKBoundary
  real, allocatable :: faceAreas(:,:,:,:) 
  real, allocatable :: rFaceAreas(:,:,:), cellVolumes(:,:,:)
  real :: bbx2, bby2, bbz2, b2, sndspd2, cfx2, cfy2, cfz2
  real :: leftFac, rghtFac
  logical :: calcOhmicHeating, calcOhmicHeatingBlk
  real :: vacDens, vacFrac
  integer :: vacSpecVar

  real :: Excent, Eycent, Ezcent
  real :: jxcent, jycent, jzcent

  !! Needed for AMR check
  logical, save :: firstCall = .true.

  !! for areas
#ifdef FIXEDBLOCKSIZE  
  real, dimension(GRID_IHI_GC) :: xCenter
#else  
  real,  allocatable,  dimension(:):: xCenter
#endif   
  real  ::  rc, rp, rm, Ap,  Am, dV

  !! for circuit BC
  real :: ILoadOld, ILoadNew
  real :: told, tnew, currold, currnew, bphiold, bphinew
 
  !! for current density  --lybbb
  real :: dxBy, dxBz, dyBx, dyBz, dzBx, dzBy
  real    :: Jx,Jy,Jz,idx,idy,idz
  real :: inv_r 
  integer :: sizex
  real,  allocatable,  dimension(:):: xCent
  !=========================================================================  

  if (.not. useDiffuse) return  
  if (.not. useDiffuseMagneticResistivity) return  

#ifdef MAGNETIC_RESISTIVITY
  vacSpecVar = res_vacSpecVar
  vacDens = res_vacDens
  vacFrac = res_vacFrac
#else
   vacSpecVar = NONEXISTENT
   vacDens = -1.
   vacFrac = 0.5
#endif

#ifdef USE_CIRCUIT
   ILoadOld = circ_ILoadOld
   ILoadNew = circ_ILoadNew
#else
   ILoadOld = 0.
   ILoadNew = 0.
#endif

  theta = diff_magThetaImplct

  cvarsRes = (/cvarPARRES, cvarPERPRES/)

  call Timers_start("diff_advanceMag")


!!calculate the current density for cylindrical geometry--lybbb
#ifdef CURZ_VAR
  do lb = 1, blockCount
     blockID = blockList(lb)

     call Grid_getBlkPtr(blockID, solnVec, CENTER)
     
     !!get x coord
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)         
     allocate(xCent(blkLimitsGC(HIGH,IAXIS)))
     call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCent, blkLimitsGC(HIGH,IAXIS))
     !!get cell size
     call Grid_getDeltas(blockID, del)
     dx = del(IAXIS)
     dy = del(JAXIS)

     dr = dx

     !!get 1/size
     idx=1./del(IAXIS)
     if (NDIM >= 2) then
        idy=1./del(JAXIS)
        if (NDIM == 3) then
           idz=1./del(KAXIS)
        endif
     endif
     
     !! calculate the curz  ,only for cylindrical --lybbb
               !! Notice that X == R, Y == Z, Z == PHI. Be aware of signs
               !! when calculating curls    
               !! 1D case : d/dy=d/dz=0
     dzBx = 0.0
     dzBy = 0.0
     dyBx = 0.0
     dyBz = 0.0
     dxBy = 0.0
     dxBz = 0.0
     Jx   = 0.0
     Jy   = 0.0
     Jz   = 0.0
      
              
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
             
              inv_r = 1.0/xCent(i)
              dxBy = (solnVec(MAGY_VAR,i+1,j,k) - solnVec(MAGY_VAR,i-1,j,k))*idx*0.5
              dxBz = (solnVec(MAGZ_VAR,i+1,j,k) - solnVec(MAGZ_VAR,i-1,j,k))*idx*0.5 &
                             +  solnVec(MAGZ_VAR,i,j,k)*inv_r

#if NDIM >= 2

               !! 2D case : d/dy .ne. 0 but d/dz=0
              dyBx = (solnVec(MAGX_VAR,i, j+1,k) - solnVec(MAGX_VAR,i,  j-1,k))*idy*0.5
              dyBz = (solnVec(MAGZ_VAR,i, j+1,k) - solnVec(MAGZ_VAR,i,  j-1,k))*idy*0.5

#endif

               !! Get Jx, Jy, Jz
              Jx = dyBz - dzBy
              Jy = dzBx - dxBz
              Jz = dxBy - dyBx
              
              solnVec(CURZ_VAR,i, j, k) = Jy
              
              !!write(*,*) "mag,idx,inv_r",solnVec(MAGZ_VAR,i,j,k),idx,inv_r
             !! write(*,*) "dxBz,Jy,Curz",dxBz,Jy,solnVec(CURZ_VAR,i, j, k)
              end do
          end do
     end do
     call Grid_releaseBlkPtr(blockID, solnVec, CENTER)
     deallocate(xCent)

  end do   !! block loop
  
#endif  


  if (diff_resistivitySolver == IMPLCT .or. diff_resistivitySolver == RES_HYBRID) then

#ifdef FLASH_GRID_PARAMESH
   call Driver_abortFlash("diff_advanceMag: implicit resistivity solver is not implemented &
                             &for AMR, switch to uniform grid (+ug) or explicit solver.")
     if (firstCall) then
        if (diff_meshMe == MASTER_PE) then
           PRINT*, "WARNING (diff_advanceMag): implicit resistivity with AMR is only implemented for &
                   &MAGZ. MAGX and MAGY should be zero."
        endif
        firstCall = .false.
     endif

     if (NDIM > 2) call Driver_abortFlash("diff_advanceMag: 3D implicit resistivity solver is not implemented &
                             &for AMR, switch to uniform grid (+ug) or explicit solver.")

     if (diff_resistivitySolver == IMPLCT) then
        if (diff_resistivityForm == FULLANISO) call Driver_abortFlash("diff_advanceMag: implicit resistivity solver for &
                                &the full anisotropic equation is not implemented for AMR, switch to uniform grid (+ug) &
                                &or explicit solver, or change resistivityForm to either PARALLEL or PERPENDICULAR.")
     endif

     if (diff_resistivitySolver == RES_HYBRID) call Driver_abortFlash("diff_advanceMag: hybrid resistivity solver is &
          &not fully implemented for AMR, suggest switching resistivitySolver to implicit or explicit.")
#endif

     bcTypes(:,1) = diff_magxDomainBC(:)
     bcTypes(:,2) = diff_magyDomainBC(:)
     bcTypes(:,3) = diff_magzDomainBC(:)
     bcValues = 0.

     !! used in CIRCUIT BC
     told = dr_simTime
     tnew = dr_simTime + dt
     currold = IloadOld
     currnew = IloadNew
     !! Calculate bphi from Ampere's Law [r in cm, curr_now in A, bphi in G normalized]
     !bphiold = 0.2*currold/sqrt(4.*PI)/diff_xmax  !! value on domain boundary
     !bphinew = 0.2*currnew/sqrt(4.*PI)/diff_xmax  !! value on domain boundary
     bphiold = 0.2*currold/diff_xmax*1.02  !! value on domain boundary
     bphinew = 0.2*currnew/diff_xmax*1.02  !! value on domain boundary
     
     where (bcTypes == PERIODIC)
        bcTypes = GRID_PDE_BND_PERIODIC
     elsewhere (bcTypes == DIRICHLET)
        bcTypes = GRID_PDE_BND_DIRICHLET
        bcValues(2,:,:) = -1.0   !! Resistivity value, not used
     elsewhere (bcTypes == OUTFLOW)
        bcTypes = GRID_PDE_BND_NEUMANN
     elsewhere (bcTypes == CIRCUIT) !! Typically for Bphi BC for a Z-pinch
        bcTypes = GRID_PDE_BND_DIRICHLET
        bcValues(1,:,:) = theta*bphinew + (1.-theta)*bphiold   !! B-field value
        bcValues(2,:,:) = -1.0   !! Resistivity value, not used
     end where

     gcmask(:) = .FALSE.
     gcmask(MAGX_VAR) = .TRUE.
     gcmask(MAGY_VAR) = .TRUE.
     gcmask(MAGZ_VAR) = .TRUE.


#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(gcmask, .FALSE., '[diff_advanceMag]', 'gcNeed')
     end if
#endif
!     call Grid_fillGuardCells(CENTER,ALLDIR,masksize=NUNK_VARS, &
!          mask=gcmask,selectBlockType=LEAF,doLogMask=.NOT.gcMaskLogged)   
     !! Make sure guard cells are filled before we do anything
     call Grid_fillGuardCells(CENTER,ALLDIR)


!! ********************************************************* !!
!!                                                           !!
!!   1. Get resistivity coefficients                         !!
!!                                                           !!
!! ********************************************************* !!
     do lb = 1, blockCount
        blockID = blockList(lb)
        call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
        call Grid_getBlkPtr(blockID, solnVec)

        do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
           do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
              do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

                 if (diff_resistivityForm == PARALLEL) then
                    call MagneticResistivity(solnVec(:,i,j,k),res_zone(1))
                    res_zone(2) = res_zone(1)
                 else
                    call MagneticResistivity(solnVec(:,i,j,k),res_zone(1),res_zone(2))
                 endif
                 solnVec(cvarPARRES,i,j,k) = res_zone(1)
                 solnVec(cvarPERPRES,i,j,k) = res_zone(2)

                 !! Set the flux limiter (use fast Alven speed):
                 bbx2 = solnVec(MAGX_VAR,i,j,k)**2/solnVec(DENS_VAR,i,j,k)
                 bby2 = solnVec(MAGY_VAR,i,j,k)**2/solnVec(DENS_VAR,i,j,k)
                 bbz2 = solnVec(MAGZ_VAR,i,j,k)**2/solnVec(DENS_VAR,i,j,k)
                 b2   = bbx2 + bby2 + bbz2
                 sndspd2 = solnVec(GAMC_VAR,i,j,k)*solnVec(PRES_VAR,i,j,k)/solnVec(DENS_VAR,i,j,k)
                 cfx2 = .5*((sndspd2+b2)+sqrt((sndspd2-b2)*(sndspd2-b2)+4.*sndspd2*(bby2+bbz2)))
                 cfy2 = .5*((sndspd2+b2)+sqrt((sndspd2-b2)*(sndspd2-b2)+4.*sndspd2*(bbx2+bbz2)))
                 cfz2 = .5*((sndspd2+b2)+sqrt((sndspd2-b2)*(sndspd2-b2)+4.*sndspd2*(bbx2+bby2)))

                 cfx2 = sqrt(cfx2) + abs(solnVec(VELX_VAR,i,j,k))
                 cfy2 = sqrt(cfy2) + abs(solnVec(VELY_VAR,i,j,k))
                 cfz2 = sqrt(cfz2) + abs(solnVec(VELZ_VAR,i,j,k))
                 maxBvel = max(abs(solnVec(MAGY_VAR,i,j,k)),abs(solnVec(MAGZ_VAR,i,j,k)))*cfx2
                 if (NDIM >= 2) then
                    maxBvel = max(maxBvel,max(abs(solnVec(MAGX_VAR,i,j,k)),abs(solnVec(MAGZ_VAR,i,j,k)))*cfy2)
                 endif
                 if (NDIM == 3) then
                    maxBvel = max(maxBvel,max(abs(solnVec(MAGX_VAR,i,j,k)),abs(solnVec(MAGY_VAR,i,j,k)))*cfz2)
                 endif

                 !!maxBvel = maxval(abs(solnVec(MAGX_VAR:MAGZ_VAR,i,j,k)))*diff_speedlt
                 solnVec(FLLM_VAR,i,j,k) = diff_magFlCoef * maxBvel

              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnVec)
     end do  !! block loop

     !! Fill guard cells with resistivity coefficients
     call Grid_fillGuardCells(CENTER,ALLDIR)


!! ********************************************************* !!
!!                                                           !!
!!   2. Time step n calculations:                            !!
!!      - calculate and store edge-centered E-fields         !!
!!      - calculate and store face-centered Poynting flux    !!
!!                                                           !!
!! ********************************************************* !!
     allocate(blkdat(1:blockCount))
     do lb = 1, blockCount

        !! In 1D, there are no edges, but we still need x-faced quantities
        call getFaceEfields()

        allocate(blkdat(lb)%oldQ(ilo:ihi,jlo:jhi,klo:khi))
        blkdat(lb)%oldQ = 0.
#if NDIM > 1
        allocate(blkdat(lb)%oldEz(ilo:ihi+1,jlo:jhi+1,klo:khi))
        if (NDIM == 3) then
           allocate(blkdat(lb)%oldEx(ilo:ihi,jlo:jhi+1,klo:khi+1))
           allocate(blkdat(lb)%oldEy(ilo:ihi+1,jlo:jhi,klo:khi+1))
        endif

        blkdat(lb)%oldEz = 0.
        if (NDIM == 3) then
           blkdat(lb)%oldEx = 0.
           blkdat(lb)%oldEy = 0.
        endif
#endif

        if (diff_geometry == CYLINDRICAL) then
#ifndef FIXEDBLOCKSIZE
           allocate(xCenter(blkLimitsGC(HIGH,IAXIS)))
#endif
           call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
        endif


        do k = klo, khi+k3D
           do j = jlo, jhi+j2D
              do i = ilo, ihi+1

                 !! each edge E-field is an average of 4 face-values
                 if (NDIM == 2 .or. (NDIM == 3 .and. k < khi+1)) then
                    blkdat(lb)%oldEz(i,j,k) = 0.25*(facexEz(i,j,k) + facexEz(i,j-1,k) + &
                                                    faceyEz(i,j,k) + faceyEz(i-1,j,k))
                 endif
#if NDIM == 3
                 if (i < ihi+1) then
                    blkdat(lb)%oldEx(i,j,k) = 0.25*(faceyEx(i,j,k) + faceyEx(i,j,k-1) + &
                                                    facezEx(i,j,k) + facezEx(i,j-1,k))
                 endif

                 if (j < jhi+1) then
                    blkdat(lb)%oldEy(i,j,k) = 0.25*(facexEy(i,j,k) + facexEy(i,j,k-1) + &
                                                    facezEy(i,j,k) + facezEy(i-1,j,k))
                 endif
#endif

                 !! oldQ = j . E
                 if (i < ihi+1 .and. j < jhi+1 .and. k < khi+1) then
                    call getCenterCurrents(solnVec,i,j,k,jxcent,jycent,jzcent)
                    call getCenterEfields(solnVec,i,j,k,jxcent,jycent,jzcent,Excent,Eycent,Ezcent)

                    blkdat(lb)%oldQ(i,j,k) = jxcent*Excent + jycent*Eycent + jzcent*Ezcent
                 endif

              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnVec)
        deallocate(facexEz, facexEy)
        if (NDIM >= 2) deallocate(faceyEz, faceyEx)
        if (NDIM == 3) deallocate(facezEx, facezEy)
#ifndef FIXEDBLOCKSIZE
        if (diff_geometry == CYLINDRICAL) deallocate(xCenter)
#endif

     end do  !! block loop


!! ********************************************************* !!
!!                                                           !!
!!   3. Diffuse zone-centered B-fields (HYPRE does an        !!
!!      implicit solve of magnetic diffusion equation)       !!
!!                                                           !!
!! ********************************************************* !!
#ifdef FLASH_GRID_PARAMESH
     call Diffuse_setContextInfo(group=0, component=MAGZ_VAR)
#else
     call Diffuse_setContextInfo(group=0, component=MAGX_VAR)
#endif

     if (diff_resistivityForm == PARALLEL) then
        ! We're only using parallel component, so limit this component.
        call Diffuse_fluxLimiter(cvarPARRES, MAGX_VAR, FLLM_VAR, &
             diff_magFlMode, blockCount, blockList)
     else
        ! Perpendicular resistivity is larger, so we limit the flux
        ! by limiting this component. Parallel resistivity is
        ! then scaled by the same factor to maintain anisotropy.
        call Diffuse_fluxLimiter(cvarPERPRES, MAGX_VAR, FLLM_VAR, &
             diff_magFlMode, blockCount, blockList)
        do lb = 1, blockCount
           blockID = blockList(lb)
           call Grid_getBlkPtr(blockID, solnVec)
           solnVec(cvarPARRES,:,:,:) = solnVec(cvarPARRES,:,:,:)*solnVec(FLLM_VAR,:,:,:)
           call Grid_releaseBlkPtr(blockID, solnVec)
        enddo
     end if

#ifdef FLASH_GRID_PARAMESH
     !! With AMR, only Bz will be solved with this call
     call Diffuse_solveScalarMag(MAGZ_VAR,cvarsRes, DFCF_VAR, bcTypes,    &
          bcValues, dt, 1.0, 1.0,           &   
          theta,pass, blockCount,blockList)
     !call Diffuse_solveScalar(MAGZ_VAR,RES2_VAR, DFCF_VAR, bcTypes(:,3),    &
     !     bcValues(:,:,3), dt, 1.0, 1.0,           &   
     !     theta,pass, blockCount,blockList)
#else
     !! With UG, all B-field components will be solved with this call
     call Diffuse_solveScalarMag(MAGX_VAR,cvarsRes, DFCF_VAR, bcTypes,    &
          bcValues, dt, 1.0, 1.0,           &   
          theta,pass, blockCount,blockList)
#endif

     !! Make sure guard cells are filled with the new B-fields
     call Grid_fillGuardCells(CENTER,ALLDIR)

     !! For testing just the diffusion, uncomment these lines
     !do lb = 1, blockCount
     !   deallocate(blkdat(lb)%oldEz)
     !   if (NDIM == 3) then
     !      deallocate(blkdat(lb)%oldEx, blkdat(lb)%oldEy)
     !   endif
     !   deallocate(blkdat(lb)%oldQ)
     !end do
     !deallocate(blkdat)
     !call Timers_stop ("diff_advanceMag")
     !return

!! ********************************************************* !!
!!                                                           !!
!!   4. Time step n+1 calculations:                          !!
!!      - calculate new edge-centered E-fields               !!
!!      - calculate new face-centered Poynting flux          !!
!!                                                           !!
!! ********************************************************* !!
     do lb = 1, blockCount
        !! In 1D, there are no edges, but we still need x-faced quantities
        call getFaceEfields()

        if (diff_geometry == CARTESIAN) then
           rc = 1.
           rp = 1.
           rm = 1.
           Ap = 1.
           Am = 1.
           dV = dx
        endif

        allocate(newQ(ilo:ihi,jlo:jhi,klo:khi))
        newQ = 0.
#if NDIM > 1
        allocate(newEz(ilo:ihi+1,jlo:jhi+1,klo:khi))
        if (NDIM == 3) then
           allocate(newEx(ilo:ihi,jlo:jhi+1,klo:khi+1))
           allocate(newEy(ilo:ihi+1,jlo:jhi,klo:khi+1))
        endif

        newEz = 0.
        if (NDIM == 3) then
           newEx = 0.
           newEy = 0.
        endif
#endif

        if (diff_geometry /= CARTESIAN) then
           !get coord info will use this for Areas and Volumes
#ifndef FIXEDBLOCKSIZE
           allocate(xCenter(blkLimitsGC(HIGH,IAXIS)))
#endif
           call Grid_getCellCoords(IAXIS,blockID, CENTER,.true.,xCenter, blkLimitsGC(HIGH,IAXIS))
        endif

        do k = klo, khi+k3D
           do j = jlo, jhi+j2D
              do i = ilo, ihi+1

                 if (NDIM == 2 .or. (NDIM == 3 .and. k < khi+1)) then
                    newEz(i,j,k) = 0.25*(facexEz(i,j,k) + facexEz(i,j-1,k) + &
                                         faceyEz(i,j,k) + faceyEz(i-1,j,k))
                 endif
#if NDIM == 3
                 if (i < ihi+1) then
                    newEx(i,j,k) = 0.25*(faceyEx(i,j,k) + faceyEx(i,j,k-1) + &
                                         facezEx(i,j,k) + facezEx(i,j-1,k))
                 endif

                 if (j < jhi+1) then
                    newEy(i,j,k) = 0.25*(facexEy(i,j,k) + facexEy(i,j,k-1) + &
                                         facezEy(i,j,k) + facezEy(i-1,j,k))
                 endif
#endif

                 !! newQ = j . E
                 if (i < ihi+1 .and. j < jhi+1 .and. k < khi+1) then
                    call getCenterCurrents(solnVec,i,j,k,jxcent,jycent,jzcent)
                    call getCenterEFields(solnVec,i,j,k,jxcent,jycent,jzcent,Excent,Eycent,Ezcent)

                    newQ(i,j,k) = jxcent*Excent + jycent*Eycent + jzcent*Ezcent
                    !!limit newQ to 5*oldQ --lllyb
                  
                 endif

              enddo
           enddo
        enddo

        deallocate(facexEz, facexEy)
        if (NDIM >= 2) deallocate(faceyEz, faceyEx)
        if (NDIM == 3) deallocate(facezEx, facezEy)


!! ********************************************************* !!
!!                                                           !!
!!   5. Update face-centered B-fields                        !!
!!      - uses old and new edge E-fields in order to do the  !!
!!        update implicitly like the diffusion solve         !!
!!      - uses Stokes' theorem                               !!
!!                                                           !!
!! ********************************************************* !!
#if NFACE_VARS > 0
        call Grid_getBlkPtr(blockID,facexData,FACEX)
        call Grid_getBlkPtr(blockID,faceyData,FACEY)
        if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)

        !! Need the face areas
        allocate(faceAreas(NDIM, blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))              
     
        datasizeGC(1:MDIM) = blkLimitsGC(HIGH,1:MDIM) - blkLimitsGC(LOW,1:MDIM)+1

        call Grid_getBlkData(blockID, CELL_FACEAREA, ILO_FACE, EXTERIOR, & 
             blkLimitsGC(LOW,:), faceAreas(IAXIS,:,:,:), datasizeGC)        
     
        call Grid_getBlkData(blockID, CELL_FACEAREA, JLO_FACE, EXTERIOR, &
             blkLimitsGC(LOW,:), faceAreas(JAXIS,:,:,:), datasizeGC)       
     
#if NDIM == 3
        call Grid_getBlkData(blockID, CELL_FACEAREA, KLO_FACE, EXTERIOR, &
             blkLimitsGC(LOW,:), faceAreas(KAXIS,:,:,:), datasizeGC)   
#endif


        do k = klo, khi+k3D
           do j = jlo, jhi+j2D
              do i = ilo, ihi+1
                 if (diff_geometry == CYLINDRICAL) then
                    rc = xCenter(i)            
                    rp = xCenter(i) + 0.5*dx
                    rm = xCenter(i) - 0.5*dx
                    Ap  = abs(rp)
                    Am  = abs(rm)
                    dV  = abs(rc)*dx
                 endif

                 if (j < jhi+1 .and. k < khi+1) then
                    facexData(MAG_FACE_VAR,i,j,k) = facexData(MAG_FACE_VAR,i,j,k) - &
                       theta*dt/dy*(newEz(i,j+1,k) - newEz(i,j,k)) - &
                       (1.-theta)*dt/dy*(blkdat(lb)%oldEz(i,j+1,k) - blkdat(lb)%oldEz(i,j,k))
#if NDIM == 3
                       facexData(MAG_FACE_VAR,i,j,k) = facexData(MAG_FACE_VAR,i,j,k) - &
                          theta*dt/dz*(newEy(i,j,k) - newEy(i,j,k+1)) - &
                          (1.-theta)*dt/dz*(blkdat(lb)%oldEy(i,j,k) - blkdat(lb)%oldEy(i,j,k+1))
#endif
                 endif

                 if (i < ihi+1 .and. k < khi+1) then
                    faceyData(MAG_FACE_VAR,i,j,k) = faceyData(MAG_FACE_VAR,i,j,k) - &
                       theta*dt/dV*(Am*newEz(i,j,k) - Ap*newEz(i+1,j,k)) - &
                       (1.-theta)*dt/dV*(Am*blkdat(lb)%oldEz(i,j,k) - Ap*blkdat(lb)%oldEz(i+1,j,k))
#if NDIM == 3
                       faceyData(MAG_FACE_VAR,i,j,k) = faceyData(MAG_FACE_VAR,i,j,k) - &
                          theta*dt/dz*(newEx(i,j,k+1) - newEx(i,j,k)) - &
                          (1.-theta)*dt/dz*(blkdat(lb)%oldEx(i,j,k+1) - blkdat(lb)%oldEx(i,j,k))
#endif
                 endif

#if NDIM == 3
                 if (i < ihi+1 .and. j < jhi+1) then
                    facezData(MAG_FACE_VAR,i,j,k) = facezData(MAG_FACE_VAR,i,j,k) - &
                       theta*dt/dy*(newEx(i,j,k) - newEx(i,j+1,k)) - &
                       (1.-theta)*dt/dy*(blkdat(lb)%oldEx(i,j,k) - blkdat(lb)%oldEx(i,j+1,k))
                    facezData(MAG_FACE_VAR,i,j,k) = facezData(MAG_FACE_VAR,i,j,k) - &
                       theta*dt/dx*(newEy(i+1,j,k) - newEy(i,j,k)) - &
                       (1.-theta)*dt/dx*(blkdat(lb)%oldEy(i+1,j,k) - blkdat(lb)%oldEy(i,j,k))
                 endif
#endif
              enddo
           enddo
        enddo
        deallocate(blkdat(lb)%oldEz, newEz, faceAreas)
        if (diff_geometry /= CARTESIAN) then
#ifndef FIXEDBLOCKSIZE
           deallocate(xCenter)
#endif
        endif
        if (NDIM == 3) then
           deallocate(blkdat(lb)%oldEx, blkdat(lb)%oldEy, newEx, newEy)
        endif
#endif


!! ********************************************************* !!
!!                                                           !!
!!   6. Update zone-centered B-fields and energies           !!
!!      - zone-centered B-fields are averages of the new     !!
!!        face-centered values                               !!
!!      - the total energy is updated from old and new       !!
!!        cell-centered Ohmic heating making it implicit   !!
!!        like the diffusion solve                           !!
!!                                                           !!
!! ********************************************************* !!

        calcOhmicHeatingBlk = .false.
	
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi

                 xdens = solnVec(DENS_VAR,i,j,k)
                 oldMAGP = solnVec(MAGP_VAR,i,j,k)

                 !! newMAGP1 is from diffused zone-centered B-fields    
                 newMAGP1 = 0.5*dot_product(solnVec(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                           solnVec(MAGX_VAR:MAGZ_VAR,i,j,k)) 

#if NFACE_VARS > 0
                 !! Update zone-centered B-fields from updated face-centered values
                 solnVec(MAGX_VAR,i,j,k) = 0.5*(facexData(MAG_FACE_VAR,i,j,k) + facexData(MAG_FACE_VAR,i+1,j,k))
                 solnVec(MAGY_VAR,i,j,k) = 0.5*(faceyData(MAG_FACE_VAR,i,j,k) + faceyData(MAG_FACE_VAR,i,j+1,k))
                 if (NDIM == 3) &
                    solnVec(MAGZ_VAR,i,j,k) = 0.5*(facezData(MAG_FACE_VAR,i,j,k) + facezData(MAG_FACE_VAR,i,j,k+1))

                 !! newMAGP2 is from the re-updated B-fields from the above averaging 
                 newMAGP2 = 0.5*dot_product(solnVec(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                           solnVec(MAGX_VAR:MAGZ_VAR,i,j,k)) 

                 solnVec(MAGP_VAR,i,j,k) = newMAGP2
#else
                 !! if no face fields, then the updated MAGP is just newMAGP1
                 newMAGP2 = newMAGP1
                 solnVec(MAGP_VAR,i,j,k) = newMAGP1
#endif

                 !! Force Ohmic heating to be zero for vacuum regions
                 !! (we can skip all of the energy calculation below)
                 calcOhmicHeating = .true.
                 if (xdens <= vacDens) then                   
                    if (vacSpecVar /= NONEXISTENT) then
                       if (solnVec(vacSpecVar,i,j,k) > vacFrac) then
                          calcOhmicHeating = .false.
                       endif
                    else
                       calcOhmicHeating = .false.
                    endif
                 endif

                 if (calcOhmicHeating) then
                    calcOhmicHeatingBlk = .true.

                    !! Energy variables in solnVec are mass weighted
                    !! Qohm = (j . E)*dt/xdens to be added to total and internal energies
                    Qohm = dt/xdens * (theta*blkdat(lb)%oldQ(i,j,k) + (1.-theta)*newQ(i,j,k))
                    solnVec(ENER_VAR,i,j,k) = solnVec(ENER_VAR,i,j,k) + Qohm

#ifdef FLASH_3T
                    solnVec(EELE_VAR,i,j,k) = solnVec(EELE_VAR,i,j,k) + Qohm
#else
                    solnVec(EINT_VAR,i,j,k) = solnVec(EINT_VAR,i,j,k) + Qohm
#endif

                    !! Abort if energy has become negative
                    !if (solnVec(ENER_VAR,i,j,k) < 0.) then
                    !   call Driver_abortFlash("ERROR: diff_advanceMag has calculated a negative &
                    !      energy! This was probably caused by a large Poynting flux and may be &
                    !      avoidable by lowering the time step.")
                    !endif
                 endif  !! if calcOhmicHeating
     
              enddo
           enddo
        enddo

        call Grid_releaseBlkPtr(blockID, solnVec)
        deallocate(blkdat(lb)%oldQ, newQ)

#if NFACE_VARS > 0
        call Grid_releaseBlkPtr(blockID, facexData)
        call Grid_releaseBlkPtr(blockID, faceyData)
        if (NDIM == 3) call Grid_releaseBlkPtr(blockID, facezData)
#endif


        !! With the new electron internal energy, a call to Eos will
        !! update the temperatures, pressures, and total internal energy
        !! We can skip this step if Ohmic heating was never calculated
        !! anywhere within the current block.
        if (calcOhmicHeatingBlk) then
#ifdef FLASH_3T
           EosMode = MODE_DENS_EI_GATHER
#else
           EosMode = MODE_DENS_EI
#endif
           call Eos_wrapped(EosMode, blkLimits, blockList(lb))         
        endif

     end do  !! block loop
     deallocate(blkdat)


  end if  !! if resistivitySolver == IMPLCT or RES_HYBRID


  if (.NOT.gcMaskLogged) then
     gcMaskLogged = .TRUE.
  end if
  
  call Timers_stop ("diff_advanceMag")
  
  return
 

contains

  !! This subroutine calculates the face-centered E-fields. Many variables
  !! are not defined here, because this is a contained subroutine that
  !! makes use of all variables defined in the calling subroutine.
  subroutine getFaceEfields()
    implicit none
    real, allocatable, dimension(:) :: xCoord  !! these are local variables
    real :: inv_dVr                            !! all other variables are global
    real :: isoEtaX, isoEtaY, isoEtaZ, &  !! face-valued resistivities
            aniEtaX, aniEtaY, aniEtaZ
    real :: facexJx, facexJy, facexJz, &  !! face-valued currents 
            faceyJx, faceyJy, faceyJz, & 
            facezJx, facezJy, facezJz 
    real :: xBx, xBy, xBz, &  !! normalized face-valued B-fields
            yBx, yBy, yBz, &
            zBx, zBy, zBz
    real :: modB
    real, PARAMETER :: TINY = 1.e-99

    blockID = blockList(lb)
    call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
    ilo = blkLimits(LOW,IAXIS)
    ihi = blkLimits(HIGH,IAXIS) 
    jlo = blkLimits(LOW,JAXIS)
    jhi = blkLimits(HIGH,JAXIS) 
    klo = blkLimits(LOW,KAXIS)
    khi = blkLimits(HIGH,KAXIS)
    j2D = 0
    k3D = 0
    if (NDIM >= 2) j2D = 1
    if (NDIM == 3) k3D = 1

    !! In 2D, we only need Ez to update face-centered Bx and By,
    !! but we need all facex and facey values for Poynting vector
    allocate(facexEz(ilo:ihi+1,jlo-j2D:jhi+j2D,klo-k3D:khi+k3D))
    allocate(facexEy(ilo:ihi+1,jlo-j2D:jhi+j2D,klo-k3D:khi+k3D))
    if (NDIM >= 2) then
       allocate(faceyEz(ilo-1:ihi+1,jlo:jhi+1,klo-k3D:khi+k3D))
       allocate(faceyEx(ilo-1:ihi+1,jlo:jhi+1,klo-k3D:khi+k3D))
    endif
    if (NDIM == 3) then
       allocate(facezEx(ilo-1:ihi+1,jlo-1:jhi+1,klo:khi+1))
       allocate(facezEy(ilo-1:ihi+1,jlo-1:jhi+1,klo:khi+1))
    endif

    call Grid_getDeltas(blockID, del)
    dx = del(IAXIS)
    dy = del(JAXIS)
    dz = del(KAXIS)

    !! Need r-position for cylindrical geometry
    if (diff_geometry == CYLINDRICAL) then
       allocate(xCoord(blkLimitsGC(HIGH,IAXIS)))
       call Grid_getCellCoords(IAXIS, blockID, CENTER, .true., xCoord, blkLimitsGC(HIGH,IAXIS))
    endif

    call Grid_getBlkBC (blockID, faces)
    call Grid_getBlkPtr(blockID, solnVec)

    facexEz = 0. ; facexEy = 0.
    if (NDIM >= 2) then
       faceyEz = 0. ; faceyEx = 0.
    endif
    if (NDIM == 3) then
       facezEx = 0. ; facezEy = 0.
    endif

    !! We need face-centered E-fields on every face within a block plus one
    !! face outside of block in every direction. For example, facexEz(1,0,1)
    !! is Ez on the first x-face below the lower j boundary. This will allow
    !! us to calculate E-field on every edge within a block later.
    do k = klo-k3D, khi+k3D
       do j = jlo-j2D, jhi+j2D
          do i = ilo-1, ihi+1

             !! Using these logicals, we can keep E-fields zero on domain
             !! boundary faces (i.e., a zero flux bc)
             notLowerIBoundary = .FALSE.
             notUpperIBoundary = .FALSE.
             notLowerJBoundary = .FALSE.
             notUpperJBoundary = .FALSE.
             notLowerKBoundary = .FALSE.
             notUpperKBoundary = .FALSE.

             if ((i /= ilo) .or. (faces(1,IAXIS) == NOT_BOUNDARY)) notLowerIBoundary = .TRUE.
             if ((i /= ihi+1) .or. (faces(2,IAXIS) == NOT_BOUNDARY)) notUpperIBoundary = .TRUE.
             if (NDIM >= 2) then
                if ((j /= jlo) .or. (faces(1,JAXIS) == NOT_BOUNDARY)) notLowerJBoundary = .TRUE.
                if ((j /= jhi+1) .or. (faces(2,JAXIS) == NOT_BOUNDARY)) notUpperJBoundary = .TRUE.
             endif
             if (NDIM == 3) then
                if ((k /= klo) .or. (faces(1,KAXIS) == NOT_BOUNDARY)) notLowerKBoundary = .TRUE.
                if ((k /= khi+1) .or. (faces(2,KAXIS) == NOT_BOUNDARY)) notUpperKBoundary = .TRUE.
             endif

             isoEtaX = 0. ; isoEtaY = 0. ; isoEtaZ = 0.
             aniEtaX = 0. ; aniEtaY = 0. ; aniEtaZ = 0.

             facexJx = 0. ; facexJy = 0. ; facexJz = 0.
             faceyJx = 0. ; faceyJy = 0. ; faceyJz = 0.
             facezJx = 0. ; facezJy = 0. ; facezJz = 0.

             !! Calculate face resistivities
             if (diff_resistivityForm == FULLANISO) then
                !! ANISOTROPIC: E = etaperp*j + (etapar-etaperp)*b*(b . j)
                isoEtaX = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i-1,j,k))
                aniEtaX = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i-1,j,k)) - isoEtaX
#if NDIM > 1
                isoEtaY = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i,j-1,k))
                aniEtaY = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i,j-1,k)) - isoEtaY
#if NDIM == 3
                isoEtaZ = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i,j,k-1))
                aniEtaZ = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i,j,k-1)) - isoEtaZ
#endif
#endif
             elseif (diff_resistivityForm == PARALLEL) then
                !! ISOTROPIC: E = etapar*j
                isoEtaX = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i-1,j,k))
#if NDIM > 1
                isoEtaY = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i,j-1,k))
#if NDIM == 3
                isoEtaZ = 0.5*(solnVec(cvarPARRES,i,j,k) + solnVec(cvarPARRES,i,j,k-1))
#endif
#endif
             elseif (diff_resistivityForm == PERPENDICULAR) then
                !! ISOTROPIC: E = etaperp*j
                isoEtaX = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i-1,j,k))
#if NDIM > 1
                isoEtaY = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i,j-1,k))
#if NDIM == 3
                isoEtaZ = 0.5*(solnVec(cvarPERPRES,i,j,k) + solnVec(cvarPERPRES,i,j,k-1))
#endif
#endif
             endif


             !! Calculate face currents
             if (i >= ilo) then
                facexJz = (solnVec(MAGY_VAR,i,j,k) - solnVec(MAGY_VAR,i-1,j,k))/dx
                if (NDIM >= 2) then
                   facexJz = facexJz &
                                    -0.25*(solnVec(MAGX_VAR,i-1,j+1,k) - solnVec(MAGX_VAR,i-1,j-1,k) &
                                            +solnVec(MAGX_VAR,i,j+1,k) - solnVec(MAGX_VAR,i,j-1,k))/dy
                endif
                if (diff_geometry == CARTESIAN) then 
                   facexJy = -(solnVec(MAGZ_VAR,i,j,k) - solnVec(MAGZ_VAR,i-1,j,k))/dx
                elseif (diff_geometry == CYLINDRICAL) then
                   !! -dBz/dx --> 1/r*d(r*Bphi)/dr (ignore sign change, it
                   !! will work out later because cross product in S = E x B) 
                   inv_dVr = xCoord(i)*xCoord(i) - xCoord(i-1)*abs(xCoord(i-1))
                   inv_dVr = 2.0/inv_dVr
                   facexJy = -(solnVec(MAGZ_VAR,i,j,k)*xCoord(i) &
                              -solnVec(MAGZ_VAR,i-1,j,k)*abs(xCoord(i-1)))*inv_dVr
                endif
                if (NDIM == 3) then
                   facexJy = facexJy &
                             +0.25*(solnVec(MAGX_VAR,i-1,j,k+1) - solnVec(MAGX_VAR,i-1,j,k-1) &
                                     +solnVec(MAGX_VAR,i,j,k+1) - solnVec(MAGX_VAR,i,j,k-1))/dz
                endif
             endif

             !! Special treatment for lower-r boundary in cylindrical geometry
             !! The Bphi derivative does not go to zero.
             !if (diff_geometry == CYLINDRICAL .and. (.not. notLowerIBoundary)) then
             !   inv_dVr = xCoord(i)*xCoord(i) - xCoord(i-1)*abs(xCoord(i-1))
             !   inv_dVr = 2.0/inv_dVr
             !   facexJy = -(solnVec(MAGZ_VAR,i,j,k)*xCoord(i) &
             !              -solnVec(MAGZ_VAR,i-1,j,k)*abs(xCoord(i-1)))*inv_dVr
             !endif

#if NDIM > 1
             if (j >= jlo) then
                faceyJx = (solnVec(MAGZ_VAR,i,j,k) - solnVec(MAGZ_VAR,i,j-1,k))/dy
                faceyJz = -(solnVec(MAGX_VAR,i,j,k) - solnVec(MAGX_VAR,i,j-1,k))/dy
                faceyJz = faceyJz + 0.25*(solnVec(MAGY_VAR,i+1,j-1,k) - solnVec(MAGY_VAR,i-1,j-1,k) &
                                  +solnVec(MAGY_VAR,i+1,j,k) - solnVec(MAGY_VAR,i-1,j,k))/dx
                if (NDIM == 3) then
                   faceyJx = faceyJx &
                             -0.25*(solnVec(MAGY_VAR,i,j-1,k+1) - solnVec(MAGY_VAR,i,j-1,k-1) &
                                     +solnVec(MAGY_VAR,i,j,k+1) - solnVec(MAGY_VAR,i,j,k-1))/dz
                endif
             endif

#if NDIM == 3
             if (k >= klo) then
                facezJx = -(solnVec(MAGY_VAR,i,j,k) - solnVec(MAGY_VAR,i,j,k-1))/dz
                facezJy = (solnVec(MAGX_VAR,i,j,k) - solnVec(MAGX_VAR,i,j,k-1))/dz
                facezJx = facezJx + 0.25*(solnVec(MAGZ_VAR,i,j+1,k-1) - solnVec(MAGZ_VAR,i,j-1,k-1) &
                                   +solnVec(MAGZ_VAR,i,j+1,k) - solnVec(MAGZ_VAR,i,j-1,k))/dy
                facezJy = facezJy - 0.25*(solnVec(MAGZ_VAR,i+1,j,k-1) - solnVec(MAGZ_VAR,i-1,j,k-1) &
                                  +solnVec(MAGZ_VAR,i+1,j,k) - solnVec(MAGZ_VAR,i-1,j,k))/dx
             endif
#endif
#endif
             if (diff_useAnisoMagRes) then
                !! Need normal face currents and normalized face B-fields
                if (i >= ilo) then
                   if (NDIM >= 2) then
                      facexJx = 0.25*(solnVec(MAGZ_VAR,i-1,j+1,k) - solnVec(MAGZ_VAR,i-1,j-1,k) &
                                       +solnVec(MAGZ_VAR,i,j+1,k) - solnVec(MAGZ_VAR,i,j-1,k))/dy
                   endif
                   if (NDIM == 3) then
                      facexJx = facexJx &
                                -0.25*(solnVec(MAGY_VAR,i-1,j,k+1) - solnVec(MAGY_VAR,i-1,j,k-1) &
                                        +solnVec(MAGY_VAR,i,j,k+1) - solnVec(MAGY_VAR,i,j,k-1))/dz
                   endif
                   xBx = 0.5*(solnVec(MAGX_VAR,i,j,k) + solnVec(MAGX_VAR,i-1,j,k))
                   xBy = 0.5*(solnVec(MAGY_VAR,i,j,k) + solnVec(MAGY_VAR,i-1,j,k))
                   xBz = 0.5*(solnVec(MAGZ_VAR,i,j,k) + solnVec(MAGZ_VAR,i-1,j,k))
                   modB = sqrt(xBx**2+xBy**2+xBz**2+TINY)
                   xBx = xBx/modB
                   xBy = xBy/modB
                   xBz = xBz/modB
                endif
#if NDIM > 1
                if (j >= jlo) then
                   if (diff_geometry == CARTESIAN) then
                      faceyJy = -0.25*(solnVec(MAGZ_VAR,i+1,j-1,k) - solnVec(MAGZ_VAR,i-1,j-1,k) &
                                       +solnVec(MAGZ_VAR,i+1,j,k) - solnVec(MAGZ_VAR,i-1,j,k))/dx
                   elseif (diff_geometry == CYLINDRICAL) then
                      inv_dVr = xCoord(i+1)*xCoord(i+1) - xCoord(i-1)*abs(xCoord(i-1))
                      inv_dVr = 2.0/inv_dVr
                      faceyJy = -0.5*(solnVec(MAGZ_VAR,i+1,j-1,k)*xCoord(i+1) - &
                                      solnVec(MAGZ_VAR,i-1,j-1,k)*abs(xCoord(i-1)) &
                                      +solnVec(MAGZ_VAR,i+1,j,k)*xCoord(i+1) - &
                                       solnVec(MAGZ_VAR,i-1,j,k)*abs(xCoord(i-1)))*inv_dVr
                   endif

                   if (NDIM == 3) then
                      faceyJy = faceyJy &
                                +0.25*(solnVec(MAGX_VAR,i,j-1,k+1) - solnVec(MAGX_VAR,i,j-1,k-1) &
                                        +solnVec(MAGX_VAR,i,j,k+1) - solnVec(MAGX_VAR,i,j,k-1))/dz
                   endif
                   yBx = 0.5*(solnVec(MAGX_VAR,i,j,k) + solnVec(MAGX_VAR,i,j-1,k))
                   yBy = 0.5*(solnVec(MAGY_VAR,i,j,k) + solnVec(MAGY_VAR,i,j-1,k))
                   yBz = 0.5*(solnVec(MAGZ_VAR,i,j,k) + solnVec(MAGZ_VAR,i,j-1,k))
                   modB = sqrt(yBx**2+yBy**2+yBz**2+TINY)
                   yBx = yBx/modB
                   yBy = yBy/modB
                   yBz = yBz/modB
                endif
#if NDIM == 3
                if (k >= klo) then
                   facezJz = 0.25*(solnVec(MAGY_VAR,i+1,j,k-1) - solnVec(MAGY_VAR,i-1,j,k-1) &
                                    +solnVec(MAGY_VAR,i+1,j,k) - solnVec(MAGY_VAR,i-1,j,k))/dx
                   facezJz = facezJz &
                             -0.25*(solnVec(MAGX_VAR,i,j+1,k-1) - solnVec(MAGX_VAR,i,j-1,k-1) &
                                     +solnVec(MAGX_VAR,i,j+1,k) - solnVec(MAGX_VAR,i,j-1,k))/dy
                   zBx = 0.5*(solnVec(MAGX_VAR,i,j,k) + solnVec(MAGX_VAR,i,j,k-1))
                   zBy = 0.5*(solnVec(MAGY_VAR,i,j,k) + solnVec(MAGY_VAR,i,j,k-1))
                   zBz = 0.5*(solnVec(MAGZ_VAR,i,j,k) + solnVec(MAGZ_VAR,i,j,k-1))
                   modB = sqrt(zBx**2+zBy**2+zBz**2+TINY)
                   zBx = zBx/modB
                   zBy = zBy/modB
                   zBz = zBz/modB
                endif
#endif
#endif
             endif   !! if diff_resistivityForm == FULLANISO


             !! Calculate face E-fields
             !! ISOTROPIC part: Eiso = isoeta*j
             if (i >= ilo) then
                facexEz(i,j,k) = isoEtaX*facexJz
                facexEy(i,j,k) = isoEtaX*facexJy
             endif

             !! Special treatment for lower-r boundary in cylindrical geometry
             !if (diff_geometry == CYLINDRICAL .and. (.not. notLowerIBoundary)) then
             !   facexEy(i,j,k) = isoEtaX*facexJy
             !endif

#if NDIM > 1
             if (j >= jlo) then
                faceyEz(i,j,k) = isoEtaY*faceyJz
                faceyEx(i,j,k) = isoEtaY*faceyJx
             endif

#if NDIM == 3
             if (k >= klo) then
                facezEx(i,j,k) = isoEtaZ*facezJx
                facezEy(i,j,k) = isoEtaZ*facezJy
             endif
#endif
#endif
             if (diff_useAnisoMagRes) then
                !! ANISOTROPIC part: E = Eiso + anieta*b*(b . j)
                if (i >= ilo) then
                   facexEz(i,j,k) = facexEz(i,j,k) + &
                                    aniEtaX*xBz*(xBx*facexJx + xBy*facexJy + xBz*facexJz)
                   facexEy(i,j,k) = facexEy(i,j,k) + &
                                    aniEtaX*xBy*(xBx*facexJx + xBy*facexJy + xBz*facexJz)
                endif
#if NDIM > 1
                if (j >= jlo) then
                   faceyEz(i,j,k) = faceyEz(i,j,k) + &
                                    aniEtaY*yBz*(yBx*faceyJx + yBy*faceyJy + yBz*faceyJz)
                   faceyEx(i,j,k) = faceyEx(i,j,k) + &
                                    aniEtaY*yBx*(yBx*faceyJx + yBy*faceyJy + yBz*faceyJz)
                endif
#if NDIM == 3
                if (k >= klo) then
                   facezEx(i,j,k) = facezEx(i,j,k) + &
                                    aniEtaZ*zBx*(zBx*facezJx + zBy*facezJy + zBz*facezJz)
                   facezEy(i,j,k) = facezEy(i,j,k) + &
                                    aniEtaZ*zBy*(zBx*facezJx + zBy*facezJy + zBz*facezJz)
                endif
#endif
#endif
             endif   !! if diff_anisoMagRes

          enddo
       enddo
    enddo

    if (diff_geometry == CYLINDRICAL) deallocate(xCoord)

    return

  end subroutine getFaceEfields


  !! Gets current at cell centers
  subroutine getCenterCurrents(U,ix,iy,iz,jx,jy,jz)
    implicit none
    real, POINTER, DIMENSION(:,:,:,:) :: U
    integer, intent(in) :: ix, iy, iz
    real, intent(out)  :: jx, jy, jz

    real :: idx, idy, idz, inv_r
    real :: dxBz, dxBy, dyBx, dyBz, dzBx, dzBy

    dzBx = 0.
    dzBy = 0.
    dyBx = 0.
    dyBz = 0.
    dxBy = 0.
    dxBz = 0.
    jx   = 0.
    jy   = 0.
    jz   = 0.

    idx = 1./dx
    if (NDIM >= 2) then
       idy = 1./dy
       if (NDIM == 3) then
          idz = 1./dz
       endif
    endif
 
    if (diff_geometry == CARTESIAN) then

       !! Compute partial derivatives to construct J
  
       !! 1D case : d/dy=d/dz=0
       dxBz = (U(MAGZ_VAR,ix+1,iy,iz) - U(MAGZ_VAR,ix-1,iy,iz))*idx*0.5
       dxBy = (U(MAGY_VAR,ix+1,iy,iz) - U(MAGY_VAR,ix-1,iy,iz))*idx*0.5
#if NDIM >= 2
       !! 2D case : d/dy .ne. 0 but d/dz=0
       dyBx = (U(MAGX_VAR,ix,iy+1,iz) - U(MAGX_VAR,ix,iy-1,iz))*idy*0.5
       dyBz = (U(MAGZ_VAR,ix,iy+1,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy*0.5
#if NDIM == 3
       dzBx = (U(MAGX_VAR,ix,iy,iz+1) - U(MAGX_VAR,ix,iy,iz-1))*idz*0.5
       dzBy = (U(MAGY_VAR,ix,iy,iz+1) - U(MAGY_VAR,ix,iy,iz-1))*idz*0.5
#endif
#endif
    endif

    if (diff_geometry == CYLINDRICAL) then

       !! Notice that X == R, Y == Z, Z == PHI. Be aware of signs
       !! when calculating curls. In CYLINDRICAL, cross products
       !! introduce another minus sign. See hy_uhd_addResistiveFluxes
       !! for more detailed notes on this point. Here, we actually end
       !! up with the wrong signs on j, but the E-fields have the same
       !! wrong sign in cylindrical, so they will cancel in j . E.
       !! 1D case : d/dy=d/dz=0
       inv_r = 1.0/xCenter(ix)
       dxBy = (U(MAGY_VAR,ix+1,iy,iz) - U(MAGY_VAR,ix-1,iy,iz))*idx*0.5
       dxBz = (U(MAGZ_VAR,ix+1,iy,iz) - U(MAGZ_VAR,ix-1,iy,iz))*idx*0.5 &
              +  U(MAGZ_VAR,ix,iy,iz)*inv_r
#if NDIM >= 2
       !! 2D case : d/dy .ne. 0 but d/dz=0
       dyBx = (U(MAGX_VAR,ix,iy+1,iz) - U(MAGX_VAR,ix,iy-1,iz))*idy*0.5
       dyBz = (U(MAGZ_VAR,ix,iy+1,iz) - U(MAGZ_VAR,ix,iy-1,iz))*idy*0.5
#endif
    endif

    !! Get jx, jy, jz
    jx = dyBz - dzBy
    jy = dzBx - dxBz
    jz = dxBy - dyBx

  end subroutine getCenterCurrents


  !! Gets E-fields at cell centers
  subroutine getCenterEFields(U,ix,iy,iz,jx,jy,jz,Ex,Ey,Ez)
    implicit none
    real, POINTER, DIMENSION(:,:,:,:) :: U
    integer, intent(in) :: ix, iy, iz
    real, intent(in)  :: jx, jy, jz
    real, intent(out)  :: Ex, Ey, Ez

    real :: isoEta, aniEta, bx, by, bz, modB, bdotj
    real, PARAMETER :: TINY = 1.e-99

    !! cell-centered B-fields (normalized)
    modB = sqrt(dot_product(U(MAGX_VAR:MAGZ_VAR,ix,iy,iz),U(MAGX_VAR:MAGZ_VAR,ix,iy,iz))) + TINY
    bx = U(MAGX_VAR,ix,iy,iz)/modB
    by = U(MAGY_VAR,ix,iy,iz)/modB
    bz = U(MAGZ_VAR,ix,iy,iz)/modB

    !! cell-centered resistivities
    if (diff_resistivityForm == FULLANISO) then
       !! ANISOTROPIC: E = etaperp*j + (etapar-etaperp)*b*(b . j)
       isoEta = U(cvarPERPRES,ix,iy,iz)
       aniEta = U(cvarPARRES,ix,iy,iz) - isoEta
    elseif (diff_resistivityForm == PARALLEL) then
       isoEta = U(cvarPARRES,ix,iy,iz)
       aniEta = 0.
    elseif (diff_resistivityForm == PERPENDICULAR) then
       isoEta = U(cvarPERPRES,ix,iy,iz)
       aniEta = 0.
    endif

    !! ISOTROPIC part of E
    Ex = isoEta*jx 
    Ey = isoEta*jy 
    Ez = isoEta*jz
 
    if (diff_useAnisoMagRes) then
       !! ANISOTROPIC part: E = Eiso + anieta*b*(b . j)
       bdotj = bx*jx + by*jy + bz*jz
       Ex = Ex + aniEta*bx*bdotj
       Ey = Ey + aniEta*by*bdotj
       Ez = Ez + aniEta*bz*bdotj
    endif

  end subroutine getCenterEFields


end subroutine diff_advanceMag
  
