!!****if* source/Simulation/SimulationMain/magnetoHD/SZP/Simulation_adjustEvolution
!!
!! NAME
!!  Simulation_adjustEvolution
!!
!!
!! SYNOPSIS
!!  Simulation_adjustEvolution( integer(IN) :: blkcnt,
!!                              integer(IN) :: blklst(blkcnt),
!!                              integer(IN) :: nstep,
!!                              real(IN) :: dt,
!!                              real(IN) :: stime )
!!
!! DESCRIPTION
!!  This routine is called every cycle. It can be used to adjust
!!  the simulation while it is running.
!!  
!! ARGUMENTS
!!  blkcnt - number of blocks
!!  blklist - block list
!!  nstep - current cycle number
!!  dt - current time step length
!!  stime - current simulation time
!!
!!***
subroutine Simulation_adjustEvolution(blkcnt, blklst, nstep, dt, stime)
  use Simulation_data
  use Driver_data,         ONLY : dr_dtOld
  use Driver_interface,  ONLY : Driver_abortFlash
  use Grid_interface,  ONLY : Grid_getSingleCellVol, Grid_getDeltas, &
                              Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_fillGuardCells, &
                              Grid_releaseBlkPtr, Grid_getCellCoords, Grid_getSingleCellVol
  use Eos_interface, ONLY: Eos_wrapped,Eos_getAbarZbar
  use MagneticResistivity_data, ONLY: res_vacFrac, res_vacDens
 ! use rt_data, ONLY : rt_useMGD
  
  use MagneticResistivity_interface, &
                        ONLY : MagneticResistivity

#ifdef FLASH_GRID_PARAMESH
  use Grid_interface,   ONLY: Grid_getBlkRefineLevel
#endif

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: blkcnt
  integer, intent(in) :: blklst(blkcnt)
  integer, intent(in) :: nstep
  real, intent(in) :: dt
  real, intent(in) :: stime

  integer :: blkLimitsGC(LOW:HIGH,MDIM)
  integer :: blkLimits(LOW:HIGH,MDIM)
  integer :: i, j, k, lb, blockID, istat, sizeX, sizeY, sizeZ
  real, pointer, dimension(:,:,:,:) :: U
  real :: delta(MDIM)
  real :: dx, dy
  real, allocatable, dimension(:) :: xCent
  real :: told, tnow, bphi, bx, by, bz, rCent, phi, dr, curr_now
  integer :: axis(MDIM)
  logical, save :: firstCall = .TRUE.

  integer :: EosMode
  real :: oldTele
  logical :: tempChanged

  real    :: res_eta
  real    :: abar,zbar

  call Grid_fillGuardCells(CENTER,ALLDIR)

  do lb = 1, blkcnt
     blockID = blklst(lb)

     call Grid_getBlkPtr(blockID, U, CENTER)
     
     !!get x coord
     call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)         
     sizeX = blkLimitsGC(HIGH,IAXIS)
     allocate(xCent(sizeX),  stat=istat)
     xCent = 0.0
     call Grid_getCellCoords(IAXIS,blockID, CENTER,    .true., xCent,   sizeX)
     !!get cell size
     call Grid_getDeltas(blockID, delta)
     dx = delta(IAXIS)
     dy = delta(JAXIS)
     if (sim_geometry == CARTESIAN) then
        dr = sqrt(dx**2+dy**2)
     elseif (sim_geometry == CYLINDRICAL) then
        dr = dx
     end if   
              
     tempChanged = .false.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
#ifdef RMGZ_VAR
              !!Calculate the vacuum current based on  2pi*(rB)=uI  ---lybbb
              U(RMGZ_VAR,i,j,k)=xCent(i)*U(MAGZ_VAR,i,j,k)*5e-3
#endif
#ifdef ZBA_VAR
              call Eos_getAbarZbar(solnVec=U(:,i,j,k),abar=abar,zbar=zbar)
              U(ZBA_VAR,i,j,k)=zbar
#endif
#ifdef RESI_VAR
	      call MagneticResistivity(U(:,i,j,k),res_eta)
              U(RESI_VAR,i,j,k)=res_eta
#endif
              if(.false. .and. U(TEMP_VAR,i,j,k) > 3.0e6) then
                 U(TEMP_VAR,i,j,k) = sim_line_maxTemp
                 tempChanged = .true.
              end if
              if(.false. .and. U(DENS_VAR,i,j,k) < res_vacDens) then
                 U(TEMP_VAR,i,j,k) = sim_line_tele
                 tempChanged = .true.
              end if
              
           end do
        end do
     end do

     !! Call Eos if the temperature has been changed by floor or ceiling
     !! This ensures that the pressures are consistent with the changed temperature
     if (tempChanged) then
        EosMode = MODE_DENS_TEMP
        call Eos_wrapped(EosMode, blkLimits, blockID)
     endif

     call Grid_releaseBlkPtr(blockID, U, CENTER)
     deallocate(xCent)

  end do   !! block loop


  call Grid_fillGuardCells(CENTER,ALLDIR)

  firstCall = .FALSE.

  return

end subroutine Simulation_adjustEvolution

