!!****if* source/Simulation/SimulationMain/magnetoHD/SZP/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(integer(IN) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!
!!
!!
!!***

subroutine Simulation_initBlock(blockId)
  use Simulation_data
  use Grid_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getBlkPtr, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer, intent(in) :: blockId
  real, dimension(MDIM) :: del
  integer :: i, j, k, n, nn, blkSeed 
  integer :: blkLimits(2, MDIM)
  integer :: blkLimitsGC(2, MDIM)
  integer :: axis(MDIM)
  integer :: sizeX, sizeY, sizeZ, datasize
  real, allocatable :: xcent(:), ycent(:), zcent(:), randoms(:)
  real :: radius, tradActual, rand, dens_pert
  real :: rho, tele, trad, tion, dx, dy, dz
  integer :: species
  real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
  real :: veloc
  real :: romax, romin, kr, r0
  real :: vacu_tradActual
  integer :: error

#ifndef AL_SPEC
  integer :: AL_SPEC = 1
#endif


  !! get the coordinate information for the current block from the database
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  datasize = sizeX*sizeY*sizeZ
  allocate(xcent(sizeX))
  allocate(ycent(sizeY))
  allocate(zcent(sizeZ))
  allocate(randoms(datasize))
  xcent = 0.0
  ycent = 0.0
  zcent = 0.0

  call Grid_getCellCoords(IAXIS, blockId, CENTER, .true., xcent, sizeX)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockId, CENTER, .true., ycent, sizeY)
  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockId, CENTER, .true., zcent, sizeZ)

  call Grid_getDeltas(blockID,del)
  dx = del(1)
  dy = del(2)
  dz = del(3)

  call Grid_getBlkPtr(blockID,solnData,CENTER)
       
#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif


  !! Loop over cells and set the initial state
  nn = 1
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           axis(IAXIS) = i
           axis(JAXIS) = j
           axis(KAXIS) = k

           radius = xcent(i)


            species = AL_SPEC
            dens_pert=sim_line_dens + 0.000*sim_line_dens*cos(2*PI*ycent(j)/sim_lamda)
            rho = dens_pert * exp(-(radius-1)**2/(2.*0.1**2))

            if (rho <= sim_vacu_dens) then
               species = VACU_SPEC

            endif

            if (species == AL_SPEC) then
               rho = rho
            else
               rho = max(rho,sim_line_minDens)  !! density floor in vacuum
            endif


           tele = sim_line_tele
           tion = tele
           trad = tele


           solnData(DENS_VAR,i,j,k) = rho
           solnData(TEMP_VAR,i,j,k) = tele
           
#ifdef FLASH_3T
           solnData(TION_VAR,i,j,k) = tion
           solnData(TELE_VAR,i,j,k) = tele

           ! Set up radiation energy density:
  !         call RadTrans_mgdEFromT(blockId, axis, trad, tradActual)
           solnData(TRAD_VAR,i,j,k) = 1e-9

#endif

           if (NSPECIES > 0) then
              ! Fill mass fractions in solution array if we have any SPECIES defined.
              do n = SPECIES_BEGIN,SPECIES_END
                 if (n==species) then
                    solnData(n,i,j,k) = 1.0e0-(NSPECIES-1)*sim_smallX
                 else
                    solnData(n,i,j,k) = sim_smallX
                 endif
              enddo
           endif

           !! velocities initially zeros
           veloc = 0.
           solnData(VELX_VAR,i,j,k) = veloc
           solnData(VELY_VAR:VELZ_VAR,i,j,k) = 0.0

           !! B-fields initially zero
           solnData(MAGX_VAR:MAGZ_VAR,i,j,k) = 0.0
#if NFACE_VARS > 0
           if (sim_killdivb) then
              facexData(MAG_FACE_VAR,i,j,k) = 0.
              if (NDIM >= 2) faceyData(MAG_FACE_VAR,i,j,k) = 0.
              if (NDIM == 3) facezData(MAG_FACE_VAR,i,j,k) = 0.
           endif
#endif

#ifdef BDRY_VAR
           solnData(BDRY_VAR,i,j,k) = -1.0
#endif

        enddo
     enddo
  enddo

  !! Make sure every proc has same sim_vacu_tradActual
  !!call MPI_AllReduce (vacu_tradActual, sim_vacu_tradActual, 1, FLASH_REAL, MPI_MAX, sim_globalComm, error)
  sim_vacu_tradActual = vacu_tradActual

  !! Make sure face B-fields are definitely zero
#if NFACE_VARS > 0
  if (sim_killdivb) then
     facexData(MAG_FACE_VAR,:,:,:) = 0.
     if (NDIM >= 2) faceyData(MAG_FACE_VAR,:,:,:) = 0.
     if (NDIM == 3) facezData(MAG_FACE_VAR,:,:,:) = 0.
     facexData(MAGI_FACE_VAR,:,:,:) = 0.
     if (NDIM >= 2) faceyData(MAGI_FACE_VAR,:,:,:) = 0.
     if (NDIM == 3) facezData(MAGI_FACE_VAR,:,:,:) = 0.
  endif
#endif

  call Grid_releaseBlkPtr(blockID,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     if (NDIM >= 2) call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif
  deallocate(xcent)
  deallocate(ycent)
  deallocate(zcent)
  deallocate(randoms)

  return

end subroutine Simulation_initBlock
