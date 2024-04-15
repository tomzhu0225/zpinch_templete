!!****if* source/Simulation/SimulationMain/magnetoHD/SZP/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for a particular simulation
!!
!!
!! ARGUMENTS
!!
!!
!!
!! PARAMETERS
!!
!!***

!! Format is easy enough to emulate, main thing is to remember to add new variables (like ones placed in config or par file) to this one, and all other relevant files.


subroutine Simulation_init()
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get, RuntimeParameters_mapStrToInt
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface, ONLY : Driver_getMype, Driver_getComm

  implicit none

#include "constants.h"
#include "Flash.h"

  character(len=MAX_STRING_LENGTH) :: str

  call Driver_getMype(MESH_COMM,sim_meshMe)
  call Driver_getMype(GLOBAL_COMM, sim_globalMe)
  call Driver_getComm(GLOBAL_COMM,sim_globalComm)

  call RuntimeParameters_get('sim_line_dens', sim_line_dens)
  call RuntimeParameters_get('sim_line_minDens', sim_line_minDens)
  call RuntimeParameters_get('sim_vacu_dens',sim_vacu_dens)
  call RuntimeParameters_get('sim_line_tele', sim_line_tele)
  call RuntimeParameters_get('sim_line_maxTemp', sim_line_maxTemp)
  call RuntimeParameters_get('sim_lamda', sim_lamda)

  call RuntimeParameters_get('smallX', sim_smallX)

#ifdef FLASH_USM_MHD
  call RuntimeParameters_get('killdivb',sim_killdivb)
#endif

  call RuntimeParameters_get('geometry', str)
  call RuntimeParameters_mapStrToInt(str, sim_geometry)

  call RuntimeParameters_get('xmin', sim_xmin)
  call RuntimeParameters_get('xmax', sim_xmax)
  call RuntimeParameters_get('ymin', sim_ymin)
  call RuntimeParameters_get('ymax', sim_ymax)

end subroutine Simulation_init

