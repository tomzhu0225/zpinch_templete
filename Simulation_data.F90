!!****if* source/Simulation/SimulationMain/magnetoHD/SZP/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!  Use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data
!!
!!
!!***
module Simulation_data

  implicit none

#include "constants.h"

  integer, save :: sim_meshMe

  !! *** Runtime Parameters *** !!

  real, save :: sim_line_dens
  real, save :: sim_line_tele
  real, save :: sim_vacu_dens
  real, save :: sim_lamda
  real, save :: sim_line_maxTemp
  real, save :: sim_line_minDens


  logical, save :: sim_killdivb = .FALSE.
  real, save :: sim_smallX
  integer, save :: sim_geometry

  real, save :: sim_xmin
  real, save :: sim_xmax
  real, save :: sim_ymin
  real, save :: sim_ymax

  integer, save :: sim_globalMe
  integer, save :: sim_globalComm

  real, save :: sim_vacu_tradActual

end module Simulation_data
