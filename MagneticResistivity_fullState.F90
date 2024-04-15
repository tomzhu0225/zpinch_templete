!!****if* source/physics/materialProperties/MagneticResistivity/MagneticResistivityMain/SpitzerHighZ/MagneticResistivity_fullState
!!
!! NAME
!!  MagneticResistivity_fullState
!!
!! SYNOPSIS
!!  call MagneticResistivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                                     real(out)   :: resPar,
!!                            OPTIONAL,real(out)   :: resPerp,
!!                            OPTIONAL,real(out)   :: resCross)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron Magnetic Resistivity for all materials,
!! including those with Z > 1. The expressions used here comes from NRL plasma form.
!!
!!  Returns Magnetic Resistivity, parallel and perpendicular components.
!!  resCross is ZERO
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   resPar   :   parallel component of Magnetic Resistivity
!!   resPerp  :   perpendicular component of Magnetic Resistivity
!!   resCross :   cross component of Magnetic Resistivity
!!
!!***

#include "Flash.h"
#include "constants.h"  


subroutine MagneticResistivity_fullState(solnVec,resPar, resPerp, resCross)
  use MagneticResistivity_interface, ONLY: MagneticResistivity
  use MagneticResistivity_data, ONLY: res_meshMe,&
       res_mele, res_qele, res_navo, res_speedlt, res_boltz
  use MagneticResistivity_data, ONLY: res_mUnit
  use MagneticResistivity_data, ONLY: res_coef
  use MagneticResistivity_data, ONLY: res_maxRes, res_vacSpecVar, res_vacDens, res_vacRes
  use Eos_interface, ONLY: Eos_getAbarZbar
  use PlasmaState_interface, ONLY: PlasmaState_logLambda

  implicit none

  real, intent(in)  :: solnVec(NUNK_VARS)
  real, intent(out) :: resPar
  real, OPTIONAL, intent(out) :: resPerp
  real, OPTIONAL, intent(out) :: resCross

  real :: dens
  real :: tele, tele_eV
  real :: tion
  real :: nele
  real :: nion
  real :: mion
  real :: abar
  real :: zbar
  real :: eqtime
  real :: resPerpLoc
  real :: ll
  logical :: firstCall = .true.
  real :: anomRes
  logical :: addedVacRes

  call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

  dens = solnVec(DENS_VAR)
  nion = dens * res_navo / abar
  nele = zbar * nion
  mion = abar / res_navo
#ifdef FLASH_3T
  tele = solnVec(TELE_VAR)
  tion = solnVec(TION_VAR)
#else
  tele = solnVec(TEMP_VAR)
  tion = solnVec(TEMP_VAR)
#endif

  tele_eV = tele/11604.5221
  call PlasmaState_logLambda(ELE_ION, tele, tion, nele, mion, zbar, ll)

  resPerpLoc = 8.21876126127e5*zbar*ll/tele_eV**1.5  !! In CGS -- Here tele has to be in eV

  ! This formula is only valid when the magnetic field is strong (and
  ! only for hydrogen):
  resPar = resPerpLoc/1.96

  !! FOR CGS--> SI unit conversion
  if (res_mUnit == "SI" .or. res_mUnit == "si" ) then
     resPerpLoc = resPerpLoc*1.0e-4
     resPar  = resPar *1.0e-4
  end if

  !! add anomalous vacuum resistivity
  addedVacRes = .false.
  if (dens <= res_vacDens) then
     if (res_vacSpecVar /= NONEXISTENT) then
        if (solnVec(res_vacSpecVar) > 0.5) then
           resPar = res_vacRes
           resPerpLoc = res_vacRes
           addedVacRes = .true.
        end if
     else
        resPar = res_vacRes
        resPerpLoc = res_vacRes
        addedVacRes = .true.
     end if
  end if

  if (present(resPerp)) then
    resPerp = resPerpLoc
    resPerp = min(resPerp, res_maxRes) ! cap resPerp to avoid unphysically cold regions
  end if

  !! cap resPar to avoid unphysically cold regions
  resPar = min(resPar, res_maxRes)
 
  !! return 0 for cross component
  if (present(resCross)) then
     resCross = 0.0
     if (res_meshMe == MASTER_PE .and. firstCall) then
        write(6,*) 'WARNING: [MagneticResistivity_fullState] The SpitzerHighZ implementation'
        write(6,*) '  does not calculate the cross component of resistivity (it is zero by default).'
        write(6,*) '  Consider using the DaviesWen implementation or set useCrossMagRes to false.'
        firstCall = .false.
     end if
  end if


end subroutine MagneticResistivity_fullState
