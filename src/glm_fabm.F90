!###############################################################################
!#                                                                             #
!# glm_fabm.F90                                                                #
!#                                                                             #
!# The interface between glm and fabm                                          #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2013 - 2018 -  The University of Western Australia                #
!#                                                                             #
!#  This file is part of GLM (General Lake Model)                              #
!#                                                                             #
!#  GLM is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  GLM is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!###############################################################################

#define _FABM_DIMENSION_COUNT_ 1
#define _FABM_DEPTH_DIMENSION_INDEX_ 1
#define _FABM_VECTORIZED_DIMENSION_INDEX_ 1
#define _FABM_VERTICAL_BOTTOM_TO_SURFACE_

#include "fabm.h"


#undef REALTYPE
#ifndef _FORTRAN_SOURCE_
#define _FORTRAN_SOURCE_ 1
#endif

#include "glm.h"

!# As of june 2014 these were no longer defined
#define varname_wind_sf  standard_variables%wind_speed
#define varname_extc     standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux
#define varname_dens     standard_variables%density
#define varname_layer_ht standard_variables%cell_thickness
#define varname_taub     standard_variables%bottom_stress
#define varname_temp     standard_variables%temperature
#define varname_salt     standard_variables%practical_salinity
#define varname_tss      standard_variables%mass_concentration_of_suspended_matter
#define varname_par      standard_variables%downwelling_photosynthetic_radiative_flux
#define varname_par_sf   standard_variables%surface_downwelling_photosynthetic_radiative_flux
#define varname_pres     standard_variables%pressure

#ifdef __GFORTRAN__
#  if __GNUC__ < 8
#    error   "You will need gfortran version 8 or better"
#  endif
#endif


!-------------------------------------------------------------------------------
MODULE glm_fabm
!
   USE ISO_C_BINDING

   USE fabm
   USE fabm_types
   USE ode_solvers

   USE glm_types
   USE glm_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private
!
#include "glm_globals.h"
#include "glm_plot.h"
#include "glm_ncdf.h"
#include "glm_csv.h"
#include "glm_mobl.h"
!
#if USE_DL_LOADER
# define _WQ_INIT_GLM        "wq_init_glm"
# define _WQ_SET_GLM_DATA    "wq_set_glm_data"
# define _WQ_DO_GLM          "wq_do_glm"
# define _WQ_CLEAN_GLM       "wq_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "wq_init_glm_output"
# define _WQ_WRITE_GLM_      "wq_write_glm"
# define _WQ_VAR_INDEX_C     "wq_var_index_c"
# define _WQ_SET_FLAGS       "wq_set_flags"
# define _WQ_IS_VAR          "wq_is_var"
#else
# define _WQ_INIT_GLM        "fabm_init_glm"
# define _WQ_SET_GLM_DATA    "fabm_set_glm_data"
# define _WQ_DO_GLM          "fabm_do_glm"
# define _WQ_CLEAN_GLM       "fabm_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "fabm_init_glm_output"
# define _WQ_WRITE_GLM_      "fabm_write_glm"
# define _WQ_VAR_INDEX_C     "fabm_var_index_c"
# define _WQ_SET_FLAGS       "fabm_set_flags"
# define _WQ_IS_VAR          "fabm_is_var"
#endif
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   AED_REAL :: lKw    !# background light attenuation (m**-1)
   LOGICAL  :: lIce = .FALSE.

   !# Namelist variables
   INTEGER :: ode_method = 1, split_factor = 1, benthic_mode = 0
   LOGICAL :: bioshade_feedback = .TRUE., repair_state = .TRUE.
   LOGICAL :: mobility_off = .FALSE.  !# flag to turn mobility off
   LOGICAL :: do_plots = .TRUE.

   !# Model
   TYPE (type_model),POINTER :: model

   !# Arrays for state and diagnostic variables
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: cc !# water quality array: nlayers, nvars
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: cc_diag
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: cc_diag_hz
!  AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: sed_zones

   !# Arrays for work, vertical movement, and cross-boundary fluxes
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: rhs_flux
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: ws
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: total
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: local

   !# Arrays for environmental variables not supplied externally.
   AED_REAL,ALLOCATABLE,DIMENSION(:) :: par,pres,tss

   !# External variables
   AED_REAL :: dt, dt_eff   ! External and internal time steps
!  INTEGER  :: w_adv_ctr    ! Scheme for vertical advection (0 IF not used)
   INTEGER  :: n_vars, n_vars_ben
   AED_REAL,POINTER,DIMENSION(:) :: rad, z, salt, temp, rho, area
   AED_REAL,POINTER,DIMENSION(:) :: extc_coef, layer_stress
   AED_REAL,POINTER              :: precip, evap, bottom_stress

   CHARACTER(len=48),ALLOCATABLE :: names(:)
#if PLOTS
   INTEGER,ALLOCATABLE,DIMENSION(:) :: plot_id_v, plot_id_sv, plot_id_d, plot_id_sd
#endif

   AED_REAL,ALLOCATABLE :: dz(:)         !# layer thickness
!===============================================================================
CONTAINS



!###############################################################################
FUNCTION MYTRIM(str) RESULT(res)
!-------------------------------------------------------------------------------
! Useful for passing string arguments to C functions
!-------------------------------------------------------------------------------
   CHARACTER(*),TARGET :: str
   CHARACTER(:),POINTER :: res
   INTEGER :: len

   len = LEN_TRIM(str)+1
   str(len:len) = achar(0)
   res => str
END FUNCTION MYTRIM
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION f_get_lun()
!-------------------------------------------------------------------------------
! Find the first free logical unit number
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER :: lun
   LOGICAL :: opened
!
!-------------------------------------------------------------------------------
!BEGIN
   DO lun = 10,99
      inquire(unit=lun, opened=opened)
      IF ( .not. opened ) THEN
         f_get_lun = lun
         RETURN
      ENDIF
   ENDDO
   f_get_lun = -1
END FUNCTION f_get_lun
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_set_flags(c_split_factor, c_mobility, c_bioshade,              &
                  c_repair_state, c_ode, c_benthic_mode, c_do_plots) BIND(C, name=_WQ_SET_FLAGS)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLOGICAL,INTENT(in) :: c_mobility, c_bioshade, c_repair_state, c_do_plots
   CINTEGER,INTENT(in) :: c_split_factor, c_ode, c_benthic_mode
!
!-------------------------------------------------------------------------------
!BEGIN
   split_factor = c_split_factor
   mobility_off = c_mobility
   bioshade_feedback = c_bioshade
   ode_method = c_ode
   repair_state = c_repair_state
   benthic_mode = c_benthic_mode
   do_plots = c_do_plots
END SUBROUTINE fabm_set_flags
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_init_glm(i_fname,len,MaxLayers,NumWQVars,NumWQBen,pKw) BIND(C, name=_WQ_INIT_GLM)
!-------------------------------------------------------------------------------
! Initialize the GLM-FABM driver by reading settings from fabm.nml.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: i_fname(*)
   CSIZET,INTENT(in)     :: len
   CINTEGER,INTENT(in)   :: MaxLayers
   CINTEGER,INTENT(out)  :: NumWQVars, NumWQBen
   AED_REAL,INTENT(in)   :: pKw
!
!LOCALS
   INTEGER :: i,rc,namlst
   CHARACTER(len=80) :: fname
   INTEGER :: TotWQVars
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(fname, i_fname, len)

   lKw = pKw

#ifdef __INTEL_COMPILER
   print *,'glm_fabm built using intel fortran version ', __INTEL_COMPILER
#else
   print *,'glm_fabm built using gfortran version ', __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
#endif
   print *,'init_glm_fabm from ', TRIM(fname)
   namlst = f_get_lun()

   !# Create model tree
   model => fabm_create_model_from_file(namlst,fname)

   !# Initialize model tree (creates metadata and assigns variable identifiers)
   CALL fabm_set_domain(model,MaxLayers)
   CALL model%set_bottom_index(1)
   CALL model%set_surface_index(MaxLayers)

   print*,'FABM : n_vars      = ', ubound(model%info%state_variables,1)
   print*,'FABM : n_vars_ben  = ', ubound(model%info%state_variables_ben,1)
   print*,'FABM : n_vars_diag = ', ubound(model%info%diagnostic_variables,1)
   print*,'FABM : n_vars_diag_sheet ', ubound(model%info%diagnostic_variables_hz,1)

   !# Report prognostic variable descriptions
   print *, 'FABM pelagic state variables:'
   DO i=1,ubound(model%info%state_variables,1)
         print *, trim(model%info%state_variables(i)%name), '  ',              &
                trim(model%info%state_variables(i)%units),'  ',                &
                trim(model%info%state_variables(i)%long_name)
   ENDDO

   print *, 'FABM benthic state variables:'
   DO i=1,ubound(model%info%state_variables_ben,1)
      print *, trim(model%info%state_variables_ben(i)%name), '  ',             &
             trim(model%info%state_variables_ben(i)%units),'  ',               &
             trim(model%info%state_variables_ben(i)%long_name)
   ENDDO

   !# Report diagnostic variable descriptions
   print *, 'FABM diagnostic variables defined on the full model domain:'
   DO i=1,ubound(model%info%diagnostic_variables,1)
      print *, trim(model%info%diagnostic_variables(i)%name), '  ',            &
             trim(model%info%diagnostic_variables(i)%units),'  ',              &
             trim(model%info%diagnostic_variables(i)%long_name)
   ENDDO

   print *, 'FABM diagnostic variables defined on a horizontal slice of the model domain:'
   DO i=1,ubound(model%info%diagnostic_variables_hz,1)
      print *, trim(model%info%diagnostic_variables_hz(i)%name), '  ',         &
             trim(model%info%diagnostic_variables_hz(i)%units),'  ',           &
             trim(model%info%diagnostic_variables_hz(i)%long_name)
   ENDDO

#if 0
   init_solver(ode_method,ode_solver)
#else
   !# Report type of solver
   print *, "Using Eulerian solver"
   SELECT CASE (ode_method)
      CASE (1)
         print *, 'Using euler_forward()'
      CASE (2)
         print *, 'Using runge_kutta_2()'
      CASE (3)
         print *, 'Using runge_kutta_4()'
      CASE (4)
         print *, 'Using patankar()'
      CASE (5)
         print *, 'Using patankar_runge_kutta_2()'
      CASE (6)
         print *, 'Using patankar_runge_kutta_4()'
      CASE (7)
         print *, 'Using modified_patankar()'
      CASE (8)
         print *, 'Using modified_patankar_2()'
      CASE (9)
         print *, 'Using modified_patankar_4()'
      CASE (10)
         print *, 'Using emp_1()'
      CASE (11)
         print *, 'Using emp_2()'
      CASE (1003)
         print *, 'Using runge_kutta_4() with pp/dd matrices'
      CASE DEFAULT
         STOP 'init_glm_fabm: no valid ode_method specified in fabm.nml!'
   END SELECT
#endif

   NumWQVars = ubound(model%info%state_variables,1)
   NumWQBen  = ubound(model%info%state_variables_ben,1)
   n_vars = NumWQVars
   n_vars_ben = NumWQBen

   TotWQVars = ubound(model%info%state_variables,1) + ubound(model%info%state_variables_ben,1)

   ALLOCATE(dz(MaxLayers))
   dz = 0.  !# initialise to zero

   !# Now that we know how many vars we need, we can allocate space for them
   ALLOCATE(cc(MaxLayers,TotWQVars))
   cc = 0.         !# initialise to zero

   CALL set_c_wqvars_ptr(cc)

!  print *,"Variable names :"
   ALLOCATE(names(1:TotWQVars),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (names)'
   DO i=1,ubound(model%info%state_variables,1)
      names(i) = trim(model%info%state_variables(i)%name)
   ENDDO
   DO i=1,ubound(model%info%state_variables_ben,1)
      names(ubound(model%info%state_variables,1)+i) = trim(model%info%state_variables_ben(i)%name)
   ENDDO
#if PLOTS
   ALLOCATE(plot_id_v(ubound(model%info%state_variables,1)))
   ALLOCATE(plot_id_sv(ubound(model%info%state_variables_ben,1)))
   ALLOCATE(plot_id_d(ubound(model%info%diagnostic_variables,1)))
   ALLOCATE(plot_id_sd(ubound(model%info%diagnostic_variables_hz,1)))
   plot_id_v = -1; plot_id_sv = -1; plot_id_d = -1; plot_id_sd = -1
#endif

!  print *,"Diag-Variable names :"
!  IF ( .not. allocated(diagnames) ) ALLOCATE(diagnames(ubound(model%info%diagnostic_variables,1)))
!  DO i=1,ubound(model%info%diagnostic_variables,1)
!     diagnames(i) = trim(model%info%diagnostic_variables(i)%name)
!     print *,trim(model%info%diagnostic_variables(i)%name)
!  ENDDO
!
!  -----------------------------------------------------------------------------
!
   !# In terms of memory use, it is a waste to allocate storage for benthic variables across the entire
   !# column (the bottom layer should suffice). However, it is important that all values at a given point
   !# in time are integrated simultaneously in multi-step algorithms. This currently can only be arranged
   !# By storing benthic values together with the pelagic, in a fully depth-explicit array.
   DO i=1,ubound(model%info%state_variables,1)
      cc(:, i) = model%info%state_variables(i)%initial_value
      CALL fabm_link_bulk_state_data(model,i,cc(:, i))
   ENDDO
   DO i=1,ubound(model%info%state_variables_ben,1)
      cc(1, ubound(model%info%state_variables,1)+i) = model%info%state_variables_ben(i)%initial_value
      CALL fabm_link_bottom_state_data(model,i,cc(1, ubound(model%info%state_variables,1)+i))
   ENDDO

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag(MaxLayers,ubound(model%info%diagnostic_variables,1)),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = _ZERO_

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag_hz(ubound(model%info%diagnostic_variables_hz,1)),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = _ZERO_

   !# Allocate array with vertical movement rates (m/s, positive for upwards),
   !# and set these to the values provided by the model.
   ALLOCATE(ws(MaxLayers,ubound(model%info%state_variables,1)),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = _ZERO_
   DO i=1,ubound(model%info%state_variables,1)
      ws(:,i) = model%info%state_variables(i)%vertical_movement
   ENDDO

   !# Allocate array for mass fluxes and initialize these to zero (no flux).
   ALLOCATE(rhs_flux(MaxLayers,TotWQVars),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (rhs_flux)'
   rhs_flux = _ZERO_

   !# Allocate array for photosynthetically active radiation (PAR).
   !# This will be calculated internally during each time step.
   ALLOCATE(par(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = _ZERO_
   CALL fabm_link_bulk_data(model,varname_par,par(1:MaxLayers))

   !# Allocate array for local pressure.
   !# This will be calculated [approximated] from layer depths internally during each time step.
   ALLOCATE(pres(1:MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   pres = _ZERO_
   CALL fabm_link_bulk_data(model,varname_pres,pres(1:MaxLayers))

   !# Allocate arrays for storing local and column-integrated values of diagnostic variables.
   !# These are used during each save.
   ALLOCATE(total(1:ubound(model%info%conserved_quantities,1)),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (total)'
   total = _ZERO_
   ALLOCATE(local(1:ubound(model%info%conserved_quantities,1)),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (local)'

   ALLOCATE(tss(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (tss)'
   tss = _ZERO_

   CALL fabm_link_bulk_data(model,varname_tss,tss)
END SUBROUTINE fabm_init_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION fabm_is_var(id,i_vname,len) BIND(C, name=_WQ_IS_VAR)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)   :: id
   CCHARACTER,INTENT(in) :: i_vname(*)
   CSIZET,INTENT(in)     :: len
!LOCALS
   CHARACTER(len=45) :: vname
   INTEGER :: i
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(vname, i_vname, len)

   DO i=1,ubound(model%info%state_variables,1)
      IF ( TRIM(model%info%state_variables(i)%name) == vname ) THEN
         fabm_is_var=i
#ifdef PLOTS
         plot_id_v(i) = id;
#endif
         RETURN
      ENDIF
   ENDDO

   DO i=1,ubound(model%info%state_variables_ben,1)
      IF ( TRIM(model%info%state_variables_ben(i)%name) == vname ) THEN
         IF (benthic_mode.EQ.1) THEN
            fabm_is_var=i
         ELSE
            fabm_is_var=-i
         ENDIF
#ifdef PLOTS
         plot_id_sv(i) = id;
#endif
         RETURN
      ENDIF
   ENDDO

   DO i=1,ubound(model%info%diagnostic_variables,1)
      IF ( TRIM(model%info%diagnostic_variables(i)%name) == vname ) THEN
         fabm_is_var=i
#ifdef PLOTS
         plot_id_d(i) = id;
#endif
         RETURN
      ENDIF
   ENDDO

   DO i=1,ubound(model%info%diagnostic_variables_hz,1)
      IF ( TRIM(model%info%diagnostic_variables_hz(i)%name) == vname ) THEN
         fabm_is_var=-i
#ifdef PLOTS
         plot_id_sd(i) = id;
#endif
         RETURN
      ENDIF
   ENDDO

   fabm_is_var = 0
END FUNCTION fabm_is_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_set_glm_data(Lake, MaxLayers, MetData, SurfData, dt_) &
                                                 BIND(C, name=_WQ_SET_GLM_DATA)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER, INTENT(in) :: MaxLayers
   TYPE(C_PTR),VALUE :: Lake
   TYPE(MetDataType),TARGET     :: MetData  !# Meteorological data
   TYPE(SurfaceDataType),TARGET :: SurfData !# Surface Data
   AED_REAL,INTENT(in)  :: dt_
!
!LOCALS
   INTEGER i
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(Lake, theLake, [MaxLayers])

   !# Save pointers to external dynamic variables that we need later (in do_glm_fabm)
   z => theLake%Height
   temp => theLake%Temp
   salt => theLake%Salinity
   rho => theLake%Density
   area => theLake%LayerArea
   rad => theLake%Light
   extc_coef => theLake%ExtcCoefSW
   layer_stress => theLake%LayerStress

   IF (benthic_mode .GT. 1) zz => z

   !# At this point we have z_cc allocated, so now we can copy the initial values
   !# from cc benthic vars to it
   DO i=1,n_zones
      z_cc(i, n_vars+1:n_vars+n_vars_ben) = cc(1, n_vars+1:n_vars+n_vars_ben)
   ENDDO

   precip => MetData%Rain
   evap   => SurfData%Evap
   bottom_stress => layer_stress(1)

   !# Copy scalars that will not change during simulation, and are needed in do_glm_fabm)
   dt = dt_

   !# Provide pointers to arrays with environmental variables to FABM.
   CALL fabm_link_bulk_data(model,varname_temp,     temp)
   CALL fabm_link_bulk_data(model,varname_salt,     salt)
   CALL fabm_link_bulk_data(model,varname_dens,     rho)
   CALL fabm_link_bulk_data(model,varname_layer_ht, dz)
   CALL fabm_link_bulk_data(model,varname_extc,     extc_coef)
   IF ( ASSOCIATED(bottom_stress) ) &
      CALL fabm_link_horizontal_data(model,varname_taub, bottom_stress)
   CALL fabm_link_horizontal_data(model,varname_wind_sf, MetData%WindSpeed)
   CALL fabm_link_horizontal_data(model,varname_par_sf,  MetData%ShortWave)

   ! Calculate and save internal time step.
   dt_eff = dt/FLOAT(split_factor)

   ! Trigger an error if FABM hasn't got all it needs from us.
   CALL fabm_check_ready(model)
END SUBROUTINE fabm_set_glm_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_do_glm(wlev, pIce) BIND(C, name=_WQ_DO_GLM)
!-------------------------------------------------------------------------------
!                           wlev is the number of levels used;
!                           nlev is the total num levels in the array
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: wlev
   CLOGICAL,INTENT(in) :: pIce
!
!LOCALS
   AED_REAL :: min_C
   INTEGER  :: i, split
   AED_REAL, POINTER :: fdhz
   AED_REAL, DIMENSION(:),POINTER :: fdiag
 ! INTEGER :: n_vars, n_vars_ben
!
!-------------------------------------------------------------------------------
!BEGIN
   lIce = pIce
   CALL model%set_surface_index(wlev)

   !# re-compute the layer heights
   dz(1) = z(1)
   DO i=2,wlev
      dz(i) = z(i) - z(i-1)
   ENDDO

   !# Calculate local pressure
   pres(1:wlev) = -z(1:wlev)
   bottom_stress => layer_stress(1)

   DO i=1,ubound(model%info%state_variables,1)
      CALL fabm_link_bulk_state_data(model,i,cc(:,i))
   ENDDO
   DO i=1,ubound(model%info%state_variables_ben,1)
      CALL fabm_link_bottom_state_data(model,i,cc(1, ubound(model%info%state_variables,1)+i))
   ENDDO
   cc_diag = 0.
   cc_diag_hz = 0.

   IF ( .NOT. mobility_off ) THEN
      !# Get updated vertical movement (m/s, positive for upwards) for biological state variables.
      CALL fabm_get_vertical_movement(model,1,wlev,ws(1:wlev,:))

      !# (3) Calculate source/sink terms due to settling rising of state
      !# variables in the water column (note that settling into benthos
      !# is done in fabm_do_benthos)

      DO i=1,ubound(model%info%state_variables,1)
         IF (model%info%state_variables(i)%vertical_movement .NE. _ZERO_) THEN
            min_C = model%info%state_variables(i)%minimum
            CALL Mobility(wlev, dt, dz, area, ws(:, i), min_C, cc(:, i))
         ENDIF
      ENDDO
   ENDIF

   !# Repair state before calling FABM
   CALL do_repair_state(wlev,'glm_fabm::do_glm_fabm, after advection/diffusion')

   DO split=1,split_factor
      IF (benthic_mode .GT. 1) THEN
         CALL copy_to_zone(cc, wlev)
         CALL calc_zone_areas(area, wlev, z(wlev))
      ENDIF

      !# Update local light field (self-shading may have changed through changes
      !# in biological state variables) changed to update_light to be inline
      !# with current aed_phyoplankton that requires only surface par then integrates over
      CALL update_light(wlev,bioshade_feedback)

      !# Time-integrate one biological time step
      CALL ode_solver(ode_method,ubound(cc,2),wlev,dt_eff,cc(:,:),right_hand_side_rhs,right_hand_side_ppdd)

      !# Provide FABM with (pointers to) updated state variables.
      DO i=1,ubound(model%info%state_variables,1)
         CALL fabm_link_bulk_state_data(model,i,cc(:, i))
      ENDDO
      DO i=1,ubound(model%info%state_variables_ben,1)
         CALL fabm_link_bottom_state_data(model,i,cc(1, ubound(model%info%state_variables,1)+i))
      ENDDO

      !# Repair state
      CALL do_repair_state(wlev,'glm_fabm::do_glm_fabm, after time integration')

      !# Time-integrate diagnostic variables defined on horizontonal slices, where needed.
      DO i=1,ubound(model%info%diagnostic_variables_hz,1)
         fdhz => fabm_get_horizontal_diagnostic_data(model,i)
         cc_diag_hz(i) = fdhz     !# Simply use last value
      ENDDO

      !# Time-integrate diagnostic variables defined on the full domain, where needed.
      DO i=1,ubound(model%info%diagnostic_variables,1)
         fdiag => fabm_get_bulk_diagnostic_data(model,i)
         cc_diag(1:wlev,i) = fdiag(1:wlev)    !# Simply use last value
      ENDDO
   ENDDO

   IF ( benthic_mode .GT. 1 ) CALL copy_from_zone(cc, cc_diag, cc_diag_hz, wlev)
 ! n_vars     = ubound(model%info%state_variables,1)
 ! n_vars_ben = ubound(model%info%state_variables_ben,1)
 ! IF ( benthic_mode .GT. 1 ) CALL copy_from_zone(cc(:,n_vars+1:n_vars+n_vars_ben), wlev)
END SUBROUTINE fabm_do_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE do_repair_state(wlev,location)
!-------------------------------------------------------------------------------
! Check the current values of all state variables
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)          :: wlev
   CHARACTER(len=*),INTENT(in) :: location
!
!LOCALS
   LOGICAL :: valid = .true., l_repair_state
!
!-------------------------------------------------------------------------------
!BEGIN
   l_repair_state = repair_state
   CALL fabm_check_state(model,1,wlev, l_repair_state, valid)

   IF (.NOT. (valid .OR. repair_state)) THEN
      write(stderr,*) 'FATAL ERROR:  State variable values are invalid and repair is not allowed.'
      write(stderr,*) 'FATAL ERROR: ', location
      STOP 'glm_fabm::do_repair_state'
   ENDIF
END SUBROUTINE do_repair_state
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE right_hand_side_rhs(first,numc,nlev,cc,rhs)
!-------------------------------------------------------------------------------
! Calculate the temporal derivatives as a derivative vector
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL, INTENT(in)     :: first
   INTEGER, INTENT(in)     :: numc,nlev
   AED_REAL, INTENT(in)    :: cc(:,:)
   AED_REAL, INTENT(out)   :: rhs(:,:)
!
!LOCALS
   INTEGER :: i,nvars,k
!
!-------------------------------------------------------------------------------
!BEGIN
   i = numc !# numc not used, this shuts up compiler warning

   !# Shortcut to the number of pelagic state variables.
   nvars = ubound(model%info%state_variables,1)

   !# Provide FABM with (pointers to) the current state.
   DO i=1,ubound(model%info%state_variables,1)
      CALL fabm_link_bulk_state_data(model,i,cc(:,i))
   ENDDO
   DO i=1,ubound(model%info%state_variables_ben,1)
      CALL fabm_link_bottom_state_data(model,i,cc(1,nvars+i))
   ENDDO

   !# If this is not the first step in the (multi-step) integration scheme,
   !# THEN first make sure that the intermediate state variable values are valid.
   IF (.not. first) CALL do_repair_state(nlev,'glm_fabm::right_hand_side_rhs')

   !# Initialization is needed because the different biogeochemical models increment or decrement
   !# the temporal derivatives, rather than setting them directly. This is needed for the simultaenous
   !# running of different coupled BGC models.
   rhs = _ZERO_
   rhs_flux = _ZERO_

   !# Start with calculating all flux terms for rhs in mass/m3/s
   !# Includes (1) surface exchange, (2) benthic flux and (3) settling/rising as calculated by glm

   !# (1) surface exchange
   !# Calculate temporal derivatives due to air water exchange.
   IF (.NOT. lIce) THEN !# no surface exchange under ice cover
      CALL fabm_get_surface_exchange(model,rhs_flux(nlev,:))
      !# Distribute surface flux into pelagic surface layer volume (i.e., divide by layer height).
      rhs(nlev,:) = rhs(nlev,:) + rhs_flux(nlev,:)/dz(nlev)
   ENDIF

   !# (2) benthic flux
   !# Calculate temporal derivatives due to exchanges at the sediment/water interface

   bottom_stress => layer_stress(1)
         !# Calculate temporal derivatives due to benthic processes.
   CALL fabm_do_benthos(model,rhs_flux(1, :),rhs(1, nvars+1:))
   !# Limit flux out of bottom layers to concentration of that layer
   !# i.e. don't flux out more than is there
   rhs_flux(1,:) = max(-1.0 * cc(1,:)*dz(1), rhs_flux(1,:))
   !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
   rhs(1,:) = rhs(1,:) + rhs_flux(1,:)/dz(1)

   IF ( benthic_mode.EQ.1 ) THEN
      DO k=2,nlev
         !# fabm wants a horizontal variable so fake it
         bottom_stress => layer_stress(k)

         !# Calculate temporal derivatives due to benthic processes.
         CALL fabm_do_benthos(model,rhs_flux(k, :),rhs(k, nvars+1:))

         !# Limit flux out of bottom layers to concentration of that layer
         !# i.e. don't flux out more than is there
         rhs_flux(k, :) = max(-1.0 * cc(k, :)*dz(k)/dt_eff, rhs_flux(k, :))

         !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
         rhs(k,:) = rhs(k,:) + rhs_flux(k,:)/dz(k) * (area(k)-area(k-1))/area(k)
      ENDDO
   ENDIF

   !# Add pelagic sink and source terms for all depth levels.
   CALL fabm_do(model,1,nlev,rhs(1:nlev,1:nvars))
END SUBROUTINE right_hand_side_rhs
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE right_hand_side_ppdd(first,numc,nlev,cc,pp,dd)
!-------------------------------------------------------------------------------
! Calculate the temporal derivatives as production/destruction matrices
!-------------------------------------------------------------------------------
!ARGUMENTS
   LOGICAL, INTENT(in)    :: first
   INTEGER, INTENT(in)    :: numc,nlev
   AED_REAL, INTENT(in)   :: cc(:,:)
   AED_REAL, INTENT(out)  :: pp(:,:,:)
   AED_REAL, INTENT(out)  :: dd(:,:,:)
!
!LOCALS
   INTEGER :: i,nvars
!
!-------------------------------------------------------------------------------
!BEGIN
   i = numc !# numc not used, this shuts up compiler warning

   !# Shortcut to the number of pelagic state variables.
   nvars = ubound(model%info%state_variables,1)

   !# Provide FABM with (pointers to) the current state.
   DO i=1,ubound(model%info%state_variables,1)
      CALL fabm_link_bulk_state_data(model,i,cc(:,i))
   ENDDO
   DO i=1,ubound(model%info%state_variables_ben,1)
      CALL fabm_link_bottom_state_data(model,i,cc(1,nvars+i))
   ENDDO

   !# If this is not the first step in the (multi-step) integration scheme,
   !# then first make sure that the intermediate state variable values are valid.
   IF (.not. first) CALL do_repair_state(nlev,'glm_fabm::right_hand_side_ppdd')

   !# Initialiaze production and destruction matrices to zero because FABM
   !# biogeochemical models increment these, rather than set these.
   pp = _ZERO_
   dd = _ZERO_

   !# Calculate temporal derivatives due to benthic processes.
   CALL fabm_do_benthos(model,pp(1,:,:),dd(1,:,:),nvars)

   !# Distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
   pp(1,1:nvars,:) = pp(1,1:nvars,:)/dz(1)
   dd(1,1:nvars,:) = dd(1,1:nvars,:)/dz(1)

   CALL fabm_do(model,1,nlev,pp(1:nlev,1:nvars,1:nvars),dd(1:nlev,1:nvars,1:nvars))
END SUBROUTINE right_hand_side_ppdd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_clean_glm() BIND(C, name=_WQ_CLEAN_GLM)
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN

   ! Deallocate internal arrays
   IF (ALLOCATED(cc_diag))    DEALLOCATE(cc_diag)
   IF (ALLOCATED(cc_diag_hz)) DEALLOCATE(cc_diag_hz)
   IF (ALLOCATED(ws))         DEALLOCATE(ws)
   IF (ALLOCATED(rhs_flux))   DEALLOCATE(rhs_flux)
   IF (ALLOCATED(total))      DEALLOCATE(total)
   IF (ALLOCATED(local))      DEALLOCATE(local)
   IF (ALLOCATED(par))        DEALLOCATE(par)
   IF (ALLOCATED(pres))       DEALLOCATE(pres)

END SUBROUTINE fabm_clean_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(nlev, bioshade_feedback)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column
! based on surface radiation, and background and biotic extinction.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: nlev
   LOGICAL,INTENT(in) :: bioshade_feedback
!
!LOCALS
   INTEGER :: i
   AED_REAL :: zz,localext
   AED_REAL :: localexts(1:nlev)

!
!-------------------------------------------------------------------------------
!BEGIN
   zz = _ZERO_

   localexts = _ZERO_
   CALL fabm_get_light_extinction(model,1,nlev,localexts)

   DO i=nlev,1,-1
      localext = localexts(i)

      zz = zz + 0.5*dz(i)

      IF (i .EQ. nlev) THEN
         par(i) = 0.45 * rad(i) * EXP( -(lKw + localext) * zz )
      ELSE
         par(i) = par(i+1) * EXP( -(lKw + localext) * zz )
      ENDIF
      zz = zz + 0.5*dz(i)

      IF (bioshade_feedback) extc_coef(i) = lKw + localext
   ENDDO
END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE fabm_init_glm_output(ncid,x_dim,y_dim,z_dim,zone_dim,time_dim) BIND(C, name=_WQ_INIT_GLM_OUTPUT)
!-------------------------------------------------------------------------------
!  Initialize the output by defining biogeochemical variables.
!-------------------------------------------------------------------------------
!
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid,x_dim,y_dim,z_dim,zone_dim,time_dim
!
!LOCALS
   INTEGER n
   INTEGER dims(4)

!
!-------------------------------------------------------------------------------
!BEGIN
   n = zone_dim !# zone_dim not used in fabm - this shuts up the compiler warning

   !# Put NetCDF library in define mode.
   CALL define_mode_on(ncid)

   !# Set up dimension indices for 3D (+ time) variables (longitude,latitude,depth,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = z_dim
   dims(4) = time_dim

   !# Add a NetCDF variable for each 3D (+ time) biogeochemical state variable.
   DO n=1,ubound(model%info%state_variables,1)
      model%info%state_variables(n)%externalid = NEW_NC_VARIABLE(ncid,         &
                            TRIM(model%info%state_variables(n)%name),          &
                            LEN_TRIM(model%info%state_variables(n)%name),      &
                            NF90_REALTYPE, 4, dims(1:4))
      CALL set_nc_attributes(ncid,model%info%state_variables(n)%externalid,    &
                            MYTRIM(model%info%state_variables(n)%units),       &
                            MYTRIM(model%info%state_variables(n)%long_name)    &
                            PARAM_FILLVALUE)
   ENDDO

   !# Add a NetCDF variable for each 3D (+ time) biogeochemical diagnostic variable.
   DO n=1,ubound(model%info%diagnostic_variables,1)
      model%info%diagnostic_variables(n)%externalid = NEW_NC_VARIABLE(ncid,    &
                       TRIM(model%info%diagnostic_variables(n)%name),          &
                       LEN_TRIM(model%info%diagnostic_variables(n)%name),      &
                       NF90_REALTYPE, 4, dims(1:4))
      CALL set_nc_attributes(ncid,model%info%diagnostic_variables(n)%externalid, &
                            MYTRIM(model%info%diagnostic_variables(n)%units),  &
                            MYTRIM(model%info%diagnostic_variables(n)%long_name)&
                            PARAM_FILLVALUE)
   ENDDO

   !# Set up dimension indices for 2D (+ time) variables (longitude,latitude,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = time_dim

   !# Add a NetCDF variable for each horizontal slice (+ time) biogeochemical state variable.
   DO n=1,ubound(model%info%state_variables_ben,1)
      model%info%state_variables_ben(n)%externalid = NEW_NC_VARIABLE(ncid,     &
                       TRIM(model%info%state_variables_ben(n)%name),           &
                       LEN_TRIM(model%info%state_variables_ben(n)%name), NF90_REALTYPE, 3, dims(1:3))
      CALL set_nc_attributes(ncid,model%info%state_variables_ben(n)%externalid,&
                            MYTRIM(model%info%state_variables_ben(n)%units),   &
                            MYTRIM(model%info%state_variables_ben(n)%long_name)&
                            PARAM_FILLVALUE)
   ENDDO

   !# Add a NetCDF variable for each 2D (longitude,latitude,time) biogeochemical diagnostic variable.
   DO n=1,ubound(model%info%diagnostic_variables_hz,1)
      model%info%diagnostic_variables_hz(n)%externalid = NEW_NC_VARIABLE(ncid, &
                TRIM(model%info%diagnostic_variables_hz(n)%name),              &
                LEN_TRIM(model%info%diagnostic_variables_hz(n)%name), NF90_REALTYPE, 3, dims(1:3))
      CALL set_nc_attributes(ncid,model%info%diagnostic_variables_hz(n)%externalid, &
                            MYTRIM(model%info%diagnostic_variables_hz(n)%units),    &
                            MYTRIM(model%info%diagnostic_variables_hz(n)%long_name) &
                            PARAM_FILLVALUE)
   ENDDO

   !# Add a variable for each conserved quantity
   DO n=1,ubound(model%info%conserved_quantities,1)
      model%info%conserved_quantities(n)%externalid = NEW_NC_VARIABLE(ncid,    &
          TRIM(TRIM(model%info%conserved_quantities(n)%name)//'_tot'),         &
          LEN_TRIM(TRIM(model%info%conserved_quantities(n)%name)//'_tot'), NF90_REALTYPE, 3, dims(1:3))
      CALL set_nc_attributes(ncid,model%info%conserved_quantities(n)%externalid, &
                         TRIM('m*'//model%info%conserved_quantities(n)%units), &
                         MYTRIM(model%info%conserved_quantities(n)%long_name)  &
                         PARAM_FILLVALUE)
   ENDDO

   !# Take NetCDF library out of define mode (ready for storing data).
   CALL define_mode_off(ncid)
END SUBROUTINE fabm_init_glm_output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!-------------------------------------------------------------------------------
SUBROUTINE fabm_write_glm(ncid,wlev,nlev,lvl,point_nlevs) BIND(C, name=_WQ_WRITE_GLM_)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid,wlev,nlev
   CINTEGER,INTENT(in) :: lvl(*),point_nlevs
!
!LOCALS
   INTEGER  :: n,i
   AED_REAL :: val_out
   CLOGICAL :: last = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,ubound(model%info%state_variables,1)
      CALL fabm_link_bulk_state_data(model,i,cc(:, i))
   ENDDO

   !# Store pelagic biogeochemical state variables.
   DO n=1,ubound(model%info%state_variables,1)
      CALL store_nc_array(ncid,model%info%state_variables(n)%externalid,XYZT_SHAPE,wlev,nlev,array=cc(:, n))
      DO i=1,point_nlevs
         IF (lvl(i) .GE. 0) THEN ; val_out = cc(lvl(i)+1, n)
         ELSE                    ; val_out = missing     ; ENDIF
         CALL write_csv_point(i,model%info%state_variables(n)%name,  &
                             len_trim(model%info%state_variables(n)%name), val_out,"",0,last=last)
      ENDDO
#ifdef PLOTS
      IF ( do_plots .AND. plot_id_v(n).GE.0 ) CALL put_glm_val(plot_id_v(n),cc(1:wlev, n))
#endif
   ENDDO

   !# Store benthic biogeochemical state variables.
   DO n=1,ubound(model%info%state_variables_ben,1)
      CALL store_nc_scalar(ncid,model%info%state_variables_ben(n)%externalid, &
                                 XYT_SHAPE,scalar=cc(1, ubound(model%info%state_variables,1)+n))
#ifdef PLOTS
      IF ( do_plots .AND. plot_id_sv(n).GE.0 ) THEN
         IF ( benthic_mode .GT. 1 ) THEN
            CALL put_glm_val(plot_id_sv(n), cc(1:wlev, ubound(model%info%state_variables,1)+n))
         ELSE
            CALL put_glm_val_s(plot_id_sv(n), cc(1, ubound(model%info%state_variables,1)+n))
         ENDIF
      ENDIF
#endif
   ENDDO

   !# Process and store diagnostic variables defined on the full domain.
   DO n=1,ubound(model%info%diagnostic_variables,1)
      !# Store diagnostic variable values.
      CALL store_nc_array(ncid,model%info%diagnostic_variables(n)%externalid,XYZT_SHAPE,wlev,nlev,array=cc_diag(:, n))
      DO i=1,point_nlevs
         IF (lvl(i) .GE. 0) THEN ; val_out = cc_diag(lvl(i)+1, n)
         ELSE                    ; val_out = missing     ; ENDIF
         CALL write_csv_point(i,model%info%diagnostic_variables(n)%name,  &
                             len_trim(model%info%diagnostic_variables(n)%name), val_out,"",0,last=last)
      ENDDO
#ifdef PLOTS
      IF ( do_plots .AND. plot_id_d(n).GE.0 ) &
         CALL put_glm_val(plot_id_d(n),cc_diag(1:wlev, n))
#endif
   ENDDO

   !# Process and store diagnostic variables defined on horizontal slices of the domain.
   DO n=1,ubound(model%info%diagnostic_variables_hz,1)
      !# Store diagnostic variable values.
      CALL store_nc_scalar(ncid,model%info%diagnostic_variables_hz(n)%externalid,XYT_SHAPE,scalar=cc_diag_hz(n))
#ifdef PLOTS
      IF ( do_plots .AND. plot_id_sd(n).GE.0 ) &
         CALL put_glm_val_s(plot_id_sd(n), cc_diag_hz(n))
#endif
   ENDDO

   !# Integrate conserved quantities over depth.
   total = _ZERO_
  !CALL fabm_get_conserved_quantities(model,1,wlev,local)
  !total = total + dz*local

   !# Store conserved quantity integrals.
   DO n=1,ubound(model%info%conserved_quantities,1)
      CALL store_nc_scalar(ncid,model%info%conserved_quantities(n)%externalid,XYT_SHAPE,scalar=total(n))
   ENDDO
END SUBROUTINE fabm_write_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION fabm_var_index_c(name, len) BIND(C, name=_WQ_VAR_INDEX_C)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: name(*)
   CSIZET,INTENT(in)     :: len
!LOCALS
   CHARACTER(len=len+1) :: tn
   CSIZET               :: i
!BEGIN
   tn = ''
   DO i=1,len
      tn=tn//' '
      tn(i:i) = name(i)
   ENDDO
   fabm_var_index_c = WQVar_Index(tn) - 1
END FUNCTION fabm_var_index_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION WQVar_Index(name)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*) :: name
!
!LOCALS
   INTEGER i
!
!-------------------------------------------------------------------------------
!BEGIN
   DO i=1,ubound(model%info%state_variables,1)+ubound(model%info%state_variables_ben,1)
      IF (name .EQ. names(i)) THEN
         WQVar_Index = i
         RETURN
      ENDIF
   ENDDO
   WQVar_Index = -1
END FUNCTION WQVar_Index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_fabm
