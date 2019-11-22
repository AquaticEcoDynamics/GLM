!###############################################################################
!#                                                                             #
!# glm_aed2.F90                                                                #
!#                                                                             #
!# The interface between glm and libaed2                                       #
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

#include "aed2.h"

#undef MISVAL
#ifndef _FORTRAN_SOURCE_
#define _FORTRAN_SOURCE_ 1
#endif

#include "glm.h"

#ifdef __GFORTRAN__
#  if __GNUC__ < 8
#    error   "You will need gfortran version 8 or better"
#  endif
#else
#  ifndef isnan
#    define isnan(x) ieee_is_nan(x)
#    define HAVE_IEEE_ARITH
#  endif
#endif


!-------------------------------------------------------------------------------
MODULE glm_aed2
!
   USE ISO_C_BINDING

   USE aed2_common
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
# define _WQ_INIT_GLM        "aed2_init_glm"
# define _WQ_SET_GLM_DATA    "aed2_set_glm_data"
# define _WQ_DO_GLM          "aed2_do_glm"
# define _WQ_CLEAN_GLM       "aed2_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "aed2_init_glm_output"
# define _WQ_WRITE_GLM_      "aed2_write_glm"
# define _WQ_VAR_INDEX_C     "aed2_var_index_c"
# define _WQ_SET_FLAGS       "aed2_set_flags"
# define _WQ_IS_VAR          "aed2_is_var"
#endif
!
!-------------------------------------------------------------------------------
!
!MODULE DATA

   AED_REAL :: lKw    !# background light attenuation (m**-1)
   LOGICAL  :: lIce = .FALSE.

   !# Namelist variables
   INTEGER :: ode_method = 1, split_factor = 1, benthic_mode
   LOGICAL :: bioshade_feedback = .TRUE., repair_state = .TRUE.
   LOGICAL :: mobility_off = .FALSE.  !# flag to turn mobility off
   LOGICAL :: do_plots = .TRUE.

   AED_REAL :: par_fraction =  0.450
   AED_REAL :: nir_fraction =  0.510
   AED_REAL :: uva_fraction =  0.035
   AED_REAL :: uvb_fraction =  0.005

   !# Arrays for state and diagnostic variables
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: cc !# water quality array: nlayers, nvars
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:,:) :: cc_diag
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: cc_diag_hz
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: tss
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: sed_zones

   !# Arrays for work, vertical movement, and cross-boundary fluxes
   AED_REAL,ALLOCATABLE,DIMENSION(:,:) :: ws
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: total
   AED_REAL,ALLOCATABLE,DIMENSION(:)   :: local
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: dz

   !# Arrays for environmental variables not supplied externally.
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: par, pres
   AED_REAL,ALLOCATABLE,TARGET,DIMENSION(:) :: uva, uvb, nir

   !# External variables
   AED_REAL :: dt, dt_eff   ! External and internal time steps
!  INTEGER  :: w_adv_ctr    ! Scheme for vertical advection (0 IF not used)
   AED_REAL,POINTER,DIMENSION(:) :: rad, z, salt, temp, rho, area
   AED_REAL,POINTER,DIMENSION(:) :: extc_coef, layer_stress
   AED_REAL,POINTER :: precip, evap, bottom_stress
   AED_REAL,POINTER :: I_0, wnd
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: depth,layer_area

   CHARACTER(len=48),ALLOCATABLE :: names(:)
   CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)
   AED_REAL,ALLOCATABLE,DIMENSION(:) :: min_, max_

   INTEGER,ALLOCATABLE,DIMENSION(:) :: externalid
#if PLOTS
   INTEGER,ALLOCATABLE,DIMENSION(:) :: plot_id_v, plot_id_sv, plot_id_d, plot_id_sd
#endif

   LOGICAL :: link_rain_loss, link_solar_shade, link_bottom_drag
   AED_REAL,POINTER :: rain_factor, sw_factor, friction

   INTEGER :: n_aed2_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
   INTEGER :: zone_var = 0
   CHARACTER(len=64) :: NULCSTR = ""
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
SUBROUTINE STOPIT(message)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(*) :: message
!-------------------------------------------------------------------------------
   PRINT *,message
   STOP
END SUBROUTINE STOPIT
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
SUBROUTINE aed2_set_flags(c_split_factor, c_mobility, c_bioshade,              &
                  c_repair_state, c_ode, c_benthic_mode, c_do_plots,           &
                  c_link_rain_loss, c_link_solar_shade, c_link_bottom_drag) BIND(C, name=_WQ_SET_FLAGS)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CLOGICAL,INTENT(in) :: c_mobility, c_bioshade, c_repair_state, c_do_plots
   CINTEGER,INTENT(in) :: c_split_factor, c_ode, c_benthic_mode
!  CLOGICAL,OPTIONAL,INTENT(in) :: c_link_rain_loss, c_link_solar_shade, c_link_bottom_drag
   CLOGICAL,INTENT(in) :: c_link_rain_loss, c_link_solar_shade, c_link_bottom_drag
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
!  IF (PRESENT(c_link_rain_loss))   link_rain_loss = c_link_rain_loss
!  IF (PRESENT(c_link_solar_shade)) link_solar_shade = c_link_solar_shade
!  IF (PRESENT(c_link_bottom_drag)) link_bottom_drag = c_link_bottom_drag
   link_rain_loss = c_link_rain_loss
   link_solar_shade = c_link_solar_shade
   link_bottom_drag = c_link_bottom_drag
END SUBROUTINE aed2_set_flags
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_init_glm(i_fname,len,MaxLayers,NumWQ_Vars,NumWQ_Ben,pKw) BIND(C, name=_WQ_INIT_GLM)
!-------------------------------------------------------------------------------
! Initialize the GLM-AED2 driver by reading settings from aed.nml.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: i_fname(*)
   CSIZET,INTENT(in)     :: len
   CINTEGER,INTENT(in)   :: MaxLayers
   CINTEGER,INTENT(out)  :: NumWQ_Vars, NumWQ_Ben
   AED_REAL,INTENT(in)   :: pKw
!
!LOCALS
   INTEGER :: i,j,namlst,status
   INTEGER :: rc, av, v, sv

   CHARACTER(len=80) :: fname
   TYPE(aed2_variable_t),POINTER :: tvar

   CHARACTER(len=64) :: models(64)
   NAMELIST /aed2_models/ models
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(fname, i_fname, len)

   lKw = pKw

#ifdef __INTEL_COMPILER
   print *,'    glm_aed2 built using intel fortran version ', __INTEL_COMPILER
#else
# ifdef __PGI
   print *,'    glm_aed2 built using pgfortran version ', __PGIC__, '.', __PGIC_MINOR__, '.', __PGIC_PATCHLEVEL__
# else
   print *,'    glm_aed2 built using gfortran version ', __GNUC__, '.', __GNUC_MINOR__, '.', __GNUC_PATCHLEVEL__
# endif
#endif
   print *,'    libaed2 enabled.... init_glm_aed2 processing: ', TRIM(fname)
   namlst = f_get_lun()

   write(*,"(/,5X,'---------- AED2 config : start ----------')")
   IF ( aed2_init_core('.') /= 0 ) STOP "     ERROR: Initialisation of aed2_core failed"
   CALL aed2_print_version

   !# Create model tree
   print *,"     Processing aed2_models config from ",TRIM(fname)
   OPEN(namlst,file=fname,action='read',status='old',iostat=status)
   IF ( status /= 0 ) CALL STOPIT("Cannot open file " // TRIM(fname))

   models = ''
   READ(namlst, nml=aed2_models, iostat=status)
   IF ( status /= 0 ) STOP "Cannot read namelist entry aed2_models"

   DO i=1,size(models)
      IF (models(i)=='') EXIT
      CALL aed2_define_model(models(i), namlst)
   ENDDO

   !# should be finished with this file
   CLOSE(namlst)
   print *,"      ... nml file parsing completed."

   n_aed2_vars = aed2_core_status(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

#if DEBUG
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         print *,"AED2 var ", i, tvar%sheet, tvar%diag, tvar%extern, TRIM(tvar%name)
      ELSE
         print *,"AED2 var ", i, " is empty"
      ENDIF
   ENDDO
#endif

   print "(/,5X,'AED2 : n_aed2_vars = ',I3,' ; MaxLayers         = ',I4)",n_aed2_vars,MaxLayers
   print "(  5X,'AED2 : n_vars      = ',I3,' ; n_vars_ben        = ',I3)",n_vars,n_vars_ben
   print "(  5X,'AED2 : n_vars_diag = ',I3,' ; n_vars_diag_sheet = ',I3,/)",n_vars_diag,n_vars_diag_sheet

   CALL check_data

   !# names = grab the names from info
   ALLOCATE(names(n_vars),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (names)'
   ALLOCATE(bennames(n_vars_ben),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (bennames)'

#if PLOTS
   ALLOCATE(plot_id_v(n_vars))
   ALLOCATE(plot_id_sv(n_vars_ben))
   ALLOCATE(plot_id_d(n_vars_diag))
   ALLOCATE(plot_id_sd(n_vars_diag_sheet))
   plot_id_v = -1; plot_id_sv = -1; plot_id_d = -1; plot_id_sd = -1
#endif

   NumWQ_Vars = n_vars
   NumWQ_Ben  = n_vars_ben
   !# Now that we know how many vars we need, we can allocate space for them
   ALLOCATE(cc(MaxLayers, (n_vars + n_vars_ben)),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (CC)'
   cc = 0.         !# initialise to zero

   CALL set_c_wqvars_ptr(cc)

   ALLOCATE(min_((n_vars + n_vars_ben))) ; ALLOCATE(max_((n_vars + n_vars_ben)))
   print "(5X,'Configured AED2 variables to simulate:')"

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%sheet .OR. tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            names(j) = TRIM(tvar%name)
            min_(j) = tvar%minimum
            max_(j) = tvar%maximum
            !print *,"     S(",j,") AED2 pelagic(3D) variable: ", TRIM(names(j))
            print "(7X,'S(',I4,') water column variable     : ',A)",j , TRIM(names(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%sheet .AND. .NOT. (tvar%diag .OR. tvar%extern) ) THEN
            j = j + 1
            bennames(j) = TRIM(tvar%name)
            min_(n_vars+j) = tvar%minimum
            max_(n_vars+j) = tvar%maximum
            !print *,"     B(",j,") AED2 benthic(2D) variable: ", TRIM(bennames(j))
            print "(7X,'B(',I4,') bottom variable           + ',A)",j , TRIM(bennames(j))
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            IF ( .NOT.  tvar%sheet ) THEN
               j = j + 1
               print "(7X,'D(',I4,') water column diagnostic   > ',A)",j , TRIM(tvar%name)
               !print *,"     D(",j,") AED2 diagnostic 3Dvariable: ", TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   j = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( tvar%diag ) THEN
            IF (tvar%sheet ) THEN
               j = j + 1
               !print *,"     D(",j,") AED2 diagnostic 2Dvariable: ", TRIM(tvar%name)
               print "(7X,'D(',I4,') bottom/surface diagnostic ~ ',A)",j , TRIM(tvar%name)
            ENDIF
         ENDIF
      ENDIF
   ENDDO

!  ALLOCATE(column(n_aed2_vars))
   ALLOCATE(externalid(n_aed2_vars))

   !----------------------------------------------------------------------------

   !# Now set initial values
   v = 0 ; sv = 0;
   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "     ERROR getting variable info"
      IF ( .NOT. ( tvar%extern .OR. tvar%diag) ) THEN  !# neither global nor diagnostic variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            cc(:, n_vars+sv) = tvar%initial
         ELSE
            v = v + 1
            cc(:, v) = tvar%initial
         ENDIF
      ENDIF
   ENDDO

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag(MaxLayers, n_vars_diag),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = zero_

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag_hz(n_vars_diag_sheet),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = zero_

   !# Allocate array with vertical movement rates (m/s, positive for upwards),
   !# and set these to the values provided by the model.
   !# allocated for all vars even though only state vars entries will be used
   ALLOCATE(ws(MaxLayers, n_aed2_vars),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (ws)'
   ws = zero_

   !# Allocate array for photosynthetically active radiation (PAR).
   !# This will be calculated internally during each time step.
   ALLOCATE(par(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = zero_

   ALLOCATE(nir(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (nir)'
   nir = zero_
   ALLOCATE(uva(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uva)'
   uva = zero_
   ALLOCATE(uvb(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (uvb)'
   uvb = zero_

   ALLOCATE(dz(MaxLayers),stat=rc)
   dz = zero_

   !# Allocate array for local pressure.
   !# This will be calculated [approximated] from layer depths internally
   !# during each time step.
   ALLOCATE(pres(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   pres = zero_

   ALLOCATE(tss(MaxLayers),stat=rc)
   IF (rc /= 0) STOP 'allocate_memory(): Error allocating (tss)'
   tss = zero_

   write(*,"(/,5X,'----------  AED2 config : end  ----------',/)")

END SUBROUTINE aed2_init_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed2_is_var(id,i_vname,len) BIND(C, name=_WQ_IS_VAR)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)   :: id
   CCHARACTER,INTENT(in) :: i_vname(*)
   CSIZET,INTENT(in)     :: len
!LOCALS
   CHARACTER(len=45) :: vname
   TYPE(aed2_variable_t),POINTER :: tvar
   INTEGER :: i, v, sv, d, sd
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(vname, i_vname, len)

   v = 0; sv = 0; d = 0; sd = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%diag .OR. tvar%extern) ) THEN
            IF ( tvar%sheet ) THEN ; sv=sv+1; ELSE ; v=v+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  IF (benthic_mode .GT. 1) THEN
                     aed2_is_var=sv
                  ELSE
                     aed2_is_var=-sv
                  ENDIF
#ifdef PLOTS
                  plot_id_sv(sv) = id;
#endif
               ELSE
                  aed2_is_var=v
#ifdef PLOTS
                  plot_id_v(v) = id;
#endif
               ENDIF
               RETURN
            ENDIF
         ELSEIF ( tvar%diag ) THEN
            IF ( tvar%sheet ) THEN ; sd=sd+1; ELSE ; d=d+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  aed2_is_var=-sd
#ifdef PLOTS
                  plot_id_sd(sd) = id;
#endif
               ELSE
                  aed2_is_var=d
#ifdef PLOTS
                  plot_id_d(d) = id;
#endif
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   aed2_is_var = 0
END FUNCTION aed2_is_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_set_glm_data(Lake, MaxLayers, MetData, SurfData, dt_,          &
                                c_rain_factor, c_sw_factor, c_friction)        &
                                                  BIND(C, name=_WQ_SET_GLM_DATA)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER, INTENT(in) :: MaxLayers
   TYPE(C_PTR),VALUE :: Lake
   TYPE(MetDataType),TARGET     :: MetData  !# Meteorological data
   TYPE(SurfaceDataType),TARGET :: SurfData !# Surface Data
   AED_REAL,INTENT(in)  :: dt_
!  AED_REAL,OPTIONAL,TARGET :: c_rain_factor, c_sw_factor, c_friction
   AED_REAL,TARGET :: c_rain_factor, c_sw_factor, c_friction
!LOCALS
   INTEGER :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(Lake, theLake, [MaxLayers])

   !# Save pointers to external dynamic variables that we need later (in do_glm_wq)
   z    => theLake%Height
   temp => theLake%Temp
   salt => theLake%Salinity
   rho  => theLake%Density
   area => theLake%LayerArea
   rad  => theLake%Light
   extc_coef => theLake%ExtcCoefSW
   layer_stress => theLake%LayerStress

   IF (benthic_mode .GT. 1) zz => z
   ALLOCATE(depth(MaxLayers))
   ALLOCATE(layer_area(MaxLayers))
   ALLOCATE(sed_zones(MaxLayers))
   sed_zones = 0.

   IF (benthic_mode .GT. 1) THEN
      ALLOCATE(z_diag(n_zones, n_vars_diag)) ; z_diag = zero_
      ALLOCATE(z_diag_hz(n_zones, n_vars_diag_sheet)) ; z_diag_hz = zero_

      !# At this point we have z_cc allocated, so now we can copy the initial values
      !# from cc benthic vars to it
      DO i=1,n_zones
         z_cc(i, n_vars+1:n_vars+n_vars_ben) = cc(1, n_vars+1:n_vars+n_vars_ben)
      ENDDO
   ENDIF

   precip => MetData%Rain
   evap   => SurfData%Evap
   bottom_stress => layer_stress(botmLayer)

   !# Copy scalars that will not change during simulation, and are needed in do_glm_wq)
   dt = dt_

   !# Provide pointers to arrays with environmental variables to aed2.
   wnd => MetData%WindSpeed
   I_0 => MetData%ShortWave

   !# Calculate and save internal time step.
   dt_eff = dt/FLOAT(split_factor)

!  IF (PRESENT(c_rain_factor)) rain_factor => c_rain_factor
!  IF (PRESENT(c_sw_factor))   sw_factor   => c_sw_factor
!  IF (PRESENT(c_friction))    friction    => c_friction
   rain_factor => c_rain_factor
   sw_factor   => c_sw_factor
   friction    => c_friction

   !# Trigger an error if WQ hasn't got all it needs from us.
   CALL check_data
END SUBROUTINE aed2_set_glm_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_data
!-------------------------------------------------------------------------------
! Check that all variable dependencies have been met
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: av !, i
   INTEGER :: v, d, sv, sd, ev, err_count
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   err_count = 0

   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; tvar%found = .true.
            CASE ( 'salinity' )    ; tvar%found = .true.
            CASE ( 'density' )     ; tvar%found = .true.
            CASE ( 'layer_ht' )    ; tvar%found = .true.
            CASE ( 'extc_coef' )   ; tvar%found = .true.
            CASE ( 'tss' )         ; tvar%found = .true.
            CASE ( 'par' )         ; tvar%found = .true.
            CASE ( 'nir' )         ; tvar%found = .true.
            CASE ( 'uva' )         ; tvar%found = .true.
            CASE ( 'uvb' )         ; tvar%found = .true.
            CASE ( 'pressure' )    ; tvar%found = .true.
            CASE ( 'depth' )       ; tvar%found = .true.
            CASE ( 'sed_zone' )    ; tvar%found = .true.
            CASE ( 'wind_speed' )  ; tvar%found = .true.
            CASE ( 'par_sf' )      ; tvar%found = .true.
            CASE ( 'taub' )        ; tvar%found = .true.
            CASE ( 'lake_depth' )  ; tvar%found = .true.
            CASE ( 'layer_area' )  ; tvar%found = .true.
            CASE ( 'rain' )        ; tvar%found = .true.
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
         ELSE
            d = d + 1
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
         ELSE
            v = v + 1
         ENDIF
      ENDIF
      IF ( .NOT. tvar%found ) THEN
         print *, "ERROR: Undefined variable ", TRIM(tvar%name)
         err_count = err_count + 1
      ENDIF
   ENDDO

   IF ( n_vars < v ) print *,"More vars than expected",v,n_vars
   IF ( n_vars_ben < sv ) print *,"More sheet vars than expected"
   IF ( n_vars_diag < d ) print *,"More diag vars than expected"
   IF ( n_vars_diag_sheet < sd ) print *,"More sheet diag vars than expected"

   IF ( err_count > 0 ) CALL STOPIT("*** Errors in configuration")
END SUBROUTINE check_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_sed_column(column, top, flux_pel, flux_atm, flux_ben)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in)  :: top
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_layers, n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av !, i
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => theZones%ztemp
            CASE ( 'salinity' )    ; column(av)%cell => theZones%zsalt
            CASE ( 'density' )     ; column(av)%cell => theZones%zrho
            CASE ( 'layer_ht' )    ; column(av)%cell => theZones%zdz
            CASE ( 'extc_coef' )   ; column(av)%cell => theZones%zextc_coef
            CASE ( 'tss' )         ; column(av)%cell => theZones%ztss
            CASE ( 'par' )         ; column(av)%cell => theZones%zpar
            CASE ( 'nir' )         ; column(av)%cell => theZones%znir
            CASE ( 'uva' )         ; column(av)%cell => theZones%zuva
            CASE ( 'uvb' )         ; column(av)%cell => theZones%zuvb
            CASE ( 'pressure' )    ; column(av)%cell => theZones%zpres
            CASE ( 'depth' )       ; column(av)%cell => theZones%zdepth
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => theZones(1)%z_sed_zones; zone_var = av
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0
            CASE ( 'taub' )        ; column(av)%cell_sheet => bottom_stress
            CASE ( 'lake_depth' )  ; column(av)%cell_sheet => depth(1)
            CASE ( 'layer_area' )  ; column(av)%cell => theZones%zarea
            CASE ( 'rain' )        ; column(av)%cell_sheet => precip
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//trim(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column(av)%cell_sheet => z_diag_hz(1,sd)
         ELSE
            d = d + 1
            column(av)%cell => z_diag(:, d)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => z_cc(1, n_vars+sv)
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => z_cc(top, n_vars+sv)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => z_cc(:, v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(:, v)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_sed_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE define_column(column, top, cc, cc_diag, flux_pel, flux_atm, flux_ben)
!-------------------------------------------------------------------------------
! Set up the current column pointers
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER, INTENT(in)  :: top
   AED_REAL, TARGET, INTENT(in) :: cc(:,:)       !# (n_layers, n_vars)
   AED_REAL, TARGET, INTENT(in) :: cc_diag(:,:)  !# (n_layers, n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_pel(:,:) !# (n_layers, n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_atm(:)   !# (n_vars)
   AED_REAL, TARGET, INTENT(inout) :: flux_ben(:)   !# (n_vars)
!
!LOCALS
   INTEGER :: av !, i
   INTEGER :: v, d, sv, sd, ev
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   v = 0 ; d = 0; sv = 0; sd = 0 ; ev = 0
   DO av=1,n_aed2_vars
      IF ( .NOT.  aed2_get_var(av, tvar) ) STOP "Error getting variable info"

      IF ( tvar%extern ) THEN !# global variable
         ev = ev + 1
         SELECT CASE (tvar%name)
            CASE ( 'temperature' ) ; column(av)%cell => temp(:)
            CASE ( 'salinity' )    ; column(av)%cell => salt(:)
            CASE ( 'density' )     ; column(av)%cell => rho(:)
            CASE ( 'layer_ht' )    ; column(av)%cell => dz(:)
            CASE ( 'extc_coef' )   ; column(av)%cell => extc_coef(:)
            CASE ( 'tss' )         ; column(av)%cell => tss(:)
            CASE ( 'par' )         ; column(av)%cell => par(:)
            CASE ( 'nir' )         ; column(av)%cell => nir(:)
            CASE ( 'uva' )         ; column(av)%cell => uva(:)
            CASE ( 'uvb' )         ; column(av)%cell => uvb(:)
            CASE ( 'pressure' )    ; column(av)%cell => pres(:)
            CASE ( 'depth' )       ; column(av)%cell => depth(:)
            CASE ( 'sed_zone' )    ; column(av)%cell_sheet => sed_zones(1)
            CASE ( 'wind_speed' )  ; column(av)%cell_sheet => wnd
            CASE ( 'par_sf' )      ; column(av)%cell_sheet => I_0
            CASE ( 'taub' )        ; column(av)%cell_sheet => bottom_stress
            CASE ( 'lake_depth' )  ; column(av)%cell_sheet => depth(1)
            CASE ( 'layer_area' )  ; column(av)%cell => layer_area(:)
            CASE ( 'rain' )        ; column(av)%cell_sheet => precip
            CASE DEFAULT ; CALL STOPIT("ERROR: external variable "//TRIM(tvar%name)//" not found.")
         END SELECT
      ELSEIF ( tvar%diag ) THEN  !# Diagnostic variable
         IF ( tvar%sheet ) THEN
            sd = sd + 1
            column(av)%cell_sheet => cc_diag_hz(sd)
         ELSE
            d = d + 1
            column(av)%cell => cc_diag(:, d)
         ENDIF
      ELSE    !# state variable
         IF ( tvar%sheet ) THEN
            sv = sv + 1
            IF ( tvar%bot ) THEN
               column(av)%cell_sheet => cc(1, n_vars+sv)
!            print *,'av',av,sv
            ELSEIF ( tvar%top ) THEN
               column(av)%cell_sheet => cc(top, n_vars+sv)
            ENDIF
            column(av)%flux_ben => flux_ben(n_vars+sv)
            column(av)%flux_atm => flux_atm(n_vars+sv)
         ELSE
            v = v + 1
            column(av)%cell => cc(:, v)
            column(av)%flux_atm => flux_atm(v)
            column(av)%flux_pel => flux_pel(:, v)
            column(av)%flux_ben => flux_ben(v)
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE define_column
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calculate_fluxes(column, wlev, column_sed, nsed, flux_pel, flux_atm, flux_ben, flux_zon)
!-------------------------------------------------------------------------------
! Checks the current values of all state variables and repairs these
!-------------------------------------------------------------------------------
#ifdef HAVE_IEEE_ARITH
!USES
   USE IEEE_ARITHMETIC
#endif
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   TYPE (aed2_column_t), INTENT(inout) :: column_sed(:)
   INTEGER, INTENT(in) :: wlev, nsed
   AED_REAL, INTENT(inout) :: flux_pel(:,:) !# (wlev, n_vars+n_vars_ben)
   AED_REAL, INTENT(inout) :: flux_atm(:)   !# (n_vars+n_vars_ben)
   AED_REAL, INTENT(inout) :: flux_ben(:)   !# (n_vars+n_vars_ben)
   AED_REAL, INTENT(inout) :: flux_zon(:,:) !# (n_zones)
!
!LOCALS
   INTEGER :: lev,zon,v_start,v_end,av,sv,sd
   AED_REAL :: scale
   AED_REAL, DIMENSION(wlev, n_vars+n_vars_ben)    :: flux_pel_pre
   AED_REAL, DIMENSION(n_zones, n_vars+n_vars_ben) :: flux_pel_z
   AED_REAL :: localrainl, localshade, localdrag
   LOGICAL :: splitZone
   TYPE(aed2_variable_t),POINTER :: tvar
!-------------------------------------------------------------------------------
!BEGIN
   flux_pel = zero_
   flux_atm = zero_
   flux_ben = zero_
   flux_zon = zero_

   !# Start with calculating all flux terms for rhs in mass/m3/s
   !# Includes (1) benthic flux, (2) surface exchange and (3) water column kinetics
   !# as calculated by glm


   !# (1) BENTHIC FLUXES
   IF ( benthic_mode .GT. 1 ) THEN
      !# Multiple static sediment zones are simulated, and therfore overlying
      !# water conditions need to be aggregated from multiple cells/layers, and output flux
      !# needs disaggregating from each zone back to the overlying cells/layers

!$OMP DO
      DO zon=1,nsed
         !# Reinitialise flux_ben to be repopulated for this zone
         flux_ben = zero_
         flux_pel_pre = zero_

         !# If multiple benthic zones, we must update the benthic variable pointer for the new zone
         IF ( zone_var .GE. 1 ) THEN
            column_sed(zone_var)%cell_sheet => theZones(zon)%z_sed_zones
    !       !MH WE NEED A COLUMN TO CC VAR MAP FOR BENTHIC GUYS
            !CAB Yes, a map (or 2 maps) would be better, but QnD since this all needs reworking
            sv = 0 ; sd = 0
            DO av=1,n_aed2_vars
               IF ( .NOT. aed2_get_var(av, tvar) ) STOP "Error getting variable info"
               IF ( .NOT. tvar%extern .AND. tvar%sheet ) THEN
                  IF ( tvar%diag ) THEN
                     sd = sd + 1
                     column(av)%cell_sheet => z_diag_hz(zon, sd)
                  ELSE
                     sv = sv + 1
                     column(av)%cell_sheet => z_cc(zon, n_vars+sv)
                  ENDIF
               ENDIF
            ENDDO
            !print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
         ENDIF
         IF ( benthic_mode .EQ. 3 ) THEN
            !# Zone is able to operated on by riparian and dry methods
            CALL aed2_calculate_riparian(column_sed, zon, theZones(zon)%z_pc_wet)
            IF (theZones(zon)%z_pc_wet .EQ. 0. ) CALL aed2_calculate_dry(column_sed, zon)

            !# update feedback arrays to host model, to reduce rain (or if -ve then add flow)
            CALL aed2_rain_loss(column, 1, localrainl);
            IF (link_rain_loss) rain_factor = localrainl

            !# update feedback arrays to shade the water (ie reduce incoming light, Io)
            CALL aed2_light_shading(column, 1, localshade)
            IF (link_solar_shade) sw_factor = localshade

            !# now the bgc updates are complete, update links to host model
            CALL aed2_bio_drag(column, 1, localdrag)
            IF (link_bottom_drag) friction = localdrag
         ENDIF
         !# Calculate temporal derivatives due to benthic processes.
         !# They are stored in flux_ben (benthic vars) and flux_pel (water vars)
         flux_pel_pre = flux_pel

!        print*,"Calling ben for zone ",zone_var,zon,z_sed_zones(zon)
         CALL aed2_calculate_benthic(column_sed, zon)

         !# Record benthic fluxes in the zone array
         flux_zon(zon, :) = flux_ben(:)

         !# Now we have to find out the water column flux that occured and
         !# disaggregate it to relevant layers
         flux_pel_z(zon,:) = flux_pel(zon,:)-flux_pel_pre(zon,:)
      ENDDO
!$OMP END DO

      !# Disaggregation of zone induced fluxes to overlying layers
      v_start = 1 ; v_end = n_vars
      zon = n_zones
      DO lev=wlev,1,-1
        IF ( zon .GT. 1 ) THEN
          IF (lev .GT. 1) THEN
            splitZone = zz(lev-1) < zone_heights(zon-1)
          ELSE
            splitZone = 0.0 < zone_heights(zon-1)
          ENDIF
        ELSE
          splitZone = .FALSE.
        ENDIF

        IF (splitZone) THEN
          IF (lev .GT. 1) THEN
            scale = (zone_heights(zon-1) - zz(lev-1)) / (zz(lev) - zz(lev-1))
          ELSE
            scale = (zone_heights(zon-1) - 0.0) / (zz(lev) - 0.0)
          ENDIF
          flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end) * scale

          zon = zon - 1

          flux_pel(lev,v_start:v_end) = flux_pel(lev,v_start:v_end) + &
                                        flux_pel_z(zon,v_start:v_end) * (1.0 - scale)
        ELSE
          flux_pel(lev,v_start:v_end) = flux_pel_z(zon,v_start:v_end)
        ENDIF
      ENDDO
      !# Limit flux out of bottom waters to concentration of that layer
      !# i.e. don't flux out more than is there & distribute
      !# bottom flux into pelagic over bottom box (i.e., divide by layer height).
      !# scaled to proportion of area that is "bottom"
      DO lev=1,wlev
         if(lev>1)flux_pel(lev, :) = flux_pel(lev, :) * (area(lev)-area(lev-1))/area(lev)
         flux_pel(lev, :) = max(-1.0 * cc(lev, :), flux_pel(lev, :)/dz(lev))
      ENDDO
   ELSE
      !# Sediment zones are not simulated and therefore just operate on the bottom-most
      !# GLM layer as the "benthos". If benthic_mode=1 then benthic fluxes will also be
      !# applied on flanks of the remaining layers, but note this is not suitable for
      !# model configurations where mass balance of benthic variables is required.

      !# Calculate temporal derivatives due to exchanges at the sediment/water interface
      IF ( zone_var .GE. 1 ) column(zone_var)%cell_sheet => theZones(1)%z_sed_zones
      CALL aed2_calculate_benthic(column, 1)

      !# Limit flux out of bottom layers to concentration of that layer
      !# i.e. don't flux out more than is there
      !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
      flux_pel(1, :) = max(-1.0 * cc(1, :), flux_pel(1, :)/dz(1))

      IF ( benthic_mode .EQ. 1 ) THEN
!$OMP DO
         DO lev=2,wlev
            !# Calculate temporal derivatives due to benthic fluxes.
            CALL aed2_calculate_benthic(column, lev)

            !# Limit flux out of bottom layers to concentration of that layer
            !# i.e. don't flux out more than is there
            !# & distribute bottom flux into pelagic over bottom box (i.e., divide by layer height).
            !# scaled to proportion of area that is "bottom"
            flux_pel(lev, :) = max(-1.0 * cc(lev, :), flux_pel(lev, :)/dz(lev))
            flux_pel(lev, :) = flux_pel(lev, :) * (area(lev)-area(lev-1))/area(lev)
         ENDDO
!$OMP END DO
      ENDIF
   ENDIF

   !# (2) SURFACE FLUXES
   !# Calculate temporal derivatives due to air-water exchange.
   IF (.NOT. lIce) THEN !# no surface exchange under ice cover
      CALL aed2_calculate_surface(column, wlev)

      !# Distribute the fluxes into pelagic surface layer
      flux_pel(wlev, :) = flux_pel(wlev, :) + flux_atm(:)/dz(wlev)
   ENDIF

   !# (3) WATER COLUMN KINETICS
   !# Add pelagic sink and source terms for all depth levels.
   DO lev=1,wlev
      CALL aed2_calculate(column, lev)
   ENDDO
END SUBROUTINE calculate_fluxes
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE check_states(column, wlev)
!-------------------------------------------------------------------------------
#ifdef HAVE_IEEE_ARITH
!USES
   USE IEEE_ARITHMETIC
#endif
!
!ARGUMENTS
   TYPE (aed2_column_t),INTENT(inout) :: column(:)
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv
   INTEGER :: i,v,lev
#if DEBUG
   INTEGER :: last_naned
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
#if DEBUG
   last_naned = -1
#endif
   DO lev=1, wlev
      CALL aed2_equilibrate(column, lev)    !MH this should be in the main do_glm routine ????!!!
      v = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%diag .OR. tv%extern) ) THEN
               v = v + 1
               IF ( repair_state ) THEN
#if DEBUG
                  IF ( isnan(cc(lev, v)) ) last_naned = i
#endif
                  IF ( .NOT. isnan(min_(v)) ) THEN
                     IF ( cc(lev, v) < min_(v) ) cc(lev, v) = min_(v)
                  ENDIF
                  IF ( .NOT. isnan(max_(v)) ) THEN
                     IF ( cc(lev, v) > max_(v) ) cc(lev, v) = max_(v)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDDO

#if DEBUG
   IF ( last_naned > -1 ) THEN
      IF ( aed2_get_var(last_naned, tv) ) THEN
         print*,"NaNs detected in CC in var ", TRIM(tv%name)
      ELSE
         print*,"NaNs detected in CC unidentified var"
      ENDIF
      STOP
   ENDIF
#endif
END SUBROUTINE check_states
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_do_glm(wlev, pIce) BIND(C, name=_WQ_DO_GLM)
!-------------------------------------------------------------------------------
!                           wlev is the number of levels used;
!-------------------------------------------------------------------------------
#ifdef HAVE_IEEE_ARITH
!USES
   USE IEEE_ARITHMETIC
#endif
!
!ARGUMENTS
   CINTEGER,INTENT(in) :: wlev
   CLOGICAL,INTENT(in) :: pIce
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv

   AED_REAL :: min_C, surf
   INTEGER  :: i, j, v, lev, split

   TYPE (aed2_column_t) :: column(n_aed2_vars)
   TYPE (aed2_column_t) :: column_sed(n_aed2_vars)
   AED_REAL :: flux_ben(n_vars+n_vars_ben), flux_atm(n_vars+n_vars_ben)
   AED_REAL :: flux(wlev, n_vars+n_vars_ben)
   AED_REAL :: flux_zone(n_zones, n_vars+n_vars_ben)
!
!-------------------------------------------------------------------------------
!BEGIN
   lIce = pIce

   surf = z(wlev)
   !# re-compute the layer heights and depths
   dz(1) = z(1)
   depth(1) = surf - z(1)
   layer_area(1) = 1
   DO i=2,wlev
      dz(i) = z(i) - z(i-1)
      depth(i) = surf - z(i)
      layer_area(i) = (area(i)-area(i-1))/area(i)
   ENDDO

   IF ( benthic_mode .GT. 1 ) THEN
      j = 1
      DO i=1,wlev
        !print *,'j',i,j
        !print *,'zone_heights',z(i),zone_heights(j)!,theZones(1)%zheight,theZones(2)%zheight
         IF (z(i) .GT. zone_heights(j)) THEN
            sed_zones(i) = j * area(i)
            j = j+1
         ELSE
            sed_zones(i) = j
         ENDIF
      ENDDO
   ENDIF

   !# Calculate local pressure
   pres(1:wlev) = -z(1:wlev)

   CALL define_column(column, wlev, cc, cc_diag, flux, flux_atm, flux_ben)
   IF (benthic_mode .GT. 1) &
      CALL define_sed_column(column_sed, n_zones, flux, flux_atm, flux_ben)

   cc_diag = 0.
   cc_diag_hz = 0.

   IF ( .NOT. mobility_off ) THEN
     v = 0
     DO i=1,n_aed2_vars
        IF ( aed2_get_var(i, tv) ) THEN
           IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
              v = v + 1
              ws(:,i) = zero_
              ! only for state_vars that are not sheet
              IF ( .NOT. isnan(tv%mobility) ) THEN
                 ! default to ws that was set during initialisation
                 ws(1:wlev,i) = tv%mobility
                 IF(i == 14) print *,'ws',i,ws(1:wlev,i)
              ENDIF
           ENDIF
        ENDIF
     ENDDO
      DO i = 1, wlev
         ! update ws for modules that use the mobility method
         CALL aed2_mobility(column, i, ws(i,:))
      ENDDO

      !# (3) Calculate source/sink terms due to the settling or rising of
      !# state variables in the water column (note that settling into benthos
      !# is done in aed2_do_benthos)
      v = 0
      DO i=1,n_aed2_vars

         IF ( aed2_get_var(i, tv) ) THEN
            IF ( .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern)   ) THEN
               v = v + 1
               !# only for state_vars that are not sheet, and also non-zero ws
               IF ( .NOT. isnan(tv%mobility) .AND. SUM(ABS(ws(1:wlev,i)))>zero_ ) THEN
                  min_C = tv%minimum
                  CALL Mobility(wlev, dt, dz, area, ws(:, i), min_C, cc(:, v))
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF
   CALL check_states(column,wlev)

   DO split=1,split_factor

      IF (benthic_mode .GT. 1) THEN
         CALL copy_to_zone(cc, wlev)
         CALL calc_zone_areas(area, wlev, z(wlev))
      ENDIF

      !# Update local light field (self-shading may have changed through
      !# changes in biological state variables). Update_light is set to
      !# be inline with current aed2_phyoplankton, which requires only
      !# surface par, then integrates over depth of a layer
      CALL update_light(column, wlev)

      !# Fudge
      nir(:) = (par(:)/par_fraction) * nir_fraction
      uva(:) = (par(:)/par_fraction) * uva_fraction
      uvb(:) = (par(:)/par_fraction) * uvb_fraction

      !# Time-integrate one biological time step
      CALL calculate_fluxes(column, wlev, column_sed, n_zones,  &
                                  flux(:,:), flux_atm, flux_ben, flux_zone(:,:))
      !# Update the water column layers
      DO v = 1, n_vars
         DO lev = 1, wlev
            cc(lev, v) = cc(lev, v) + dt_eff*flux(lev, v)
         ENDDO
      ENDDO
      !# Now update benthic variables, depending on whether zones are simulated
      IF ( benthic_mode .GT. 1 ) THEN
         ! Loop through benthic state variables to update their mass
         DO v = n_vars+1, n_vars+n_vars_ben
            ! Loop through each sediment zone
            DO lev = 1, n_zones
               ! Update the main cc_sed data array with the
               z_cc(lev, v) = z_cc(lev, v)+ dt_eff*flux_zone(lev, v)
            ENDDO
         ENDDO
      ELSE
         DO v = n_vars+1, n_vars+n_vars_ben
            cc(1, v) = cc(1, v) + dt_eff*flux_ben(v)
         ENDDO
      ENDIF

      !# Distribute cc-sed benthic properties back into main cc array
      IF ( benthic_mode .GT. 1 ) &
         CALL copy_from_zone(cc, cc_diag, cc_diag_hz, wlev)

      CALL check_states(column, wlev)
   ENDDO
END SUBROUTINE aed2_do_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_clean_glm() BIND(C, name=_WQ_CLEAN_GLM)
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN
   ! Deallocate internal arrays
   IF (ALLOCATED(cc_diag))    DEALLOCATE(cc_diag)
   IF (ALLOCATED(cc_diag_hz)) DEALLOCATE(cc_diag_hz)
   IF (ALLOCATED(ws))         DEALLOCATE(ws)
   IF (ALLOCATED(total))      DEALLOCATE(total)
   IF (ALLOCATED(local))      DEALLOCATE(local)
   IF (ALLOCATED(par))        DEALLOCATE(par)
   IF (ALLOCATED(nir))        DEALLOCATE(nir)
   IF (ALLOCATED(uva))        DEALLOCATE(uva)
   IF (ALLOCATED(uvb))        DEALLOCATE(uvb)
   IF (ALLOCATED(pres))       DEALLOCATE(pres)
   IF (ALLOCATED(dz))         DEALLOCATE(dz)
END SUBROUTINE aed2_clean_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE update_light(column, nlev)
!-------------------------------------------------------------------------------
! Calculate photosynthetically active radiation over entire column based
! on surface radiation, attenuated based on background & biotic extinction
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE (aed2_column_t), INTENT(inout) :: column(:)
   INTEGER,INTENT(in)  :: nlev
!
!LOCALS
   INTEGER :: i
   AED_REAL :: localext, localext_up
!
!-------------------------------------------------------------------------------
!BEGIN

   localext = zero_; localext_up = zero_

   ! Surface Kd
   CALL aed2_light_extinction(column, nlev, localext)

   ! Surface PAR
   par(nlev) = par_fraction * rad(nlev) * EXP( -(lKw+localext)*1e-6*dz(nlev) )

   ! Now set the top of subsequent layers, down to the bottom
   DO i = (nlev-1),1,-1

      localext_up = localext
      CALL aed2_light_extinction(column, i, localext)

      par(i) = par(i+1) * EXP( -(lKw + localext_up) * dz(i+1) )

      IF (bioshade_feedback) extc_coef(i) = lKw + localext

   ENDDO

END SUBROUTINE update_light
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE aed2_init_glm_output(ncid,x_dim,y_dim,z_dim,zone_dim,time_dim) BIND(C, name=_WQ_INIT_GLM_OUTPUT)
!-------------------------------------------------------------------------------
!  Initialize the output by defining biogeochemical variables.
!-------------------------------------------------------------------------------
!
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid,x_dim,y_dim,z_dim,zone_dim,time_dim
!
!LOCALS
   INTEGER i !, v, d
   INTEGER dims(4)

   TYPE(aed2_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Put NetCDF library in define mode.
   CALL define_mode_on(ncid)

   !# Set up dimension indices for 3D (+ time) variables (longitude,latitude,depth,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = z_dim
   dims(4) = time_dim

!  v = 0; d = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tv) ) THEN
         IF ( .NOT. (tv%sheet .OR. tv%extern) ) THEN
            !# only for state and diag vars that are not sheet
            externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 4, dims(1:4))
            CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
         ENDIF
      ENDIF
   ENDDO

   IF ( n_zones .GT. 0 ) THEN
      !# Set up dimension indices for 3D (+ time) variables (longitude,latitude,zone,time).
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = zone_dim
      dims(4) = time_dim

!print*,x_dim,y_dim,zone_dim,time_dim
!     v = 0; d = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( tv%sheet .AND. .NOT. tv%extern ) THEN
               !# only for state and diag sheet vars
               externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 4, dims(1:4))
               CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
            ENDIF
         ENDIF
      ENDDO
   ELSE
      !# Set up dimension indices for 2D (+ time) variables (longitude,latitude,time).
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = time_dim

!     v = 0; d = 0
      DO i=1,n_aed2_vars
         IF ( aed2_get_var(i, tv) ) THEN
            IF ( tv%sheet .AND. .NOT. tv%extern ) THEN
               !# only for state and diag sheet vars
               externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 3, dims(1:3))
               CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   !# Take NetCDF library out of define mode (ready for storing data).
   CALL define_mode_off(ncid)
END SUBROUTINE aed2_init_glm_output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!-------------------------------------------------------------------------------
SUBROUTINE aed2_write_glm(ncid,wlev,nlev,lvl,point_nlevs) BIND(C, name=_WQ_WRITE_GLM_)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid, wlev, nlev
   CINTEGER,INTENT(in) :: lvl(*), point_nlevs
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv

   INTEGER  :: i, j, v, d, sv, sd
   AED_REAL :: val_out
   CLOGICAL :: last = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0; d = 0; sv = 0; sd = 0
   DO i=1,n_aed2_vars
      IF ( aed2_get_var(i, tv) ) THEN
         IF ( tv%diag ) THEN
            !# Process and store diagnostic variables.
            IF ( tv%sheet ) THEN
               sd = sd + 1
               !# Process and store diagnostic variables defined on horizontal slices of the domain.
               IF ( n_zones .GT. 0 ) THEN
                  CALL store_nc_array(ncid, externalid(i), XYNT_SHAPE, n_zones, n_zones, array=cc_diag(1:n_zones, sd))
               ELSE
                  CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc_diag_hz(sd))
               ENDIF
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_sd(sd).GE.0 ) CALL put_glm_val_s(plot_id_sd(sd),cc_diag_hz(sd))
!IF ( do_plots .AND. plot_id_d(sd).GE.0 ) print*,"PLOT ",d,plot_id_d(sd),cc_diag_hz(sd)
#endif
            ELSE
               d = d + 1
               !# Store diagnostic variable values defined on the full domain.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc_diag(:, d))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_d(d).GE.0 ) CALL put_glm_val(plot_id_d(d), cc_diag(1:wlev, d))
!IF ( do_plots .AND. plot_id_d(d).GE.0 ) print*,"PLOT ",d,plot_id_d(d),cc_diag(1:3, d)
#endif
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc_diag(lvl(j)+1, d)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ELSE IF ( .NOT. tv%extern ) THEN
            IF ( tv%sheet ) THEN
               sv = sv + 1
               !# Store benthic biogeochemical state variables.
               IF ( n_zones .GT. 0 ) THEN
                  CALL store_nc_array(ncid, externalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_cc(1:n_zones, n_vars+sv))
               ELSE
                  CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc(1, n_vars+sv))
               ENDIF
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_sv(sv).GE.0 ) THEN
                  IF (benthic_mode .GT. 1) THEN
                     CALL put_glm_val(plot_id_sv(sv), cc(1:n_zones, n_vars+sv))
                  ELSE
                     CALL put_glm_val_s(plot_id_sv(sv), cc(1, n_vars+sv))
                  ENDIF
               ENDIF
#endif
            ELSE
               v = v + 1
               !# Store pelagic biogeochemical state variables.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc(:, v))
#ifdef PLOTS
               IF ( do_plots .AND. plot_id_v(v).GE.0 ) CALL put_glm_val(plot_id_v(v), cc(1:wlev, v))
#endif
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc(lvl(j)+1, v)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE aed2_write_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed2_var_index_c(name, len) BIND(C, name=_WQ_VAR_INDEX_C)
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
   aed2_var_index_c = WQVar_Index(tn) - 1
END FUNCTION aed2_var_index_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
INTEGER FUNCTION WQVar_Index(name)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CHARACTER(len=*) :: name
!
!LOCALS
   TYPE(aed2_variable_t),POINTER :: tv
   INTEGER i,v
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0
   DO i=1, n_aed2_vars
      IF ( aed2_get_var(i, tv) .AND. .NOT. (tv%sheet .OR. tv%diag .OR. tv%extern) ) THEN
         v = v + 1
         IF ( name .EQ. tv%name ) THEN
            WQVar_Index = v
            RETURN
         ENDIF
      ENDIF
   ENDDO
   WQVar_Index = -1
END FUNCTION WQVar_Index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_aed2
