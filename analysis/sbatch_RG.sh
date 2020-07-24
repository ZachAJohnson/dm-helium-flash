#!/bin/bash
# Call the job this name
#SBATCH --job-name MESA_RG
# Number of nodes you want
#SBATCH --nodes=1
# Number of cpus you want
#SBATCH --cpus-per-task=5
# Memory per cpu
#SBATCH --mem-per-cpu=5g
# Job time
#SBATCH --time=19:00:00
# Always do this unless you get an account
#SBATCH --account=bsafdi1
# standard, or do debug if job less than 4 hours and only single-core
#SBATCH --partition=standard
# Don't email
#SBATCH --mail-type=NONE
# Can delete the following if not job array.
# If job array, the array is indexed by bash variable ${SLURM_ARRAY_TASK_ID}
#SBATCH --array=1-23
# Now do what you actually want to do
# Load modules

# Run code

# Setup
export OMP_NUM_THREADS=5
cd /scratch/bsafdi_root/bsafdi/dessert/MESA_models/RG/null

# Read from file
data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /nfs/turbo/bsafdi/dessert/dm-helium-flash/data/MESA_RG/RG_masses.dat)
initial_mass=${data}
# initial_mass=1.3
echo "Running a MESA simulation with initial mass ${initial_mass} solar. The task ID is ${SLURM_ARRAY_TASK_ID}."

star_folder=RG${SLURM_ARRAY_TASK_ID}_no_injection
cp -r $MESA_DIR/star/work $star_folder
# cp -r $MESA_DIR/star/test_suite/he_core_flash $star_folder # Could check how this looks, to make sure we're using everything we need. Although we have imported from there to below.
cd $star_folder
./mk

# Create project inlist
rm -r inlist_project

echo "! inlist to evolve a ${initial_mass} solar mass star

! For the sake of future readers of this file (yourself included),
! ONLY include the controls you are actually using.  DO NOT include
! all of the other controls that simply have their default values.

&star_job

  ! hack starting from pre-ms
    create_pre_main_sequence_model = .true.

  ! save a model at the end of the run
    save_model_when_terminate = .true.
    save_model_filename = '${star_folder}.mod'

  ! display on-screen plots
    pgstar_flag = .false.

/ !end of star_job namelist


&controls

    use_gold_tolerances = .true.
    use_eosDT2 = .true.
    use_eosELM = .true.
    
    use_eps_mdot = .false.

    use_dedt_form_of_energy_eqn = .true.
    min_cell_energy_fraction_for_dedt_form = 0
    use_eps_correction_for_KE_plus_PE_in_dLdm_eqn = .true.

    convergence_ignore_equL_residuals = .true. ! needed during flash
    ! note: can get rid of this if turn on conv_vel's.  just takes more steps.

  ! starting specifications
    initial_mass = ${initial_mass} ! in Msun units
    initial_y = 0.28 ! What about uncertainties on this?
    initial_z = 0.02

  ! use C/O enhanced opacities
  ! important for He-burning onwards
    use_Type2_opacities = .true.
    Zbase = 0.02

  use_ledoux_criterion = .true.

  cool_wind_RGB_scheme = 'Reimers'
  cool_wind_AGB_scheme = 'Blocker'
  RGB_to_AGB_wind_switch = 0.75d0
  Reimers_scaling_factor = 0.1d0  
  Blocker_scaling_factor = 0.2d0  

  ! same parameters
  mixing_length_alpha = 1.8

  ! output parameters
      photo_interval = 100
      profile_interval = 10
      history_interval = 10
      terminal_interval = -1
      write_header_frequency = -1
      max_num_profile_models = -1

  ! accuracy
  !mesh_delta_coeff = 0.8d0
  !mesh_dlog_3alf_dlogP_extra = 0.20d0
  !delta_lgRho_limit = 0.1d0
  !delta_lgT_limit = 0.1d0
  !min_timestep_limit = 1d-6 ! (seconds)

  ! Inject extra erg/s, only correct for 1.3 solar
  !start_time_for_inject_extra_ergs_sec = 1.39E+17 ! s
  !base_of_inject_extra_ergs_sec = 0 !9.09E-02 ! Msun
  !total_mass_for_inject_extra_ergs_sec = 6.10E-04 ! Msun
  !duration_for_inject_extra_ergs_sec = 2.89E+09 ! s
  !inject_extra_ergs_sec = 9.94E+${SLURM_ARRAY_TASK_ID} ! ergs/s ! 46 total

  ! extra_power_source = 0.5

  !inject_uniform_extra_heat = 1E+5
  !min_q_for_uniform_extra_heat = 0
  !max_q_for_uniform_extra_heat = 1E-01

  ! when to stop
    HB_limit = 1d-4

/ ! end of controls namelist" > inlist_project


./rn

