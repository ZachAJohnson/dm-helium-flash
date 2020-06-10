#!/bin/bash
# Call the job this name
#SBATCH --job-name MESA_RG
# Number of nodes you want
#SBATCH --nodes=1
# Number of cpus you want
#SBATCH --cpus-per-task=3
# Memory per cpu
#SBATCH --mem-per-cpu=5g
# Job time
#SBATCH --time=4:00:00
# Always do this unless you get an account
#SBATCH --account=bsafdi1
# standard, or do debug if job less than 4 hours and only single-core
#SBATCH --partition=standard
# Don't email
#SBATCH --mail-type=NONE
# Can delete the following if not job array.
# If job array, the array is indexed by bash variable ${SLURM_ARRAY_TASK_ID}
#SBATCH --array=19-23
# Now do what you actually want to do
# Load modules

# Run code

# Setup
export OMP_NUM_THREADS=5
cd /scratch/bsafdi_root/bsafdi/dessert/MESA_models/RG

# Read from file
data=$(sed "${SLURM_ARRAY_TASK_ID}q;d" /nfs/turbo/bsafdi/dessert/helium-flash/data/MESA_RG/RG_masses.dat)
initial_mass=${data}
echo "Running a MESA simulation with initial mass ${initial_mass} solar. The task ID is ${SLURM_ARRAY_TASK_ID}."

star_folder=RG${SLURM_ARRAY_TASK_ID}_v4
cp -r $MESA_DIR/star/work $star_folder
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

  ! starting specifications
    initial_mass = ${initial_mass} ! in Msun units
    initial_y = 0.28
    initial_z = 0.02

  ! use C/O enhanced opacities
  ! important for He-burning onwards
    use_Type2_opacities = .true.
    Zbase = 0.02

  ! configure mass loss on RGB & AGB
  !  cool_wind_RGB_scheme = 'Dutch'
  !  cool_wind_AGB_scheme = 'Dutch'
  !  RGB_to_AGB_wind_switch = 1d-4
  !  Dutch_scaling_factor = 0.8

  ! same parameters
  mixing_length_alpha = 1.8

  ! output parameters
      photo_interval = -1
      profile_interval = 10
      history_interval = 10
      terminal_interval = -1
      write_header_frequency = -1
      max_num_profile_models = 10000

  ! accuracy
  mesh_delta_coeff = 0.8d0
  mesh_dlog_3alf_dlogP_extra = 0.20d0
  delta_lgRho_limit = 0.1d0
  delta_lgT_limit = 0.1d0

  ! when to stop
  !  max_age = 11.75d9
    HB_limit = 1d-4

/ ! end of controls namelist" > inlist_project


./rn

