#!/bin/bash

#SBATCH -A geomechanics
#SBATCH --job-name="Main_Periodic_Ocean"
#SBATCH --output="./RunOutputs/Main_Periodic_Ocean.%j"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00
#SBATCH --export=ALL
#SBATCH --constraint=skylake
#SBATCH --partition=expansion

#Step 1 Process
export OMP_NUM_THREADS=1
#srun Main_Periodic_Ocean 1 1 $1 $2;
srun Main_Periodic_Ocean 1 1 $1 $2 $3; #For seed based runs

#Step 2 Call python and then run post-processing
####source /groups/geomechanics/rigo/peri-ice/peri-ice/env/bin/activate
#python3 plotConfigurationOcean.py $1 Test_qv_nb_postrev_$1 $2
####python3 plotConfigurationOcean.py $1 Test_qv_nb_latmelt_$1 $2
##python3 plotFSDvsConcSlope2018v2.py
module load python/3.11.6-gcc-11.3.1-rphh4kv;
module load openmpi/5.0.1-gcc-11.3.1-j4o6ryt;
module load eigen/3.4.0-gcc-11.3.1-jvlgwff;
module load lsmlib/1.0.1;
source /groups/geomechanics/rigo/peri-ice/peri-ice/env/bin/activate;

python3 plotConfigurationOcean.py $1 Test_qv_nb_$1 $2 $3;

##python3 plotConfigurationOceanNP.py $1 Test_qv_nb_$1 $2
##python3 plotConfigurationOceanNP.py $1 Test_qv_nb_alterqv_$1 $2
#python3 plotConfigurationOcean_ML.py $1 Test_qv_nb_ML_alter_$1 $2
##python3 plotConfigurationOcean_F.py $1 Test_qv_nb_FF_$1 $2
##python3 plotConfigurationOcean_ML.py $1 Test_qv_nb_ML_alter_$1 $2
#python3 plotConfigurationOcean_ML.py $1 Test_qv_nb_ML_$1 $2
##python3 plotConfigurationOcean_Thickness.py $1 Test_qv_nb_Thickness_$1 $2

#make -f makefile_fluid;
