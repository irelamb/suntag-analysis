#!/bin/bash

# to be modified by the user #
input_path="..." # path to input data (tif/ometif)
env_path="..." # path to python environment
job_name="test"
#

model_name="gaussian_plane_v1"
n_files=$(($(ls ${input_path} | wc -l)))
echo $n_files

cp 02_template_ch0.sh 03_slurm_ch0.sh
sed -i "s#PATH_IN#$input_path#" 03_slurm_ch0.sh
sed -i "s#PATH_ENV#$env_path#" 03_slurm_ch0.sh
sed -i "s/N_FILES/$(($n_files-1))/" 03_slurm_ch0.sh
sed -i "s#JOB_NAME#$job_name#" 03_slurm_ch0.sh
sed -i "s#MODEL_NAME#$model_name#" 03_slurm_ch0.sh

cp 02_template_ch1.sh 03_slurm_ch1.sh
sed -i "s#PATH_IN#$input_path#" 03_slurm_ch1.sh
sed -i "s#PATH_ENV#$env_path#" 03_slurm_ch1.sh
sed -i "s/N_FILES/$(($n_files-1))/" 03_slurm_ch1.sh
sed -i "s#JOB_NAME#$job_name#" 03_slurm_ch1.sh
sed -i "s#MODEL_NAME#$model_name#" 03_slurm_ch1.sh

sbatch --qos=serial 03_slurm_ch0.sh
sbatch --qos=serial 03_slurm_ch1.sh
