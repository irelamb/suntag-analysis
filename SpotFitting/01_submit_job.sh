#!/bin/bash

input_path="/scratch/ilambert/SunTag/Data/Images/20240725/COL-H-AG_01" #20230920/*
job_name="test"
model_name="gaussian_plane_v1"
#chan0_avg="140"
#chan1_avg="70"
#$(echo ${input_path} | cut -d "/" -f 8)
#n_files=$(($(ls ${input_path} | wc -l)*2))
n_files=$(($(ls ${input_path} | wc -l)))
echo $n_files

cp 02_template_ch0.sh 03_slurm_ch0.sh
sed -i "s#PATH_IN#$input_path#" 03_slurm_ch0.sh
sed -i "s/N_FILES/$(($n_files-1))/" 03_slurm_ch0.sh
sed -i "s#JOB_NAME#$job_name#" 03_slurm_ch0.sh
#sed -i "s#CHAN0_AVG#$chan0_avg#" 03_slurm_ch0.sh
sed -i "s#MODEL_NAME#$model_name#" 03_slurm_ch0.sh

cp 02_template_ch1.sh 03_slurm_ch1.sh
sed -i "s#PATH_IN#$input_path#" 03_slurm_ch1.sh
sed -i "s/N_FILES/$(($n_files-1))/" 03_slurm_ch1.sh
sed -i "s#JOB_NAME#$job_name#" 03_slurm_ch1.sh
#sed -i "s#CHAN1_AVG#$chan1_avg#" 03_slurm_ch1.sh
sed -i "s#MODEL_NAME#$model_name#" 03_slurm_ch1.sh

sbatch --qos=serial 03_slurm_ch0.sh
sbatch --qos=serial 03_slurm_ch1.sh
