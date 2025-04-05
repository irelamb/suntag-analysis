#!/bin/bash -l

#SBATCH --job-name=JOB_NAME
#SBATCH --output=/scratch/ilambert/SunTag/JOB_NAME_%A_%a.out
#SBATCH --error=/scratch/ilambert/SunTag/JOB_NAME_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 24
#SBATCH --mem 20G
#SBATCH --array=0-N_FILES

# In this version PATH_IN and N_FILES are just "placeholders" that will be substituted by the input path and the number of files in the input path, respectively, by a bash script. 

# ----------------- #
dx=0.2752 # microns
dy=0.2752
dt=20 # seconds
# ----------------- #

model_name="MODEL_NAME"

#chan0_avg=180 # average intensity green channel
#chan0_std=55 # fixed

#chan1_avg=CHAN1_AVG #200 # average intensity red channel
#chan1_std=55 # fixed

path="PATH_IN"
file_names=($(ls $path)) # file names of input images (either .tif or .ome.tif)

file_name=${file_names[${SLURM_ARRAY_TASK_ID}]}

source activate /home/ilambert/miniconda3/envs/microscopy_analysis 

#python3 main/makeFit_distributed.py $file_name "0" $chan0_avg $chan0_std $model_name $dx $dy $dt # channel 0
python3 main/makeFit_distributed.py $file_name "1" $model_name $dx $dy $dt # channel 1
