#!/bin/sh

# ========================================

#SBATCH --partition=high_mem           #Selecting “shared” Queue
#SBATCH --job-name="ROS_Subj17to24"              #Name of Jobs (displayed in squeue)
#SBATCH --nodes=2                     #No of nodes to run job #2
#SBATCH --ntasks-per-node=8           #No of cores to use per node #8
#SBATCH --time=23:00:00            #Maximum time for job to run
#SBATCH --mem=8G                      #Amount of memory per node #8G
#SBATCH --output=ROS_Subj17to24.%N.%j.out        #Output file for stdout (optional)
#SBATCH -e ROS_Subj17to24.%N.%j.err              #Error file

# 
#  Script for extracting subject-level timeseries from data-driven ROIs active during Switch>Stay
#  ====================================
# 
#  Danielle Kurtin, July 2022
#

#######  STEPS ############
#### 1.  Prepare tools

module load fsl
module load matlab/R2019b

# ===============================
#### 2.  Define data and give permissions

## NOTE: CHANGE THE STUDY TO THE ONE YOU WANT TO RUN
export AtlasDir=/users/psychology01/parallel_scratch/projects/ROS/WatershedROIs
export OutputDir=/users/psychology01/parallel_scratch/projects/ROS/WatershedROIs/Timecourses
export DataDir=/users/psychology01/parallel_scratch/projects/ROS/Study_2
#export DataDir=/users/psychology01/parallel_scratch/projects/ROS/Study_2/NewData

#SubjNum=(1 2 3 4 5 6 7 9 10 12 13 14 15 16) 
#SubjNum=(19 20 21 22 23 24 25 26 27)
SubjNum=(8)
Run=(1 2 3)


# ==============================
#### 3. Extract and move the data

for ((s=0;s<${#SubjNum[@]};++s)); do
for ((r=0;r<${#Run[@]};++r)); do
#for ((c=0;c<${#Cope[@]};++c)); do

	echo Running Subject ${SubjNum[$s]}, Run${Run[$r]} timecourse extraction  

################
# Uncomment below for subjs 1-16
################  
#   flirt -in "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data.nii.gz -ref $FSLDIR/data/standard/MNI152_T1_2mm -applyxfm -usesqform -out "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data_2mm_v2.nii.gz

#   fslmeants -i "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data_2mm_v2.nii.gz --label="$AtlasDir"/FWS_S2_M2_c22_forLEiDA.nii -o "$OutputDir"/SwitchStayROIs/S${SubjNum[$s]}_R${Run[$r]}.csv

#   fslmeants -i "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data_2mm_v2.nii.gz --label=/users/psychology01/software/standard/Schaefer/MNI/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz  -o "$OutputDir"/SchaeferROIs/S${SubjNum[$s]}_R${Run[$r]}.csv


################
# Uncomment below for subjs 17-24
################  

# Get the transform from standard to subject
convert_xfm -inverse "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/reg/example_func2standard.mat -omat "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/reg/std2func.mat

# Split the 4D timeseries to the first volume so FLIRT works
fslsplit "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data.nii.gz "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/Split -t

# Flirt the atlas to subject space
flirt -in "$AtlasDir"/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii.gz -ref "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/Split0038 -applyxfm -init "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/reg/std2func.mat -out "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/TransAtlas

# Extract timeseries
fslmeants -i "$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/filtered_func_data.nii.gz --label="$DataDir"/ROS_DEPRES_${SubjNum[$s]}/1stlevel/Model1/run${Run[$r]}.feat/TransAtlas -o "$OutputDir"/SchaeferROIs/S${SubjNum[$s]}_R${Run[$r]}.csv


done
done
#done
