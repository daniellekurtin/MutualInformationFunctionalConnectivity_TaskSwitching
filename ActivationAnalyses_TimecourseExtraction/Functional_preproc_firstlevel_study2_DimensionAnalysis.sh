#!/bin/sh

Usage() {
  cat <<EOF

usage:  Functional_preproc_study2.sh -i <Subject>

e.g. sbatch Functional_preproc_study2.sh  -i Subject_1

EOF
  exit 1
}

if [  $# -le 1 ];
  then
  Usage
  exit 1
 fi

while [ $# -ge 1 ];
do
  case "$1" in

        -i)
                    subject=$2;
                    shift;;

  esac
  shift
done


echo " Input job = $subject "

# ========================================

#SBATCH --partition=shared            #Selecting “shared” Queue
#SBATCH --job-name="ROS"              #Name of Jobs (displayed in squeue)
#SBATCH --nodes=2                     #No of nodes to run job #2
#SBATCH --ntasks-per-node=8           #No of cores to use per node #8
#SBATCH --time=01-00:00:00            #Maximum time for job to run
#SBATCH --mem=15G                      #Amount of memory per node #8G
#SBATCH --output=ROS.%N.%j.out        #Output file for stdout (optional)
#SBATCH -e ROS.%N.%j.err              #Error file

echo ${subject}
# 
#  Script for pre-processing data ROS
#  ====================================
# 
#  Ines Violante, Sep 2021
#
#
#	Adapted by Danielle Kurtin, Sep 2021
#	NOTE: SOME SUBJECTS NEEDED TO BE RUN DIFFERENTLY. This is because not everyone's T1 was as high resolution as it needed to be. Those subjects were rerun using secondary or tertiary MPRAGE files to ensure their T1 was in order.

#######  STEPS ############
#### 1.  prepare tools

#source /users/psychology01/bin/conda_bash /users/psychology01/software/conda/imaging
module load fsl
module load matlab/R2019b


# ===============================
#### 2.  define data and give permissions

export MainDir=/users/dk00549/psychology01/parallel_scratch/projects/ROS/Study_2/NewData
export ScriptsDir=/users/psychology01/projects/ROS/Scripts


# # # # # # =========================================
# # # # # # Run GLM Models

export FSFdir="$ScriptsDir"/FSFfiles
export SubjectFSFdir="$FSFdir"/SubjectFSF 
Run=(run1 run2 run3)
Model=(DimensionAnalysis)


for ((j=0;j<${#Model[@]};++j)); do

        # If all feats need to be replaced 
        # echo ===== Removing old files ${subject} ${Model[$j]} =====
        # rm -rf "$MainDir"/${subject}/1stlevel/${Model[$j]}

   for ((k=0;k<${#Run[@]};++k)); do

			 if [ -f "$MainDir"/${subject}/preproc/${Run[$k]}.nii.gz ]      #make sure the niftis are there
			 then
			 	if [ ! -d "$MainDir"/${subject}/1stlevel/${Model[$j]}/${Run[$k]}.feat ]  #check if Feat's been run; if not, make the directory for it, find and replace the subject, model, and run variables in the template fsf
			 	then

			 		echo ==== Running FEAT ${subject} ${Model[$j]} ${Run[$k]} =====

			 		nvols="$(fslnvols "$MainDir"/${subject}/preproc/${Run[$k]}.nii.gz)"

                	subj=$(ls "$MainDir"/${subject}/Behaviour/)
              
					sed -e "s/ROS_DEPRES_19/${subject}/g" -e"s/set fmri(npts) 378/set fmri(npts) ${nvols}/g" -e"s/run1/${Run[$k]}/g" -e"s/DimensionAnalysis/${Model[$j]}/g" -e"s/Subject_19/${subj}/g" "$FSFdir"/NewStudy2_${Model[$j]}_run1.fsf > "$SubjectFSFdir"/${subject}_NewStudy2_${Model[$j]}_${Run[$k]}.fsf

			         feat "$SubjectFSFdir"/${subject}_NewStudy2_${Model[$j]}_${Run[$k]}.fsf
			 	else

					echo ===== ${subject}_NewStudy2_${Model[$j]}_${Run[$k]}.feat, already exists, skipping this step ====
				fi
			else
     			echo ===== No ${Run[$k]} image for ${subject}, cannot run FEAT ====
			fi
    done
done



# echo ==== Changing permissions ======
# chgrp -R psychology01 "$MainDir"
# chmod -R 2770 "$MainDir"


# #source deactivate /users/iv0004/psychology01/bin/conda_bash /users/iv0004/psychology01/software/conda/imaging

