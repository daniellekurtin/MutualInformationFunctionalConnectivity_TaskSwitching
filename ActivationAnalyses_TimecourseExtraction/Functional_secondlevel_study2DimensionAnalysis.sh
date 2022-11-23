#!/bin/sh

Usage() {
  cat <<EOF

usage:  Functional_higherlevel.sh -i <Subject>

e.g. sbatch Functional_higherlevel.sh  -i Subject_1

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
#SBATCH --mem=8G                      #Amount of memory per node #8G
#SBATCH --output=ROS.%N.%j.out        #Output file for stdout (optional)
#SBATCH -e ROS.%N.%j.err              #Error file

echo ${subject}
# 
#        Script for running second-level data ROS
#        ====================================
# 
#  Danielle Kurtin and Ines Violante, Sep 2021
#

#######  STEPS ############
#### 1.  Prepare tools

#source /users/psychology01/bin/conda_bash /users/psychology01/software/conda/imaging
module load fsl
module load matlab/R2019b


# ===============================
#### 2.  Define data and give permissions

export MainDir=/users/dk00549/psychology01/parallel_scratch/projects/ROS/Study_2/NewData
export ScriptsDir=/users/psychology01/projects/ROS/Scripts


# ===============================
#### 3.  Run GLM Models

export FSFdir="$ScriptsDir"/FSFfiles
export SubjectFSFdir="$FSFdir"/SubjectFSF 
Model=(DimensionAnalysis)


for ((j=0;j<${#Model[@]};++j)); do

	#if [ -f "$MainDir"/${subject}/1stlevel/${Model[$j]}/run1.feat ]   #make sure the first level FEATs are there
	#then
		if [ ! -d "$MainDir"/${subject}/2ndlevel/${Model[$j]}.feat ]  #check if Feat's been run; if not, make the directory for it, find and replace the subject and model variables in the template fsf
		then

#			mkdir -p "$MainDir"/${subject}/2ndlevel/

			echo ==== Running FEAT ${subject} ${Model[$j]} =====

			sed -e"s/ROS_DEPRES_19/${subject}/g" -e"s/DimensionAnalysis/${Model[$j]}/g" "$FSFdir"/NewStudy2_${Model[$j]}_2ndlevel.fsf > "$SubjectFSFdir"/${subject}_NewStudy2_${Model[$j]}_2ndlevel.fsf

			feat "$SubjectFSFdir"/${subject}_NewStudy2_${Model[$j]}_2ndlevel.fsf

		else
			echo ===== ${Model[$j]}.feat, already exists, skipping this step ====
		fi
	#else
     		#echo ===== No first-level FEATs for ${subject}, cannot run higher-level FEATs ====
	#fi

done


echo ==== Changing permissions ======
#chgrp -R psychology01 "$MainDir"
#chmod -R 2770 "$MainDir"


# #source deactivate /users/iv0004/psychology01/bin/conda_bash /users/iv0004/psychology01/software/conda/imaging

