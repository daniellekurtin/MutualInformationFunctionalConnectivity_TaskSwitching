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

export MainDir=/users/psychology01/parallel_scratch/projects/ROS/Study_2
export ScriptsDir=/users/psychology01/parallel_scratch/projects/ROS/Scripts


# ===============================
#### 3.  Run GLM Models

export FSFdir="$ScriptsDir"/FSFfiles
export SubjectFSFdir="$FSFdir"/SubjectFSF 
Model=(Model3)


for ((j=0;j<${#Model[@]};++j)); do
	
			echo ==== Running FEAT ${subject} ${Model[$j]} =====

			sed -e"s/ROS_DEPRES_1/${subject}/g" -e"s/Model3/${Model[$j]}/g" "$FSFdir"/Study2_${Model[$j]}_2ndlevel_NoInCorr.fsf > "$SubjectFSFdir"/${subject}_Study2_${Model[$j]}_2ndlevel_NoInCorr.fsf

			feat "$SubjectFSFdir"/${subject}_Study2_${Model[$j]}_2ndlevel_NoInCorr.fsf


done


