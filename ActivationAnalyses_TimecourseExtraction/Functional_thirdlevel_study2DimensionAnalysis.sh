#!/bin/sh

Usage() {
  cat <<EOF

usage:  Functional_thirdlevel_v3.sh -i <cope>

e.g. sbatch Functional_thirdlevel_v3.sh  -i cope1

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
                    cope=$2;
                    shift;;

  esac
  shift
done


echo " Input job = $cope "
# ========================================

#SBATCH --partition=shared            #Selecting “shared” Queue
#SBATCH --job-name="ROS"              #Name of Jobs (displayed in squeue)
#SBATCH --nodes=2                     #No of nodes to run job #2
#SBATCH --ntasks-per-node=8           #No of cores to use per node #8
#SBATCH --time=01-00:00:00            #Maximum time for job to run
#SBATCH --mem=8G                      #Amount of memory per node #8G
#SBATCH --output=ROS.%N.%j.out        #Output file for stdout (optional)
#SBATCH -e ROS.%N.%j.err              #Error file


# 
#  Script for running third-level data ROS
#  ====================================
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


export ScriptsDir=/users/psychology01/projects/ROS/Scripts


# ===============================
#### 3.  Run GLM Models

export FSFdir="$ScriptsDir"/FSFfiles
export GROUPFSFDir="$FSFdir"/GroupLevelFSF 
#Model=(Model1 Model2)
Model=(DimensionAnalysis)


#if [ -f "$MainDir"/${subject}/1stlevel/${Model[$j]}/run1.feat ]   #make sure the first level FEATs are there
#then
#if [ ! -d "$MainDir"/GroupLevel/${Model[$j]}/${Model[$j]}_3rdlevel_${Cope[c]}.feat] 
#then

echo ==== Running FEAT ${Model} ${cope} =====

sed -e"s/cope1/${cope}/g" "$FSFdir"/NewStudy2_${Model}_3rdlevel.fsf > "$GROUPFSFDir"/NewStudy2_${Model}_3rdlevel_${cope}.fsf

feat "$GROUPFSFDir"/NewStudy2_${Model}_3rdlevel_${cope}.fsf



