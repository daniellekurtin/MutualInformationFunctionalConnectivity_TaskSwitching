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

export MainDir=/users/psychology01/projects/ROS/Data/Study_2
export RawDataDir=/users/psychology01/projects/ROS/Data/Raw/Study_2
export BehDir=/users/psychology01/projects/ROS/DMN_Data/source_data/Study_2/Study_2_Behavioural_data
export ScriptsDir=/users/psychology01/projects/ROS/Scripts


# ===============================
#### 3.  Convert to dcm


if [ ! -f "$MainDir"/${subject}/preproc/*.nii.gz ]
then 
  echo ==== Organizing folders ${subject} =======
    
  mkdir -p "$MainDir"/${subject}/preproc

  echo ==== Converting to dcm ${subject} =======

  /users/psychology01/bin/dcm2niix -o "$MainDir"/${subject}/preproc/ -z y -d y -f %n_%p "$RawDataDir"/${subject}/

else
echo ==== Preproc folder with nifti files already present for ${subject}, skipping dcm2niix =======

fi


# ===============================
#### 4. Reorient to standard for T1 MPRAGE images


if [ ! -f "$MainDir"/${subject}/preproc/T1.nii.gz ]
then 
    echo ==== Converting to standard ${subject} =======
    cd "$MainDir"/${subject}/preproc
        
    find . -type f -name "*MPRAGE_ADNI_P2.nii.gz" | xargs -n 1 -I '{}' mv {} T1_native.nii.gz

    fslreorient2std T1_native.nii.gz T1.nii.gz
else
    echo ==== T1 file already present for ${subject}, skipping reorient to standard =======
fi


# ===============================
#### 5.  BET for T1 MPRAGE images

if [ ! -f "$MainDir"/${subject}/preproc/T1_brain.nii.gz ]
then
   cd "$MainDir"/${subject}/preproc
    
   echo ==== BET ${subject} =====
    
   bet T1 T1_brain
   $FSLDIR/bin/standard_space_roi T1_brain T1_cut -roiNONE -ssref ${FSLDIR}/data/standard/MNI152_T1_2mm_brain -altinput T1.nii.gz

   $FSLDIR/bin/bet T1_cut T1_brain -f 0.1 -B   

else
    echo ==== T1_brain file already present for ${subject}, skipping BET for T1 image =======
fi



# ===============================
#### 6. Fsl Motion outliers

# # # Usage: fsl_motion_outliers -i <input 4D image> -o <output confound file> [options]


if [ ! -f "$MainDir"/${subject}/preproc/MotionOutliers/*DVARS.txt ]
then 

   cd "$MainDir"/${subject}/preproc

   find . -type f -name "*_fMRI_Switch.nii.gz" | xargs -n 1 -I '{}' cp -R {} run1.nii.gz
   find . -type f -name "*_fMRI_Switch_run1.nii.gz" | xargs -n 1 -I '{}' cp -R {} run1.nii.gz
   find . -type f -name "*_fMRI_Switch_run2.nii.gz" | xargs -n 1 -I '{}' cp -R {} run2.nii.gz
   find . -type f -name "*_fMRI_Switch_run3.nii.gz" | xargs -n 1 -I '{}' cp -R {} run3.nii.gz

   mkdir "$MainDir"/${subject}/preproc/MotionOutliers

   for file in run*; do
    

       filename_aux2="$file"
       filename="${filename_aux2%.nii.gz}"


       echo ==== runnning motion outliers ${subject} "$filename" DVARS ====

       fsl_motion_outliers -i "${filename_aux2}" -o MotionOutliers/"$filename"_ConfMatrix_DVARS --dvars --thresh=50 -s MotionOutliers/"$filename"_DVARS.txt -p MotionOutliers/"$filename"_DVARS -v > MotionOutliers/"$filename"_DVARS_output.txt

       echo ==== runnning motion outliers ${subject} "$filename" FD ====

       fsl_motion_outliers -i "${filename_aux2}" -o MotionOutliers/"$filename"_ConfMatrix_FD --fd --thresh=0.2 -s MotionOutliers/"$filename"_FD.txt -p MotionOutliers/"$filename"_FD -v > MotionOutliers/"$filename"_FD_output.txt

        
   done

else
 echo ==== MotionOutliers folder with DVARS file already present for ${subject}, skiping this step =======

fi



# # # =========================================
# # # Copy EV files (assuming easy matching) 


if [ ! -d "$MainDir"/${subject}/Behaviour/Subject_* ]
then
	mkdir "$MainDir"/${subject}/Behaviour
    echo === copying ${subject} to "$MainDir"/${subject}/Behaviour/
    cp -R "$BehDir"/${subject} "$MainDir"/${subject}/Behaviour/

    subj=$(ls "$MainDir"/${subject}/Behaviour/)
    mv  "$MainDir"/${subject}/Behaviour/${subj}/EVs/1 "$MainDir"/${subject}/Behaviour/${subj}/EVs/run1
    mv  "$MainDir"/${subject}/Behaviour/${subj}/EVs/2 "$MainDir"/${subject}/Behaviour/${subj}/EVs/run2
    mv  "$MainDir"/${subject}/Behaviour/${subj}/EVs/3 "$MainDir"/${subject}/Behaviour/${subj}/EVs/run3

else
    echo ==== Behaviour directory for ${subject} already exists ===
fi

# # # # # # =========================================
# # # # # # Run GLM Models

export FSFdir="$ScriptsDir"/FSFfiles
export SubjectFSFdir="$FSFdir"/SubjectFSF 
Run=(run1 run2 run3)
Model=(Model1 Model2)

if [ ! -d "$SubjectFSFdir" ]
then
	mkdir "$SubjectFSFdir"
fi

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
              
					sed -e "s/ROS_DEPRES_2/${subject}/g" -e"s/set fmri(npts) 378/set fmri(npts) ${nvols}/g" -e"s/run1/${Run[$k]}/g" -e"s/Model1/${Model[$j]}/g" -e"s/Subject_2/${subj}/g" "$FSFdir"/Study2_${Model[$j]}_run1.fsf > "$SubjectFSFdir"/${subject}_Study2_${Model[$j]}_${Run[$k]}.fsf

			         feat "$SubjectFSFdir"/${subject}_Study2_${Model[$j]}_${Run[$k]}.fsf
			 	else

					echo ===== ${subject}_Study2_${Model[$j]}_${Run[$k]}.feat, already exists, skipping this step ====
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

