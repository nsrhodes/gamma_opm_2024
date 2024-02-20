#!/bin/bash
#extension="_2.5-3bl"
#extension="_4mm_2.5-3bl"
#extension="_4mm_3-3.5bl"
extension=""
if [[ ( $@ == "--help") ||  $@ == "-h" || $# == 0 ]]
then
  	echo " "
	echo "Usage: $0 <projectDirectory> <subjectNumber> <norm>"
	echo " "
	echo "projectDirectory...BIDS directory containing /sub-xxx/ses-yyy/anat/sub-xxx_anat.nii"
        echo "subjectNumber......value of xxx above. For example 001,002 etc."
	echo "norm...............set 1 if using Tstats with _norm suffix"
        echo " "
        echo "Needs: projectDirectory/derivatives/sourcespace/sub-xxx/sub-xxx_brain.nii (run fsl brain extraction tool)"
        echo "       projectDirectory/derivatives/Tstats${extension}/sub-xxx/D1_peak.txt (peak desync location for index)"
        echo "       projectDirectory/derivatives/Tstats${extension}/sub-xxx/D4_peak.txt (peak desync location for pinky)"
        echo "       projectDirectory/derivatives/sourcespace/MNI152_T1_1mm_brain.nii.gz (MNI template)"
	echo "	     projectDirectory/derivatives/sourcespace/sub-xxx/sub-xxx_anat2mni.mat (flirt transform from anatomical to MNI space)"
	exit 0
fi

project_dir=$1 #Bids directory containing /sub-001/ses-001/anat/sub-001_anat.nii
sub=$2 #001,002....
echo $sub
ANAT="${project_dir}/derivatives/sourcespace/sub-${sub}/sub-${sub}_brain"
ANAT2MNI="${project_dir}/derivatives/sourcespace/sub-${sub}/sub-${sub}_anat2mni.mat"
MNI_brain="${project_dir}/derivatives/sourcespace/MNI152_T1_1mm_brain.nii.gz"
peak="${project_dir}/derivatives/Tstats${extension}/sub-${sub}/circles_peak_vm_SKchans.txt"
peak_mni="${project_dir}/derivatives/Tstats${extension}/sub-${sub}/circles_peak_vm__SKchans_mni.txt"


# Check files exist
if ! [[ -f "$ANAT.nii" || -f "$ANAT.nii.gz" ]]; then
	   echo "$ANAT doesn't exist!"
	      exit 1
fi
if ! [ -f "$MNI_brain" ]; then
	   echo "$MNI_brain doesn't exist!"
	      exit 1
fi
if ! [ -f "$peak" ]; then
	   echo "$peak doesn't exist!"
	      exit 1
fi
if ! [ -f "$ANAT2MNI" ]; then
	   echo "$ANAT2MNI doesn't exist!"
	      exit 1
fi

img2imgcoord -src ${ANAT} -dest ${MNI_brain} -xfm ${ANAT2MNI} -mm ${peak} > ${peak_mni}
