#!/bin/bash
#extension="_4mm_2.5-3bl"
extension="_4mm_3-3.5bl"
extension="_4mm"
if [[ ( $@ == "--help") ||  $@ == "-h" || $# == 0 ]]
then
  	echo " "
	echo "Usage: $0 <projectDirectory> <subjectNumber> <norm>"
	echo " "
	echo "projectDirectory...BIDS directory containing /sub-xxx/ses-yyy/anat/sub-xxx_anat.nii"
        echo "subjectNumber......value of xxx above. For example 001,002 etc."
	echo "norm...............set 1 if using Tstats with _norm suffix"
        echo " "
        echo "Needs: projectDirectory/derivatives/sourcespace/sub-xxx/sub-xxx_brain.nii (run fsl brain extraction tool or use fieldtrip)"
        echo "       projectDirectory/derivatives/sourcespace/sub-xxx/sub-xxx_brain_4mm.nii (run fsl brain extraction tool or use fieldtrip)"
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
ANAT_4mm="${project_dir}/derivatives/sourcespace/sub-${sub}/sub-${sub}_brain_4mm_mm"
ANAT2MNI_4mm="${project_dir}/derivatives/sourcespace/sub-${sub}/sub-${sub}_anat2mni_4mm"
MNI_brain="${project_dir}/derivatives/sourcespace/MNI152_T1_4mm_brain.nii.gz"

ses="001"
task="faces_circles"
run="001"

tstat="${project_dir}/derivatives/Tstats/sub-${sub}/sub-${sub}_ses-${ses}_task-${task}_run-${run}_pseudoT_circles_adj_mm"


# Check files exist

if ! [[ -f "$ANAT_4mm.nii" || -f "$ANAT_4mm.nii.gz" ]]; then
	   echo "$ANAT_4mm doesn't exist!"
	      exit 1
fi

if ! [ -f "$MNI_brain" ]; then
	   echo "$MNI_brain doesn't exist!"
	      exit 1
fi
if ! [ -f "$tstat.nii" ]; then
	echo "$tstat doesn't exist!"
	exit 1
fi

flirt -in ${ANAT_4mm}.nii -ref ${MNI_brain} -omat ${ANAT2MNI_4mm}.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear
# final xfm from 4mm to 4mm to MNI space
echo applying flirt xfm to vlumes
fslchpixdim ${tstat}.nii 4 4 4 

flirt -in ${tstat}.nii -applyxfm -init ${ANAT2MNI_4mm}.mat -out ${tstat}_MNI_4mm -paddingsize 0.0 -interp trilinear -ref ${MNI_brain}
# Run all subs like this:
# echo {001..027} | xargs -n 1 | xargs -I {} ~/braillekids/desync_peaks2mni_4mm.sh ~/bids_braillekids {} 0
