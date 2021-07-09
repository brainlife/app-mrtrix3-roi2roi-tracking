#!/bin/bash

set -x
set -e

NCORE=8

# make top dirs
mkdir -p track
mkdir -p brainmask
mkdir -p wmc
mkdir -p 5tt
mkdir -p raw

# set variables
anat=`jq -r '.anat' config.json`
dtiinit=`jq -r '.dtiinit' config.json`
dwi=`jq -r '.dwi' config.json`
bvecs=`jq -r '.bvecs' config.json`
bvals=`jq -r '.bvals' config.json`
mask=`jq -r '.mask' config.json`
brainmask=`jq -r '.brainmask' config.json`
white_matter_mask=`jq -r '.wm_mask' config.json`
max_lmax=`jq -r '.lmax' config.json`
rois=`jq -r '.rois' config.json`
count=`jq -r '.count' config.json`
roipair=`jq -r '.roiPair' config.json`
min_fod_amp=`jq -r '.minfodamp' config.json`
curvatures=`jq -r '.curvatures' config.json`
seed_max_trials=`jq -r '.maxtrials' config.json`
lmax2=`jq -r '.lmax2' config.json`
lmax4=`jq -r '.lmax4' config.json`
lmax6=`jq -r '.lmax6' config.json`
lmax8=`jq -r '.lmax8' config.json`
lmax10=`jq -r '.lmax10' config.json`
lmax12=`jq -r '.lmax12' config.json`
lmax14=`jq -r '.lmax14' config.json`
response=`jq -r '.response' config.json`
single_lmax=`jq -r '.single_lmax' config.json`
multiple_seed=`jq -r '.multiple_seed' config.json`
reslice=`jq -r '.reslice' config.json`
step_size=`jq -r '.stepsize' config.json`
min_length=`jq -r '.min_length' config.json`
max_length=`jq -r '.max_length' config.json`
act=`jq -r '.act' config.json`
exclusion=`jq -r '.exclusion' config.json`

if [ "$reslice" == "true" ]; then
    echo "using resliced_rois"
    rois=resliced_rois
fi

# parse whether dtiinit or dwi input
if [[ ! ${dtiinit} == "null" ]]; then
        input_nii_gz=$dtiinit/*dwi_aligned*.nii.gz
        BVALS=$dtiinit/*.bvals
        BVECS=$dtiinit/*.bvecs
        brainmask=$dtiinit/dti/bin/brainMask.nii.gz
        [ ! -f mask.mif ] && mrconvert ${brainmask} mask.mif -force -nthreads $NCORE
else
	input_nii_gz=${dwi}
fi

# convert input diffusion nifti to mrtrix format
[ ! -f dwi.b ] && mrconvert -fslgrad $bvecs $bvals ${input_nii_gz} dwi.mif --export_grad_mrtrix dwi.b -nthreads $NCORE

# create mask of dwi
if [[ ${brainmask} == 'null' ]]; then
	[ ! -f mask.mif ] && dwi2mask dwi.mif mask.mif -nthreads $NCORE
else
	echo "brainmask input exists. converting to mrtrix format"
	[ ! -f mask.mif ] && mrconvert ${brainmask} -stride 1,2,3,4 mask.mif -force -nthreads $NCORE
fi

# generate 5-tissue-type (5TT) tracking mask
if [[ ${mask} == 'null' ]]; then
	[ ! -f anat.mif ] && mrconvert ${anat} anat.mif -force -nthreads $NCORE
	[ ! -f 5tt.mif ] && 5ttgen fsl anat.mif 5tt.mif -nocrop -sgm_amyg_hipp -scratch ./tmp -force -nthreads $NCORE
else
	echo "input 5tt mask exists. converting to mrtrix format"
	[ ! -f 5tt.mif ] && mrconvert ${mask} -stride 1,2,3,4 5tt.mif -force -nthreads $NCORE
fi

## generate csf,gm,wm masks
[ ! -f gm.mif ] && mrconvert -coord 3 0 5tt.mif gm.mif -force -nthreads $NCORE
[ ! -f csf.mif ] && mrconvert -coord 3 3 5tt.mif csf.mif -force -nthreads $NCORE
[ ! -f csf_bin.nii.gz ] && mrconvert csf.mif -stride 1,2,3,4 csf.nii.gz -force -nthreads $NCORE && fslmaths csf.nii.gz -thr 0.5 -bin csf_bin.nii.gz

if [[ ${white_matter_mask} == "null" ]]; then

	[ ! -f wm.mif ] && mrconvert -coord 3 2 5tt.mif wm.mif -force -nthreads $NCORE
	[ ! -f wm.nii.gz ] && mrconvert wm.mif -stride 1,2,3,4 wm.nii.gz -force -nthreads $NCORE
else
	cp ${white_matter_mask} wm.nii.gz
fi

[ ! -f ./5tt/mask.nii.gz ] && mrconvert 5tt.mif -stride 1,2,3,4 ./5tt/mask.nii.gz -force -nthreads $NCORE

# brainmask
[ ! -f ./brainmask/mask.nii.gz ] && mrconvert mask.mif -stride 1,2,3,4 ./brainmask/mask.nii.gz -force -nthreads $NCORE

# generate sequence of lmax spherical harmonic order for single or ensemble
if [[ ${single_lmax} == true ]]; then
	lmaxs=$(seq ${max_lmax} ${max_lmax})
else
	lmaxs=$(seq 2 2 ${max_lmax})
fi

# Run mrtrix3
if [[ ${act} == true ]]; then
	act_line="-act 5tt.mif"
	backtrack_line="-backtrack"
else
	act_line="-mask total_mask.nii.gz"
	backtrack_line=""
fi

pairs=($roipair)
range=` expr ${#pairs[@]}`
nTracts=` expr ${range} / 2`

if [[ ! ${exclusion} == "" ]]; then
	exclus=($exclusion)
	exclude='true'
else
	exclus=""
	exclude='false'
fi

for (( i=0; i<$nTracts; i+=1 )); do
	[ -f track$((i+1)).tck ] && continue

	echo "creating seed for tract $((i+1))"
	if [ ! -f $rois/ROI${pairs[$((i*2))]}.nii.gz ]; then
		roi1=$rois/${pairs[$((i*2))]}.nii.gz
		if [[ ${exclude} == true ]]; then
			Exclusion=$rois/${exclus[$((i-1))]}.nii.gz
			exclusion_line="-exclude ${Exclusion}"
		else
			exclusion_line=""
		fi

	else
		roi1=$rois/ROI${pairs[$((i*2))]}.nii.gz
		if [[ ${exclude} == true ]]; then
			Exclusion=$rois/ROI${exclus[$((i-1))]}.nii.gz
			exclusion_line="-exclude ${Exclusion}"
		else
			exclusion_line=""
		fi
	fi

	if [ ! -f $rois/ROI${pairs[$((i*2+1))]}.nii.gz ]; then
		roi2=$rois/${pairs[$((i*2+1))]}.nii.gz
	else
		roi2=$rois/ROI${pairs[$((i*2+1))]}.nii.gz
	fi

	if [[ ${multiple_seed} == true ]]; then
		seed=seed_${pairs[$((i*2))]}_${pairs[$((i*2+1))]}.nii.gz
		[ ! -f $seed ] && mrcalc $roi1 $roi2 -add $seed -force -quiet -nthreads $NCORE && fslmaths $seed -bin $seed
		l1="-include ${roi1}"
		l2="-include ${roi2}"
		l3=""
		if [[ ${act} == false ]];  then
			[ ! -f total_mask.nii.gz ] && mrtransform wm.nii.gz wm_dwi.nii.gz -template dwi.mif -interp nearest -force -nthreads $NCORE -quiet &&  mrcalc $seed wm_dwi.nii.gz -add total_mask.nii.gz -force -quiet -nthreads $NCORE && fslmaths total_mask.nii.gz -bin total_mask.nii.gz
		fi
	else
		seed=$roi1
		l1="-include ${roi2}"
		l2=""
		l3="-seed_unidirectional"
		if [[ ${act} == false ]]; then
			[ ! -f total_mask.nii.gz ] && mrtransform wm.nii.gz wm_dwi.nii.gz -template dwi.mif -interp nearest -force -nthreads $NCORE -quiet && mrcalc $roi1 $roi2 -add temp_mask.nii.gz -force -quiet -nthreads $NCORE && mrcalc temp_mask.nii.gz wm_dwi.nii.gz -add total_mask.nii.gz -force -quiet -nthreads $NCORE && fslmaths total_mask.nii.gz -bin total_mask.nii.gz
		fi
	fi


	for LMAXS in ${lmaxs}; do
		input_csd=$(eval "echo \$lmax${LMAXS}")
		mrconvert ${input_csd} csd${LMAXS}.mif -force -nthreads $NCORE -quiet
		echo "running tracking with mrtrix3 iFOD2 tracking on lmax ${LMAXS}"
		for CURV in ${curvatures}; do
			echo "curvature ${CURV}"
			for STEP in ${step_size}; do
				echo "step size ${STEP}"
				for FOD in ${min_fod_amp}; do
					echo "FOD amplitude ${FOD}"
					if [ ! -f track$((i+1))_lmax${LMAXS}_curv${CURV}_step${STEP}_amp${FOD}.vtk ]; then
							output="track$((i+1))_lmax${LMAXS}_curv${CURV}_step${STEP}_amp${FOD}.tck"
							tckgen ${input_csd} \
								-algorithm iFOD2 \
								${act_line} \
								${backtrack_line} \
								-select ${count} \
								-seed_image ${seed} \
								${l1} \
								${l2} \
								-minlength ${min_length} \
								-maxlength ${max_length} \
								-step ${STEP} \
								-angle ${CURV} \
								-cutoff ${FOD} \
								-trials ${seed_max_trials} \
								${l3} \
								${exclusion_line} \
								$output \
								-quiet \
								-force \
								-nthreads $NCORE 
							fi
					done
				done
			done
		done
		
		output=track$((i+1)).tck
		tcks=(track$((i+1))*.tck)
		if [ ${#tcks[@]} == 1 ]; then
			mv ${tcks[0]} $output
		else
			tckedit ${tcks[*]} $output
			mv ${tcks[*]} ./raw/
		fi
		tckinfo $output > track_info$((i+1)).txt
	done

if [ -f track1.tck ]; then
	mv *.mif *.b* *.nii.gz ./raw/
	holder=(track*.tck)
	if [ ${#holder[@]} == 1 ]; then
		cp -v ${holder[0]} ./track/track.tck
	else
		tckedit ${holder[*]} ./track/track.tck
	fi
	tckinfo ./track/track.tck > ./track/track_info.txt
else
	echo "tracking did not generate. please check derivatives and log files for debugging"
	exit 1
fi

