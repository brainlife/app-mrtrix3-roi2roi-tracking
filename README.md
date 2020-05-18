[![Abcdspec-compliant](https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg)](https://github.com/brain-life/abcd-spec)
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-brainlife.app.337-blue.svg)](https://doi.org/10.25663/brainlife.app.337)

# Mrtrix3 iFOD2 ROI Tracking (dwi) 

This app will perform ensemble tracking between 2 or more cortical regions of interest (ROIs). This app uses MrTrix3.0's tckgen iFOD2 probabilistic tracking functionality to track between pairs of ROIs. If the user desires to use the Anatomically Constrained Tractography (ACT) framework, the user can set the 'act' field to true. If not desired, or if not providing adequate results, this can be set to false. This app takes in the following required inputs: DWI, rois, CSD, T1, Freesurfer. Optionally, a user can input precomputed brainmask, 5-tissue type segmentation, and white matter masks. If these inputs are not provided, the app will generate these datatypes internally and output the derivatives. This is why the Freesurfer and T1 are required inputs.

To specify which ROIs are desired, the user can specify each pair as the ROI numbers found in the ROI datatype seperated by a space. For example, if the user wanted to track between ROIs 0001 and 0002 from their ROIs datatype, the user should input 0001 0002 in the roiPair field. The first value will generally be treated as the seed ROI and the second as the termination ROI. However, the user can specify seeding in both ROIs by setting the 'multiple_seeds' field to true. If the user wants to make multiple tracks, enter the user can input multiple pairs by creating a new line in the roiPair field. The output classification structure will then contain the same number of tracks as ROI pairs.

This app provides the user with a large number of exposed parameters to guide and shape tractography. These include a maximum spherical harmonic order (lmax), number of repetitions, minimum and maximum length of streamlines, step size, maximum number of attempts, streamline count, FOD cutoff, and maximum angle of curvature. For lmax, the user can specify whether or not to track on a single lmax or 'ensemble' across lmax's. If the user wants to track in just a single lmax, set the 'single_lmax' field to true. Else, leave as false. If the user does not know which lmax to use, they can leave the 'lmax' field empty and the app will automatically calculate the maximum lmax based on the number of unique directions in the non-b0 weighted volumes of the DWI. For FOD cutoff, maximum angle of curvature, and step size, the user can input multiple values to perform 'ensemble tracking'. If this is desired, the user can input each value separated by a space in the respective fields (example: 0.25 0.5). The outputs of each iteration will be merged together in the final output.

This app requires the ROIs and DWI datatypes to have the same dimensions. If the ROIs were generated with the 'Generate ROIs in dMRI space' app, then the ROIs and DWI will have the same dimensions. If another app was used to generate the ROIs input, the user will need to set the 'reslice' field to true. This will reslice the ROIs to have the same dimensions as the DWI image. This also assumes proper alignment between the DWI image and the ROIs.
This app uses multiple software packages, including Freesurfer, MrTrix3.0, Matlab, and python. 

### Authors 

- Brad Caron (bacaron@iu.edu) 

### Contributors 

- Soichi Hayashi (hayashis@iu.edu)
- Brent McPherson (bcmcpher@iu.edu) 

### Funding 

[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-1636893](https://img.shields.io/badge/NSF_BCS-1636893-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-ACI-1916518](https://img.shields.io/badge/NSF_ACI-1916518-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1916518)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://grantome.com/grant/NIH/R01-EB029272-01)

### Citations 

Please cite the following articles when publishing papers that used data, code or other resources created by the brainlife.io community. 

1. Tournier, J.-D.; Calamante, F. & Connelly, A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. Proceedings of the International Society for Magnetic Resonance in Medicine, 2010, 1670 

## Running the App 

### On Brainlife.io 

You can submit this App online at [https://doi.org/10.25663/brainlife.app.337](https://doi.org/10.25663/brainlife.app.337) via the 'Execute' tab. 

### Running Locally (on your machine) 

1. git clone this repo 

2. Inside the cloned directory, create `config.json` with something like the following content with paths to your input files. 

```json 
{
   "dwi":    "testdata/dwi/dwi.nii.gz",
   "bvals":    "testdata/dwi/dwi.bvals",
   "bvecs":    "testdata/dwi/dwi.bvals",
   "lmax2":    "testdata/csd/lmax2.nii.gz/",
   "rois":    "testdata/rois/rois/",
   "anat":    "testdata/anat/t1.nii.gz",
   "min_length":    10,
   "max_length":    200,
   "lmax":    2,
   "minfodamp":    "0.025",
   "roiPair":    "001 0002",
   "stepsize":    "0.25",
   "count":    500,
   "curvatures":    "45",
   "maxtrials":    100000,
   "single_lmax":    true,
   "reslice":    false,
   "multiple_seed":    false,
   "act":    true
} 
``` 

### Sample Datasets 

You can download sample datasets from Brainlife using [Brainlife CLI](https://github.com/brain-life/cli). 

```
npm install -g brainlife 
bl login 
mkdir input 
bl dataset download 
``` 

3. Launch the App by executing 'main' 

```bash 
./main 
``` 

## Output 

The main output of this App is contains the whole-brain tractogram (tck) and the internally computed masks, a white-matter classification structure (wmc), and all of the other derivatives generated. If masks were inputted, the output is simply copies of the inputs. 

#### Product.json 

The secondary output of this app is `product.json`. This file allows web interfaces, DB and API calls on the results of the processing. 

### Dependencies 

This App requires the following libraries when run locally. 

- MrTrix3.0: 
- Matlab: 
- jsonlab: 
- singularity: 
- FSL: 
