<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

**Tracks the human optic tract between LGN ROIs and the optic chiasm with MRtrix3 tckgen, using iFOD2 probabilistic or SD_STREAM deterministic tracking, and outputs a tractogram and tract classification.**

Same as the iFOD2 optic tract App but with tracking_type choosing MRtrix3's IFOD2 (probabilistic) or SD_STREAM (deterministic, no backtracking) algorithm. Uses the LGN ROI(s) in roiPair and the optic chiasm ROI in oc; ROIs must be on the DWI grid and aligned with it.

Run it on brainlife.io: [https://doi.org/10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407)

| Input | Datatype | Required | Description |
|---|---|---|---|
| `dwi` | `neuro/dwi` | yes | This is the path's to the DWI image and it's associated bvals and bvecs. |
| `csd` | `neuro/csd` | yes | This is the path to the constrained spherical deconvolution (CSD) outputs, including all lmax csd images generated and the response.txt file. |
| `rois` | `neuro/rois` | yes | This is the path to the rois directory including both the ROIs listed in the 'roiPair' config. ROIs must be on the same grid as the DWI image (or set reslice where available) and aligned to it. |
| `mask` | `neuro/mask` (5tt_masks) | yes | This is the path to the five tissue type brain image generated from mrtrix3's 5ttgen function. |
| `brainmask` | `neuro/mask` | no | This is the path to the brain mask of the DWI image. This is optional. If left empty, app will generate brainmask automatically. |
| `wm_mask` | `neuro/mask` (white_matter) | no | This is the path to the precomputed white matter mask. if left empty, one will be generated automatically |

| Output | Datatype | Description |
|---|---|---|
| `track` | `neuro/track/tck` | This is the output tractogram containing all of the tracks generated |
| `wmc` | `neuro/wmc` | This is the classification structure that refers to the generated tractogram |
| `brainmask` | `neuro/mask` | This is the brainmask output. If input was used, this is just a copy and doesn't need to be archived. |

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_length` | number | `10` | The minimum allowable streamline length to be generated (number) example: 10 |
| `max_length` | number | `200` | The maximum allowable streamline length to be generated (number) example: 10 |
| `lmax` | number | `8` | The maxmium spherical harmonic order (lmax) to use. If left empty, this will be identified automatically. (number) |
| `minfodamp` | string | `0.025` | The minimum FOD amplitude to use as a cutoff (mrtrix3 flag -cutoff) (number) ex: 0.025 |
| `roiPair` | string |  | this field is here in order to allow for the option for only tracking one hemisphere or both hemispheres. If you want to track one hemisphere, enter the exact name of the LGN roi you want to track from (example: '008109'). If both hemispheres are desired, enter both LGN ROI numbers (example: '008109 008209'). Please see your ROI datatype for the exact ROI numbers/names. |
| `stepsize` | string | `0.25` | Step size to use for tracking (number) ex: 0.25 |
| `count` | number | `500` | total number of streamlines |
| `oc` | string |  | this is the exact ROI number/name of the optic chiasm from your ROIs datatype that you would like to track to. Example: '085'. Please see your ROI datatype for the exact ROI numbers/names. |
| `maxtrials` | number | `1000000` | Sets the maximum number of attempts to generate streamline from the seed point. Default=1000000. |
| `curvatures` | string | `15` | maximum angle of curvature (in degrees). can be multiple (i.e. ensemble) or one value. default: 45 |
| `single_lmax` | boolean | `true` | Determines whether single lmax (true) or ensemble lmax (false) will be performed |
| `multiple_seed` | boolean | `false` | If true, this will set the seeds of tractography to be present in both the seeding and termination ROIs. This will also set the tracking to be bidirectional. |
| `act` | boolean | `true` | If true, this will use the anatomically constrained tractography framework for tracking. This will include backtracking. If false, the white matter mask will be used as a tracking mask. The ROIs will be added to the white matter mask. If the ROIs haven't been resliced to diffusion space, set reslice to true in order to get this to work. |
| `tracking_type` | enum |  | Tracking algorithm passed to tckgen: ifod2 (IFOD2, probabilistic) or deterministic (SD_STREAM, backtracking disabled). Options: `ifod2`, `deterministic`. |

Runtime: container `brainlife/mrtrix3:3.0.0`; 8 CPUs; 30 GB memory; 8 h walltime.

<!-- brainlife-app:end -->

# app-mrtrix3-roi2roi-tracking