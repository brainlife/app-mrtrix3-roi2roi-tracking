<h1 align="center">MRtrix3 ROI-to-ROI tracking</h1>

<p align="center">
  Tractograms and tract classifications between pairs of ROIs, including the human optic tract and optic radiation, tracked with MRtrix3 from DWI and CSD fibre orientation distributions.
</p>

<p align="center">
  <a href="https://github.com/brainlife/abcd-spec"><img src="https://img.shields.io/badge/ABCD_Spec-v1.1-green.svg" alt="ABCD Spec v1.1"></a>
  <a href="https://doi.org/10.25663/brainlife.app.407"><img src="https://img.shields.io/badge/Brainlife-brainlife.app.407-blue.svg" alt="Run on brainlife.io"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-lightgrey.svg" alt="MIT license"></a>
</p>

> 📍 **This branch, `optic-tract-v1.1`, is the brainlife.io App Track the Human Optic TRAct (THOTRA) - mrtrix3 (dwi)** ([doi:10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407)). The other Apps in this repository are listed under **Versions on brainlife.io** below.

## 🧠 At a glance

<!-- brainlife-app:begin (generated from the App manifest; edit the manifest, not this block) -->

> **Tracks the human optic tract between LGN ROIs and the optic chiasm with MRtrix3 tckgen, using iFOD2 probabilistic or SD_STREAM deterministic tracking, and outputs a tractogram and tract classification.**

Same as the iFOD2 optic tract App but with tracking_type choosing MRtrix3's IFOD2 (probabilistic) or SD_STREAM (deterministic, no backtracking) algorithm. Uses the LGN ROI(s) in roiPair and the optic chiasm ROI in oc; ROIs must be on the DWI grid and aligned with it.

| | |
|---|---|
| 📥 **Inputs** | `neuro/dwi` · `neuro/csd` · `neuro/rois` · `neuro/mask` · `neuro/mask` (optional) ×2 |
| 📤 **Outputs** | `neuro/track/tck` · `neuro/wmc` · `neuro/mask` |
| 🐳 **Container** | `brainlife/mrtrix3:3.0.0` |
| ⚙️ **Resources** | 8 CPUs · 30 GB memory · 8 h walltime |
| 🔗 **Run it** | [doi.org/10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407) |

#### 📥 Inputs

| Input | Datatype | Required | Description |
|---|---|---|---|
| `dwi` | `neuro/dwi` | yes | This is the path's to the DWI image and it's associated bvals and bvecs. |
| `csd` | `neuro/csd` | yes | This is the path to the constrained spherical deconvolution (CSD) outputs, including all lmax csd images generated and the response.txt file. |
| `rois` | `neuro/rois` | yes | This is the path to the rois directory including both the ROIs listed in the 'roiPair' config. ROIs must be on the same grid as the DWI image (or set reslice where available) and aligned to it. |
| `mask` | `neuro/mask` (5tt_masks) | yes | This is the path to the five tissue type brain image generated from mrtrix3's 5ttgen function. |
| `brainmask` | `neuro/mask` | no | This is the path to the brain mask of the DWI image. This is optional. If left empty, app will generate brainmask automatically. |
| `wm_mask` | `neuro/mask` (white_matter) | no | This is the path to the precomputed white matter mask. if left empty, one will be generated automatically |

#### 📤 Outputs

| Output | Datatype | Description |
|---|---|---|
| `track` | `neuro/track/tck` | This is the output tractogram containing all of the tracks generated |
| `wmc` | `neuro/wmc` | This is the classification structure that refers to the generated tractogram |
| `brainmask` | `neuro/mask` | This is the brainmask output. If input was used, this is just a copy and doesn't need to be archived. |

<details>
<summary><b>⚙️ Parameters</b> (14)</summary>

| Parameter | Type | Default | Description |
|---|---|---|---|
| `min_length` | number | `10` | The minimum allowable streamline length to be generated (number) example: 10 |
| `max_length` | number | `200` | The maximum allowable streamline length to be generated (number) example: 10 |
| `lmax` | number | `8` | The maximum spherical harmonic order (lmax) to use. If left empty, this will be identified automatically. (number) |
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

</details>

<!-- brainlife-app:end -->

## 🏷️ Versions on brainlife.io

Each App is built from its own branch of this repository.

| App | Tracks | Branch | Run it | Status |
|---|---|---|---|---|
| **Track the Human Optic TRAct (THOTRA) - mrtrix3 (dwi)** (this branch) | LGN to optic chiasm, iFOD2 or SD_STREAM | [`optic-tract-v1.1`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/optic-tract-v1.1) | [doi:10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407) | current |
| Track the Human Optic TRAct (THOTRA) - mrtrix3 iFOD2 (dwi) | LGN to optic chiasm, iFOD2 | [`optic-tract-v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/optic-tract-v1.0) | [doi:10.25663/brainlife.app.352](https://doi.org/10.25663/brainlife.app.352) | available |
| Track the Human Optic RAdiation (THORA) - mrtrix3 iFOD2 (dwi) | LGN to V1, iFOD2 | [`thora-v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/thora-v1.0) | [doi:10.25663/brainlife.app.347](https://doi.org/10.25663/brainlife.app.347) | available |
| Mrtrix3 iFOD2 ROI Tracking (dwi) | any seed/termination ROI pairs, iFOD2 | [`v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/v1.0) | [doi:10.25663/brainlife.app.337](https://doi.org/10.25663/brainlife.app.337) | available |
| Track the Human Optic TRAct (THOTRA) - mrtrix3 IFOD2 (dwi) | LGN to optic chiasm, iFOD2 | [`optic-tract-v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/optic-tract-v1.0) | [doi:10.25663/brainlife.app.406](https://doi.org/10.25663/brainlife.app.406) | deprecated |
| Track the Human Optic RAdiation - mrtrix3 iFOD2 (dtiinit) | LGN to V1, dtiInit input | [`thora-v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/thora-v1.0) | [doi:10.25663/brainlife.app.348](https://doi.org/10.25663/brainlife.app.348) | deprecated |
| Mrtrix3 iFOD2 ROI Tracking (dtiinit) | ROI pairs, dtiInit input | [`v1.0`](https://github.com/brainlife/app-mrtrix3-roi2roi-tracking/tree/v1.0) | [doi:10.25663/brainlife.app.339](https://doi.org/10.25663/brainlife.app.339) | deprecated |

## 🚀 Running the App

### ☁️ On brainlife.io

Open the App at [doi.org/10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407) (or another App from the table) and choose, from your project, a `neuro/dwi` dataset, the matching `neuro/csd` dataset (CSD images for each lmax and the response), a `neuro/rois` dataset and a 5-tissue-type `neuro/mask`. A brain mask and a white-matter mask are optional; if left empty they are generated. Then set the ROIs and press **Execute**.

- **ROIs.** The values are the ROI names as they appear in the `rois` dataset (after the `ROI` prefix and before `.nii.gz`). For the optic tract Apps, `roiPair` holds the LGN ROI(s), one per hemisphere you want to track, and `oc` the optic chiasm ROI (for example `085`); for the optic radiation App, `v1` names the V1 ROI. For the generic ROI tracking App, enter each seed/termination pair separated by a space (for example `0001 0002`), one pair per line; the classification will contain one tract per pair. `exclusion` lists exclusion ROIs in the same format.
- **Algorithm.** `tracking_type` selects MRtrix3 `iFOD2` (probabilistic) or `SD_STREAM` (deterministic, no backtracking). With `act` on, tracking uses the Anatomically Constrained Tractography framework with backtracking; with it off, the white-matter mask is used as the tracking mask.
- **Ensembles.** `curvatures`, `stepsize` and `minfodamp` (FOD cutoff) accept several space-separated values; every combination is tracked and the results merged. Set `single_lmax` to false to also ensemble over lmax; leave `lmax` empty to compute it from the number of diffusion directions.
- **Other options.** `count` (streamlines per tract, default 500), `min_length`, `max_length`, `maxtrials` and `multiple_seed` (seed in both ROIs and track bidirectionally).

ROIs must be on the DWI grid and aligned with it. ROIs from the *Generate ROIs in dMRI space* App already are; for the generic ROI tracking App, set `reslice` when the ROIs come from elsewhere. The outputs are a `neuro/track/tck` tractogram, a `neuro/wmc` classification and the brain mask.

### 💻 Locally

You need [Singularity/Apptainer](https://apptainer.org/) and [`jq`](https://jqlang.github.io/jq/). MRtrix3 3.0.0 and the compiled MATLAB classification step (MATLAB Runtime R2017a) run inside containers.

1. Clone this repository and check out the App you want (for example `git checkout optic-tract-v1.1`).
2. Create `config.json` with the paths to your inputs and the options you want:

   ```json
   {
       "dwi": "/input/dwi/dwi.nii.gz",
       "bvecs": "/input/dwi/dwi.bvecs",
       "bvals": "/input/dwi/dwi.bvals",
       "lmax2": "/input/csd/lmax2.nii.gz",
       "lmax4": "/input/csd/lmax4.nii.gz",
       "lmax6": "/input/csd/lmax6.nii.gz",
       "lmax8": "/input/csd/lmax8.nii.gz",
       "response": "/input/csd/response.txt",
       "rois": "/input/rois/rois",
       "mask": "/input/5tt/mask.nii.gz",
       "roiPair": "<LGN ROI names>",
       "oc": "<optic chiasm ROI name>",
       "exclusion": "",
       "tracking_type": "ifod2",
       "act": true,
       "lmax": 8,
       "single_lmax": true,
       "curvatures": "15",
       "stepsize": "0.25",
       "minfodamp": "0.025",
       "count": 500,
       "min_length": 10,
       "max_length": 200,
       "maxtrials": 1000000,
       "multiple_seed": false
   }
   ```

3. Run the App:

   ```bash
   ./main
   ```

## 📦 Outputs

```text
track/
└── track.tck               all tracts merged into one tractogram
wmc/
├── classification.mat      assignment of streamlines to tracts
└── tracts/                 per-tract JSON files for visualisation
brainmask/                  brain mask (a copy of the input if one was given)
raw/                        per-run tractograms, track_info files and intermediate images
```

## 🩺 Troubleshooting

- **`tracking did not generate. please check derivatives and log files for debugging`**: no streamlines were produced. Check that the ROI names match files in the `rois` dataset and that the ROIs are on the DWI grid; try turning `act` off or relaxing the curvature.

## 📚 Citing

> [!TIP]
> If you use this App in your work, please cite the App, the software it builds on and brainlife.io (below). GitHub's **Cite this repository** button gives the App citation in APA and BibTeX.

- Tournier J-D, Calamante F, Connelly A. Improved probabilistic streamlines tractography by 2nd order integration over fibre orientation distributions. *Proc Intl Soc Mag Reson Med* 18:1670 (2010).
- Track the Human Optic TRAct (THOTRA) - mrtrix3 (dwi). brainlife.io App. [doi:10.25663/brainlife.app.407](https://doi.org/10.25663/brainlife.app.407)
- Fischl. FreeSurfer. *NeuroImage* 62:774-781 (2012). [doi:10.1016/j.neuroimage.2012.01.021](https://doi.org/10.1016/j.neuroimage.2012.01.021)
- Tournier et al. MRtrix3: A fast, flexible and open software framework for medical image processing and visualisation. *NeuroImage* 202:116137 (2019). [doi:10.1016/j.neuroimage.2019.116137](https://doi.org/10.1016/j.neuroimage.2019.116137)
- Hayashi S, Caron BA, Heinsfeld AS, et al. brainlife.io: a decentralized and open-source cloud platform to support neuroscience research. *Nature Methods* 21:809–813 (2024). [doi:10.1038/s41592-024-02237-2](https://doi.org/10.1038/s41592-024-02237-2)
- Avesani P, McPherson B, Hayashi S, et al. The open diffusion data derivatives, brain data upcycling via integrated publishing of derivatives and reproducible open cloud services. *Scientific Data* 6:69 (2019). [doi:10.1038/s41597-019-0073-y](https://doi.org/10.1038/s41597-019-0073-y)

## 👥 Authors

- Brad Caron

**Contributors:** Soichi Hayashi, Brent McPherson

## 💰 Funding

[![Wellcome-226486/Z/22/Z](https://img.shields.io/badge/Wellcome-226486%2FZ%2F22%2FZ-purple.svg)](https://wellcome.org/)
[![NIH-NINDS-U24NS140384](https://img.shields.io/badge/NIH_NINDS-U24NS140384-green.svg)](https://reporter.nih.gov/project-details/11261037)
[![NIH-NINDS-UM1NS132207](https://img.shields.io/badge/NIH_NINDS-UM1NS132207-green.svg)](https://reporter.nih.gov/project-details/11411328)
[![NIH-NIBIB-R01EB029272](https://img.shields.io/badge/NIH_NIBIB-R01EB029272-green.svg)](https://reporter.nih.gov/project-details/10262925)
[![NIH-NIBIB-R01EB030896](https://img.shields.io/badge/NIH_NIBIB-R01EB030896-green.svg)](https://reporter.nih.gov/project-details/10428625)
[![NIH-NIMH-R01MH126699](https://img.shields.io/badge/NIH_NIMH-R01MH126699-green.svg)](https://reporter.nih.gov/project-details/10460628)
[![NIH-NIDCD-R01DC020088](https://img.shields.io/badge/NIH_NIDCD-R01DC020088-green.svg)](https://reporter.nih.gov/project-details/11140284)
[![NSF-BCS-1734853](https://img.shields.io/badge/NSF_BCS-1734853-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1734853)
[![NSF-BCS-2203524](https://img.shields.io/badge/NSF_BCS-2203524-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2203524)
[![NSF-IIS-1636893](https://img.shields.io/badge/NSF_IIS-1636893-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1636893)
[![NSF-OAC-2148729](https://img.shields.io/badge/NSF_OAC-2148729-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148729)
[![NSF-IIS-1912270](https://img.shields.io/badge/NSF_IIS-1912270-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1912270)
[![NSF-IIS-2148700](https://img.shields.io/badge/NSF_IIS-2148700-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2148700)
[![NSF-OAC-1916518](https://img.shields.io/badge/NSF_OAC-1916518-blue.svg)](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1916518)

<sub>brainlife.io is publicly funded, and acknowledging that support helps keep it running. This App was developed with support from Wellcome, the U.S. National Institutes of Health (including the BRAIN Initiative BRAIN CONNECTS program: The Axonal Projectome EXchange, U24NS140384, and the Center for Mesoscale Connectomics, UM1NS132207) and the U.S. National Science Foundation.</sub>

## ⚖️ License

[MIT](LICENSE). MRtrix3 itself is distributed under its own [license](https://www.mrtrix.org/).
