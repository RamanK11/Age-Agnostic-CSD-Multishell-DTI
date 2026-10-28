# MRtrix3 diffusion and ACT workflow

A single age-selectable workflow runs AP/PA DWI preprocessing, TOPUP/EDDY, Dhollander three-tissue response estimation, MSMT-CSD, five-tissue segmentation, and ACT tractography. Bash and Python runners are provided. They are research tools; review outputs and QC before interpreting results.

## Choose the age-specific route

Use the age flag in the top-level configuration block at the top of the Bash or Python script. Enter age at scan in months.

| Age at scan | Bash settings | Python settings | Route |
|---|---|---|---|
| Newborn (0 months) | `age_less_than_60m=true`, `age_months=0` | `AGE_LESS_THAN_60M = True`, `AGE_MONTHS = 0` | Infant FreeSurfer newborn mode (`--newborn`) |
| Infant / paediatric, under 60 months | `age_less_than_60m=true`, set `age_months` | `AGE_LESS_THAN_60M = True`, set `AGE_MONTHS` | Age-specific Infant FreeSurfer (`--age`) |
| 60 months or older | `age_less_than_60m=false`; age may be blank or >=60 | `AGE_LESS_THAN_60M = False`; `AGE_MONTHS = None` or >=60 | FSL `5ttgen` adult OOB/BATMAN-style route |

For newborns, use age `0`: the script selects Infant FreeSurfer's dedicated newborn mode. Ages greater than `0` and less than `60` use Infant FreeSurfer's age-specific mode. At `60` months and above, select the FSL route. The Infant FreeSurfer branch requires its SIF container and a FreeSurfer license; the FSL branch requires BET, FAST, and FIRST. Both branches share the same DWI processing stages.

## Install the environment

On an internet-connected Ubuntu x86_64 machine with Conda installed, run the setup script from the repository root (the command assumes it is stored in `conda_installation/`):

```bash
bash conda_installation/setup_neonatal_mrtrix3_conda.sh install
```

By default, the installer creates `$HOME/conda_envs/mrtrix3_v1`, includes Apptainer, and is configured to pack an archive under `$HOME/conda_env_archives/`. It refuses to overwrite an existing environment. Set `ENV_PREFIX` if you want a different path; the path printed by Conda is the actual environment location.

For direct use on the machine where Conda is installed, activate and load FSL with:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/conda_envs/mrtrix3_v1"
export FSLDIR="$CONDA_PREFIX"
export FSLOUTPUTTYPE=NIFTI_GZ
source "$FSLDIR/etc/fslconf/fsl.sh"
```

## Infant FreeSurfer files

The setup script does not download the Infant FreeSurfer image or license. Keep both outside the repository. Download the SIF and obtain a FreeSurfer license:

```bash
mkdir -p "$HOME/containers"
wget -c -O "$HOME/containers/freesurfer_infant-8.1.0.sif" \
  https://ftp.nmr.mgh.harvard.edu/pub/dist/lcnpublic/dist/infantfs/freesurfer_infant-8.1.0.sif
```

Set `INFANT_SIF` and `FS_LICENSE_FILE` near the top of the selected script to those paths. The license is available from [FreeSurfer registration](https://surfer.nmr.mgh.harvard.edu/registration.html).

## Inputs

Set `PROJECT`, `SUBJECT`, and `SESSION` at the top of the script. Inputs go under `PROJECT/in/` in this layout:

```text
in/sub-XXXX/ses-YYYY/anat/sub-XXXX_ses-YYYY_T1w.nii.gz
in/sub-XXXX/ses-YYYY/dwi/sub-XXXX_ses-YYYY_dir-AP_dwi.nii.gz + .bval + .bvec + .json
in/sub-XXXX/ses-YYYY/dwi/sub-XXXX_ses-YYYY_dir-PA_dwi.nii.gz + .bval + .bvec + .json
```

The JSON sidecars need valid signed `PhaseEncodingDirection` and `TotalReadoutTime` values. The PA series must contain a genuine b=0 volume. Check AP/PA phase encoding, image grids, readout times, and reference contrast before running TOPUP.

## Top-level configuration before running

Review the configuration block at the top of the selected script (Section 0 in Bash) and set:

- `PROJECT`, `SUBJECT`, and `SESSION` to locate the participant's input files.
- `OUT` to a fresh output directory for this subject and configuration.
- The age flag and age as described above.
- `INFANT_SIF` and `FS_LICENSE_FILE` when using Infant FreeSurfer.
- `THREADS` and `interactive` / `INTERACTIVE` for your machine and whether MRView should open.

The script builds AP/PA DWI and T1 paths from these values using the layout below. If your BIDS filenames differ, edit the corresponding path settings at the top of the script.

## Run

After completing the top-level configuration above and activating the environment, run one version:

```bash
bash -n mrtrix3_full_workflow_age_agnostic.sh
bash mrtrix3_full_workflow_age_agnostic.sh
```

Or use Python:

```bash
python3 -m py_compile mrtrix3_full_workflow_age_agnostic.py
python3 mrtrix3_full_workflow_age_agnostic.py
```

Set `interactive` in Bash or `INTERACTIVE` in Python to `true` / `True` to open MRView QC windows, or false for headless processing. Use a separate output directory per participant/configuration. A successful run is not a substitute for visual QC of EDDY, segmentation, registration, FODs, and the ACT pilot tractogram.

Do not commit participant data, outputs, the `.sif` image, or the FreeSurfer license to Git.

## Appendix: Offline deployment

After the installer completes successfully and creates the archive, transfer the archive and setup script to the offline target. Run:

```bash
bash setup_neonatal_mrtrix3_conda.sh unpack /path/to/mrtrix3_v1.tar.gz
bash setup_neonatal_mrtrix3_conda.sh check
```

## Appendix: Top-level configuration reference

This section explains the settings in the Bash configuration block. The Python runner uses the corresponding uppercase constants where applicable. Values such as `WORK`, `QC`, `COMMON`, and `TRACK_OPTIONS` are derived by the scripts; normally edit their source settings rather than the derived values.

### Participant and file paths

| Bash setting | Python setting | Purpose |
|---|---|---|
| `PROJECT` | `PROJECT` | Project root. The `in/` and output folders are located under this directory. |
| `SUBJECT` | `SUBJECT` | BIDS subject label, for example `sub-XXXX`. |
| `SESSION` | `SESSION` | BIDS session label, for example `ses-YYYY`. |
| `INPUT` | `INPUT` | Derived subject/session input directory: `PROJECT/in/SUBJECT/SESSION`. |
| `DWI`, `BVAL`, `BVEC`, `JSON_AP` | Same names | AP diffusion NIfTI, its gradient table, and its JSON metadata. |
| `REVERSE`, `REVERSE_BVAL`, `REVERSE_BVEC`, `JSON_PA` | Same names | Reverse-PE (PA) series and its matching gradient and JSON files. PA is used to form the TOPUP pair. |
| `T1` | `T1` | Structural T1-weighted image used for five-tissue segmentation. |
| `OUT` | `OUT` | Output directory for this participant/session and configuration. Use a new directory for a new run. |

The current Bash and Python scripts construct exact filenames from `SUBJECT` and `SESSION`. A `*` written inside a quoted path is treated literally; these scripts do not expand wildcards. If your filenames differ, edit the path construction so each variable points to one real file.

### Route and execution

| Bash setting | Python setting | Purpose |
|---|---|---|
| `age_less_than_60m` | `AGE_LESS_THAN_60M` | `true` / `True` selects Infant FreeSurfer; `false` / `False` selects FSL `5ttgen`. |
| `age_months` | `AGE_MONTHS` | Age at scan in months. Required below 60 months; use `0` for newborn mode. Leave blank / `None` for the FSL route, or enter an age of at least 60. |
| `interactive` | `INTERACTIVE` | Enables MRView checkpoints when true. False skips viewers for headless runs. `show_mrview` is the Bash helper that applies this setting. |
| `INFANT_SIF` | `INFANT_SIF` | Path to the Infant FreeSurfer `.sif` image; used only in the InfantFS branch. |
| `FS_LICENSE_FILE` | `FS_LICENSE_FILE` | Path to the FreeSurfer license; used only in the InfantFS branch and mounted read-only into the container. |
| `THREADS` | `THREADS` | Thread count for MRtrix3 and the ITK/ANTs settings. The scripts set `OMP_NUM_THREADS` and `ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS` from it. |

`FSLOUTPUTTYPE=NIFTI_GZ` selects compressed NIfTI for FSL outputs. `FSLDIR`, Conda activation, and `PATH` examples in the Bash file are optional environment setup instructions; use the values for your actual installation. Do not set `FSLDIR` to `CONDA_PREFIX` unless that prefix really contains the FSL installation. `CUDA_VISIBLE_DEVICES` is an optional GPU selection and is not needed for the default CPU workflow.

### Diffusion and tractography parameters

| Bash setting | Python setting | Purpose |
|---|---|---|
| `BZERO_THRESHOLD` | `BZERO_THRESHOLD` | MRtrix b=0 classification threshold. The current value, `10`, excludes the reverse-series b=50 volume from the TOPUP b=0 pair. |
| `GIBBS_AXES` | `GIBBS_AXES` | Image-axis pair passed to `mrdegibbs`. Confirm these axes match the acquisition plane. |
| `EDDY_OPTIONS` | `EDDY_OPTIONS` | Additional EDDY options; current defaults enable outlier replacement and shelled-data handling. Slice-to-volume correction is not enabled. |
| `DHOLLANDER_FA` | `DHOLLANDER_FA` | FA threshold passed to Dhollander response estimation. |
| `DTI_SHELLS` | `DTI_SHELLS` | Shells used only for the optional tensor metrics. This does not control MSMT-CSD. |
| `TRACK_PILOT` | `TRACK_PILOT` | Number of streamlines requested for the pilot ACT tractogram, which should be reviewed before the final run. |
| `TRACK_FINAL` | `TRACK_FINAL` | Number of streamlines requested for the final ACT tractogram. |
| `FOD_CUTOFF` | `FOD_CUTOFF` | FOD cutoff for tractography. MRtrix applies its ACT-specific factor to this value. |
| `MIN_LENGTH_MM`, `MAX_LENGTH_MM`, `STEP_MM` | Same names | Optional tractography length and step settings in mm. Leave blank / `None` to use MRtrix defaults. |

`TRACK_OPTIONS` is assembled from the cutoff and optional length/step settings. `COMMON` carries the MRtrix thread count and b=0 threshold. `BIAS_COMMON` passes the threshold and ITK thread configuration to bias correction; it intentionally does not pass MRtrix `-nthreads` to the N4 wrapper.

### Output and Infant FreeSurfer paths

`WORK`, `QC`, `METRICS`, `FOD`, `ACT`, `REG`, `TRACKS`, and `SCRATCH` are subdirectories derived from `OUT` for intermediate files, QC, metrics, FODs, anatomy, registration, tractography, and temporary processing data. `BIAS`, `MASK`, and `WM_FOD` name key derived outputs used by later stages.

`INFANT_OUT` is the Infant FreeSurfer output directory. `INFANT_SEG` and `INFANT_REF` identify its expected segmentation and structural reference. `CONTAINER_LUT` is the FreeSurfer color lookup table inside the SIF. `LUT_RAW` stores a copy of that table, and `LUT_MRTRIX` is the normalized copy used by MRtrix. `FIVE_TT` is selected automatically: `5tt_infantfs.mif` for Infant FreeSurfer or `5tt_fsl.mif` for the FSL route.
