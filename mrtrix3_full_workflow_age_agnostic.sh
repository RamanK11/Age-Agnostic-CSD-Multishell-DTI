#!/usr/bin/env bash
# MRTRIX3 + AGE-SELECTED INFANT/FSL 5TT + ACT — COMMAND-BY-COMMAND WORKFLOW
# age_less_than_60m=true: Infant FreeSurfer with age-specific segmentation.
# age_less_than_60m=false: standard MRtrix 5ttgen fsl route.
# Edit section 0, then copy/run numbered sections in order in the SAME Bash session.
# This is a manual runbook: inspection commands do not automatically approve QC.
# No commands are hidden inside a Python processing wrapper. Existing outputs are
# not overwritten by default. Use a new OUT for each parameter/model experiment.
# Assumes one fixed-PE AP DWI series and a reverse-PE series containing true b0.
# Uses standard FSL 5TT generation for the >=60-month branch.
# Retains the same DWI preprocessing, three-tissue Dhollander/MSMT-CSD, and ACT stages.
# If executing the entire file, inspect QC windows before closing them to continue.
set -euo pipefail

# ============================================================================
# 0. EDIT PATHS, ENVIRONMENT AND ANALYSIS SETTINGS
# ============================================================================
PROJECT="${HOME}/projects/mrtrix3"
SUBJECT="sub-XXXX"
SESSION="ses-YYYY"
INPUT="$PROJECT/in/$SUBJECT/$SESSION"
DWI="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-AP_dwi.nii.gz"
BVAL="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-AP_dwi.bval"
BVEC="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-AP_dwi.bvec"
JSON_AP="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-AP_dwi.json"
REVERSE="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-PA_dwi.nii.gz"
REVERSE_BVAL="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-PA_dwi.bval"
REVERSE_BVEC="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-PA_dwi.bvec"
JSON_PA="$INPUT/dwi/${SUBJECT}_${SESSION}_dir-PA_dwi.json"
T1="$INPUT/anat/${SUBJECT}_${SESSION}_T1w.nii.gz"
OUT="$PROJECT/out_age_switch/$SUBJECT/$SESSION"
# Set these before running. InfantFS requires scan age in months.
# Choose the branch explicitly. For InfantFS, set age_months (age at scan) below 60.
age_less_than_60m=false
age_months=""
interactive=false  # true opens MRView at QC checkpoints; false runs headless.
INFANT_SIF="${HOME}/containers/freesurfer_infant-8.1.0.sif"
FS_LICENSE_FILE="${HOME}/license.txt"
THREADS=8
if [[ "$interactive" != true && "$interactive" != false ]]; then
    echo "interactive must be exactly true or false (got: $interactive)" >&2
    exit 2
fi
show_mrview() {
    if [[ "$interactive" == true ]]; then
        mrview "$@"
    fi
}
export FSLOUTPUTTYPE=NIFTI_GZ
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS="$THREADS"
export OMP_NUM_THREADS="$THREADS"
# Load your installed environments BEFORE running commands. Examples only:
# source /path/to/conda/etc/profile.d/conda.sh
# conda activate your_mrtrix_environment
# export FSLDIR=/path/to/fsl
# source "$FSLDIR/etc/fslconf/fsl.sh"
# export PATH="/path/to/mrtrix3/bin:/path/to/ants/bin:$PATH"
# Optional GPU selection, only if appropriate for your machine:
# export CUDA_VISIBLE_DEVICES=0
# Do not infer FSLDIR from CONDA_PREFIX unless that really is your FSL installation.

BZERO_THRESHOLD=10       # Consistent MRtrix policy; excludes reverse b=50.
GIBBS_AXES="0,1"         # Verify acquisition plane in MRtrix image axes.
EDDY_OPTIONS=" --repol --data_is_shelled"
# No slice-to-volume correction is enabled by this string. If required, configure
# --mporder and correct slice timing/group information for YOUR acquisition.
DHOLLANDER_FA=0.2        # Preserves original setting; 0.1/0.15 are candidates, not mandates.
# Optional DTI branch: preserves original all-shell fit. For a separately justified
# conventional DTI analysis change this to "0,1000" or another prespecified subset.
DTI_SHELLS="0,500,1000,2000,3000"
TRACK_PILOT=200000
TRACK_FINAL=10000000
FOD_CUTOFF=0.1          # MRtrix multiplies this by 0.5 with ACT.
# Leave blank to retain voxel-size-dependent MRtrix defaults. If filled, units mm.
MIN_LENGTH_MM=""
MAX_LENGTH_MM=""
STEP_MM=""
TRACK_OPTIONS=(-cutoff "$FOD_CUTOFF")
[[ -z "$MIN_LENGTH_MM" ]] || TRACK_OPTIONS+=(-minlength "$MIN_LENGTH_MM")
[[ -z "$MAX_LENGTH_MM" ]] || TRACK_OPTIONS+=(-maxlength "$MAX_LENGTH_MM")
[[ -z "$STEP_MM" ]] || TRACK_OPTIONS+=(-step "$STEP_MM")
COMMON=(-nthreads "$THREADS" -config BZeroThreshold "$BZERO_THRESHOLD")
# Keep the N4 wrapper's command line free of -nthreads; its ITK thread count is
# controlled by ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS, set to THREADS above.
BIAS_COMMON=(-config BZeroThreshold "$BZERO_THRESHOLD" -config NumberOfThreads "$THREADS")
WORK="$OUT/work"
QC="$OUT/qc"
METRICS="$OUT/metrics"
FOD="$OUT/fod"
ACT="$OUT/act"
INFANT_OUT="$OUT/infantfs"
INFANT_SEG="$INFANT_OUT/mri/aparc+aseg.mgz"
INFANT_REF="$INFANT_OUT/mri/norm.nii.gz"
CONTAINER_LUT="/usr/local/freesurfer/luts/FreeSurferColorLUT.txt"
LUT_RAW="$OUT/FreeSurferColorLUT.txt"
LUT_MRTRIX="$ACT/FreeSurferColorLUT_mrtrix.txt"
if [[ "$age_less_than_60m" == true ]]; then
    FIVE_TT="$ACT/5tt_infantfs.mif"
elif [[ "$age_less_than_60m" == false ]]; then
    FIVE_TT="$ACT/5tt_fsl.mif"
else
    echo "age_less_than_60m must be exactly true or false (got: $age_less_than_60m)" >&2
    exit 2
fi
REG="$ACT/registration"
TRACKS="$OUT/tractography"
SCRATCH="$OUT/scratch"
BIAS="$OUT/dwi_preprocessed_biascorrected.mif"
MASK="$OUT/brain_mask.mif"
WM_FOD="$FOD/wm_fod_normalised.mif"

# ============================================================================
# 1. INPUTS, DIRECTORIES AND PROVENANCE
# ============================================================================
if [[ "$age_less_than_60m" == true ]]; then
    [[ -n "$age_months" ]] || { echo "Set age_months when age_less_than_60m=true." >&2; exit 2; }
    [[ "$age_months" =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "age_months must be a non-negative number." >&2; exit 2; }
    awk -v age="$age_months" 'BEGIN { exit !(age < 60) }' || {
        echo "age_months must be less than 60 when age_less_than_60m=true." >&2; exit 2;
    }
else
    if [[ -n "$age_months" ]]; then
        [[ "$age_months" =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo "age_months must be a non-negative number." >&2; exit 2; }
        awk -v age="$age_months" 'BEGIN { exit !(age >= 60) }' || {
            echo "Age is under 60 months; set age_less_than_60m=true and provide the age-specific InfantFS route." >&2; exit 2;
        }
    fi
fi
for f in "$DWI" "$BVAL" "$BVEC" "$JSON_AP" "$REVERSE" "$REVERSE_BVAL" \
         "$REVERSE_BVEC" "$JSON_PA" "$T1"; do
    [[ -f "$f" ]] || { echo "Missing input: $f" >&2; exit 1; }
done
mkdir -p "$WORK" "$QC" "$METRICS" "$FOD" "$REG" "$TRACKS" "$SCRATCH" "$ACT"
if [[ "$interactive" == true ]]; then
    command -v mrview >/dev/null || { echo "interactive=true but mrview was not found in PATH" >&2; exit 1; }
fi
mrconvert -version > "$QC/mrtrix_version.txt" 2>&1
5ttgen -version > "$QC/5ttgen_version.txt" 2>&1
command -v flirt >/dev/null || { echo "Missing FSL command: flirt" >&2; exit 1; }
if [[ "$age_less_than_60m" == true ]]; then
    for f in "$INFANT_SIF" "$FS_LICENSE_FILE"; do
        [[ -f "$f" ]] || { echo "Missing InfantFS dependency: $f" >&2; exit 1; }
    done
    command -v apptainer >/dev/null || { echo "Missing command: apptainer" >&2; exit 1; }
    mkdir -p "$INFANT_OUT"
    # Catch a container/API mismatch before DWI preprocessing takes hours.
    apptainer exec "$INFANT_SIF" infant_recon_all --help > "$QC/infant_recon_all_help.txt" 2>&1 || true
    infantfs_age_flag='--age'
    if awk -v age="$age_months" 'BEGIN { exit !(age == 0) }'; then
        infantfs_age_flag='--newborn'
    fi
    if ! grep -q -- "$infantfs_age_flag" "$QC/infant_recon_all_help.txt"; then
        echo "This infant_recon_all does not advertise $infantfs_age_flag; check its CLI before running the long pipeline." >&2
        exit 1
    fi
else
    for fslcmd in bet fast run_first_all; do
        command -v "$fslcmd" >/dev/null || { echo "Missing FSL command: $fslcmd" >&2; exit 1; }
    done
fi
if [[ -n "${FSLDIR:-}" && -f "$FSLDIR/etc/fslversion" ]]; then
    cp "$FSLDIR/etc/fslversion" "$QC/fsl_version.txt"
fi
# Save the configured script alongside results for reproducibility.
# cp /absolute/path/to/this_edited_script.sh "$QC/workflow_used.sh"

# ============================================================================
# 2. IMPORT ORIGINAL NIFTI + MATCHING JSON TOGETHER (PRESERVES PE AXIS HANDLING)
# ============================================================================
# Both JSONs must have signed PhaseEncodingDirection and valid TotalReadoutTime
# in seconds. If only EstimatedTotalReadoutTime exists, validate it and prepare
# a working JSON copy with TotalReadoutTime set BEFORE these imports.
# Likewise, put any approved readout override into that working JSON, so the
# DWI and its extracted TOPUP reference inherit the SAME value.
mrconvert "$DWI" "$WORK/ap_imported.mif" -fslgrad "$BVEC" "$BVAL" \
    -json_import "$JSON_AP" "${COMMON[@]}"
mrconvert "$REVERSE" "$WORK/pa_imported.mif" -fslgrad "$REVERSE_BVEC" "$REVERSE_BVAL" \
    -json_import "$JSON_PA" "${COMMON[@]}"
mrinfo "$WORK/ap_imported.mif" -shell_bvalues -shell_sizes "${COMMON[@]}"
mrinfo "$WORK/ap_imported.mif" -petable > "$QC/ap_pe_table.txt"
mrinfo "$WORK/pa_imported.mif" -petable > "$QC/pa_pe_table.txt"
mrinfo "$WORK/ap_imported.mif" -size -spacing -transform > "$QC/ap_geometry.txt"
mrinfo "$WORK/pa_imported.mif" -size -spacing -transform > "$QC/pa_geometry.txt"
# CHECK: compatible spatial grids/orientations, opposite PE directions after
# import, correct readout times. AP and PA readout times need not be identical.
# Do not fix real acquisition/grid differences by simply editing headers.

# ============================================================================
# 3. DENOISE AND REMOVE GIBBS RINGING BEFORE MOTION CORRECTION
# ============================================================================
dwidenoise "$WORK/ap_imported.mif" "$WORK/ap_denoised.mif" \
    -noise "$METRICS/noise.mif" "${COMMON[@]}"
mrcalc "$WORK/ap_imported.mif" "$WORK/ap_denoised.mif" -subtract \
    "$QC/denoise_residual.mif" "${COMMON[@]}"
mrdegibbs "$WORK/ap_denoised.mif" "$WORK/ap_degibbs.mif" \
    -axes "$GIBBS_AXES" "${COMMON[@]}"
# Do not run MP-PCA denoising on a reverse series containing only 1–2 volumes.
# Apply Gibbs correction to reverse images only with verified acquisition axes.
mrdegibbs "$WORK/pa_imported.mif" "$WORK/pa_degibbs.mif" \
    -axes "$GIBBS_AXES" "${COMMON[@]}"
show_mrview "$METRICS/noise.mif"
show_mrview "$QC/denoise_residual.mif"

# ============================================================================
# 4. TRUE-B0 TOPUP PAIR (FIRST ELIGIBLE B0 FROM EACH SERIES, NOT b=50)
# ============================================================================
dwiextract "$WORK/ap_degibbs.mif" "$WORK/ap_bzeros.mif" -bzero "${COMMON[@]}"
dwiextract "$WORK/pa_degibbs.mif" "$WORK/pa_bzeros.mif" -bzero "${COMMON[@]}"
mrconvert "$WORK/ap_bzeros.mif" "$WORK/ap_b0.mif" -coord 3 0 "${COMMON[@]}"
mrconvert "$WORK/pa_bzeros.mif" "$WORK/pa_b0.mif" -coord 3 0 "${COMMON[@]}"
mrcat "$WORK/ap_b0.mif" "$WORK/pa_b0.mif" "$WORK/se_epi_pair.mif" \
    -axis 3 "${COMMON[@]}"
mrinfo "$WORK/se_epi_pair.mif" -petable > "$QC/topup_pair_pe_table.txt"
show_mrview "$WORK/se_epi_pair.mif"
# CHECK: both references have usable signal, comparable contrast/TE/TR and
# opposite PE. Do not average unaligned b0s from a moving infant.

# ============================================================================
# 5. TOPUP + EDDY WITH QC OUTPUTS
# ============================================================================
# Header mode supports different AP/PA readout times. -align_seepi explicitly
# anchors TOPUP to the first DWI b0; use only with matching reference contrast.
dwifslpreproc "$WORK/ap_degibbs.mif" "$WORK/dwi_preprocessed.mif" \
    -rpe_header \
    -se_epi "$WORK/se_epi_pair.mif" \
    -align_seepi \
    -topup_options " --nthr=8" \
    -eddy_options "$EDDY_OPTIONS" \
    -eddyqc_all "$QC/eddy" \
    -scratch "$SCRATCH" \
    -config BZeroThreshold "$BZERO_THRESHOLD"
# CHECK saved EDDY motion/outlier outputs and corrected images. eddy_quad reports
# require eddy_quad and its dependencies to be installed; retaining files alone
# does not certify correction quality.
show_mrview "$WORK/dwi_preprocessed.mif"

# ============================================================================
# 6. BIAS CORRECTION, MASK, MEAN B0 AND FINAL GRADIENT EXPORT
# ============================================================================
dwibiascorrect ants "$WORK/dwi_preprocessed.mif" "$BIAS" \
    -bias "$METRICS/bias_field.mif" -scratch "$SCRATCH" "${BIAS_COMMON[@]}" -force
dwi2mask "$BIAS" "$MASK" "${COMMON[@]}"
dwiextract "$BIAS" "$WORK/bzeros_corrected.mif" -bzero "${COMMON[@]}"
mrmath "$WORK/bzeros_corrected.mif" mean "$REG/mean_b0.mif" -axis 3 "${COMMON[@]}"
mrconvert "$REG/mean_b0.mif" "$REG/mean_b0.nii.gz" "${COMMON[@]}"
mrconvert "$BIAS" "$OUT/dwi_preprocessed_biascorrected.nii.gz" \
    -export_grad_fsl "$OUT/final_rotated.bvec" "$OUT/final.bval" "${COMMON[@]}"
show_mrview "$REG/mean_b0.mif" -overlay.load "$MASK"
# Gradient QC is diagnostic here: do not blindly replace the corrected gradients.
dwigradcheck "$BIAS" -mask "$MASK" \
    -export_grad_fsl "$QC/gradcheck_candidate.bvec" "$QC/gradcheck_candidate.bval" \
    "${COMMON[@]}"
# If a real axis/sign error is confirmed, correct the input convention and rerun
# affected processing. The candidate files above are NOT final rotated gradients.

# ============================================================================
# 7. OPTIONAL AUXILIARY DTI METRICS (NOT USED BY THE CSD TRACTOGRAPHY)
# ============================================================================
# Skip this entire section if FA/MD/AD/RD are not needed. DTI_SHELLS above preserves
# the original all-shell analysis; selecting 0,1000 is a distinct analysis choice.
dwiextract "$BIAS" "$METRICS/dwi_dti.mif" -shells "$DTI_SHELLS" "${COMMON[@]}"
dwi2tensor "$METRICS/dwi_dti.mif" "$METRICS/tensor.mif" -mask "$MASK" "${COMMON[@]}"
tensor2metric "$METRICS/tensor.mif" -fa "$METRICS/fa.mif" -adc "$METRICS/md.mif" \
    -ad "$METRICS/ad.mif" -rd "$METRICS/rd.mif" "${COMMON[@]}"

# ============================================================================
# 8. THREE-TISSUE RESPONSE ESTIMATION, MSMT-CSD, NORMALIZATION
# ============================================================================
dwi2response dhollander "$BIAS" "$FOD/wm_response.txt" "$FOD/gm_response.txt" \
    "$FOD/csf_response.txt" -mask "$MASK" -fa "$DHOLLANDER_FA" \
    -voxels "$QC/response_voxels.mif" -scratch "$SCRATCH" "${COMMON[@]}"
show_mrview "$REG/mean_b0.mif" -overlay.load "$QC/response_voxels.mif"
# CHECK response selections and curves before proceeding. Three-tissue signal
# separation is not guaranteed by the shell count or by successful execution.
dwi2fod msmt_csd "$BIAS" \
    "$FOD/wm_response.txt" "$FOD/wm_fod.mif" \
    "$FOD/gm_response.txt" "$FOD/gm.mif" \
    "$FOD/csf_response.txt" "$FOD/csf.mif" -mask "$MASK" "${COMMON[@]}"
mtnormalise "$FOD/wm_fod.mif" "$WM_FOD" \
    "$FOD/gm.mif" "$FOD/gm_normalised.mif" \
    "$FOD/csf.mif" "$FOD/csf_normalised.mif" -mask "$MASK" "${COMMON[@]}"
show_mrview "$REG/mean_b0.mif" -odf.load_sh "$WM_FOD"
# For cohort FOD-amplitude/fixel comparisons, specify a shared response-function
# strategy before processing the cohort. This runbook retains subject responses.

# ============================================================================
# 9. AGE-SELECTED ANATOMICAL SEGMENTATION / 5TT GENERATION
# ============================================================================
if [[ "$age_less_than_60m" == true ]]; then
    # InfantFS age is in months; the age argument selects the age-specific model.
    # Bind the license and data in, and retain the established output layout.
    INFANTFS_AGE_ARGS=(--age "$age_months")
    if awk -v age="$age_months" 'BEGIN { exit !(age == 0) }'; then
        # Stable InfantFS 8.1.0 mishandles --age 0; its documented newborn mode
        # selects the dedicated newborn training set instead.
        INFANTFS_AGE_ARGS=(--newborn)
    fi
    apptainer exec \
        --bind "$FS_LICENSE_FILE:/usr/local/freesurfer/license.txt:ro" \
        --bind "$INPUT:$INPUT:ro" --bind "$INFANT_OUT:$INFANT_OUT:rw" \
        "$INFANT_SIF" infant_recon_all -s "$SUBJECT" -i "$T1" \
        -o "$INFANT_OUT" "${INFANTFS_AGE_ARGS[@]}"
    apptainer exec "$INFANT_SIF" cat "$CONTAINER_LUT" > "$LUT_RAW"

    # InfantFS LUTs may add a seventh sclimbic class field. Normalize a copy for
    # MRtrix labelconvert while retaining the source LUT unchanged.
    python3 - "$LUT_RAW" "$LUT_MRTRIX" <<'PYLUT'
from pathlib import Path
import sys
source, destination = map(Path, sys.argv[1:])
normalized = []
removed = 0
for line_number, line in enumerate(source.read_text().splitlines(keepends=True), 1):
    body, marker, comment = line.partition('#')
    fields = body.split()
    if not fields or len(fields) == 6:
        normalized.append(line)
        continue
    if len(fields) == 7 and fields[-1] in {'2', '3'}:
        ending = '\r\n' if line.endswith('\r\n') else ('\n' if line.endswith('\n') else '')
        suffix = f' #{comment.rstrip(chr(13)+chr(10))}' if marker else ''
        normalized.append(' '.join(fields[:6]) + suffix + ending)
        removed += 1
        continue
    raise SystemExit(f"Unsupported LUT row at {source}:{line_number}: {len(fields)} fields: {line.rstrip()}")
destination.parent.mkdir(parents=True, exist_ok=True)
destination.write_text(''.join(normalized))
print(f"Validated LUT: {destination} ({removed} appended class column(s) removed)")
PYLUT
    [[ -s "$INFANT_SEG" && -s "$INFANT_REF" ]] || {
        echo "InfantFS expected outputs missing; verify INFANT_SEG and INFANT_REF for this container." >&2
        exit 1
    }
    mrconvert "$INFANT_REF" "$ACT/structural_reference.nii.gz" "${COMMON[@]}"
    5ttgen freesurfer "$INFANT_SEG" "$FIVE_TT" \
        -lut "$LUT_MRTRIX" -nocrop -sgm_amyg_hipp \
        -scratch "$SCRATCH" "${COMMON[@]}"
else
    mrconvert "$T1" "$ACT/T1_raw.mif" "${COMMON[@]}"
    5ttgen fsl "$ACT/T1_raw.mif" "$FIVE_TT" \
        -nocrop -scratch "$SCRATCH" "${COMMON[@]}"
    mrconvert "$ACT/T1_raw.mif" "$ACT/structural_reference.nii.gz" "${COMMON[@]}"
fi

# ============================================================================
# 10. STRUCTURAL FIVE-TISSUE IMAGE AND QC
# ============================================================================
mrinfo "$FIVE_TT" -size
5ttcheck "$FIVE_TT" "${COMMON[@]}"
5tt2vis "$FIVE_TT" "$ACT/5tt_vis.mif" "${COMMON[@]}"
show_mrview "$ACT/structural_reference.nii.gz" -overlay.load "$ACT/5tt_vis.mif"
# Structural QC must check anatomy, not merely five volumes / successful 5ttcheck.
mrconvert "$FIVE_TT" "$ACT/cgm.mif" -coord 3 0 -axes 0,1,2 "${COMMON[@]}"
mrconvert "$FIVE_TT" "$ACT/sgm.mif" -coord 3 1 -axes 0,1,2 "${COMMON[@]}"
mrconvert "$FIVE_TT" "$ACT/wm.mif" -coord 3 2 -axes 0,1,2 "${COMMON[@]}"
mrconvert "$FIVE_TT" "$ACT/csf.mif" -coord 3 3 -axes 0,1,2 "${COMMON[@]}"
mrconvert "$FIVE_TT" "$ACT/pathology.mif" -coord 3 4 -axes 0,1,2 "${COMMON[@]}"

# ============================================================================
# 11. RIGID REGISTRATION + MATRIX CONVERSION + 5TT INTO DWI SPACE
# ============================================================================
flirt -in "$REG/mean_b0.nii.gz" -ref "$ACT/structural_reference.nii.gz" \
    -out "$REG/mean_b0_in_struct.nii.gz" -omat "$REG/diff2struct_fsl.mat" \
    -dof 6 -cost mutualinfo -interp trilinear
show_mrview "$ACT/structural_reference.nii.gz" -overlay.load "$REG/mean_b0_in_struct.nii.gz"
transformconvert "$REG/diff2struct_fsl.mat" "$REG/mean_b0.nii.gz" \
    "$ACT/structural_reference.nii.gz" flirt_import "$REG/diff2struct_mrtrix.txt"
# Header transform only: retain structural voxel resolution. No -template needed.
mrtransform "$FIVE_TT" "$REG/5tt_coreg.mif" \
    -linear "$REG/diff2struct_mrtrix.txt" -inverse "${COMMON[@]}"
5ttcheck "$REG/5tt_coreg.mif" "${COMMON[@]}"
5tt2vis "$REG/5tt_coreg.mif" "$REG/5tt_coreg_vis.mif" "${COMMON[@]}"
mrconvert "$REG/5tt_coreg.mif" "$REG/wm_coreg.mif" -coord 3 2 -axes 0,1,2 "${COMMON[@]}"
show_mrview "$REG/mean_b0.mif" -overlay.load "$REG/5tt_coreg_vis.mif"
5tt2gmwmi "$REG/5tt_coreg.mif" "$REG/gmwm_seed_coreg.mif" "${COMMON[@]}"
show_mrview "$REG/mean_b0.mif" -overlay.load "$REG/gmwm_seed_coreg.mif"

# ============================================================================
# 12. ACT PILOT — INSPECT BEFORE THE LARGE RUN
# ============================================================================
tckgen "$WM_FOD" "$TRACKS/tracks_ACT_pilot.tck" -algorithm iFOD2 \
    -act "$REG/5tt_coreg.mif" -backtrack \
    -seed_gmwmi "$REG/gmwm_seed_coreg.mif" -select "$TRACK_PILOT" \
    "${TRACK_OPTIONS[@]}" "${COMMON[@]}"
tckstats "$TRACKS/tracks_ACT_pilot.tck" > "$QC/ACT_pilot_stats.txt"
show_mrview "$REG/mean_b0.mif" -tractography.load "$TRACKS/tracks_ACT_pilot.tck"
# CHECK coverage, spurious pathways, endpoints, and actual accepted count.
# Do not treat streamline count or attractive rendering as biological validation.

# ============================================================================
# 13. FINAL ACT TRACTOGRAM — SAME PARAMETERS AS ACCEPTED PILOT
# ============================================================================
tckgen "$WM_FOD" "$TRACKS/tracks_ACT_final.tck" -algorithm iFOD2 \
    -act "$REG/5tt_coreg.mif" -backtrack \
    -seed_gmwmi "$REG/gmwm_seed_coreg.mif" -select "$TRACK_FINAL" \
    "${TRACK_OPTIONS[@]}" "${COMMON[@]}"
tckstats "$TRACKS/tracks_ACT_final.tck" > "$QC/ACT_final_stats.txt"

# ============================================================================
# 14. OPTIONAL NON-ACT EXPLORATORY TRACTOGRAM (RUN MANUALLY IF NEEDED)
# ============================================================================
# tckgen "$WM_FOD" "$TRACKS/whole_brain_nonACT.tck" -algorithm iFOD2 \
#     -seed_dynamic "$WM_FOD" -mask "$MASK" -select "$TRACK_PILOT" \
#     "${TRACK_OPTIONS[@]}" "${COMMON[@]}"
# tckstats "$TRACKS/whole_brain_nonACT.tck"
# This comparison also changes seeding and ACT-dependent tracking defaults;
# it does not isolate the effect of ACT.

# Optional future quantitative connectivity stage (not automatically executed):
# tcksift2 "$TRACKS/tracks_ACT_final.tck" "$WM_FOD" "$TRACKS/sift2_weights.txt" \
#     -act "$REG/5tt_coreg.mif" -out_mu "$TRACKS/sift2_mu.txt" "${COMMON[@]}"
# A connectome additionally needs a validated neonatal parcellation aligned to
# DWI space and an explicit edge definition. Streamline counts are not axon counts.
# This FSL 5TT workflow does not generate a parcellation atlas.
