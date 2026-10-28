#!/usr/bin/env python3
"""Single-command Python runner matching mrtrix3_full_workflow_age_agnostic.sh.

Edit the CONFIGURATION section, activate the MRtrix3/FSL/ANTs environment, then run:
    python3 mrtrix3_full_workflow_age_agnostic.py
This orchestrates the same external MRtrix3, FSL, ANTs, and (for InfantFS) Apptainer
commands as the Bash CLI. It does not replace visual or methodological QC.
"""
from __future__ import annotations

import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
from typing import Sequence

# ============================================================================
# CONFIGURATION — edit these values to match the Bash script's section 0
# ============================================================================
PROJECT = Path.home() / "projects" / "mrtrix3"
SUBJECT = "sub-XXXX"
SESSION = "ses-YYYY"

# Branch choice: True = Infant FreeSurfer for an age below 60 months;
# False = standard FSL 5ttgen route (normally age 60 months or older).
AGE_LESS_THAN_60M = False
AGE_MONTHS: float | None = None  # Required when InfantFS is selected; 0 means newborn.
INTERACTIVE = False  # True opens MRView QC checkpoints; False is headless.

OUT = PROJECT / "out_age_switch" / SUBJECT / SESSION
INFANT_SIF = Path.home() / "containers" / "freesurfer_infant-8.1.0.sif"
FS_LICENSE_FILE = Path.home() / "license.txt"
THREADS = 8

BZERO_THRESHOLD = 10
GIBBS_AXES = "0,1"  # Verify the acquisition plane in MRtrix image axes.
EDDY_OPTIONS = "--repol --data_is_shelled"
DHOLLANDER_FA = 0.2
DTI_SHELLS = "0,500,1000,2000,3000"
TRACK_PILOT = 200_000
TRACK_FINAL = 10_000_000
FOD_CUTOFF = 0.1
MIN_LENGTH_MM: float | None = None
MAX_LENGTH_MM: float | None = None
STEP_MM: float | None = None
CONTAINER_LUT = "/usr/local/freesurfer/luts/FreeSurferColorLUT.txt"

# Input files follow the same BIDS-like names as the Bash CLI.
INPUT = PROJECT / "in" / SUBJECT / SESSION
DWI = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-AP_dwi.nii.gz"
BVAL = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-AP_dwi.bval"
BVEC = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-AP_dwi.bvec"
JSON_AP = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-AP_dwi.json"
REVERSE = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-PA_dwi.nii.gz"
REVERSE_BVAL = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-PA_dwi.bval"
REVERSE_BVEC = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-PA_dwi.bvec"
JSON_PA = INPUT / "dwi" / f"{SUBJECT}_{SESSION}_dir-PA_dwi.json"
T1 = INPUT / "anat" / f"{SUBJECT}_{SESSION}_T1w.nii.gz"

# ============================================================================
# Runner
# ============================================================================

def normalize_freesurfer_lut(source: Path, destination: Path) -> None:
    """Normalize only known 7-column InfantFS rows for MRtrix labelconvert."""
    normalized: list[str] = []
    removed = 0
    for line_number, line in enumerate(source.read_text().splitlines(keepends=True), 1):
        body, marker, comment = line.partition("#")
        fields = body.split()
        if not fields or len(fields) == 6:
            normalized.append(line)
            continue
        if len(fields) == 7 and fields[-1] in {"2", "3"}:
            ending = "\r\n" if line.endswith("\r\n") else ("\n" if line.endswith("\n") else "")
            suffix = f" #{comment.rstrip(chr(13) + chr(10))}" if marker else ""
            normalized.append(" ".join(fields[:6]) + suffix + ending)
            removed += 1
            continue
        raise RuntimeError(
            f"Unsupported LUT row at {source}:{line_number}: {len(fields)} fields: {line.rstrip()}"
        )
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text("".join(normalized))
    print(f"Validated LUT: {destination} ({removed} appended class column(s) removed)", flush=True)


class Workflow:
    def __init__(self) -> None:
        self.env = os.environ.copy()
        self.env["FSLOUTPUTTYPE"] = "NIFTI_GZ"
        self.env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(THREADS)
        self.env["OMP_NUM_THREADS"] = str(THREADS)
        self.log_path = OUT / "workflow.log"
        self.common = ["-nthreads", str(THREADS), "-config", "BZeroThreshold", str(BZERO_THRESHOLD)]
        # Keep -nthreads out of dwibiascorrect's N4 invocation; use MRtrix config
        # and the ITK environment thread setting instead, matching the Bash script.
        self.bias_common = ["-config", "BZeroThreshold", str(BZERO_THRESHOLD),
                            "-config", "NumberOfThreads", str(THREADS)]
        self.work = OUT / "work"
        self.qc = OUT / "qc"
        self.metrics = OUT / "metrics"
        self.fod = OUT / "fod"
        self.act = OUT / "act"
        self.infant_out = OUT / "infantfs"
        self.reg = self.act / "registration"
        self.tracks = OUT / "tractography"
        self.scratch = OUT / "scratch"
        self.bias = OUT / "dwi_preprocessed_biascorrected.mif"
        self.mask = OUT / "brain_mask.mif"
        self.wm_fod = self.fod / "wm_fod_normalised.mif"
        self.five_tt = self.act / ("5tt_infantfs.mif" if AGE_LESS_THAN_60M else "5tt_fsl.mif")

    def run(self, args: Sequence[object], *, capture: Path | None = None,
            check: bool = True, append_log: bool = True) -> str:
        command = [str(arg) for arg in args]
        print("Command:  " + shlex.join(command), flush=True)
        if append_log:
            self.log_path.parent.mkdir(parents=True, exist_ok=True)
            log = self.log_path.open("a", encoding="utf-8")
        else:
            log = None
        capture_stream = capture.open("w", encoding="utf-8") if capture else None
        try:
            if log:
                from datetime import datetime, timezone
                log.write(f"\n[{datetime.now(timezone.utc).isoformat()}] {shlex.join(command)}\n")
                log.flush()
            proc = subprocess.Popen(command, env=self.env, stdout=subprocess.PIPE,
                                    stderr=subprocess.STDOUT, text=True, bufsize=1)
            assert proc.stdout is not None
            for line in proc.stdout:
                sys.stdout.write(line)
                sys.stdout.flush()
                if log:
                    log.write(line)
                    log.flush()
                if capture_stream:
                    capture_stream.write(line)
                    capture_stream.flush()
            code = proc.wait()
            if code and check:
                raise subprocess.CalledProcessError(code, command)
            return ""
        finally:
            if capture_stream:
                capture_stream.close()
            if log:
                log.close()

    def view(self, *args: object) -> None:
        if INTERACTIVE:
            self.run(["mrview", *args])

    def preflight(self) -> None:
        if not isinstance(INTERACTIVE, bool):
            raise ValueError("INTERACTIVE must be True or False")
        if THREADS < 1:
            raise ValueError("THREADS must be a positive integer")
        if AGE_LESS_THAN_60M:
            if AGE_MONTHS is None or not math.isfinite(float(AGE_MONTHS)) or not 0 <= AGE_MONTHS < 60:
                raise ValueError("Set AGE_MONTHS to a numeric age >=0 and <60 for Infant FreeSurfer")
        elif AGE_MONTHS is not None and (not math.isfinite(float(AGE_MONTHS)) or AGE_MONTHS < 60):
            raise ValueError("For FSL 5ttgen, omit AGE_MONTHS or set it to >=60")
        for path in (DWI, BVAL, BVEC, JSON_AP, REVERSE, REVERSE_BVAL, REVERSE_BVEC, JSON_PA, T1):
            if not path.is_file():
                raise FileNotFoundError(f"Missing input: {path}")
        self.out_dirs = (self.work, self.qc, self.metrics, self.fod, self.reg,
                         self.tracks, self.scratch, self.act)
        for directory in self.out_dirs:
            directory.mkdir(parents=True, exist_ok=True)
        if INTERACTIVE and not shutil.which("mrview"):
            raise RuntimeError("INTERACTIVE=True but mrview was not found in PATH")
        required = ["mrconvert", "5ttgen", "flirt", "dwidenoise", "mrcalc", "mrdegibbs",
                    "dwiextract", "mrinfo", "mrcat", "dwifslpreproc", "dwibiascorrect",
                    "dwi2mask", "mrmath", "dwi2response", "dwi2fod", "mtnormalise",
                    "5ttcheck", "5tt2vis", "5tt2gmwmi", "transformconvert", "mrtransform",
                    "tckgen", "tckstats", "dwigradcheck", "dwi2tensor", "tensor2metric"]
        if AGE_LESS_THAN_60M:
            required.append("apptainer")
            for path in (INFANT_SIF, FS_LICENSE_FILE):
                if not path.is_file():
                    raise FileNotFoundError(f"Missing Infant FreeSurfer dependency: {path}")
        else:
            required.extend(["bet", "fast", "run_first_all"])
        missing = [name for name in required if shutil.which(name) is None]
        if missing:
            raise RuntimeError("Required commands not found in active environment: " + ", ".join(missing))

        if AGE_LESS_THAN_60M:
            (self.infant_out).mkdir(parents=True, exist_ok=True)
            help_file = self.qc / "infant_recon_all_help.txt"
            self.run(["apptainer", "exec", INFANT_SIF, "infant_recon_all", "--help"],
                     capture=help_file, check=False)
            flag = "--newborn" if AGE_MONTHS == 0 else "--age"
            if flag not in help_file.read_text(errors="replace"):
                raise RuntimeError(f"InfantFS container help does not advertise {flag}; check its CLI")

    def provenance(self) -> None:
        for command, destination in ((["mrconvert", "-version"], self.qc / "mrtrix_version.txt"),
                                     (["5ttgen", "-version"], self.qc / "5ttgen_version.txt")):
            self.run(command, capture=destination)
        fsldir = self.env.get("FSLDIR")
        if fsldir and (Path(fsldir) / "etc/fslversion").is_file():
            shutil.copy2(Path(fsldir) / "etc/fslversion", self.qc / "fsl_version.txt")

    def pipeline(self) -> None:
        w, q, m, f, act, reg, tr, scratch = (
            self.work, self.qc, self.metrics, self.fod, self.act, self.reg, self.tracks, self.scratch
        )
        common = self.common
        # 2. Import original NIfTI with its matching JSON and FSL gradients.
        self.run(["mrconvert", DWI, w / "ap_imported.mif", "-fslgrad", BVEC, BVAL,
                  "-json_import", JSON_AP, *common])
        self.run(["mrconvert", REVERSE, w / "pa_imported.mif", "-fslgrad", REVERSE_BVEC,
                  REVERSE_BVAL, "-json_import", JSON_PA, *common])
        self.run(["mrinfo", w / "ap_imported.mif", "-shell_bvalues", "-shell_sizes", *common])
        self.run(["mrinfo", w / "ap_imported.mif", "-petable"], capture=q / "ap_pe_table.txt")
        self.run(["mrinfo", w / "pa_imported.mif", "-petable"], capture=q / "pa_pe_table.txt")
        self.run(["mrinfo", w / "ap_imported.mif", "-size", "-spacing", "-transform"], capture=q / "ap_geometry.txt")
        self.run(["mrinfo", w / "pa_imported.mif", "-size", "-spacing", "-transform"], capture=q / "pa_geometry.txt")

        # 3. Denoise and remove Gibbs ringing before motion correction.
        self.run(["dwidenoise", w / "ap_imported.mif", w / "ap_denoised.mif",
                  "-noise", m / "noise.mif", *common])
        self.run(["mrcalc", w / "ap_imported.mif", w / "ap_denoised.mif", "-subtract",
                  q / "denoise_residual.mif", *common])
        self.run(["mrdegibbs", w / "ap_denoised.mif", w / "ap_degibbs.mif",
                  "-axes", GIBBS_AXES, *common])
        self.run(["mrdegibbs", w / "pa_imported.mif", w / "pa_degibbs.mif",
                  "-axes", GIBBS_AXES, *common])
        self.view(m / "noise.mif")
        self.view(q / "denoise_residual.mif")

        # 4. Use the first eligible true b=0 from each series for TOPUP.
        self.run(["dwiextract", w / "ap_degibbs.mif", w / "ap_bzeros.mif", "-bzero", *common])
        self.run(["dwiextract", w / "pa_degibbs.mif", w / "pa_bzeros.mif", "-bzero", *common])
        self.run(["mrconvert", w / "ap_bzeros.mif", w / "ap_b0.mif", "-coord", 3, 0, *common])
        self.run(["mrconvert", w / "pa_bzeros.mif", w / "pa_b0.mif", "-coord", 3, 0, *common])
        self.run(["mrcat", w / "ap_b0.mif", w / "pa_b0.mif", w / "se_epi_pair.mif",
                  "-axis", 3, *common])
        self.run(["mrinfo", w / "se_epi_pair.mif", "-petable"], capture=q / "topup_pair_pe_table.txt")
        self.view(w / "se_epi_pair.mif")

        # 5. TOPUP + EDDY. Preserve the known-working Bash option placement:
        # -topup_options gets --nthr=8; do NOT pass MRtrix -nthreads here.
        self.run(["dwifslpreproc", w / "ap_degibbs.mif", w / "dwi_preprocessed.mif",
                  "-rpe_header", "-se_epi", w / "se_epi_pair.mif", "-align_seepi",
                  "-topup_options", " --nthr=8", "-eddy_options", f" {EDDY_OPTIONS.strip()}",
                  "-eddyqc_all", q / "eddy", "-scratch", scratch,
                  "-config", "BZeroThreshold", str(BZERO_THRESHOLD)])
        self.view(w / "dwi_preprocessed.mif")

        # 6. Bias correction, mask, mean b0 and corrected gradient export.
        self.run(["dwibiascorrect", "ants", w / "dwi_preprocessed.mif", self.bias,
                  "-bias", m / "bias_field.mif", "-scratch", scratch, *self.bias_common, "-force"])
        self.run(["dwi2mask", self.bias, self.mask, *common])
        self.run(["dwiextract", self.bias, w / "bzeros_corrected.mif", "-bzero", *common])
        self.run(["mrmath", w / "bzeros_corrected.mif", "mean", reg / "mean_b0.mif", "-axis", 3, *common])
        self.run(["mrconvert", reg / "mean_b0.mif", reg / "mean_b0.nii.gz", *common])
        self.run(["mrconvert", self.bias, OUT / "dwi_preprocessed_biascorrected.nii.gz",
                  "-export_grad_fsl", OUT / "final_rotated.bvec", OUT / "final.bval", *common])
        self.view(reg / "mean_b0.mif", "-overlay.load", self.mask)
        self.run(["dwigradcheck", self.bias, "-mask", self.mask,
                  "-export_grad_fsl", q / "gradcheck_candidate.bvec", q / "gradcheck_candidate.bval", *common])

        # 7. Optional auxiliary tensor metrics (same enabled default as Bash).
        self.run(["dwiextract", self.bias, m / "dwi_dti.mif", "-shells", DTI_SHELLS, *common])
        self.run(["dwi2tensor", m / "dwi_dti.mif", m / "tensor.mif", "-mask", self.mask, *common])
        self.run(["tensor2metric", m / "tensor.mif", "-fa", m / "fa.mif", "-adc", m / "md.mif",
                  "-ad", m / "ad.mif", "-rd", m / "rd.mif", *common])

        # 8. Dhollander three-tissue response, MSMT-CSD and normalisation.
        self.run(["dwi2response", "dhollander", self.bias, f / "wm_response.txt", f / "gm_response.txt",
                  f / "csf_response.txt", "-mask", self.mask, "-fa", DHOLLANDER_FA,
                  "-voxels", q / "response_voxels.mif", "-scratch", scratch, *common])
        self.view(reg / "mean_b0.mif", "-overlay.load", q / "response_voxels.mif")
        self.run(["dwi2fod", "msmt_csd", self.bias, f / "wm_response.txt", f / "wm_fod.mif",
                  f / "gm_response.txt", f / "gm.mif", f / "csf_response.txt", f / "csf.mif",
                  "-mask", self.mask, *common])
        self.run(["mtnormalise", f / "wm_fod.mif", self.wm_fod,
                  f / "gm.mif", f / "gm_normalised.mif", f / "csf.mif", f / "csf_normalised.mif",
                  "-mask", self.mask, *common])
        self.view(reg / "mean_b0.mif", "-odf.load_sh", self.wm_fod)

        # 9. Age-selected anatomical segmentation / 5TT generation.
        if AGE_LESS_THAN_60M:
            infant_seg = self.infant_out / "mri" / "aparc+aseg.mgz"
            infant_ref = self.infant_out / "mri" / "norm.nii.gz"
            age_args: list[object] = ["--newborn"] if AGE_MONTHS == 0 else ["--age", AGE_MONTHS]
            self.run(["apptainer", "exec", "--bind",
                      f"{FS_LICENSE_FILE}:/usr/local/freesurfer/license.txt:ro",
                      "--bind", f"{INPUT}:{INPUT}:ro", "--bind", f"{self.infant_out}:{self.infant_out}:rw",
                      INFANT_SIF, "infant_recon_all", "-s", SUBJECT, "-i", T1, "-o", self.infant_out,
                      *age_args])
            lut_raw = OUT / "FreeSurferColorLUT.txt"
            self.run(["apptainer", "exec", INFANT_SIF, "cat", CONTAINER_LUT], capture=lut_raw)
            lut_mrtrix = act / "FreeSurferColorLUT_mrtrix.txt"
            normalize_freesurfer_lut(lut_raw, lut_mrtrix)
            if not infant_seg.is_file() or not infant_ref.is_file():
                raise RuntimeError(f"InfantFS expected outputs missing: {infant_seg} or {infant_ref}")
            self.run(["mrconvert", infant_ref, act / "structural_reference.nii.gz", *common])
            self.run(["5ttgen", "freesurfer", infant_seg, self.five_tt, "-lut", lut_mrtrix,
                      "-nocrop", "-sgm_amyg_hipp", "-scratch", scratch, *common])
        else:
            self.run(["mrconvert", T1, act / "T1_raw.mif", *common])
            self.run(["5ttgen", "fsl", act / "T1_raw.mif", self.five_tt,
                      "-nocrop", "-scratch", scratch, *common])
            self.run(["mrconvert", act / "T1_raw.mif", act / "structural_reference.nii.gz", *common])

        # 10. Five-tissue image QC outputs.
        self.run(["mrinfo", self.five_tt, "-size"])
        self.run(["5ttcheck", self.five_tt, *common])
        self.run(["5tt2vis", self.five_tt, act / "5tt_vis.mif", *common])
        self.view(act / "structural_reference.nii.gz", "-overlay.load", act / "5tt_vis.mif")
        for tissue, index in (("cgm", 0), ("sgm", 1), ("wm", 2), ("csf", 3), ("pathology", 4)):
            self.run(["mrconvert", self.five_tt, act / f"{tissue}.mif", "-coord", 3, index,
                      "-axes", "0,1,2", *common])

        # 11. Rigid registration and 5TT to DWI space.
        self.run(["flirt", "-in", reg / "mean_b0.nii.gz", "-ref", act / "structural_reference.nii.gz",
                  "-out", reg / "mean_b0_in_struct.nii.gz", "-omat", reg / "diff2struct_fsl.mat",
                  "-dof", 6, "-cost", "mutualinfo", "-interp", "trilinear"])
        self.view(act / "structural_reference.nii.gz", "-overlay.load", reg / "mean_b0_in_struct.nii.gz")
        self.run(["transformconvert", reg / "diff2struct_fsl.mat", reg / "mean_b0.nii.gz",
                  act / "structural_reference.nii.gz", "flirt_import", reg / "diff2struct_mrtrix.txt"])
        self.run(["mrtransform", self.five_tt, reg / "5tt_coreg.mif", "-linear",
                  reg / "diff2struct_mrtrix.txt", "-inverse", *common])
        self.run(["5ttcheck", reg / "5tt_coreg.mif", *common])
        self.run(["5tt2vis", reg / "5tt_coreg.mif", reg / "5tt_coreg_vis.mif", *common])
        self.run(["mrconvert", reg / "5tt_coreg.mif", reg / "wm_coreg.mif", "-coord", 3, 2,
                  "-axes", "0,1,2", *common])
        self.view(reg / "mean_b0.mif", "-overlay.load", reg / "5tt_coreg_vis.mif")
        self.run(["5tt2gmwmi", reg / "5tt_coreg.mif", reg / "gmwm_seed_coreg.mif", *common])
        self.view(reg / "mean_b0.mif", "-overlay.load", reg / "gmwm_seed_coreg.mif")

        # 12–13. ACT pilot and final tractograms, with the Bash script parameters.
        track_options: list[object] = ["-cutoff", FOD_CUTOFF]
        if MIN_LENGTH_MM is not None:
            track_options += ["-minlength", MIN_LENGTH_MM]
        if MAX_LENGTH_MM is not None:
            track_options += ["-maxlength", MAX_LENGTH_MM]
        if STEP_MM is not None:
            track_options += ["-step", STEP_MM]
        for label, count in (("pilot", TRACK_PILOT), ("final", TRACK_FINAL)):
            track_path = self.tracks / f"tracks_ACT_{label}.tck"
            self.run(["tckgen", self.wm_fod, track_path, "-algorithm", "iFOD2", "-act",
                      reg / "5tt_coreg.mif", "-backtrack", "-seed_gmwmi", reg / "gmwm_seed_coreg.mif",
                      "-select", count, *track_options, *common])
            self.run(["tckstats", track_path], capture=self.qc / f"ACT_{label}_stats.txt")
            if label == "pilot":
                self.view(reg / "mean_b0.mif", "-tractography.load", track_path)

    def run_all(self) -> None:
        self.preflight()
        self.provenance()
        self.pipeline()


def main() -> int:
    try:
        Workflow().run_all()
    except subprocess.CalledProcessError as exc:
        print(f"\nFAILED (exit {exc.returncode}): {shlex.join([str(x) for x in exc.cmd])}", file=sys.stderr)
        return exc.returncode or 1
    except (Exception, KeyboardInterrupt) as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        return 130 if isinstance(exc, KeyboardInterrupt) else 1
    print(f"\nProcessing complete; review QC before interpretation: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
