#!/usr/bin/env bash
# ============================================================================
# Neonatal diffusion MRI + ACT environment installer
#
# Scope: local_mrtrix3_AP_PA_topup_fixed.py + Infant FreeSurfer .sif workflow.
# Includes: MRtrix3, FSL TOPUP / CPU EDDY / FLIRT / utilities / EDDY QC,
#           ANTs N4BiasFieldCorrection, Python QC utilities, conda-pack,
#           Apptainer (optional; enabled by default), optional CUDA EDDY.
#
# Internet-connected Ubuntu laptop:
#   bash setup_neonatal_mrtrix3_conda.sh install
#
# Transfer the resulting archive and this script to the OFFLINE Debian server:
#   bash setup_neonatal_mrtrix3_conda.sh unpack /path/to/mrtrix3_fsl_av_v2.tar.gz
#   bash setup_neonatal_mrtrix3_conda.sh check
#
# Customize with environment variables (examples):
#   ENV_PREFIX="$HOME/conda_envs/mrtrix3_fsl_av_v2"
#   ARCHIVE_DIR="$HOME/conda_env_archives"
#   INCLUDE_APPTAINER=0    # existing host apptainer will be used instead
#   ENABLE_EDDY_CUDA=1     # requires compatible NVIDIA driver on target
#   MAKE_ARCHIVE=0         # install only; pack later with 'pack' subcommand
#   MRTRIX_VERSION=3.0.8
#
# IMPORTANT:
# - Creates a NEW environment. Never modifies or deletes an existing prefix.
# - Does NOT install Infant FreeSurfer itself (already supplied as a .sif file).
# - Does NOT install the host's NVIDIA driver or configure user namespaces.
# - Building Linux x86_64 on Ubuntu and unpacking on Debian x86_64 requires
#   runtime compatibility testing on the actual server.
# - FSL and MRtrix3 use their respective licenses; review them before use.
# ============================================================================
set -Eeuo pipefail
IFS=$'\n\t'

ENV_PREFIX="${ENV_PREFIX:-${HOME}/conda_envs/mrtrix3_v1}"
ARCHIVE_DIR="${ARCHIVE_DIR:-${HOME}/conda_env_archives}"
ARCHIVE_NAME="${ARCHIVE_NAME:-mrtrix3_v1.tar.gz}"
MRTRIX_VERSION="${MRTRIX_VERSION:-3.0.8}"
INCLUDE_APPTAINER="${INCLUDE_APPTAINER:-1}"
ENABLE_EDDY_CUDA="${ENABLE_EDDY_CUDA:-0}"
MAKE_ARCHIVE="${MAKE_ARCHIVE:-1}"
FSL_CHANNEL='https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/'

log() { printf '\n[INFO] %s\n' "$*" >&2; }
warn() { printf '\n[WARN] %s\n' "$*" >&2; }
die() { printf '\n[ERROR] %s\n' "$*" >&2; exit 1; }

usage() {
  cat <<USAGE
Neonatal MRtrix3 / FSL / ANTs conda environment

Usage:
  bash $0 install               Create a new online conda environment; pack by default
  bash $0 pack                  Pack an already installed environment
  bash $0 unpack ARCHIVE.tar.gz Unpack on the offline Debian server; run conda-unpack
  bash $0 check                 Check an installed/unpacked environment

Configuration (environment variables):
  ENV_PREFIX        ${ENV_PREFIX}
  ARCHIVE_DIR       ${ARCHIVE_DIR}
  ARCHIVE_NAME      ${ARCHIVE_NAME}
  MRTRIX_VERSION    ${MRTRIX_VERSION}
  INCLUDE_APPTAINER ${INCLUDE_APPTAINER} (1=yes; 0=use host installation)
  ENABLE_EDDY_CUDA  ${ENABLE_EDDY_CUDA} (1=install optional fsl-eddy-cuda)
  MAKE_ARCHIVE      ${MAKE_ARCHIVE} (1=pack after install)

To use the environment after installing/unpacking:
  source "\${ENV_PREFIX}/activate_neonatal_mrtrix3.sh"

This installs the diffusion + InfantFS/ACT workflow, not the entire Aesclepian
Vision Streamlit/Celery/model-training application.
USAGE
}

assert_linux_x86_64() {
  [[ "$(uname -s)" == Linux ]] || die 'This installer supports Linux only.'
  [[ "$(uname -m)" == x86_64 ]] || die 'Build/target must be x86_64.'
}

find_conda() {
  if command -v conda >/dev/null 2>&1; then
    command -v conda
    return
  fi
  local path
  for path in \
    "${HOME}/miniforge3/bin/conda" \
    "${HOME}/anaconda3/bin/conda" \
    "${HOME}/miniconda3/bin/conda" \
    "${HOME}/mambaforge/bin/conda"; do
    if [[ -x "$path" ]]; then
      printf '%s\n' "$path"
      return
    fi
  done
  die 'Conda not found. Install Miniforge on the ONLINE build machine first.'
}

write_activation_helper() {
  local helper="${ENV_PREFIX}/activate_neonatal_mrtrix3.sh"
  cat > "$helper" <<'ACTIVATE'
#!/usr/bin/env bash
# Source, do not execute: source /path/to/env/activate_neonatal_mrtrix3.sh
_neonatal_env_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
if [[ -f "${_neonatal_env_dir}/bin/activate" ]]; then
  # conda-pack supplies bin/activate in the unpacked environment.
  source "${_neonatal_env_dir}/bin/activate"
else
  # A freshly created Conda prefix does not necessarily contain bin/activate.
  _neonatal_conda="${CONDA_EXE:-}"
  if [[ -z "${_neonatal_conda}" ]] && command -v conda >/dev/null 2>&1; then
    _neonatal_conda="$(command -v conda)"
  fi
  if [[ -z "${_neonatal_conda}" ]]; then
    for _candidate in \
      "${HOME}/miniforge3/bin/conda" \
      "${HOME}/anaconda3/bin/conda" \
      "${HOME}/miniconda3/bin/conda" \
      "${HOME}/mambaforge/bin/conda"; do
      if [[ -x "${_candidate}" ]]; then _neonatal_conda="${_candidate}"; break; fi
    done
  fi
  if [[ -z "${_neonatal_conda}" || ! -x "${_neonatal_conda}" ]]; then
    echo "Conda not found; cannot activate ${_neonatal_env_dir}" >&2
    return 1
  fi
  eval "$("${_neonatal_conda}" shell.bash hook)"
  conda activate "${_neonatal_env_dir}"
fi
export FSLDIR="${_neonatal_env_dir}"
export FSLOUTPUTTYPE=NIFTI_GZ
if [[ -f "${FSLDIR}/etc/fslconf/fsl.sh" ]]; then
  source "${FSLDIR}/etc/fslconf/fsl.sh"
fi
# Prioritize FSL's public command wrappers, without modifying system files.
if [[ -d "${FSLDIR}/share/fsl/bin" ]]; then
  export PATH="${FSLDIR}/share/fsl/bin:${PATH}"
fi
unset _neonatal_env_dir _neonatal_conda _candidate
ACTIVATE
  chmod 644 "$helper"
}

validate_env() (
  [[ -d "$ENV_PREFIX" ]] || die "Environment not found: $ENV_PREFIX"
  local helper="$ENV_PREFIX/activate_neonatal_mrtrix3.sh"
  [[ -f "$helper" ]] || die "Missing activation helper: $helper"
  # Activate only inside a subshell, never change the caller's environment.
  # shellcheck disable=SC1090
  set +u  # Third-party FSL/Conda activation scripts may refer to unset vars.
  source "$helper"
  set -u

  log "Environment: $ENV_PREFIX"
  log "FSLDIR: $FSLDIR"
  [[ -f "$FSLDIR/etc/fslconf/fsl.sh" ]] || die 'Missing FSL activation configuration (fsl-base).'
  [[ -f "$FSLDIR/etc/flirtsch/b02b0.cnf" ]] || die 'Missing FSL TOPUP b02b0.cnf configuration.'
  local failed=0 cmd
  local needed=(
    python mrconvert mrinfo mrview dwigradcheck dwidenoise mrdegibbs
    dwifslpreproc dwibiascorrect dwi2mask dwi2tensor tensor2metric
    dwi2response dwi2fod mtnormalise dwiextract mrmath mrcat
    tckgen tckstats 5ttgen 5ttcheck 5tt2vis 5tt2gmwmi
    mrtransform transformconvert mrstats mrcalc
    topup applytopup flirt N4BiasFieldCorrection
  )
  for cmd in "${needed[@]}"; do
    if command -v "$cmd" >/dev/null 2>&1; then
      printf '  [OK]      %-23s %s\n' "$cmd" "$(command -v "$cmd")"
    else
      printf '  [MISSING] %s\n' "$cmd" >&2
      failed=1
    fi
  done
  if command -v eddy_cpu >/dev/null 2>&1; then
    printf '  [OK]      eddy_cpu              %s\n' "$(command -v eddy_cpu)"
  elif command -v eddy_cuda >/dev/null 2>&1; then
    printf '  [OK]      eddy_cuda             %s\n' "$(command -v eddy_cuda)"
    warn 'CPU eddy is not visible; CUDA eddy will require GPU compatibility.'
  else
    printf '  [MISSING] Neither eddy_cpu nor eddy_cuda is visible.\n' >&2
    failed=1
  fi
  if command -v eddy_cuda >/dev/null 2>&1; then
    printf '  [OPTIONAL] eddy_cuda available\n'
  fi
  if command -v apptainer >/dev/null 2>&1; then
    printf '  [OK]      apptainer             %s\n' "$(command -v apptainer)"
  else
    warn 'apptainer is missing. Infant FreeSurfer .sif cannot run until it is installed.'
    failed=1
  fi
  python -c 'import numpy, scipy, nibabel, pandas, matplotlib; print("[OK] Python QC packages")' || failed=1
  if [[ $failed -ne 0 ]]; then
    die 'One or more requirements are missing; fix these before running the pipeline.'
  fi
  log 'All command-presence checks passed.'
  log 'NOTE: This does not validate TOPUP/EDDY output, registrations, anatomical QC, or GPU runtime.'
)

make_pack() {
  assert_linux_x86_64
  [[ -d "$ENV_PREFIX" ]] || die "Missing environment: $ENV_PREFIX"
  mkdir -p "$ARCHIVE_DIR"
  local archive="$ARCHIVE_DIR/$ARCHIVE_NAME"
  [[ ! -e "$archive" ]] || die "Archive already exists; refusing overwrite: $archive"
  [[ -x "$ENV_PREFIX/bin/conda-pack" ]] || die 'conda-pack missing from environment.'
  log "Packing environment into $archive"
  "$ENV_PREFIX/bin/conda-pack" --prefix "$ENV_PREFIX" --output "$archive" --format tar.gz
  (cd "$ARCHIVE_DIR" && sha256sum "$ARCHIVE_NAME" > "$ARCHIVE_NAME.sha256")
  log "Archive complete: $archive"
  log "Checksum: ${archive}.sha256"
  log 'Transfer the tar.gz, its .sha256, and this .sh file to the offline Debian server.'
}

install_env() {
  assert_linux_x86_64
  [[ ! -e "$ENV_PREFIX" ]] || die "Refusing to modify existing environment: $ENV_PREFIX"
  local conda_bin
  conda_bin="$(find_conda)"
  mkdir -p "$(dirname "$ENV_PREFIX")" "$ARCHIVE_DIR"
  log "Using Conda: $conda_bin"
  log "Creating isolated environment: $ENV_PREFIX"

  # Official FSL channel should come before conda-forge; MRtrix3 comes from
  # Bioconda. Using --override-channels avoids an unintended defaults channel.
  local channels=(
    --override-channels
    -c "$FSL_CHANNEL"
    -c conda-forge
    -c bioconda
  )
  local packages=(
    "python=3.12"
    "mrtrix3=${MRTRIX_VERSION}"
    ants
    fsl-base
    fsl-topup
    fsl-eddy
    fsl-flirt
    fsl-avwutils
    fsl-utils
    fsl-misc_scripts
    fsl-eddy_qc
    numpy scipy nibabel pandas matplotlib pyyaml
    conda-pack
  )
  if [[ "$ENABLE_EDDY_CUDA" == 1 ]]; then
    packages+=(fsl-eddy-cuda)
    warn 'CUDA EDDY is included; test against the actual target GPU and NVIDIA driver.'
  fi
  if [[ "$INCLUDE_APPTAINER" == 1 ]]; then
    packages+=(apptainer)
  else
    warn 'Conda Apptainer disabled; you must have a working host apptainer.'
  fi

  "$conda_bin" create -y --prefix "$ENV_PREFIX" "${channels[@]}" "${packages[@]}"
  [[ -d "$ENV_PREFIX" ]] || die 'Conda reported success but environment does not exist.'

  write_activation_helper
  # Store a package inventory in the archive for reproducibility.
  "$conda_bin" list --prefix "$ENV_PREFIX" --explicit > "$ENV_PREFIX/packages-explicit.txt"
  "$conda_bin" list --prefix "$ENV_PREFIX" > "$ENV_PREFIX/packages-readable.txt"
  printf 'created_on=%s\nbuild_os=%s\nbuild_arch=%s\n' \
    "$(date -Is)" "$(uname -sr)" "$(uname -m)" \
    > "$ENV_PREFIX/build-metadata.txt"

  validate_env
  if [[ "$MAKE_ARCHIVE" == 1 ]]; then
    make_pack
  else
    log 'To create an archive later: bash setup_neonatal_mrtrix3_conda.sh pack'
  fi
  log 'Installation finished without modifying other Conda environments.'
}

unpack_env() {
  assert_linux_x86_64
  local archive="${1:-}"
  [[ -n "$archive" ]] || die 'Usage: unpack /path/to/archive.tar.gz'
  [[ -f "$archive" ]] || die "Archive not found: $archive"
  [[ ! -e "$ENV_PREFIX" ]] || die "Refusing to overwrite existing environment: $ENV_PREFIX"

  local checksum="${archive}.sha256"
  if [[ -f "$checksum" ]]; then
    log "Verifying checksum from $checksum"
    (cd "$(dirname "$archive")" && sha256sum --check "$(basename "$checksum")")
  else
    warn "No checksum found at $checksum; cannot verify archive integrity."
  fi

  mkdir -p "$ENV_PREFIX"
  log "Unpacking $archive to $ENV_PREFIX"
  tar -xzf "$archive" -C "$ENV_PREFIX"
  if [[ ! -x "$ENV_PREFIX/bin/conda-unpack" ]]; then
    die 'conda-unpack not found in archive; not a valid conda-pack archive.'
  fi
  # Repair absolute prefixes to the new path. Must be done only after unpack.
  # Use the pack-provided activation script; no Conda solver/internet required.
  (
    # shellcheck disable=SC1090
    set +u
    source "$ENV_PREFIX/bin/activate"
    set -u
    conda-unpack
  )
  log 'Conda prefixes repaired.'
  validate_env
  log "Offline deployment completed: $ENV_PREFIX"
  log "Use: source '$ENV_PREFIX/activate_neonatal_mrtrix3.sh'"
}

case "${1:-help}" in
  install) [[ $# -eq 1 ]] || die 'Usage: install'; install_env ;;
  pack)    [[ $# -eq 1 ]] || die 'Usage: pack'; make_pack ;;
  unpack)  [[ $# -eq 2 ]] || die 'Usage: unpack ARCHIVE.tar.gz'; unpack_env "$2" ;;
  check)   [[ $# -eq 1 ]] || die 'Usage: check'; validate_env ;;
  help|-h|--help) usage ;;
  *) usage; die "Unknown command: $1" ;;
esac
