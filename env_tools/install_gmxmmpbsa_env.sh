#!/bin/bash
# Rebuild the gmxMMPBSA conda env from scratch -- e.g. after a scratch
# purge wipes it and no Acacia backup exists to restore from (check
# restore_env_from_acacia.sh first: `rclone lsd acacia:volet-envs` -- a
# restore is much faster than reinstalling if a backup exists).
#
# Mirrors gmx_MMPBSA's own official conda-based install method, plus the
# exact package pins this project validated working (2026-08-14, GSTA1-1 +
# lig_03506 ran clean, 0 errors/warnings). Builds on scratch, not
# /software -- the file-count quota there is what caused the original
# install saga (see docs/implementation_notes.md). Do NOT install AmberTools
# via /software or add `pip install --user` -- both were real bugs hit and
# fixed during the original install.
#
# STAGE 2 (AmberTools via conda) is the step most likely to need attention
# if package versions have moved on since this was last validated -- watch
# it and paste back any errors.
#
# Usage:
#   bash install_gmxmmpbsa_env.sh [env_dir]
# Default env_dir: $MYSCRATCH/conda_envs/gmxMMPBSA
set -euo pipefail

ENV_DIR=${1:-$MYSCRATCH/conda_envs/gmxMMPBSA}

if [ -d "$ENV_DIR" ] && [ -n "$(ls -A "$ENV_DIR" 2>/dev/null)" ]; then
    echo "ERROR: $ENV_DIR already has content -- refusing to overwrite."
    echo "Remove it first if you really want to rebuild from scratch, or use"
    echo "restore_env_from_acacia.sh instead if a backup exists."
    exit 1
fi

source $MYSOFTWARE/miniconda3/etc/profile.d/conda.sh

echo "======================================================================"
echo "STAGE 1: create conda env at $ENV_DIR (python 3.11)"
echo "======================================================================"
conda create -p "$ENV_DIR" python=3.11 -y

conda activate "$ENV_DIR"

echo ""
echo "======================================================================"
echo "STAGE 2: mpi4py + AmberTools via conda-forge"
echo "======================================================================"
# ambertools<24 pin: this project's validated working install used this
# constraint. gmx_MMPBSA's compatibility with newer AmberTools hasn't been
# re-tested -- if this fails or pulls something clearly broken, that's the
# first thing to investigate.
conda install -c conda-forge "mpi4py=4.0.1" "ambertools<24" -y

echo ""
echo "======================================================================"
echo "STAGE 3: gmx_MMPBSA via pip (into the active env)"
echo "======================================================================"
# No --user here: it installs outside the active env (into PYTHONUSERBASE)
# and breaks `gmx_MMPBSA` on PATH -- a real bug hit and fixed during the
# original install (see docs/implementation_notes.md).
python3 -m pip install gmx_MMPBSA

echo ""
echo "======================================================================"
echo "STAGE 4: sanity check"
echo "======================================================================"
python3 -c "import GMXMMPBSA; print('[OK] GMXMMPBSA importable')"
echo "[OK] gmx_MMPBSA resolves to: $(command -v gmx_MMPBSA)"

echo ""
echo "[OK] gmxMMPBSA env rebuilt at $ENV_DIR"
echo ""
echo "This does NOT prove it works end-to-end yet -- re-run the validated"
echo "smoke test before trusting it for real work:"
echo "  bash md/05_run_mmpbsa.sh --run-dir md/runs/gsta1_lig03506 --start-frame 500 --end-frame 520"
echo ""
echo "Once confirmed, back it up so the next purge doesn't cost you this again:"
echo "  bash env_tools/backup_env_to_acacia.sh \"$ENV_DIR\" acacia:volet-envs/gmxMMPBSA --do-it"
