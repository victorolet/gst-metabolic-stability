#!/bin/bash --login
# MM-PBSA rescoring of a completed 03_run_md.sh production run, via
# gmx_MMPBSA (Single Trajectory protocol). Receptor = Protein + GSH,
# ligand = docked candidate (UNL) -- a modelling choice, see
# docs/implementation_notes.md. Uses PB (not GB) with CHARMM radii
# (PBRadii=7, radiopt=0), per gmx_MMPBSA's own CHARMM guidance.
#
# What this does: builds a Protein_LIG index group (receptor), a full-
# density PBC-corrected trajectory, writes mmpbsa.in, stages topol.top +
# itp/prm/posre files + the FF dir alongside the output dir (gmx_MMPBSA
# re-parses topol.top via ParmEd), then runs gmx_MMPBSA.
#
# Frame range defaults to a small validation slice (20 frames) -- PB
# electrostatics is expensive per frame. Widen --start-frame/--end-frame/
# --interval once this runs clean (50-200+ frames is typical for a
# converged estimate). For real parallel scale-out, MPI mode is blocked on
# Setonix (see notes doc) -- use slurm/submit_mmpbsa_array.sh instead.
#
# Validated 2026-08-14: 20-frame slice (500-520), 0 errors/warnings,
# ~4.5 min single-core.
#
# Requires: the gmxMMPBSA conda env active (gmx_MMPBSA on PATH), gmx_mpi
# (module load gromacs/2024.3-mixed).
#
# Usage:
#   ./05_run_mmpbsa.sh --run-dir md/runs/gsta1_lig03506 \
#       [--out-dir md/runs/gsta1_lig03506/mmpbsa] \
#       [--start-frame 500] [--end-frame 520] [--interval 1] \
#       [--istrng 0.15] [--fillratio 4.0] [--nranks 1]
set -e

STARTFRAME=300
ENDFRAME=1000
INTERVAL=5
ISTRNG=0.15
FILLRATIO=4.0
NRANKS=1

while [ $# -gt 0 ]; do
    case "$1" in
        --run-dir) RUN_DIR=$2; shift 2 ;;
        --out-dir) OUT_DIR=$2; shift 2 ;;
        --start-frame) STARTFRAME=$2; shift 2 ;;
        --end-frame) ENDFRAME=$2; shift 2 ;;
        --interval) INTERVAL=$2; shift 2 ;;
        --istrng) ISTRNG=$2; shift 2 ;;
        --fillratio) FILLRATIO=$2; shift 2 ;;
        --nranks) NRANKS=$2; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ -z "$RUN_DIR" ]; then
    echo "ERROR: --run-dir is required"
    echo "Usage: $0 --run-dir D [--out-dir D] [--start-frame N] [--end-frame N] [--interval N] [--istrng F] [--fillratio F] [--nranks N]"
    exit 1
fi

for f in "$RUN_DIR/md_0_5.xtc" "$RUN_DIR/md_0_5.tpr" "$RUN_DIR/index.ndx" "$RUN_DIR/topol.top"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: expected input not found: $f"
        exit 1
    fi
done

if ! command -v gmx_MMPBSA >/dev/null 2>&1; then
    echo "ERROR: gmx_MMPBSA not found on PATH -- activate the gmxMMPBSA conda env first"
    echo "  (e.g. conda activate \$MYSCRATCH/conda_envs/gmxMMPBSA)"
    exit 1
fi

OUT_DIR=${OUT_DIR:-"$RUN_DIR/mmpbsa"}

# Resolve to absolute paths before cd'ing.
RUN_DIR=$(cd "$RUN_DIR" && pwd)
case "$OUT_DIR" in
    /*) : ;;
    *) OUT_DIR="$(pwd)/$OUT_DIR" ;;
esac
mkdir -p "$OUT_DIR"

cd "$RUN_DIR"

echo "======================================================================"
echo "STAGE 1: Protein_LIG (protein + GSH, no ligand) index group"
echo "======================================================================"
cp index.ndx "$OUT_DIR/mmpbsa.ndx"
NDX_LISTING=$(printf "q\n" | gmx_mpi make_ndx -f md_0_5.tpr -n "$OUT_DIR/mmpbsa.ndx" -o "$OUT_DIR/mmpbsa.ndx" 2>&1)

get_group_idx() {
    echo "$NDX_LISTING" | grep -E "^\s*[0-9]+\s+$1\s" | head -1 | awk '{print $1}'
}

PROT_IDX=$(get_group_idx "Protein")
LIG_IDX=$(get_group_idx "LIG")
UNL_IDX=$(get_group_idx "UNL")
for v in PROT_IDX LIG_IDX UNL_IDX; do
    if [ -z "${!v}" ]; then
        echo "ERROR: could not find the group behind $v in md_0_5.tpr's index listing:"
        echo "$NDX_LISTING"
        exit 1
    fi
done

MAX_IDX=$(echo "$NDX_LISTING" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)
PROT_LIG_IDX=$(( MAX_IDX + 1 ))
printf "%s | %s\nname %s Protein_LIG\nq\n" "$PROT_IDX" "$LIG_IDX" "$PROT_LIG_IDX" \
    | gmx_mpi make_ndx -f md_0_5.tpr -n "$OUT_DIR/mmpbsa.ndx" -o "$OUT_DIR/mmpbsa.ndx" >/dev/null 2>&1
echo "  Protein=$PROT_IDX  GSH(LIG)=$LIG_IDX  -> Protein_LIG (receptor)=$PROT_LIG_IDX  |  ligand UNL=$UNL_IDX"

echo ""
echo "======================================================================"
echo "STAGE 2: PBC-corrected, water/ion-stripped complex trajectory"
echo "======================================================================"
# Full complex (protein+GSH+ligand) -- -cg in stage 4 splits receptor vs
# ligand, the trajectory itself stays whole.
PROT_LIG_UNL_IDX=$(get_group_idx "Protein_LIG_UNL")
if [ -z "$PROT_LIG_UNL_IDX" ]; then
    echo "ERROR: could not find Protein_LIG_UNL group -- was this run built by the current 03_run_md.sh?"
    exit 1
fi
printf "%s\n%s\n" "$PROT_LIG_UNL_IDX" "$PROT_LIG_UNL_IDX" \
    | gmx_mpi trjconv -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/mmpbsa.ndx" \
        -pbc mol -center -ur compact -o "$OUT_DIR/mmpbsa_traj.xtc"
printf "%s\n%s\n" "$PROT_LIG_UNL_IDX" "$PROT_LIG_UNL_IDX" \
    | gmx_mpi trjconv -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/mmpbsa.ndx" \
        -pbc mol -center -ur compact -dump 0 -o "$OUT_DIR/mmpbsa_ref.pdb"

echo ""
echo "======================================================================"
echo "STAGE 3: stage topology + mmpbsa.in"
echo "======================================================================"
cp topol.top lig.itp lig.prm unl.itp unl.prm posre*.itp "$OUT_DIR/" 2>/dev/null
FFDIR_NAME=$(basename "$(find . -maxdepth 1 -name '*.ff' -print -quit)")
if [ -z "$FFDIR_NAME" ]; then
    echo "ERROR: no *.ff force field directory found in $RUN_DIR"
    exit 1
fi
ln -sfn "$RUN_DIR/$FFDIR_NAME" "$OUT_DIR/$FFDIR_NAME"
echo "[staged] $(ls "$OUT_DIR"/*.itp "$OUT_DIR"/*.prm "$OUT_DIR"/topol.top 2>/dev/null | wc -l) topology files + $FFDIR_NAME symlink"

cat > "$OUT_DIR/mmpbsa.in" <<EOF
Rescoring input for GST + GSH + docked ligand (Single Trajectory, CHARMM36/CGenFF).
Receptor = Protein + GSH, Ligand = docked candidate. PB recommended over GB
for CHARMM-prepped systems per gmx_MMPBSA's own CHARMMff documentation.

&general
sys_name="GST_rescoring",
startframe=$STARTFRAME,
endframe=$ENDFRAME,
interval=$INTERVAL,
PBRadii=7,
/
&pb
istrng=$ISTRNG,
fillratio=$FILLRATIO,
radiopt=0,
/
EOF
echo "[mmpbsa.in] frames $STARTFRAME-$ENDFRAME (interval $INTERVAL), istrng=$ISTRNG, fillratio=$FILLRATIO, PBRadii=7 (charmm_radii), radiopt=0"

echo ""
echo "======================================================================"
echo "STAGE 4: gmx_MMPBSA"
echo "======================================================================"
cd "$OUT_DIR"
# -cs needs the absolute path -- md_0_5.tpr lives in $RUN_DIR, not $OUT_DIR.
#
# --nranks: parallelises over frames via MPI (each PB grid solve is
# embarrassingly parallel). Defaults to 1 (serial, the validated config).
# On Setonix, real MPI parallelism is blocked -- see docs/implementation_notes.md
# and slurm/run_mmpbsa_array.sh for the working parallel path.
GMX_MMPBSA_CMD="gmx_MMPBSA"
if [ "$NRANKS" -gt 1 ]; then
    GMX_MMPBSA_CMD="mpirun -np $NRANKS gmx_MMPBSA"
fi
$GMX_MMPBSA_CMD -O -i mmpbsa.in -cs "$RUN_DIR/md_0_5.tpr" -ci mmpbsa.ndx -cg "$PROT_LIG_IDX" "$UNL_IDX" \
    -ct mmpbsa_traj.xtc -cp topol.top -o FINAL_RESULTS_MMPBSA.dat -eo FINAL_RESULTS_MMPBSA.csv -nogui

echo ""
echo "[OK] MM-PBSA complete: $OUT_DIR/FINAL_RESULTS_MMPBSA.dat / .csv"
echo "Once this small validation slice runs clean, widen --start-frame/--end-frame/--interval"
echo "for a statistically meaningful estimate (50-200+ frames is typical)."
