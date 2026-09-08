#!/bin/bash --login
# Post-processing for a completed 03_run_md.sh production run: QC energies,
# RMSD/RMSF, warhead-to-GSH-thiolate distance, plus a VMD visualization
# trajectory and automatic snapshot + movie rendering (VMD headless +
# ffmpeg, no manual follow-up unless --skip-render is given).
#
# Produces (under --out-dir, default <run-dir>/analysis):
#   energy_{potential,temperature,pressure,density}.xvg/.png
#   rmsd_protein_backbone / rmsd_ligand_fit_on_protein / rmsd_gsh_fit_on_protein
#   rmsf_protein_calpha
#   warhead_sg2_distance -- ligand warhead to GSH thiolate (SG2) over time
#   viz_protein_lig_gsh.pdb/.xtc, vmd_load_trajectory.tcl,
#   vmd_render_snapshot_and_frames.tcl, snapshot + frames + md_movie.mp4
#
# Warhead-atom lookup and the VMD/ffmpeg rendering approach are explained in
# docs/implementation_notes.md.
#
# Requires: gmx_mpi (module load gromacs/2024.3-mixed), python3 with rdkit +
# matplotlib (gst_ml env), vmd 1.9.3 (spack load vmd@1.9.3, or --vmd-bin),
# ffmpeg on PATH. --skip-render stops after the .pdb/.xtc/.tcl files.
#
# Usage:
#   ./04_analyze_md.sh --run-dir md/runs/gsta1_lig03506 \
#       --ligand-mol2 md/system/gsta1_lig03506_prepped/3506_A_fixed.mol2 \
#       [--out-dir md/runs/gsta1_lig03506/analysis] [--stride 20] \
#       [--sg-resname LIG] [--sg-atomname SG2] [--vmd-bin /path/to/vmd] \
#       [--skip-render]
set -e

TOOLS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/tools"

STRIDE=20
SG_RESNAME=LIG
SG_ATOMNAME=SG2
VMD_BIN=""
SKIP_RENDER=0
# Known-good spack install (confirmed on Setonix, 2026-08-14); fallback only
# if `vmd` isn't on PATH and --vmd-bin wasn't given.
SPACK_VMD="/software/projects/pawsey1376/volet/setonix/2025.08/software/linux-sles15-zen3/gcc-14.2.0/vmd-1.9.3-cokkao5hjqplpzh5zjwr7h4irrcb647h/bin/vmd"

while [ $# -gt 0 ]; do
    case "$1" in
        --run-dir) RUN_DIR=$2; shift 2 ;;
        --ligand-mol2) LIGAND_MOL2=$2; shift 2 ;;
        --out-dir) OUT_DIR=$2; shift 2 ;;
        --stride) STRIDE=$2; shift 2 ;;
        --sg-resname) SG_RESNAME=$2; shift 2 ;;
        --sg-atomname) SG_ATOMNAME=$2; shift 2 ;;
        --vmd-bin) VMD_BIN=$2; shift 2 ;;
        --skip-render) SKIP_RENDER=1; shift 1 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

for v in RUN_DIR LIGAND_MOL2; do
    if [ -z "${!v}" ]; then
        flag=$(echo "$v" | tr '[:upper:]_' '[:lower:]-')
        echo "ERROR: --${flag} is required"
        echo "Usage: $0 --run-dir D --ligand-mol2 F [--out-dir D] [--stride N] [--sg-resname R] [--sg-atomname A]"
        exit 1
    fi
done

for f in "$RUN_DIR/md_0_5.xtc" "$RUN_DIR/md_0_5.tpr" "$RUN_DIR/md_0_5.edr" \
         "$RUN_DIR/index.ndx" "$RUN_DIR/topol.top" "$LIGAND_MOL2"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: expected input not found: $f"
        exit 1
    fi
done

OUT_DIR=${OUT_DIR:-"$RUN_DIR/analysis"}

# Resolve everything to absolute before cd'ing into RUN_DIR.
RUN_DIR=$(cd "$RUN_DIR" && pwd)
LIGAND_MOL2=$(cd "$(dirname "$LIGAND_MOL2")" && pwd)/$(basename "$LIGAND_MOL2")
case "$OUT_DIR" in
    /*) : ;;
    *) OUT_DIR="$(pwd)/$OUT_DIR" ;;
esac
mkdir -p "$OUT_DIR"

cd "$RUN_DIR"

echo "======================================================================"
echo "STAGE 1: QC -- energy stability"
echo "======================================================================"
for term in Potential Temperature Pressure Density; do
    echo "[energy] $term"
    printf "%s\n0\n" "$term" | gmx_mpi energy -f md_0_5.edr -o "$OUT_DIR/energy_${term,,}.xvg"
done

echo ""
echo "======================================================================"
echo "STAGE 2: build analysis index groups"
echo "======================================================================"
echo "[warhead] identifying the ligand's electrophilic atom from $LIGAND_MOL2"
WARHEAD_LOCAL_NUM=$(python3 "$TOOLS_DIR/find_warhead_atom.py" "$LIGAND_MOL2")
echo "  ligand-local atom number: $WARHEAD_LOCAL_NUM (see stderr above for element/name)"

cp index.ndx "$OUT_DIR/analysis.ndx"
NDX_LISTING0=$(printf "q\n" | gmx_mpi make_ndx -f md_0_5.tpr -n "$OUT_DIR/analysis.ndx" -o "$OUT_DIR/analysis.ndx" 2>&1)

get_group_idx() {
    # $1 = group name to match (word-bounded); echoes the group number
    echo "$NDX_LISTING0" | grep -E "^\s*[0-9]+\s+$1\s" | head -1 | awk '{print $1}'
}
get_group_natoms() {
    # $1 = group name to match; echoes the atom count reported on that line
    echo "$NDX_LISTING0" | grep -E "^\s*[0-9]+\s+$1\s" | head -1 | grep -oE '[0-9]+ atoms' | grep -oE '^[0-9]+'
}

LIG_IDX=$(get_group_idx "LIG")
UNL_IDX=$(get_group_idx "UNL")
BACKBONE_IDX=$(get_group_idx "Backbone")
CALPHA_IDX=$(get_group_idx "C-alpha")
PROT_LIG_UNL_IDX=$(get_group_idx "Protein_LIG_UNL")
PROTEIN_NATOMS=$(get_group_natoms "Protein")
LIG_NATOMS=$(get_group_natoms "LIG")
for v in LIG_IDX UNL_IDX BACKBONE_IDX CALPHA_IDX PROT_LIG_UNL_IDX PROTEIN_NATOMS LIG_NATOMS; do
    if [ -z "${!v}" ]; then
        echo "ERROR: could not find/parse the group behind $v in md_0_5.tpr's index listing:"
        echo "$NDX_LISTING0"
        exit 1
    fi
done

# Whole-system atom number for the warhead = protein + GSH atom counts
# (complex.gro's fixed ordering) + its local position. See implementation_notes.md.
WARHEAD_GLOBAL_NUM=$(( PROTEIN_NATOMS + LIG_NATOMS + WARHEAD_LOCAL_NUM ))
echo "  sanity check -- unl.gro line for local atom $WARHEAD_LOCAL_NUM (should show a real atom, not blank):"
sed -n "$((WARHEAD_LOCAL_NUM + 2))p" unl.gro | sed 's/^/    /'
echo "  -> whole-system atom number: $PROTEIN_NATOMS (protein) + $LIG_NATOMS (GSH) + $WARHEAD_LOCAL_NUM (local) = $WARHEAD_GLOBAL_NUM"

add_ndx_group() {
    # $1 = selection expression piped to make_ndx; sets $NEW_GROUP_IDX to
    # the resulting (always-highest-numbered) new group.
    local listing
    listing=$(printf "%s\nq\n" "$1" | gmx_mpi make_ndx -f md_0_5.tpr -n "$OUT_DIR/analysis.ndx" -o "$OUT_DIR/analysis.ndx" 2>&1)
    NEW_GROUP_IDX=$(echo "$listing" | grep -E "^\s*[0-9]+\s+\S+\s*:" | awk '{print $1}' | sort -n | tail -1)
    if [ -z "$NEW_GROUP_IDX" ]; then
        echo "ERROR: make_ndx group creation failed for expression: $1"
        echo "$listing"
        exit 1
    fi
}

add_ndx_group "$LIG_IDX & ! a H*";                    LIG_NOH_IDX=$NEW_GROUP_IDX
add_ndx_group "$UNL_IDX & ! a H*";                    UNL_NOH_IDX=$NEW_GROUP_IDX
add_ndx_group "$LIG_IDX & a $SG_ATOMNAME";            SG2_IDX=$NEW_GROUP_IDX
add_ndx_group "a $WARHEAD_GLOBAL_NUM";                WARHEAD_IDX=$NEW_GROUP_IDX
echo "  LIG(noH)=$LIG_NOH_IDX  UNL(noH)=$UNL_NOH_IDX  SG2=$SG2_IDX  warhead=$WARHEAD_IDX  Backbone=$BACKBONE_IDX  C-alpha=$CALPHA_IDX"

echo ""
echo "======================================================================"
echo "STAGE 3: RMSD / RMSF"
echo "======================================================================"
echo "[rmsd] protein backbone (self-fit)"
printf "%s\n%s\n" "$BACKBONE_IDX" "$BACKBONE_IDX" \
    | gmx_mpi rms -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -tu ns -o "$OUT_DIR/rmsd_protein_backbone.xvg"

echo "[rmsd] ligand (UNL), fit on protein backbone"
printf "%s\n%s\n" "$BACKBONE_IDX" "$UNL_NOH_IDX" \
    | gmx_mpi rms -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -tu ns -o "$OUT_DIR/rmsd_ligand_fit_on_protein.xvg"

echo "[rmsd] GSH (LIG), fit on protein backbone"
printf "%s\n%s\n" "$BACKBONE_IDX" "$LIG_NOH_IDX" \
    | gmx_mpi rms -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -tu ns -o "$OUT_DIR/rmsd_gsh_fit_on_protein.xvg"

echo "[rmsf] protein C-alpha, per-residue"
printf "%s\n" "$CALPHA_IDX" \
    | gmx_mpi rmsf -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -res -o "$OUT_DIR/rmsf_protein_calpha.xvg"

echo ""
echo "======================================================================"
echo "STAGE 4: warhead -> GSH thiolate (SG2) distance over time"
echo "======================================================================"
printf "%s\n" "$WARHEAD_IDX" \
    | gmx_mpi traj -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -com -tu ns -ox "$OUT_DIR/_warhead_pos.xvg"
printf "%s\n" "$SG2_IDX" \
    | gmx_mpi traj -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" -com -tu ns -ox "$OUT_DIR/_sg2_pos.xvg"
python3 "$TOOLS_DIR/compute_distance_xvg.py" "$OUT_DIR/_warhead_pos.xvg" "$OUT_DIR/_sg2_pos.xvg" "$OUT_DIR/warhead_sg2_distance.xvg"
rm -f "$OUT_DIR/_warhead_pos.xvg" "$OUT_DIR/_sg2_pos.xvg"

echo ""
echo "======================================================================"
echo "STAGE 5: plots"
echo "======================================================================"
for f in energy_potential energy_temperature energy_pressure energy_density \
         rmsd_protein_backbone rmsd_ligand_fit_on_protein rmsd_gsh_fit_on_protein \
         warhead_sg2_distance; do
    [ -f "$OUT_DIR/$f.xvg" ] && python3 "$TOOLS_DIR/plot_xvg.py" "$OUT_DIR/$f.xvg" "$OUT_DIR/$f.png"
done
python3 "$TOOLS_DIR/plot_xvg.py" "$OUT_DIR/rmsf_protein_calpha.xvg" "$OUT_DIR/rmsf_protein_calpha.png" \
    --xlabel "Residue" --ylabel "RMSF (nm)" --title "Per-residue RMSF (protein C-alpha)"

echo ""
echo "======================================================================"
echo "STAGE 6: VMD-ready visualization trajectory"
echo "======================================================================"
echo "[trjconv] PBC-corrected, water/ion-stripped, stride=$STRIDE"
printf "%s\n%s\n" "$PROT_LIG_UNL_IDX" "$PROT_LIG_UNL_IDX" \
    | gmx_mpi trjconv -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" \
        -pbc mol -center -ur compact -skip "$STRIDE" -o "$OUT_DIR/viz_protein_lig_gsh.xtc"
printf "%s\n%s\n" "$PROT_LIG_UNL_IDX" "$PROT_LIG_UNL_IDX" \
    | gmx_mpi trjconv -s md_0_5.tpr -f md_0_5.xtc -n "$OUT_DIR/analysis.ndx" \
        -pbc mol -center -ur compact -dump 0 -o "$OUT_DIR/viz_protein_lig_gsh.pdb"

echo "[vmd script] writing $OUT_DIR/vmd_load_trajectory.tcl"
cat > "$OUT_DIR/vmd_load_trajectory.tcl" <<'TCL_EOF'
# Loads structure + trajectory + representations only -- no rendering, no
# `quit`. Use this to watch the animation live in VMD's own GUI:
#   vmd -e vmd_load_trajectory.tcl
# vmd_render_snapshot_and_frames.tcl sources this same file for the headless
# path, so edit here if you want the view to look different.

mol new viz_protein_lig_gsh.pdb
mol addfile viz_protein_lig_gsh.xtc type xtc waitfor all

mol delrep 0 top
mol representation NewCartoon
mol color Structure
mol selection "protein"
mol addrep top

mol representation Licorice 0.3 12 12
mol color Name
mol selection "resname LIG"
mol addrep top

mol representation Licorice 0.3 12 12
mol color Name
mol selection "resname UNL"
mol addrep top

display projection Orthographic
display depthcue off
color Display Background white
axes location Off
TCL_EOF

echo "[vmd script] writing $OUT_DIR/vmd_render_snapshot_and_frames.tcl"
cat > "$OUT_DIR/vmd_render_snapshot_and_frames.tcl" <<'TCL_EOF'
# Headless rendering: one hero snapshot (final frame) + a numbered frame
# sequence, via TachyonInternal. Sources vmd_load_trajectory.tcl for the
# structure/representation setup. Run automatically by 04_analyze_md.sh;
# only run by hand to re-render without repeating the whole analysis:
#   vmd -dispdev text -e vmd_render_snapshot_and_frames.tcl

source vmd_load_trajectory.tcl

animate goto end
render TachyonInternal snapshot_final_frame.tga

set nf [molinfo top get numframes]
for {set i 0} {$i < $nf} {incr i} {
    animate goto $i
    render TachyonInternal [format "frame_%04d.tga" $i]
}
puts "Rendered $nf frames + snapshot_final_frame.tga"
quit
TCL_EOF

if [ "$SKIP_RENDER" = "1" ]; then
    echo ""
    echo "[OK] Analysis complete: $OUT_DIR (--skip-render set, VMD/ffmpeg not run)"
    echo "Render later with: vmd -dispdev text -e vmd_render_snapshot_and_frames.tcl"
    echo "                   then: ffmpeg -framerate 10 -i frame_%04d.tga -pix_fmt yuv420p md_movie.mp4"
    exit 0
fi

echo ""
echo "======================================================================"
echo "STAGE 7: render snapshot + movie (VMD, headless, then ffmpeg)"
echo "======================================================================"
if [ -n "$VMD_BIN" ]; then
    :
elif command -v vmd >/dev/null 2>&1; then
    VMD_BIN=vmd
elif [ -x "$SPACK_VMD" ]; then
    VMD_BIN="$SPACK_VMD"
else
    echo "ERROR: vmd not found. Run 'spack load vmd@1.9.3' first, add vmd to PATH,"
    echo "pass --vmd-bin /path/to/vmd, or re-run with --skip-render and do this part later."
    exit 1
fi
if ! command -v ffmpeg >/dev/null 2>&1; then
    echo "ERROR: ffmpeg not found on PATH -- needed to stitch rendered frames into md_movie.mp4."
    echo "Re-run with --skip-render and do this part later once ffmpeg is available."
    exit 1
fi

echo "[vmd] $VMD_BIN -dispdev text -e vmd_render_snapshot_and_frames.tcl"
( cd "$OUT_DIR" && "$VMD_BIN" -dispdev text -e vmd_render_snapshot_and_frames.tcl )

echo "[ffmpeg] stitching frames -> md_movie.mp4"
( cd "$OUT_DIR" && ffmpeg -y -framerate 10 -i frame_%04d.tga -pix_fmt yuv420p md_movie.mp4 )

echo ""
echo "[OK] Analysis complete: $OUT_DIR"
echo "  snapshot_final_frame.tga, frame_0000.tga.., md_movie.mp4 -- all produced by this run."
echo "  To watch it live instead: vmd -e vmd_load_trajectory.tcl (from inside $OUT_DIR)"
