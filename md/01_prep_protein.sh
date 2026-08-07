#!/bin/bash --login
# =============================================================================
# 01_prep_protein.sh -- one-time protein preparation
# =============================================================================
# Wraps Harry's documented step:
#   gmx_mpi pdb2gmx -f protein.pdb -o protein_processed.gro -ignh
#   (manual: SPCE water, N-terminus NH3+, C-terminus COO-)
#
# Run this ONCE per protein structure (GSTA1-1 / GSTM1-1 / GSTP1-1) -- the
# output is reused across every ligand run against that receptor.
#
# WHY THIS DOESN'T USE -ter (changed after hitting a real hang on Setonix):
# -ter makes pdb2gmx prompt interactively, once per PROTEIN CHAIN, for the
# N- and C-terminus type. On Setonix's gmx_mpi build, piped stdin answers
# worked for chain 1 but the terminal sat waiting for real keyboard input
# on chain 2 -- piped multi-prompt stdin isn't reaching the process
# reliably across chains. Rather than fight that, this drops -ter
# entirely: GROMACS's documented DEFAULT (no -ter flag at all) already
# applies charged termini automatically and non-interactively -- NH3+ on
# the N-terminus, COO- on the C-terminus -- for every chain, which is
# exactly what Harry's protocol calls for. No prompts, no hang.
#
# --nterm/--cterm are kept below only for the rare case you deliberately
# want NON-default termini; supplying either re-enables -ter and pipes
# your answer, so be aware that reintroduces the per-chain prompt risk
# above. For the normal charged-termini case (the common one), leave
# them unset.
#
# After running, still worth opening topol.top once to confirm both
# chains got NH3+ / COO- as expected.
#
# IMPORTANT: the input PDB must NOT contain GSH (or any other cofactor/
# ligand HETATM) -- GSH gets its own separate mol2 -> CGenFF -> itp
# pipeline via 02_prep_molecule.sh and is merged back in later by
# 03_run_md.sh. Leaving GSH in here makes pdb2gmx try to treat it as
# another chain needing its own terminus assignment (GSH is chemically a
# tripeptide, so pdb2gmx doesn't reject it outright -- it just gets
# confused). Strip all non-protein HETATM records (GSH included) before
# running this script.
#
# --allow-missing passes -missing through to pdb2gmx, letting it geometry-
# idealize side-chain atoms the crystal structure didn't resolve (check the
# PDB's REMARK 470 section first -- if it lists missing atoms, you need this
# flag or pdb2gmx errors out instead of silently guessing anything).
#
# Usage:
#   ./01_prep_protein.sh protein.pdb outdir/ [ffdir] [--nterm N] [--cterm N] [--allow-missing]
# =============================================================================
set -e

PROTEIN_PDB=$1
OUTDIR=$2
FFDIR=${3:-charmm36-jul2022.ff}
NTERM=""
CTERM=""
ALLOW_MISSING=0
shift 3 2>/dev/null || true
while [ $# -gt 0 ]; do
    case "$1" in
        --nterm) NTERM=$2; shift 2 ;;
        --cterm) CTERM=$2; shift 2 ;;
        --allow-missing) ALLOW_MISSING=1; shift 1 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done
MISSING_FLAG=""
if [ "$ALLOW_MISSING" = "1" ]; then
    MISSING_FLAG="-missing"
fi
# -ter only gets turned on if the caller explicitly asked for non-default
# termini via --nterm/--cterm (both must be given together).
USE_TER=0
if [ -n "$NTERM" ] || [ -n "$CTERM" ]; then
    if [ -z "$NTERM" ] || [ -z "$CTERM" ]; then
        echo "ERROR: --nterm and --cterm must both be given together"
        exit 1
    fi
    USE_TER=1
fi

if [ -z "$PROTEIN_PDB" ] || [ -z "$OUTDIR" ]; then
    echo "Usage: $0 protein.pdb outdir/ [ffdir] [--nterm N] [--cterm N]"
    exit 1
fi
if [ ! -f "$PROTEIN_PDB" ]; then
    echo "ERROR: $PROTEIN_PDB not found"
    exit 1
fi

mkdir -p "$OUTDIR"
cp "$PROTEIN_PDB" "$OUTDIR/protein_input.pdb"
# pdb2gmx looks for the force field directory locally (or via $GMXLIB) --
# symlink it in so -ff resolves regardless of where FFDIR actually lives.
FFDIR_ABS=$(cd "$(dirname "$FFDIR")" && pwd)/$(basename "$FFDIR")
ln -sfn "$FFDIR_ABS" "$OUTDIR/$(basename "$FFDIR")"
cd "$OUTDIR"

if [ "$USE_TER" = "1" ]; then
    echo "[1/1] gmx_mpi pdb2gmx (water=spce, ignoring existing H, -ter with piped answers)"
    echo "NOTE: piped multi-chain -ter answers have hung before on Setonix's"
    echo "gmx_mpi build past chain 1. If this hangs, Ctrl+C and re-run without"
    echo "--nterm/--cterm to use default charged termini instead."
    TERM_ANSWERS=$(for _ in $(seq 1 20); do printf "%s\n%s\n" "$NTERM" "$CTERM"; done)
    printf "%s" "$TERM_ANSWERS" | gmx_mpi pdb2gmx \
        -f protein_input.pdb \
        -o protein_processed.gro \
        -p topol.top \
        -water spce \
        -ff "$(basename "$FFDIR" .ff)" \
        -ignh -ter $MISSING_FLAG
else
    echo "[1/1] gmx_mpi pdb2gmx (water=spce, ignoring existing H, default charged termini)"
    gmx_mpi pdb2gmx \
        -f protein_input.pdb \
        -o protein_processed.gro \
        -p topol.top \
        -water spce \
        -ff "$(basename "$FFDIR" .ff)" \
        -ignh $MISSING_FLAG
fi

echo ""
echo "[OK] protein_processed.gro + topol.top written to $OUTDIR"
echo ""
echo "VERIFY: open topol.top and confirm the first and last residue entries"
echo "carry NH3+ / COO- (GROMACS's default charged termini) on every chain."
