#!/bin/bash --login
# One-time protein prep, once per isoform: pdb2gmx (CHARMM36, SPC/E water,
# default charged termini, no -ter). Input PDB must have GSH/HETATM already
# stripped. Rationale + gotchas: docs/implementation_notes.md.
#
# Usage:
#   ./01_prep_protein.sh protein.pdb outdir/ [ffdir] [--nterm N] [--cterm N] [--allow-missing]
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
# -ter only turns on if the caller explicitly asked for non-default termini.
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
# Symlink the FF dir in so -ff resolves regardless of where FFDIR lives.
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
