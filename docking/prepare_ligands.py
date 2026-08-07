#!/usr/bin/env python3
# =============================================================================
# Ligand preparation: SMILES -> 3D conformer (RDKit) -> docking-ready PDBQT (meeko)
# =============================================================================
#
# Reads the same training_data.csv used by the ML pipeline and produces the
# PDBQT ligand files that Vina_rigid_A.pl / Vina_rigid_M.pl / Vina_rigid_P.pl
# expect, plus a ligand.txt manifest in the format those scripts read.
#
# This removes the need to prepare ligands by hand: everything from SMILES
# to a docking-ready PDBQT happens on the command line via RDKit + meeko,
# no external webserver or GUI tool involved.
#
# Outputs (written to the current working directory):
#   ligand_pdbqt/lig_<ID>.pdbqt   one file per successfully prepared ligand
#   ligand.txt                     manifest of basenames, for Vina_rigid_*.pl
#   ligand_prep_report.csv         per-ligand status (ok / failed + reason)
#   ligand_prep.log                run log
#
# Usage (run from inside docking/, sibling to ml_pipeline/):
#   python3 prepare_ligands.py --outdir ligand_pdbqt
# (defaults to reading ../ml_pipeline/training_data.csv; override with --input
# if you're running from somewhere else)
#
# Requires: rdkit, meeko, pandas, scipy, gemmi
#   pip install rdkit meeko pandas scipy gemmi
#   (rdkit and pandas are already in requirements.txt; meeko/scipy/gemmi are new)
#
# Tested on the actual training_data.csv (150-ligand and 900-ligand subsets,
# 0 failures) before being added to this pipeline.
# =============================================================================

import argparse
import csv
import logging
import multiprocessing as mp
import os
import sys
import time

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")  # silence RDKit's warning spam on messy SMILES

# Same seed used throughout the ML pipeline (pipeline_HAL20260712_ensemble.py)
# so conformer generation here is reproducible and consistent with it.
GLOBAL_SEED = 42


def prepare_one(job):
    """Embed one SMILES in 3D, optimise, and write a PDBQT via meeko.

    Returns a (ligand_id, smiles, status, detail) tuple. `detail` is the
    output filename on success, or a short failure reason on failure.
    Never raises -- failures are caught and reported so one bad SMILES
    can't take down a batch run.
    """
    lig_id, smiles, outdir = job
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return (lig_id, smiles, "failed", "invalid SMILES")

        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = GLOBAL_SEED
        cid = AllChem.EmbedMolecule(mol, params)
        if cid < 0:
            # retry once with random coordinates for awkward/strained molecules
            params.useRandomCoords = True
            cid = AllChem.EmbedMolecule(mol, params)
        if cid < 0:
            return (lig_id, smiles, "failed", "3D embedding failed")

        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
        except Exception:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=2000)
            except Exception:
                pass  # keep the embedded (unoptimised) geometry rather than dropping the ligand

        mk_prep = MoleculePreparation()  # default charge_model="gasteiger"
        mol_setups = mk_prep.prepare(mol)
        if not mol_setups:
            return (lig_id, smiles, "failed", "meeko produced no setup")

        pdbqt_string, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
        if not is_ok:
            return (lig_id, smiles, "failed", f"meeko write failed: {err}")

        fname = f"lig_{int(lig_id):05d}.pdbqt"
        with open(os.path.join(outdir, fname), "w") as fh:
            fh.write(pdbqt_string)
        return (lig_id, smiles, "ok", fname)

    except Exception as e:
        return (lig_id, smiles, "failed", f"exception: {e}")


def load_ligands(filepath: str) -> pd.DataFrame:
    df = pd.read_csv(filepath)
    col_map = {}
    for canon, cands in [
        ("ID", ["ID", "id", "ligand_id", "Ligand_ID"]),
        ("SMILES", ["SMILES", "Smiles", "smiles", "SMILE"]),
    ]:
        col = next((c for c in cands if c in df.columns), None)
        if col:
            col_map[canon] = col
    missing = [c for c in ("ID", "SMILES") if c not in col_map]
    if missing:
        raise ValueError(f"Could not find columns {missing} in {filepath}")
    return df.rename(columns={v: k for k, v in col_map.items()})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="../ml_pipeline/training_data.csv",
                     help="CSV with ID + SMILES columns (default: ../ml_pipeline/training_data.csv)")
    ap.add_argument("--outdir", default="ligand_pdbqt", help="Directory for output PDBQT files")
    ap.add_argument("--nprocs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, default=None, help="Only process the first N rows (for testing)")
    ap.add_argument("--ids", type=str, default=None,
                     help="Comma-separated list of specific ligand IDs to include, e.g. for "
                          "topping up a --limit test batch with IDs you already have CGenFF "
                          ".str files for. Combined (union) with --limit if both are given.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler("ligand_prep.log"), logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger(__name__)

    df = load_ligands(args.input)
    selected = None
    if args.ids:
        wanted = {int(x) for x in args.ids.split(",") if x.strip()}
        selected = df[df["ID"].astype(int).isin(wanted)]
    if args.limit:
        head_df = df.head(args.limit)
        selected = head_df if selected is None else pd.concat([selected, head_df]).drop_duplicates(subset="ID")
    if selected is not None:
        df = selected.sort_values("ID")

    jobs = [(row["ID"], row["SMILES"], args.outdir) for _, row in df.iterrows()]
    log.info("Preparing %d ligands with %d processes -> %s/", len(jobs), args.nprocs, args.outdir)

    def checkpoint(results):
        """Overwrite the report + manifest with progress so far. Cheap enough
        to call periodically so a walltime cutoff or preemption on the
        cluster doesn't lose completed work."""
        with open("ligand_prep_report.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Ligand_ID", "SMILES", "status", "detail"])
            for lig_id, smiles, status, detail in sorted(results, key=lambda r: int(r[0])):
                w.writerow([lig_id, smiles, status, detail])
        ok_files = sorted(detail for _, _, status, detail in results if status == "ok")
        with open("ligand.txt", "w") as fh:
            for fname in ok_files:
                fh.write(fname + "\n")
        return len(ok_files)

    t0 = time.time()
    results = []
    with mp.Pool(args.nprocs) as pool:
        for i, res in enumerate(pool.imap_unordered(prepare_one, jobs, chunksize=8), 1):
            results.append(res)
            if i % 200 == 0 or i == len(jobs):
                log.info("  %d/%d processed", i, len(jobs))
                checkpoint(results)

    n_ok = sum(1 for r in results if r[2] == "ok")
    n_fail = len(results) - n_ok
    elapsed = time.time() - t0
    log.info("Done in %.1fs (%.2fs/ligand): %d ok, %d failed",
              elapsed, elapsed / max(1, len(jobs)), n_ok, n_fail)

    n_written = checkpoint(results)
    log.info("Wrote ligand_prep_report.csv and ligand.txt (%d entries)", n_written)

    if n_fail:
        log.warning("%d/%d ligands failed prep -- see ligand_prep_report.csv for reasons",
                     n_fail, len(jobs))


if __name__ == "__main__":
    main()
