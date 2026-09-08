# GST Metabolic Stability Pipeline

This project predicts and investigates the metabolic stability of covalent drug
candidates against glutathione S-transferase (GST) — the Phase II detox enzyme
that conjugates glutathione (GSH) onto electrophilic ("warhead"-bearing)
xenobiotics via its cysteine thiolate. A drug candidate that GST clears too
quickly has poor metabolic stability (short half-life, therapeutic failure),
so the goal across both halves of this project is to flag which candidates
are likely to be reactive GST substrates.

There are two largely independent approaches living side by side in this
repo, which is why it's split into `ml_pipeline/` and a structure-based
`docking/` + `md/` pair:

- **`ml_pipeline/`** — a ligand-only machine-learning classifier. Fast,
  works directly from SMILES, no 3D structure or HPC docking/MD required.
  Already trained and validated on 3604 ligands.
- **`docking/` + `md/`** — a structure-based pipeline: dock each candidate
  against the real GST isoform structure with AutoDock Vina, then run full
  molecular dynamics (GROMACS) on the best candidates to see how the
  ligand actually behaves in the binding site over time. Slower, mechanistic,
  complements the ML pipeline rather than replacing it.

The two aren't fully siloed — `docking/check_warhead_proximity.py` reuses the
same SMARTS-based electrophile ("warhead") detection logic that the ML
pipeline uses for its `alert_score` feature, so both approaches agree on what
counts as a reactive substructure.

A visual architecture diagram of the docking + MD side is at
[`docs/pipeline_workflow_diagram.svg`](docs/pipeline_workflow_diagram.svg).
Script comments are kept short on purpose; the engineering detail behind
them (bug fixes, environment quirks, modelling choices) lives in
[`docs/implementation_notes.md`](docs/implementation_notes.md).

## Repository layout

```
pipeline_development/
├── ml_pipeline/            # ML stability classifier (SMILES -> stable/unstable)
├── docking/                 # AutoDock Vina docking (structure-based, per isoform)
│   ├── receptors/            # per-isoform receptor PDBQT + grid configs
│   ├── configs/               #   (raw-PDB-coordinate frame, one grid box per isoform)
│   ├── scripts/               # Vina_rigid_{A,M,P}.pl -- per-isoform Perl docking wrappers
│   ├── slurm/                 # SLURM array job scripts (Setonix)
│   ├── prepare_ligands.py     # SMILES -> 3D -> PDBQT (meeko)
│   ├── check_warhead_proximity.py  # candidate selection: score + warhead-distance ranking
│   ├── ligand_pdbqt/, results/, logs/  # generated, not tracked in git
├── md/                       # GROMACS molecular dynamics (structure-based, per candidate)
│   ├── 01_prep_protein.sh     # protein-only PDB -> GROMACS topology (pdb2gmx, CHARMM36)
│   ├── 02_prep_molecule.sh    # GSH or ligand mol2 -> CGenFF .str -> GROMACS itp/prm
│   ├── 03_run_md.sh           # assemble system -> solvate -> EM -> NVT -> NPT -> production MD
│   ├── 04_analyze_md.sh       # post-MD QC, RMSD/RMSF, warhead-distance, VMD snapshot + MP4
│   ├── 05_run_mmpbsa.sh       # MM-PBSA rescoring (gmx_MMPBSA), single-trajectory ST protocol
│   ├── mdp/                   # ions.mdp / minim.mdp / nvt.mdp / npt.mdp / md.mdp
│   ├── tools/                 # helper + third-party scripts 03/04/05 depend on
│   ├── slurm/                 # batch wrappers: run_md.sh, run_mmpbsa.sh, run_mmpbsa_array.sh + submit_mmpbsa_array.sh
│   ├── system/                 # per-isoform prepped protein/GSH directories (generated)
│   ├── ligand_mol2_cache/, ligand_str_cache/  # cached CGenFF inputs/outputs per isoform
├── docs/                     # reference diagrams (this file's companion)
└── README.md                 # this file
```

---

## Part 1 — ML pipeline (`ml_pipeline/pipeline_HAL20260712_ensemble.py`)

**What it predicts:** whether a covalent ligand is metabolically *unstable*
(i.e. reactive) to GST-mediated GSH conjugation, directly from its SMILES —
no docking or MD needed.

**Model:** a soft-voting ensemble of three classifiers (Random Forest +
SVM + Gradient Boosting), trained with SMOTE class balancing on 8 features:

- 5 RDKit molecular descriptors (`SMR_VSA7`, `SMR_VSA9`, `PEOE_VSA9`,
  `SLogP_VSA1`, `CSP3`)
- 3 mechanism-based electrophilicity descriptors: `omega_global` (global
  electrophilicity index), `omega_local_warhead` (electrophilicity localised
  to the detected warhead atom), and `alert_score` (weighted count of ~50
  structural-alert SMARTS patterns matched — SNAr leaving groups, Michael
  acceptors, epoxides, alkyl halides, etc.)

**Pipeline stages (`main()`):**

1. **`stage1_build_features`** — compute RDKit + electrophilicity descriptors
   for every ligand in `training_data.csv` (3604 ligands: `smiles, ID,
   GROUND_TRUTH`). Electrophilicity descriptors use xtb if available on
   `PATH`, otherwise fall back to an RDKit-based surrogate.
2. **`stage2_train_ml_model`** — 80/20 train/test split (721 held-out test
   ligands), 10-fold CV, trains the ensemble.
3. **`stage3_benchmark_and_report`** — evaluates on the held-out test set,
   runs validation controls (bootstrap 95% CIs, Williams-leverage
   applicability domain, y-randomisation label-shuffling control, SHAP,
   permutation importance), and writes everything to `pipeline_outputs/`.

**Validated results (`pipeline_outputs/performance_report.txt`, most recent
run):**

| | Train (n=2882) | Test (n=721) |
|---|---|---|
| MCC | 0.735 | 0.753 |
| ROC-AUC | 0.985 | 0.974 |
| Precision | 0.806 | 0.857 |
| Recall | 0.676 | 0.667 |
| Balanced accuracy | 0.837 | 0.833 |

Y-randomisation control: true test MCC 0.753 vs. a null distribution of
-0.005 ± 0.057 over 100 label-shuffled permutations (empirical p = 0.0099) —
the model is learning genuine signal, not fitting noise. Applicability
domain: only 2.5% of the test set falls outside the training descriptor
space (Williams leverage).

**Running it:** no CLI arguments — reads `training_data.csv` in the working
directory, writes to `pipeline_outputs/`. On Setonix, submit via `submit.sh`
(SLURM, `gst_ml` conda env, 16 CPUs). Needs `rdkit`, `scikit-learn`,
`imbalanced-learn` (SMOTE), `shap`, `matplotlib` — see `requirements.txt`.

---

## Part 2 — Docking pipeline (`docking/`)

**What it does:** docks the full candidate ligand library against each GST
isoform's real 3D structure using AutoDock Vina, then ranks/filters the
results down to one strong candidate per isoform to carry forward into MD.

**Inputs:**

- Receptors: `1PKW_GSH.pdbqt` (GSTA1-1), `3GUR_GSH.pdbqt` (GSTM1-1),
  `1AQW_GSH.pdbqt` (GSTP1-1) — protein+GSH co-crystal structures, each a
  single-chain monomer, converted to PDBQT. (`legacy_overlaid/` holds an
  earlier shared-coordinate-frame receptor set, superseded once it was
  confirmed to only work correctly for 1 of the 3 proteins — kept for
  reference only.)
- Matching per-isoform grid-box configs (`configs/config_{A1,M1,P1}.txt`) —
  the config's own `receptor=` line is authoritative for Vina, so receptor
  and config filenames must stay paired.
- Ligand library: `prepare_ligands.py` converts SMILES to 3D to PDBQT via
  meeko (`--limit` for a quick test batch, `--ids` to top up with specific
  IDs).

**Running docking:** one SLURM array script per isoform
(`slurm/gsta_docking.sh`, `gstm_docking.sh`, `gstp_docking.sh`), each
wrapping the matching `Vina_rigid_{A,M,P}.pl` script. Ligands are chunked
across array tasks (`CHUNK_SIZE`, default 50) to fit inside the `work`
partition's walltime cap. Submit from `docking/`:

```
N=$(wc -l < ligand.txt); CH=50
sbatch --array=0-$(( (N + CH - 1) / CH - 1 )) slurm/gsta_docking.sh
```

**Candidate selection (`check_warhead_proximity.py`):** ranks docked poses
by two independent signals — the Vina affinity score, and the 3D distance
from the ligand's detected electrophilic warhead atom (via the same SMARTS
patterns the ML pipeline uses) to the GSH thiolate sulfur (`SG2` by
default). A good pose should score well *and* place its warhead near the
reactive sulfur; visual inspection is still worth doing for the shortlist,
since distance alone won't catch steric clashes or strained geometry.

**Why docking doesn't need CGenFF:** Vina treats the receptor as rigid and
scores with an empirical function, so it only needs fixed PDBQT atom
types/charges — no bonded force-field parameters. That requirement only
shows up once a candidate moves on to MD (see below).

---

## Part 3 — Molecular dynamics pipeline (`md/`)

**What it does:** takes one selected ligand (chosen from docking results)
and runs full explicit-solvent MD with GROMACS to see how it behaves in the
GST binding site over real simulated time, using the CHARMM36 protein force
field plus CGenFF parameters for GSH and the ligand (neither is a standard
amino acid, so CHARMM36's built-in library doesn't cover them).

**Step 1 — protein prep (`01_prep_protein.sh`), once per isoform:**
strips GSH/HETATM records, runs `gmx_mpi pdb2gmx` (CHARMM36, SPC/E water,
default charged termini — NH3+/COO-, applied automatically and
non-interactively per chain). `--allow-missing` reconstructs any
crystallographically-disordered side-chain atoms flagged in the PDB's
`REMARK 470`. Produces `protein_processed.gro`, `topol.top`, and `posre.itp`
(the protein's own position-restraint file).

**Step 2 — GSH / ligand CGenFF prep (`02_prep_molecule.sh`), once per
molecule:** takes a "fixed" mol2 (explicit hydrogens, correct protonation —
GSH specifically carries a net -3 charge: three deprotonated carboxylates, a
protonated amine, and a deprotonated Cys thiolate, matching GST's catalytic
mechanism) through `sort_mol2_bonds.pl` and residue renaming, then pauses for
manual CGenFF webserver submission (no local license yet), and finishes by
running `cgenff_charmm2gmx.py` (Lemkul-Lab) to produce GROMACS `itp`/`prm`
files. **The internal CGenFF residue name (`RESI` line inside the `.str`
file) must exactly match the `RESNAME` argument** — GSH's stream files all
declare `RESI LIG` (not `GSH`), docked-ligand stream files declare
`RESI UNL`; check with `grep '^RESI' file.str` before running.

**Step 3 — system build & simulation (`03_run_md.sh`):** assembles the
combined topology (protein + GSH + ligand), merges coordinates, solvates in
a cubic water box, neutralises with `genion` (auto-detects the solvent
group, no hardcoded group index), runs energy minimisation, generates
position restraints for all three components, then NVT → NPT → 10 ns
production MD. Validated end-to-end on 2026-08-10 (GSTA1-1 + GSH +
`lig_03506`): every stage ran clean through NVT/NPT, production confirmed
running at ~43 ns/day before being manually stopped (10 ns interactively
isn't appropriate — see `slurm/run_md.sh`).

```
bash md/03_run_md.sh \
    --protein-dir md/system/gsta1_prepped --gsh-dir md/system/gsta1_gsh_prepped \
    --ligand-dir md/system/gsta1_lig03506_prepped --mdp-dir md/mdp \
    --ffdir ./topology/charmm36-feb2026_cgenff-5.0.ff --outdir md/runs/gsta1_lig03506
```

**Running production MD as a batch job:** `slurm/run_md.sh` wraps the same
call for `sbatch` (needed since a 10 ns run takes hours, not something to
run under an interactive `salloc`):

```
mkdir -p md/logs
sbatch md/slurm/run_md.sh md/system/gsta1_prepped md/system/gsta1_gsh_prepped \
    md/system/gsta1_lig03506_prepped md/runs/gsta1_lig03506
```

Requires `module load gromacs/2024.3-mixed` and the `gst_ml` conda env
(supplies `python3` for the helper scripts in `md/tools/`) on Setonix.

**Two real bugs worth knowing about if you're debugging a new system:**
the protein's own `posre.itp` (written by `pdb2gmx`, not by `03_run_md.sh`)
has to be copied into each run directory or NVT fails immediately; and
`nvt.mdp`/`npt.mdp`/`md.mdp`'s `tc-grps` must reference `Water_and_ions`,
not `Water` — `Water` alone excludes the neutralising ion from any
temperature-coupling group.

**Step 4 — post-MD QC & analysis (`04_analyze_md.sh`):** runs on a
completed production trajectory. Computes QC energy checks, RMSD/RMSF
(protein backbone self-fit, plus ligand and GSH each fit onto the protein),
and the distance from the ligand's detected warhead atom to the GSH
thiolate sulfur (`SG2`) over time — the key mechanistic readout for whether
the covalent geometry is staying plausible during dynamics. A single
invocation also builds a PBC-corrected, water/ion-stripped visualization
trajectory and drives VMD 1.9.3 headlessly (`-dispdev text` + Tachyon
software rendering, no display needed) to produce a snapshot image and a
rendered MP4 of the trajectory via `ffmpeg` — no separate manual step:

```
bash md/04_analyze_md.sh --run-dir md/runs/gsta1_lig03506 \
    --ligand-mol2 md/ligand_mol2_cache/GSTA/3506_A_fixed.mol2
```

**Step 5 — MM-PBSA rescoring (`05_run_mmpbsa.sh`):** rescoring via
`gmx_MMPBSA`, using the single-trajectory (ST) protocol — receptor
(protein+GSH) and ligand are both sliced from the one production
trajectory, no separate apo/free-ligand simulations needed. Runs Poisson-
Boltzmann (PB, not GB — recommended for CHARMM-parameterised systems) with
the dedicated CHARMM radii set (`PBRadii=7`) and `radiopt=0`. As configured
(no `&entropy` block) the result is an *enthalpic* ΔG_bind estimate — useful
for relative ranking between candidates, not a full Gibbs free energy.
**Validated 2026-08-14** on GSTA1-1 + `lig_03506`: a 20-frame slice ran
clean, 0 errors/warnings.

Runs in its own conda env (`gmxMMPBSA`, built on Pawsey scratch —
AmberTools' pinned dependencies are permanently incompatible with the
`gst_ml` env's newer numpy/shap). gmx_MMPBSA's MPI-parallel mode
(`mpirun -np N`) turned out to be unusable on Setonix: it explicitly refuses
`gmx_mpi`/`gmx_mpi_d` when running in parallel, and Setonix's
`gromacs/2024.3-mixed` module only ships `gmx_mpi` (built against Cray
MPICH) — incompatible with the conda env's `mpi4py` (its own bundled MPI).
Worked around with a **SLURM job array**: `slurm/submit_mmpbsa_array.sh`
splits the frame range into independent chunks and submits one fully serial
`gmx_MMPBSA` run per chunk (the exact validated config, `--nranks 1`) as
array tasks — embarrassingly parallel, zero MPI-compatibility risk. Once all
tasks finish, `tools/merge_mmpbsa_chunks.py` pools the per-frame `DELTA
TOTAL` values across chunks into one mean ± SEM (not an average of
per-chunk averages):

```
bash md/slurm/submit_mmpbsa_array.sh md/runs/gsta1_lig03506 300 1000 5 20 8
python3 md/tools/merge_mmpbsa_chunks.py md/runs/gsta1_lig03506/mmpbsa_chunks \
    --out md/runs/gsta1_lig03506/mmpbsa_merged.csv
```

**Reading MM-PBSA against the Vina score:** neither number is a complete
Gibbs free energy — Vina's score is an empirical function regression-fit to
reproduce known experimental affinities (not derived from a physical free
energy calculation), and MM-PBSA here omits the entropy term, so it's an
effective/enthalpic estimate rather than a full ΔG. The two aren't
numerically comparable even though both report kcal/mol. MM-PBSA is used to
check whether Vina's *relative ranking* across candidates holds up once the
pose is allowed to relax and sample dynamically — agreement in trend is the
useful signal, not agreement in magnitude.

---

## Current status (2026-08-18)

- ML pipeline: trained and validated (see Part 1 results above).
- Docking: all three isoforms wired up and run successfully on real
  receptor/config data; candidate selection tooling in place.
- MD: fully automated end-to-end and validated for one system
  (GSTA1-1 + GSH + `lig_03506`); production MD, post-MD QC/RMSD/warhead-
  distance analysis with VMD visualization, and MM-PBSA rescoring are all
  scripted and batch-submittable via SLURM.
- MM-PBSA: validated on the one completed system (20-frame slice, 0
  errors/warnings); now parallelised via a SLURM job array (see Step 5
  above) to scale up the frame count beyond that validation slice.
- Open items: confirm the monomer-vs-dimer MD modeling choice; most ligands
  are still waiting on the CGenFF webserver's daily generation cap (GSTA/
  GSTM have 13/105 parameterised, GSTP has 6/105); repeat docking + MD +
  MM-PBSA for additional candidates as more `.str` files clear; compare the
  MM-PBSA result against the original Vina score as a ranking sanity check.
