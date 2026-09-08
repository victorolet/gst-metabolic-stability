# Implementation notes

Engineering detail behind the docking + MD scripts — bug fixes, environment
quirks, and modelling choices that don't belong in the script comments
themselves. See `pipeline_workflow_diagram.svg` for the pipeline overview
and `../README.md` for the narrative walkthrough; this doc is the "why"
reference for anyone maintaining or debugging the scripts.

## Docking (`docking/slurm/gst{a,m,p}_docking.sh`)

- Split into SLURM array jobs (one array task per `CHUNK_SIZE` ligands, default 50) because Setonix's `work` partition caps walltime at 24h, too short to dock ~3600 ligands serially in one job.
- Ligand PDBQTs are symlinked into scratch per task rather than copied, since there can be thousands.
- `CHUNK_SIZE`/`--time` are first guesses — submit a small array (`--array=0-2`) and check real per-chunk time via `sacct` before scaling up.
- Exit-code gotcha: the scripts exit with `$EXIT_CODE` directly rather than `[ $EXIT_CODE -ne 0 ] && exit $EXIT_CODE`. The latter, when `$EXIT_CODE` is 0, makes the `[ ]` test itself the script's last command — which evaluates false/exit-1 — so SLURM reports a spurious FAILED status even though the run succeeded.
- Receptor/config pairing: Vina reads the receptor path from inside the config file (`receptor=` line), not from the submit script directly, so a receptor and its config must always be swapped together.

## Protein prep (`md/01_prep_protein.sh`)

- Runs `pdb2gmx` **without** `-ter`. `-ter` prompts interactively per chain for N-/C-terminus type; piped stdin answers worked for chain 1 but hung waiting for real input on chain 2 on Setonix's `gmx_mpi` build. Omitting `-ter` uses GROMACS's documented default (charged termini, NH3+/COO-) non-interactively for every chain — which is what the project protocol calls for anyway. `--nterm`/`--cterm` re-enable `-ter` for the rare non-default case; be aware the hang risk returns.
- Input PDB must have GSH (and any other HETATM) already stripped — GSH is chemically a tripeptide, so `pdb2gmx` doesn't reject it, it just tries to assign it its own terminus and gets confused. GSH gets its own prep via `02_prep_molecule.sh`.
- `--allow-missing` passes `-missing` through to let `pdb2gmx` idealize side-chain atoms flagged in the PDB's `REMARK 470`.

## Molecule prep (`md/02_prep_molecule.sh`)

- Same procedure for GSH and any docked ligand. Converts to mol2 (obabel — not in either conda env, resolved via a spack path fallback), sorts bonds (`sort_mol2_bonds.pl`), renames the residue, then pauses for manual CGenFF webserver submission (no local license) before finishing with `cgenff_charmm2gmx.py`.
- **RESNAME must exactly match the `RESI` line inside the paired `.str` file**, not just be a convenient label — confirmed from real delivered files: GSH's stream files all say `RESI LIG` (a generic placeholder, not `GSH`), docked-ligand stream files say `RESI UNL`. Check with `grep '^RESI' file.str` before picking RESNAME. A mismatch fails with a natoms-mismatch error, not a silent wrong answer.
- PyMOL-exported mol2 files carry a leading `# created with PyMOL...` comment line that `sort_mol2_bonds.pl` chokes on (it hard-requires the file to start with `TRIPOS`) — stripped automatically before the perl step.
- LP/LPH lines are **not** stripped by default (unlike the original manual doc), since this build of `cgenff_charmm2gmx.py` has built-in lone-pair support for CGenFF≥4.0 halogens that stripping would disable — relevant to the halogenated warheads this project's `alert_score` descriptor flags. `STRIP_LP=1` forces the old behaviour if ever needed.

## MD system build (`md/03_run_md.sh`)

Two real bugs found and fixed during the first full validation run (GSTA1-1 + GSH + `lig_03506`, 2026-08-10):

- The protein's own `posre.itp` (written by `pdb2gmx`, not by this script) has to be copied into the run directory or NVT fails immediately — `topol.top`'s `#ifdef POSRES` block expects it locally.
- `nvt.mdp`/`npt.mdp`/`md.mdp`'s `tc-grps` must reference `Water_and_ions`, not `Water` — `Water` alone excludes the neutralising ion from any temperature-coupling group.

Other notes: `genion`'s SOL group index is auto-detected from `make_ndx` output rather than hardcoded, since it isn't guaranteed stable across systems. GSH's GROMACS residue/group name is `LIG` throughout (per the RESI-naming note above), not `GSH` — anything grepping index-group listings for `GSH` will find nothing.

Validated end-to-end 2026-08-10: all stages ran clean through NVT/NPT, production confirmed running at ~43 ns/day.

## Post-MD analysis (`md/04_analyze_md.sh`)

- The ligand mol2 used for CGenFF has non-unique atom names (bare `C`, `N`, `H` repeated), so the warhead atom can't be selected by name. `tools/find_warhead_atom.py` returns the warhead's local (1-based) atom number within the ligand's own mol2/`unl.gro`; the script converts that to a whole-system atom number by adding the protein and GSH atom counts that precede it in `complex.gro`'s fixed ordering (protein, then GSH, then ligand — set by `merge_gro.py` in `03_run_md.sh`, never reordered by `solvate`/`genion`, which only append). GSH's `SG2` stays name-based since GSH's mol2 names are unique.
- Visualization is fully automatic: builds a PBC-corrected, water/ion-stripped, stride-reduced trajectory, then drives VMD 1.9.3 headlessly (`-dispdev text` + `TachyonInternal` software rendering) for a snapshot + frame sequence, then `ffmpeg` stitches the frames into an MP4 — one invocation, no manual follow-up. `vmd_load_trajectory.tcl` is also written standalone (no rendering, no `quit`) so the trajectory can be watched live in VMD's own GUI separately.
- `--skip-render` stops after the `.pdb`/`.xtc`/`.tcl` files if VMD or ffmpeg aren't available in a given session.

## MM-PBSA rescoring (`md/05_run_mmpbsa.sh`)

**Modelling choice (not obviously correct either way — flag if this should change):** receptor = Protein + GSH, ligand = the docked candidate (UNL). This mirrors how docking treats the system throughout the project (protein+GSH as one fixed unit, the xenobiotic as the binding partner). A 3-way protein/GSH/ligand split would need gmx_MMPBSA's more involved multicomponent protocol.

**Why PB, not GB:** gmx_MMPBSA's own CHARMM documentation recommends PB for CHARMM-prepped systems — GB/PB radii sets built for AMBER atom types aren't well tested against CHARMM. Uses the dedicated `charmm_radii` set (`PBRadii=7`) and `radiopt=0` (radii from the topology, not recomputed), both per gmx_MMPBSA's own recommendation.

**Environment:** runs in its own conda env (`gmxMMPBSA`, built on Pawsey scratch), separate from `gst_ml` — AmberTools' pinned numpy/pandas/scipy/matplotlib versions are permanently incompatible with `gst_ml`'s newer numpy/shap requirement. Built on `$MYSCRATCH` rather than `/software` after hitting a **file-count** quota there (241k/250k files; space itself was trivial, 60GB/16TB) — conda environments with heavy cheminformatics stacks (RDKit, AmberTools) are inherently file-count-expensive.

**MPI incompatibility (why the SLURM job-array version exists):** gmx_MMPBSA's parallel mode (`mpirun -np N`) requires `mpi4py` and GROMACS to share the same underlying MPI library. Setonix's `gromacs/2024.3-mixed` module only ships `gmx_mpi` (built against Cray MPICH), and gmx_MMPBSA explicitly refuses `gmx_mpi`/`gmx_mpi_d` for `-np>1` — confirmed no plain `gmx` binary exists in that module directory either. Real MPI parallelism isn't available here without a rebuild. `md/slurm/run_mmpbsa_array.sh` (submitted via `submit_mmpbsa_array.sh`) works around this by splitting the frame range into independent chunks and running one fully serial `gmx_MMPBSA` per chunk as SLURM array tasks — embarrassingly parallel, zero MPI risk. `tools/merge_mmpbsa_chunks.py` then pools the per-frame `DELTA TOTAL` values across chunks into one mean ± SEM (not an average of per-chunk averages).

**Path bugs fixed during the first real run:** `-cs` needs `md_0_5.tpr`'s absolute path (it lives in `$RUN_DIR`, one level above `$OUT_DIR`, where the script `cd`s to run gmx_MMPBSA). Separately, the CHARMM force-field directory (`*.ff`) has to be symlinked into `$OUT_DIR` — gmx_MMPBSA re-parses `topol.top` via ParmEd rather than just reading the `.tpr`, so every `#include` (the ff directory's own `forcefield.itp` included) must resolve from wherever `topol.top` sits.

**Validated** 2026-08-14 on GSTA1-1 + `lig_03506`: 20-frame slice (500–520), 0 errors/warnings, ~4.5 min single-core. This is a deliberately small validation slice, not a converged estimate — published MM-PBSA studies typically use 50–200+ frames.

**Reading the result against Vina's score:** neither is a complete Gibbs free energy. Vina's score is an empirical function regression-fit to reproduce known experimental affinities, not derived from a physical free energy calculation. MM-PBSA here omits the entropy term (no `&entropy` block), so it's an effective/enthalpic estimate. The two aren't numerically comparable even though both report kcal/mol — use MM-PBSA to check whether Vina's *relative ranking* across candidates holds up, not to match magnitudes.

## Environment storage (`env_tools/`)

Pawsey's `/software` quota is a **file-count** limit (100k files), not a space limit -- a single mature conda environment (RDKit, AmberTools, etc.) can be tens of thousands of files.

**SquashFS/FUSE does NOT work on Setonix for this -- confirmed 2026-09-07, do not retry without a different approach.** Pawsey's own documentation describes packing such environments into a SquashFS image, but their documented path loads it *inside an Apptainer/Singularity container* (`singularity exec --overlay`), which uses Apptainer's own privileged, admin-installed mount helper. A personal, spack-built `squashfuse` (the standalone route, which would have let the rest of the pipeline keep using plain `conda activate` unchanged) builds fine but fails to actually mount: `fusermount3: mount failed: Operation not permitted`. Setonix blocks unprivileged FUSE mounts at the platform level for regular users -- this is a site policy, not a config problem on our end. `squash_conda_env.sh` and `mount_conda_env.sh` are kept in this directory for reference but are **not used** by any script.

**Current approach: scratch + periodic Acacia backup**, applied uniformly to both `gst_ml` and `gmxMMPBSA`:
- Both envs live on `$MYSCRATCH/conda_envs/<name>` (no file-count quota there, only a 21-day purge policy).
- `env_tools/backup_env_to_acacia.sh` syncs a live env to Acacia (Pawsey's S3-compatible object storage, `rclone`-based, space-quota not file-count-quota, see `acacia_rclone.conf.example` for the one-time rclone setup) as insurance against that purge. Defaults to a dry run; `--do-it` actually uploads.
- `env_tools/restore_env_from_acacia.sh` restores a backup, always to the *same* path the env was built at (conda environments bake absolute paths into scripts/shebangs/conda-meta, so restoring to a different path can silently break things) and refuses to overwrite a directory that still has real content.
- Not yet automated/scheduled -- run `backup_env_to_acacia.sh` manually after a meaningful env change (e.g. once `gmxMMPBSA` or `gst_ml` gets new packages installed).

**rclone drops symlinks by default.** `backup_env_to_acacia.sh` uses `--copy-links` for this reason -- conda envs are full of symlinks (`bin/python`, shared libraries), and a plain `rclone sync`/`copy` silently skips them, producing an incomplete backup with no error. This also caught out a restore of `md/runs/gsta1_lig03506` from an older, unrelated `volet-pipeline-development` Acacia backup (2026-09-09): `03_run_md.sh` sets up the CHARMM force-field directory *inside* each run directory as a symlink back to `topology/`, which never made it into that backup. Restoring run data from Acacia will need that symlink recreated by hand (`ln -sfn <topology dir> <run dir>/<ffname>.ff`) regardless of what gets pulled down.

**Env locations aren't where `.condarc`'s `envs_dirs` might suggest.** `envs_dirs` only affects where *new* environments get created/searched by bare name -- conda always additionally searches `<install prefix>/envs` regardless. `gst_ml` was originally created (before the scratch-cache `.condarc` change made during the gmxMMPBSA saga) at `/software/projects/pawsey1376/volet/miniconda3/envs/gst_ml`, not under the `conda_envs/` directory `envs_dirs` points at -- `conda_envs/` only ever ended up holding `vina`. Check with `ls` before assuming a path. `gst_ml` has since been relocated to `$MYSCRATCH/conda_envs/gst_ml` per the above.

- `env_tools/squash_conda_env.sh <env_dir> <out_sqfs_path>` builds the image and smoke-tests it via a temporary mount. It never touches the original environment -- swapping the live directory for the squashed one is a separate, manual, printed step, so you can verify first.
- `env_tools/mount_conda_env.sh` (source it, call `mount_conda_env <env_dir> <sqfs_path>`) is the idempotent mount-before-activate helper wired into `run_md.sh` (`gst_ml`) and `run_mmpbsa.sh`/`run_mmpbsa_array.sh` (`gmxMMPBSA`). It only mounts when the `.sqfs` exists **and** `env_dir` is currently empty/missing -- a live directory with real content is never shadowed. This means it's safe in every script regardless of whether an env has actually been squashed yet, and it doubles as disaster recovery: `gmxMMPBSA` lives on scratch (21-day purge), so if the live copy is ever purged, the next job run automatically restores it from the `/software`-backed `.sqfs` instead of failing.
- SquashFS images are read-only, so this suits an environment that's largely done changing (`gmxMMPBSA` was fully validated before being squashed) rather than one still being actively `pip install`-ed into.
- Squashfuse mounts don't persist across a new login shell or a new SLURM job's node allocation -- every session/job re-mounts via `mount_conda_env`, which is why it's called at the top of each SLURM script rather than assumed to already be in place.

## Recurring conventions worth knowing

- **Resolve paths to absolute before any `cd`.** Every script that changes directory (`02_prep_molecule.sh`, `03_run_md.sh`, `04_analyze_md.sh`, `05_run_mmpbsa.sh`) resolves its input/output directory arguments to absolute paths first. A relative path silently breaks once the working directory has moved — this exact bug has recurred enough times across the pipeline to be worth calling out generally rather than per-script.
- **SLURM exit codes:** every batch wrapper exits with `$EXIT_CODE` captured right after the underlying command, not a conditional test on it (see the docking note above for why the conditional form produces false FAILED statuses).
- **SLURM array wrappers needing a computed `--array` range** (`run_mmpbsa_array.sh`) are paired with a plain-bash `submit_*.sh` launcher, since `--array=X-Y` has to be known at `sbatch` submit time and can't be computed inside the batch script itself from positional args.
