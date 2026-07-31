#!/usr/bin/env python3
# =============================================================================
# ML Pipeline for GST Stability Prediction of Covalent Ligands
# =============================================================================
#
# Predicts whether a covalent ligand is metabolically UNSTABLE to
# glutathione S-transferase (GST)-mediated glutathione (GSH) conjugation.
#
# Predictive model: SOFT-VOTING ENSEMBLE (Random Forest + SVM + Gradient
# Boosting) trained on RDKit molecular descriptors PLUS
# mechanism-based electrophilicity descriptors (global electrophilicity
# index, warhead-localised local electrophilicity, structural-alert warhead
# count).
#
# Pipeline stages:
# 1. Compute RDKit + electrophilicity descriptors for every ligand
# 2. Train/test split + train the soft-voting ensemble ML model (with SMOTE)
# 3. Evaluate on test set, run controls, and produce benchmarking outputs
#
# =============================================================================

import csv
import os
import sys
import shutil
import subprocess
import tempfile
import time
import numpy as np
import pandas as pd
import pickle
from typing import Dict, List, Tuple, Optional, Set
from functools import lru_cache

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from datetime import datetime

# sklearn
from sklearn.model_selection import cross_val_score, train_test_split, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier, VotingClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_auc_score, matthews_corrcoef,
                              roc_curve, auc, make_scorer,
                              precision_recall_curve)
from sklearn.inspection import permutation_importance

# RDKit
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem import rdPartialCharges

from tqdm import tqdm

# Optional SMOTE (depending on input data)
try:
    from imblearn.over_sampling import SMOTE
    SMOTE_AVAILABLE = True
except ImportError:
    SMOTE = None
    SMOTE_AVAILABLE = False
    print("WARNING: imbalanced-learn not installed. SMOTE will be disabled.")

import shap

RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]

# ===============
# GLOBAL SETTINGS
# ===============

BASE_DIR = os.path.abspath(os.getcwd())

GLOBAL_SEED = 42

# =======================
# ENSEMBLE MODEL REGISTRY
# =======================

ENSEMBLE_CONFIG = {
    "random_forest": {
        "estimator": RandomForestClassifier,
        "params": {
            "n_estimators": 300,
            "max_depth": 3,
            "min_samples_split": 500,
            "min_samples_leaf": 250,
            "class_weight": "balanced",
            "random_state": GLOBAL_SEED,
            "n_jobs": -1,
        },
        "weight": 1.0,
    },
    "svm": {
        "estimator": SVC,
        "params": {
            "kernel": "rbf",
            "C": 0.3,
            "gamma": 0.05,
            "class_weight": "balanced",
            "probability": True,
            "random_state": GLOBAL_SEED,
        },
        "weight": 1.0,
    },
    "gradient_boosting": {
        "estimator": GradientBoostingClassifier,
        "params": {
            "n_estimators": 500,
            "max_depth": 3,
            "learning_rate": 0.02,
            "subsample": 0.6,
            "min_samples_split": 500,
            "min_samples_leaf": 250,
            "n_iter_no_change": 20,
            "random_state": GLOBAL_SEED,
        },
        "weight": 1.0,
    },
}


def build_ensemble_model(config: dict = None) -> VotingClassifier:
    if config is None:
        config = ENSEMBLE_CONFIG
    estimators = [(name, spec["estimator"](**spec["params"]))
                  for name, spec in config.items()]
    weights = [spec["weight"] for spec in config.values()]
    return VotingClassifier(estimators=estimators, voting="soft",
                             weights=weights, n_jobs=-1)


# Feature set: RDKit descriptors + mechanism-based electrophilicity descriptors.
DESCRIPTOR_NAMES = ["SMR_VSA7", "SMR_VSA9", "PEOE_VSA9", "SLogP_VSA1", "CSP3"]
ELECTRO_FEATURE_NAMES = ["omega_global", "omega_local_warhead", "alert_score"]
ALL_FEATURE_NAMES = DESCRIPTOR_NAMES + ELECTRO_FEATURE_NAMES

# SMOTE
USE_SMOTE = SMOTE_AVAILABLE
SMOTE_RATIO = 0.4

# Train/test split + CV
TEST_SPLIT = 0.2
CV_FOLDS = 10

# Backend selection for the electrophilicity descriptors: "auto" -> xtb if
# available (shutil.which("xtb")) else the RDKit surrogate.
ELECTRO_BACKEND = "auto"

# Decision-threshold grid (MCC-tuned on train).
THRESHOLD_GRID = np.linspace(0.0, 1.0, 101)

# Number of label-permutation runs for the y-randomisation control.
Y_RANDOM_RUNS = 100

# SHAP settings
SHAP_BACKGROUND_SIZE = 50
SHAP_SAMPLE_SIZE = 200
SHAP_NSAMPLES = 100

# Permutation importance settings.
PERM_IMPORTANCE_REPEATS = 30

# Directories.
OUTPUT_DIR = os.path.join(BASE_DIR, "pipeline_outputs")
TRAINING_DATA_FILE = os.path.join(BASE_DIR, "training_data.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ========================
# IN-MEMORY COMPUTE CACHES
# ========================

_ELECTRO_CACHE: Dict[str, dict] = {}  # canonical SMILES -> electrophilicity record

# ============
# DATA LOADING
# ============

def load_training_data(filepath: str) -> pd.DataFrame:
    # Load labelled SMILES; flexible column detection; binarise label.
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Training data file not found: {filepath}")
    df = pd.read_csv(filepath)
    print(f"[OK] File loaded: {len(df)} rows")

    column_map = {}
    for canon, cands in [
        ("ID", ["ID", "id", "ligand_id", "Ligand_ID"]),
        ("SMILES", ["SMILES", "Smiles", "smiles", "SMILE"]),
        ("GROUND_TRUTH", ["GROUND_TRUTH", "Label", "label", "Stability"]),
    ]:
        col = next((c for c in cands if c in df.columns), None)
        if col:
            column_map[canon] = col
    missing = [c for c in ("ID", "SMILES", "GROUND_TRUTH") if c not in column_map]
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df.rename(columns={v: k for k, v in column_map.items()})
    df["Ligand_ID"] = df["ID"].astype(int)

    if df["GROUND_TRUTH"].dtype == "object":
        stable_vals = {"stable", "Stable", "yes", "true", "active", "1"}
        unstable_vals = {"unstable", "Unstable", "no", "false", "inactive", "0"}

        def map_label(v):
            s = str(v).strip()
            return 1 if s in stable_vals else (0 if s in unstable_vals else None)

        df["Label"] = df["GROUND_TRUTH"].apply(map_label)
    else:
        df["Label"] = df["GROUND_TRUTH"].astype(int)

    df = df.dropna(subset=["SMILES", "Label"])
    n = len(df)
    n_stable = int((df["Label"] == 1).sum())
    print(f"[OK] Loaded {n} ligands | Stable={n_stable} | Unstable={n - n_stable}")
    return df


def load_ground_truth(filepath: str = TRAINING_DATA_FILE) -> Dict[int, str]:
    # Map ligand id -> "stable"/"unstable" for evaluation.
    df = pd.read_csv(filepath)
    id_col = next((c for c in ["ID", "id", "ligand_id", "Ligand_ID"] if c in df.columns), None)
    gt_col = next((c for c in ["GROUND_TRUTH", "Label", "label", "Stability"] if c in df.columns), None)
    if not id_col or not gt_col:
        raise ValueError("Cannot find ID and GROUND_TRUTH columns")
    gt = {}
    for _, row in df.iterrows():
        lid = int(row[id_col])
        val = str(row[gt_col]).strip().lower()
        if val in {"stable", "1", "yes", "true", "active"}:
            gt[lid] = "stable"
        elif val in {"unstable", "0", "no", "false", "inactive"}:
            gt[lid] = "unstable"
    return gt


def filter_ground_truth(gt: Dict[int, str], ids: Set[int]) -> Dict[int, str]:
    return {lid: lab for lid, lab in gt.items() if lid in ids}


# =============================
# PART 2: MOLECULAR DESCRIPTORS
# =============================

@lru_cache(maxsize=4000)
def compute_molecular_descriptors_cached(smiles: str) -> Optional[tuple]:
    # Compute the descriptor tuple for one SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    try:
        return (
            Descriptors.SMR_VSA7(mol),
            Descriptors.SMR_VSA9(mol),
            Descriptors.PEOE_VSA9(mol),
            Descriptors.SlogP_VSA1(mol),
            Descriptors.FractionCSP3(mol),
        )
    except Exception:
        return None


def compute_descriptors(smiles: str) -> Optional[Dict[str, float]]:
    res = compute_molecular_descriptors_cached(smiles)
    return dict(zip(DESCRIPTOR_NAMES, res)) if res is not None else None


# ==================================
# PART 3: ELECTROPHILICITY FEATURES
# ==================================
# Mechanism-based soft-electrophile reactivity toward the GSH thiolate.
# 1. omega_global        : global electrophilicity index
#                           omega = mu^2 / (2*eta),
#                           mu = (eLUMO + eHOMO) / 2, eta = (eLUMO - eHOMO)
#                           (Parr, Szentpaly, Liu 1999).
# 2. omega_local_warhead  : local electrophilicity omega * f+ (condensed
#                           Fukui function; Wondrousch et al. 2010)
# 3. alert_score          : total count of GSH-reactive structural-alert
#                           warhead centres in the molecule (RDKit SMARTS)
#
# Backend: xTB GFN2 frontier orbitals + Fukui function (--vfukui) when the
# `xtb` binary is on PATH; otherwise a deterministic RDKit Gasteiger-charge
# surrogate.
#
# Literature basis:
#   Parr, Szentpaly, Liu 1999. doi:10.1021/ja983494x
#   Wondrousch, Bohme, Thaens, Ost, Schuurmann 2010. doi:10.1021/jz100247x
#   Mayer & Ofial 2019. doi:10.1002/anie.201909803
#   Bohme, Thaens, Paschke, Schuurmann 2009. doi:10.1021/tx900044e
#   Miller, Hughes, Swamidass 2015. doi:10.1021/acs.chemrestox.5b00017
#   Bannwarth, Ehlert, Grimme 2019 (GFN2-xTB). doi:10.1021/acs.jctc.8b01176

# --------
# Backends
# --------

def _xtb_available() -> bool:
    return shutil.which("xtb") is not None


def _embed_3d(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()  # type: ignore[attr-defined]
        params.randomSeed = GLOBAL_SEED
        cid = AllChem.EmbedMolecule(mol, params)  # type: ignore[attr-defined]
        if cid < 0:
            params.useRandomCoords = True
            cid = AllChem.EmbedMolecule(mol, params)  # type: ignore[attr-defined]
            if cid < 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(mol)  # type: ignore[attr-defined]
        except Exception:
            pass
        if mol.GetNumConformers() == 0:
            return None
        return mol
    except Exception:
        return None


def _xtb_frontier_and_fukui(mol) -> Optional[dict]:
    if not _xtb_available():
        return None
    heavy_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    n_heavy = len(heavy_idx)
    if n_heavy == 0:
        return None
    tmpdir = tempfile.mkdtemp(prefix="xtb_")
    try:
        xyz_path = os.path.join(tmpdir, "mol.xyz")
        Chem.MolToXYZFile(mol, xyz_path)
        try:
            proc = subprocess.run(
                ["xtb", "mol.xyz", "--gfn", "2", "--vfukui", "--chrg", "0", "--uhf", "0"],
                cwd=tmpdir, timeout=120, capture_output=True, text=True, check=False)
        except Exception:
            return None
        out = proc.stdout or ""
        if not out:
            return None
        eHOMO = eLUMO = None
        for line in out.splitlines():
            if "(HOMO)" in line:
                for p in line.split():
                    try:
                        eHOMO = float(p)
                    except ValueError:
                        continue
            elif "(LUMO)" in line:
                for p in line.split():
                    try:
                        eLUMO = float(p)
                    except ValueError:
                        continue
        if eHOMO is None or eLUMO is None:
            return None
        fplus = None
        fukui_file = os.path.join(tmpdir, "fukui")
        if os.path.exists(fukui_file):
            try:
                vals = []
                with open(fukui_file) as fh:
                    for line in fh:
                        toks = line.split()
                        if len(toks) >= 2:
                            try:
                                vals.append(float(toks[1]))
                            except ValueError:
                                continue
                if len(vals) >= n_heavy:
                    fplus = np.array(vals[:n_heavy], dtype=float)
            except Exception:
                fplus = None
        if fplus is None:
            in_block = False
            vals = []
            for line in out.splitlines():
                if "f(+)" in line and "f(-)" in line:
                    in_block = True
                    continue
                if in_block:
                    toks = line.split()
                    if len(toks) >= 3 and toks[0].isdigit():
                        try:
                            vals.append(float(toks[2]))
                        except (ValueError, IndexError):
                            continue
                    elif toks and not toks[0].isdigit():
                        if vals:
                            break
            if len(vals) >= n_heavy:
                fplus = np.array(vals[:n_heavy], dtype=float)
        if fplus is None or not np.all(np.isfinite(fplus)):
            return None
        return {"eHOMO": float(eHOMO), "eLUMO": float(eLUMO), "fplus": fplus}
    except Exception:
        return None
    finally:
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass


def _rdkit_surrogate_frontier_and_fukui(mol, electrophilic_atoms=None) -> dict:
    try:
        m = Chem.Mol(mol)
        rdPartialCharges.ComputeGasteigerCharges(m)
        heavy_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
        pos_full = np.zeros(m.GetNumAtoms(), dtype=float)
        for a in m.GetAtoms():
            try:
                c = float(a.GetProp("_GasteigerCharge"))
            except (KeyError, ValueError):
                c = 0.0
            if not np.isfinite(c):
                c = 0.0
            pos_full[a.GetIdx()] = max(c, 0.0)

        pos_heavy = np.array([pos_full[i] for i in heavy_idx], dtype=float) if heavy_idx else np.zeros(1)
        omega_global = float(np.sum(pos_heavy) * 2.0)  # deterministic scale
        if not np.isfinite(omega_global):
            omega_global = 0.0

        elec_set = set(electrophilic_atoms) if electrophilic_atoms else set()
        raw = np.zeros(len(heavy_idx), dtype=float)
        if elec_set:
            for h, idx in enumerate(heavy_idx):
                if idx in elec_set:
                    env = pos_full[idx]
                    atom = m.GetAtomWithIdx(idx)
                    for nb in atom.GetNeighbors():
                        env += pos_full[nb.GetIdx()]
                    raw[h] = max(env, 0.10)
        else:
            raw = pos_heavy.copy()

        max_raw = float(np.max(raw)) if raw.size else 0.0
        fplus = raw / max_raw if max_raw > 0 else np.zeros_like(raw)
        fplus = np.nan_to_num(fplus, nan=0.0, posinf=0.0, neginf=0.0)
        return {"omega_global": omega_global, "fplus": fplus}
    except Exception:
        n_heavy = max(1, sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1))
        return {"omega_global": 0.0, "fplus": np.zeros(n_heavy)}


# ------------------------------------------------------
# Structural-alert SMARTS (GSH-reactive warhead classes)
# ------------------------------------------------------

_ELECTRO_WARHEAD_ATOM_IDX = 0

ELECTRO_SMARTS = {
    "snar_p_no2": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[N+](=O)[O-]"), 2.0),
    "snar_p_cn": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC#N"), 1.0),
    "snar_p_cf3": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC(F)(F)F"), 1.0),
    "snar_p_sulfone": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cS(=O)(=O)[C,c,N]"), 1.0),
    "snar_p_carbonyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:cC(=O)"), 1.0),
    "snar_p_ammonium": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[N+]"), 1.0),
    "snar_p_noxide": (Chem.MolFromSmarts("[n+]1([O-])ccc([F,Cl,Br])cc1"), 1.0),
    "snar_o_no2": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[N+](=O)[O-]"), 2.0),
    "snar_o_cn": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC#N"), 1.0),
    "snar_o_cf3": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC(F)(F)F"), 1.0),
    "snar_o_sulfone": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cS(=O)(=O)[C,c,N]"), 1.0),
    "snar_o_carbonyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):cC(=O)"), 1.0),
    "snar_o_ammonium": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[N+]"), 1.0),
    "snar_o_noxide": (Chem.MolFromSmarts("[n+]1([O-])c([F,Cl,Br])cccc1"), 1.0),
    "snar_2_pyridyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n"), 1.5),
    "snar_3_pyridyl": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:n"), 1.5),
    "snar_pyridazine": (Chem.MolFromSmarts("[c]([F,Cl,Br])(:n):n"), 1.5),
    "snar_pyrimidine": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n:c:n"), 1.5),
    "snar_trazine": (Chem.MolFromSmarts("[c]([F,Cl,Br]):n:c:n:c:n"), 1.5),
    "snar_o_trihalo": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c[F,Cl,Br,I]"), 1.0),
    "snar_p_trihalo": (Chem.MolFromSmarts("[c]([F,Cl,Br]):c:c:c[F,Cl,Br,I]"), 1.0),
    "snar_sulfonate_LG1": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfonate_LG2": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC#N"), 1.0),
    "snar_sulfonate_LG3": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC(F)(F)F"), 1.0),
    "snar_sulfonate_LG4": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:cC(=O)"), 1.0),
    "snar_sulfonate_LG5": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:c[N+]"), 1.0),
    "snar_sulfonate_LG6": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c[N+](=O)[O-]"), 1.0),
    "snar_sulfonate_LG7": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):cC#N"), 1.0),
    "snar_sulfonate_LG8": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):cC(=O)"), 1.0),
    "snar_sulfonate_LG9": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):n"), 1.0),
    "snar_sulfonate_LG10": (Chem.MolFromSmarts("[c](OS(=O)(=O)[C,c]):c:c:n"), 1.0),
    "snar_sulfone_LG1": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG2": (Chem.MolFromSmarts("[c](S(=O)(=O)c):c:c:c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG3": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:cC#N"), 1.0),
    "snar_sulfone_LG4": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:cC(=O)"), 1.0),
    "snar_sulfone_LG5": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c[N+](=O)[O-]"), 1.0),
    "snar_sulfone_LG6": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):cC#N"), 1.0),
    "snar_sulfone_LG7": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):cC(=O)"), 1.0),
    "snar_sulfone_LG8": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):n"), 1.0),
    "snar_sulfone_LG9": (Chem.MolFromSmarts("[c](S(=O)(=O)[C,c]):c:c:n"), 1.0),
    "ab_unsat_carbonyl": (Chem.MolFromSmarts("[CX3:1]=[CX3]C=O"), 1.0),
    "acrylate_ester": (Chem.MolFromSmarts("[CX3:1]=[CX3]C(=O)[O,N]"), 1.0),
    "acrylonitrile": (Chem.MolFromSmarts("[CX3:1]=[CX3]C#N"), 1.0),
    "vinyl_sulfone": (Chem.MolFromSmarts("[CX3:1]=[CX3]S(=O)(=O)[C,c]"), 1.0),
    "nitroalkene": (Chem.MolFromSmarts("[CX3:1]=[CX3][N+](=O)[O-]"), 1.0),
    "divinyl_enone": (Chem.MolFromSmarts("[CX3:1]=[CX3]C(=O)C"), 1.0),
    "prim_ahalide": (Chem.MolFromSmarts("[CX4;!$(C(F)(F)F):1][Cl,Br,I]"), 1.0),
    "asulfonate": (Chem.MolFromSmarts("[CX4:1]OS(=O)(=O)[C,c]"), 1.0),
    "atriflate": (Chem.MolFromSmarts("[CX4:1]OS(=O)(=O)C(F)(F)F"), 1.0),
    "aammonium": (Chem.MolFromSmarts("[CH2X4:1][N+](C)(C)C"), 1.0),
    "epoxide": (Chem.MolFromSmarts("[CX4:1]1[OX2][CX4]1"), 1.5),
}


def _electrophilic_atom_indices(mol) -> List[int]:
    idxs: Set[int] = set()
    for patt, _weight in ELECTRO_SMARTS.values():
        if patt is None:
            continue
        try:
            for match in mol.GetSubstructMatches(patt):
                if match:
                    idxs.add(match[0])
        except Exception:
            continue
    if idxs:
        return sorted(idxs)
    try:
        m = Chem.Mol(mol)
        rdPartialCharges.ComputeGasteigerCharges(m)
        fallback = []
        for a in m.GetAtoms():
            if a.GetAtomicNum() <= 1:
                continue
            try:
                c = float(a.GetProp("_GasteigerCharge"))
            except (KeyError, ValueError):
                continue
            if np.isfinite(c) and c > 0.05:
                fallback.append(a.GetIdx())
        return fallback
    except Exception:
        return []


def _alert_score(mol) -> float:
    centres: Set[int] = set()
    for patt, _weight in ELECTRO_SMARTS.values():
        if patt is None:
            continue
        try:
            for match in mol.GetSubstructMatches(patt):
                if match:
                    centres.add(match[_ELECTRO_WARHEAD_ATOM_IDX])
        except Exception:
            continue
    return float(len(centres))


# -------------------
# Descriptor assembly
# -------------------

def compute_electrophilicity_for_ligand(smiles: str) -> Optional[dict]:
    mol = _embed_3d(smiles)
    if mol is None:
        return None

    backend = "rdkit_surrogate"
    omega_global = 0.0
    fplus_heavy = None

    mol_noH = Chem.RemoveHs(mol)
    n_heavy_noH = mol_noH.GetNumAtoms()
    electrophilic_atoms = _electrophilic_atom_indices(mol_noH)

    if ELECTRO_BACKEND in ("auto", "xtb_gfn2") and _xtb_available():
        xtb_res = _xtb_frontier_and_fukui(mol)
        if xtb_res is not None:
            eHOMO, eLUMO = xtb_res["eHOMO"], xtb_res["eLUMO"]
            mu = (eLUMO + eHOMO) / 2.0
            eta = (eLUMO - eHOMO)
            if np.isfinite(mu) and np.isfinite(eta) and abs(eta) > 1e-9:
                omega_global = float((mu ** 2) / (2.0 * eta))
            fplus_heavy = xtb_res["fplus"]
            backend = "xtb_gfn2"

    if fplus_heavy is None:
        surrogate = _rdkit_surrogate_frontier_and_fukui(
            mol, electrophilic_atoms=electrophilic_atoms)
        omega_global = float(surrogate["omega_global"])
        fplus_heavy = surrogate["fplus"]
        backend = "rdkit_surrogate"

    if not np.isfinite(omega_global):
        omega_global = 0.0

    fplus_heavy = np.asarray(fplus_heavy, dtype=float)
    fplus_heavy = np.nan_to_num(fplus_heavy, nan=0.0, posinf=0.0, neginf=0.0)
    if fplus_heavy.size != n_heavy_noH:
        if fplus_heavy.size > n_heavy_noH:
            fplus_heavy = fplus_heavy[:n_heavy_noH]
        else:
            fplus_heavy = np.pad(fplus_heavy, (0, n_heavy_noH - fplus_heavy.size))

    omega_local = omega_global * fplus_heavy

    warhead_atoms: Set[int] = set()
    for patt, _weight in ELECTRO_SMARTS.values():
        if patt is None:
            continue
        try:
            for match in mol_noH.GetSubstructMatches(patt):
                if match:
                    warhead_atoms.add(match[_ELECTRO_WARHEAD_ATOM_IDX])
        except Exception:
            continue
    warhead_atoms &= set(electrophilic_atoms) if electrophilic_atoms else warhead_atoms
    valid_warhead = [a for a in warhead_atoms if a < omega_local.size]
    omega_local_warhead = float(np.max(omega_local[valid_warhead])) if valid_warhead else 0.0
    if not np.isfinite(omega_local_warhead):
        omega_local_warhead = 0.0

    alert_score = _alert_score(mol_noH)

    feats = [omega_global, omega_local_warhead, alert_score]
    all_finite = all(np.isfinite(v) for v in feats)
    has_reactive = bool(all_finite and (alert_score > 0 or omega_local_warhead > 0))

    return {
        "omega_global": float(omega_global),
        "omega_local_warhead": float(omega_local_warhead),
        "alert_score": float(alert_score),
        "has_reactive": has_reactive,
        "backend": backend,
    }


def compute_electrophilicity_all_ligands(smiles_dict: Dict[int, str]) -> Dict[int, dict]:
    records: Dict[int, dict] = {}
    for lid, smi in tqdm(sorted(smiles_dict.items()), desc="Electrophilicity",
                          unit="lig", file=sys.stdout):
        mol = Chem.MolFromSmiles(smi)
        canon = Chem.MolToSmiles(mol) if mol is not None else smi
        if canon in _ELECTRO_CACHE:
            rec = _ELECTRO_CACHE[canon]
        else:
            rec = compute_electrophilicity_for_ligand(smi)
            _ELECTRO_CACHE[canon] = rec
        if rec is not None:
            records[lid] = rec
    return records


def build_full_feature_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[int], Dict[int, dict]]:
    # Build the combined RDKit-descriptor + electrophilicity feature matrix.
    smiles_dict = {int(r["Ligand_ID"]): r["SMILES"] for _, r in df.iterrows()}
    electro_records = compute_electrophilicity_all_ligands(smiles_dict)

    X_list, valid_indices, failed = [], [], 0
    for idx, row in df.iterrows():
        lid = int(row["Ligand_ID"])
        desc = compute_descriptors(row["SMILES"])
        erec = electro_records.get(lid)
        if desc is None or erec is None:
            failed += 1
            continue
        feats = [desc[name] for name in DESCRIPTOR_NAMES]
        feats += [erec[name] for name in ELECTRO_FEATURE_NAMES]
        if not all(np.isfinite(v) for v in feats):
            failed += 1
            continue
        X_list.append(feats)
        valid_indices.append(idx)
    X = np.array(X_list)
    print(f"[OK] Combined descriptors: {len(valid_indices)}/{len(df)} ligands "
          f"(failed: {failed}), shape {X.shape}")
    return X, valid_indices, electro_records


# ==============
# SHARED METRICS
# ==============

def recall_score(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred == 1)); fn = np.sum((y_true == 1) & (y_pred == 0))
    return tp / (tp + fn) if (tp + fn) else 0.0


def specificity_score(y_true, y_pred):
    tn = np.sum((y_true == 0) & (y_pred == 0)); fp = np.sum((y_true == 0) & (y_pred == 1))
    return tn / (tn + fp) if (tn + fp) else 0.0


def precision_score_custom(y_true, y_pred):
    tp = np.sum((y_true == 1) & (y_pred == 1)); fp = np.sum((y_true == 0) & (y_pred == 1))
    return tp / (tp + fp) if (tp + fp) else 0.0


def f1_score_custom(y_true, y_pred):
    r, p = recall_score(y_true, y_pred), precision_score_custom(y_true, y_pred)
    return 2 * p * r / (p + r) if (p + r) else 0.0


def safe_roc_auc(y_true, scores):
    return roc_auc_score(y_true, scores) if len(np.unique(y_true)) > 1 else float("nan")


def safe_pr_auc(y_true, scores):
    if len(np.unique(y_true)) < 2:
        return float("nan")
    prec, rec, _ = precision_recall_curve(y_true, scores)
    return auc(rec, prec)


LOGAUC_LAMBDA = 0.001
RANDOM_LOGAUC = (1.0 - LOGAUC_LAMBDA) / (np.log(10) * np.log10(1.0 / LOGAUC_LAMBDA))


def build_log_auc_curve(ytrue, scores, min_frac=LOGAUC_LAMBDA, max_frac=1.0):
    ytrue = np.asarray(ytrue); scores = np.asarray(scores, dtype=float)
    n = len(ytrue); n_pos = int(np.sum(ytrue == 1)); n_neg = n - n_pos
    if n_pos == 0 or n_neg == 0:
        return np.array([min_frac, max_frac]), np.array([0.0, 1.0]), float("nan")
    order = np.argsort(scores)[::-1]
    y_sorted = ytrue[order]
    fpr = np.cumsum(y_sorted == 0) / n_neg
    tpr = np.cumsum(y_sorted == 1) / n_pos
    fpr_full = np.concatenate([[0.0], fpr]); tpr_full = np.concatenate([[0.0], tpr])
    tpr_at_min = float(np.interp(min_frac, fpr_full, tpr_full))
    tpr_at_max = float(np.interp(max_frac, fpr_full, tpr_full))
    mask = (fpr > min_frac) & (fpr < max_frac)
    x = np.concatenate([[min_frac], fpr[mask], [max_frac]])
    y = np.concatenate([[tpr_at_min], tpr[mask], [tpr_at_max]])
    norm = np.log10(max_frac / min_frac)
    _trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    raw_log_auc = float(_trapz(y, np.log10(x)) / norm)
    adjusted_log_auc = (raw_log_auc - RANDOM_LOGAUC) * 100.0
    return x, y, adjusted_log_auc


def enrichment_factor(y_true, scores, fraction):
    y_true = np.asarray(y_true); scores = np.asarray(scores, dtype=float)
    n = len(y_true); n_pos = int(np.sum(y_true == 1))
    if n == 0 or n_pos == 0:
        return float("nan")
    k = max(1, int(np.ceil(fraction * n)))
    order = np.argsort(scores)[::-1]
    top_hits = int(np.sum(y_true[order[:k]] == 1))
    return float((top_hits / k) / (n_pos / n))


def bedroc_score(y_true, scores, alpha=20.0):
    y_true = np.asarray(y_true); scores = np.asarray(scores, dtype=float)
    n = len(y_true); n_pos = int(np.sum(y_true == 1))
    if n == 0 or n_pos == 0 or n_pos == n:
        return float("nan")
    Ra = n_pos / n
    order = np.argsort(scores)[::-1]
    ranks = np.where(y_true[order] == 1)[0] + 1
    rie_num = np.sum(np.exp(-alpha * ranks / n))
    rie_den = (n_pos / n) * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = rie_num / rie_den
    factor = Ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * Ra))
    bedroc = rie * factor + 1.0 / (1.0 - np.exp(alpha * (1 - Ra)))
    return float(bedroc)


def balanced_accuracy_custom(y_true, y_pred):
    return 0.5 * (recall_score(y_true, y_pred) + specificity_score(y_true, y_pred))


def cohen_kappa_custom(y_true, y_pred):
    y_true = np.asarray(y_true); y_pred = np.asarray(y_pred)
    n = len(y_true)
    if n == 0:
        return float("nan")
    po = np.mean(y_true == y_pred)
    p_yes = np.mean(y_true == 1) * np.mean(y_pred == 1)
    p_no = np.mean(y_true == 0) * np.mean(y_pred == 0)
    pe = p_yes + p_no
    return float((po - pe) / (1 - pe)) if (1 - pe) else 0.0


def brier_score(y_true, probs):
    y_true = np.asarray(y_true, dtype=float); probs = np.asarray(probs, dtype=float)
    if len(y_true) == 0:
        return float("nan")
    return float(np.mean((probs - y_true) ** 2))


def classification_metrics(y_true, scores, threshold):
    y_pred = (scores >= threshold).astype(int)
    return {
        "mcc": matthews_corrcoef(y_true, y_pred),
        "recall": recall_score(y_true, y_pred),
        "specificity": specificity_score(y_true, y_pred),
        "precision": precision_score_custom(y_true, y_pred),
        "f1": f1_score_custom(y_true, y_pred),
        "balanced_acc": balanced_accuracy_custom(y_true, y_pred),
        "kappa": cohen_kappa_custom(y_true, y_pred),
        "roc_auc": safe_roc_auc(y_true, scores),
        "pr_auc": safe_pr_auc(y_true, scores),
        "log_auc": build_log_auc_curve(y_true, scores)[2],
        "ef1": enrichment_factor(y_true, scores, 0.01),
        "ef5": enrichment_factor(y_true, scores, 0.05),
        "bedroc": bedroc_score(y_true, scores, alpha=20.0),
    }


def optimize_mcc_threshold(y_true: np.ndarray, scores: np.ndarray,
                            threshold_grid: np.ndarray) -> dict:
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    best = {"threshold": float(threshold_grid[0]), "mcc": -2.0, "recall": float("nan")}
    if len(y_true) == 0 or len(np.unique(y_true)) < 2:
        return best
    for thr in threshold_grid:
        y_pred = (scores >= thr).astype(int)
        m = matthews_corrcoef(y_true, y_pred)
        if m > best["mcc"]:
            best = {"threshold": float(thr), "mcc": float(m),
                    "recall": float(recall_score(y_true, y_pred))}
    return best


def bootstrap_cis(y_true, scores, threshold, n_boot=1000, alpha=0.05, seed=GLOBAL_SEED):
    rng = np.random.default_rng(seed)
    y_true = np.asarray(y_true); scores = np.asarray(scores, dtype=float)
    n = len(y_true)
    score_keys = ["mcc", "recall", "specificity", "precision", "f1",
                  "balanced_acc", "kappa", "roc_auc", "pr_auc", "log_auc",
                  "ef1", "ef5", "bedroc"]
    raw = {k: [] for k in score_keys}
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt = y_true[idx]
        if len(np.unique(yt)) < 2:
            continue
        m = classification_metrics(yt, scores[idx], threshold)
        for k in score_keys:
            raw[k].append(m[k])
    lo, hi = 100 * alpha / 2, 100 * (1 - alpha / 2)
    ci = {}
    for k, vals in raw.items():
        if vals:
            v = np.array(vals)
            ci[k] = {"mean": float(np.mean(v)), "lower": float(np.percentile(v, lo)),
                      "upper": float(np.percentile(v, hi)), "n_valid": len(v)}
        else:
            ci[k] = {"mean": np.nan, "lower": np.nan, "upper": np.nan, "n_valid": 0}
    return ci


# ==================================================
# STAGE 1: FEATURE MATRIX (RDKit + electrophilicity)
# ==================================================

def stage1_build_features(df: pd.DataFrame) -> Tuple[pd.DataFrame, np.ndarray, Dict[int, dict]]:
    print("\n" + "=" * 70 + "\nSTAGE 1: FEATURE MATRIX (RDKit + electrophilicity)\n" + "=" * 70)
    backend = "xtb_gfn2" if (_xtb_available() and ELECTRO_BACKEND in ("auto", "xtb_gfn2")) else "rdkit_surrogate"
    print(f"  Electrophilicity backend: {backend}")
    X, valid_indices, electro_records = build_full_feature_matrix(df)
    df_valid = df.loc[valid_indices].reset_index(drop=True)
    print("[OK] Stage 1 complete")
    return df_valid, X, electro_records


# ==================================================
# STAGE 2: ML MODEL TRAINING (soft-voting ensemble)
# ==================================================

def train_ml_model_with_smote(X, y, ensemble_config=None, use_smote=True,
                               smote_ratio=SMOTE_RATIO, verbose=True):
    if ensemble_config is None:
        ensemble_config = ENSEMBLE_CONFIG
    if verbose:
        names = ", ".join(ensemble_config.keys())
        print(f"\nTRAINING ENSEMBLE MODEL ({'with SMOTE' if use_smote else 'standard'}) "
              f"[{names}] | n={len(X)}")
    n_stable = int(np.sum(y == 1)); n_unstable = int(np.sum(y == 0))
    if verbose:
        print(f"  Original: Stable={n_stable}, Unstable={n_unstable}")

    if use_smote and SMOTE_AVAILABLE and n_stable > 1:
        smote = SMOTE(sampling_strategy=smote_ratio, random_state=GLOBAL_SEED,
                       k_neighbors=min(5, n_stable - 1))
        X_res, y_res = smote.fit_resample(X, y)
        if verbose:
            print(f"  After SMOTE: Stable={int(np.sum(y_res == 1))}, "
                  f"Unstable={int(np.sum(y_res == 0))}")
    else:
        X_res, y_res = X, y

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_res)
    ml_model = build_ensemble_model(ensemble_config)
    ml_model.fit(X_scaled, y_res)

    if verbose:
        cv = cross_val_score(ml_model, scaler.transform(X), y,
                              cv=CV_FOLDS, scoring=make_scorer(matthews_corrcoef))
        print(f"  {CV_FOLDS}-fold CV MCC: {cv.mean():.3f} +/- {cv.std():.3f}")
    return ml_model, scaler


def stage2_train_ml_model(df_valid: pd.DataFrame, X: np.ndarray) -> dict:
    print("\n" + "=" * 70 + "\nSTAGE 2: ML MODEL TRAINING & TRAIN/TEST SPLIT\n" + "=" * 70)
    y = df_valid["Label"].values

    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
        X, y, np.arange(len(df_valid)),
        test_size=TEST_SPLIT, stratify=y, random_state=GLOBAL_SEED)
    df_train = df_valid.iloc[idx_train].reset_index(drop=True)
    df_test = df_valid.iloc[idx_test].reset_index(drop=True)
    train_ligand_ids = set(df_train["Ligand_ID"].values)
    test_ligand_ids = set(df_test["Ligand_ID"].values)
    print(f"  Train: {len(X_train)} | Test: {len(X_test)} (held out)")

    ml_model, scaler_ml = train_ml_model_with_smote(
        X_train, y_train, ENSEMBLE_CONFIG, use_smote=USE_SMOTE, smote_ratio=SMOTE_RATIO)

    oof_train_proba = {}
    train_ids_arr = df_train["Ligand_ID"].values
    skf = StratifiedKFold(n_splits=CV_FOLDS, shuffle=True, random_state=GLOBAL_SEED)
    for tr_idx, te_idx in skf.split(X_train, y_train):
        fold_model, fold_scaler = train_ml_model_with_smote(
            X_train[tr_idx], y_train[tr_idx], ENSEMBLE_CONFIG,
            use_smote=USE_SMOTE, smote_ratio=SMOTE_RATIO, verbose=False)
        proba = fold_model.predict_proba(fold_scaler.transform(X_train[te_idx]))[:, 1]
        for j, p_val in zip(te_idx, proba):
            oof_train_proba[int(train_ids_arr[j])] = float(p_val)

    ml_artifacts = {
        "ml_model": ml_model, "scaler_ml": scaler_ml, "features": ALL_FEATURE_NAMES,
        "ensemble_config": ENSEMBLE_CONFIG,
        "df_train": df_train, "df_test": df_test,
        "X_train": X_train, "X_test": X_test, "y_train": y_train, "y_test": y_test,
        "train_ligand_ids": train_ligand_ids, "test_ligand_ids": test_ligand_ids,
        "oof_train_proba": oof_train_proba,
        "id_to_row": {int(lid): i for i, lid in enumerate(df_valid["Ligand_ID"].values)},
        "X_full": X, "df_valid": df_valid,
    }
    print(f"[OK] Out-of-fold ML predictions for {len(oof_train_proba)} train ligands")
    print("[OK] Stage 2 complete")
    return ml_artifacts


def ml_probs_all_ligands(ml_artifacts) -> Dict[int, float]:
    ml_model, scaler_ml = ml_artifacts["ml_model"], ml_artifacts["scaler_ml"]
    id_to_row, X_full = ml_artifacts["id_to_row"], ml_artifacts["X_full"]
    Xs = scaler_ml.transform(X_full)
    proba = ml_model.predict_proba(Xs)[:, 1]
    return {lid: float(proba[row]) for lid, row in id_to_row.items()}


# ===========================================
# STAGE 3: BENCHMARKING, CONTROLS & REPORTING
# ===========================================

def evaluate_ml_model(ml_probs, ground_truth, ligand_ids, threshold,
                       set_name="Set", compute_ci=False, n_boot=1000) -> dict:
    lids = [l for l in ligand_ids if l in ml_probs and l in ground_truth]
    y_true = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids])
    scores = np.array([ml_probs[l] for l in lids])
    metrics = classification_metrics(y_true, scores, threshold)
    metrics["n_samples"] = len(y_true)
    print(f"  [{set_name}] n={len(y_true)} MCC={metrics['mcc']:.3f} "
          f"recall={metrics['recall']:.3f} ROC-AUC={metrics['roc_auc']:.3f}")
    if compute_ci:
        metrics["bootstrap_cis"] = bootstrap_cis(y_true, scores, threshold, n_boot=n_boot)
    return metrics


def y_randomisation_control(ml_artifacts, ml_probs, ground_truth,
                             n_runs=Y_RANDOM_RUNS, seed=GLOBAL_SEED) -> dict:
    print("\n" + "=" * 70 + "\nY-RANDOMISATION CONTROL (label shuffling)\n" + "=" * 70)
    rng = np.random.default_rng(seed)
    train_ids = ml_artifacts["train_ligand_ids"]
    test_ids = ml_artifacts["test_ligand_ids"]
    X_train, y_train = ml_artifacts["X_train"], ml_artifacts["y_train"]
    X_test = ml_artifacts["X_test"]
    lids_te = [l for l in ml_artifacts["df_test"]["Ligand_ID"].values if l in ground_truth]
    y_te = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids_te])

    null_mcc = []
    for _ in range(n_runs):
        y_perm = rng.permutation(y_train)
        if len(np.unique(y_perm)) < 2:
            continue
        model, scaler = train_ml_model_with_smote(
            X_train, y_perm, ENSEMBLE_CONFIG, use_smote=USE_SMOTE,
            smote_ratio=SMOTE_RATIO, verbose=False)
        proba_tr = model.predict_proba(scaler.transform(X_train))[:, 1]
        thr = optimize_mcc_threshold(y_perm, proba_tr, THRESHOLD_GRID)["threshold"]
        proba_te = model.predict_proba(scaler.transform(X_test))[:, 1]
        y_pred_te = (proba_te >= thr).astype(int)
        if len(np.unique(y_te)) > 1:
            null_mcc.append(matthews_corrcoef(y_te, y_pred_te))
    null_mcc = np.array(null_mcc) if null_mcc else np.array([np.nan])

    true_scores = np.array([ml_probs[l] for l in lids_te])
    true_thr = ml_artifacts["ml_threshold"]
    true_pred = (true_scores >= true_thr).astype(int)
    true_mcc = (matthews_corrcoef(y_te, true_pred)
                if len(np.unique(y_te)) > 1 else float("nan"))
    p_emp = float((np.sum(null_mcc >= true_mcc) + 1) / (len(null_mcc) + 1))
    res = {"true_test_mcc": float(true_mcc),
           "null_mean": float(np.nanmean(null_mcc)),
           "null_std": float(np.nanstd(null_mcc)),
           "null_max": float(np.nanmax(null_mcc)),
           "p_empirical": p_emp, "n_runs": int(len(null_mcc)),
           "null_mcc": null_mcc.tolist()}
    print(f"  True test MCC={true_mcc:.3f} | null mean={res['null_mean']:.3f} "
          f"+/- {res['null_std']:.3f} (max={res['null_max']:.3f}) | p={p_emp:.4f}")
    return res


def applicability_domain(ml_artifacts) -> dict:
    print("\n" + "=" * 70 + "\nAPPLICABILITY DOMAIN (descriptor-space leverage)\n" + "=" * 70)
    scaler = ml_artifacts["scaler_ml"]
    Xtr = scaler.transform(ml_artifacts["X_train"])
    Xte = scaler.transform(ml_artifacts["X_test"])
    n, p = Xtr.shape
    gram = Xtr.T @ Xtr + 1e-8 * np.eye(p)
    gram_inv = np.linalg.inv(gram)
    h_train = np.einsum("ij,jk,ik->i", Xtr, gram_inv, Xtr)
    h_test = np.einsum("ij,jk,ik->i", Xte, gram_inv, Xte)
    h_star = 3.0 * p / n
    n_out_tr = int(np.sum(h_train > h_star))
    n_out_te = int(np.sum(h_test > h_star))
    res = {"h_star": float(h_star), "p": int(p), "n_train": int(n),
           "h_train": h_train.tolist(), "h_test": h_test.tolist(),
           "n_outside_train": n_out_tr, "n_outside_test": n_out_te,
           "frac_outside_test": float(n_out_te / len(h_test)) if len(h_test) else float("nan")}
    print(f"  Warning leverage h*={h_star:.4f} | test outside AD: "
          f"{n_out_te}/{len(h_test)} ({100 * res['frac_outside_test']:.1f}%)")
    return res


def compute_permutation_importance(ml_model, scaler, X, y, feature_names,
                                    n_repeats=PERM_IMPORTANCE_REPEATS, seed=GLOBAL_SEED):
    print("\n" + "=" * 70 + "\nPERMUTATION FEATURE IMPORTANCE (ensemble, model-agnostic)\n" + "=" * 70)
    Xs = scaler.transform(X)
    result = permutation_importance(ml_model, Xs, y, scoring=make_scorer(matthews_corrcoef),
                                      n_repeats=n_repeats, random_state=seed, n_jobs=-1)
    out = {name: {"mean": float(m), "std": float(s)}
           for name, m, s in zip(feature_names, result.importances_mean, result.importances_std)}
    for name, vals in sorted(out.items(), key=lambda x: x[1]["mean"], reverse=True):
        print(f"  {name:>20s}: dMCC={vals['mean']:.4f} +/- {vals['std']:.4f}")
    return out


# -----
# Plots
# -----

def plot_roc_curves(ml_probs, ground_truth, train_ids, test_ids, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, (ids, name) in zip(axes, [(train_ids, "Training Set"), (test_ids, "Test Set")]):
        lids = [l for l in ids if l in ml_probs and l in ground_truth]
        y_true = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids])
        scores = np.array([ml_probs[l] for l in lids])
        if len(np.unique(y_true)) < 2:
            ax.set_title(f"ROC - {name} (insufficient classes)"); continue
        fpr, tpr, _ = roc_curve(y_true, scores)
        ax.plot(fpr, tpr, "g-", lw=2, label=f"Ensemble model (AUC={auc(fpr, tpr):.3f})")
        ax.plot([0, 1], [0, 1], "k--", lw=1.5, label="Random")
        ax.set_xlabel("False Positive Rate"); ax.set_ylabel("True Positive Rate")
        ax.set_title(f"ROC Curve - {name}", fontweight="bold")
        ax.legend(loc="lower right"); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "roc_curves_train_test.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] ROC curves saved: {out}")


def plot_pr_curves(ml_probs, ground_truth, train_ids, test_ids, output_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ax, (ids, name) in zip(axes, [(train_ids, "Training Set"), (test_ids, "Test Set")]):
        lids = [l for l in ids if l in ml_probs and l in ground_truth]
        y_true = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids])
        scores = np.array([ml_probs[l] for l in lids])
        if len(np.unique(y_true)) < 2:
            ax.set_title(f"PR - {name} (insufficient classes)"); continue
        prec, rec, _ = precision_recall_curve(y_true, scores)
        ax.plot(rec, prec, "g-", lw=2, label=f"Ensemble model (PR-AUC={auc(rec, prec):.3f})")
        ax.axhline(y=y_true.mean(), color="k", linestyle="--", lw=1.5, label="Random")
        ax.set_xlabel("Recall"); ax.set_ylabel("Precision")
        ax.set_title(f"PR Curve - {name}", fontweight="bold")
        ax.legend(loc="upper right"); ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "pr_curves_train_test.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] PR curves saved: {out}")


def plot_electrophilicity_distribution(electro_records, ground_truth, output_dir):
    stable_w, unstable_w = [], []
    for lid, r in electro_records.items():
        if lid not in ground_truth or not r["has_reactive"] or not np.isfinite(r["omega_local_warhead"]):
            continue
        (stable_w if ground_truth[lid] == "stable" else unstable_w).append(r["omega_local_warhead"])
    plt.figure(figsize=(8, 5))
    bins = np.linspace(0, max([0.5] + stable_w + unstable_w), 25)
    plt.hist(unstable_w, bins=bins, alpha=0.6, label=f"Unstable (n={len(unstable_w)})", color="crimson")
    plt.hist(stable_w, bins=bins, alpha=0.6, label=f"Stable (n={len(stable_w)})", color="seagreen")
    plt.xlabel("Warhead-gated local electrophilicity, omega_local_warhead (eV)"); plt.ylabel("Count")
    plt.title("Electrophilicity Feature (Warhead-Gated) by Class", fontweight="bold")
    plt.legend(); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "electrophilicity_distribution.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Electrophilicity distribution plot saved: {out}")


def plot_calibration_curve(ml_probs, ground_truth, test_ids, output_dir, n_bins=10):
    lids = [l for l in test_ids if l in ml_probs and l in ground_truth]
    y_true = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids])
    probs = np.array([ml_probs[l] for l in lids])
    brier = brier_score(y_true, probs)
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    plt.figure(figsize=(6, 6))
    plt.plot([0, 1], [0, 1], "k--", lw=1.5, label="Perfectly calibrated")
    xs, ys = [], []
    for i in range(n_bins):
        m = (probs >= edges[i]) & (probs < edges[i + 1] if i < n_bins - 1
                                    else probs <= edges[i + 1])
        if np.sum(m) > 0:
            xs.append(float(np.mean(probs[m]))); ys.append(float(np.mean(y_true[m])))
    plt.plot(xs, ys, "gs-", lw=2, label=f"Ensemble model (Brier={brier:.3f})")
    plt.xlabel("Mean predicted P(stable)"); plt.ylabel("Observed fraction stable")
    plt.title("Calibration (Test Set)", fontweight="bold")
    plt.xlim(0, 1); plt.ylim(0, 1); plt.legend(loc="upper left"); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "ml_calibration.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Calibration curve saved: {out}")
    return {"brier_ml": float(brier)}


def plot_confusion_matrix(ml_probs, ground_truth, test_ids, threshold, output_dir):
    lids = [l for l in test_ids if l in ml_probs and l in ground_truth]
    y_true = np.array([1 if ground_truth[l] == "stable" else 0 for l in lids])
    y_pred = np.array([int(ml_probs[l] >= threshold) for l in lids])
    fig, ax = plt.subplots(figsize=(5, 5))
    labels = ["Unstable", "Stable"]
    cm = np.zeros((2, 2), dtype=int)
    for t, pp in zip(y_true, y_pred):
        cm[t, pp] += 1
    ax.imshow(cm, cmap="Blues")
    ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
    ax.set_xticklabels(labels); ax.set_yticklabels(labels)
    ax.set_xlabel("Predicted"); ax.set_ylabel("True")
    ax.set_title("Ensemble Model (Test Set)", fontweight="bold")
    for i in range(2):
        for j in range(2):
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    color="white" if cm[i, j] > cm.max() / 2 else "black",
                    fontsize=14, fontweight="bold")
    plt.tight_layout()
    out = os.path.join(output_dir, "confusion_matrix_test.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Confusion matrix saved: {out}")


def plot_applicability_domain(ad, output_dir):
    plt.figure(figsize=(9, 5))
    h_tr = np.array(ad["h_train"]); h_te = np.array(ad["h_test"])
    plt.scatter(np.arange(len(h_tr)), h_tr, s=12, c="steelblue", alpha=0.6, label="Train")
    plt.scatter(np.arange(len(h_tr), len(h_tr) + len(h_te)), h_te, s=18,
                c="crimson", alpha=0.7, label="Test")
    plt.axhline(ad["h_star"], color="k", ls="--", lw=1.5,
                label=f"h* = {ad['h_star']:.3f}")
    plt.xlabel("Ligand index"); plt.ylabel("Leverage (h)")
    plt.title("Applicability Domain (descriptor-space leverage)", fontweight="bold")
    plt.legend(); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "applicability_domain.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Applicability-domain plot saved: {out}")


def plot_y_randomisation(yrand, output_dir):
    plt.figure(figsize=(8, 5))
    null = np.array(yrand["null_mcc"])
    plt.hist(null, bins=20, color="lightgray", edgecolor="k", alpha=0.8,
              label=f"Permuted-label null (n={yrand['n_runs']})")
    plt.axvline(yrand["true_test_mcc"], color="crimson", lw=2.5,
                label=f"True model (MCC={yrand['true_test_mcc']:.3f}, p={yrand['p_empirical']:.3f})")
    plt.xlabel("Test-set MCC"); plt.ylabel("Count")
    plt.title("Y-Randomisation Control", fontweight="bold")
    plt.legend(); plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "y_randomisation.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Y-randomisation plot saved: {out}")


def plot_permutation_importance(importances, output_dir):
    names = list(importances.keys())
    means = [importances[n]["mean"] for n in names]
    stds = [importances[n]["std"] for n in names]
    order = np.argsort(means)
    plt.figure(figsize=(8, 5))
    plt.barh([names[i] for i in order], [means[i] for i in order],
              xerr=[stds[i] for i in order], color="seagreen")
    plt.xlabel("Decrease in MCC when feature is permuted")
    plt.title("Permutation Feature Importance (Ensemble)", fontweight="bold")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(output_dir, "permutation_importance.png")
    plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
    print(f"[OK] Permutation importance plot saved: {out}")


def run_shap_analysis_ml(ml_model, scaler_ml, X_train, X_test, feature_names, output_dir,
                          background_size=SHAP_BACKGROUND_SIZE,
                          sample_size=SHAP_SAMPLE_SIZE, nsamples=SHAP_NSAMPLES,
                          seed=GLOBAL_SEED):
    # Model-agnostic SHAP via KernelExplainer
    print("RUNNING SHAP ANALYSIS FOR ML MODEL (KernelExplainer, ensemble-compatible)")
    rng = np.random.RandomState(seed)
    Xtr_s = scaler_ml.transform(X_train)
    n_bg = min(background_size, len(Xtr_s))
    try:
        background = shap.kmeans(Xtr_s, n_bg)
    except Exception:
        bg_idx = rng.choice(len(Xtr_s), size=n_bg, replace=False)
        background = Xtr_s[bg_idx]
    explainer = shap.KernelExplainer(ml_model.predict_proba, background)

    for tag, X in [("train", X_train), ("test", X_test)]:
        Xs = scaler_ml.transform(X)
        n_sample = min(sample_size, len(Xs))
        sample_idx = rng.choice(len(Xs), size=n_sample, replace=False)
        Xs_sample = Xs[sample_idx]
        raw = explainer.shap_values(Xs_sample, nsamples=nsamples)
        sv = raw[1] if isinstance(raw, list) else (raw[:, :, 1] if getattr(raw, "ndim", 0) == 3 else raw)
        plt.figure()
        shap.summary_plot(sv, Xs_sample, feature_names=feature_names, show=False)  # type: ignore
        plt.tight_layout()
        out = os.path.join(output_dir, f"shap_summary_ml_{tag}.png")
        plt.savefig(out, dpi=150, bbox_inches="tight"); plt.close()
        print(f"[OK] SHAP summary ({tag}) saved: {out}")


def save_final_predictions(ml_probs, ground_truth, train_ids, test_ids, threshold, output_dir):
    out = os.path.join(output_dir, "final_predictions.csv")
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ligand_ID", "Set", "Ground_Truth", "ML_Proba", "ML_Pred", "ML_Correct"])
        for lid in sorted(ml_probs):
            if lid not in ground_truth:
                continue
            s = "train" if lid in train_ids else ("test" if lid in test_ids else "unknown")
            gt = 1 if ground_truth[lid] == "stable" else 0
            proba = ml_probs[lid]
            pred = int(proba >= threshold)
            w.writerow([lid, s, gt, f"{proba:.4f}", pred, int(pred == gt)])
    print(f"[OK] Final predictions saved: {out}")


def save_complete_model(ml_artifacts, output_dir):
    model = {
        "ml_model": ml_artifacts["ml_model"], "scaler_ml": ml_artifacts["scaler_ml"],
        "features": ml_artifacts["features"],
        "ensemble_config": ml_artifacts["ensemble_config"],
        "electro_feature_names": list(ELECTRO_FEATURE_NAMES),
        "electro_backend": "xtb_gfn2" if (_xtb_available() and ELECTRO_BACKEND in ("auto", "xtb_gfn2")) else "rdkit_surrogate",
        "ml_threshold": ml_artifacts["ml_threshold"],
        "training_date": datetime.now().isoformat(),
    }
    out = os.path.join(output_dir, "complete_ml_model.pkl")
    with open(out, "wb") as f:
        pickle.dump(model, f)
    print(f"[OK] Complete model saved: {out}")


def generate_performance_report(ml_artifacts, train_metrics, test_metrics,
                                 calibration, ad, yrand, perm_importance, output_dir):
    rep = os.path.join(output_dir, "performance_report.txt")
    score_keys = ["mcc", "recall", "specificity", "precision", "f1",
                  "balanced_acc", "kappa", "roc_auc", "pr_auc", "log_auc",
                  "ef1", "ef5", "bedroc"]
    hdr = {"mcc": "MCC", "recall": "RECALL", "specificity": "SPEC",
           "precision": "PREC", "f1": "F1", "balanced_acc": "BAL_ACC",
           "kappa": "KAPPA", "roc_auc": "ROC_AUC", "pr_auc": "PR_AUC",
           "log_auc": "LOGAUC", "ef1": "EF1%", "ef5": "EF5%", "bedroc": "BEDROC"}

    def line(c="="):
        return c * 78 + "\n"

    with open(rep, "w", encoding="utf-8") as f:
        f.write(line() + "ML-ONLY MODEL PERFORMANCE REPORT\n" + line())
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Train samples: {len(ml_artifacts['X_train'])} | "
                f"Test samples: {len(ml_artifacts['X_test'])} (held out until Stage 3)\n\n")

        f.write(line("-") + "1. MODEL CONFIGURATION\n" + line("-"))
        f.write("Model: soft-voting ensemble classifier (VotingClassifier, voting='soft')\n")
        f.write("Base estimators:\n")
        for name, spec in ml_artifacts["ensemble_config"].items():
            f.write(f"  - {name} ({spec['estimator'].__name__}), weight={spec['weight']}\n")
            for pname, pval in spec["params"].items():
                f.write(f"      {pname} = {pval}\n")
        f.write("Features: RDKit descriptors + electrophilicity descriptors\n")
        for name in ml_artifacts["features"]:
            f.write(f"  - {name}\n")
        f.write(f"Decision threshold (train-tuned, MCC-optimal): {ml_artifacts['ml_threshold']:.4f}\n\n")

        f.write(line("-") + "2. PERFORMANCE METRICS\n" + line("-"))
        for set_name, metrics in [("TRAINING SET", train_metrics), ("TEST SET", test_metrics)]:
            f.write(f"  {set_name} (n={metrics['n_samples']})\n")
            f.write("  " + "".join(f"{hdr[k]:>9}" for k in score_keys) + "\n")
            f.write("  " + "".join(f"{metrics[k]:>9.3f}" for k in score_keys) + "\n\n")

        if "bootstrap_cis" in test_metrics:
            f.write(line("-") + "3. BOOTSTRAP 95% CONFIDENCE INTERVALS (TEST SET, n=1000)\n" + line("-"))
            for k in score_keys:
                ci = test_metrics["bootstrap_cis"][k]
                f.write(f"  {hdr[k]:<10} {test_metrics[k]:>8.4f} "
                        f"[{ci['lower']:>7.4f}, {ci['upper']:>7.4f}]\n")
            f.write("\n")

        f.write(line("-") + "4. PROBABILITY CALIBRATION (TEST SET)\n" + line("-"))
        f.write(f"Ensemble model Brier score : {calibration['brier_ml']:.4f}\n")
        f.write("(0 = perfect, lower is better). Reliability curve: ml_calibration.png\n\n")

        f.write(line("-") + "5. APPLICABILITY DOMAIN (descriptor-space leverage)\n" + line("-"))
        f.write(f"Descriptors p={ad['p']} | train n={ad['n_train']} | "
                f"warning leverage h*={ad['h_star']:.4f}\n")
        f.write(f"Train ligands outside AD: {ad['n_outside_train']}\n")
        f.write(f"Test ligands outside AD : {ad['n_outside_test']}/"
                f"{len(ad['h_test'])} ({100 * ad['frac_outside_test']:.1f}%)\n")
        f.write("Williams leverage plot saved to applicability_domain.png\n\n")

        f.write(line("-") + "6. Y-RANDOMISATION CONTROL (label shuffling)\n" + line("-"))
        f.write(f"Permutations: {yrand['n_runs']}\n")
        f.write(f"True test MCC       : {yrand['true_test_mcc']:.4f}\n")
        f.write(f"Null MCC mean+/-sd  : {yrand['null_mean']:.4f} +/- {yrand['null_std']:.4f}\n")
        f.write(f"Null MCC max        : {yrand['null_max']:.4f}\n")
        f.write(f"Empirical p-value   : {yrand['p_empirical']:.4f}\n")
        f.write("  p << 0.05 indicates the model has learnt genuine signal, not chance.\n\n")

        f.write(line("-") + "7. PERMUTATION FEATURE IMPORTANCE (TEST SET, dMCC)\n" + line("-"))
        f.write("(Model-agnostic; replaces RF-only feature_importances_ since the\n"
                " ensemble mixes RF, SVM, and GB estimators.)\n")
        for name, vals in sorted(perm_importance.items(), key=lambda x: x[1]["mean"], reverse=True):
            f.write(f"  {name:<20} {vals['mean']:.4f} +/- {vals['std']:.4f}\n")
        f.write("Bar chart saved to permutation_importance.png\n\n")

        f.write(line("-") + "8. SHAP ANALYSIS\n" + line("-"))
        f.write("Model-agnostic SHAP via KernelExplainer (background="
                f"{SHAP_BACKGROUND_SIZE} pts, sample={SHAP_SAMPLE_SIZE} pts, "
                f"nsamples={SHAP_NSAMPLES}).\n")
        f.write("Plots saved to shap_summary_ml_train.png / shap_summary_ml_test.png\n")

    print(f"[OK] Performance report saved: {rep}")


def stage3_benchmark_and_report(ml_artifacts, electro_records) -> dict:
    print("\n" + "=" * 70 + "\nSTAGE 3: BENCHMARKING, CONTROLS & REPORTING\n" + "=" * 70)
    gt = load_ground_truth(TRAINING_DATA_FILE)
    train_ids, test_ids = ml_artifacts["train_ligand_ids"], ml_artifacts["test_ligand_ids"]

    ml_probs = ml_probs_all_ligands(ml_artifacts)
    lids_tr = [l for l in train_ids if l in ml_probs and l in gt]
    y_tr = np.array([1 if gt[l] == "stable" else 0 for l in lids_tr])
    scores_tr = np.array([ml_probs[l] for l in lids_tr])
    thr = optimize_mcc_threshold(y_tr, scores_tr, THRESHOLD_GRID)
    ml_artifacts["ml_threshold"] = thr["threshold"]
    print(f"[OK] ML decision threshold={thr['threshold']:.3f} | "
          f"train MCC={thr['mcc']:.3f} recall={thr['recall']:.3f}")

    train_metrics = evaluate_ml_model(ml_probs, gt, train_ids, thr["threshold"], "Training Set")
    test_metrics = evaluate_ml_model(ml_probs, gt, test_ids, thr["threshold"], "Test Set",
                                      compute_ci=True, n_boot=1000)
    ad = applicability_domain(ml_artifacts)
    yrand = y_randomisation_control(ml_artifacts, ml_probs, gt)

    perm_importance = compute_permutation_importance(
        ml_artifacts["ml_model"], ml_artifacts["scaler_ml"],
        ml_artifacts["X_test"], ml_artifacts["y_test"], ml_artifacts["features"])
    plot_permutation_importance(perm_importance, OUTPUT_DIR)

    run_shap_analysis_ml(ml_artifacts["ml_model"], ml_artifacts["scaler_ml"],
                          ml_artifacts["X_train"], ml_artifacts["X_test"],
                          ml_artifacts["features"], OUTPUT_DIR)
    plot_roc_curves(ml_probs, gt, train_ids, test_ids, OUTPUT_DIR)
    plot_pr_curves(ml_probs, gt, train_ids, test_ids, OUTPUT_DIR)
    plot_electrophilicity_distribution(electro_records, gt, OUTPUT_DIR)
    calibration = plot_calibration_curve(ml_probs, gt, test_ids, OUTPUT_DIR)
    plot_confusion_matrix(ml_probs, gt, test_ids, thr["threshold"], OUTPUT_DIR)
    plot_applicability_domain(ad, OUTPUT_DIR)
    plot_y_randomisation(yrand, OUTPUT_DIR)
    save_final_predictions(ml_probs, gt, train_ids, test_ids, thr["threshold"], OUTPUT_DIR)
    save_complete_model(ml_artifacts, OUTPUT_DIR)
    generate_performance_report(ml_artifacts, train_metrics, test_metrics,
                                 calibration, ad, yrand, perm_importance, OUTPUT_DIR)

    print(f"[OK] Stage 3 complete - outputs in {OUTPUT_DIR}")
    return {"train_metrics": train_metrics, "test_metrics": test_metrics,
            "calibration": calibration, "applicability_domain": ad,
            "y_randomisation": yrand, "permutation_importance": perm_importance}


def main():
    start = time.perf_counter()
    print("=" * 70 + "\nML-ONLY PIPELINE (RDKit + electrophilicity features, ensemble model)\n" + "=" * 70)

    df = load_training_data(TRAINING_DATA_FILE)

    df_valid, X, electro_records = stage1_build_features(df)

    ml_artifacts = stage2_train_ml_model(df_valid, X)

    output_results = stage3_benchmark_and_report(ml_artifacts, electro_records)

    print("\n" + "=" * 70 + "\nPIPELINE COMPLETE\n" + "=" * 70)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Outputs in: {OUTPUT_DIR}")
    print(f"Elapsed: {time.perf_counter() - start:.2f} s")
    return {"ml_artifacts": ml_artifacts, "output_results": output_results}


if __name__ == "__main__":
    try:
        main()
        print("\n[OK] Pipeline executed successfully!")
        sys.exit(0)
    except KeyboardInterrupt:
        print("\n[!] Interrupted by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\n[ERROR] Pipeline failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
