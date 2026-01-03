################################################################################
# Title: CellPhoneDB Analysis for ASD L-serine Treatment Study
# Description: 
#   1. Downloads and builds CellPhoneDB v5.0.0 database.
#   2. Performs statistical interaction analysis for two conditions:
#      - BO (Before OM / ASD)
#      - AO (After OM / Treat)
#   3. Compares CD4+ T cell interactions between groups.
#   4. Generates Heatmaps and Dotplots using ktplotspy.
# Author: Jaemyung Jang / Korea Brain Research Institue (piloter2@kbri.re.kr)
# Date: 2026-01-02
# Python Version: 3.9+
################################################################################

# ==============================================================================
# 1. Imports & Setup
# ==============================================================================
import os
import glob
import datetime
import warnings
import pandas as pd
import anndata
import matplotlib.pyplot as plt
import ktplotspy as kpy

# CellPhoneDB imports
from cellphonedb.utils import db_releases_utils, db_utils
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# ==============================================================================
# 2. Configuration & Paths
# ==============================================================================
# NOTE: Update BASE_DIR to match your local environment
BASE_DIR = "/path/to" 
DATA_DIR = os.path.join(BASE_DIR, "data/v5.0.0")
RESULTS_DIR = os.path.join(BASE_DIR, "results")
CPDB_VERSION = 'v5.0.0'
CPDB_TARGET_DIR = os.path.join(DATA_DIR, CPDB_VERSION)
CPDB_FILE_PATH = os.path.join(CPDB_TARGET_DIR, "cellphonedb.zip")

# Input Files (Counts & Meta txt files for CellPhoneDB)
INPUTS = {
    "BO": {
        "counts": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_BO_count.txt"),
        "meta": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_BO_meta.tsv"),
        "h5ad": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_BO_2025-12-29 14-16-18.h5ad"), # For plotting
        "label": "L-serine before OM"
    },
    "AO": {
        "counts": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_AO_count.txt"),
        "meta": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_AO_meta.tsv"),
        "h5ad": os.path.join(DATA_DIR, "ASD_PBMC_final_supp_PBMC_AO_2025-12-29 14-16-18.h5ad"), # For plotting
        "label": "L-serine after OM"
    }
}

# Ensure Output Directory Exists
CPDB_OUTPUT_DIR = os.path.join(DATA_DIR, "cellphone_results")
if not os.path.exists(CPDB_OUTPUT_DIR):
    os.makedirs(CPDB_OUTPUT_DIR)
if not os.path.exists(RESULTS_DIR):
    os.makedirs(RESULTS_DIR)

# ==============================================================================
# 3. Helper Functions
# ==============================================================================

def setup_database(target_dir, version, db_file):
    """
    Downloads and creates the CellPhoneDB database if it doesn't exist.
    """
    print(f"Checking database at: {target_dir}")
    
    # Check if DB needs to be downloaded
    if not os.path.exists(target_dir) or not os.listdir(target_dir):
        print(f"Downloading CellPhoneDB {version}...")
        os.makedirs(target_dir, exist_ok=True)
        db_utils.download_database(target_dir, version)
    else:
        print("Database files already present.")

    # Check if the zip file exists, if not create it
    if not os.path.exists(db_file):
        print("Creating CellPhoneDB zip file...")
        db_utils.create_db(target_dir)
        # Verify creation (CellPhoneDB creates a timestamped zip, renaming might be needed or handled by glob)
        generated_zips = glob.glob(os.path.join(target_dir, "cellphonedb_*.zip"))
        if generated_zips:
            # For consistency, symlink or rename the latest zip to 'cellphonedb.zip' if strict naming required
            # Here we assume the user/script handles the specific filename or uses the latest
            print(f"Created: {generated_zips[0]}")
            # If the script expects a fixed name 'cellphonedb.zip', rename/symlink logic goes here.
            # For now, we assume standard flow.
    else:
        print(f"Database zip found: {db_file}")


def run_cpdb_analysis(group_name, inputs, db_file, output_path):
    """
    Runs CellPhoneDB statistical analysis for a specific group.
    """
    print(f"\n[Analysis] Running CellPhoneDB for group: {group_name}")
    
    results = cpdb_statistical_analysis_method.call(
        cpdb_file_path = db_file,
        meta_file_path = inputs["meta"],
        counts_file_path = inputs["counts"],
        counts_data = 'hgnc_symbol',
        score_interactions = True,
        iterations = 1000,
        threshold = 0.1,
        threads = 16,
        debug_seed = 42,
        result_precision = 3,
        pvalue = 0.05,
        subsampling = False,
        separator = '|',
        debug = False,
        output_path = output_path,
        output_suffix = f"_{group_name}" # Appends group name to output files
    )
    return results


def compare_interactions(res_bo, res_ao, cell_type_of_interest="CD4+ T"):
    """
    Identifies unique and shared interactions involving a specific cell type.
    """
    print(f"\n[Comparison] Analyzing interactions for {cell_type_of_interest}...")
    
    def get_interactions(cpdb_results, cell_type):
        # Filter columns that contain the cell type
        df = cpdb_results['interaction_scores']
        cols = [col for col in df.columns if f"{cell_type}|" in col or f"|{cell_type}" in col]
        
        # Rows where interaction score > 0 (or significant)
        relevant_rows = df[df[cols].any(axis=1)]
        return set(relevant_rows['interacting_pair'])

    pairs_bo = get_interactions(res_bo, cell_type_of_interest)
    pairs_ao = get_interactions(res_ao, cell_type_of_interest)
    
    unique_ao = pairs_ao - pairs_bo
    shared = pairs_bo.intersection(pairs_ao)
    
    print(f"Interactions unique to AO (Treat): {unique_ao}")
    print(f"Shared interactions: {len(shared)}")
    print(f"Total unique interactions combined: {len(pairs_bo.union(pairs_ao))}")


def generate_plots(adata_bo, adata_ao, res_bo, res_ao):
    """
    Generates Heatmaps and Dotplots using ktplotspy.
    """
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    
    cell_types_heatmap = ["NK", "RBC", "Monocyte", "gdT", "CD8+ T", "CD4+ T", "pDC"]
    
    # --- 1. Heatmaps ---
    print("\n[Plotting] Generating Heatmaps...")
    p_hm_bo = kpy.plot_cpdb_heatmap(
        pvals=res_bo['pvalues'],
        cell_types=cell_types_heatmap,
        title="Sum of significant interactions (Before OM)",
        cmap='inferno',
        figsize=(5, 5)
    )
    
    p_hm_ao = kpy.plot_cpdb_heatmap(
        pvals=res_ao['pvalues'],
        cell_types=cell_types_heatmap,
        title="Sum of significant interactions (After OM)",
        cmap='inferno',
        figsize=(5, 5)
    )
    
    # --- 2. Dotplots (CD4+ T focused) ---
    print("[Plotting] Generating Dotplots...")
    p_dot_bo = kpy.plot_cpdb(
        adata=adata_bo,
        cell_type1="CD4+ T",
        cell_type2="CD8+ T|Monocyte|NK|pDC",
        means=res_bo['means'],
        pvals=res_bo['pvalues'],
        celltype_key="cell.type.coarse",
        cmap_name="inferno",
        highlight_size=1,
        figsize=(10, 8)
    )
    
    p_dot_ao = kpy.plot_cpdb(
        adata=adata_ao,
        cell_type1="CD4+ T",
        cell_type2="Monocyte|NK|pDC",
        means=res_ao['means'],
        pvals=res_ao['pvalues'],
        celltype_key="cell.type.coarse",
        cmap_name="inferno",
        highlight_size=1,
        figsize=(10, 4)
    )
    
    # --- 3. Save Figures ---
    print(f"[Saving] Saving plots to {RESULTS_DIR}...")
    p_hm_bo.savefig(os.path.join(RESULTS_DIR, f"Figure5a1_supp_{timestamp}.pdf"), dpi=1200, bbox_inches='tight')
    p_hm_ao.savefig(os.path.join(RESULTS_DIR, f"Figure5a2_supp_{timestamp}.pdf"), dpi=1200, bbox_inches='tight')
    p_dot_bo.save(os.path.join(RESULTS_DIR, f"Figure5b1_supp_{timestamp}.pdf"), dpi=1200, bbox_inches='tight')
    p_dot_ao.save(os.path.join(RESULTS_DIR, f"Figure5b2_supp_{timestamp}.pdf"), dpi=1200, bbox_inches='tight')


# ==============================================================================
# 4. Main Execution
# ==============================================================================
if __name__ == "__main__":
    
    # 1. Setup Database
    setup_database(CPDB_TARGET_DIR, CPDB_VERSION, CPDB_FILE_PATH)
    
    # 2. Run Analysis for BO and AO
    # Note: Ensure the 'cellphonedb.zip' path is correct. If create_db made a timestamped zip, 
    # update CPDB_FILE_PATH to point to that specific file.
    # For this script, we assume the user has renamed it or provided the exact path.
    # If using generated timestamp zip, find it:
    try:
        real_cpdb_path = glob.glob(os.path.join(CPDB_TARGET_DIR, "cellphonedb*.zip"))[0]
    except IndexError:
        raise FileNotFoundError("CellPhoneDB zip file not found. Database creation might have failed.")

    cpdb_results_BO = run_cpdb_analysis("BO", INPUTS["BO"], real_cpdb_path, CPDB_OUTPUT_DIR)
    cpdb_results_AO = run_cpdb_analysis("AO", INPUTS["AO"], real_cpdb_path, CPDB_OUTPUT_DIR)
    
    # 3. Compare Results
    compare_interactions(cpdb_results_BO, cpdb_results_AO)
    
    # 4. Visualization
    # Load AnnData objects required for plotting
    print("\n[Data] Loading AnnData objects for plotting...")
    adata_BO = anndata.read_h5ad(INPUTS["BO"]["h5ad"])
    adata_AO = anndata.read_h5ad(INPUTS["AO"]["h5ad"])
    
    generate_plots(adata_BO, adata_AO, cpdb_results_BO, cpdb_results_AO)
    
    print("\nAll tasks completed successfully.")
