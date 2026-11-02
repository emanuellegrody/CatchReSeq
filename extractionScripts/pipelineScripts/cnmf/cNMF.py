"""
cNMF Analysis Pipeline with SLURM Support

Complete pipeline for consensus NMF analysis:
  1. preprocess: Load h5ad, filter genes, prepare counts matrix
  2. prepare: Initialize cNMF and set up factorization jobs  
  3. factorize: Run NMF (parallelized via SLURM array jobs)
  4. combine: Merge factorization results from all workers
  5. k_select: Generate K selection plot and export stability metrics
     ** USER INSPECTS PLOT AND CHOOSES K **
  6. postprocess: Run consensus clustering and export results for R

Usage (two-phase workflow):
    # Phase 1: run via cNMF_workflow.sh (submits SLURM chain)
    submissionScripts/Figure3/cNMF_workflow.sh <sample>

    # Phase 1 produces: {RUN_NAME}.k_selection.png
    # User inspects plot, picks K

    # Phase 2: run postprocess with chosen K
    submissionScripts/Figure3/cNMF_postprocess.sh <run_name> <K> [density]

Results are exported as CSV files for downstream analysis in R.
"""

import os
import sys
import pickle
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc

# Use non-interactive backend before importing pyplot.
# Required on headless SLURM nodes where no display server is available.
import matplotlib
matplotlib.use('Agg')

from cnmf import cNMF

warnings.filterwarnings('ignore')

# =============================================================================
# Configuration (set via environment variables from shell script)
# =============================================================================
OUTPUT_DIR = os.environ.get('CNMF_OUTPUT_DIR', '.')
RUN_NAME = os.environ.get('CNMF_RUN_NAME', 'cnmf_run')
INPUT_H5AD = os.environ.get('CNMF_INPUT_H5AD', '')
COUNTS_FILE = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_adata_for_cnmf.h5ad")

# Parse K_RANGE from comma-separated string (e.g., "5,6,7,8,9,10")
K_RANGE_STR = os.environ.get('CNMF_K_RANGE', '5,6,7,8,9,10')
K_RANGE = [int(k) for k in K_RANGE_STR.split(',')]

N_ITER = int(os.environ.get('CNMF_N_ITER', 100))
N_HVG = int(os.environ.get('CNMF_N_HVG', 2000))
SEED = int(os.environ.get('CNMF_SEED', 42))
N_WORKERS = int(os.environ.get('CNMF_N_WORKERS', 4))

CNMF_PICKLE_PATH = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_cnmf_object.pkl")


def save_cnmf_object(cnmf_obj, path):
    """Serialize cNMF object for persistence between steps."""
    with open(path, 'wb') as f:
        pickle.dump(cnmf_obj, f)
    print(f"cNMF object saved: {path}")


def load_cnmf_object(path):
    """Load serialized cNMF object."""
    with open(path, 'rb') as f:
        return pickle.load(f)


def preprocess():
    """
    Load raw h5ad, ensure raw counts in .X, filter genes, and save prepared data.

    cNMF requires raw counts in .X. This step ensures:
    - Raw counts are in the main matrix (not normalized data)
    - Genes expressed in too few cells are removed

    Note: HVG selection is handled internally by cNMF during prepare().
    """
    print(f"Loading data: {INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)
    print(f"  Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    # Ensure raw counts are in .X for cNMF
    # cNMF expects integer counts, not normalized data
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
        print("  Using counts from layers['counts']")

    # Save prepared data
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    adata.write(COUNTS_FILE)
    print(f"  Saved prepared data: {COUNTS_FILE}")
    print("Preprocessing complete.")


def prepare():
    """
    Initialize cNMF analysis and prepare factorization jobs.

    This step:
    - Normalizes counts per cell (TPM-like)
    - Selects highly variable genes
    - Creates iteration files for parallel factorization
    """
    print(f"Preparing cNMF analysis: {RUN_NAME}")
    print(f"  K range: {K_RANGE}")
    print(f"  Iterations per K: {N_ITER}")
    print(f"  HVG count: {N_HVG}")

    cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)

    cnmf_obj.prepare(
        counts_fn=COUNTS_FILE,
        components=K_RANGE,
        n_iter=N_ITER,
        seed=SEED,
        num_highvar_genes=N_HVG,
        genes_file=None
    )

    save_cnmf_object(cnmf_obj, CNMF_PICKLE_PATH)
    print("Preparation complete.")


def factorize():
    """
    Run NMF factorization for a single SLURM array task.

    Each worker processes a subset of the total iterations.
    Worker index is read from SLURM_ARRAY_TASK_ID environment variable.
    """
    worker_i = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
    print(f"Factorization worker {worker_i}/{N_WORKERS}")

    # Reinitialize cNMF object (required for each worker)
    cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)
    cnmf_obj.factorize(worker_i=worker_i, total_workers=N_WORKERS)

    print(f"Worker {worker_i} complete.")


def combine():
    """
    Combine factorization results from all parallel workers.

    Merges the NMF results from each worker into unified result files.
    """
    print("Combining factorization results...")

    cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)
    cnmf_obj.combine()

    save_cnmf_object(cnmf_obj, CNMF_PICKLE_PATH)
    print("Results combined.")


def k_select():
    """
    Generate the K selection plot and export stability/error metrics.

    This calls cNMF's built-in k_selection_plot(), which:
    - Computes mean stability (silhouette score) for each K
    - Computes mean reconstruction error (Frobenius norm) for each K
    - Saves a diagnostic PNG and a tab-separated stats file

    Output files (inside {OUTPUT_DIR}/{RUN_NAME}/):
    - {RUN_NAME}.k_selection.png: Two-panel plot (stability vs K, error vs K)
    - {RUN_NAME}.k_selection_stats.txt: Tab-separated metrics per K

    Additionally exports a CSV copy to the top-level output directory.

    The user should inspect the PNG to choose K. Per Kotliar et al. 2019
    (eLife 8:e43803), select the largest K that remains reasonably stable
    and/or corresponds to a local maximum in stability.
    """
    print("Generating K selection plot...")

    cnmf_obj = cNMF(output_dir=OUTPUT_DIR, name=RUN_NAME)
    cnmf_obj.k_selection_plot()

    # The library writes the plot to {output_dir}/{name}/{name}.k_selection.png
    plot_path = os.path.join(OUTPUT_DIR, RUN_NAME, f"{RUN_NAME}.k_selection.png")
    stats_path = os.path.join(OUTPUT_DIR, RUN_NAME, f"{RUN_NAME}.k_selection_stats.txt")

    if os.path.exists(plot_path):
        print(f"  K selection plot: {plot_path}")
    else:
        print(f"  WARNING: Expected plot not found at {plot_path}")

    # Export stats to CSV in the top-level output directory for convenience
    if os.path.exists(stats_path):
        stats_df = pd.read_csv(stats_path, sep='\t')
        csv_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_k_selection.csv")
        stats_df.to_csv(csv_path, index=False)
        print(f"  K selection stats: {stats_path}")
        print(f"  K selection CSV:   {csv_path}")
        print("")
        print("  K selection summary:")
        print(stats_df.to_string(index=False))
    else:
        print(f"  WARNING: Stats file not found at {stats_path}")

    save_cnmf_object(cnmf_obj, CNMF_PICKLE_PATH)

    print("")
    print("==============================================")
    print("Phase 1 complete. Inspect the K selection plot:")
    print(f"  {plot_path}")
    print("")
    print("Then run Phase 2 with your chosen K:")
    print("  submissionScripts/Figure3/cNMF_postprocess.sh <run_name> <K> [density]")
    print("==============================================")


def postprocess(k, density_threshold):
    """
    Run consensus clustering and export results for R analysis.

    Parameters
    ----------
    k : int
        Number of patterns (components) to use. Select based on k_selection_plot.
    density_threshold : float
        Density threshold for filtering outlier patterns (typically 0.05-0.2).
        Lower values are more stringent.

    Exports
    -------
    - {RUN_NAME}_usage_K{k}.csv: Cell x Pattern usage matrix
    - {RUN_NAME}_spectra_K{k}.csv: Gene x Pattern loadings (z-scores)
    - {RUN_NAME}_spectra_tpm_K{k}.csv: Gene x Pattern loadings (TPM-normalized)
    - {RUN_NAME}_top_genes_K{k}.csv: Top genes per pattern
    """
    print(f"Postprocessing: K={k}, density_threshold={density_threshold}")

    cnmf_obj = load_cnmf_object(CNMF_PICKLE_PATH)

    # Run consensus clustering
    # This filters outlier iterations and computes final consensus matrices
    cnmf_obj.consensus(k=k, density_threshold=density_threshold, show_clustering=True)

    # Load results
    # usage: cells x patterns (how much each cell uses each pattern)
    # spectra_scores: genes x patterns (z-scored gene loadings)
    # spectra_tpm: genes x patterns (TPM-normalized gene loadings)
    # top_genes: genes x patterns (sorted by loading)
    usage, spectra_scores, spectra_tpm, top_genes = cnmf_obj.load_results(
        K=k, density_threshold=density_threshold
    )

    print(f"  Usage matrix: {usage.shape[0]} cells x {usage.shape[1]} patterns")
    print(f"  Spectra matrix: {spectra_scores.shape[0]} genes x {spectra_scores.shape[1]} patterns")

    # Rename pattern columns for clarity
    pattern_names = [f"Pattern_{i+1}" for i in range(k)]
    usage.columns = pattern_names
    spectra_scores.columns = pattern_names
    spectra_tpm.columns = pattern_names
    top_genes.columns = pattern_names

    # Export to CSV for R
    usage_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_usage_K{k}.csv")
    spectra_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_spectra_K{k}.csv")
    spectra_tpm_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_spectra_tpm_K{k}.csv")
    top_genes_path = os.path.join(OUTPUT_DIR, f"{RUN_NAME}_top_genes_K{k}.csv")

    usage.to_csv(usage_path)
    spectra_scores.to_csv(spectra_path)
    spectra_tpm.to_csv(spectra_tpm_path)
    top_genes.to_csv(top_genes_path)

    print(f"  Exported: {usage_path}")
    print(f"  Exported: {spectra_path}")
    print(f"  Exported: {spectra_tpm_path}")
    print(f"  Exported: {top_genes_path}")

    # Print top genes summary
    print("\nTop 10 genes per pattern:")
    for col in top_genes.columns:
        genes = top_genes[col].head(10).index.tolist()
        print(f"  {col}: {', '.join(genes)}")

    print("\nPostprocessing complete. Results ready for R analysis.")


def main():
    parser = argparse.ArgumentParser(
        description='cNMF analysis pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Steps (two-phase workflow):
  Phase 1 (automated via cNMF_workflow.sh):
    1. preprocess  - Prepare h5ad file (filter genes, ensure raw counts)
    2. prepare     - Initialize cNMF (HVG selection, create iteration files)
    3. factorize   - Run NMF (submit via SLURM array job)
    4. combine     - Merge results from parallel workers
    5. k_select    - Generate K selection plot and stability metrics

  ** User inspects K selection plot and chooses K **

  Phase 2 (run via cNMF_postprocess.sh <K> [density]):
    6. postprocess - Consensus clustering and export CSVs for R

Environment variables (set in shell script):
  CNMF_OUTPUT_DIR   - Output directory for results
  CNMF_RUN_NAME     - Name prefix for output files
  CNMF_INPUT_H5AD   - Path to input h5ad file
  CNMF_K_RANGE      - Comma-separated K values (e.g., "5,6,7,8,9,10")
  CNMF_N_ITER       - Iterations per K (default: 100)
  CNMF_N_HVG        - Number of highly variable genes (default: 2000)
  CNMF_SEED         - Random seed (default: 42)
  CNMF_N_WORKERS    - Number of parallel workers (default: 4)
        """
    )

    subparsers = parser.add_subparsers(dest='command', help='Pipeline step to run')

    subparsers.add_parser('preprocess', help='Prepare input h5ad file')
    subparsers.add_parser('prepare', help='Initialize cNMF analysis')
    subparsers.add_parser('factorize', help='Run factorization (use via SLURM)')
    subparsers.add_parser('combine', help='Combine factorization results')
    subparsers.add_parser('k_select', help='Generate K selection plot (end of Phase 1)')

    postprocess_parser = subparsers.add_parser('postprocess', help='Export results for R (Phase 2)')
    postprocess_parser.add_argument('--k', type=int, required=True,
                                    help='Number of patterns to use')
    postprocess_parser.add_argument('--density', type=float, default=0.1,
                                    help='Density threshold (default: 0.1)')

    subparsers.add_parser('config', help='Print current configuration')

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == 'config':
        print("Current configuration:")
        print(f"  OUTPUT_DIR: {OUTPUT_DIR}")
        print(f"  RUN_NAME: {RUN_NAME}")
        print(f"  INPUT_H5AD: {INPUT_H5AD}")
        print(f"  COUNTS_FILE: {COUNTS_FILE}")
        print(f"  K_RANGE: {K_RANGE}")
        print(f"  N_ITER: {N_ITER}")
        print(f"  N_HVG: {N_HVG}")
        print(f"  SEED: {SEED}")
        print(f"  N_WORKERS: {N_WORKERS}")
        sys.exit(0)

    if args.command == 'preprocess' and not INPUT_H5AD:
        print("ERROR: CNMF_INPUT_H5AD environment variable not set")
        sys.exit(1)

    commands = {
        'preprocess': preprocess,
        'prepare': prepare,
        'factorize': factorize,
        'combine': combine,
        'k_select': k_select,
    }

    if args.command in commands:
        commands[args.command]()
    elif args.command == 'postprocess':
        postprocess(args.k, args.density)


if __name__ == "__main__":
    main()