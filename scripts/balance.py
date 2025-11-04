#!/usr/bin/env python3
"""
Class Balancing

Balance perturbation classes via downsampling to prevent model bias.
Strategy: Downsample overrepresented perturbations to target cell count.
"""

import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from collections import Counter

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def balance_perturbations(adata, perturbation_key, target_cells, min_cells, random_state=42):
    """
    Balance perturbation classes via intelligent sampling
    
    Strategy:
    - If perturbation has > target_cells: downsample to target_cells
    - If perturbation has < min_cells: remove (too few for training)
    - Otherwise: keep all cells
    
    Parameters:
    -----------
    adata : AnnData
        Input dataset
    perturbation_key : str
        Column name for perturbation labels
    target_cells : int
        Target number of cells per perturbation
    min_cells : int
        Minimum cells required to keep perturbation
    random_state : int
        Random seed for reproducibility
    
    Returns:
    --------
    adata_balanced : AnnData
        Balanced dataset
    balance_report : dict
        Report with balancing statistics
    """
    logger.info("="*60)
    logger.info("CLASS BALANCING")
    logger.info("="*60)
    logger.info(f"Input: {adata.n_obs:,} cells")
    
    # Check if perturbation key exists
    if perturbation_key not in adata.obs.columns:
        raise ValueError(f"Perturbation key '{perturbation_key}' not found in adata.obs")
    
    # Get perturbation counts
    pert_counts = adata.obs[perturbation_key].value_counts()
    n_perturbations_before = len(pert_counts)
    
    logger.info(f"Perturbations: {n_perturbations_before}")
    logger.info(f"Cell distribution: {pert_counts.min()}-{pert_counts.max()} cells per perturbation")
    logger.info(f"Mean cells per perturbation: {pert_counts.mean():.1f} ± {pert_counts.std():.1f}")
    
    logger.info(f"\nBalancing strategy:")
    logger.info(f"  Target cells: {target_cells}")
    logger.info(f"  Minimum cells: {min_cells}")
    
    # Categorize perturbations
    oversized = pert_counts[pert_counts > target_cells]
    undersized = pert_counts[pert_counts < min_cells]
    balanced = pert_counts[(pert_counts >= min_cells) & (pert_counts <= target_cells)]
    
    logger.info(f"\nPerturbation categories:")
    logger.info(f"  Oversized (>{target_cells} cells): {len(oversized)} perturbations")
    logger.info(f"  Balanced ({min_cells}-{target_cells} cells): {len(balanced)} perturbations")
    logger.info(f"  Undersized (<{min_cells} cells): {len(undersized)} perturbations - WILL BE REMOVED")
    
    # Remove undersized perturbations
    if len(undersized) > 0:
        logger.warning(f"\nRemoving {len(undersized)} perturbations with <{min_cells} cells:")
        for pert in undersized.index[:10]:  # Show first 10
            logger.warning(f"  - {pert}: {undersized[pert]} cells")
        if len(undersized) > 10:
            logger.warning(f"  ... and {len(undersized)-10} more")
    
    # Filter to keep only perturbations with sufficient cells
    valid_perturbations = pert_counts[pert_counts >= min_cells].index
    adata = adata[adata.obs[perturbation_key].isin(valid_perturbations)].copy()
    
    logger.info(f"\nAfter removing undersized: {adata.n_obs:,} cells")
    
    # Balance by downsampling
    logger.info(f"\nDownsampling oversized perturbations...")
    np.random.seed(random_state)
    
    balanced_cells = []
    balance_stats = []
    
    for pert in adata.obs[perturbation_key].unique():
        pert_data = adata[adata.obs[perturbation_key] == pert].copy()
        n_cells = pert_data.n_obs
        
        if n_cells > target_cells:
            # Downsample
            indices = np.random.choice(n_cells, size=target_cells, replace=False)
            pert_data = pert_data[indices].copy()
            action = 'downsampled'
        else:
            # Keep all
            action = 'kept'
        
        balanced_cells.append(pert_data)
        balance_stats.append({
            'perturbation': pert,
            'original_cells': n_cells,
            'final_cells': pert_data.n_obs,
            'action': action
        })
    
    # Concatenate balanced cells
    logger.info("Concatenating balanced datasets...")
    adata_balanced = sc.concat(balanced_cells, join='outer')
    
    # Create balance report
    balance_report = {
        'n_perturbations_before': n_perturbations_before,
        'n_perturbations_after': len(adata_balanced.obs[perturbation_key].unique()),
        'n_cells_before': adata.n_obs,
        'n_cells_after': adata_balanced.n_obs,
        'n_removed_perturbations': len(undersized),
        'n_downsampled_perturbations': len(oversized),
        'n_kept_unchanged': len(balanced),
        'stats': balance_stats
    }
    
    # Summary
    logger.info("\n" + "="*60)
    logger.info("BALANCING COMPLETE")
    logger.info("="*60)
    logger.info(f"Cells: {adata.n_obs:,} → {adata_balanced.n_obs:,}")
    logger.info(f"Perturbations: {n_perturbations_before} → {balance_report['n_perturbations_after']}")
    logger.info(f"Retention: {100*adata_balanced.n_obs/adata.n_obs:.1f}% of cells")
    
    # Final distribution
    final_counts = adata_balanced.obs[perturbation_key].value_counts()
    logger.info(f"\nFinal distribution:")
    logger.info(f"  Min cells: {final_counts.min()}")
    logger.info(f"  Max cells: {final_counts.max()}")
    logger.info(f"  Mean cells: {final_counts.mean():.1f} ± {final_counts.std():.1f}")
    
    return adata_balanced, balance_report


def write_balance_report(report, output_path):
    """Write detailed balance report to text file"""
    
    with open(output_path, 'w') as f:
        f.write("="*60 + "\n")
        f.write("CLASS BALANCING REPORT\n")
        f.write("="*60 + "\n\n")
        
        f.write("SUMMARY\n")
        f.write("-"*60 + "\n")
        f.write(f"Perturbations before: {report['n_perturbations_before']}\n")
        f.write(f"Perturbations after:  {report['n_perturbations_after']}\n")
        f.write(f"Removed:              {report['n_removed_perturbations']}\n\n")
        
        f.write(f"Cells before:         {report['n_cells_before']:,}\n")
        f.write(f"Cells after:          {report['n_cells_after']:,}\n")
        f.write(f"Retention:            {100*report['n_cells_after']/report['n_cells_before']:.1f}%\n\n")
        
        f.write(f"Downsampled:          {report['n_downsampled_perturbations']} perturbations\n")
        f.write(f"Kept unchanged:       {report['n_kept_unchanged']} perturbations\n\n")
        
        f.write("\nPER-PERTURBATION DETAILS\n")
        f.write("-"*60 + "\n")
        f.write(f"{'Perturbation':<30} {'Original':>10} {'Final':>10} {'Action':>15}\n")
        f.write("-"*60 + "\n")
        
        for stat in report['stats']:
            f.write(f"{stat['perturbation']:<30} "
                   f"{stat['original_cells']:>10} "
                   f"{stat['final_cells']:>10} "
                   f"{stat['action']:>15}\n")
    
    logger.info(f"✓ Balance report saved to {output_path}")


def main():
    """Main execution function"""
    
    if 'snakemake' not in globals():
        logger.error("This script must be run via Snakemake")
        sys.exit(1)
    
    # Get parameters
    input_path = snakemake.input[0]
    output_h5ad = snakemake.output.h5ad
    output_report = snakemake.output.report
    
    perturbation_key = snakemake.params.perturbation_key
    target_cells = snakemake.params.target_cells
    min_cells = snakemake.params.min_cells
    random_state = snakemake.params.random_state
    
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_h5ad}")
    
    # Load data
    logger.info("\nLoading data...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Balance
    adata_balanced, balance_report = balance_perturbations(
        adata,
        perturbation_key=perturbation_key,
        target_cells=target_cells,
        min_cells=min_cells,
        random_state=random_state
    )
    
    # Save output
    logger.info(f"\nSaving balanced dataset to {output_h5ad}")
    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata_balanced.write_h5ad(output_h5ad, compression='gzip')
    logger.info("✓ Saved successfully")
    
    # Write report
    Path(output_report).parent.mkdir(parents=True, exist_ok=True)
    write_balance_report(balance_report, output_report)


if __name__ == "__main__":
    main()