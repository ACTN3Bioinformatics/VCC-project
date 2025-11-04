#!/usr/bin/env python3
"""
Data Splitting

Create train/validation/test splits using leave-genes-out strategy.
This ensures test perturbations are completely unseen during training.
"""

import sys
import logging
import json
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
from sklearn.model_selection import train_test_split

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def leave_genes_out_split(adata, perturbation_key, train_ratio, val_ratio, test_ratio, random_state=42):
    """
    Split data by perturbations (leave-genes-out strategy)
    
    This ensures that test perturbations are completely unseen during training,
    which is the most realistic evaluation for perturbation prediction.
    
    Parameters:
    -----------
    adata : AnnData
        Input dataset
    perturbation_key : str
        Column with perturbation labels
    train_ratio : float
        Proportion of perturbations for training
    val_ratio : float
        Proportion of perturbations for validation
    test_ratio : float
        Proportion of perturbations for testing
    random_state : int
        Random seed
    
    Returns:
    --------
    splits : dict
        Dictionary with train, val, test AnnData objects and split info
    """
    logger.info("="*60)
    logger.info("LEAVE-GENES-OUT DATA SPLITTING")
    logger.info("="*60)
    logger.info(f"Strategy: Split by {perturbation_key}")
    logger.info(f"Ratios: {train_ratio:.0%} train, {val_ratio:.0%} val, {test_ratio:.0%} test")
    
    # Check ratios sum to 1
    if not np.isclose(train_ratio + val_ratio + test_ratio, 1.0):
        raise ValueError(f"Ratios must sum to 1.0, got {train_ratio + val_ratio + test_ratio}")
    
    # Get unique perturbations
    perturbations = adata.obs[perturbation_key].unique()
    n_perturbations = len(perturbations)
    
    logger.info(f"\nTotal perturbations: {n_perturbations}")
    logger.info(f"Total cells: {adata.n_obs:,}")
    
    # Shuffle perturbations
    np.random.seed(random_state)
    perturbations_shuffled = np.random.permutation(perturbations)
    
    # Split perturbations
    n_train = int(n_perturbations * train_ratio)
    n_val = int(n_perturbations * val_ratio)
    # n_test = remainder
    
    train_perts = perturbations_shuffled[:n_train]
    val_perts = perturbations_shuffled[n_train:n_train+n_val]
    test_perts = perturbations_shuffled[n_train+n_val:]
    
    logger.info(f"\nPerturbations per split:")
    logger.info(f"  Train: {len(train_perts)} perturbations")
    logger.info(f"  Val:   {len(val_perts)} perturbations")
    logger.info(f"  Test:  {len(test_perts)} perturbations")
    
    # Split cells by perturbation
    train_mask = adata.obs[perturbation_key].isin(train_perts)
    val_mask = adata.obs[perturbation_key].isin(val_perts)
    test_mask = adata.obs[perturbation_key].isin(test_perts)
    
    adata_train = adata[train_mask].copy()
    adata_val = adata[val_mask].copy()
    adata_test = adata[test_mask].copy()
    
    logger.info(f"\nCells per split:")
    logger.info(f"  Train: {adata_train.n_obs:,} cells ({100*adata_train.n_obs/adata.n_obs:.1f}%)")
    logger.info(f"  Val:   {adata_val.n_obs:,} cells ({100*adata_val.n_obs/adata.n_obs:.1f}%)")
    logger.info(f"  Test:  {adata_test.n_obs:,} cells ({100*adata_test.n_obs/adata.n_obs:.1f}%)")
    
    # Verify no overlap
    assert len(set(train_perts) & set(val_perts)) == 0, "Train-Val overlap!"
    assert len(set(train_perts) & set(test_perts)) == 0, "Train-Test overlap!"
    assert len(set(val_perts) & set(test_perts)) == 0, "Val-Test overlap!"
    
    logger.info("\n✓ No perturbation overlap between splits")
    
    # Create split info
    split_info = {
        'strategy': 'leave_genes_out',
        'perturbation_key': perturbation_key,
        'train_ratio': train_ratio,
        'val_ratio': val_ratio,
        'test_ratio': test_ratio,
        'random_state': random_state,
        'n_perturbations_total': n_perturbations,
        'n_perturbations_train': len(train_perts),
        'n_perturbations_val': len(val_perts),
        'n_perturbations_test': len(test_perts),
        'n_cells_total': adata.n_obs,
        'n_cells_train': adata_train.n_obs,
        'n_cells_val': adata_val.n_obs,
        'n_cells_test': adata_test.n_obs,
        'train_perturbations': train_perts.tolist(),
        'val_perturbations': val_perts.tolist(),
        'test_perturbations': test_perts.tolist()
    }
    
    # Add split labels to obs
    adata_train.obs['split'] = 'train'
    adata_val.obs['split'] = 'val'
    adata_test.obs['split'] = 'test'
    
    return {
        'train': adata_train,
        'val': adata_val,
        'test': adata_test,
        'info': split_info
    }


def random_cell_split(adata, train_ratio, val_ratio, test_ratio, random_state=42):
    """
    Split data randomly by cells (NOT RECOMMENDED for perturbation prediction)
    
    This allows same perturbation in train and test, which inflates performance.
    Only use for baseline comparison.
    """
    logger.info("="*60)
    logger.info("RANDOM CELL SPLITTING")
    logger.info("="*60)
    logger.warning("⚠ This strategy is NOT recommended for perturbation prediction!")
    logger.warning("⚠ It allows data leakage and inflates performance metrics")
    
    # Random split
    train_size = train_ratio
    temp_size = val_ratio + test_ratio
    
    # First split: train vs (val+test)
    train_idx, temp_idx = train_test_split(
        np.arange(adata.n_obs),
        train_size=train_size,
        random_state=random_state,
        shuffle=True
    )
    
    # Second split: val vs test
    val_size = val_ratio / temp_size
    val_idx, test_idx = train_test_split(
        temp_idx,
        train_size=val_size,
        random_state=random_state,
        shuffle=True
    )
    
    # Create splits
    adata_train = adata[train_idx].copy()
    adata_val = adata[val_idx].copy()
    adata_test = adata[test_idx].copy()
    
    # Add split labels
    adata_train.obs['split'] = 'train'
    adata_val.obs['split'] = 'val'
    adata_test.obs['split'] = 'test'
    
    split_info = {
        'strategy': 'random_cells',
        'train_ratio': train_ratio,
        'val_ratio': val_ratio,
        'test_ratio': test_ratio,
        'random_state': random_state,
        'n_cells_total': adata.n_obs,
        'n_cells_train': adata_train.n_obs,
        'n_cells_val': adata_val.n_obs,
        'n_cells_test': adata_test.n_obs
    }
    
    return {
        'train': adata_train,
        'val': adata_val,
        'test': adata_test,
        'info': split_info
    }


def main():
    """Main execution function"""
    
    if 'snakemake' not in globals():
        logger.error("This script must be run via Snakemake")
        sys.exit(1)
    
    # Get parameters
    input_path = snakemake.input.h5ad
    output_train = snakemake.output.train
    output_val = snakemake.output.val
    output_test = snakemake.output.test
    output_info = snakemake.output.split_info
    
    train_ratio = snakemake.params.train_ratio
    val_ratio = snakemake.params.val_ratio
    test_ratio = snakemake.params.test_ratio
    gene_split = snakemake.params.gene_split
    random_state = snakemake.params.random_state
    perturbation_key = snakemake.params.perturbation_key
    
    logger.info(f"Input: {input_path}")
    
    # Load data
    logger.info("\nLoading data...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Check perturbation key
    if gene_split and perturbation_key not in adata.obs.columns:
        logger.error(f"Perturbation key '{perturbation_key}' not found!")
        logger.error(f"Available columns: {adata.obs.columns.tolist()}")
        sys.exit(1)
    
    # Perform split
    if gene_split:
        splits = leave_genes_out_split(
            adata,
            perturbation_key=perturbation_key,
            train_ratio=train_ratio,
            val_ratio=val_ratio,
            test_ratio=test_ratio,
            random_state=random_state
        )
    else:
        splits = random_cell_split(
            adata,
            train_ratio=train_ratio,
            val_ratio=val_ratio,
            test_ratio=test_ratio,
            random_state=random_state
        )
    
    # Save splits
    logger.info("\n" + "="*60)
    logger.info("SAVING SPLITS")
    logger.info("="*60)
    
    for split_name, output_path in [('train', output_train), ('val', output_val), ('test', output_test)]:
        logger.info(f"Saving {split_name} split to {output_path}")
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        splits[split_name].write_h5ad(output_path, compression='gzip')
    
    # Save split info
    logger.info(f"Saving split info to {output_info}")
    with open(output_info, 'w') as f:
        json.dump(splits['info'], f, indent=2)
    
    logger.info("\n" + "="*60)
    logger.info("SPLITTING COMPLETE")
    logger.info("="*60)
    logger.info("✓ All files saved successfully")


if __name__ == "__main__":
    main()