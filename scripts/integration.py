#!/usr/bin/env python3
"""
Batch Integration

Correct batch effects using Harmony or BBKNN while preserving biological variation.
Includes checks for batch-perturbation confounding.
"""

import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def check_batch_perturbation_confounding(adata, batch_key, perturbation_key):
    """
    Check if batch and perturbation are confounded
    Returns confounding score (0-1, higher = more confounded)
    """
    if batch_key not in adata.obs.columns or perturbation_key not in adata.obs.columns:
        return 0.0
    
    # Create contingency table
    contingency = pd.crosstab(adata.obs[batch_key], adata.obs[perturbation_key])
    
    # Calculate confounding score (normalized mutual information)
    from sklearn.metrics import normalized_mutual_info_score
    
    score = normalized_mutual_info_score(
        adata.obs[batch_key].astype(str),
        adata.obs[perturbation_key].astype(str)
    )
    
    return score


def integrate_harmony(adata, batch_key, perturbation_key=None):
    """
    Batch correction using Harmony
    
    Parameters:
    -----------
    adata : AnnData
        Input dataset
    batch_key : str
        Column with batch information
    perturbation_key : str, optional
        Column with perturbation information (to preserve)
    
    Returns:
    --------
    adata : AnnData
        Dataset with batch-corrected PCA
    """
    logger.info("Using Harmony for batch correction")
    
    # Compute PCA if not present
    if 'X_pca' not in adata.obsm.keys():
        logger.info("Computing PCA (50 components)...")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
        logger.info("✓ PCA computed")
    
    # Run Harmony
    try:
        import harmonypy as hm
        
        logger.info("Running Harmony...")
        
        # Prepare data
        data_mat = adata.obsm['X_pca']
        meta_data = adata.obs[[batch_key]].copy()
        
        # Add perturbation as covariate if specified
        vars_use = [batch_key]
        if perturbation_key and perturbation_key in adata.obs.columns:
            meta_data[perturbation_key] = adata.obs[perturbation_key]
            # Don't correct perturbation signal, only batch
            logger.info(f"Preserving {perturbation_key} signal")
        
        # Run Harmony
        ho = hm.run_harmony(
            data_mat,
            meta_data,
            vars_use,
            max_iter_harmony=20,
            verbose=False
        )
        
        # Store corrected PCA
        adata.obsm['X_pca_harmony'] = ho.Z_corr.T
        adata.obsm['X_pca_original'] = adata.obsm['X_pca'].copy()
        adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']
        
        logger.info("✓ Harmony integration complete")
        
    except ImportError:
        logger.error("Harmony not installed! Install with: pip install harmonypy")
        logger.error("Skipping batch correction")
        return adata
    
    return adata


def integrate_bbknn(adata, batch_key):
    """
    Batch correction using BBKNN (Batch Balanced KNN)
    
    Parameters:
    -----------
    adata : AnnData
        Input dataset
    batch_key : str
        Column with batch information
    
    Returns:
    --------
    adata : AnnData
        Dataset with batch-corrected neighbors
    """
    logger.info("Using BBKNN for batch correction")
    
    # Compute PCA if not present
    if 'X_pca' not in adata.obsm.keys():
        logger.info("Computing PCA (50 components)...")
        sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
        logger.info("✓ PCA computed")
    
    # Run BBKNN
    try:
        import bbknn
        
        logger.info("Running BBKNN...")
        bbknn.bbknn(
            adata,
            batch_key=batch_key,
            n_pcs=30,
            neighbors_within_batch=3
        )
        logger.info("✓ BBKNN integration complete")
        
    except ImportError:
        logger.error("BBKNN not installed! Install with: pip install bbknn")
        logger.error("Skipping batch correction")
        return adata
    
    return adata


def main():
    """Main execution function"""
    
    if 'snakemake' not in globals():
        logger.error("This script must be run via Snakemake")
        sys.exit(1)
    
    # Get parameters
    input_path = snakemake.input[0]
    output_path = snakemake.output[0]
    
    batch_key = snakemake.params.batch_key
    method = snakemake.params.method
    run_integration = snakemake.params.run
    
    logger.info("="*60)
    logger.info("BATCH INTEGRATION")
    logger.info("="*60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")
    logger.info(f"Method: {method}")
    logger.info(f"Run integration: {run_integration}")
    
    # Load data
    logger.info("\nLoading data...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Check if batch correction should run
    if not run_integration:
        logger.info("\nBatch correction disabled in config - copying input to output")
        adata.write_h5ad(output_path, compression='gzip')
        logger.info("✓ Done")
        return
    
    # Check if batch key exists
    if batch_key not in adata.obs.columns:
        logger.warning(f"\nBatch key '{batch_key}' not found in data!")
        logger.warning("Available columns: " + ", ".join(adata.obs.columns[:10]))
        logger.warning("Skipping batch correction - copying input to output")
        adata.write_h5ad(output_path, compression='gzip')
        return
    
    # Check number of batches
    n_batches = adata.obs[batch_key].nunique()
    logger.info(f"\nNumber of batches: {n_batches}")
    
    if n_batches <= 1:
        logger.warning("Only one batch found - batch correction not needed")
        logger.warning("Copying input to output")
        adata.write_h5ad(output_path, compression='gzip')
        return
    
    # Check for perturbation key
    perturbation_key = None
    possible_pert_keys = ['perturbation', 'target_gene', 'target_gene_name', 'gene']
    for key in possible_pert_keys:
        if key in adata.obs.columns:
            perturbation_key = key
            break
    
    # Check confounding if perturbation key found
    if perturbation_key:
        confounding_score = check_batch_perturbation_confounding(
            adata, batch_key, perturbation_key
        )
        logger.info(f"\nBatch-perturbation confounding score: {confounding_score:.3f}")
        
        if confounding_score > 0.7:
            logger.warning("⚠ HIGH CONFOUNDING between batch and perturbation!")
            logger.warning("⚠ Batch correction may remove perturbation signal!")
            logger.warning("⚠ Proceed with caution")
    
    # Run integration
    logger.info(f"\nRunning {method} integration...")
    
    if method.lower() == 'harmony':
        adata = integrate_harmony(adata, batch_key, perturbation_key)
    elif method.lower() == 'bbknn':
        adata = integrate_bbknn(adata, batch_key)
    else:
        logger.error(f"Unknown integration method: {method}")
        logger.error("Supported methods: harmony, bbknn")
        sys.exit(1)
    
    # Recompute UMAP with integrated data
    logger.info("\nRecomputing UMAP with integrated data...")
    
    if 'X_pca_harmony' in adata.obsm.keys():
        # Use Harmony-corrected PCA
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=30)
    else:
        # BBKNN already computed neighbors
        pass
    
    sc.tl.umap(adata)
    logger.info("✓ UMAP computed")
    
    # Save output
    logger.info(f"\nSaving integrated dataset to {output_path}")
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression='gzip')
    
    logger.info("\n" + "="*60)
    logger.info("INTEGRATION COMPLETE")
    logger.info("="*60)
    logger.info("✓ Saved successfully")


if __name__ == "__main__":
    main()