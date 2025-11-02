#!/usr/bin/env python3
"""
Verify integrity of downloaded demo data
"""

import sys
import logging
import scanpy as sc
import json

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def verify_h5ad(filepath):
    """Verify h5ad file can be loaded and has expected structure"""
    try:
        logger.info(f"Verifying {filepath}")
        adata = sc.read_h5ad(filepath)
        
        # Basic checks
        assert adata.n_obs > 0, "No cells in dataset"
        assert adata.n_vars > 0, "No genes in dataset"
        assert adata.X is not None, "No expression matrix"
        
        # Check for perturbation column
        possible_keys = ["gene", "target_gene", "target_gene_name", "perturbation"]
        has_pert = any(key in adata.obs.columns for key in possible_keys)
        assert has_pert, "No perturbation column found"
        
        logger.info("✓ Data verification passed")
        return True
        
    except Exception as e:
        logger.error(f"✗ Verification failed: {e}")
        return False


def main():
    if 'snakemake' in globals():
        input_file = snakemake.input[0]
    else:
        input_file = sys.argv[1]
    
    if not verify_h5ad(input_file):
        sys.exit(1)


if __name__ == "__main__":
    main()