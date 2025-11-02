#!/usr/bin/env python3
"""
Combine features with train/val/test splits
"""

import logging
import scanpy as sc

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def combine_features_splits(features_path, train_path, val_path, test_path, output_path):
    """Combine feature-engineered data with splits"""
    
    logger.info("Loading data...")
    adata_features = sc.read_h5ad(features_path)
    adata_train = sc.read_h5ad(train_path)
    adata_val = sc.read_h5ad(val_path)
    adata_test = sc.read_h5ad(test_path)
    
    # Add split annotations
    adata_train.obs['split'] = 'train'
    adata_val.obs['split'] = 'val'
    adata_test.obs['split'] = 'test'
    
    # Concatenate
    adata_combined = sc.concat([adata_train, adata_val, adata_test], join='outer')
    
    # Transfer features
    if 'varm' in dir(adata_features) and len(adata_features.varm.keys()) > 0:
        for key in adata_features.varm.keys():
            adata_combined.varm[key] = adata_features.varm[key]
    
    # Save
    logger.info(f"Saving combined data to {output_path}")
    adata_combined.write_h5ad(output_path, compression='gzip')
    logger.info("✓ Complete")


def main():
    if 'snakemake' in globals():
        combine_features_splits(
            snakemake.input.features,
            snakemake.input.train,
            snakemake.input.val,
            snakemake.input.test,
            snakemake.output.final
        )
    else:
        import sys
        combine_features_splits(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])


if __name__ == "__main__":
    main()