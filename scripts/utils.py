#!/usr/bin/env python3
"""
Utility Functions

Common helper functions used across scripts.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import logging

logger = logging.getLogger(__name__)


def find_perturbation_key(adata, possible_keys=None):
    """
    Automatically find perturbation column in adata.obs
    
    Parameters:
    -----------
    adata : AnnData
        Dataset to search
    possible_keys : list of str, optional
        List of possible column names
    
    Returns:
    --------
    key : str
        Found perturbation key
    
    Raises:
    -------
    ValueError if no perturbation key found
    """
    if possible_keys is None:
        possible_keys = [
            'perturbation', 'target_gene', 'target_gene_name', 
            'gene', 'gene_symbol', 'sgRNA_target', 'guide_id'
        ]
    
    for key in possible_keys:
        if key in adata.obs.columns:
            logger.info(f"Found perturbation key: '{key}'")
            return key
    
    raise ValueError(
        f"Could not find perturbation column. "
        f"Tried: {possible_keys}. "
        f"Available columns: {adata.obs.columns.tolist()[:20]}"
    )


def check_data_integrity(adata):
    """
    Check basic data integrity
    
    Returns dictionary with issues found
    """
    issues = {}
    
    # Check for NaN values
    if hasattr(adata.X, 'toarray'):
        X_array = adata.X.toarray()
    else:
        X_array = adata.X
    
    n_nans = np.isnan(X_array).sum()
    if n_nans > 0:
        issues['nans'] = f"Found {n_nans} NaN values in expression matrix"
    
    # Check for infinite values
    n_infs = np.isinf(X_array).sum()
    if n_infs > 0:
        issues['infs'] = f"Found {n_infs} infinite values in expression matrix"
    
    # Check for negative values (after log transform, shouldn't have any)
    n_negative = (X_array < 0).sum()
    if n_negative > 0:
        issues['negative'] = f"Found {n_negative} negative values"
    
    # Check gene names are unique
    if adata.var_names.duplicated().any():
        n_dup = adata.var_names.duplicated().sum()
        issues['duplicate_genes'] = f"Found {n_dup} duplicate gene names"
    
    # Check cell names are unique
    if adata.obs_names.duplicated().any():
        n_dup = adata.obs_names.duplicated().sum()
        issues['duplicate_cells'] = f"Found {n_dup} duplicate cell names"
    
    return issues


def make_unique_names(names):
    """
    Make names unique by appending numbers to duplicates
    
    Parameters:
    -----------
    names : list or Index
        Names to make unique
    
    Returns:
    --------
    unique_names : list
        Unique names
    """
    seen = {}
    unique = []
    
    for name in names:
        if name in seen:
            seen[name] += 1
            unique.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            unique.append(name)
    
    return unique


def compute_percent_mt(adata):
    """
    Compute mitochondrial gene percentage
    
    Parameters:
    -----------
    adata : AnnData
        Dataset
    
    Returns:
    --------
    adata : AnnData
        Dataset with 'pct_counts_mt' in obs
    """
    # Identify mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    
    # Calculate percentage
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt'],
        percent_top=None,
        log1p=False,
        inplace=True
    )
    
    return adata


def safe_divide(a, b, fill_value=0):
    """
    Safely divide arrays, replacing division by zero with fill_value
    
    Parameters:
    -----------
    a : array
        Numerator
    b : array
        Denominator
    fill_value : float
        Value to use when dividing by zero
    
    Returns:
    --------
    result : array
        Division result with safe handling of zeros
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.true_divide(a, b)
        result[~np.isfinite(result)] = fill_value
    return result


def get_memory_usage(adata):
    """
    Estimate memory usage of AnnData object
    
    Returns size in MB
    """
    size_bytes = 0
    
    # X matrix
    if hasattr(adata.X, 'data'):
        # Sparse matrix
        size_bytes += adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes
    else:
        # Dense matrix
        size_bytes += adata.X.nbytes
    
    # Layers
    for layer_name, layer_data in adata.layers.items():
        if hasattr(layer_data, 'data'):
            size_bytes += layer_data.data.nbytes + layer_data.indices.nbytes + layer_data.indptr.nbytes
        else:
            size_bytes += layer_data.nbytes
    
    # Rough estimate for obs, var, etc.
    size_bytes += adata.obs.memory_usage(deep=True).sum()
    size_bytes += adata.var.memory_usage(deep=True).sum()
    
    return size_bytes / 1024 / 1024  # Convert to MB


def log_dataset_info(adata, title="Dataset Info"):
    """
    Log comprehensive dataset information
    """
    logger.info("=" * 60)
    logger.info(title)
    logger.info("=" * 60)
    logger.info(f"Cells: {adata.n_obs:,}")
    logger.info(f"Genes: {adata.n_vars:,}")
    logger.info(f"Memory: {get_memory_usage(adata):.1f} MB")
    
    if adata.layers:
        logger.info(f"Layers: {', '.join(adata.layers.keys())}")
    
    if adata.obsm:
        logger.info(f"Obsm keys: {', '.join(adata.obsm.keys())}")
    
    if adata.varm:
        logger.info(f"Varm keys: {', '.join(adata.varm.keys())}")
    
    # Check for perturbation info
    try:
        pert_key = find_perturbation_key(adata)
        n_perts = adata.obs[pert_key].nunique()
        logger.info(f"Perturbations: {n_perts} ({pert_key})")
    except ValueError:
        logger.info("Perturbations: Not found")
    
    logger.info("=" * 60)


def save_metadata(adata, output_path, additional_info=None):
    """
    Save dataset metadata to JSON file
    
    Parameters:
    -----------
    adata : AnnData
        Dataset
    output_path : str or Path
        Path to save metadata JSON
    additional_info : dict, optional
        Additional information to include
    """
    import json
    from datetime import datetime
    
    metadata = {
        'n_cells': int(adata.n_obs),
        'n_genes': int(adata.n_vars),
        'layers': list(adata.layers.keys()) if adata.layers else [],
        'obsm_keys': list(adata.obsm.keys()) if adata.obsm else [],
        'varm_keys': list(adata.varm.keys()) if adata.varm else [],
        'obs_columns': list(adata.obs.columns),
        'var_columns': list(adata.var.columns),
        'memory_mb': get_memory_usage(adata),
        'timestamp': datetime.now().isoformat()
    }
    
    # Add perturbation info if available
    try:
        pert_key = find_perturbation_key(adata)
        metadata['perturbation_key'] = pert_key
        metadata['n_perturbations'] = int(adata.obs[pert_key].nunique())
        metadata['perturbations'] = adata.obs[pert_key].unique().tolist()
    except ValueError:
        pass
    
    # Add additional info
    if additional_info:
        metadata.update(additional_info)
    
    # Save
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"✓ Metadata saved to {output_path}")