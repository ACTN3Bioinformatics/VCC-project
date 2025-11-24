#!/usr/bin/env python3
"""
Quality Control, Filtering, Normalization and Scaling - MEMORY EFFICIENT

This script handles three Snakemake rules:
1. filter_cells - QC and filtering
2. normalize - Count normalization and log transform
3. scale - Z-score scaling

The script detects which operation to perform based on input/output paths.

The script supports backed mode for large files:
- Automatically detects if file is too large for RAM
- Uses backed='r' mode to avoid loading entire matrix
- Converts to in-memory only after filtering reduces size
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


def estimate_memory_requirement(file_path):
    """
    Estimate memory required to load file
    Returns size in GB
    """
    file_size_gb = Path(file_path).stat().st_size / 1e9
    # Rough estimate: file size * 2 (uncompressed + working memory)
    estimated_ram_gb = file_size_gb * 2
    return estimated_ram_gb


def load_data_smart(file_path, max_ram_gb=8.0):
    """
    Smart data loading based on file size
    
    - Small files (< 4GB estimated RAM): load directly
    - Large files (> 4GB estimated RAM): use backed mode
    
    Parameters:
    -----------
    file_path : str
        Path to h5ad file
    max_ram_gb : float
        Maximum RAM to use for loading (default: 8GB for 16GB system)
    
    Returns:
    --------
    adata : AnnData
        Loaded data (may be in backed mode)
    is_backed : bool
        Whether data was loaded in backed mode
    """
    estimated_ram = estimate_memory_requirement(file_path)
    logger.info(f"File size: {Path(file_path).stat().st_size / 1e9:.2f} GB")
    logger.info(f"Estimated RAM needed: {estimated_ram:.2f} GB")
    
    if estimated_ram > max_ram_gb:
        logger.warning(f"File requires > {max_ram_gb:.1f} GB RAM")
        logger.info("Using BACKED MODE (read from disk, minimal RAM)")
        adata = sc.read_h5ad(file_path, backed='r')
        return adata, True
    else:
        logger.info("Loading into memory...")
        adata = sc.read_h5ad(file_path)
        return adata, False


def detect_operation(input_path, output_path):
    """Detect which operation to perform based on file paths"""
    input_path = Path(input_path)
    output_path = Path(output_path)
    
    if 'filtered' in output_path.name:
        return 'filter'
    elif 'normalized' in output_path.name:
        return 'normalize'
    elif 'scaled' in output_path.name:
        return 'scale'
    else:
        raise ValueError(f"Cannot determine operation from paths: {input_path} -> {output_path}")


def filter_cells(adata, min_genes=200, max_genes=6000, max_pct_mt=15, 
                 min_cells_per_gene=3, qc_plots_dir=None):
    """
    Filter low-quality cells and genes
    MEMORY EFFICIENT: Works with backed mode

    Parameters:
    -----------
    adata : AnnData
        Input dataset
    min_genes : int
        Minimum number of genes per cell
    max_genes : int
        Maximum number of genes per cell (doublet threshold)
    max_pct_mt : float
        Maximum mitochondrial percentage
    min_cells_per_gene : int
        Minimum cells expressing each gene
    qc_plots_dir : str or Path
        Directory to save QC plots
    
    Returns:
    --------
    adata : AnnData
        Filtered dataset
    """
    logger.info("="*60)
    logger.info("QUALITY CONTROL AND FILTERING")
    logger.info("="*60)
    
    # Check if data is backed
    is_backed = hasattr(adata, 'isbacked') and adata.isbacked
    logger.info(f"Input: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    logger.info(f"Mode: {'BACKED (disk-based)' if is_backed else 'IN-MEMORY'}")
    
    # Calculate QC metrics if not already present
    if 'n_genes_by_counts' not in adata.obs.columns:
        logger.info("Computing QC metrics...")
        
        # Identify mitochondrial genes
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        
        sc.pp.calculate_qc_metrics(
            adata, 
            qc_vars=['mt'], 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        logger.info("✓ QC metrics computed")
    else:
        logger.info("✓ QC metrics already present")
    
    # Print summary statistics BEFORE filtering
    logger.info("\nQC Statistics (before filtering):")
    logger.info(f"  Genes per cell: {adata.obs['n_genes_by_counts'].median():.0f} (median)")
    logger.info(f"  Total counts: {adata.obs['total_counts'].median():.0f} (median)")
    logger.info(f"  MT content: {adata.obs['pct_counts_mt'].median():.2f}% (median)")
    
    # Count cells before filtering
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars
    
    # Filter cells
    logger.info("\nApplying cell filters:")
    logger.info(f"  min_genes: {min_genes}")
    logger.info(f"  max_genes: {max_genes}")
    logger.info(f"  max_pct_mt: {max_pct_mt}%")
    
    # Create filter masks
    gene_filter = (adata.obs['n_genes_by_counts'] >= min_genes) & \
                  (adata.obs['n_genes_by_counts'] <= max_genes)
    mt_filter = adata.obs['pct_counts_mt'] <= max_pct_mt
    
    # Apply filters
    if is_backed:
        # For backed mode: select indices, then load to memory
        logger.info("Converting to in-memory after filtering...")
        cell_indices = np.where(gene_filter & mt_filter)[0]
        adata = adata[cell_indices, :].to_memory()
        logger.info("✓ Converted to in-memory mode")
    else:
        # For in-memory: direct filtering
        adata = adata[gene_filter & mt_filter, :].copy()
    
    logger.info(f"  Cells passing filters: {adata.n_obs:,} ({100*adata.n_obs/n_cells_before:.1f}%)")
    logger.info(f"  Cells removed: {n_cells_before - adata.n_obs:,}")
    
    # Filter genes
    logger.info(f"\nApplying gene filter (min_cells: {min_cells_per_gene}):")
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)
    
    logger.info(f"  Genes passing filter: {adata.n_vars:,} ({100*adata.n_vars/n_genes_before:.1f}%)")
    logger.info(f"  Genes removed: {n_genes_before - adata.n_vars:,}")
    
    # Final summary
    logger.info("\n" + "="*60)
    logger.info("FILTERING COMPLETE")
    logger.info("="*60)
    logger.info(f"Final dataset: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    logger.info(f"Retention rate: {100*adata.n_obs/n_cells_before:.1f}% cells, "
                f"{100*adata.n_vars/n_genes_before:.1f}% genes")
    
    # Estimate memory usage
    if hasattr(adata.X, 'data'):
        # Sparse
        mem_mb = (adata.X.data.nbytes + adata.X.indices.nbytes + adata.X.indptr.nbytes) / 1e6
    else:
        # Dense
        mem_mb = adata.X.nbytes / 1e6
    logger.info(f"Memory usage: ~{mem_mb:.1f} MB")
    
    # Save QC plots if directory specified
    if qc_plots_dir:
        qc_plots_dir = Path(qc_plots_dir)
        qc_plots_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info(f"\nSaving QC plots to {qc_plots_dir}")
        
        try:
            # Violin plots
            sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
                         jitter=0.4, multi_panel=True, save='_after_filter.png', show=False)
            
            # Scatter plots
            sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', 
                          save='_counts_vs_genes.png', show=False)
            sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
                          save='_counts_vs_mt.png', show=False)
            
            logger.info("✓ QC plots saved")
        except Exception as e:
            logger.warning(f"Could not save QC plots: {e}")
    
    return adata


def normalize_counts(adata, target_sum=1e4, log_transform=True, regress_out=None):
    """
    Normalize counts and optionally log-transform

    Parameters:
    -----------
    adata : AnnData
        Filtered dataset
    target_sum : float
        Target sum for total-count normalization
    log_transform : bool
        Whether to apply log1p transformation
    regress_out : list of str
        Variables to regress out (optional)
    
    Returns:
    --------
    adata : AnnData
        Normalized dataset
    """
    logger.info("="*60)
    logger.info("NORMALIZATION")
    logger.info("="*60)
    
    # Store raw counts
    if 'counts' not in adata.layers:
        logger.info("Storing raw counts in layer 'counts'")
        adata.layers['counts'] = adata.X.copy()
    
    # Total-count normalization
    logger.info(f"Normalizing to {target_sum:,.0f} counts per cell...")
    sc.pp.normalize_total(adata, target_sum=target_sum)
    logger.info("✓ Total-count normalization complete")
    
    # Log transformation
    if log_transform:
        logger.info("Applying log1p transformation...")
        sc.pp.log1p(adata)
        logger.info("✓ Log transformation complete")
    
    # Store normalized counts
    adata.layers['normalized'] = adata.X.copy()
    
    # Regress out technical variables (optional)
    if regress_out and len(regress_out) > 0:
        logger.info(f"Regressing out: {', '.join(regress_out)}")
        
        # Check if variables exist
        missing_vars = [v for v in regress_out if v not in adata.obs.columns]
        if missing_vars:
            logger.warning(f"Variables not found, skipping: {missing_vars}")
            regress_out = [v for v in regress_out if v in adata.obs.columns]
        
        if regress_out:
            sc.pp.regress_out(adata, regress_out)
            logger.info("✓ Regression complete")
    
    logger.info("\n" + "="*60)
    logger.info("NORMALIZATION COMPLETE")
    logger.info("="*60)
    
    return adata


def scale_data(adata, max_value=10):
    """
    Scale data to unit variance and zero mean

    Parameters:
    -----------
    adata : AnnData
        Normalized dataset
    max_value : float
        Clip values to [-max_value, max_value]
    
    Returns:
    --------
    adata : AnnData
        Scaled dataset
    """
    logger.info("="*60)
    logger.info("SCALING")
    logger.info("="*60)
    
    # Scale
    logger.info(f"Scaling to zero mean, unit variance (max_value={max_value})...")
    sc.pp.scale(adata, max_value=max_value)
    logger.info("✓ Scaling complete")
    
    # Store scaled data
    adata.layers['scaled'] = adata.X.copy()
    
    logger.info("\n" + "="*60)
    logger.info("SCALING COMPLETE")
    logger.info("="*60)
    
    return adata


def main():
    """Main execution function"""
    
    # Get parameters from snakemake
    if 'snakemake' in globals():
        input_path = snakemake.input.h5ad if hasattr(snakemake.input, 'h5ad') else snakemake.input[0]
        output_path = snakemake.output.h5ad if hasattr(snakemake.output, 'h5ad') else snakemake.output[0]
        
        # Detect operation
        operation = detect_operation(input_path, output_path)
        
        logger.info(f"Operation: {operation}")
        logger.info(f"Input: {input_path}")
        logger.info(f"Output: {output_path}")
        
        # Load data with smart mode selection
        logger.info("\nLoading data...")
        adata, is_backed = load_data_smart(input_path, max_ram_gb=8.0)
        logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        
        # Execute appropriate operation
        if operation == 'filter':
            qc_plots_dir = snakemake.output.qc_plots if hasattr(snakemake.output, 'qc_plots') else None
            adata = filter_cells(
                adata,
                min_genes=snakemake.params.min_genes,
                max_genes=snakemake.params.max_genes,
                max_pct_mt=snakemake.params.max_pct_mt,
                min_cells_per_gene=snakemake.params.min_cells_per_gene,
                qc_plots_dir=qc_plots_dir
            )
        
        elif operation == 'normalize':
            adata = normalize_counts(
                adata,
                target_sum=snakemake.params.target_sum,
                log_transform=snakemake.params.log_transform,
                regress_out=snakemake.params.get('regress_out', [])
            )
        
        elif operation == 'scale':
            adata = scale_data(
                adata,
                max_value=snakemake.params.max_value
            )
        
        # Save output
        logger.info(f"\nSaving to {output_path}")
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path, compression='gzip')
        logger.info("✓ Saved successfully")
        
    else:
        logger.error("This script must be run via Snakemake")
        sys.exit(1)


if __name__ == "__main__":
    main()