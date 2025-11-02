#!/usr/bin/env python3
"""
Generate HTML QC report comparing raw vs filtered data
"""

import logging
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def generate_qc_report(raw_path, filtered_path, output_path):
    """Generate comprehensive QC report"""
    
    logger.info("Loading data...")
    adata_raw = sc.read_h5ad(raw_path)
    adata_filtered = sc.read_h5ad(filtered_path)
    
    # Create output directory
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    
    # Generate plots
    fig, axes = plt.subplots(3, 2, figsize=(12, 15))
    
    # Plot 1: Cell counts
    ax = axes[0, 0]
    counts = [adata_raw.n_obs, adata_filtered.n_obs]
    ax.bar(['Raw', 'Filtered'], counts, color=['lightcoral', 'lightgreen'])
    ax.set_ylabel('Number of cells')
    ax.set_title('Cell Counts')
    
    # Plot 2: Gene counts
    ax = axes[0, 1]
    counts = [adata_raw.n_vars, adata_filtered.n_vars]
    ax.bar(['Raw', 'Filtered'], counts, color=['lightcoral', 'lightgreen'])
    ax.set_ylabel('Number of genes')
    ax.set_title('Gene Counts')
    
    # Plot 3: Genes per cell
    ax = axes[1, 0]
    ax.hist(adata_raw.obs['n_genes_by_counts'], bins=50, alpha=0.5, label='Raw')
    ax.hist(adata_filtered.obs['n_genes_by_counts'], bins=50, alpha=0.5, label='Filtered')
    ax.set_xlabel('Genes per cell')
    ax.set_ylabel('Frequency')
    ax.legend()
    ax.set_title('Genes per Cell Distribution')
    
    # Plot 4: UMI counts
    ax = axes[1, 1]
    ax.hist(adata_raw.obs['total_counts'], bins=50, alpha=0.5, label='Raw')
    ax.hist(adata_filtered.obs['total_counts'], bins=50, alpha=0.5, label='Filtered')
    ax.set_xlabel('Total UMI counts')
    ax.set_ylabel('Frequency')
    ax.legend()
    ax.set_title('UMI Counts Distribution')
    
    # Plot 5: Mitochondrial percentage
    ax = axes[2, 0]
    ax.hist(adata_raw.obs['pct_counts_mt'], bins=50, alpha=0.5, label='Raw')
    ax.hist(adata_filtered.obs['pct_counts_mt'], bins=50, alpha=0.5, label='Filtered')
    ax.set_xlabel('Mitochondrial %')
    ax.set_ylabel('Frequency')
    ax.legend()
    ax.set_title('Mitochondrial Content')
    
    # Plot 6: Summary stats
    ax = axes[2, 1]
    ax.axis('off')
    summary_text = f"""
    QC Summary:
    
    Raw Data:
    - Cells: {adata_raw.n_obs:,}
    - Genes: {adata_raw.n_vars:,}
    
    Filtered Data:
    - Cells: {adata_filtered.n_obs:,} ({100*adata_filtered.n_obs/adata_raw.n_obs:.1f}%)
    - Genes: {adata_filtered.n_vars:,} ({100*adata_filtered.n_vars/adata_raw.n_vars:.1f}%)
    
    Cells removed: {adata_raw.n_obs - adata_filtered.n_obs:,}
    Genes removed: {adata_raw.n_vars - adata_filtered.n_vars:,}
    """
    ax.text(0.1, 0.5, summary_text, fontsize=10, family='monospace')
    
    plt.tight_layout()
    
    # Save as PNG
    png_path = output_path.replace('.html', '.png')
    plt.savefig(png_path, dpi=150, bbox_inches='tight')
    logger.info(f"Saved QC plots to {png_path}")
    
    # Generate HTML
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>QC Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333; }}
            img {{ max-width: 100%; height: auto; }}
        </style>
    </head>
    <body>
        <h1>Quality Control Report</h1>
        <p>Generated: {pd.Timestamp.now()}</p>
        <img src="{Path(png_path).name}" alt="QC Plots">
    </body>
    </html>
    """
    
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    logger.info(f"✓ QC report saved to {output_path}")


def main():
    if 'snakemake' in globals():
        raw_path = snakemake.input.raw
        filtered_path = snakemake.input.filtered
        output_path = snakemake.output.report
    else:
        import sys
        raw_path = sys.argv[1]
        filtered_path = sys.argv[2]
        output_path = sys.argv[3]
    
    generate_qc_report(raw_path, filtered_path, output_path)


if __name__ == "__main__":
    main()