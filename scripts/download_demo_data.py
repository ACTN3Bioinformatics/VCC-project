#!/usr/bin/env python3
"""
Download and prepare demonstration dataset
Replogle et al. 2022 - K562 essential Perturb-seq from Zenodo/Figshare

Source: 
- Primary: Zenodo (scPerturb database)
- Backup: Figshare (original authors)

Original size: 1.55 GB (~188k cells, ~2000 perturbations)
After subsetting: ~500MB (~10k cells, ~150 perturbations)

Optimized for 16GB RAM systems (AMD Ryzen 5 7535HS)
"""

import os
import sys
import logging
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def download_file_with_progress(url, output_path, chunk_size=8192):
    """Download file with progress indicator"""
    logger.info(f"Downloading from:\n  {url}")
    
    try:
        # Try to open URL
        req = urllib.request.Request(url)
        req.add_header('User-Agent', 'Mozilla/5.0')  # Some servers require this
        
        with urllib.request.urlopen(req) as response:
            total_size = int(response.headers.get('Content-Length', 0))
            downloaded = 0
            
            logger.info(f"File size: {total_size / 1e9:.2f} GB")
            logger.info("Downloading...")
            
            with open(output_path, 'wb') as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    
                    if total_size > 0:
                        progress = (downloaded / total_size) * 100
                        mb_downloaded = downloaded / 1e6
                        mb_total = total_size / 1e6
                        print(f"\rProgress: {progress:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", 
                              end='', flush=True)
        
        print()  # New line after progress
        logger.info(f"✓ Download complete: {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"✗ Download failed: {e}")
        return False


def subset_data(adata, max_cells=None, max_perturbations=None, 
                perturbation_key="gene", random_state=42):
    """
    Subset data to specified limits
    
    Strategy:
    1. Select top N perturbations by frequency
    2. Randomly sample cells from selected perturbations
    3. Balance if needed
    """
    logger.info("="*60)
    logger.info("DATA SUBSETTING")
    logger.info("="*60)
    logger.info(f"Input: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Find perturbation key
    possible_keys = ["gene", "target", "target_gene", "target_gene_name", 
                     "perturbation", "gene_symbol", "sgRNA_target"]
    
    if perturbation_key not in adata.obs.columns:
        for key in possible_keys:
            if key in adata.obs.columns:
                perturbation_key = key
                logger.info(f"Found perturbation key: '{perturbation_key}'")
                break
        else:
            logger.error("Could not find perturbation column!")
            logger.info(f"Available columns: {adata.obs.columns.tolist()[:20]}")
            sys.exit(1)
    
    # Get perturbation counts
    pert_counts = adata.obs[perturbation_key].value_counts()
    logger.info(f"Found {len(pert_counts):,} unique perturbations")
    logger.info(f"Cell distribution: {pert_counts.min()}-{pert_counts.max()} cells per perturbation")
    
    # Subset perturbations if needed
    if max_perturbations and len(pert_counts) > max_perturbations:
        logger.info(f"Selecting top {max_perturbations} perturbations by frequency...")
        # Select top perturbations
        top_perts = pert_counts.head(max_perturbations).index.tolist()
        adata = adata[adata.obs[perturbation_key].isin(top_perts)].copy()
        pert_counts = adata.obs[perturbation_key].value_counts()
        logger.info(f"After perturbation filter: {adata.n_obs:,} cells")
    
    # Subset cells if needed
    if max_cells and adata.n_obs > max_cells:
        logger.info(f"Randomly sampling {max_cells:,} cells...")
        np.random.seed(random_state)
        sample_idx = np.random.choice(
            adata.n_obs, 
            size=max_cells, 
            replace=False
        )
        adata = adata[sample_idx].copy()
    
    logger.info(f"Output: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    logger.info(f"Output: {adata.obs[perturbation_key].nunique()} perturbations")
    
    return adata, perturbation_key


def prepare_metadata(adata, perturbation_key="gene", source_url=""):
    """Prepare metadata JSON for demo dataset"""
    
    # Get perturbation list (truncate if too many)
    pert_list = adata.obs[perturbation_key].unique().tolist()
    if len(pert_list) > 50:
        pert_list_truncated = pert_list[:50] + [f"... and {len(pert_list)-50} more"]
    else:
        pert_list_truncated = pert_list
    
    metadata = {
        "source": "Replogle et al. 2022 (via scPerturb/Zenodo)",
        "doi": "10.1016/j.cell.2022.05.013",
        "scperturb_doi": "10.1038/s41592-023-02144-y",
        "zenodo_record": "7041849",
        "source_url": source_url,
        "cell_type": "K562 (chronic myelogenous leukemia)",
        "cell_line_source": "ATCC CCL-243",
        "technology": "Perturb-seq (CRISPRi)",
        "perturbation_type": "Essential gene knockdowns",
        "day_post_transduction": 6,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_perturbations": int(adata.obs[perturbation_key].nunique()),
        "perturbations_sample": pert_list_truncated,
        "processed_date": pd.Timestamp.now().isoformat(),
        "pipeline_version": "1.0.0",
        "note": "Subset created for VCC-project demonstration purposes",
        "hardware_optimized": "AMD Ryzen 5 7535HS, 16GB RAM",
        "citation": "Replogle JM et al. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. Cell. 2022;185(14):2559-2575."
    }
    return metadata


def main():
    """Main execution function"""
    
    # Get parameters from snakemake or defaults
    if 'snakemake' in globals():
        output_h5ad = snakemake.output.h5ad
        output_metadata = snakemake.output.metadata
        url = snakemake.params.url
        alternative_url = snakemake.params.get('alternative_url', '')
        max_cells = snakemake.params.get('max_cells', None)
        max_perts = snakemake.params.get('target_perturbations', None)
    else:
        # Standalone execution
        output_h5ad = "data_local/demo/replogle_subset.h5ad"
        output_metadata = "data_local/demo/metadata.json"
        # Zenodo - scPerturb curated dataset (PRIMARY)
        url = "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad"
        # Figshare - Original authors (BACKUP)
        alternative_url = "https://plus.figshare.com/ndownloader/files/42444534"
        max_cells = 10000
        max_perts = 150
    
    # Create output directory
    output_dir = Path(output_h5ad).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Temporary download location
    temp_h5ad = output_dir / "temp_download.h5ad"
    
    try:
        # ================================================================
        # STEP 1: Download
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 1: DOWNLOADING DEMO DATA")
        logger.info("="*60)
        logger.info("Source: Replogle et al. 2022 - K562 Essential Perturb-seq")
        logger.info("File size: ~1.55 GB (will be subsetted to ~500MB)")
        
        download_success = download_file_with_progress(url, temp_h5ad)
        
        if not download_success and alternative_url:
            logger.warning("Primary download failed, trying backup URL...")
            logger.info(f"Backup: {alternative_url}")
            download_success = download_file_with_progress(alternative_url, temp_h5ad)
        
        if not download_success:
            logger.error("All download attempts failed!")
            logger.error("\nTroubleshooting:")
            logger.error("1. Check internet connection")
            logger.error("2. Try manual download from:")
            logger.error(f"   {url}")
            logger.error(f"   or: {alternative_url}")
            logger.error("3. Place downloaded file at: data_local/demo/")
            sys.exit(1)
        
        # ================================================================
        # STEP 2: Load data
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 2: LOADING DATA")
        logger.info("="*60)
        
        logger.info("Reading h5ad file (this may take 1-2 minutes for 1.44GB file)...")
        adata = sc.read_h5ad(temp_h5ad)
        logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        
        # ================================================================
        # STEP 3: Subset data
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 3: CREATING OPTIMIZED SUBSET")
        logger.info("="*60)
        
        adata_subset, pert_key = subset_data(
            adata,
            max_cells=max_cells,
            max_perturbations=max_perts,
            random_state=42
        )
        
        # Free memory
        del adata
        
        # ================================================================
        # STEP 4: QC metrics
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 4: COMPUTING QC METRICS")
        logger.info("="*60)
        
        if 'n_genes_by_counts' not in adata_subset.obs.columns:
            logger.info("Computing QC metrics...")
            sc.pp.calculate_qc_metrics(
                adata_subset, 
                qc_vars=['mt'], 
                percent_top=None, 
                log1p=False, 
                inplace=True
            )
            logger.info("✓ QC metrics computed")
        else:
            logger.info("✓ QC metrics already present")
        
        # ================================================================
        # STEP 5: Save
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 5: SAVING PROCESSED DATA")
        logger.info("="*60)
        
        logger.info(f"Saving subset to: {output_h5ad}")
        adata_subset.write_h5ad(output_h5ad, compression='gzip')
        
        final_size_mb = Path(output_h5ad).stat().st_size / 1e6
        logger.info(f"✓ Saved ({final_size_mb:.1f} MB)")
        
        # Save metadata
        logger.info(f"Creating metadata...")
        metadata = prepare_metadata(adata_subset, pert_key, url)
        import json
        with open(output_metadata, 'w') as f:
            json.dump(metadata, f, indent=2)
        logger.info(f"✓ Metadata saved to: {output_metadata}")
        
        # ================================================================
        # Cleanup
        # ================================================================
        logger.info("\nCleaning up temporary files...")
        temp_h5ad.unlink(missing_ok=True)
        
        # ================================================================
        # SUCCESS SUMMARY
        # ================================================================
        print("\n" + "="*60)
        print("✓ DEMO DATA PREPARATION COMPLETE!")
        print("="*60)
        print(f"Cells:         {adata_subset.n_obs:,}")
        print(f"Genes:         {adata_subset.n_vars:,}")
        print(f"Perturbations: {adata_subset.obs[pert_key].nunique()}")
        print(f"File size:     {final_size_mb:.1f} MB")
        print(f"Location:      {output_h5ad}")
        print("="*60)
        print("\n📊 Data Summary:")
        print(f"  - Original dataset: Replogle et al. 2022")
        print(f"  - Cell type: K562 (myeloid leukemia)")
        print(f"  - Technology: Perturb-seq (CRISPRi)")
        print(f"  - Perturbation targets: Essential genes")
        print("\n🚀 Next Steps:")
        print("  1. Explore data:")
        print("     jupyter notebook notebooks/demo_exploration.ipynb")
        print("  2. Run QC report:")
        print("     snakemake reports/demo/qc_report.html --cores 2")
        print("  3. Run full pipeline:")
        print("     snakemake --cores 4")
        print("="*60)
        
    except Exception as e:
        logger.error(f"\n✗ ERROR during processing: {e}")
        import traceback
        traceback.print_exc()
        
        # Cleanup on error
        temp_h5ad.unlink(missing_ok=True)
        
        logger.error("\nIf download keeps failing, try manual download:")
        logger.error(f"1. Download from: {url}")
        logger.error(f"   (or backup: {alternative_url})")
        logger.error(f"2. Save as: {output_h5ad}")
        logger.error(f"3. Re-run this script")
        
        sys.exit(1)


if __name__ == "__main__":
    main()