#!/usr/bin/env python3
"""
Download and prepare demonstration dataset - MEMORY EFFICIENT VERSION
Replogle et al. 2022 - K562 essential Perturb-seq from Zenodo/Figshare

CRITICAL: Uses backed mode to avoid loading 10GB file into RAM
Original size: 1.44 GB (~310k cells, ~8500 genes)
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
        
        with urllib.request.urlopen(req, timeout=60) as response:
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


def get_perturbation_key(adata):
    """Auto-detect perturbation column"""
    possible_keys = [
        "gene", "target", "target_gene", "target_gene_name", 
        "perturbation", "gene_symbol", "sgRNA_target", "guide_id"
    ]
    
    for key in possible_keys:
        if key in adata.obs.columns:
            logger.info(f"Found perturbation key: '{key}'")
            return key
    
    # If not found, show available columns
    logger.error("Could not find perturbation column!")
    logger.info(f"Available columns: {adata.obs.columns.tolist()[:20]}")
    return None


def subset_data_memory_efficient(input_path, output_path, max_cells=10000, 
                                 max_perturbations=150, random_state=42):
    """
    Create optimized subset WITHOUT loading full file into memory
    
    Strategy:
    1. Open in backed mode (read from disk, don't load to RAM)
    2. Sample cells based on perturbation
    3. Only load selected cells into memory
    4. Save subset
    
    This allows processing 10GB files on 16GB RAM laptops!
    """
    logger.info("="*60)
    logger.info("MEMORY-EFFICIENT SUBSETTING")
    logger.info("="*60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_path}")
    
    # ========================================================================
    # STEP 1: Open in BACKED mode (read-only, no RAM loading)
    # ========================================================================
    logger.info("\nStep 1: Opening file in backed mode (minimal RAM usage)...")
    try:
        adata = sc.read_h5ad(input_path, backed='r')
        logger.info(f"✓ Opened: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
        logger.info(f"  Memory footprint: MINIMAL (backed mode)")
    except Exception as e:
        logger.error(f"Failed to open file: {e}")
        return False
    
    # ========================================================================
    # STEP 2: Detect perturbation key
    # ========================================================================
    logger.info("\nStep 2: Detecting perturbation column...")
    perturbation_key = get_perturbation_key(adata)
    if not perturbation_key:
        return False
    
    # ========================================================================
    # STEP 3: Analyze perturbation distribution (in backed mode)
    # ========================================================================
    logger.info("\nStep 3: Analyzing perturbation distribution...")
    pert_counts = adata.obs[perturbation_key].value_counts()
    logger.info(f"Found {len(pert_counts):,} unique perturbations")
    logger.info(f"Cell distribution: {pert_counts.min()}-{pert_counts.max()} cells per perturbation")
    
    # ========================================================================
    # STEP 4: Select perturbations to keep
    # ========================================================================
    logger.info(f"\nStep 4: Selecting top {max_perturbations} perturbations...")
    
    # Select top perturbations by frequency
    top_perts = pert_counts.head(max_perturbations).index.tolist()
    logger.info(f"Selected {len(top_perts)} perturbations")
    
    # ========================================================================
    # STEP 5: Select cell indices (without loading data)
    # ========================================================================
    logger.info(f"\nStep 5: Selecting up to {max_cells:,} cells...")
    
    # Get boolean mask for cells with selected perturbations
    cell_mask = adata.obs[perturbation_key].isin(top_perts)
    selected_indices = np.where(cell_mask)[0]
    
    logger.info(f"Found {len(selected_indices):,} cells with selected perturbations")
    
    # Randomly sample if too many cells
    if len(selected_indices) > max_cells:
        logger.info(f"Randomly sampling {max_cells:,} cells...")
        np.random.seed(random_state)
        selected_indices = np.random.choice(
            selected_indices, 
            size=max_cells, 
            replace=False
        )
        selected_indices = np.sort(selected_indices)  # Keep sorted for efficient access
    
    logger.info(f"Final selection: {len(selected_indices):,} cells")
    
    # ========================================================================
    # STEP 6: Load ONLY selected cells into memory
    # ========================================================================
    logger.info(f"\nStep 6: Loading selected cells into memory...")
    logger.info("  This may take 1-2 minutes for large selections...")
    
    try:
        # Create new AnnData with only selected cells
        # This is the ONLY time we load data into RAM
        adata_subset = adata[selected_indices, :].to_memory()
        
        logger.info(f"✓ Loaded: {adata_subset.n_obs:,} cells, {adata_subset.n_vars:,} genes")
        logger.info(f"  Perturbations: {adata_subset.obs[perturbation_key].nunique()}")
        
    except Exception as e:
        logger.error(f"Failed to load subset: {e}")
        logger.error("Tip: Try reducing max_cells further")
        return False
    
    # Close backed file
    adata.file.close()
    
    # ========================================================================
    # STEP 7: Compute basic QC metrics if not present
    # ========================================================================
    logger.info("\nStep 7: Computing QC metrics...")
    
    if 'n_genes_by_counts' not in adata_subset.obs.columns:
        # Identify mitochondrial genes
        adata_subset.var['mt'] = adata_subset.var_names.str.startswith('MT-')
        
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
    
    # ========================================================================
    # STEP 8: Save subset
    # ========================================================================
    logger.info(f"\nStep 8: Saving subset to {output_path}...")
    
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata_subset.write_h5ad(output_path, compression='gzip')
    
    final_size_mb = Path(output_path).stat().st_size / 1e6
    logger.info(f"✓ Saved ({final_size_mb:.1f} MB)")
    
    # ========================================================================
    # Summary
    # ========================================================================
    logger.info("\n" + "="*60)
    logger.info("SUBSETTING COMPLETE")
    logger.info("="*60)
    logger.info(f"Original: {adata.n_obs:,} cells")
    logger.info(f"Subset:   {adata_subset.n_obs:,} cells ({100*adata_subset.n_obs/adata.n_obs:.1f}%)")
    logger.info(f"Perturbations: {adata_subset.obs[perturbation_key].nunique()}")
    logger.info(f"File size: {final_size_mb:.1f} MB")
    logger.info("="*60)
    
    return True


def prepare_metadata(adata_path, perturbation_key, source_url=""):
    """Prepare metadata JSON for demo dataset"""
    
    # Load just obs (minimal memory)
    adata = sc.read_h5ad(adata_path, backed='r')
    
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
    
    adata.file.close()
    
    return metadata


def main():
    """Main execution function"""
    
    # Get parameters from snakemake or defaults
    if 'snakemake' in globals():
        output_h5ad = snakemake.output.h5ad
        output_metadata = snakemake.output.metadata
        url = snakemake.params.url
        alternative_url = snakemake.params.get('alternative_url', '')
        max_cells = snakemake.params.get('max_cells', 10000)
        max_perts = snakemake.params.get('target_perturbations', 150)
    else:
        # Standalone execution
        output_h5ad = "data_local/demo/replogle_subset.h5ad"
        output_metadata = "data_local/demo/metadata.json"
        url = "https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad"
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
        # STEP 1: Check if file already exists and is correct size
        # ================================================================
        if Path(output_h5ad).exists():
            size_mb = Path(output_h5ad).stat().st_size / 1e6
            logger.info(f"Output file already exists: {output_h5ad}")
            logger.info(f"Size: {size_mb:.1f} MB")
            
            # If file is small (< 800MB), it's likely already a subset
            if size_mb < 800:
                logger.info("✓ File appears to be a subset already")
                
                # Generate metadata if missing
                if not Path(output_metadata).exists():
                    logger.info("Generating metadata...")
                    adata_test = sc.read_h5ad(output_h5ad, backed='r')
                    pert_key = get_perturbation_key(adata_test)
                    adata_test.file.close()
                    
                    if pert_key:
                        metadata = prepare_metadata(output_h5ad, pert_key, url)
                        import json
                        with open(output_metadata, 'w') as f:
                            json.dump(metadata, f, indent=2)
                        logger.info(f"✓ Metadata saved to {output_metadata}")
                
                logger.info("✓ Demo data preparation complete!")
                return
            else:
                logger.info(f"File is large ({size_mb:.1f} MB), needs subsetting...")
                temp_h5ad = Path(output_h5ad)  # Use existing file
                download_needed = False
        else:
            download_needed = True
        
        # ================================================================
        # STEP 2: Download if needed
        # ================================================================
        if download_needed:
            logger.info("\n" + "="*60)
            logger.info("STEP 1: DOWNLOADING DEMO DATA")
            logger.info("="*60)
            logger.info("Source: Replogle et al. 2022 - K562 Essential Perturb-seq")
            logger.info("File size: ~1.44 GB (will be subsetted to ~500MB)")
            
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
        # STEP 3: Subset data (MEMORY EFFICIENT)
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 2: CREATING MEMORY-EFFICIENT SUBSET")
        logger.info("="*60)
        
        success = subset_data_memory_efficient(
            input_path=temp_h5ad,
            output_path=output_h5ad,
            max_cells=max_cells,
            max_perturbations=max_perts,
            random_state=42
        )
        
        if not success:
            sys.exit(1)
        
        # ================================================================
        # STEP 4: Generate metadata
        # ================================================================
        logger.info("\n" + "="*60)
        logger.info("STEP 3: GENERATING METADATA")
        logger.info("="*60)
        
        # Detect perturbation key from subset
        adata_test = sc.read_h5ad(output_h5ad, backed='r')
        pert_key = get_perturbation_key(adata_test)
        adata_test.file.close()
        
        if pert_key:
            metadata = prepare_metadata(output_h5ad, pert_key, url)
            import json
            with open(output_metadata, 'w') as f:
                json.dump(metadata, f, indent=2)
            logger.info(f"✓ Metadata saved to: {output_metadata}")
        
        # ================================================================
        # STEP 5: Cleanup
        # ================================================================
        if download_needed and temp_h5ad != Path(output_h5ad):
            logger.info("\nCleaning up temporary files...")
            temp_h5ad.unlink(missing_ok=True)
        
        # ================================================================
        # SUCCESS SUMMARY
        # ================================================================
        print("\n" + "="*60)
        print("✓ DEMO DATA PREPARATION COMPLETE!")
        print("="*60)
        print(f"Location:      {output_h5ad}")
        final_size_mb = Path(output_h5ad).stat().st_size / 1e6
        print(f"File size:     {final_size_mb:.1f} MB")
        print("="*60)
        print("\n🚀 Next Steps:")
        print("  1. Run QC: snakemake results/demo/filtered.h5ad --cores 4")
        print("  2. Full pipeline: snakemake results/demo/final.h5ad --cores 4")
        print("="*60)
        
    except Exception as e:
        logger.error(f"\n✗ ERROR during processing: {e}")
        import traceback
        traceback.print_exc()
        
        # Cleanup on error
        if temp_h5ad.exists() and download_needed:
            temp_h5ad.unlink(missing_ok=True)
        
        logger.error("\nIf download keeps failing, try manual download:")
        logger.error(f"1. Download from: {url}")
        logger.error(f"2. Save as: {temp_h5ad}")
        logger.error(f"3. Re-run: snakemake --cores 4")
        
        sys.exit(1)


if __name__ == "__main__":
    main()