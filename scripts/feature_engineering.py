#!/usr/bin/env python3
"""
Feature Engineering

Extract biological features to enhance model performance:
- Pathway membership (KEGG, Reactome, GO)
- TF target information
- Gene-gene interactions
- Cell state signatures
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


def extract_pathway_features(adata, databases=['KEGG', 'Reactome', 'GO_Biological_Process']):
    """
    Extract pathway membership features using gene sets
    
    Creates binary matrix: genes x pathways
    """
    logger.info("="*60)
    logger.info("PATHWAY FEATURES")
    logger.info("="*60)
    
    try:
        import gseapy as gp
        
        pathway_features = {}
        gene_names = adata.var_names.tolist()
        
        for db in databases:
            logger.info(f"\nProcessing {db}...")
            
            try:
                # Get gene sets for this database
                db_name_mapping = {
                    'KEGG': 'KEGG_2021_Human',
                    'Reactome': 'Reactome_2022',
                    'GO_Biological_Process': 'GO_Biological_Process_2023'
                }
                
                gseapy_db = db_name_mapping.get(db, db)
                
                # Download gene sets
                gene_sets = gp.get_library(name=gseapy_db, organism='Human')
                
                if not gene_sets:
                    logger.warning(f"No gene sets found for {db}")
                    continue
                
                logger.info(f"  Found {len(gene_sets)} gene sets")
                
                # Create binary matrix: genes x pathways
                pathway_matrix = []
                pathway_names = []
                
                for pathway_name, pathway_genes in gene_sets.items():
                    # Check which genes in our dataset are in this pathway
                    is_in_pathway = [1 if gene in pathway_genes else 0 for gene in gene_names]
                    
                    # Only keep pathways with at least 5 genes
                    if sum(is_in_pathway) >= 5:
                        pathway_matrix.append(is_in_pathway)
                        pathway_names.append(pathway_name)
                
                if pathway_matrix:
                    pathway_matrix = np.array(pathway_matrix).T  # Transpose: genes x pathways
                    pathway_features[db] = {
                        'matrix': pathway_matrix,
                        'names': pathway_names
                    }
                    logger.info(f"  ✓ Extracted {len(pathway_names)} pathways")
                else:
                    logger.warning(f"  No valid pathways found for {db}")
                
            except Exception as e:
                logger.warning(f"  Error processing {db}: {e}")
                continue
        
        # Store in adata
        for db, features in pathway_features.items():
            adata.varm[f'pathway_{db}'] = features['matrix']
            adata.uns[f'pathway_{db}_names'] = features['names']
        
        logger.info(f"\n✓ Pathway features extracted for {len(pathway_features)} databases")
        
        return adata
        
    except ImportError:
        logger.warning("gseapy not installed - skipping pathway features")
        logger.warning("Install with: pip install gseapy")
        return adata


def extract_tf_targets(adata):
    """
    Extract TF-target relationships
    
    Uses dorothea database of TF regulons
    """
    logger.info("\n" + "="*60)
    logger.info("TRANSCRIPTION FACTOR TARGETS")
    logger.info("="*60)
    
    try:
        # Try to load dorothea (DecoupleR database)
        import requests
        
        logger.info("Downloading dorothea TF database...")
        
        # Download dorothea regulons
        url = "https://raw.githubusercontent.com/saezlab/dorothea/master/data/dorothea_hs.csv"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            # Parse CSV
            from io import StringIO
            df = pd.read_csv(StringIO(response.text))
            
            # Filter for high confidence (A, B, C levels)
            df = df[df['confidence'].isin(['A', 'B', 'C'])]
            
            logger.info(f"  Found {df['tf'].nunique()} TFs")
            logger.info(f"  Found {len(df)} TF-target interactions")
            
            # Create TF-target matrix
            gene_names = adata.var_names.tolist()
            tfs = df['tf'].unique()
            
            tf_matrix = []
            tf_names = []
            
            for tf in tfs:
                tf_targets = df[df['tf'] == tf]['target'].tolist()
                
                # Binary vector: is each gene a target of this TF?
                is_target = [1 if gene in tf_targets else 0 for gene in gene_names]
                
                # Only keep TFs with at least 5 targets in our dataset
                if sum(is_target) >= 5:
                    tf_matrix.append(is_target)
                    tf_names.append(tf)
            
            if tf_matrix:
                tf_matrix = np.array(tf_matrix).T  # genes x TFs
                adata.varm['tf_targets'] = tf_matrix
                adata.uns['tf_names'] = tf_names
                logger.info(f"  ✓ Extracted {len(tf_names)} TF regulons")
            else:
                logger.warning("  No valid TF regulons found")
        else:
            logger.warning(f"  Failed to download dorothea: HTTP {response.status_code}")
        
        return adata
        
    except Exception as e:
        logger.warning(f"Error extracting TF targets: {e}")
        return adata


def extract_cell_signatures(adata):
    """
    Extract cell state signatures (cell cycle, stress, etc.)
    """
    logger.info("\n" + "="*60)
    logger.info("CELL STATE SIGNATURES")
    logger.info("="*60)
    
    # Define signature gene sets
    signatures = {
        'cell_cycle_G1S': ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6'],
        'cell_cycle_G2M': ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2'],
        'stress': ['HSP90AA1', 'HSPA1A', 'HSPA1B', 'DNAJB1', 'HSPH1', 'HSPA6', 'DDIT3', 'ATF4'],
        'apoptosis': ['CASP3', 'CASP7', 'CASP9', 'BAX', 'BAK1', 'BID', 'BBC3', 'PMAIP1'],
        'proliferation': ['MKI67', 'TOP2A', 'PCNA', 'CCNB1', 'CCNA2', 'CDK1', 'AURKA', 'AURKB']
    }
    
    # Compute signature scores
    for sig_name, sig_genes in signatures.items():
        # Find genes present in dataset
        sig_genes_present = [g for g in sig_genes if g in adata.var_names]
        
        if len(sig_genes_present) >= 3:  # Need at least 3 genes
            logger.info(f"  {sig_name}: {len(sig_genes_present)}/{len(sig_genes)} genes found")
            
            # Compute mean expression as signature score
            sig_data = adata[:, sig_genes_present].X
            if hasattr(sig_data, 'toarray'):
                sig_data = sig_data.toarray()
            
            adata.obs[f'{sig_name}_score'] = sig_data.mean(axis=1)
        else:
            logger.warning(f"  {sig_name}: too few genes found ({len(sig_genes_present)})")
    
    logger.info("\n✓ Cell state signatures computed")
    
    return adata


def extract_gene_statistics(adata):
    """
    Extract gene-level statistics as features
    """
    logger.info("\n" + "="*60)
    logger.info("GENE STATISTICS")
    logger.info("="*60)
    
    # Mean expression
    if 'mean_expression' not in adata.var.columns:
        mean_expr = np.array(adata.X.mean(axis=0)).flatten()
        adata.var['mean_expression'] = mean_expr
        logger.info("  ✓ Mean expression")
    
    # Variance
    if 'variance' not in adata.var.columns:
        var_expr = np.array(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X).var(axis=0)
        adata.var['variance'] = var_expr
        logger.info("  ✓ Variance")
    
    # Dropout rate
    if 'dropout_rate' not in adata.var.columns:
        X_array = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        dropout = (X_array == 0).sum(axis=0) / adata.n_obs
        adata.var['dropout_rate'] = dropout
        logger.info("  ✓ Dropout rate")
    
    # Coefficient of variation
    if 'cv' not in adata.var.columns:
        mean = adata.var['mean_expression']
        std = np.sqrt(adata.var['variance'])
        cv = std / (mean + 1e-10)  # Add small constant to avoid division by zero
        adata.var['cv'] = cv
        logger.info("  ✓ Coefficient of variation")
    
    logger.info("\n✓ Gene statistics computed")
    
    return adata


def main():
    """Main execution function"""
    
    if 'snakemake' not in globals():
        logger.error("This script must be run via Snakemake")
        sys.exit(1)
    
    # Get parameters
    input_path = snakemake.input.h5ad
    output_h5ad = snakemake.output.h5ad
    output_list = snakemake.output.feature_list
    
    use_pathways = snakemake.params.use_pathways
    use_tf_targets = snakemake.params.use_tf_targets
    use_regulatory_networks = snakemake.params.use_regulatory_networks
    databases = snakemake.params.databases
    
    logger.info("="*60)
    logger.info("FEATURE ENGINEERING")
    logger.info("="*60)
    logger.info(f"Input: {input_path}")
    logger.info(f"Output: {output_h5ad}")
    
    # Load data
    logger.info("\nLoading data...")
    adata = sc.read_h5ad(input_path)
    logger.info(f"✓ Loaded: {adata.n_obs:,} cells, {adata.n_vars:,} genes")
    
    # Track extracted features
    extracted_features = []
    
    # Extract pathway features
    if use_pathways:
        logger.info("\n" + "="*60)
        logger.info("EXTRACTING PATHWAY FEATURES")
        logger.info("="*60)
        adata = extract_pathway_features(adata, databases=databases)
        extracted_features.append('pathways')
    
    # Extract TF targets
    if use_tf_targets:
        adata = extract_tf_targets(adata)
        extracted_features.append('tf_targets')
    
    # Extract cell state signatures
    adata = extract_cell_signatures(adata)
    extracted_features.append('cell_signatures')
    
    # Extract gene statistics
    adata = extract_gene_statistics(adata)
    extracted_features.append('gene_statistics')
    
    # Regulatory networks (placeholder - requires more complex implementation)
    if use_regulatory_networks:
        logger.info("\n" + "="*60)
        logger.info("REGULATORY NETWORKS")
        logger.info("="*60)
        logger.warning("Regulatory network features not yet implemented")
        logger.warning("This would require GRN inference (e.g., SCENIC, CellOracle)")
    
    # Save feature list
    logger.info(f"\nSaving feature list to {output_list}")
    Path(output_list).parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_list, 'w') as f:
        f.write("EXTRACTED FEATURES\n")
        f.write("="*60 + "\n\n")
        
        for feature in extracted_features:
            f.write(f"- {feature}\n")
        
        f.write("\nVARM KEYS:\n")
        f.write("-"*60 + "\n")
        for key in adata.varm.keys():
            shape = adata.varm[key].shape
            f.write(f"  {key}: {shape}\n")
        
        f.write("\nVAR COLUMNS:\n")
        f.write("-"*60 + "\n")
        for col in adata.var.columns:
            f.write(f"  {col}\n")
        
        f.write("\nOBS COLUMNS:\n")
        f.write("-"*60 + "\n")
        for col in adata.obs.columns:
            f.write(f"  {col}\n")
    
    # Save output
    logger.info(f"\nSaving feature-enriched dataset to {output_h5ad}")
    Path(output_h5ad).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad, compression='gzip')
    
    logger.info("\n" + "="*60)
    logger.info("FEATURE ENGINEERING COMPLETE")
    logger.info("="*60)
    logger.info(f"✓ Extracted features: {', '.join(extracted_features)}")
    logger.info("✓ Saved successfully")


if __name__ == "__main__":
    main()