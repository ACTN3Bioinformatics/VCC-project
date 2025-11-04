"""
VCC-project Main Snakemake Workflow
Single-cell CRISPR perturbation data processing pipeline

IMPORTANT: This file should be in the ROOT directory of the project
Run from root: snakemake --cores 4
"""

import os
from pathlib import Path
import yaml

# ============================================================================
# Configuration - Use absolute paths relative to Snakefile location
# ============================================================================

# Get the directory where this Snakefile is located
SNAKEFILE_DIR = Path(workflow.basedir)

configfile: str(SNAKEFILE_DIR / "config" / "config.yaml")

# Load datasets configuration
datasets_path = SNAKEFILE_DIR / "config" / "datasets.yaml"
with open(datasets_path, 'r') as f:
    datasets = yaml.safe_load(f)

# Merge into main config
config['datasets'] = datasets

# Get list of datasets to process (only those with existing input files)
ALL_DATASETS = list(config['datasets'].keys())

# Filter to only datasets where input file exists
DATASETS = []
for dataset in ALL_DATASETS:
    input_path = config['datasets'][dataset].get('input_path')
    if input_path and Path(input_path).exists():
        DATASETS.append(dataset)

if not DATASETS:
    # If no datasets have existing input files, include 'demo' anyway for download
    if 'demo' in ALL_DATASETS:
        DATASETS = ['demo']

print(f"Available datasets: {ALL_DATASETS}")
print(f"Processing datasets: {DATASETS}")

# ============================================================================
# Global variables
# ============================================================================

DATA_LOCAL = config.get("data_local_dir", "data_local")
RESULTS = config.get("results_dir", "results")
REPORTS = config.get("reports_dir", "reports")
LOGS = config.get("logs_dir", "logs")

# ============================================================================
# Include rules - use absolute paths
# ============================================================================

include: str(SNAKEFILE_DIR / "workflows" / "rules" / "download.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "qc.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "normalize.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "balance.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "integrate.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "split.smk")
include: str(SNAKEFILE_DIR / "workflows" / "rules" / "features.smk")

# ============================================================================
# Target rules
# ============================================================================

rule all:
    """
    Default target: process all configured datasets through complete pipeline
    """
    input:
        # Final outputs for each dataset
        expand("results/{dataset}/final.h5ad", dataset=DATASETS),
        # QC reports
        expand("reports/{dataset}/qc_report.html", dataset=DATASETS)


rule download_all:
    """
    Download all demo datasets
    """
    input:
        expand("{data_local}/demo/replogle_subset.h5ad", data_local=DATA_LOCAL)


rule qc_all:
    """
    Run QC on all datasets
    """
    input:
        expand("results/{dataset}/filtered.h5ad", dataset=DATASETS)


rule balance_all:
    """
    Balance all datasets
    """
    input:
        expand("results/{dataset}/balanced.h5ad", dataset=DATASETS)


rule features_all:
    """
    Extract features for all datasets
    """
    input:
        expand("results/{dataset}/features.h5ad", dataset=DATASETS)


# ============================================================================
# Utility rules
# ============================================================================

rule clean:
    """
    Remove all generated files (use with caution!)
    """
    shell:
        """
        rm -rf results/*
        rm -rf reports/*
        rm -rf logs/*
        rm -rf .snakemake/*
        echo "✓ Cleaned all results and logs"
        """


rule clean_temp:
    """
    Remove temporary files only
    """
    shell:
        """
        find results -name "*.tmp" -delete
        find {DATA_LOCAL} -name "*.tmp" -delete
        echo "✓ Cleaned temporary files"
        """


rule test:
    """
    Test rule to verify setup
    """
    run:
        print("="*60)
        print("VCC-PROJECT CONFIGURATION TEST")
        print("="*60)
        print(f"Snakefile directory: {SNAKEFILE_DIR}")
        print(f"Working directory: {os.getcwd()}")
        print(f"Available datasets: {ALL_DATASETS}")
        print(f"Datasets to process: {DATASETS}")
        print(f"Data directory: {DATA_LOCAL}")
        print(f"Results directory: {RESULTS}")
        print("="*60)
        
        # Check if files exist
        for dataset in ALL_DATASETS:
            input_path = config['datasets'][dataset].get('input_path')
            exists = "✓" if Path(input_path).exists() else "✗"
            print(f"{exists} {dataset}: {input_path}")
        
        print("="*60)
