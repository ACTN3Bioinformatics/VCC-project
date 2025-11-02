# VCC-project Pipeline Guide

Complete documentation for the single-cell CRISPR perturbation data processing pipeline.

------------------------------------------------------------------------

## Table of Contents

1.  [Overview](#overview)
2.  [Data Flow](#data-flow)
3.  [Pipeline Stages](#pipeline-stages)
4.  [Configuration](#configuration)
5.  [Hardware Requirements](#hardware-requirements)
6.  [Troubleshooting](#troubleshooting)
7.  [Advanced Usage](#advanced-usage)

------------------------------------------------------------------------

## Overview {#overview}

### Purpose

The VCC-project pipeline transforms raw single-cell RNA-seq data with CRISPR perturbations into high-quality, ML-ready datasets. It implements best practices for:

-   Quality control and filtering
-   Count normalization
-   Class balancing to prevent model bias
-   Batch effect correction
-   Feature engineering with biological annotations
-   Train/validation/test splitting

### Design Principles

1.  **Reproducibility**: All steps documented, versioned, and deterministic
2.  **Scalability**: Handles datasets from 10k to 1M+ cells
3.  **Modularity**: Each stage is independent and configurable
4.  **Efficiency**: Optimized for laptop and cluster environments
5.  **FAIR Compliance**: Findable, Accessible, Interoperable, Reusable

------------------------------------------------------------------------

## Data Flow {#data-flow}

```         
Raw Data (h5ad)
    ↓
[1. Quality Control]  ← Filter low-quality cells/genes
    ↓
[2. Normalization]    ← Count normalization + log transform
    ↓
[3. Scaling]          ← Zero-mean, unit variance
    ↓
[4. Balancing]        ← Downsample overrepresented perturbations
    ↓
[5. Integration*]     ← Batch correction (optional)
    ↓
[6. Feature Engineering] ← Add pathway/TF annotations
    ↓
[7. Splitting]        ← Train/val/test splits
    ↓
Final Dataset (h5ad)
```

\*Optional steps marked with asterisk

------------------------------------------------------------------------

## Pipeline Stages {#pipeline-stages}

### Stage 1: Data Acquisition

**Purpose**: Download and prepare demonstration datasets

**Input**: URL or GEO accession **Output**: `data_local/demo/replogle_subset.h5ad`

**What it does**: 1. Downloads Replogle et al. 2022 K562 essential Perturb-seq data 2. Subsets to \~10,000 cells and \~100-150 perturbations 3. Optimizes for 16GB RAM systems 4. Validates data integrity

**Snakemake rule**:

``` bash
snakemake download_demo_data --cores 1
```

**Manual execution**:

``` python
python scripts/download_demo_data.py
```

**Interactive exploration**: After download, explore the data with our demo notebook:

``` bash
jupyter notebook notebooks/demo_exploration.ipynb
```

The notebook provides: - Interactive data inspection - QC metrics visualization - Perturbation distribution analysis - UMAP and PCA plots - Gene expression patterns

**Configuration** (`config/config.yaml`):

``` yaml
demo_data_url: "https://huggingface.co/datasets/scEvalsJam/datasets/resolve/main/1gene-replogle-essential-split.h5ad"
demo_max_cells: 10000
demo_target_perturbations: 150
```

------------------------------------------------------------------------

### Stage 2: Quality Control

**Purpose**: Remove low-quality cells and uninformative genes

**Input**: Raw h5ad file **Output**: `results/{dataset}/filtered.h5ad`

**Quality metrics calculated**:

1.  **Per-cell metrics**:
    -   `n_genes_by_counts`: Number of genes detected
    -   `total_counts`: Total UMI counts
    -   `pct_counts_mt`: Mitochondrial gene percentage
    -   `doublet_score`: Doublet detection score (Scrublet)
2.  **Per-gene metrics**:
    -   `n_cells_by_counts`: Number of cells expressing gene
    -   `mean_counts`: Average expression level
    -   `pct_dropout_by_counts`: Percentage of zeros

**Filtering thresholds** (configurable per dataset):

``` yaml
demo:
  min_genes: 200          # Minimum genes per cell
  max_genes: 6000         # Maximum genes (doublet threshold)
  max_pct_mt: 15          # Maximum mitochondrial %
  min_cells_per_gene: 3   # Minimum cells per gene
```

**What gets removed**: - Dead/dying cells (high mt%, low gene count) - Doublets (very high gene count) - Empty droplets (very low UMI count) - Uninformative genes (expressed in \<3 cells)

**Snakemake rule**:

``` bash
snakemake results/demo/filtered.h5ad --cores 4
```

**Expected cell retention**: 70-85% of input cells

**QC report**: `reports/demo/qc_report.html`

------------------------------------------------------------------------

### Stage 3: Normalization

**Purpose**: Make counts comparable across cells with different sequencing depth

**Input**: Filtered h5ad **Output**: `results/{dataset}/normalized.h5ad`

**Steps**:

1.  **Total count normalization**:

    ``` python
    # Normalize each cell to same total count
    sc.pp.normalize_total(adata, target_sum=10000)
    ```

    -   Corrects for sequencing depth differences
    -   Default: 10,000 counts per cell
    -   Preserves relative expression within cells

2.  **Log transformation**:

    ``` python
    # Log1p transform for variance stabilization
    sc.pp.log1p(adata)
    ```

    -   Reduces impact of highly expressed genes
    -   Makes distribution more normal
    -   `log1p(x) = log(1 + x)` prevents log(0)

3.  **Regression** (optional):

    ``` python
    # Remove technical variation
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    ```

    -   Removes unwanted technical effects
    -   Use cautiously - can remove biological signal

**Configuration**:

``` yaml
demo:
  normalize: true
  target_sum: 10000
  log_transform: true
  regress_out:
    - "total_counts"
    - "pct_counts_mt"
```

**Snakemake rule**:

``` bash
snakemake results/demo/normalized.h5ad --cores 4
```

------------------------------------------------------------------------

### Stage 4: Scaling

**Purpose**: Standardize gene expression to zero mean and unit variance

**Input**: Normalized h5ad **Output**: `results/{dataset}/scaled.h5ad`

**What it does**:

``` python
# For each gene:
# scaled = (expression - mean) / std_dev
sc.pp.scale(adata, max_value=10)
```

**Why scale**: - Makes genes with different expression levels comparable - Required for PCA and many ML methods - Prevents high-expression genes from dominating

**max_value parameter**: - Clips extreme values to ±10 - Prevents outliers from distorting analysis - Recommended for scRNA-seq data

**Configuration**:

``` yaml
demo:
  scale: true
  max_value: 10
```

**Snakemake rule**:

``` bash
snakemake results/demo/scaled.h5ad --cores 4
```

------------------------------------------------------------------------

### Stage 5: Class Balancing

**Purpose**: Balance perturbation classes to prevent model bias

**Input**: Scaled h5ad **Output**: `results/{dataset}/balanced.h5ad`

**The problem**: - Some perturbations have 1000+ cells - Others have only 50 cells - ML models will bias toward overrepresented classes

**The solution - Downsampling**:

``` python
# For each perturbation:
target_cells = 100
if n_cells > target_cells:
    randomly_sample(target_cells)
else:
    keep_all()
```

**Configuration**:

``` yaml
demo:
  balance: true
  target_cells_per_perturbation: 100
  min_cells_per_perturbation: 10
  perturbation_key: "target_gene_name"
```

**Strategy by dataset**:

| Dataset    | Strategy         | Target Cells | Rationale                        |
|------------|------------------|--------------|----------------------------------|
| Demo       | Downsample       | 100          | Fair comparison, fast processing |
| Training   | Downsample       | 500          | More data for training           |
| Validation | Downsample       | 500          | Match training distribution      |
| Test       | **No balancing** | N/A          | Preserve real distribution       |

**Why not balance test set**: - Test should reflect real-world distribution - Balancing would artificially inflate performance - Want to see how model handles class imbalance

**Snakemake rule**:

``` bash
snakemake results/demo/balanced.h5ad --cores 2
```

**Output report**: `reports/demo/balance_report.txt`

------------------------------------------------------------------------

### Stage 6: Batch Integration (Optional)

**Purpose**: Remove technical batch effects while preserving biological variation

**Input**: Balanced h5ad **Output**: `results/{dataset}/integrated.h5ad`

**When to use**: - Data from multiple experiments - Different sequencing runs - Different labs/protocols

**When NOT to use**: - Single experiment (like demo) - Batch = biological condition of interest - Very small datasets

**Methods available**:

1.  **Harmony** (default, recommended):

    ``` python
    import harmonypy as hm
    ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    ```

    -   Fast, scalable
    -   Works well with large datasets
    -   Preserves biological structure

2.  **BBKNN** (alternative):

    ``` python
    import bbknn
    bbknn.bbknn(adata, batch_key='batch')
    ```

    -   Graph-based integration
    -   Good for small-medium datasets
    -   Can over-correct if not careful

**Configuration**:

``` yaml
demo:
  batch_correction: false  # Not needed for demo
  batch_key: "batch"
  integration_method: "harmony"

training:
  batch_correction: true   # Multiple experiments
  batch_key: "experiment_id"
  integration_method: "harmony"
```

**Snakemake rule**:

``` bash
snakemake results/training/integrated.h5ad --cores 8
```

**Validation**: - Check UMAP before/after integration - Verify biological signal preserved - Ensure batches are mixed

------------------------------------------------------------------------

### Stage 7: Feature Engineering

**Purpose**: Add biological annotations to enhance model performance

**Input**: Balanced or integrated h5ad **Output**: `results/{dataset}/features.h5ad`

**Features extracted**:

1.  **Pathway membership**:
    -   KEGG pathways
    -   Reactome pathways
    -   GO Biological Process
    -   GO Molecular Function
2.  **Transcription factor targets**:
    -   ChIP-seq data (ENCODE)
    -   Predicted TF binding sites
    -   Regulatory networks
3.  **H1-hESC specific features** (for VCC data):
    -   Pluripotency markers (OCT4, SOX2, NANOG)
    -   Differentiation markers
    -   Cell-cycle genes
4.  **Protein-protein interactions**:
    -   STRING database
    -   BioGRID interactions
    -   Pathway crosstalk

**Output format**:

``` python
# Added to adata.varm (gene-level features)
adata.varm['pathways']      # Binary matrix: gene x pathway
adata.varm['tf_targets']    # Binary matrix: gene x TF
adata.varm['ppi']           # Interaction scores

# Added to adata.uns (global annotations)
adata.uns['pathway_names']  # List of pathway names
adata.uns['tf_names']       # List of TF names
```

**Configuration**:

``` yaml
demo:
  extract_features: true
  use_pathways: true
  use_tf_targets: true
  use_regulatory_networks: false  # Too large for demo

training:
  extract_features: true
  use_pathways: true
  use_tf_targets: true
  use_regulatory_networks: true
```

**Databases used**: - KEGG: https://www.genome.jp/kegg/ - Reactome: https://reactome.org/ - GO: http://geneontology.org/ - ENCODE: https://www.encodeproject.org/

**Snakemake rule**:

``` bash
snakemake results/demo/features.h5ad --cores 4
```

**Memory note**: Feature engineering can be memory-intensive for large gene sets

------------------------------------------------------------------------

### Stage 8: Data Splitting

**Purpose**: Create train/validation/test splits for ML model development

**Input**: Feature-engineered h5ad **Output**: - `results/{dataset}/splits/train.h5ad` - `results/{dataset}/splits/val.h5ad` - `results/{dataset}/splits/test.h5ad` - `results/{dataset}/splits/split_info.json`

**Splitting strategy**:

**Leave-genes-out approach** (recommended for perturbation data):

``` python
# Split by perturbation, not by cells
# Ensures test genes are completely unseen
genes = unique_perturbations
train_genes = 70% of genes
val_genes = 15% of genes
test_genes = 15% of genes

# All cells for each gene go to same split
```

**Why this matters**: - Mimics real challenge: predicting unseen perturbations - Prevents data leakage - More realistic performance estimates

**Alternative: Random cell split** (not recommended):

``` python
# Would allow same gene in train and test
# Inflates performance metrics
# Doesn't test generalization
```

**Configuration**:

``` yaml
# In workflows/rules/split.smk
train_ratio: 0.7
val_ratio: 0.15
test_ratio: 0.15
gene_split: true  # Leave-genes-out
random_state: 42  # For reproducibility
```

**Split distribution example**:

```         
Dataset: 10,000 cells, 150 perturbations

Train:  7,000 cells (105 genes)
Val:    1,500 cells (23 genes)
Test:   1,500 cells (22 genes)
```

**Snakemake rule**:

``` bash
snakemake results/demo/splits/train.h5ad --cores 2
```

**Validation**: Check `split_info.json` for gene assignments

------------------------------------------------------------------------

### Stage 9: Final Assembly

**Purpose**: Combine all processed data into single h5ad file

**Input**: - Feature-engineered data - Train/val/test splits

**Output**: `results/{dataset}/final.h5ad`

**Final dataset structure**:

``` python
adata
├── X: scaled expression matrix
├── obs (cell metadata):
│   ├── n_genes_by_counts
│   ├── total_counts
│   ├── pct_counts_mt
│   ├── target_gene (perturbation)
│   ├── split (train/val/test)
│   └── batch (if applicable)
├── var (gene metadata):
│   ├── n_cells_by_counts
│   ├── mean_counts
│   └── highly_variable
├── varm (gene features):
│   ├── pathways
│   ├── tf_targets
│   └── ppi
├── uns (global metadata):
│   ├── pipeline_version
│   ├── processing_date
│   ├── config_used
│   └── split_info
└── layers:
    ├── counts (raw counts)
    ├── normalized (log-normalized)
    └── scaled (z-scored)
```

**Snakemake rule**:

``` bash
snakemake results/demo/final.h5ad --cores 2
```

------------------------------------------------------------------------

## Configuration {#configuration}

### Global Configuration (`config/config.yaml`)

``` yaml
# Data sources
demo_data_url: "https://huggingface.co/datasets/scEvalsJam/datasets/resolve/main/1gene-replogle-essential-split.h5ad"
demo_max_cells: 10000
demo_target_perturbations: 150

# Compute resources (AMD Ryzen 5 7535HS, 16GB RAM)
max_threads: 8
default_threads: 4
max_memory_mb: 14000

# Random seeds
random_seed: 42

# Feature engineering
pathway_databases:
  - "KEGG"
  - "Reactome"
  - "GO_Biological_Process"

tf_databases:
  - "ENCODE"
  - "ChIP-Atlas"
```

### Dataset-Specific Configuration (`config/datasets.yaml`)

``` yaml
demo:
  # Paths
  input_path: "data_local/demo/replogle_subset.h5ad"
  output_dir: "results/demo"
  
  # QC
  min_genes: 200
  max_genes: 6000
  max_pct_mt: 15
  min_cells_per_gene: 3
  
  # Normalization
  normalize: true
  target_sum: 10000
  log_transform: true
  scale: true
  
  # Balancing
  balance: true
  target_cells_per_perturbation: 100
  
  # Integration
  batch_correction: false
  
  # Features
  extract_features: true
```

### Customizing for Your Data

1.  **Add new dataset**:

``` yaml
my_dataset:
  input_path: "data_local/raw/my_data.h5ad"
  output_dir: "results/my_dataset"
  # ... copy other params from demo
```

2.  **Run pipeline**:

``` bash
snakemake results/my_dataset/final.h5ad --cores 4
```

------------------------------------------------------------------------

## Hardware Requirements {#hardware-requirements}

### Demo Data (10k cells, 150 perturbations)

**Minimum**: - CPU: 4 cores - RAM: 8GB - Storage: 20GB - Time: \~30 minutes

**Recommended** (AMD Ryzen 5 7535HS): - CPU: 8 cores \@ 3.55 GHz - RAM: 16GB LPDDR5x-6400 - GPU: AMD Radeon 660M (optional) - Storage: 50GB SSD - Time: \~15 minutes

### Full VCC 2025 Dataset (300k cells, 7000+ perturbations)

**Minimum**: - CPU: 16 cores - RAM: 64GB - Storage: 200GB - Time: \~4 hours

**Recommended** (HPC cluster): - CPU: 32+ cores - RAM: 128GB+ - Storage: 500GB SSD - Time: \~1-2 hours

### Memory Optimization Tips

1.  **Use backed mode** for very large files:

``` python
adata = sc.read_h5ad(file, backed='r')
```

2.  **Process in chunks**:

``` python
# In Snakemake rules
resources:
    mem_mb = lambda wildcards, attempt: 8000 * attempt
```

3.  **Reduce subset size**:

``` yaml
demo_max_cells: 5000  # Instead of 10000
```

------------------------------------------------------------------------

## Troubleshooting {#troubleshooting}

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for detailed solutions.

**Common issues**:

1.  **Out of memory**: Reduce `demo_max_cells` or `max_genes`
2.  **Download fails**: Check internet connection, use manual download
3.  **Perturbation column not found**: Check `perturbation_key` in config
4.  **Pipeline stuck**: Check logs in `logs/`

------------------------------------------------------------------------

## Advanced Usage {#advanced-usage}

### Running Specific Stages

``` bash
# Only QC
snakemake results/demo/filtered.h5ad --cores 4

# Through normalization
snakemake results/demo/normalized.h5ad --cores 4

# Skip feature engineering
snakemake results/demo/balanced.h5ad --cores 4
```

### Parallel Processing

``` bash
# Use all cores
snakemake --cores all

# Specific number
snakemake --cores 8
```

### Cluster Execution

``` bash
# Submit to SLURM
snakemake --cluster "sbatch -n {threads} --mem {resources.mem_mb}" \
          --jobs 10

# With custom profile
snakemake --profile slurm/
```

### Custom Rules

Add to `workflows/rules/custom.smk`:

``` python
rule my_custom_analysis:
    input:
        "results/demo/final.h5ad"
    output:
        "results/demo/my_analysis.pdf"
    script:
        "../../scripts/my_analysis.py"
```

------------------------------------------------------------------------

## Interactive Data Exploration

### Jupyter Notebooks

The project includes interactive Jupyter notebooks for hands-on data exploration and visualization.

#### Demo Exploration Notebook

**Location**: `notebooks/demo_exploration.ipynb`

**Purpose**: Comprehensive introduction to the demo dataset with interactive visualizations.

**What's included**:

1.  **Setup and Data Loading**
    -   Environment configuration
    -   Loading processed data
    -   Inspecting data structure
2.  **Quality Control Visualization**
    -   Genes per cell distributions
    -   UMI count distributions
    -   Mitochondrial content analysis
    -   QC threshold visualization
3.  **Perturbation Analysis**
    -   Perturbation distribution plots
    -   Cell counts per perturbation
    -   Top perturbations identification
    -   Imbalance visualization
4.  **Dimensionality Reduction**
    -   PCA variance explained
    -   UMAP visualization
    -   Coloring by QC metrics
    -   Highlighting specific perturbations
5.  **Gene Expression Analysis**
    -   Marker gene expression
    -   Expression patterns on UMAP
    -   Control vs perturbed cells
    -   Differential expression visualization
6.  **Pipeline Output Comparison**
    -   Comparing filtered vs balanced data
    -   Tracking cell/gene counts through pipeline
    -   File size comparisons

**Launch notebook**:

``` bash
# Start Jupyter and open notebook
jupyter notebook notebooks/demo_exploration.ipynb

# Or launch Jupyter Lab
jupyter lab
```

**Learning outcomes**: After working through the notebook, you'll understand: - How to load and inspect h5ad files - Quality control metrics and their interpretation - Perturbation class imbalance issues - Dimensionality reduction for visualization - How pipeline stages transform the data

#### Creating Custom Notebooks

Want to create your own analysis notebooks?

1.  **Copy the template**:

``` bash
cp notebooks/demo_exploration.ipynb notebooks/my_analysis.ipynb
```

2.  **Customize for your data**:

``` python
# Load your processed data
adata = sc.read_h5ad('results/your_dataset/final.h5ad')
```

3.  **Share your notebook**:

-   Notebooks are tracked in Git (not in `.gitignore`)
-   Great for reproducible analyses
-   Include in publications as supplementary material

**Best practices for notebooks**: - Keep cells focused and well-commented - Use markdown cells for explanations - Save outputs (figures) for documentation - Clear outputs before committing to Git - Use version control for notebook history

------------------------------------------------------------------------

## Best Practices

### For Demo/Testing

1.  Use small subset (5-10k cells)
2.  Disable expensive steps (integration, extensive features)
3.  Monitor memory usage
4.  Save intermediate results

### For Production

1.  Full QC with multiple thresholds
2.  Extensive feature engineering
3.  Cross-validation with multiple folds
4.  Document all parameters
5.  Version control configs

### For Reproducibility

1.  Fix all random seeds
2.  Save exact package versions
3.  Archive raw data with checksums
4.  Document hardware used
5.  Share configs in publications

------------------------------------------------------------------------

## Citation

If you use this pipeline, please cite:

``` bibtex
@software{vcc_project_2025,
  author = {Szymon Myrta},
  title = {VCC-project: Single-Cell CRISPR Perturbation Pipeline},
  year = {2025},
  url = {https://github.com/ACTN3Bioinformatics/VCC-project},
  doi = {10.5281/zenodo.17503792}
}
```

**Demo data citation**:

``` bibtex
@article{replogle2022mapping,
  title={Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq},
  author={Replogle, Joseph M and Saunders, Reuben A and others},
  journal={Cell},
  volume={185},
  number={14},
  pages={2559--2575},
  year={2022}
}
```

------------------------------------------------------------------------

## Support

-   **Issues**: https://github.com/ACTN3Bioinformatics/VCC-project/issues
-   **Discussions**: https://github.com/ACTN3Bioinformatics/VCC-project/discussions
-   **Email**: your.email\@example.com

------------------------------------------------------------------------

**Last updated**: October 27, 2025\
**Pipeline version**: 1.0.0