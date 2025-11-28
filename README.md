# VCC-project: Single-Cell CRISPR Perturbation Pipeline

[![CI](https://github.com/ACTN3Bioinformatics/VCC-project/actions/workflows/ci.yml/badge.svg)](https://github.com/ACTN3Bioinformatics/VCC-project/actions/workflows/ci.yml) [![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Snakemake](https://img.shields.io/badge/snakemake-≥7.0-brightgreen.svg)](https://snakemake.readthedocs.io) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17508488.svg)](https://doi.org/10.5281/zenodo.17508488)

A reproducible computational pipeline for processing and analyzing single-cell RNA-seq data with CRISPR perturbations, designed for the Virtual Cell Challenge 2025.

## 🎯 Overview

This pipeline transforms raw single-cell RNA-seq data with genetic perturbations into high-quality, balanced datasets suitable for training ML/AI models to predict perturbation effects. It supports multiple dataset types and provides comprehensive quality control, normalization, and integration capabilities.

### Key Features

-   **🔄 Automated Workflow**: Snakemake-based pipeline with dependency management
-   **✅ Quality Control**: Comprehensive filtering of low-quality cells and genes
-   **⚖️ Class Balancing**: Smart downsampling to prevent model bias
-   **🔗 Batch Integration**: Harmonization of multiple datasets using Harmony/BBKNN
-   **🧬 Feature Engineering**: Biological feature extraction (pathways, TFs, regulatory networks)
-   **📊 Cross-Validation**: Leave-genes-out CV strategy for unseen perturbations
-   **🔬 Reproducibility**: Fully containerized environment with version control
-   **📓 Interactive Notebooks**: Jupyter notebooks for data exploration and visualization

## 📋 Quick Start

### Prerequisites

-   Python 3.9+
-   Conda/Mamba
-   **Hardware**: 4+ cores, 8GB+ RAM (16GB recommended)
-   **Storage**: \~20GB for demo data, \~100GB for full datasets

### Installation

``` bash
# Clone repository
git clone https://github.com/ACTN3Bioinformatics/VCC-project.git
cd VCC-project

# Create mamba environment
mamba env create -f environment.yml
mamba activate vcc2025
```

### Demo Data Setup

Download demonstration data (optimized subset of Replogle et al. 2022):

``` bash
# Automatic download and preparation
snakemake download_demo_data --cores 1

# This creates: data_local/demo/replogle_subset.h5ad (~500MB, 10k cells)
```

### Running the Pipeline

``` bash
# Run complete pipeline on demo data
snakemake --cores 4 --configfile config/datasets.yaml

# Run specific stages
snakemake results/demo/filtered.h5ad --cores 4          # QC only
snakemake results/demo/balanced.h5ad --cores 4          # Through balancing
snakemake results/demo/final.h5ad --cores 4             # Complete pipeline

# Dry-run to see execution plan
snakemake -n

# Test configuration
snakemake test --cores 1

# Generate workflow visualization
python scripts/generate_workflow_diagram.py
# Creates: docs/workflow_diagram.png (if graphviz installed)
#      or: docs/workflow_dag.txt (always works)
```

### Explore with Jupyter Notebook

``` bash
# Launch demo exploration notebook
jupyter notebook notebooks/demo_exploration.ipynb

# Or explore processed results
jupyter notebook
```

## 📁 Project Structure

```         
VCC-project/
├── workflows/              # Snakemake workflow definitions
│   ├── Snakefile          # Main workflow entry point
│   └── rules/             # Individual pipeline rules
│       ├── download.smk   # Data acquisition
│       ├── qc.smk         # Quality control & filtering
│       ├── normalize.smk  # Normalization & scaling
│       ├── balance.smk    # Class balancing
│       ├── integrate.smk  # Batch integration
│       ├── split.smk      # Train/val/test splits
│       └── features.smk   # Feature engineering
├── scripts/               # Core Python modules
│   ├── download_demo_data.py
│   ├── filter_normalize.py
│   ├── balance.py
│   ├── integration.py
│   ├── split_data.py
│   ├── feature_engineering.py
│   └── utils.py
├── config/                # Configuration files
│   ├── datasets.yaml      # Dataset-specific parameters
│   └── config.yaml        # Global pipeline settings
├── data_local/            # Local data storage (NOT tracked in Git)
│   ├── demo/             # Demonstration datasets
│   ├── raw/              # Raw input data
│   └── processed/        # Intermediate outputs
├── results/               # Final pipeline outputs
│   └── demo/             # Demo results
├── reports/               # QC reports and visualizations
├── logs/                  # Snakemake and script logs
├── notebooks/             # Jupyter notebooks
│   └── demo_exploration.ipynb  # Interactive demo notebook
├── docs/                  # Extended documentation
│   ├── PIPELINE_GUIDE.md  # Detailed pipeline guide
│   ├── QUICKSTART.md     # 5-minute tutorial
│   └── TROUBLESHOOTING.md # Common issues
├── tests/                 # Unit tests
│   ├── test_qc.py
│   └── test_balance.py
├── environment.yml        # Conda environment specification
├── LICENSE               # MIT License
├── CITATION.cff          # Citation metadata
├── CONTRIBUTING.md       # Contribution guidelines
└── README.md             # This file
```

## 🔬 Pipeline Overview

The pipeline consists of modular stages executed by Snakemake:

1.  **📥 Data Acquisition** - Download and prepare demo data
2.  **🔍 Quality Control** - Filter low-quality cells and genes
3.  **📊 Normalization** - Count normalization and log transformation
4.  **⚖️ Class Balancing** - Balance perturbation classes
5.  **🔗 Batch Integration** - Harmonize datasets (optional)
6.  **🧬 Feature Engineering** - Extract biological features
7.  **✂️ Data Splitting** - Create train/val/test splits
8.  **📈 Benchmarking** - Evaluate baseline models

**For detailed information**, see [docs/PIPELINE_GUIDE.md](docs/PIPELINE_GUIDE.md).

## 📊 Dataset Types

| Dataset | Description | Size | Purpose | Processing |
|---------------|---------------|---------------|---------------|---------------|
| **Demo** | Replogle K562 subset | \~10k cells | Testing/Learning | Full pipeline |
| **Training** | H1-hESC CRISPRi | \~300k cells | Model training | Full QC + balancing |
| **Validation** | H1-hESC validation | \~50k cells | Model selection | Same as training |
| **Test** | Unseen perturbations | \~50k cells | Final evaluation | Minimal processing |
| **Public** | External datasets | Variable | Pre-training/augmentation | Full integration |

## 🔧 Configuration

Customize processing via `config/datasets.yaml`:

``` yaml
demo:
  input_path: "data_local/demo/replogle_subset.h5ad"
  output_dir: "results/demo"
  
  # QC thresholds
  min_genes: 200
  max_genes: 6000
  max_pct_mt: 15
  min_cells_per_gene: 3
  
  # Processing options
  normalize: true
  log_transform: true
  scale: true
  balance: true
  target_cells_per_perturbation: 100
  
  # Integration (for multi-batch data)
  batch_correction: false
  batch_key: "batch"
```

See [docs/PIPELINE_GUIDE.md#configuration](docs/PIPELINE_GUIDE.md#configuration) for all options.

## 💻 System Requirements

### Minimum (Demo Data)

-   **CPU**: 4 cores
-   **RAM**: 8GB
-   **Storage**: 20GB SSD
-   **Time**: \~30 minutes

### Recommended (Demo Data)

-   **CPU**: AMD Ryzen 5 7535HS or equivalent (8 cores \@ 3.55 GHz)
-   **RAM**: 16GB LPDDR5x-6400
-   **GPU**: AMD Radeon 660M (optional, for ML training)
-   **Storage**: 50GB SSD
-   **Time**: \~15 minutes

### Full VCC 2025 Dataset

-   **CPU**: 16+ cores
-   **RAM**: 64GB+
-   **GPU**: 16GB+ VRAM for deep learning
-   **Storage**: 200GB+ SSD
-   **Time**: \~2-4 hours

**Note**: Demo data is specifically optimized for laptop processing on AMD Ryzen 5 7535HS system (16GB RAM).

## 📚 Documentation

-   [**Quick Start Guide**](docs/QUICKSTART.md) - Get running in 5 minutes
-   [**Pipeline Guide**](docs/PIPELINE_GUIDE.md) - Complete pipeline documentation
-   [**Troubleshooting**](docs/TROUBLESHOOTING.md) - Common issues and solutions
-   [**Demo Notebook**](notebooks/demo_exploration.ipynb) - Interactive demo data exploration
-   [**Demo Report**](docs/demo_report.html) - HTML report for demo data
-   [**Contributing Guide**](CONTRIBUTING.md) - How to contribute

## 📓 Jupyter Notebooks

### Demo Exploration Notebook

The `notebooks/demo_exploration.ipynb` provides an interactive introduction to:

-   Loading and inspecting processed data
-   Visualizing QC metrics
-   Exploring perturbation effects
-   Dimensionality reduction (PCA, UMAP)
-   Comparing pipeline stages

**Launch notebook**:

``` bash
jupyter notebook notebooks/demo_exploration.ipynb
```

## 🧪 Testing

Run unit tests:

``` bash
# All tests
pytest tests/

# Specific test
pytest tests/test_qc.py -v

# With coverage
pytest --cov=scripts tests/

# Test CI locally
./test_ci_locally.sh
```

## 🤝 Contributing

Contributions welcome! Please:

1.  Fork the repository
2.  Create a feature branch (`git checkout -b feature/amazing-feature`)
3.  Commit changes (`git commit -m 'Add amazing feature'`)
4.  Push to branch (`git push origin feature/amazing-feature`)
5.  Open a Pull Request

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

## 👤 Author

### Szymon Myrta (kontakt@actn3.pl)

I'm passionate about applying data science and bioinformatics to advance pharmaceutical research and precision medicine. This repository reflects my commitment to continuous learning, knowledge sharing, and contributing to the open-source R community in pharma.

-   8+ years of experience in bioinformatics in pharmaceutical and biotech settings
-   NGS data analysis (RNAseq, scRNA-seq, ChIP-seq, TCR/BCR-seq, WES/WGS etc.) \[Bioconductor, Seurat, edgeR, DESeq2, limma\]
-   Functional genomics data analysis (CRISPR / ORF overexpression screens)
-   Strong R programming skills, including development of R packages and Shiny apps
-   Developer of NGS data analysis pipelines and reproducible research workflows and documentation (RMarkdown, Quarto, Snakemake, Bash, Git, CI/CD, Docker)
-   Data visualization and interpretation of results \[ggplot2, Shiny, ComplexHeatmap\]
-   Background in computational biology, cancer genomics, immuno-oncology
-   Co-author of multiple peer-reviewed scientific publications in top-tier journals
-   **Interested in multi-omics data integration, precision medicine, AI-powered analyses**

## 📝 Citation

If you use this pipeline in your research, please cite:

``` bibtex
@software{vcc_project_2025,
  author = {Szymon Myrta},
  title = {VCC-project: Single-Cell CRISPR Perturbation Pipeline},
  year = {2025},
  url = {https://github.com/ACTN3Bioinformatics/VCC-project},
  doi = {10.5281/zenodo.17503792}
}
```

**Demo data**: If using the demo dataset, please also cite: - Replogle et al. (2022). "Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq." Cell. DOI: [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

-   Virtual Cell Challenge 2025 organizers
-   Replogle et al. for public Perturb-seq data
-   scPerturb database for curated datasets
-   Scanpy and AnnData developers
-   Snakemake community

## 📮 Contact

-   **Issues**: [GitHub Issues](https://github.com/ACTN3Bioinformatics/VCC-project/issues)
-   **Discussions**: [GitHub Discussions](https://github.com/ACTN3Bioinformatics/VCC-project/discussions)
-   **Email**: kontakt\@actn3.pl

------------------------------------------------------------------------

**Status**: 🚀 Active Development \| **Version**: 1.0.0 \| **Last Updated**: November 2025