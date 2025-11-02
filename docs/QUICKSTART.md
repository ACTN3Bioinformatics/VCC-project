# Quick Start Guide

Get up and running with VCC-project in 5 minutes!

## Prerequisites
```bash
# Check Python version (need 3.9+)
python --version

# Check conda
conda --version
```

## Installation (2 minutes)
```bash
# 1. Clone repository
git clone https://github.com/ACTN3Bioinformatics/VCC-project.git
cd VCC-project

# 2. Create environment
conda env create -f environment.yml

# 3. Activate environment
conda activate vcc2025
```

## Run Demo Pipeline (3 minutes)
```bash
# Download demo data
snakemake download_demo_data --cores 1

# Run complete pipeline
snakemake --cores 4 --configfile config/datasets.yaml

# View results
ls results/demo/
```

## What Just Happened?

The pipeline:
1. ✅ Downloaded ~10k cells from public Perturb-seq data
2. ✅ Filtered low-quality cells (QC)
3. ✅ Normalized and log-transformed counts
4. ✅ Balanced perturbation classes
5. ✅ Created train/val/test splits
6. ✅ Generated QC reports

## Explore Results
```python
import scanpy as sc

# Load processed data
adata = sc.read_h5ad('results/demo/final.h5ad')

# Check it out
print(adata)
print(adata.obs['split'].value_counts())
```

## Next Steps

- 📖 Read [PIPELINE_GUIDE.md](PIPELINE_GUIDE.md) for detailed documentation
- 📓 Explore data interactively: `jupyter notebook notebooks/demo_exploration.ipynb`
- 🔧 Customize parameters in `config/datasets.yaml`
- 📊 Check QC report: `reports/demo/qc_report.html`
- 🚀 Process your own data by adding to `config/datasets.yaml`

## Interactive Data Exploration

Want to visualize and explore your data? Launch the demo notebook:
```bash
# Start Jupyter
jupyter notebook notebooks/demo_exploration.ipynb

# The notebook covers:
# - Loading processed data
# - QC metrics visualization
# - Perturbation analysis
# - PCA and UMAP plots
# - Gene expression patterns
```

## Troubleshooting

**Out of memory?**
```yaml
# In config/datasets.yaml, reduce:
demo:
  max_genes: 4000
  target_cells_per_perturbation: 50
```

**Download failed?**
```bash
# Manual download
wget https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad
mv ReplogleWeissman2022_K562_essential.h5ad data_local/demo/replogle_subset.h5ad
```

Need more help? See [TROUBLESHOOTING.md](TROUBLESHOOTING.md)