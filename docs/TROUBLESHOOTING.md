# Troubleshooting Guide

Common issues and quick solutions for VCC-project.

---

## Installation Issues

### Conda environment fails
```bash
# Try mamba (faster)
conda install mamba -c conda-forge
mamba env create -f environment.yml

# Or update conda
conda update -n base conda
conda env create -f environment.yml
```

### Snakemake not found
```bash
conda activate vcc2025
conda install -c conda-forge -c bioconda snakemake
```

---

## Data Download Issues

### Download timeout
```bash
# Use wget with resume
wget -c https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad
mv ReplogleWeissman2022_K562_essential.h5ad data_local/demo/

# Or use backup URL
wget -c https://plus.figshare.com/ndownloader/files/42444534
```

### Demo data too large
```yaml
# Edit config/config.yaml
demo_max_cells: 5000  # Reduce from 10000
demo_target_perturbations: 75
```

---

## Memory Issues

### Out of memory errors
**Solution 1** - Reduce data size (config/datasets.yaml):
```yaml
demo:
  max_genes: 4000  # Lower from 6000
  target_cells_per_perturbation: 50
```

**Solution 2** - Reduce threads:
```bash
snakemake --cores 2  # Instead of 4
```

**Solution 3** - Use backed mode for large files:
```python
adata = sc.read_h5ad(file, backed='r')
```

---

## Pipeline Execution Issues

### Snakemake rules fail
```bash
# Check logs
cat logs/qc/demo_filter.log

# Run with verbose output
snakemake --cores 4 -p

# Unlock if locked
snakemake --unlock
```

### Missing input files
```bash
# Re-download demo data
snakemake download_demo_data --cores 1 --forcerun

# Check config paths
cat config/datasets.yaml
```

### Import errors
```bash
# Verify environment
conda activate vcc2025
python -c "import scanpy, snakemake, pandas"

# Reinstall if needed
conda env update -f environment.yml --prune
```

---

## Data Processing Issues

### h5ad file corrupted
```python
import h5py

# Check integrity
with h5py.File('file.h5ad', 'r') as f:
    print(f.keys())
```

### Perturbation column not found
```python
import scanpy as sc

# Check available columns
adata = sc.read_h5ad('data.h5ad')
print(adata.obs.columns)

# Update config with correct column name
# perturbation_key: "gene" or "target_gene_name"
```

---

## Performance Issues

### Pipeline too slow
**Speed up**:
```bash
# Use more cores
snakemake --cores 8

# Process specific target
snakemake results/demo/balanced.h5ad --cores 4

# Disable expensive steps
# In config/datasets.yaml:
#   extract_features: false
#   batch_correction: false
```

---

## AMD Ryzen 5 7535HS Specific

### GPU not utilized
```bash
# Check ROCm support
rocminfo

# For PyTorch with ROCm
pip install torch --index-url https://download.pytorch.org/whl/rocm5.4.2
```

### Thermal throttling
```bash
# Monitor temperatures
watch -n 1 sensors

# Reduce threads
snakemake --cores 4  # Instead of 8
```

---

## File System Issues

### data_local tracked by Git
```bash
# Remove from tracking
git rm -r --cached data_local/
git commit -m "Remove data_local from tracking"

# Ensure .gitignore is correct
echo "data_local/" >> .gitignore
```

### Disk space full
```bash
# Clean caches
rm -rf .snakemake/
conda clean --all

# Remove old results
rm -rf results/old_*
```

---

## Getting Help

1. **Check logs**: `cat logs/[step]/[dataset].log`
2. **Dry run**: `snakemake -n` shows execution plan
3. **GitHub Issues**: https://github.com/ACTN3Bioinformatics/VCC-project/issues
4. **Snakemake docs**: https://snakemake.readthedocs.io/
5. **Scanpy docs**: https://scanpy.readthedocs.io/

---

## Reporting Bugs

Include in your report:
- Error message and traceback
- Log file contents
- System: `uname -a`
- Python: `python --version`
- Packages: `conda list | grep -E "(scanpy|snakemake)"`
- Config: `cat config/datasets.yaml`

---

**Last updated**: October 27, 2025