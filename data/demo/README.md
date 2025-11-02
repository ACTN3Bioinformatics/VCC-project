# Demo Dataset

This directory contains demonstration data for VCC-project.

## Dataset Information

**Source**: Replogle et al. (2022)  
**Title**: Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq  
**DOI**: [10.1016/j.cell.2022.05.013](https://doi.org/10.1016/j.cell.2022.05.013)  
**GEO Accession**: GSE140023

## Contents

- `replogle_subset.h5ad` - Subset of K562 Perturb-seq data (~10k cells, ~150 perturbations)
- `metadata.json` - Dataset metadata and provenance

## Download

Data is automatically downloaded when running:
```bash
snakemake download_demo_data --cores 1
```

Or manually from:
- **Primary**: https://zenodo.org/records/7041849/files/ReplogleWeissman2022_K562_essential.h5ad
- **Backup**: https://plus.figshare.com/ndownloader/files/42444534

## Data Characteristics

- **Cell type**: K562 (chronic myelogenous leukemia cell line)
- **Technology**: Perturb-seq (CRISPR perturbations + scRNA-seq)
- **Perturbation type**: CRISPRi knockdowns
- **Size**: ~500MB (subsetted from 1.44 GB original)
- **Cells**: ~10,000 (subset from ~188k)
- **Genes**: ~18,000
- **Perturbations**: ~150 genes (subset from ~2000)

## Citation

If you use this demo data, please cite:
```bibtex
@article{replogle2022mapping,
  title={Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq},
  author={Replogle, Joseph M and Saunders, Reuben A and others},
  journal={Cell},
  volume={185},
  number={14},
  pages={2559--2575},
  year={2022},
  publisher={Elsevier}
}
```

## Notes

- This is a **subset** optimized for demonstration and testing
- Data has been **downsampled** to be processable on laptops (16GB RAM)
- Not suitable for production model training - use full VCC 2025 datasets instead
- For full dataset, visit [GEO GSE140023](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE140023)