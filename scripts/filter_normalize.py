import scanpy as sc
import argparse
import yaml

def load_config(path: str) -> dict:
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def filter_and_normalize(adata, config):
    # Load QC params
    qc = config.get('qc', {})
    min_genes = qc.get('min_genes', 200)
    max_mt_pct = qc.get('max_mt_pct', 15)

    # Calculate mitochondrial metrics and overall QC
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Filter cells
    adata = adata[(adata.obs['n_genes_by_counts'] >= min_genes) &
                  (adata.obs['pct_counts_mt'] <= max_mt_pct)].copy()

    # Doublet detection (optional)
    if config.get('doublet_detection', {}).get('enabled', True):
        import scrublet as scr
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        adata = adata[~adata.obs['predicted_doublet']].copy()

    # Normalize and log transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Cell cycle adjustment: none/regress/covariate
    cell_cycle_cfg = config.get('cell_cycle', {})
    mode = cell_cycle_cfg.get('mode', 'none')

    # Load provided gene lists or fallback to Scanpy defaults
    s_genes = cell_cycle_cfg.get('s_phase_genes')
    g2m_genes = cell_cycle_cfg.get('g2m_phase_genes')

    if s_genes is None or g2m_genes is None:
        import scanpy as sc_inner
        # Use Scanpy default cell cycle genes
        from scanpy.preprocess._utils import _get_default_s_genes, _get_default_g2m_genes
        s_genes = s_genes if s_genes else _get_default_s_genes()
        g2m_genes = g2m_genes if g2m_genes else _get_default_g2m_genes()

    if mode != 'none':
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        if mode == 'regress':
            # regress out cell cycle scores then scale
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
            if config.get('scaling', {}).get('enabled', True):
                sc.pp.scale(adata)

        elif mode == 'covariate':
            # Do not regress, keep scores as covariates (no scaling here)
            pass
    else:
        # Default scaling if enabled and no cell cycle regression
        if config.get('scaling', {}).get('enabled', True):
            sc.pp.scale(adata)

    return adata

def main():
    parser = argparse.ArgumentParser(description="Filter and normalize VCC dataset with cell cycle options")
    parser.add_argument("--dataset_name", required=True, help="Name of dataset")
    parser.add_argument("--config", required=True, help="Path to YAML config with parameters")
    args = parser.parse_args()

    config_all = load_config(args.config)
    dataset_cfg = config_all['datasets'][args.dataset_name]

    adata = sc.read_h5ad(dataset_cfg['input_file'])
    adata_filtered = filter_and_normalize(adata, dataset_cfg)

    output_file = dataset_cfg.get('output_file', f"output/{args.dataset_name}_filtered.h5ad")
    adata_filtered.write_h5ad(output_file)
    print(f"Filtered and normalized data saved to: {output_file}")

if __name__ == "__main__":
    main()
