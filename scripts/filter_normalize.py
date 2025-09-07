import scanpy as sc
import anndata
import yaml
import argparse

def load_config(path):
    with open(path, 'r') as stream:
        config = yaml.safe_load(stream)
    return config

def run_qc_filtering(adata, config):
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Basic filtering
    adata = adata[adata.obs['n_genes_by_counts'] >= config['qc']['min_genes'], :]
    adata = adata[adata.obs['pct_counts_mt'] <= config['qc']['max_mt_pct'], :]

    # Filter based on total counts (optional)
    if 'min_counts' in config['qc']:
        adata = adata[adata.obs['total_counts'] >= config['qc']['min_counts'], :]
    if 'max_counts' in config['qc']:
        adata = adata[adata.obs['total_counts'] <= config['qc']['max_counts'], :]

    # Filter based on mitochondrial counts (optional)
    if 'min_mt_pct' in config['qc']:
        adata = adata[adata.obs['pct_counts_mt'] >= config['qc']['min_mt_pct'], :]
    return adata

def run_normalization(adata, config):
    sc.pp.normalize_total(adata, target_sum=config['normalization']['target_sum'])
    sc.pp.log1p(adata)
    if config['normalization']['regress_vars']:
        sc.pp.regress_out(adata, keys=config['normalization']['regress_vars'])
    sc.pp.scale(adata)
    return adata
  
def run_doublet_detection(adata, config):
    import scrublet as scr
    if config['doublet_detection']['enable']:
        scrub = scr.Scrublet(adata.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
        threshold = config['doublet_detection']['score_threshold']
        adata = adata[adata.obs['doublet_score'] < threshold]
    return adata

def regress_cell_cycle(adata, config):
    if config['cell_cycle']['enable']:
        # Load gene lists for cell cycle from config or default sets
        s_genes = config['cell_cycle']['s_phase_genes']
        g2m_genes = config['cell_cycle']['g2m_phase_genes']

        # Score cell cycle phases
        sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

        # Regress out cell cycle score if enabled
        if config['cell_cycle']['regress']:
            sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
            sc.pp.scale(adata)
    return adata

def main(config_path):
    config = load_config(config_path)

    # Load data
    adata = sc.read_h5ad(config['data']['input_file'])
    print(f"Loaded data: {adata}")

    # QC filtering
    adata = run_qc_filtering(adata, config)
    print(f"After basic QC filtering: {adata}")

    # Normalization
    adata = run_normalization(adata, config)
    print(f"After normalization: {adata}")

    # Doublet detection and filtering
    adata = run_doublet_detection(adata, config)
    print(f"After doublet filtering: {adata}")

    # Cell cycle scoring and optional regression
    adata = regress_cell_cycle(adata, config)
    print(f"After cell cycle processing: {adata}")    

    # Save processed data
    adata.write_h5ad(config['data']['output_file'])
    print(f"Processed data saved to {config['data']['output_file']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter cells, normalize, detect doublets and regress cell cycle")
    parser.add_argument("--config", type=str, required=True, help="Path to config.yaml")
    args = parser.parse_args()
    main(args.config)
