import scanpy as sc
import argparse
import yaml
import numpy as np

def load_config(path: str) -> dict:
    """
    Load YAML configuration file.
    
    Parameters:
    - path (str): Path to the YAML configuration file.
    
    Returns:
    - dict: Loaded configuration dictionary.
    """
    with open(path, 'r') as stream:
        return yaml.safe_load(stream)

def balance_classes(
    adata: sc.AnnData,
    target_key: str,
    mode: str = 'min_class',
    target_count: int = None,
    random_state: int = 42
) -> sc.AnnData:
    """
    Balance perturbation classes to avoid class imbalance bias.
    Performs downsampling of larger classes to either the minimum class size or a fixed target count.
    
    Parameters:
    - adata (AnnData): Annotated data matrix with cell metadata.
    - target_key (str): Column name in adata.obs containing class labels.
    - mode (str): Balancing strategy ('min_class' or 'fixed').
    - target_count (int): Target cell count per class when mode='fixed'.
    - random_state (int): Seed for reproducibility.
    
    Returns:
    - AnnData: Subset AnnData object with balanced classes.
    """
    np.random.seed(random_state)

    # Print class counts before balancing
    counts = adata.obs[target_key].value_counts()
    print("Sample counts per class before balancing:")
    print(counts)

    # Determine downsampling target count
    if mode == 'min_class':
        target_count = counts.min()
        print(f"Balancing all classes to minimum class size: {target_count}")
    elif mode == 'fixed':
        if target_count is None:
            raise ValueError("In 'fixed' mode, target_count must be provided.")
        print(f"Balancing all classes to fixed size: {target_count}")
    elif mode == 'none':
        print("Balancing disabled for this dataset.")
        return adata
    else:
        raise ValueError(f"Unknown mode '{mode}'. Supported modes: 'min_class', 'fixed', 'none'.")

    selected_indices = []

    # For each perturbation class, sample or keep all depending on class size
    for cls, group in adata.obs.groupby(target_key):
        n_samples = len(group)
        if n_samples > target_count:
            sampled = group.sample(n=target_count, random_state=random_state).index
        else:
            sampled = group.index
        selected_indices.extend(sampled)

    # Return balanced subset
    balanced_adata = adata[selected_indices].copy()

    # Print class counts after balancing
    counts_after = balanced_adata.obs[target_key].value_counts()
    print("Sample counts per class after balancing:")
    print(counts_after)

    return balanced_adata

def main(config_path: str):
    """
    Main function to run class balancing based on configuration from file.

    Parameters:
    - config_path (str): Path to configuration YAML file.
    """
    config = load_config(config_path)

    # Load AnnData object specified in config
    adata = sc.read_h5ad(config['data']['input_file'])

    # Extract balancing parameters from config
    target_key = config['balance']['target_key']
    mode = config['balance'].get('mode', 'min_class')
    target_count = config['balance'].get('target_count', None)
    random_state = config['balance'].get('random_state', 42)

    # Perform balancing
    balanced_adata = balance_classes(
        adata,
        target_key=target_key,
        mode=mode,
        target_count=target_count,
        random_state=random_state,
    )

    # Save balanced dataset
    balanced_adata.write_h5ad(config['data']['output_file'])
    print(f"Balanced dataset saved to {config['data']['output_file']}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Balance perturbation classes in single-cell dataset")
    parser.add_argument("--config", required=True, help="Path to YAML config file")
    args = parser.parse_args()
    main(args.config)
