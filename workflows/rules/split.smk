"""
Data splitting rules for train/val/test
"""

rule split_data:
    """
    Create train/validation/test splits
    Leave-genes-out strategy for CV
    """
    input:
        h5ad = lambda wildcards: (
            f"results/{wildcards.dataset}/integrated.h5ad" 
            if config['datasets'][wildcards.dataset].get("batch_correction", False)
            else f"results/{wildcards.dataset}/balanced.h5ad"
        )
    output:
        train = "results/{dataset}/splits/train.h5ad",
        val = "results/{dataset}/splits/val.h5ad",
        test = "results/{dataset}/splits/test.h5ad",
        split_info = "results/{dataset}/splits/split_info.json"
    params:
        train_ratio = 0.7,
        val_ratio = 0.15,
        test_ratio = 0.15,
        gene_split = True,  # Leave-genes-out CV
        random_state = 42,
        perturbation_key = lambda wildcards: config['datasets'][wildcards.dataset].get("perturbation_key", "target_gene")
    log:
        "logs/split/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../../scripts/split_data.py"