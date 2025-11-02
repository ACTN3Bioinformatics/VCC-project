"""
Data splitting rules
"""

rule split_data:
    input:
        h5ad = lambda wildcards: (
            f"{{results}}/{{dataset}}/integrated.h5ad" 
            if datasets_config[wildcards.dataset].get("batch_correction", False)
            else f"{{results}}/{{dataset}}/balanced.h5ad"
        )
    output:
        train = "{results}/{dataset}/splits/train.h5ad",
        val = "{results}/{dataset}/splits/val.h5ad",
        test = "{results}/{dataset}/splits/test.h5ad",
        split_info = "{results}/{dataset}/splits/split_info.json"
    params:
        train_ratio = 0.7,
        val_ratio = 0.15,
        test_ratio = 0.15,
        gene_split = True,
        random_state = 42
    log:
        "logs/split/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 2
    resources:
        mem_mb = 8000
    script:
        "../../scripts/split_data.py"