"""
Class balancing rules
"""

rule balance_classes:
    input:
        "{results}/{dataset}/scaled.h5ad"
    output:
        h5ad = "{results}/{dataset}/balanced.h5ad",
        report = "{reports}/{dataset}/balance_report.txt"
    params:
        target_cells = lambda wildcards: datasets_config[wildcards.dataset].get("target_cells_per_perturbation", 100),
        min_cells = lambda wildcards: datasets_config[wildcards.dataset].get("min_cells_per_perturbation", 10),
        perturbation_key = lambda wildcards: datasets_config[wildcards.dataset].get("perturbation_key", "target_gene"),
        random_state = 42
    log:
        "logs/balance/{dataset}.log"
    conda:
        "../../environment.yml"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 30
    script:
        "../../scripts/balance.py"