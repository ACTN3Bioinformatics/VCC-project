"""
Class balancing rules
Balance perturbation classes to prevent model bias
"""

rule balance_classes:
    """
    Balance perturbation classes via downsampling
    """
    input:
        "results/{dataset}/scaled.h5ad"
    output:
        h5ad = "results/{dataset}/balanced.h5ad",
        report = "reports/{dataset}/balance_report.txt"
    params:
        target_cells = lambda wildcards: config['datasets'][wildcards.dataset].get("target_cells_per_perturbation", 100),
        min_cells = lambda wildcards: config['datasets'][wildcards.dataset].get("min_cells_per_perturbation", 10),
        perturbation_key = lambda wildcards: config['datasets'][wildcards.dataset].get("perturbation_key", "target_gene"),
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